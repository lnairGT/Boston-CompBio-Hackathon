"""The bounded Claude tool-use loop.

This is an agent-directed investigation, not a fixed pipeline with a chat summary bolted on.
The loop is: observe persisted state -> select the next evidence-gathering or design action
-> execute an allowed typed tool -> validate the result -> update the evidence graph and
uncertainty -> continue, revise, ask a targeted question, or stop.

What the model does and does not do:

* It **chooses** which tool to call next, and it **explains** results using retrieved
  evidence IDs.
* It does **not** compute scores, invent rubric weights, or decide whether a design job may
  launch. Those are deterministic services (``pipeline.ranking``, ``DesignRequest.assert_launchable``).

The LLM is injected rather than imported, so the same loop runs three ways: inside Claude
Science via ``host.llm``, against the Anthropic API from a Modal container, and in tests
against a scripted client with no network at all.
"""

from __future__ import annotations

import json
import os
import uuid
from typing import Any, Protocol

from ..contracts import AgentOutcome, Budgets, ResearchBrief, utcnow
from .state import RunStore
from .tools import ToolRegistry

MAX_STEPS_HARD_CEILING = 60

SYSTEM_PROMPT = """You are the investigation agent in an evidence-to-binder workflow for drug \
discovery. A scientist gives you an indication and a desired mechanism of action; you decide \
which target to recommend, why, and whether a protein-binder design job is justified.

You work by calling typed tools. Each turn: look at the state summary, decide the single most \
useful next action, call one tool, read the result, and revise your plan. Give a one-sentence \
reason for every call.

Rules that are not negotiable:

1. Resolve the indication first. Never assume an ontology ID string is valid.
2. Retrieve a bounded candidate pool. Never claim exhaustive discovery.
3. Screen modality suitability EARLY with assess_structure, before spending deep-evidence \
calls. A soluble binder can only engage secreted proteins or extracellular domains. If a \
strongly-associated target is intracellular, say so plainly, keep its biology visible, and \
investigate accessible candidates instead of pretending it is reachable.
4. Association is not causation, not a probability of success, and says nothing about \
whether inhibition or activation helps. Verify direction separately and keep unknown \
direction explicit.
5. You never compute a score. compare_targets runs a deterministic versioned rubric; your \
job is to explain its output, name the evidence IDs behind it, and identify which missing \
evidence would change the decision. If it reports a tie, report the tie. Do not invent a \
tiebreaker.
6. A tool result of capability_unavailable or source_unavailable means the CAPABILITY or \
SOURCE is missing. It never means the biological thing is absent. Record the gap and \
qualify the report.
7. Never manufacture expression values, citations, sequences or design metrics. If you \
cannot support a claim with a retrieved evidence ID, do not make the claim.
8. Stop when you have enough to recommend a target and a binding region, or when budget is \
exhausted, or when a specific limitation blocks progress. Returning an honest partial result \
is a success; a fabricated complete one is a failure.

When you are done investigating, call no further tools and reply with a final report as \
JSON in a ```json fenced block with keys: outcome (one of completed_candidate_report, \
evidence_report_only, needs_clarification, insufficient_evidence, no_binder_feasible_target, \
source_unavailable, budget_exhausted, design_failed, embedding_unavailable), \
recommended_target_id, recommended_gene_symbol, rationale, binding_region_proposal, \
evidence_ids, remaining_uncertainty, next_action."""


class LLMClient(Protocol):
    def complete(self, messages: list[dict[str, Any]], tools: list[dict[str, Any]],
                 system: str) -> dict[str, Any]:
        """Return {"content": [blocks], "stop_reason": str}. Blocks are Anthropic-shaped."""
        ...


class HostLLMClient:
    """Runs the loop inside Claude Science using the in-kernel ``host.llm`` bridge."""

    def __init__(self, host: Any, model: str | None = None):
        self.host = host
        self.model = model or host.reasoning_model()

    def complete(self, messages, tools, system):
        resp = self.host.llm({
            "messages": messages, "tools": tools, "system": system,
            "model": self.model, "max_tokens": 4096,
        })
        content = resp.get("content")
        if content is None:
            content = [{"type": "text", "text": resp.get("text", "")}]
        # Two shape fixes for the host bridge, kept here so the loop stays provider-agnostic.
        #
        # 1. It omits the `id` on tool_use blocks, but a tool_result must carry a matching
        #    tool_use_id or the following turn is rejected. Mint stable ids.
        # 2. It strips the cryptographic signature from `thinking` blocks, and the API
        #    rejects an unsigned thinking block echoed back into the conversation. Drop them
        #    entirely -- which is also what we want on the record: the spec asks for brief
        #    action rationales and source references, not private chain-of-thought.
        cleaned = []
        for block in content:
            if block.get("type") in ("thinking", "redacted_thinking"):
                continue
            if block.get("type") == "tool_use" and not block.get("id"):
                block["id"] = f"toolu_{uuid.uuid4().hex[:20]}"
            cleaned.append(block)
        return {"content": cleaned, "stop_reason": resp.get("stop_reason", "end_turn")}


class AnthropicClient:
    """Runs the loop from a Modal container against the Anthropic API.

    The key comes from the environment, which on Modal is populated from a Secret. It is
    never passed through a browser, logged, or written into a run manifest.
    """

    def __init__(self, model: str = "claude-sonnet-4-5", api_key: str | None = None):
        import anthropic  # imported lazily so the package is optional
        key = api_key or os.environ.get("ANTHROPIC_API_KEY")
        if not key:
            raise RuntimeError("ANTHROPIC_API_KEY is not set; provide it via a Modal Secret.")
        self.client = anthropic.Anthropic(api_key=key)
        self.model = model

    def complete(self, messages, tools, system):
        resp = self.client.messages.create(
            model=self.model, max_tokens=4096, system=system, tools=tools, messages=messages)
        return {"content": [b.model_dump() for b in resp.content], "stop_reason": resp.stop_reason}


class AgentLoop:
    def __init__(self, llm: LLMClient, store: RunStore, registry: ToolRegistry | None = None):
        self.llm = llm
        self.store = store
        self.registry = registry or ToolRegistry(store=store)
        self.transcript: list[dict[str, Any]] = []

    # -- observation ----------------------------------------------------------------

    def observe(self) -> str:
        """The state summary the model reasons over. Persisted state, not conversation history."""
        s = self.store
        caps = self.registry.capabilities()
        lines = [
            f"BRIEF: indication='{s.brief.raw_indication}' desired_effect='{s.brief.desired_effect}' "
            f"moa='{s.brief.moa}' modality='{s.brief.modality}' species={s.brief.species}",
            f"RESOLVED: {s.brief.resolved_disease_id or 'NOT YET RESOLVED'} "
            f"({s.brief.resolved_disease_label or '-'})",
            f"STAGE: {s.run.current_stage}",
            f"BUDGET: {s.run.budgets.tool_decisions_used}/{s.run.budgets.max_tool_decisions} tool "
            f"decisions used; paid design batches "
            f"{s.run.budgets.paid_design_batches_used}/{s.run.budgets.max_paid_design_batches}",
            f"POOL: {len(s.targets)} target(s) retrieved"
            + (f" of {s.pool_scope.get('total_associated_in_release')} associated in release" if s.pool_scope else ""),
        ]
        if s.structures:
            acc_ok = [a for a, v in s.structures.items() if v["accessibility"]["accessible"]]
            lines.append(
                f"ACCESSIBILITY SCREENED: {len(s.structures)} protein(s); {len(acc_ok)} accessible to a "
                f"soluble binder, {len(s.structures) - len(acc_ok)} not")
        if s.evidence:
            lines.append("DEEP EVIDENCE: " + ", ".join(
                f"{s.targets[t].gene_symbol}({len(e)} items)" for t, e in s.evidence.items() if t in s.targets))
        else:
            lines.append("DEEP EVIDENCE: none retrieved yet")
        lines.append(f"CELL SUMMARIES: {len(s.cell_summaries)}")
        if s.comparison:
            lines.append(
                f"COMPARISON ({s.comparison.rubric_version}): order="
                f"{[s.targets[t].gene_symbol for t in s.comparison.ordered_target_ids[:5] if t in s.targets]} "
                f"ties={[[s.targets[t].gene_symbol for t in g if t in s.targets] for g in s.comparison.ties]} "
                f"excluded={len(s.comparison.excluded)} insufficient={s.comparison.insufficient_evidence}")
            gaps = {t: c for t, c in s.comparison.unresolved_criteria.items()
                    if t in s.comparison.ordered_target_ids[:3]}
            if gaps:
                lines.append("UNRESOLVED CRITERIA on shortlist: " + json.dumps(
                    {s.targets[t].gene_symbol: c for t, c in gaps.items() if t in s.targets}))
        unavailable = [k for k, v in caps.items() if isinstance(v, dict) and not v.get("available", True)]
        if unavailable:
            lines.append("CAPABILITIES UNAVAILABLE: " + ", ".join(
                f"{k} ({caps[k]['detail']})" for k in unavailable))
        if s.design_runs:
            lines.append("DESIGN RUNS: " + ", ".join(
                f"{r.run_id}={r.status.value}" for r in s.design_runs.values()))
        return "\n".join(lines)

    # -- the loop -------------------------------------------------------------------

    def run(self, max_steps: int | None = None, verbose: bool = True) -> dict[str, Any]:
        budget = self.store.run.budgets
        limit = min(max_steps or budget.max_tool_decisions, MAX_STEPS_HARD_CEILING)

        messages: list[dict[str, Any]] = [{
            "role": "user",
            "content": (
                "Investigate this brief and recommend a target and binding region.\n\n"
                f"=== PERSISTED STATE ===\n{self.observe()}\n\n"
                "Begin. Call one tool at a time and give a one-sentence reason for each."
            ),
        }]
        tools = self.registry.anthropic_tools()
        final_report: dict[str, Any] | None = None
        step = 0

        while step < limit:
            if budget.exhausted():
                self.store.finish("budget_exhausted",
                                  f"Tool-decision budget of {budget.max_tool_decisions} spent.")
                return self._result("budget_exhausted", final_report)

            resp = self.llm.complete(messages, tools, SYSTEM_PROMPT)
            content = resp.get("content") or []
            messages.append({"role": "assistant", "content": content})

            tool_uses = [b for b in content if b.get("type") == "tool_use"]
            texts = [b.get("text", "") for b in content if b.get("type") == "text"]

            if texts and verbose:
                for t in texts:
                    if t.strip():
                        self.transcript.append({"step": step, "kind": "reasoning", "text": t.strip()[:1200]})

            if not tool_uses:
                final_report = _parse_report("\n".join(texts))
                break

            results = []
            for tu in tool_uses:
                step += 1
                name, args = tu["name"], (tu.get("input") or {})
                reason = _nearest_reason(texts) or f"step {step}: {name}"
                result, action = self.registry.call(name, args, reason=reason, step=step)
                self.store.log_action(action)
                self.store.set_stage(name)
                if verbose:
                    self.transcript.append({
                        "step": step, "kind": "tool_call", "tool": name,
                        "arguments": action.arguments, "reason": reason,
                        "status": result.get("status"), "error": action.error,
                        "elapsed_ms": action.elapsed_ms,
                    })
                results.append({
                    "type": "tool_result", "tool_use_id": tu["id"],
                    "content": json.dumps(_trim(result), default=str)[:14000],
                    "is_error": result.get("status") in ("error",),
                })

            results.append({
                "type": "text",
                "text": f"=== PERSISTED STATE (refreshed) ===\n{self.observe()}",
            })
            messages.append({"role": "user", "content": results})

        outcome: AgentOutcome = (final_report or {}).get("outcome") or "evidence_report_only"
        valid = set(AgentOutcome.__args__)  # type: ignore[attr-defined]
        if outcome not in valid:
            outcome = "evidence_report_only"

        # The model proposes an outcome; the store validates it against what actually happened.
        outcome = self._reconcile(outcome)
        self.store.finish(outcome, (final_report or {}).get("rationale"))
        return self._result(outcome, final_report)

    def _reconcile(self, claimed: AgentOutcome) -> AgentOutcome:
        """Never let the model claim a completed candidate report without real candidates."""
        s = self.store
        has_real_candidates = any(
            r.status.value == "succeeded" and r.artifact_refs for r in s.design_runs.values())
        if claimed == "completed_candidate_report" and not has_real_candidates:
            return "evidence_report_only"
        if not s.targets and claimed not in ("needs_clarification", "source_unavailable"):
            return "source_unavailable"
        if s.comparison and s.comparison.insufficient_evidence and claimed == "evidence_report_only":
            return "insufficient_evidence"
        if s.comparison is not None and not s.comparison.ordered_target_ids and s.comparison.excluded:
            return "no_binder_feasible_target"
        return claimed

    def _result(self, outcome: str, report: dict[str, Any] | None) -> dict[str, Any]:
        paths = self.store.save()
        return {
            "outcome": outcome,
            "report": report,
            "manifest": self.store.manifest(),
            "transcript": self.transcript,
            "artifact_paths": paths,
        }


def _parse_report(text: str) -> dict[str, Any] | None:
    if "```json" in text:
        blob = text.split("```json", 1)[1].split("```", 1)[0]
    elif "{" in text:
        blob = text[text.index("{"): text.rindex("}") + 1]
    else:
        return {"rationale": text.strip()[:4000]} if text.strip() else None
    try:
        return json.loads(blob)
    except json.JSONDecodeError:
        return {"rationale": text.strip()[:4000]}


def _nearest_reason(texts: list[str]) -> str | None:
    for t in reversed(texts):
        cleaned = " ".join(t.split())
        if cleaned:
            return cleaned[:300]
    return None


def _trim(result: dict[str, Any]) -> dict[str, Any]:
    """Keep full matrices and long vectors out of the model's context."""
    out: dict[str, Any] = {}
    for k, v in result.items():
        if isinstance(v, list) and len(v) > 40:
            out[k] = v[:40] + [f"<{len(v) - 40} more omitted; reference the artifact>"]
        else:
            out[k] = v
    return out


def new_brief(raw_indication: str, *, desired_effect: str | None = None, moa: str | None = None,
              modality: str = "de novo minibinder", species: str = "homo_sapiens") -> ResearchBrief:
    return ResearchBrief(
        brief_id=f"brief-{uuid.uuid4().hex[:10]}", raw_indication=raw_indication,
        desired_effect=desired_effect, moa=moa, modality=modality, species=species, status="pending")


def start_run(raw_indication: str, *, llm: LLMClient, desired_effect: str | None = None,
              moa: str | None = None, modality: str = "de novo minibinder",
              budgets: Budgets | None = None, root: str = "runs") -> tuple[AgentLoop, RunStore]:
    brief = new_brief(raw_indication, desired_effect=desired_effect, moa=moa, modality=modality)
    store = RunStore(run_id=f"run-{uuid.uuid4().hex[:10]}", brief=brief, root=root, budgets=budgets)
    return AgentLoop(llm=llm, store=store), store
