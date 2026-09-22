# Decision log

Chronological, with who decided and on what basis. Decisions that were *mine* rather than the
team's are marked, so they can be overridden.

## Scope and product

**D1 · Indication and MoA — atopic dermatitis, antagonise type 2 inflammatory signalling.**
Delegated to me by the team ("use your best judgment"). Basis: it has a verified CELLxGENE
catalog entry so Census coverage is plausible, a dense Open Targets evidence base (3,207
associated targets), and its biology is dominated by secreted and cell-surface proteins, so a
binder modality is genuinely appropriate rather than forced onto an intracellular target. NSCLC
chosen as the second indication for the generality test because it is immunologically distinct
and produces a disjoint target pool. **The target was not chosen by me** — the agent discovers
it from Open Targets at runtime.

**D2 · Modality — de novo minibinder, human.** Follows from the brief. Recorded because it is
what the accessibility screen is screening *for*: a soluble binder reaches secreted proteins and
extracellular domains and nothing else.

**D3 · Three-way tie is broken by data, not preference.** Sara's call, 12:50. The rubric
returned IL4R, TNFSF4 and IL31RA all at 0.613 with `cellular_context` unknown for all seven
ranked targets. Rather than break it on structural convenience, the decision was to wait for
CELLxGENE donor-aware expression and let the evidence separate them. Consequence, accepted
knowingly: design is the longest pole, so this puts the one paid GPU run behind the Census lane.
Contingency is mine to trigger — if Census has not resolved it in time to finish a real design
batch before the 15:45 freeze, I escalate back to Sara with a feasibility-based choice.

## Architecture

**D4 · Empty repo means build the foundation.** `lnairGT/Boston-CompBio-Hackathon` cloned with
zero commits. Spec section 5 prescribes a small Python app with a server-rendered interface and
typed records for exactly this case, so there was no existing stack to inherit and no teammate
work to preserve.

**D5 · Seven lanes with hard file ownership.** Each lane creates and edits only its own paths,
reads everything, and never touches `contracts.py` / `config.py` / `pipeline/` / `agent/`.
Lanes hand off through typed records rather than calling into each other's modules, so a lane
that has not landed cannot block one that has. Recorded in `docs/workstreams.md`.

**D6 · Capability degradation instead of stubs.** When a lane's adapter is absent its tools
return `{"status": "capability_unavailable"}` — an honest result the agent branches on — rather
than raising, and emphatically rather than returning a plausible fake value. This is what makes
`embedding_unavailable` and `source_unavailable` real outcomes instead of error paths.

**D7 · The LLM is injected, not imported.** `HostLLMClient` runs the loop inside Claude Science,
`AnthropicClient` runs it from a Modal container against the API, and a scripted client runs it
in tests with no network. Same loop in all three.

## Scientific method

**D8 · Accessibility is decided on topology, not location strings.** Changed after the first
live run: a keyword screen over Open Targets `subcellularLocations` passed PDE4A, PDE4B and
PDE4D as binder-accessible because they are annotated "plasma membrane" while being cytosolic
enzymes docked on the inner face. The screen now requires a UniProt signal peptide without a
transmembrane segment (secreted) or an explicit `Extracellular` topological domain (ectodomain).

**D9 · Minimum epitope surface of 40 residues.** HRH1 passes the topology test with a 29-residue
extracellular N-terminus — real, but not a design surface for a minibinder. Targets whose largest
engageable region is below the threshold are excluded with that specific reason rather than
silently passing. The threshold is a judgement call and is a named constant
(`MIN_EPITOPE_RESIDUES`) so it can be argued with.

**D10 · Direction of modulation is verified separately, and precedent does not transfer.**
Association says nothing about whether inhibition or activation helps. Direction fit is rated
from annotated agent mechanisms, partitioned into those pursued **for this indication** and those
pursued for others. Cintredekin besudotox is an IL4R-directed agent at Phase 3 — for brain
cancer — and it is recorded as off-indication precedent, not as support for atopic dermatitis.

**D11 · Text-mined literature is capped at `moderate`.** Literature co-occurrence scales with how
much a gene has been written about. When it is the only datatype contributing to a criterion, the
rating cannot reach `strong`. Implements "literature volume alone should not determine the winner".

**D12 · Single-donor observations cannot be rated `strong`.** Thousands of cells from one donor are
not thousands of independent replicates, so the cellular-context criterion is capped at `weak`
when `n_donors < 2` regardless of the expression value.

**D13 · Rubric version bumped v1 → v2** when the datatype mapping was corrected and the literature
cap added. Rankings recorded under v1 are not comparable to v2, which is why `rubric_version`
travels on every `Assessment` and `ComparisonResult`.

**D14 · Unscreened ≠ not-excluded.** See `capabilities.md` bug 6. `pending_screen` exists because
a target with no adverse evidence *yet* must not be presented as a finalist.

**D15 · The model may not declare a completed candidate report.** `AgentLoop._reconcile` downgrades
a claimed `completed_candidate_report` to `evidence_report_only` unless a design run actually
succeeded and produced artifact references. The outcome is checked against what happened, not
against what the model asserts.
