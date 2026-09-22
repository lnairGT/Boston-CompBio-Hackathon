# Three-minute demo script — evidence to binder

Boston Computational Biology Hackathon · 22 September 2026 · Lane 5 (experience)

**Total 3:00.** Times are cumulative. Everything in *italics* is what you say; everything
in `code` is what you click or type. Two rules for the whole demo: **read the origin badge
out loud the first time one appears**, and **never say "affinity" about a model score**.

---

## Before you start (not on the clock)

```bash
# from the repo root
pip install fastapi jinja2 'pydantic>=2' uvicorn python-multipart
uvicorn e2b.ui.app:app --port 8000
```

Open `http://127.0.0.1:8000`. Have **two browser tabs** ready:

| Tab | URL | Why |
|---|---|---|
| A | `http://127.0.0.1:8000` | the live brief you will type into |
| B | `http://127.0.0.1:8000/fixtures/walkthrough` | the guaranteed complete flow |

Tab B is the safety net. If the live run stalls, switch to it and keep talking — it is
labelled `fixture` on every card, so showing it is honest, not a fudge.

---

## 0:00 – 0:30 · Brief — free text in, resolved identity out

Tab A. Type into the two boxes (do **not** use a preset — typing proves nothing is
hardcoded):

- Indication: `atopic dermatitis`
- Desired effect: `suppress type 2 inflammatory signalling`
- Direction: `antagonise`   Modality: `de novo minibinder`

Click `Resolve this brief`.

> *"No dropdown, no allowlist. We hand the free text to the resolver, and it comes back
> with candidate ontology terms — atopic dermatitis actually resolves to MONDO_0004980,
> 'atopic eczema'. We show the ambiguity instead of silently picking one."*

Point at the **Agent stage** strip: stage name, timestamp, and one line of what happens
next.

> *"There is no progress bar on this page anywhere. The agent reports a stage, not a
> percentage, so we show the stage and the time. A bar would be a number we invented."*

If several terms are offered, select one and click `Use the selected term`.

## 0:30 – 1:15 · Compare — evidence, gaps, and feasibility kept apart

Scroll to **Compare** (or switch to tab B if the live run has not reached it).

> *"Up to three targets. Each card has the four rubric criteria with a rating and the
> evidence behind it."*

Point at a criterion rated **no evidence retrieved**.

> *"That one is a gap. It says what we looked for and did not find. It scores zero
> weight — missing evidence never becomes a favourable score."*

Now point at the **Design feasibility** block, below the dashed rule.

> *"This is deliberately separate from biological support. A target can have excellent
> genetics and still be undesignable for this modality — the third card is annotated on
> the cytoplasmic face, so a soluble binder cannot reach it. We do not let feasibility
> quietly inflate a biology score."*

If the run shows ties or *insufficient evidence*, say so explicitly:

> *"When the rubric cannot separate two targets, it says they are tied. It does not
> invent a tiebreaker."*

## 1:15 – 1:45 · Inspect — the evidence trail

Scroll to **Inspect**.

> *"Every claim on the cards is here with its source, its locator and its disease scope —
> whether the item is about this disease or was propagated from a related term."*

Open one **Provenance** disclosure.

> *"Provenance is one click away, not in your face: endpoint, query, release, retrieval
> time."*

Point at the **Source status** table (present in the `timed_out` fixture):

> *"And when a source times out we record it as a source failure. A timeout is not
> evidence that the data is absent, and this table says so in those words."*

## 1:45 – 2:20 · Design — the gate

Scroll to **Design**. Show the binding site, hotspots and the MoA rationale.

> *"Selected target, chosen site, and why engaging here would produce the effect asked
> for. The launch button is greyed out."*

Tick `I have read the target, site and rationale above and confirm these inputs.`

> *"It only enables once a human confirms the inputs — and the request itself has to be
> approved and complete: structure, hotspots, rationale, evidence links, sequence
> checksum. The button carries a one-shot token, so a double-click or a refresh cannot
> start a second paid run."*

Click `Launch design run`. Read the **Job status** panel aloud: state, stage, route,
timestamps.

> *"Real job state from the executor. If nothing is attached, it says so and spends
> nothing — which is what you are seeing right now."*

## 2:20 – 3:00 · Results — and what we refuse to claim

Switch to tab B (`walkthrough`) if the live job has not returned.

> *"Candidate sequence, structure reference, and the model scores."*

Point at the metric table and then at the null metric.

> *"These are model scores. Interface PAE, pLDDT. None of them is a measured binding
> affinity and we never print one under a name like Kd. The third row is null, not zero —
> it was not computed, and zero would be a lie."*

Point at the fixture badges.

> *"Everything on this tab is labelled fixture. It is here to prove the interface, and it
> is not a biological result — we would rather show you a labelled stand-in than an
> unlabelled guess."*

Close on the next-test list.

> *"And the honest end of the pipeline: the next experimental test. Biolayer
> interferometry for a real Kd, and a cell-based assay to show the binder antagonises
> rather than merely binds. The computer has proposed a hypothesis; it has not validated
> one."*

---

## If the design job is still running at demo time

Do **not** wait on it and do not narrate over a spinner.

1. Say: *"the design job is still running — here is its real state"*, and read the Job
   status panel: stage name and last-update timestamp.
2. Switch to tab B and finish the Results section on the `walkthrough` fixture, saying
   out loud that it is a fixture.
3. If asked when the real one lands, quote the timestamp on the panel and say you do not
   know — the executor does not report an ETA and inventing one would be the exact
   failure mode this UI is built to avoid.

## If a source is down

Open `http://127.0.0.1:8000/fixtures/timed_out` and use it deliberately:

> *"This is the state we care most about getting right. The source timed out. We show a
> source failure with the query we issued, and we explicitly do not convert it into
> 'no evidence found'."*

## Fixture URLs (all labelled, all reachable)

| State | URL |
|---|---|
| complete run | `/fixtures/walkthrough` |
| replayed cached result | `/fixtures/cached` |
| ambiguous indication | `/fixtures/ambiguous` |
| job pending | `/fixtures/pending` |
| nothing retrieved | `/fixtures/no_data` |
| job failed | `/fixtures/failed` |
| retrieval + job timeout | `/fixtures/timed_out` |

## Serving it on Modal

`e2b/modal_app/web.py` defines the ASGI entrypoint (`modal deploy e2b/modal_app/web.py`).
It needs the `e2b` package on the image, `fastapi`/`jinja2`/`pydantic`, and **no secret**
— the UI holds no credential. Pin to one container so the in-memory run store stays
coherent for one audience. As of this writing the entrypoint has been import-checked but
**not deployed**; do not quote a Modal URL in the demo unless a deploy has actually
printed one.
