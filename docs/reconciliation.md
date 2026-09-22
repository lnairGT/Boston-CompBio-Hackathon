# Reconciling the agent work with `indication2binder`

22 September 2026 · written by the integration lane after reading `main` at `489d289`.

## What happened, plainly

I cloned this repository at about 11:45 and it had **zero commits**. I recorded "the repo is
empty, so we are building the foundation" in `docs/workstreams.md` and `docs/decision_log.md`,
and built on that assumption for roughly ninety minutes without re-checking.

Lakshmi seeded the repository at 12:24 and merged four pull requests between 13:03 and 13:07.
By the time I pushed `lane/1-integration` at 13:17, `main` already carried a complete staged
pipeline solving the same problem.

The two trees have **no common ancestor**. That is my error, not a merge accident: I should
have re-fetched before pushing, and certainly before telling the team the repo was empty.

**`main` is the foundation.** What follows re-homes the agent work onto it.

## What `indication2binder` already is

A staged CLI pipeline, `src/ind2b/`, stages 0-5 all implemented: resolve an indication to an
ontology term, assemble evidence layers, score and rank, keep only targets with an
experimentally determined complex, extract interface residues and hotspots, emit a BindCraft2
specification. It deliberately **stops before running a design**, because that decision costs
GPU time and belongs to a human reading the epitope shortlist.

Infrastructure worth noting, because it is better than mine: every HTTP response is disk-cached
and keyed by method, URL, params and body, so a populated cache replays a whole run offline;
every tunable lives in `config.py` or a weights JSON rather than inline; and `ind2b schema-check`
introspects the live Open Targets schema and fails loudly when a queried field disappears.

## Independent convergence

We built these separately and arrived at the same place. That is mutual validation, and it is
the strongest evidence either of us has that these choices are right.

| Decision | `ind2b` | agent spine |
|---|---|---|
| Accessibility from **topology**, not location keywords | `sources/uniprot.py` | `adapters/structures.py` |
| Accessibility is a **filter**, not a ranking tiebreak | stage 1 gate | `pending_screen` + exclusion |
| Interface contact cutoff | **5.0 Å** | **5.0 Å** |
| Minimum usable target chain / epitope | **40 residues** | **40 residues** |
| Minimum ectodomain span | 30 residues | signal-peptide / topological-domain test |
| PDB author → UniProt numbering via SIFTS | `stage4_interface.auth_to_uniprot` | `structures.to_uniprot_numbering` |
| Antibody/Fab partners are **not** a usable interface | description regex, stage 3 | unmapped-chain detection |
| Weighted, versioned, re-runnable rubric | `DEFAULT_WEIGHTS` + weights JSON | `CRITERIA` + `RUBRIC_VERSION` |
| Rank ≠ causality; model score ≠ affinity | README "does not establish" | contract validators |
| Open Targets schema drift must be checked | `ind2b schema-check` | `verify_schema()` |
| Parent association already aggregates descendants — don't double-count | `association_scope` | ontology `descendant` never folded in |

Two of these deserve emphasis. Both of us independently rejected keyword matching on subcellular
location, and both of us independently excluded antibody complexes from epitope derivation. Those
are the two errors most likely to silently produce a binder aimed at the wrong place.

## Correction: `ind2b`'s partner filter is better than mine, not complementary to it

I first wrote that our two antibody filters were complementary and should be combined. Reading
`stage3_complexes.py` properly, that is wrong and worth correcting rather than quietly dropping.

- **Mine is a denylist.** A partner chain is rejected when it has no UniProt mapping, or is
  another copy of the target. Everything else is admitted, so it depends on me having
  enumerated the ways a chain can be wrong.
- **`ind2b`'s is an allowlist.** A chain becomes a `native_partner` only if its accession
  appears in that target's **Open Targets interaction-partner list**. An antibody, a nanobody,
  a DARPin, a crystallisation chaperone and an unnamed engineered binder all fail that test by
  construction, without anyone having to think of them in advance. The immunoglobulin and
  chaperone regexes are only used to *label* the rejected chains for the audit trail.

An allowlist keyed to known interactors is the stronger design and mine adds nothing to it.

> **Unreviewed.** This comparison is my own reading of her code. The verification lane declined to
> check it and said so explicitly, so it carries no independent confirmation — unlike the
> `target_chains[0]` finding below, where review corrected me.
`ind2b` also has an insight I did not: **one epitope per partner chain, not per entry**, because
a target contacting two partners at once yields a composite interface of several thousand square
ångström that no single mini-binder reproduces. That is a real constraint on what is designable
and my per-entry contact set quietly ignored it.

## One real gap, verified on live data

`stage4_interface.analyse_target` picks the target chain with:

```python
target_chain = target_chains[0]
```

When an entry contains several copies of the target this silently takes the first. Two entries
in our own candidate set do:

```
8K4Q  P24394 -> chains A, C   offsets {A: +25, C: +25}   agree -> picks A silently
3L5W  P35225 -> chains I, J   offsets {I: +33, J: +33}   agree -> picks I silently
1IAR  P24394 -> chain  B      offsets {B: +25}           single chain, no ambiguity
```

**Corrected after independent review — the verdict held, my mechanism did not.**

I first wrote that the offsets "happen to agree", so her output is correct by luck, and that a
disagreeing offset would shift every reported residue position. That is wrong.
`auth_to_uniprot` builds the map **per chain**, from that chain's own author-id list, for the
chain that is then analysed — so a disagreeing offset **could not** shift the reported positions.
Her code is *structurally* safe here, not accidentally safe. Describing a colleague's code as
luckily correct when it is correct by construction is a costly kind of wrong.

The gap at `target_chains[0]` is still real, but it is a different gap:

1. **Copy selection.** If the chosen copy is not the one a given partner engages,
   `interface_residues` returns empty and the record is rejected with a stated reason — a
   *visible false negative*, an epitope lost rather than a wrong epitope.
2. **Under-sampling**, the quieter case: a less complete or differently packed copy yields a
   narrower accepted epitope that still clears `MIN_BURIED_AREA`, and nothing flags it. This one
   is **reasoned from her gates, not demonstrated** — no real entry was found where copy `[0]`
   gives a materially narrower interface than its sibling.
3. **Unmapped chains fail safe**: an empty map yields `status='rejected'` rather than author
   positions presented as canonical.

Suggested fix, unchanged and small: record `target_chains_available` alongside `target_chain` so
a reader can see the copy choice was arbitrary. **No change to the numbering path is warranted.**

Suggested minimal change, in the spirit of the existing code: when `len(target_chains) > 1`,
build the numbering map for each, and either record `alternative_chains` with their offsets on
the epitope row, or reject with `ambiguous_target_chain` when the offsets disagree. A caller
should not be able to miss it.

## Disposition of every piece of agent work

**Adopt `ind2b`, drop mine** — duplicated, and theirs has better infrastructure:

| Mine | Superseded by |
|---|---|
| `adapters/open_targets.py` | `sources/opentargets.py` (disk-cached, `schema-check` CLI) |
| topology accessibility in `adapters/structures.py` | `sources/uniprot.py` |
| `interface_residues` contact counting | `stage4_interface.py` (real SASA via biopython/gemmi; mine is a dependency-free approximation) |
| `pipeline/ranking.py` rubric | `stage2_score.py` (weights in a file, re-rank without re-fetching) |
| `config.py` constants | `config.py` |

**Contribute into `ind2b` — genuinely additive, nothing there covers it:**

| Capability | Proposed home | Why it is new |
|---|---|---|
| Single-cell disease-context expression | `sources/cellxgene.py` + a stage 1 layer | `ind2b` has no single-cell layer at all. Donor-aware, disease-vs-reference. |
| Cross-vocabulary disease mapping (MONDO) | extend `stage0_resolve` | stage 0 resolves a term; this maps that term onto *another corpus's* vocabulary, classifying exact / xref / synonym / descendant / ancestor. Required to reach CELLxGENE at all. |
| Design **execution** | new `stage6_design` (opt-in) | `ind2b` stops at spec by design. Proven route on Modal: RFdiffusion → ProteinMPNN → Boltz-2 → metrics. **Must preserve the stage-5 stop**: off by default, explicit human approval, and the approval bound to the exact parameters. |
| Candidate triage | new `stage7_triage` | developability, interface-metric interpretation, paralogue epitope cross-reactivity |
| Wet-lab plan | new `stage8_experiments` | assays with controls and pre-registered go/no-go criteria |
| Pathway alternatives | `sources/reactome.py` + stage | other members of implicated pathways as intervention hypotheses |
| Agent loop + web UI | separate optional front-end over the same stage outputs | `ind2b` is a CLI; the agent is an alternative driver, not a replacement |

**Contribute as fixes to `ind2b`:**

1. **Chain multiplicity in `stage4_interface`** — the `target_chains[0]` gap above. This is the
   one code-level fix I am confident `ind2b` needs, and it is verified on a real entry.
2. **Bind approval to parameters.** If a stage-6 gate is added, the approval must carry a hash
   of the exact parameters and refuse when they change — approving a 4-design run and then
   setting it to 512 must not launch. Nested in-place mutation is the case that evades a naive
   check, because pydantic's `validate_assignment` does not fire on it.
4. **An empty result is not a negative result.** A source timeout, an unscreened target and a
   disease absent from a corpus must each be distinguishable from a measured zero.

## Things I got wrong that are worth recording

- **The repo was not empty.** Stated as fact for ninety minutes without re-checking.
- **I accused a sub-agent of string-matching** a disease label. I never observed it; I inferred
  it. Its conclusion was correct: `MONDO:0004980` has zero cells in all six Census builds
  probed. What I found was the CELLxGENE **Discover** corpus, which the Census lags by about six
  months. Two corpora, not one.
- **I mis-labelled a metric column** when relaying it between lanes, calling per-side interface
  area the total. Caught downstream only because that lane read the convention from each
  record's own method text rather than from my label.

## Suggested sequence

1. Open a PR from `lane/1-integration` **for review, not merge** — it is a parallel tree and
   should not land on `main` as-is.
2. Land the two filter improvements against `main` first; they are small and independently useful.
3. Port the single-cell layer as `sources/cellxgene.py`, the highest-value additive piece.
4. Decide as a team whether `ind2b` should ever execute a design. If yes, stage 6 goes in
   behind the gate described above. If no, the Modal route stays a separate tool that consumes
   stage 5 specs — which is arguably the cleaner boundary anyway.
