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

## The two antibody filters are complementary — combine them

This is the clearest immediate improvement, and neither approach is sufficient alone.

- **`ind2b`** matches the *entity description* against an immunoglobulin regex (`fab`, `scfv`,
  `vhh`, `darpin`, …) plus a chaperone regex. Catches an antibody chain that **has** a UniProt
  accession. Misses a novel binder scaffold whose description uses none of those words.
- **agent spine** flags any partner chain with **no UniProt mapping**, plus chains that are
  another copy of the target. Catches unnamed and engineered binders. Misses a UniProt-mapped
  antibody, and false-positives on any legitimately unmapped partner.

Running both, and requiring a partner to pass each, is strictly stronger than either.

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

1. **Unmapped-chain partner check** added alongside the antibody regex (above).
2. **Chain multiplicity signalling.** 8K4Q maps `P24394` to chains A *and* C. Returning one
   chain without flagging the other lets a caller design against an arbitrary copy.
3. **Bind approval to parameters.** If a stage-6 gate is added, the approval must carry a hash
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
