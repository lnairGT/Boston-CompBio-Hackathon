"""indication2binder: indication -> ranked targets -> epitope -> BindCraft2 spec.

The pipeline is deliberately staged. Each stage reads the previous stage's
JSON/CSV output from the run directory and writes its own, so any stage can be
re-run in isolation without repeating network work.

Stages
------
0 ``stage0_resolve``    indication string -> disease ontology id + subtree
1 ``stage1_evidence``   disease id -> per-target evidence payloads
2 ``stage2_score``      evidence -> ranked target table
3 ``stage3_complexes``  ranked targets -> qualifying experimental complexes
4 ``stage4_interface``  complexes -> target-side interface residues
5 ``stage5_specs``      epitopes -> BindCraft2 target specifications

Nothing in this package submits a design job. Stage 5 writes specification
files and stops; running BindCraft2 is a separate, explicit step.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
