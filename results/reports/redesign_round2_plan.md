# Round-2 Focused Redesign Plan

Phase 7A prepares deterministic ProteinMPNN-ready campaign inputs for a narrow redesign round centered on the three metadynamics-validated seed scaffolds.
No ProteinMPNN, LigandMPNN, Rosetta, MD, QM, or quantum execution is started here.

## Why QM Is Not Started Yet

- Phase 6C1 still classifies all three selected seed scaffolds as `candidate_keep_for_qm=FALSE`, so there is no metadynamics-backed signal yet to escalate any design directly into QM.
- The immediate failure mode is shared across the seeds: `Dy` remains `retained_bound`, while both `Al` and `Fe` remain `persistent_capture`, so the next efficient step is to narrow sequence space around the tolerated local neighborhoods before spending on QM.
- This phase therefore only exports a deterministic round-2 redesign plan and leaves all expensive downstream methods for later phases.

## Why These Three Seed Scaffolds

| seed | topology | panel_role | Dy | Al | Fe | seed_mutations | why_kept |
| --- | --- | --- | --- | --- | --- | --- | --- |
| am1_mex_ss_only_u02 | am1_monomer | am1_mex_designed_candidate | retained_bound | persistent_capture | persistent_capture | 18,65,69,100 | Retained as the best AM1/Mex design by integrated Phase 5B ranking; LigandMPNN ligand_confidence=0.3307 and Rosetta within-topology rank=1. |
| hans_pocket_ss_only_u02 | hans_monomer | hans_pocket_designed_candidate | retained_bound | persistent_capture | persistent_capture | 18,83,93,98 | Retained as the best Hans pocket-focused design; LigandMPNN ligand_confidence=0.4251 and Rosetta within-topology rank=1. |
| hans_interface_ss_plus_if_u04 | hans_interface_multichain | hans_interface_aware_candidate | retained_bound | persistent_capture | persistent_capture | 18,38,93,98 | Retained as the best Hans interface-aware design to preserve the multichain hypothesis; LigandMPNN ligand_confidence=0.3853 and Rosetta within-topology rank=1. |

These three seeds span the AM1/Mex monomer, the Hans pocket-focused monomer, and the Hans interface-aware multichain context while preserving the exact metadynamics phenotype that motivates round 2: Dy retention without relief of persistent Al/Fe capture.

## Why The Scope Is Intentionally Narrow

- `fixed_first_shell` and `protected_positions` stay fixed for every seed, so the Dy-retaining first-shell architecture and numbering-protected positions are not reopened.
- Redesignable sites must already be marked `mutable_second_sphere` or `mutable_interface`, must remain inside mature AM1 numbering, and must lie within +/-3 canonical positions of a seed mutation or of a known interface-aware mutable site.
- This keeps round 2 focused on local neighborhoods most likely to tune competitor capture without destabilizing the Dy-compatible scaffold cores.

## Planned Round-2 Campaigns

| seed | designable_count | designable_canonical_positions | designable_residue_ids |
| --- | --- | --- | --- |
| am1_mex_ss_only_u02 | 18 | 18,24,27,38,39,40,44,65,69,73,75,76,79,80,83,97,98,100 | A38,A44,A47,A56,A57,A58,A62,A83,A87,A91,A93,A94,A97,A98,A101,A115,A116,A118 |
| hans_pocket_ss_only_u02 | 17 | 18,24,27,38,39,40,44,73,75,76,79,80,83,93,97,98,100 | A35,A41,A44,A55,A56,A57,A61,A90,A92,A93,A96,A97,A100,A110,A114,A115,A117 |
| hans_interface_ss_plus_if_u04 | 17 | 18,24,27,38,39,40,44,73,75,76,79,80,83,93,97,98,100 | A35,A41,A44,A55,A56,A57,A61,A90,A92,A93,A96,A97,A100,A110,A114,A115,A117 |

## am1_mex_ss_only_u02

- Starting structure: `results/ligandmpnn_smoke/am1_mex_ss_only_u02/packed/am1_mex_ss_only_u02_packed_2_1.pdb`
- Seed mutation canonical positions: `18,65,69,100`
- Redesignable canonical positions: `18,24,27,38,39,40,44,65,69,73,75,76,79,80,83,97,98,100`

| residue_id | canonical | am1 | seed_aa | mutable_second_sphere | mutable_interface | reason |
| --- | --- | --- | --- | --- | --- | --- |
| A38 | 18 | 16 | D | True | True | seed_mutation_window,known_interface_window |
| A44 | 24 | 22 | L | True | False | known_interface_window |
| A47 | 27 | 25 | A | False | True | known_interface_window |
| A56 | 38 | 34 | D | False | True | known_interface_window |
| A57 | 39 | 35 | K | False | True | known_interface_window |
| A58 | 40 | 36 | L | True | False | known_interface_window |
| A62 | 44 | 40 | K | True | True | known_interface_window |
| A83 | 65 | 61 | Y | True | False | seed_mutation_window |
| A87 | 69 | 65 | K | True | False | seed_mutation_window |
| A91 | 73 | 69 | L | True | False | known_interface_window |
| A93 | 75 | 71 | K | False | True | known_interface_window |
| A94 | 76 | 72 | K | False | True | known_interface_window |
| A97 | 79 | 75 | L | False | True | known_interface_window |
| A98 | 80 | 76 | A | False | True | known_interface_window |
| A101 | 83 | 79 | E | True | True | known_interface_window |
| A115 | 97 | 93 | I | True | False | seed_mutation_window |
| A116 | 98 | 94 | D | True | False | seed_mutation_window |
| A118 | 100 | 96 | A | True | False | seed_mutation_window |

## hans_pocket_ss_only_u02

- Starting structure: `results/ligandmpnn_smoke/hans_pocket_ss_only_u02/packed/hans_pocket_ss_only_u02_packed_2_1.pdb`
- Seed mutation canonical positions: `18,83,93,98`
- Redesignable canonical positions: `18,24,27,38,39,40,44,73,75,76,79,80,83,93,97,98,100`

| residue_id | canonical | am1 | seed_aa | mutable_second_sphere | mutable_interface | reason |
| --- | --- | --- | --- | --- | --- | --- |
| A35 | 18 | 16 | T | True | True | seed_mutation_window,known_interface_window |
| A41 | 24 | 22 | L | True | False | known_interface_window |
| A44 | 27 | 25 | A | False | True | known_interface_window |
| A55 | 38 | 34 | T | False | True | known_interface_window |
| A56 | 39 | 35 | A | False | True | known_interface_window |
| A57 | 40 | 36 | I | True | False | known_interface_window |
| A61 | 44 | 40 | G | True | True | known_interface_window |
| A90 | 73 | 69 | L | True | False | known_interface_window |
| A92 | 75 | 71 | M | False | True | known_interface_window |
| A93 | 76 | 72 | D | False | True | known_interface_window |
| A96 | 79 | 75 | L | False | True | known_interface_window |
| A97 | 80 | 76 | K | False | True | seed_mutation_window,known_interface_window |
| A100 | 83 | 79 | K | True | True | seed_mutation_window,known_interface_window |
| A110 | 93 | 89 | G | True | False | seed_mutation_window |
| A114 | 97 | 93 | L | True | False | seed_mutation_window |
| A115 | 98 | 94 | D | True | False | seed_mutation_window |
| A117 | 100 | 96 | A | True | False | seed_mutation_window |

## hans_interface_ss_plus_if_u04

- Starting structure: `results/ligandmpnn_smoke/hans_interface_ss_plus_if_u04/packed/hans_interface_ss_plus_if_u04_packed_2_1.pdb`
- Seed mutation canonical positions: `18,38,93,98`
- Redesignable canonical positions: `18,24,27,38,39,40,44,73,75,76,79,80,83,93,97,98,100`

| residue_id | canonical | am1 | seed_aa | mutable_second_sphere | mutable_interface | reason |
| --- | --- | --- | --- | --- | --- | --- |
| A35 | 18 | 16 | S | True | True | seed_mutation_window,known_interface_window |
| A41 | 24 | 22 | L | True | False | known_interface_window |
| A44 | 27 | 25 | A | False | True | known_interface_window |
| A55 | 38 | 34 | W | False | True | seed_mutation_window,known_interface_window |
| A56 | 39 | 35 | A | False | True | seed_mutation_window,known_interface_window |
| A57 | 40 | 36 | I | True | False | seed_mutation_window,known_interface_window |
| A61 | 44 | 40 | G | True | True | known_interface_window |
| A90 | 73 | 69 | L | True | False | known_interface_window |
| A92 | 75 | 71 | M | False | True | known_interface_window |
| A93 | 76 | 72 | D | False | True | known_interface_window |
| A96 | 79 | 75 | L | False | True | known_interface_window |
| A97 | 80 | 76 | K | False | True | known_interface_window |
| A100 | 83 | 79 | R | True | True | known_interface_window |
| A110 | 93 | 89 | G | True | False | seed_mutation_window |
| A114 | 97 | 93 | L | True | False | seed_mutation_window |
| A115 | 98 | 94 | E | True | False | seed_mutation_window |
| A117 | 100 | 96 | A | True | False | seed_mutation_window |
