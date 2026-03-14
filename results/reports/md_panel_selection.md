# MD Panel Selection

Phase 5B integrates ProteinMPNN prioritization, LigandMPNN confidence, and Rosetta score_jd2 screening to nominate a deterministic MD validation panel.

## Scope

- target metal: `Dy`
- competitors for later MD validation: `Nd, Y, Al, Fe`
- temperature: `298 K`
- excluded in this phase: `MD`, `QM`, and quantum steps

## Why raw Rosetta total_score is not used across topologies

- Raw Rosetta total_score is not a valid cross-topology comparator; only within-backbone/topology rank and normalized score are used across AM1 monomer, Hans monomer, and Hans multichain assemblies.
- `am1_monomer` LigandMPNN sequence lengths in the current set: `105` residues.
- `hans_interface_multichain` LigandMPNN sequence lengths in the current set: `435` residues.
- `hans_monomer` LigandMPNN sequence lengths in the current set: `110` residues.
- Phase 5B therefore compares Rosetta scores only within the same backbone/topology cohort and uses within-backbone rank plus normalized score in the integrated ranking.

## Integrated ranking

| integrated_rank | candidate_id | topology_class | proteinmpnn_rank | ligand_confidence | rosetta_rank_within_topology | rosetta_total_score |
| ---: | --- | --- | ---: | ---: | ---: | ---: |
| 1 | hans_pocket_ss_only_u02 | hans_monomer | 9 | 0.4251 | 1 | 150.783 |
| 2 | am1_mex_ss_only_u02 | am1_monomer | 7 | 0.3307 | 1 | 111.315 |
| 3 | hans_interface_ss_plus_if_u04 | hans_interface_multichain | 10 | 0.3853 | 1 | 531.324 |
| 4 | hans_interface_ss_plus_if_u02 | hans_interface_multichain | 2 | 0.3879 | 2 | 599.480 |
| 5 | hans_interface_ss_plus_if_u01 | hans_interface_multichain | 1 | 0.3621 | 3 | 603.399 |
| 6 | hans_pocket_ss_only_u01 | hans_monomer | 5 | 0.4244 | 2 | 155.946 |
| 7 | hans_interface_if_only_u01 | hans_interface_multichain | 4 | 0.2383 | 4 | 604.090 |
| 8 | am1_mex_ss_only_u03 | am1_monomer | 11 | 0.3526 | 2 | 111.951 |

## MD panel

The selected panel keeps scaffold coverage, includes explicit wild-type baselines, and preserves candidates that are best-supported by LigandMPNN confidence plus within-topology Rosetta standing.

- `am1_mex_ss_only_u02` (am1_mex_designed_candidate): Retained as the best AM1/Mex design by integrated Phase 5B ranking; LigandMPNN ligand_confidence=0.3307 and Rosetta within-topology rank=1.
- `am1_mex_wt_reference` (am1_mex_wild_type_reference): Kept as the AM1/Mex wild-type baseline so later MD can measure whether redesigned AM1 variants improve Dy selectivity relative to the unmodified monomer scaffold.
- `hans_pocket_ss_only_u02` (hans_pocket_designed_candidate): Retained as the best Hans pocket-focused design; LigandMPNN ligand_confidence=0.4251 and Rosetta within-topology rank=1.
- `hans_interface_ss_plus_if_u04` (hans_interface_aware_candidate): Retained as the best Hans interface-aware design to preserve the multichain hypothesis; LigandMPNN ligand_confidence=0.3853 and Rosetta within-topology rank=1.
- `hans_interface_wt_reference` (hans_wild_type_reference): Kept as the Hans wild-type baseline so later MD can separate mutation-driven effects from the native Hans scaffold and interface context.
- `hans_interface_ss_plus_if_u02` (best_remaining_designed_candidate): Retained as the strongest remaining designed candidate after mandatory scaffold coverage; LigandMPNN ligand_confidence=0.3879 and Rosetta within-topology rank=2.

## Dy-vs-Nd/Y/Al/Fe objective support

- The panel spans both AM1/Mex and Hans scaffolds, so later MD can test whether `Dy` preference emerges from scaffold identity or from specific redesigned residues.
- The Hans pocket-focused and interface-aware designs decouple local pocket tuning from multichain/interface effects, which is important for later Dy-versus-Nd/Y discrimination and for checking whether interface stabilization helps or hurts selectivity.
- The two wild-type references provide baseline behavior for off-target Al/Fe coordination and for Dy-versus-Nd/Y comparison before mutations are introduced.
- Because this phase keeps one best remaining designed candidate after satisfying coverage constraints, the panel still allocates one slot to the strongest extra design signal rather than stopping at categorical minimums.
