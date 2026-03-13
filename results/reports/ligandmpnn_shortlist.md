# LigandMPNN Candidate Shortlist

Deterministic Phase 4C deduplication of existing LigandMPNN smoke outputs into a Rosetta-ready shortlist.

- Per-candidate prefilter: top `2` unique designed-chain sequences by ligand_confidence then overall_confidence.
- Final shortlist limit: `8` designs.
- `hans_interface_ss_plus_if` contributes at least `2` designs when available.
- `hans_pocket_ss_only` contributes at least `1` design when available.
- `am1_mex*` campaigns contribute at least `1` design when available.

## Candidate Counts

| candidate_id | campaign_id | raw_designs | unique_designed_chain_sequences | top2_considered | retained_shortlist |
| --- | --- | ---: | ---: | ---: | ---: |
| hans_interface_ss_plus_if_u01 | hans_interface_ss_plus_if | 2 | 2 | 2 | 1 |
| hans_interface_ss_plus_if_u02 | hans_interface_ss_plus_if | 2 | 2 | 2 | 1 |
| am1_mex_ss_only_u01 | am1_mex_ss_only | 2 | 1 | 1 | 0 |
| hans_interface_if_only_u01 | hans_interface_if_only | 2 | 2 | 2 | 1 |
| hans_pocket_ss_only_u01 | hans_pocket_ss_only | 2 | 1 | 1 | 1 |
| hans_interface_ss_plus_if_u03 | hans_interface_ss_plus_if | 2 | 1 | 1 | 0 |
| am1_mex_ss_only_u02 | am1_mex_ss_only | 2 | 2 | 2 | 1 |
| hans_interface_if_only_u02 | hans_interface_if_only | 2 | 2 | 2 | 0 |
| hans_pocket_ss_only_u02 | hans_pocket_ss_only | 2 | 1 | 1 | 1 |
| hans_interface_ss_plus_if_u04 | hans_interface_ss_plus_if | 2 | 2 | 2 | 1 |
| am1_mex_ss_only_u03 | am1_mex_ss_only | 2 | 1 | 1 | 1 |
| hans_interface_if_only_u03 | hans_interface_if_only | 2 | 2 | 2 | 0 |

## Duplicate Counts

| scope | count |
| --- | ---: |
| raw generated designs | 24 |
| within-candidate unique designed-chain sequences | 19 |
| globally unique designed-chain sequences | 16 |
| within-candidate duplicates collapsed | 5 |
| cross-candidate duplicates collapsed | 3 |
| cross-candidate duplicate groups | 3 |

## Retained Shortlist

| rank | sequence_id | candidate_id | campaign_id | ligand_confidence | overall_confidence | seq_recovery | mutations_vs_candidate_input | redesigned_residue_ids | retained_because |
| ---: | --- | --- | --- | ---: | ---: | ---: | --- | --- | --- |
| 1 | hans_interface_ss_plus_if_u02_design_02 | hans_interface_ss_plus_if_u02 | hans_interface_ss_plus_if | 0.3879 | 0.4410 | 0.4286 | T12S,T16D,D21P,Y59A,K77R | A35,A39,A44,A82,A100 | required interface campaign coverage |
| 2 | hans_interface_ss_plus_if_u04_design_02 | hans_interface_ss_plus_if_u04 | hans_interface_ss_plus_if | 0.3853 | 0.4239 | 0.3333 | T12S,D21A,V32W,Y59A,D92E | A35,A44,A55,A82,A115 | required interface campaign coverage |
| 3 | hans_pocket_ss_only_u02_design_02 | hans_pocket_ss_only_u02 | hans_pocket_ss_only | 0.4251 | 0.4261 | 0.3333 | T16D,Y59A,E77K | A39,A82,A100 | required pocket campaign coverage |
| 4 | am1_mex_ss_only_u03_design_01 | am1_mex_ss_only_u03 | am1_mex_ss_only | 0.3526 | 0.3547 | 0.1667 | G10D,V87I | A38,A115 | required AM1/Mex coverage |
| 5 | hans_interface_ss_plus_if_u01_design_02 | hans_interface_ss_plus_if_u01 | hans_interface_ss_plus_if | 0.3621 | 0.4121 | 0.2857 | T12S,D21P,V32W,Y59A,K77R,D92E | A35,A44,A55,A82,A100,A115 | balanced round-robin fill (round 1) |
| 6 | hans_pocket_ss_only_u01_design_01 | hans_pocket_ss_only_u01 | hans_pocket_ss_only | 0.4244 | 0.4220 | 0.2000 | Y59A,E77K,D92E | A82,A100,A115 | balanced round-robin fill (round 1) |
| 7 | am1_mex_ss_only_u02_design_02 | am1_mex_ss_only_u02 | am1_mex_ss_only | 0.3307 | 0.3364 | 0.0000 | N10D | A38 | balanced round-robin fill (round 1) |
| 8 | hans_interface_if_only_u01_design_02 | hans_interface_if_only_u01 | hans_interface_if_only | 0.2383 | 0.2920 | 0.3333 | T12S,D21P,K77R | A35,A44,A100 | balanced round-robin fill (round 1) |
