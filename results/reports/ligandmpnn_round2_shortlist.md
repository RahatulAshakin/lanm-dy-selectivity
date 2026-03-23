# LigandMPNN Round-2 Shortlist

Phase 7E deterministic deduplication of Phase 7D LigandMPNN smoke outputs into a reduced Rosetta score-only round-2 shortlist.

- Final shortlist limit: `4` unique designs.
- Top `am1_mex*` round-2 candidate retained when available.
- Top `hans_pocket_ss_only` round-2 candidate retained when available.
- Top `hans_interface_ss_plus_if` round-2 candidate retained when available.
- One additional best remaining unique candidate retained across all remaining candidates when available.
- Ranking uses mean ligand confidence, then mean overall confidence, then smallest round-2 shortlist rank.
- No Rosetta, MD, QM, or quantum steps were run in this phase.

## Seed Counts

| seed_candidate_id | campaign_id | topology_class | raw_designs | unique_designed_chain_sequences | duplicates_collapsed | retained_shortlist |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| am1_mex_ss_only_u02 | am1_mex_ss_only | am1_monomer | 4 | 3 | 1 | 1 |
| hans_pocket_ss_only_u02 | hans_pocket_ss_only | hans_monomer | 4 | 1 | 3 | 1 |
| hans_interface_ss_plus_if_u04 | hans_interface_ss_plus_if | hans_interface_multichain | 4 | 2 | 2 | 2 |

## Duplicate Counts

| scope | count |
| --- | ---: |
| raw generated designs | 12 |
| within-candidate unique designed-chain sequences | 8 |
| globally unique designed-chain sequences | 6 |
| within-candidate duplicates collapsed | 4 |
| cross-candidate duplicates collapsed | 2 |

## Retained Shortlist

| rank | sequence_id | candidate_id | seed_scaffold | campaign_id | mean_ligand_confidence | mean_overall_confidence | best_seq_recovery | mutations_vs_round2_seed | redesigned_residue_ids | retained_because |
| ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- | --- |
| 1 | hans_pocket_ss_only_u02_r2u01_design_01 | hans_pocket_ss_only_u02_r2u01 | hans_pocket_ss_only_u02 | hans_pocket_ss_only | 1.0000 | 0.3380 | 0.0000 | A44:A>D,A55:T>W,A92:M>L,A97:K>A | A44,A55,A92,A97 | top Hans pocket round-2 candidate |
| 2 | hans_interface_ss_plus_if_u04_r2u01_design_02 | hans_interface_ss_plus_if_u04_r2u01 | hans_interface_ss_plus_if_u04 | hans_interface_ss_plus_if | 0.3400 | 0.3331 | 1.0000 | native | native | top Hans interface-aware round-2 candidate |
| 3 | hans_interface_ss_plus_if_u04_r2u02_design_01 | hans_interface_ss_plus_if_u04_r2u02 | hans_interface_ss_plus_if_u04 | hans_interface_ss_plus_if | 0.3300 | 0.3361 | 0.8000 | A100:R>K | A100 | best remaining unique candidate |
| 4 | am1_mex_ss_only_u02_r2u02_design_01 | am1_mex_ss_only_u02_r2u02 | am1_mex_ss_only_u02 | am1_mex_ss_only | 0.2097 | 0.2994 | 0.5000 | A57:K>A,A93:K>R,A94:K>A | A57,A93,A94 | top AM1/Mex round-2 candidate |
