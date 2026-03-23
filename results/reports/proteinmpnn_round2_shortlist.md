# ProteinMPNN Round-2 Shortlist

Phase 7C deterministic deduplication and shortlist selection from Phase 7B round-2 ProteinMPNN smoke outputs.

- Per-seed prefilter: top `2` unique sequences by best_score then best_global_score.
- Global shortlist limit: `6` candidates.
- At least one candidate per seed scaffold is retained when a non-duplicate sequence is available.
- Exact designed-sequence duplicates across seed scaffolds were excluded from the final shortlist: `0` skipped candidates.

## Seed Counts

| seed_candidate_id | campaign_id | raw_sequences | unique_sequences | top2_considered | backbone_id |
| --- | --- | ---: | ---: | ---: | --- |
| am1_mex_ss_only_u02 | am1_mex_ss_only | 40 | 29 | 2 | am1_mex_8fns_chain_a |
| hans_pocket_ss_only_u02 | hans_pocket_ss_only | 40 | 33 | 2 | hans_pocket_8fnr_chain_a |
| hans_interface_ss_plus_if_u04 | hans_interface_ss_plus_if | 40 | 12 | 2 | hans_interface_8fnr_a_b_c_d |

## Final Shortlist

| rank | candidate_id | seed_scaffold | campaign_rank | best_score | best_global_score | best_seq_recovery | redesigned_residue_ids | retained_because |
| ---: | --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| 1 | am1_mex_ss_only_u02_r2u01 | am1_mex_ss_only_u02 | 1 | 0.7243 | 0.9420 | 0.8333 | A38,A56,A93 | required seed scaffold coverage |
| 2 | hans_pocket_ss_only_u02_r2u01 | hans_pocket_ss_only_u02 | 1 | 0.8822 | 1.4167 | 0.7647 | A44,A55,A92,A97 | required seed scaffold coverage |
| 3 | hans_interface_ss_plus_if_u04_r2u01 | hans_interface_ss_plus_if_u04 | 1 | 0.5580 | 1.3632 | 0.7647 | A35,A44,A55,A115 | required seed scaffold coverage |
| 4 | am1_mex_ss_only_u02_r2u02 | am1_mex_ss_only_u02 | 2 | 0.7261 | 0.9915 | 0.6667 | A38,A56,A57,A83,A93,A94 | balanced round-robin fill (round 2) |
| 5 | hans_pocket_ss_only_u02_r2u02 | hans_pocket_ss_only_u02 | 2 | 0.8902 | 1.4679 | 0.7647 | A44,A55,A92,A97 | balanced round-robin fill (round 2) |
| 6 | hans_interface_ss_plus_if_u04_r2u02 | hans_interface_ss_plus_if_u04 | 2 | 0.5630 | 1.3628 | 0.7059 | A35,A44,A55,A100,A115 | balanced round-robin fill (round 2) |
