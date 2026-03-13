# ProteinMPNN Smoke Run

Deterministic Phase 3B ProteinMPNN smoke run across the exported Phase 3A campaigns.

- ProteinMPNN root: `/home/ashak/apps/ProteinMPNN`
- num_seq_per_target: `20`
- sampling_temp: `0.1 0.15`
- seed: `37`
- batch_size: `1`

| campaign_id | sequences_generated | designed_chains | backbone_used | diversity_summary |
| --- | ---: | --- | --- | --- |
| am1_mex_ss_only | 40 | A | am1_mex_8fns_chain_a (8FNS) | 9/40 unique; pairwise identity mean 0.989 (min 0.962, max 1.000) |
| hans_interface_if_only | 40 | A | hans_interface_8fnr_a_b_c_d (8FNR) | 10/40 unique; pairwise identity mean 0.988 (min 0.964, max 1.000) |
| hans_interface_ss_plus_if | 40 | A | hans_interface_8fnr_a_b_c_d (8FNR) | 16/40 unique; pairwise identity mean 0.985 (min 0.945, max 1.000) |
| hans_pocket_ss_only | 40 | A | hans_pocket_8fnr_chain_a (8FNR) | 9/40 unique; pairwise identity mean 0.994 (min 0.982, max 1.000) |
