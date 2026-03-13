# ProteinMPNN Candidate Shortlist

Deterministic Phase 3C deduplication of Phase 3B ProteinMPNN smoke outputs.

- Per-campaign prefilter: top `5` unique sequences by best score.
- Global shortlist limit: `12` sequences.
- Priority campaign `hans_interface_ss_plus_if` contributes at least `2` sequences when available.
- Exact designed-sequence duplicates across campaigns were excluded from the final shortlist: `0` skipped candidates.

## Campaign Counts

| campaign_id | raw_sequences | unique_sequences | top5_considered | backbone_id | design_set_name |
| --- | ---: | ---: | ---: | --- | --- |
| am1_mex_ss_only | 40 | 9 | 5 | am1_mex_8fns_chain_a | campaign_ss_only |
| hans_interface_if_only | 40 | 10 | 5 | hans_interface_8fnr_a_b_c_d | campaign_if_only |
| hans_interface_ss_plus_if | 40 | 16 | 5 | hans_interface_8fnr_a_b_c_d | campaign_ss_plus_if |
| hans_pocket_ss_only | 40 | 9 | 5 | hans_pocket_8fnr_chain_a | campaign_ss_only |

## Mutation Count Distribution

| mutation_count | unique_sequences |
| --- | ---: |
| 2 | 1 |
| 3 | 5 |
| 4 | 6 |
| 5 | 12 |
| 6 | 11 |
| 7 | 8 |
| 8 | 1 |

## Final Shortlist

| rank | candidate_id | campaign_id | best_score | best_global_score | best_seq_recovery | mutations | retained_because |
| ---: | --- | --- | ---: | ---: | ---: | --- | --- |
| 1 | hans_interface_ss_plus_if_u01 | hans_interface_ss_plus_if | 0.5185 | 1.3563 | 0.6957 | K12T,A21D,T32V,A59Y,R77K,K87G,T92D | required interface campaign coverage |
| 2 | hans_interface_ss_plus_if_u02 | hans_interface_ss_plus_if | 0.5251 | 1.3791 | 0.6957 | K12T,D16T,A21D,A59Y,R77K,K87G,T92D | required interface campaign coverage |
| 3 | am1_mex_ss_only_u01 | am1_mex_ss_only | 0.5021 | 0.9404 | 0.6875 | K10G,L55Y,N59K,N83G,R90A | required campaign coverage |
| 4 | hans_interface_if_only_u01 | hans_interface_if_only | 0.6158 | 1.3814 | 0.7000 | K12T,A21D,R77K | required campaign coverage |
| 5 | hans_pocket_ss_only_u01 | hans_pocket_ss_only | 0.5847 | 1.4511 | 0.6875 | K12T,A59Y,R77E,K87G,T92D | required campaign coverage |
| 6 | hans_interface_ss_plus_if_u03 | hans_interface_ss_plus_if | 0.5423 | 1.3616 | 0.7391 | K12T,A21D,A59Y,R77K,K87G,T92D | balanced round-robin fill (round 1) |
| 7 | am1_mex_ss_only_u02 | am1_mex_ss_only | 0.5134 | 0.9620 | 0.7500 | K10N,L55Y,N59K,R90A | balanced round-robin fill (round 1) |
| 8 | hans_interface_if_only_u02 | hans_interface_if_only | 0.6293 | 1.3662 | 0.6000 | K12T,A21D,T32V,R77K | balanced round-robin fill (round 1) |
| 9 | hans_pocket_ss_only_u02 | hans_pocket_ss_only | 0.6053 | 1.4756 | 0.6250 | K12T,D16T,A59Y,R77E,K87G,T92D | balanced round-robin fill (round 1) |
| 10 | hans_interface_ss_plus_if_u04 | hans_interface_ss_plus_if | 0.5450 | 1.3557 | 0.7391 | K12T,A21D,T32V,A59Y,K87G,T92D | balanced round-robin fill (round 2) |
| 11 | am1_mex_ss_only_u03 | am1_mex_ss_only | 0.5408 | 0.9529 | 0.6250 | K10G,L55Y,N59K,N83G,I87V,R90A | balanced round-robin fill (round 2) |
| 12 | hans_interface_if_only_u03 | hans_interface_if_only | 0.6306 | 1.3681 | 0.6000 | K12T,A21D,T32R,R77K | balanced round-robin fill (round 2) |
