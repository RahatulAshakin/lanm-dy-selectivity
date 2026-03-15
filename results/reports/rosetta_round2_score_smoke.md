# Rosetta Round-2 Score Smoke

Deterministic Phase 7F Rosetta `score_jd2` baseline scoring for the 4 round-2 shortlisted candidates.

- rosetta_bin_dir: `/home/ashak/apps/rosetta.source.release-408/main/source/bin`
- rosetta_database: `/home/ashak/apps/rosetta.source.release-408/main/database`
- score application: `score_jd2`
- seed: `37`
- nstruct: `1`
- fa_elec weight override: `0.0`
- representative selection: `highest ligand_confidence`, then `highest overall_confidence`, then `smallest design_id`
- ranking scope: `within topology class only`
- cross-topology note: Raw Rosetta total_score is not a valid cross-topology comparator; Phase 7F ranks candidates only within am1_monomer, hans_monomer, and hans_interface_multichain topology classes.
- excluded in this phase: `relax`, `MD`, `QM`, and quantum steps

## Within-Topology Ranking

| topology_class | rosetta_rank_within_topology | candidate_id | representative_design_id | ligand_confidence | total_score |
| --- | ---: | --- | ---: | ---: | ---: |
| am1_monomer | 1 | am1_mex_ss_only_u02_r2u02 | 1 | 0.2213 | 218.234 |
| hans_monomer | 1 | hans_pocket_ss_only_u02_r2u01 | 1 | 1.0000 | 144.070 |
| hans_interface_multichain | 1 | hans_interface_ss_plus_if_u04_r2u01 | 2 | 0.3454 | 530.712 |
| hans_interface_multichain | 2 | hans_interface_ss_plus_if_u04_r2u02 | 1 | 0.3300 | 538.078 |

## Reduced MD Rescreen Panel

| panel_rank | panel_role | candidate_id | topology_class | rosetta_rank_within_topology | ligand_confidence | total_score |
| ---: | --- | --- | --- | ---: | ---: | ---: |
| 1 | top_am1_mex_round2_candidate | am1_mex_ss_only_u02_r2u02 | am1_monomer | 1 | 0.2213 | 218.234 |
| 2 | top_hans_pocket_round2_candidate | hans_pocket_ss_only_u02_r2u01 | hans_monomer | 1 | 1.0000 | 144.070 |
| 3 | best_hans_interface_aware_round2_candidate | hans_interface_ss_plus_if_u04_r2u01 | hans_interface_multichain | 1 | 0.3454 | 530.712 |
