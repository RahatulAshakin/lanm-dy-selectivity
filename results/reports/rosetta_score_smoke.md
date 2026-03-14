# Rosetta Score Smoke

Deterministic Phase 5A1 Rosetta `score_jd2` baseline scoring for the LigandMPNN shortlist.

- rosetta_bin_dir: `/home/ashak/apps/rosetta.source.release-408/main/source/bin`
- rosetta_database: `/home/ashak/apps/rosetta.source.release-408/main/database`
- score application: `score_jd2`
- seed: `37`
- nstruct: `1`
- fa_elec weight override: `0.0`
- representative selection: `highest ligand_confidence`, then `highest overall_confidence`, then `smallest design_id`
- preserved input handling: original packed PDB scored directly with waters retained, auto metal setup enabled, headers preserved, and PDB renumbering disabled
- excluded in this phase: `relax`, `fixbb`, `MD`, `QM`, and quantum steps

| rosetta_rank | candidate_id | campaign_id | backbone_id | metal | representative_design_id | total_score |
| ---: | --- | --- | --- | --- | ---: | ---: |
| 1 | am1_mex_ss_only_u02 | am1_mex_ss_only | am1_mex_8fns_chain_a | ND | 2 | 111.315 |
| 2 | am1_mex_ss_only_u03 | am1_mex_ss_only | am1_mex_8fns_chain_a | ND | 1 | 111.951 |
| 3 | hans_pocket_ss_only_u02 | hans_pocket_ss_only | hans_pocket_8fnr_chain_a | DY | 2 | 150.783 |
| 4 | hans_pocket_ss_only_u01 | hans_pocket_ss_only | hans_pocket_8fnr_chain_a | DY | 1 | 155.946 |
| 5 | hans_interface_ss_plus_if_u04 | hans_interface_ss_plus_if | hans_interface_8fnr_a_b_c_d | DY | 2 | 531.324 |
| 6 | hans_interface_ss_plus_if_u02 | hans_interface_ss_plus_if | hans_interface_8fnr_a_b_c_d | DY | 2 | 599.480 |
| 7 | hans_interface_ss_plus_if_u01 | hans_interface_ss_plus_if | hans_interface_8fnr_a_b_c_d | DY | 2 | 603.399 |
| 8 | hans_interface_if_only_u01 | hans_interface_if_only | hans_interface_8fnr_a_b_c_d | DY | 2 | 604.090 |
