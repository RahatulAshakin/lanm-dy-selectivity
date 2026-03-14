# Rosetta Relax Smoke

Deterministic Phase 5A2 restrained Rosetta `relax` smoke refinement for the top Phase 5A1 score-ranked shortlist candidates.

- rosetta_bin_dir: `/home/ashak/apps/rosetta.source.release-408/main/source/bin`
- rosetta_database: `/home/ashak/apps/rosetta.source.release-408/main/database`
- relax application: `relax`
- selected candidates: `am1_mex_ss_only_u02, am1_mex_ss_only_u03, hans_pocket_ss_only_u02`
- selected top-N: `3`
- seed: `37`
- nstruct: `1`
- fa_elec weight override: `0.0`
- restraint mode: `-relax:constrain_relax_to_start_coords` with constraint ramping disabled
- single-threading: `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `NUMEXPR_NUM_THREADS=1`, `VECLIB_MAXIMUM_THREADS=1`
- preserved input handling: representative packed PDBs reused from Phase 5A1 with waters retained, auto metal setup enabled, headers preserved, and PDB renumbering disabled
- excluded in this phase: `fixbb`, `MD`, `QM`, and quantum steps

| rosetta_relax_rank | candidate_id | campaign_id | backbone_id | metal | pre_relax_total_score | post_relax_total_score | delta_total_score |
| ---: | --- | --- | --- | --- | ---: | ---: | ---: |
| 1 | am1_mex_ss_only_u02 | am1_mex_ss_only | am1_mex_8fns_chain_a | ND | 111.315 | -114.657 | -225.972 |
| 2 | am1_mex_ss_only_u03 | am1_mex_ss_only | am1_mex_8fns_chain_a | ND | 111.951 | -113.531 | -225.482 |
| 3 | hans_pocket_ss_only_u02 | hans_pocket_ss_only | hans_pocket_8fnr_chain_a | DY | 150.783 | -100.899 | -251.682 |
