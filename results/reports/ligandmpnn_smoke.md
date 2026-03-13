# LigandMPNN Smoke Run

Deterministic Phase 4B LigandMPNN smoke pass across the 12 shortlisted candidates from Phase 4A.

- LigandMPNN root: `/home/ashak/apps/LigandMPNN`
- model_type: `ligand_mpnn`
- seed: `37`
- batch_size: `2`
- number_of_batches: `1`
- temperature: `0.1`
- ligand atom context: `True`
- side-chain context: `True`
- pack_side_chains: `True`
- pack_with_ligand_context: `True`

| candidate_id | source_backbone | preserved_metal_identity | redesigned_residues | sequence_count | confidence_diversity_summary |
| --- | --- | --- | --- | ---: | --- |
| hans_interface_ss_plus_if_u01 | hans_interface_8fnr_a_b_c_d (8FNR) | DY | A35,A44,A55,A82,A100,A110,A115 | 2 | overall mean 0.407 (min 0.402, max 0.412); ligand mean 0.359 (min 0.357, max 0.362); 2/2 unique; pairwise identity mean 0.998 (min 0.998, max 0.998) |
| hans_interface_ss_plus_if_u02 | hans_interface_8fnr_a_b_c_d (8FNR) | DY | A35,A39,A44,A82,A100,A110,A115 | 2 | overall mean 0.430 (min 0.419, max 0.441); ligand mean 0.380 (min 0.373, max 0.388); 2/2 unique; pairwise identity mean 0.998 (min 0.998, max 0.998) |
| am1_mex_ss_only_u01 | am1_mex_8fns_chain_a (8FNS) | ND | A38,A83,A87,A111,A118 | 2 | overall mean 0.312 (min 0.305, max 0.318); ligand mean 0.289 (min 0.281, max 0.298); 1/2 unique; pairwise identity mean 1.000 (min 1.000, max 1.000) |
| hans_interface_if_only_u01 | hans_interface_8fnr_a_b_c_d (8FNR) | DY | A35,A44,A100 | 2 | overall mean 0.279 (min 0.267, max 0.292); ligand mean 0.227 (min 0.216, max 0.238); 2/2 unique; pairwise identity mean 0.995 (min 0.995, max 0.995) |
| hans_pocket_ss_only_u01 | hans_pocket_8fnr_chain_a (8FNR) | DY | A35,A82,A100,A110,A115 | 2 | overall mean 0.423 (min 0.422, max 0.424); ligand mean 0.423 (min 0.421, max 0.424); 1/2 unique; pairwise identity mean 1.000 (min 1.000, max 1.000) |
| hans_interface_ss_plus_if_u03 | hans_interface_8fnr_a_b_c_d (8FNR) | DY | A35,A44,A82,A100,A110,A115 | 2 | overall mean 0.417 (min 0.414, max 0.420); ligand mean 0.348 (min 0.346, max 0.350); 1/2 unique; pairwise identity mean 1.000 (min 1.000, max 1.000) |
| am1_mex_ss_only_u02 | am1_mex_8fns_chain_a (8FNS) | ND | A38,A83,A87,A118 | 2 | overall mean 0.335 (min 0.334, max 0.336); ligand mean 0.325 (min 0.319, max 0.331); 2/2 unique; pairwise identity mean 0.990 (min 0.990, max 0.990) |
| hans_interface_if_only_u02 | hans_interface_8fnr_a_b_c_d (8FNR) | DY | A35,A44,A55,A100 | 2 | overall mean 0.294 (min 0.290, max 0.298); ligand mean 0.234 (min 0.230, max 0.238); 2/2 unique; pairwise identity mean 0.998 (min 0.998, max 0.998) |
| hans_pocket_ss_only_u02 | hans_pocket_8fnr_chain_a (8FNR) | DY | A35,A39,A82,A100,A110,A115 | 2 | overall mean 0.421 (min 0.416, max 0.426); ligand mean 0.420 (min 0.415, max 0.425); 1/2 unique; pairwise identity mean 1.000 (min 1.000, max 1.000) |
| hans_interface_ss_plus_if_u04 | hans_interface_8fnr_a_b_c_d (8FNR) | DY | A35,A44,A55,A82,A110,A115 | 2 | overall mean 0.420 (min 0.416, max 0.424); ligand mean 0.379 (min 0.372, max 0.385); 2/2 unique; pairwise identity mean 0.998 (min 0.998, max 0.998) |
| am1_mex_ss_only_u03 | am1_mex_8fns_chain_a (8FNS) | ND | A38,A83,A87,A111,A115,A118 | 2 | overall mean 0.349 (min 0.343, max 0.355); ligand mean 0.345 (min 0.337, max 0.353); 1/2 unique; pairwise identity mean 1.000 (min 1.000, max 1.000) |
| hans_interface_if_only_u03 | hans_interface_8fnr_a_b_c_d (8FNR) | DY | A35,A44,A55,A100 | 2 | overall mean 0.294 (min 0.290, max 0.298); ligand mean 0.234 (min 0.230, max 0.238); 2/2 unique; pairwise identity mean 0.998 (min 0.998, max 0.998) |
