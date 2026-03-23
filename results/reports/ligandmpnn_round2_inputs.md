# LigandMPNN Round-2 Inputs

Phase 7C deterministic LigandMPNN-ready input preparation from the shortlisted round-2 ProteinMPNN candidates.

- Source PDBs were copied directly from the Phase 7A round-2 campaign exports.
- No LigandMPNN, Rosetta, MD, QM, or quantum execution steps were run in this phase.

| rank | candidate_id | seed_scaffold | source_pdb | redesigned_residue_ids | mutation_string | exported_pdb |
| ---: | --- | --- | --- | --- | --- | --- |
| 1 | am1_mex_ss_only_u02_r2u01 | am1_mex_ss_only_u02 | results/design_inputs/proteinmpnn_round2/campaigns/am1_mex_ss_only_u02/am1_mex_ss_only_u02.pdb | A38,A56,A93 | A38:D>N,A56:D>K,A93:K>E | results/design_inputs/ligandmpnn_round2/am1_mex_ss_only_u02_r2u01/am1_mex_ss_only_u02_r2u01.pdb |
| 2 | hans_pocket_ss_only_u02_r2u01 | hans_pocket_ss_only_u02 | results/design_inputs/proteinmpnn_round2/campaigns/hans_pocket_ss_only_u02/hans_pocket_ss_only_u02.pdb | A44,A55,A92,A97 | A44:A>D,A55:T>K,A92:M>L,A97:K>A | results/design_inputs/ligandmpnn_round2/hans_pocket_ss_only_u02_r2u01/hans_pocket_ss_only_u02_r2u01.pdb |
| 3 | hans_interface_ss_plus_if_u04_r2u01 | hans_interface_ss_plus_if_u04 | results/design_inputs/proteinmpnn_round2/campaigns/hans_interface_ss_plus_if_u04/hans_interface_ss_plus_if_u04.pdb | A35,A44,A55,A115 | A35:S>T,A44:A>D,A55:W>V,A115:E>D | results/design_inputs/ligandmpnn_round2/hans_interface_ss_plus_if_u04_r2u01/hans_interface_ss_plus_if_u04_r2u01.pdb |
| 4 | am1_mex_ss_only_u02_r2u02 | am1_mex_ss_only_u02 | results/design_inputs/proteinmpnn_round2/campaigns/am1_mex_ss_only_u02/am1_mex_ss_only_u02.pdb | A38,A56,A57,A83,A93,A94 | A38:D>N,A56:D>K,A57:K>A,A83:Y>W,A93:K>E,A94:K>A | results/design_inputs/ligandmpnn_round2/am1_mex_ss_only_u02_r2u02/am1_mex_ss_only_u02_r2u02.pdb |
| 5 | hans_pocket_ss_only_u02_r2u02 | hans_pocket_ss_only_u02 | results/design_inputs/proteinmpnn_round2/campaigns/hans_pocket_ss_only_u02/hans_pocket_ss_only_u02.pdb | A44,A55,A92,A97 | A44:A>D,A55:T>W,A92:M>L,A97:K>A | results/design_inputs/ligandmpnn_round2/hans_pocket_ss_only_u02_r2u02/hans_pocket_ss_only_u02_r2u02.pdb |
| 6 | hans_interface_ss_plus_if_u04_r2u02 | hans_interface_ss_plus_if_u04 | results/design_inputs/proteinmpnn_round2/campaigns/hans_interface_ss_plus_if_u04/hans_interface_ss_plus_if_u04.pdb | A35,A44,A55,A100,A115 | A35:S>T,A44:A>D,A55:W>T,A100:R>K,A115:E>D | results/design_inputs/ligandmpnn_round2/hans_interface_ss_plus_if_u04_r2u02/hans_interface_ss_plus_if_u04_r2u02.pdb |
