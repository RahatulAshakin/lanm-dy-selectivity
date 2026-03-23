# LigandMPNN Round-2 Smoke

Deterministic Phase 7D LigandMPNN smoke pass across the 6 shortlisted round-2 candidates from Phase 7C.

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
- No Rosetta, MD, QM, or quantum steps were run in this phase.

| candidate_id | seed_scaffold | preserved_metal_identity | redesigned_residues | sequence_count | overall_confidence | ligand_confidence | unique_sequence_count |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| am1_mex_ss_only_u02_r2u01 | am1_mex_ss_only_u02 | ND | A38,A56,A93 | 2 | 0.2380 | 0.2009 | 2 |
| hans_pocket_ss_only_u02_r2u01 | hans_pocket_ss_only_u02 | DY | A44,A55,A92,A97 | 2 | 0.3380 | 1.0000 | 1 |
| hans_interface_ss_plus_if_u04_r2u01 | hans_interface_ss_plus_if_u04 | DY | A35,A44,A55,A115 | 2 | 0.3311 | 0.3453 | 1 |
| am1_mex_ss_only_u02_r2u02 | am1_mex_ss_only_u02 | ND | A38,A56,A57,A83,A93,A94 | 2 | 0.2994 | 0.2097 | 1 |
| hans_pocket_ss_only_u02_r2u02 | hans_pocket_ss_only_u02 | DY | A44,A55,A92,A97 | 2 | 0.3380 | 1.0000 | 1 |
| hans_interface_ss_plus_if_u04_r2u02 | hans_interface_ss_plus_if_u04 | DY | A35,A44,A55,A100,A115 | 2 | 0.3366 | 0.3297 | 2 |

## am1_mex_ss_only_u02_r2u01

- Seed scaffold: `am1_mex_ss_only_u02`
- Backbone/topology: `am1_mex_8fns_chain_a` / `am1_monomer`
- Preserved metal identity: `ND`
- Preserved solvent residue count: `0`
- Redesigned residues: `A38,A56,A93`
- Sequence count: `2` with `2` unique sequences
- Overall confidence: `mean 0.2380; min 0.2352; max 0.2408`
- Ligand confidence: `mean 0.2009; min 0.1989; max 0.2028`

## hans_pocket_ss_only_u02_r2u01

- Seed scaffold: `hans_pocket_ss_only_u02`
- Backbone/topology: `hans_pocket_8fnr_chain_a` / `hans_monomer`
- Preserved metal identity: `DY`
- Preserved solvent residue count: `0`
- Redesigned residues: `A44,A55,A92,A97`
- Sequence count: `2` with `1` unique sequences
- Overall confidence: `mean 0.3380; min 0.3361; max 0.3399`
- Ligand confidence: `mean 1.0000; min 1.0000; max 1.0000`

## hans_interface_ss_plus_if_u04_r2u01

- Seed scaffold: `hans_interface_ss_plus_if_u04`
- Backbone/topology: `hans_interface_8fnr_a_b_c_d` / `hans_interface_multichain`
- Preserved metal identity: `DY`
- Preserved solvent residue count: `0`
- Redesigned residues: `A35,A44,A55,A115`
- Sequence count: `2` with `1` unique sequences
- Overall confidence: `mean 0.3311; min 0.3303; max 0.3320`
- Ligand confidence: `mean 0.3453; min 0.3451; max 0.3454`

## am1_mex_ss_only_u02_r2u02

- Seed scaffold: `am1_mex_ss_only_u02`
- Backbone/topology: `am1_mex_8fns_chain_a` / `am1_monomer`
- Preserved metal identity: `ND`
- Preserved solvent residue count: `0`
- Redesigned residues: `A38,A56,A57,A83,A93,A94`
- Sequence count: `2` with `1` unique sequences
- Overall confidence: `mean 0.2994; min 0.2964; max 0.3023`
- Ligand confidence: `mean 0.2097; min 0.1980; max 0.2213`

## hans_pocket_ss_only_u02_r2u02

- Seed scaffold: `hans_pocket_ss_only_u02`
- Backbone/topology: `hans_pocket_8fnr_chain_a` / `hans_monomer`
- Preserved metal identity: `DY`
- Preserved solvent residue count: `0`
- Redesigned residues: `A44,A55,A92,A97`
- Sequence count: `2` with `1` unique sequences
- Overall confidence: `mean 0.3380; min 0.3361; max 0.3399`
- Ligand confidence: `mean 1.0000; min 1.0000; max 1.0000`

## hans_interface_ss_plus_if_u04_r2u02

- Seed scaffold: `hans_interface_ss_plus_if_u04`
- Backbone/topology: `hans_interface_8fnr_a_b_c_d` / `hans_interface_multichain`
- Preserved metal identity: `DY`
- Preserved solvent residue count: `0`
- Redesigned residues: `A35,A44,A55,A100,A115`
- Sequence count: `2` with `2` unique sequences
- Overall confidence: `mean 0.3366; min 0.3361; max 0.3370`
- Ligand confidence: `mean 0.3297; min 0.3294; max 0.3300`
