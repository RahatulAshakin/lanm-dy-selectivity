# MD Input Preparation

Phase 6A1 prepares deterministic OpenMM-ready starting structures and metadata for the selected MD validation panel.

## Scope

- panel members prepared: `6`
- target metals prepared per panel member: `Dy, Nd, Y, Al, Fe`
- total systems written: `30`
- replicate_count recorded for later production MD: `3`
- target temperature recorded for later production MD: `298 K`
- suggested collective variables recorded for later metadynamics: `coordination_number, mean_metal_oxygen_distance`
- excluded in this phase: `OpenMM System building`, `MD`, `metadynamics`, `QM`, and quantum steps

## Prepared Systems

| panel_member_id | target_metal | topology_class | metal_site_count | preserved_solvent_residue_count | source_structure |
| --- | --- | --- | ---: | ---: | --- |
| am1_mex_ss_only_u02 | Dy | am1_monomer | 4 | 0 | results/ligandmpnn_smoke/am1_mex_ss_only_u02/packed/am1_mex_ss_only_u02_packed_2_1.pdb |
| am1_mex_ss_only_u02 | Nd | am1_monomer | 4 | 0 | results/ligandmpnn_smoke/am1_mex_ss_only_u02/packed/am1_mex_ss_only_u02_packed_2_1.pdb |
| am1_mex_ss_only_u02 | Y | am1_monomer | 4 | 0 | results/ligandmpnn_smoke/am1_mex_ss_only_u02/packed/am1_mex_ss_only_u02_packed_2_1.pdb |
| am1_mex_ss_only_u02 | Al | am1_monomer | 4 | 0 | results/ligandmpnn_smoke/am1_mex_ss_only_u02/packed/am1_mex_ss_only_u02_packed_2_1.pdb |
| am1_mex_ss_only_u02 | Fe | am1_monomer | 4 | 0 | results/ligandmpnn_smoke/am1_mex_ss_only_u02/packed/am1_mex_ss_only_u02_packed_2_1.pdb |
| am1_mex_wt_reference | Dy | am1_monomer | 4 | 25 | results/design_inputs/proteinmpnn/backbones/am1_mex_8fns_chain_a.pdb |
| am1_mex_wt_reference | Nd | am1_monomer | 4 | 25 | results/design_inputs/proteinmpnn/backbones/am1_mex_8fns_chain_a.pdb |
| am1_mex_wt_reference | Y | am1_monomer | 4 | 25 | results/design_inputs/proteinmpnn/backbones/am1_mex_8fns_chain_a.pdb |
| am1_mex_wt_reference | Al | am1_monomer | 4 | 25 | results/design_inputs/proteinmpnn/backbones/am1_mex_8fns_chain_a.pdb |
| am1_mex_wt_reference | Fe | am1_monomer | 4 | 25 | results/design_inputs/proteinmpnn/backbones/am1_mex_8fns_chain_a.pdb |
| hans_pocket_ss_only_u02 | Dy | hans_monomer | 4 | 0 | results/ligandmpnn_smoke/hans_pocket_ss_only_u02/packed/hans_pocket_ss_only_u02_packed_2_1.pdb |
| hans_pocket_ss_only_u02 | Nd | hans_monomer | 4 | 0 | results/ligandmpnn_smoke/hans_pocket_ss_only_u02/packed/hans_pocket_ss_only_u02_packed_2_1.pdb |
| hans_pocket_ss_only_u02 | Y | hans_monomer | 4 | 0 | results/ligandmpnn_smoke/hans_pocket_ss_only_u02/packed/hans_pocket_ss_only_u02_packed_2_1.pdb |
| hans_pocket_ss_only_u02 | Al | hans_monomer | 4 | 0 | results/ligandmpnn_smoke/hans_pocket_ss_only_u02/packed/hans_pocket_ss_only_u02_packed_2_1.pdb |
| hans_pocket_ss_only_u02 | Fe | hans_monomer | 4 | 0 | results/ligandmpnn_smoke/hans_pocket_ss_only_u02/packed/hans_pocket_ss_only_u02_packed_2_1.pdb |
| hans_interface_ss_plus_if_u04 | Dy | hans_interface_multichain | 14 | 0 | results/ligandmpnn_smoke/hans_interface_ss_plus_if_u04/packed/hans_interface_ss_plus_if_u04_packed_2_1.pdb |
| hans_interface_ss_plus_if_u04 | Nd | hans_interface_multichain | 14 | 0 | results/ligandmpnn_smoke/hans_interface_ss_plus_if_u04/packed/hans_interface_ss_plus_if_u04_packed_2_1.pdb |
| hans_interface_ss_plus_if_u04 | Y | hans_interface_multichain | 14 | 0 | results/ligandmpnn_smoke/hans_interface_ss_plus_if_u04/packed/hans_interface_ss_plus_if_u04_packed_2_1.pdb |
| hans_interface_ss_plus_if_u04 | Al | hans_interface_multichain | 14 | 0 | results/ligandmpnn_smoke/hans_interface_ss_plus_if_u04/packed/hans_interface_ss_plus_if_u04_packed_2_1.pdb |
| hans_interface_ss_plus_if_u04 | Fe | hans_interface_multichain | 14 | 0 | results/ligandmpnn_smoke/hans_interface_ss_plus_if_u04/packed/hans_interface_ss_plus_if_u04_packed_2_1.pdb |
| hans_interface_wt_reference | Dy | hans_interface_multichain | 14 | 40 | results/design_inputs/proteinmpnn/backbones/hans_interface_8fnr_a_b_c_d.pdb |
| hans_interface_wt_reference | Nd | hans_interface_multichain | 14 | 40 | results/design_inputs/proteinmpnn/backbones/hans_interface_8fnr_a_b_c_d.pdb |
| hans_interface_wt_reference | Y | hans_interface_multichain | 14 | 40 | results/design_inputs/proteinmpnn/backbones/hans_interface_8fnr_a_b_c_d.pdb |
| hans_interface_wt_reference | Al | hans_interface_multichain | 14 | 40 | results/design_inputs/proteinmpnn/backbones/hans_interface_8fnr_a_b_c_d.pdb |
| hans_interface_wt_reference | Fe | hans_interface_multichain | 14 | 40 | results/design_inputs/proteinmpnn/backbones/hans_interface_8fnr_a_b_c_d.pdb |
| hans_interface_ss_plus_if_u02 | Dy | hans_interface_multichain | 14 | 0 | results/ligandmpnn_smoke/hans_interface_ss_plus_if_u02/packed/hans_interface_ss_plus_if_u02_packed_2_1.pdb |
| hans_interface_ss_plus_if_u02 | Nd | hans_interface_multichain | 14 | 0 | results/ligandmpnn_smoke/hans_interface_ss_plus_if_u02/packed/hans_interface_ss_plus_if_u02_packed_2_1.pdb |
| hans_interface_ss_plus_if_u02 | Y | hans_interface_multichain | 14 | 0 | results/ligandmpnn_smoke/hans_interface_ss_plus_if_u02/packed/hans_interface_ss_plus_if_u02_packed_2_1.pdb |
| hans_interface_ss_plus_if_u02 | Al | hans_interface_multichain | 14 | 0 | results/ligandmpnn_smoke/hans_interface_ss_plus_if_u02/packed/hans_interface_ss_plus_if_u02_packed_2_1.pdb |
| hans_interface_ss_plus_if_u02 | Fe | hans_interface_multichain | 14 | 0 | results/ligandmpnn_smoke/hans_interface_ss_plus_if_u02/packed/hans_interface_ss_plus_if_u02_packed_2_1.pdb |

## Wild-Type Reference Handling

- Wild-type reference backbones remain the panel entry points, but missing metal-site context is restored deterministically from the original backbone source structure when needed.
- `am1_mex_wt_reference` uses `data/raw/local_bundle/8fns_atoms.csv` as the metal/water template while preserving the protein backbone from `results/design_inputs/proteinmpnn/backbones/am1_mex_8fns_chain_a.pdb`.
- `hans_interface_wt_reference` uses `data/raw/public/structures/8FNR.cif` as the metal/water template while preserving the protein backbone from `results/design_inputs/proteinmpnn/backbones/hans_interface_8fnr_a_b_c_d.pdb`.
