# MD System Building Plan

Phase 6A2 creates deterministic OpenMM build requests and a metal-model registry for the selected MD validation panel.

## Build Defaults

- protein force field scaffold: `amber19-all.xml`
- water model scaffold: `amber19/opc3.xml`
- box padding target: `1.20 nm`
- ionic strength target: `0.15 M`
- neutralization ion policy: `monovalent_background_ions_only_excluding_panel_metals`

## Why These Metals Stay Custom

- Dy, Nd, Y, Al, and Fe are treated as custom bound species because they are carried in the starting structures as pre-bound LanM-site cofactors, not as freely exchanged bulk-solvent ions.
- The project brief explicitly calls for a custom 12-6-4 / chelator-tuned treatment of the metal center, and the repo does not yet contain numeric bound-site parameters for any panel metal.
- Phase 6A2 therefore keeps generic protein and water setup separate from the missing metal-site parameterization step needed for full OpenMM `System` creation.

## Registry Summary

| metal | formal_charge | intended_model_family | standard_forcefield_supported | custom_required | parameter_source_status |
| --- | ---: | --- | --- | --- | --- |
| Dy | +3 | custom_bound_site_12_6_4_lj_chelator_tuned | False | True | custom_parameters_not_present_in_repo |
| Nd | +3 | custom_bound_site_12_6_4_lj_chelator_tuned | False | True | custom_parameters_not_present_in_repo |
| Y | +3 | custom_bound_site_12_6_4_lj_chelator_tuned | False | True | custom_parameters_not_present_in_repo |
| Al | +3 | custom_bound_site_12_6_4_lj_chelator_tuned | False | True | custom_parameters_not_present_in_repo |
| Fe | +3 | custom_bound_site_12_6_4_lj_chelator_tuned | False | True | custom_parameters_not_present_in_repo |

## System Readiness

- systems planned: `30`
- systems that can proceed directly with generic protein/water preparation as written: `10`
- systems with optional template-water restoration flagged before solvation: `20`
- systems still blocked from full OpenMM `System` creation pending custom metal parameters: `30`

| panel_member_id | target_metal | topology_class | generic_protein_water_preparation_ready | restore_template_waters | custom_parameters_still_required | template_water_reference |
| --- | --- | --- | --- | --- | --- | --- |
| am1_mex_ss_only_u02 | Al | am1_monomer | True | True | True | am1_mex_wt_reference |
| am1_mex_ss_only_u02 | Dy | am1_monomer | True | True | True | am1_mex_wt_reference |
| am1_mex_ss_only_u02 | Fe | am1_monomer | True | True | True | am1_mex_wt_reference |
| am1_mex_ss_only_u02 | Nd | am1_monomer | True | True | True | am1_mex_wt_reference |
| am1_mex_ss_only_u02 | Y | am1_monomer | True | True | True | am1_mex_wt_reference |
| am1_mex_wt_reference | Al | am1_monomer | True | False | True | am1_mex_wt_reference |
| am1_mex_wt_reference | Dy | am1_monomer | True | False | True | am1_mex_wt_reference |
| am1_mex_wt_reference | Fe | am1_monomer | True | False | True | am1_mex_wt_reference |
| am1_mex_wt_reference | Nd | am1_monomer | True | False | True | am1_mex_wt_reference |
| am1_mex_wt_reference | Y | am1_monomer | True | False | True | am1_mex_wt_reference |
| hans_pocket_ss_only_u02 | Al | hans_monomer | True | True | True | hans_interface_wt_reference |
| hans_pocket_ss_only_u02 | Dy | hans_monomer | True | True | True | hans_interface_wt_reference |
| hans_pocket_ss_only_u02 | Fe | hans_monomer | True | True | True | hans_interface_wt_reference |
| hans_pocket_ss_only_u02 | Nd | hans_monomer | True | True | True | hans_interface_wt_reference |
| hans_pocket_ss_only_u02 | Y | hans_monomer | True | True | True | hans_interface_wt_reference |
| hans_interface_ss_plus_if_u04 | Al | hans_interface_multichain | True | True | True | hans_interface_wt_reference |
| hans_interface_ss_plus_if_u04 | Dy | hans_interface_multichain | True | True | True | hans_interface_wt_reference |
| hans_interface_ss_plus_if_u04 | Fe | hans_interface_multichain | True | True | True | hans_interface_wt_reference |
| hans_interface_ss_plus_if_u04 | Nd | hans_interface_multichain | True | True | True | hans_interface_wt_reference |
| hans_interface_ss_plus_if_u04 | Y | hans_interface_multichain | True | True | True | hans_interface_wt_reference |
| hans_interface_wt_reference | Al | hans_interface_multichain | True | False | True | hans_interface_wt_reference |
| hans_interface_wt_reference | Dy | hans_interface_multichain | True | False | True | hans_interface_wt_reference |
| hans_interface_wt_reference | Fe | hans_interface_multichain | True | False | True | hans_interface_wt_reference |
| hans_interface_wt_reference | Nd | hans_interface_multichain | True | False | True | hans_interface_wt_reference |
| hans_interface_wt_reference | Y | hans_interface_multichain | True | False | True | hans_interface_wt_reference |
| hans_interface_ss_plus_if_u02 | Al | hans_interface_multichain | True | True | True | hans_interface_wt_reference |
| hans_interface_ss_plus_if_u02 | Dy | hans_interface_multichain | True | True | True | hans_interface_wt_reference |
| hans_interface_ss_plus_if_u02 | Fe | hans_interface_multichain | True | True | True | hans_interface_wt_reference |
| hans_interface_ss_plus_if_u02 | Nd | hans_interface_multichain | True | True | True | hans_interface_wt_reference |
| hans_interface_ss_plus_if_u02 | Y | hans_interface_multichain | True | True | True | hans_interface_wt_reference |

## Interpretation

- Wild-type reference systems already preserve their metal-proximal waters and can move straight into generic protein/water setup planning.
- Designed systems with zero preserved solvent are compared against the matching wild-type/template family so optional metal-proximal waters can be restored deterministically before bulk solvation.
- All systems still require custom metal parameters before an actual OpenMM `System` object can be created.
