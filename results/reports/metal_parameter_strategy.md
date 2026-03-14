# Metal Parameter Strategy

Phase 6A2b converts the Phase 6A2 metal-model registry and system-build manifest into an explicit parameter-source mapping.
It does not create numeric metal parameters and it does not build OpenMM `System` objects.

## Why Phase 6A2 Marked Every System Custom

- All `30` Phase 6A2 panel systems still carry `custom_metal_parameters_still_required=True` in `results/tables/md_system_build_manifest.csv`.
- Dy, Nd, Y, Al, and Fe are modeled as pre-bound LanM-site metals rather than ordinary bulk-solvent ions.
- Phase 6A2 recorded only an intended custom bound-site family in `config/metal_models.yaml`; it did not add numeric Dy/Nd/Y/Al/Fe parameter files to the repo.
- Result: the protein/water scaffolding is ready to audit, but full OpenMM `System` creation remains blocked until the metal-parameter source is turned into actual numeric inputs.

## Source Family Overview

- Metals with only generic baseline support: `Al, Fe`
- Metals with LanM-adjacent chelator-tuned support: `Dy, Nd, Y`
- Metals that still require derivation before OpenMM System creation: `Dy, Nd, Y`
- Metals that still require an explicit proxy decision before OpenMM System creation: `Al, Fe`

| source_family_label | family_scope | covered_metals | repo_numeric_parameters_present | water_model_compatibility |
| --- | --- | --- | --- | --- |
| generic_12_6_4_highly_charged | generic_baseline | Al, Fe | False | amber19/opc3.xml |
| chelator_tuned_12_6_4_lanmodulin | lanm_adjacent_chelator_tuned | Dy, Nd, Y | False | amber19/opc3.xml |

## Per-Metal Mapping

| metal | formal_charge | intended_parameter_family | source_family_label | direct_support_status | water_model_compatibility |
| --- | ---: | --- | --- | --- | --- |
| Dy | +3 | chelator_tuned_12_6_4_lanmodulin | chelator_tuned_12_6_4_lanmodulin | direct | amber19/opc3.xml |
| Nd | +3 | chelator_tuned_12_6_4_lanmodulin | chelator_tuned_12_6_4_lanmodulin | direct | amber19/opc3.xml |
| Y | +3 | chelator_tuned_12_6_4_lanmodulin | chelator_tuned_12_6_4_lanmodulin | direct | amber19/opc3.xml |
| Al | +3 | chelator_tuned_12_6_4_lanmodulin | generic_12_6_4_highly_charged | proxy | amber19/opc3.xml |
| Fe | +3 | chelator_tuned_12_6_4_lanmodulin | generic_12_6_4_highly_charged | proxy | amber19/opc3.xml |

## Per-System Mapping

- systems mapped: `30`
- systems with direct family mappings: `18`
- systems with proxy family mappings: `12`
- systems ready for OpenMM `System` build: `0`

| panel_member_id | target_metal | chosen_parameter_family | direct_support_status | ready_for_openmm_system_build | rationale |
| --- | --- | --- | --- | --- | --- |
| am1_mex_ss_only_u02 | Dy | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| am1_mex_ss_only_u02 | Nd | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| am1_mex_ss_only_u02 | Y | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| am1_mex_ss_only_u02 | Al | generic_12_6_4_highly_charged | proxy | False | Phase 6A2 kept this system in the custom-parameter-required state because the current mapping for Al is a generic 12-6-4 proxy rather than a direct LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed. |
| am1_mex_ss_only_u02 | Fe | generic_12_6_4_highly_charged | proxy | False | Phase 6A2 kept this system in the custom-parameter-required state because the current mapping for Fe is a generic 12-6-4 proxy rather than a direct LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed. |
| am1_mex_wt_reference | Dy | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| am1_mex_wt_reference | Nd | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| am1_mex_wt_reference | Y | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| am1_mex_wt_reference | Al | generic_12_6_4_highly_charged | proxy | False | Phase 6A2 kept this system in the custom-parameter-required state because the current mapping for Al is a generic 12-6-4 proxy rather than a direct LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed. |
| am1_mex_wt_reference | Fe | generic_12_6_4_highly_charged | proxy | False | Phase 6A2 kept this system in the custom-parameter-required state because the current mapping for Fe is a generic 12-6-4 proxy rather than a direct LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed. |
| hans_pocket_ss_only_u02 | Dy | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| hans_pocket_ss_only_u02 | Nd | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| hans_pocket_ss_only_u02 | Y | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| hans_pocket_ss_only_u02 | Al | generic_12_6_4_highly_charged | proxy | False | Phase 6A2 kept this system in the custom-parameter-required state because the current mapping for Al is a generic 12-6-4 proxy rather than a direct LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed. |
| hans_pocket_ss_only_u02 | Fe | generic_12_6_4_highly_charged | proxy | False | Phase 6A2 kept this system in the custom-parameter-required state because the current mapping for Fe is a generic 12-6-4 proxy rather than a direct LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed. |
| hans_interface_ss_plus_if_u04 | Dy | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| hans_interface_ss_plus_if_u04 | Nd | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| hans_interface_ss_plus_if_u04 | Y | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| hans_interface_ss_plus_if_u04 | Al | generic_12_6_4_highly_charged | proxy | False | Phase 6A2 kept this system in the custom-parameter-required state because the current mapping for Al is a generic 12-6-4 proxy rather than a direct LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed. |
| hans_interface_ss_plus_if_u04 | Fe | generic_12_6_4_highly_charged | proxy | False | Phase 6A2 kept this system in the custom-parameter-required state because the current mapping for Fe is a generic 12-6-4 proxy rather than a direct LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed. |
| hans_interface_wt_reference | Dy | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| hans_interface_wt_reference | Nd | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| hans_interface_wt_reference | Y | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| hans_interface_wt_reference | Al | generic_12_6_4_highly_charged | proxy | False | Phase 6A2 kept this system in the custom-parameter-required state because the current mapping for Al is a generic 12-6-4 proxy rather than a direct LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed. |
| hans_interface_wt_reference | Fe | generic_12_6_4_highly_charged | proxy | False | Phase 6A2 kept this system in the custom-parameter-required state because the current mapping for Fe is a generic 12-6-4 proxy rather than a direct LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed. |
| hans_interface_ss_plus_if_u02 | Dy | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| hans_interface_ss_plus_if_u02 | Nd | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| hans_interface_ss_plus_if_u02 | Y | chelator_tuned_12_6_4_lanmodulin | direct | False | Phase 6A2 kept this system in the custom-parameter-required state because the chosen 'chelator_tuned_12_6_4_lanmodulin' family is only a conceptual LanM-adjacent mapping in the repo; numeric parameters still need to be derived. |
| hans_interface_ss_plus_if_u02 | Al | generic_12_6_4_highly_charged | proxy | False | Phase 6A2 kept this system in the custom-parameter-required state because the current mapping for Al is a generic 12-6-4 proxy rather than a direct LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed. |
| hans_interface_ss_plus_if_u02 | Fe | generic_12_6_4_highly_charged | proxy | False | Phase 6A2 kept this system in the custom-parameter-required state because the current mapping for Fe is a generic 12-6-4 proxy rather than a direct LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed. |
