# Metal Parameter Strategy

Phase 6A2c preserves the Phase 6A2b conceptual family mapping, ingests published OPC3 12-6-4 baseline coefficients, and updates per-system baseline-versus-tuned build readiness.
It does not build OpenMM `System` objects, run minimization, or run MD.

## Why Phase 6A2b Was Conservative

- All `30` Phase 6A2 panel systems still carry `custom_metal_parameters_still_required=True` in `results/tables/md_system_build_manifest.csv` because that artifact predated any audited numeric metal registry.
- Phase 6A2b intentionally stopped at conceptual family mapping in `config/metal_parameters.yaml`; it did not ingest explicit Dy/Nd/Y/Al/Fe coefficient sets into the repo.
- With only family labels and no auditable coefficient registry, the conservative and correct status in Phase 6A2b was to keep every system out of build-ready state.

## Published OPC3 Baseline Registry

- Published generic OPC3 12-6-4 baseline coefficients are now present for `Dy, Nd, Y, Al, Fe` in `config/metal_parameter_values.yaml`.
- The conceptual LanM-tuned family is still the direct refinement path for `Dy, Nd, Y`, while `Al, Fe` remain conceptually generic/proxy metals.
- Result: the generic baseline family is now numerically build-ready for Dy, Nd, Y, Al, and Fe under OPC3, while the LanM-tuned family remains a future refinement path until explicit tuned coefficients are added.

| metal | source_family_label | water_model | rmin_half_A | epsilon_kcal_per_mol | c4_kcal_per_mol_A4 |
| --- | --- | --- | ---: | ---: | ---: |
| Dy | generic_12_6_4_highly_charged | amber19/opc3.xml | 1.632 | 0.09620220 | 183 |
| Nd | generic_12_6_4_highly_charged | amber19/opc3.xml | 1.712 | 0.14640930 | 184 |
| Y | generic_12_6_4_highly_charged | amber19/opc3.xml | 1.626 | 0.09289608 | 192 |
| Al | generic_12_6_4_highly_charged | amber19/opc3.xml | 1.361 | 0.01031847 | 363 |
| Fe | generic_12_6_4_highly_charged | amber19/opc3.xml | 1.455 | 0.02662782 | 429 |

## Family Status

| source_family_label | family_scope | covered_metals | repo_numeric_parameters_present | parameter_provenance | water_model_compatibility |
| --- | --- | --- | --- | --- | --- |
| generic_12_6_4_highly_charged | generic_baseline | Dy, Nd, Y, Al, Fe | True | published_opc3_12_6_4_baseline | amber19/opc3.xml |
| chelator_tuned_12_6_4_lanmodulin | lanm_adjacent_chelator_tuned | Dy, Nd, Y | False | explicit_tuned_numeric_coefficients_in_repo | amber19/opc3.xml |

## Per-Metal Mapping

| metal | intended_parameter_family | conceptual_source_family | direct_support_status | baseline_parameter_family | baseline_numeric_values_present | tuned_parameter_family | tuned_numeric_values_present |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Dy | chelator_tuned_12_6_4_lanmodulin | chelator_tuned_12_6_4_lanmodulin | direct | generic_12_6_4_highly_charged | True | chelator_tuned_12_6_4_lanmodulin | False |
| Nd | chelator_tuned_12_6_4_lanmodulin | chelator_tuned_12_6_4_lanmodulin | direct | generic_12_6_4_highly_charged | True | chelator_tuned_12_6_4_lanmodulin | False |
| Y | chelator_tuned_12_6_4_lanmodulin | chelator_tuned_12_6_4_lanmodulin | direct | generic_12_6_4_highly_charged | True | chelator_tuned_12_6_4_lanmodulin | False |
| Al | chelator_tuned_12_6_4_lanmodulin | generic_12_6_4_highly_charged | proxy | generic_12_6_4_highly_charged | True | chelator_tuned_12_6_4_lanmodulin | False |
| Fe | chelator_tuned_12_6_4_lanmodulin | generic_12_6_4_highly_charged | proxy | generic_12_6_4_highly_charged | True | chelator_tuned_12_6_4_lanmodulin | False |

## Per-System Readiness

- systems mapped: `30`
- systems baseline-ready for OpenMM `System` build: `30`
- systems tuned-ready for OpenMM `System` build: `0`
- per-system readiness audit written to `results/tables/metal_build_readiness.csv`.

| panel_member_id | target_metal | chosen_parameter_family | baseline_numeric_values_present | ready_for_openmm_system_build_baseline | tuned_numeric_values_present | ready_for_openmm_system_build_tuned | readiness_summary |
| --- | --- | --- | --- | --- | --- | --- | --- |
| am1_mex_ss_only_u02 | Dy | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| am1_mex_ss_only_u02 | Nd | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| am1_mex_ss_only_u02 | Y | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| am1_mex_ss_only_u02 | Al | generic_12_6_4_highly_charged | True | True | False | False | baseline_ready_tuned_pending |
| am1_mex_ss_only_u02 | Fe | generic_12_6_4_highly_charged | True | True | False | False | baseline_ready_tuned_pending |
| am1_mex_wt_reference | Dy | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| am1_mex_wt_reference | Nd | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| am1_mex_wt_reference | Y | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| am1_mex_wt_reference | Al | generic_12_6_4_highly_charged | True | True | False | False | baseline_ready_tuned_pending |
| am1_mex_wt_reference | Fe | generic_12_6_4_highly_charged | True | True | False | False | baseline_ready_tuned_pending |
| hans_pocket_ss_only_u02 | Dy | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| hans_pocket_ss_only_u02 | Nd | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| hans_pocket_ss_only_u02 | Y | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| hans_pocket_ss_only_u02 | Al | generic_12_6_4_highly_charged | True | True | False | False | baseline_ready_tuned_pending |
| hans_pocket_ss_only_u02 | Fe | generic_12_6_4_highly_charged | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_ss_plus_if_u04 | Dy | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_ss_plus_if_u04 | Nd | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_ss_plus_if_u04 | Y | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_ss_plus_if_u04 | Al | generic_12_6_4_highly_charged | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_ss_plus_if_u04 | Fe | generic_12_6_4_highly_charged | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_wt_reference | Dy | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_wt_reference | Nd | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_wt_reference | Y | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_wt_reference | Al | generic_12_6_4_highly_charged | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_wt_reference | Fe | generic_12_6_4_highly_charged | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_ss_plus_if_u02 | Dy | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_ss_plus_if_u02 | Nd | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_ss_plus_if_u02 | Y | chelator_tuned_12_6_4_lanmodulin | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_ss_plus_if_u02 | Al | generic_12_6_4_highly_charged | True | True | False | False | baseline_ready_tuned_pending |
| hans_interface_ss_plus_if_u02 | Fe | generic_12_6_4_highly_charged | True | True | False | False | baseline_ready_tuned_pending |

## Interpretation

- The published generic OPC3 12-6-4 baseline family `generic_12_6_4_highly_charged` is now numerically ready for all panel metals under `amber19/opc3.xml`.
- The LanM-tuned family `chelator_tuned_12_6_4_lanmodulin` remains intentionally not build-ready because `explicit_tuned_numeric_coefficients_in_repo` are absent from the repo.
- The generic baseline values come from `published_opc3_12_6_4_baseline` and are being handed off as the auditable baseline path; tuned coefficients remain a future refinement path rather than a claimed ready state.
