# OpenMM Metadynamics Smoke

Phase 6C1 runs one short deterministic well-tempered metadynamics triage trajectory per reduced-panel system to test whether local pocket binding is retained under bias.

## Protocol

- reduced-panel systems attempted: `25`
- successful systems: `25`
- failed systems: `0`
- panel members: `am1_mex_ss_only_u02, am1_mex_wt_reference, hans_pocket_ss_only_u02, hans_interface_ss_plus_if_u04, hans_interface_wt_reference`
- target metals: `Dy, Nd, Y, Al, Fe`
- duration per system: `500.0 ps`
- temperature: `298 K`
- timestep: `2.0 fs`
- reporter stride: `500` steps (`1.0 ps`)
- Gaussian deposition stride: `250` steps (`0.5 ps`)
- bias factor: `8.0`
- Gaussian height: `1.00 kJ/mol`
- local-pocket donor cutoff: `3.2 A`
- escape heuristic: `coordination <= 5.0` and `mean distance >= 3.2 A`

## Status Counts

- `retained_bound`: `15`
- `escape_or_weak_binding`: `0`
- `persistent_capture`: `10`
- `escaped_under_bias`: `0`
- `failed`: `0`

## Per-System Summary

| panel_member_id | target_metal | success_status | minimum_coordination_number | maximum_mean_metal_oxygen_distance_A | escape_event_detected | final_coordination_number | final_mean_metal_oxygen_distance_A | metadynamics_status |
| --- | --- | --- | ---: | ---: | --- | ---: | ---: | --- |
| am1_mex_ss_only_u02 | Dy | success | 7.579774 | 2.369706 | FALSE | 7.716559 | 2.323844 | retained_bound |
| am1_mex_ss_only_u02 | Nd | success | 8.033158 | 2.473268 | FALSE | 8.470017 | 2.435236 | retained_bound |
| am1_mex_ss_only_u02 | Y | success | 7.927117 | 2.351134 | FALSE | 8.341104 | 2.309657 | retained_bound |
| am1_mex_ss_only_u02 | Al | success | 6.219741 | 1.953113 | FALSE | 6.466359 | 1.901380 | persistent_capture |
| am1_mex_ss_only_u02 | Fe | success | 5.733933 | 2.110943 | FALSE | 6.648948 | 2.055731 | persistent_capture |
| am1_mex_wt_reference | Dy | success | 7.570040 | 2.374114 | FALSE | 7.688422 | 2.345194 | retained_bound |
| am1_mex_wt_reference | Nd | success | 8.083591 | 2.479514 | FALSE | 8.654737 | 2.438195 | retained_bound |
| am1_mex_wt_reference | Y | success | 7.912337 | 2.370116 | FALSE | 8.115236 | 2.319582 | retained_bound |
| am1_mex_wt_reference | Al | success | 6.215909 | 1.904032 | FALSE | 6.321990 | 1.881529 | persistent_capture |
| am1_mex_wt_reference | Fe | success | 6.570415 | 2.138048 | FALSE | 6.734467 | 2.094306 | persistent_capture |
| hans_pocket_ss_only_u02 | Dy | success | 8.450577 | 2.394182 | FALSE | 8.662168 | 2.354380 | retained_bound |
| hans_pocket_ss_only_u02 | Nd | success | 8.685505 | 2.462077 | FALSE | 9.361835 | 2.435706 | retained_bound |
| hans_pocket_ss_only_u02 | Y | success | 8.055975 | 2.361767 | FALSE | 8.241793 | 2.342883 | retained_bound |
| hans_pocket_ss_only_u02 | Al | success | 5.827908 | 1.892757 | FALSE | 5.896461 | 1.838681 | persistent_capture |
| hans_pocket_ss_only_u02 | Fe | success | 6.793803 | 2.127505 | FALSE | 7.253453 | 2.059346 | persistent_capture |
| hans_interface_ss_plus_if_u04 | Dy | success | 8.497619 | 2.347378 | FALSE | 8.592645 | 2.333689 | retained_bound |
| hans_interface_ss_plus_if_u04 | Nd | success | 8.861396 | 2.444889 | FALSE | 9.020901 | 2.426476 | retained_bound |
| hans_interface_ss_plus_if_u04 | Y | success | 8.261168 | 2.335582 | FALSE | 8.339267 | 2.318138 | retained_bound |
| hans_interface_ss_plus_if_u04 | Al | success | 6.328784 | 1.901782 | FALSE | 6.407392 | 1.878182 | persistent_capture |
| hans_interface_ss_plus_if_u04 | Fe | success | 7.062326 | 2.088424 | FALSE | 7.304961 | 2.071761 | persistent_capture |
| hans_interface_wt_reference | Dy | success | 8.005960 | 2.355922 | FALSE | 8.069297 | 2.336961 | retained_bound |
| hans_interface_wt_reference | Nd | success | 8.900506 | 2.448768 | FALSE | 9.048840 | 2.431462 | retained_bound |
| hans_interface_wt_reference | Y | success | 8.287275 | 2.338644 | FALSE | 8.375759 | 2.324034 | retained_bound |
| hans_interface_wt_reference | Al | success | 6.332740 | 1.890833 | FALSE | 6.425909 | 1.875584 | persistent_capture |
| hans_interface_wt_reference | Fe | success | 7.005078 | 2.092269 | FALSE | 7.051648 | 2.087419 | persistent_capture |

## Panel Status

| panel_member_id | Dy | Nd | Y | Al | Fe | candidate_keep_for_qm |
| --- | --- | --- | --- | --- | --- | --- |
| am1_mex_ss_only_u02 | retained_bound | retained_bound | retained_bound | persistent_capture | persistent_capture | FALSE |
| am1_mex_wt_reference | retained_bound | retained_bound | retained_bound | persistent_capture | persistent_capture | FALSE |
| hans_pocket_ss_only_u02 | retained_bound | retained_bound | retained_bound | persistent_capture | persistent_capture | FALSE |
| hans_interface_ss_plus_if_u04 | retained_bound | retained_bound | retained_bound | persistent_capture | persistent_capture | FALSE |
| hans_interface_wt_reference | retained_bound | retained_bound | retained_bound | persistent_capture | persistent_capture | FALSE |
