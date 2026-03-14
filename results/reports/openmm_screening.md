# OpenMM Screening

Phase 6B2 runs short deterministic 3-replicate NPT screening trajectories for every successful Phase 6B1 system to score stable rare-earth binding and persistent Al/Fe capture.

## Screening Protocol

- screened systems: `30`
- replicate attempts: `90`
- successful replicates: `90`
- failed replicates: `0`
- replicate schedule: `1->101`, `2->102`, `3->103`
- target temperature: `298 K`
- pressure target: `1.0 atm`
- timestep: `2.0 fs`
- production duration per replicate: `100.0 ps`
- reporter stride: `500` steps (`1.0 ps`)
- inner-sphere cutoff: `3.2 A`

## Status Counts

- `stable_bound`: `18`
- `weak_binding_flag`: `0`
- `persistent_capture_flag`: `12`
- `no_persistent_capture`: `0`
- `failed`: `0`

## Per-System Summary

| panel_member_id | target_metal | topology_class | replicate_success_count | median_mean_min_metal_oxygen_distance_A | median_inner_sphere_occupancy_fraction | screening_status |
| --- | --- | --- | ---: | ---: | ---: | --- |
| am1_mex_ss_only_u02 | Al | am1_monomer | 3 | 1.721077 | 1.000000 | persistent_capture_flag |
| am1_mex_ss_only_u02 | Dy | am1_monomer | 3 | 2.201175 | 1.000000 | stable_bound |
| am1_mex_ss_only_u02 | Fe | am1_monomer | 3 | 1.920004 | 1.000000 | persistent_capture_flag |
| am1_mex_ss_only_u02 | Nd | am1_monomer | 3 | 2.302743 | 1.000000 | stable_bound |
| am1_mex_ss_only_u02 | Y | am1_monomer | 3 | 2.197013 | 1.000000 | stable_bound |
| am1_mex_wt_reference | Al | am1_monomer | 3 | 1.724824 | 1.000000 | persistent_capture_flag |
| am1_mex_wt_reference | Dy | am1_monomer | 3 | 2.203971 | 1.000000 | stable_bound |
| am1_mex_wt_reference | Fe | am1_monomer | 3 | 1.935318 | 1.000000 | persistent_capture_flag |
| am1_mex_wt_reference | Nd | am1_monomer | 3 | 2.303163 | 1.000000 | stable_bound |
| am1_mex_wt_reference | Y | am1_monomer | 3 | 2.194080 | 1.000000 | stable_bound |
| hans_interface_ss_plus_if_u02 | Al | hans_interface_multichain | 3 | 1.726544 | 1.000000 | persistent_capture_flag |
| hans_interface_ss_plus_if_u02 | Dy | hans_interface_multichain | 3 | 2.209440 | 1.000000 | stable_bound |
| hans_interface_ss_plus_if_u02 | Fe | hans_interface_multichain | 3 | 1.928160 | 1.000000 | persistent_capture_flag |
| hans_interface_ss_plus_if_u02 | Nd | hans_interface_multichain | 3 | 2.315071 | 1.000000 | stable_bound |
| hans_interface_ss_plus_if_u02 | Y | hans_interface_multichain | 3 | 2.203336 | 1.000000 | stable_bound |
| hans_interface_ss_plus_if_u04 | Al | hans_interface_multichain | 3 | 1.723012 | 1.000000 | persistent_capture_flag |
| hans_interface_ss_plus_if_u04 | Dy | hans_interface_multichain | 3 | 2.207140 | 1.000000 | stable_bound |
| hans_interface_ss_plus_if_u04 | Fe | hans_interface_multichain | 3 | 1.927614 | 1.000000 | persistent_capture_flag |
| hans_interface_ss_plus_if_u04 | Nd | hans_interface_multichain | 3 | 2.312553 | 1.000000 | stable_bound |
| hans_interface_ss_plus_if_u04 | Y | hans_interface_multichain | 3 | 2.190667 | 1.000000 | stable_bound |
| hans_interface_wt_reference | Al | hans_interface_multichain | 3 | 1.725407 | 1.000000 | persistent_capture_flag |
| hans_interface_wt_reference | Dy | hans_interface_multichain | 3 | 2.203732 | 1.000000 | stable_bound |
| hans_interface_wt_reference | Fe | hans_interface_multichain | 3 | 1.927193 | 1.000000 | persistent_capture_flag |
| hans_interface_wt_reference | Nd | hans_interface_multichain | 3 | 2.310992 | 1.000000 | stable_bound |
| hans_interface_wt_reference | Y | hans_interface_multichain | 3 | 2.195130 | 1.000000 | stable_bound |
| hans_pocket_ss_only_u02 | Al | hans_monomer | 3 | 1.722567 | 1.000000 | persistent_capture_flag |
| hans_pocket_ss_only_u02 | Dy | hans_monomer | 3 | 2.216472 | 1.000000 | stable_bound |
| hans_pocket_ss_only_u02 | Fe | hans_monomer | 3 | 1.928160 | 1.000000 | persistent_capture_flag |
| hans_pocket_ss_only_u02 | Nd | hans_monomer | 3 | 2.319641 | 1.000000 | stable_bound |
| hans_pocket_ss_only_u02 | Y | hans_monomer | 3 | 2.208915 | 1.000000 | stable_bound |

## Panel Status

| panel_member_id | topology_class | Dy | Nd | Y | Al | Fe |
| --- | --- | --- | --- | --- | --- | --- |
| am1_mex_ss_only_u02 | am1_monomer | stable_bound | stable_bound | stable_bound | persistent_capture_flag | persistent_capture_flag |
| am1_mex_wt_reference | am1_monomer | stable_bound | stable_bound | stable_bound | persistent_capture_flag | persistent_capture_flag |
| hans_interface_ss_plus_if_u02 | hans_interface_multichain | stable_bound | stable_bound | stable_bound | persistent_capture_flag | persistent_capture_flag |
| hans_interface_ss_plus_if_u04 | hans_interface_multichain | stable_bound | stable_bound | stable_bound | persistent_capture_flag | persistent_capture_flag |
| hans_interface_wt_reference | hans_interface_multichain | stable_bound | stable_bound | stable_bound | persistent_capture_flag | persistent_capture_flag |
| hans_pocket_ss_only_u02 | hans_monomer | stable_bound | stable_bound | stable_bound | persistent_capture_flag | persistent_capture_flag |
