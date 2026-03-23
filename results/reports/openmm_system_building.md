# OpenMM System Building

Phase 6A4 builds deterministic baseline OpenMM `System` objects and simulation configs for all baseline-ready MD panel systems.

## Build Summary

- baseline-ready systems attempted: `30`
- successful OpenMM `System` serializations: `30`
- failures recorded: `0`
- force-field files: `amber19-all.xml`, `amber19/opc3.xml`
- baseline metal parameter family: `generic_12_6_4_highly_charged`
- integrator defaults: `LangevinMiddleIntegrator`, `298 K`, `1.0 ps^-1`, `2.0 fs`
- nonbonded defaults: `NoCutoff`, cutoff `1.0 nm`, `hydrogen_mass_repartitioning=False`

| panel_member_id | target_metal | topology_class | success_status | atom_count | residue_count | failure_reason |
| --- | --- | --- | --- | ---: | ---: | --- |
| am1_mex_ss_only_u02 | Al | am1_monomer | success | 1608 | 134 | - |
| am1_mex_ss_only_u02 | Dy | am1_monomer | success | 1608 | 134 | - |
| am1_mex_ss_only_u02 | Fe | am1_monomer | success | 1608 | 134 | - |
| am1_mex_ss_only_u02 | Nd | am1_monomer | success | 1608 | 134 | - |
| am1_mex_ss_only_u02 | Y | am1_monomer | success | 1608 | 134 | - |
| am1_mex_wt_reference | Al | am1_monomer | success | 1623 | 134 | - |
| am1_mex_wt_reference | Dy | am1_monomer | success | 1623 | 134 | - |
| am1_mex_wt_reference | Fe | am1_monomer | success | 1623 | 134 | - |
| am1_mex_wt_reference | Nd | am1_monomer | success | 1623 | 134 | - |
| am1_mex_wt_reference | Y | am1_monomer | success | 1623 | 134 | - |
| hans_pocket_ss_only_u02 | Al | hans_monomer | success | 1676 | 126 | - |
| hans_pocket_ss_only_u02 | Dy | hans_monomer | success | 1676 | 126 | - |
| hans_pocket_ss_only_u02 | Fe | hans_monomer | success | 1676 | 126 | - |
| hans_pocket_ss_only_u02 | Nd | hans_monomer | success | 1676 | 126 | - |
| hans_pocket_ss_only_u02 | Y | hans_monomer | success | 1676 | 126 | - |
| hans_interface_ss_plus_if_u04 | Al | hans_interface_multichain | success | 6697 | 489 | - |
| hans_interface_ss_plus_if_u04 | Dy | hans_interface_multichain | success | 6697 | 489 | - |
| hans_interface_ss_plus_if_u04 | Fe | hans_interface_multichain | success | 6697 | 489 | - |
| hans_interface_ss_plus_if_u04 | Nd | hans_interface_multichain | success | 6697 | 489 | - |
| hans_interface_ss_plus_if_u04 | Y | hans_interface_multichain | success | 6697 | 489 | - |
| hans_interface_wt_reference | Al | hans_interface_multichain | success | 6716 | 489 | - |
| hans_interface_wt_reference | Dy | hans_interface_multichain | success | 6716 | 489 | - |
| hans_interface_wt_reference | Fe | hans_interface_multichain | success | 6716 | 489 | - |
| hans_interface_wt_reference | Nd | hans_interface_multichain | success | 6716 | 489 | - |
| hans_interface_wt_reference | Y | hans_interface_multichain | success | 6716 | 489 | - |
| hans_interface_ss_plus_if_u02 | Al | hans_interface_multichain | success | 6688 | 489 | - |
| hans_interface_ss_plus_if_u02 | Dy | hans_interface_multichain | success | 6688 | 489 | - |
| hans_interface_ss_plus_if_u02 | Fe | hans_interface_multichain | success | 6688 | 489 | - |
| hans_interface_ss_plus_if_u02 | Nd | hans_interface_multichain | success | 6688 | 489 | - |
| hans_interface_ss_plus_if_u02 | Y | hans_interface_multichain | success | 6688 | 489 | - |
