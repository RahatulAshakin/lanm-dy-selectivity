# OpenMM System Smoke

Phase 6A3 validates one Dy baseline-ready representative per topology class using the generic OPC3 12-6-4 family path.

## Build Summary

- representative systems attempted: `3`
- successful OpenMM `System` serializations: `3`
- failures recorded: `0`
- force-field files: `amber19-all.xml`, `amber19/opc3.xml`
- custom metal family audited for this phase: `generic_12_6_4_highly_charged`

| panel_member_id | target_metal | topology_class | success_status | restored_template_waters | atom_count | residue_count | failure_reason |
| --- | --- | --- | --- | ---: | ---: | ---: | --- |
| am1_mex_ss_only_u02 | Dy | am1_monomer | success | 25 | 1608 | 134 | - |
| hans_pocket_ss_only_u02 | Dy | hans_monomer | success | 12 | 1676 | 126 | - |
| hans_interface_ss_plus_if_u04 | Dy | hans_interface_multichain | success | 40 | 6697 | 489 | - |
