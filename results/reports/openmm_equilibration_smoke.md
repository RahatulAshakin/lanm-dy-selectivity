# OpenMM Equilibration Smoke

Phase 6B1 runs one deterministic minimization plus short NVT/NPT equilibration smoke replicate for every successful Phase 6A4 baseline OpenMM build.

## Smoke Protocol

- baseline-built systems attempted: `30`
- successful smoke runs: `30`
- failures recorded: `0`
- deterministic smoke replicate count per system: `1`
- target temperature: `298 K`
- pressure target: `1.0 atm`
- timestep: `2.0 fs`
- NVT duration: `20.0 ps`
- NPT duration: `20.0 ps`
- NPT note: baseline systems are serialized as non-periodic `NoCutoff` systems, so Phase 6B1 promotes only the pressure-coupled leg to a deterministic orthorhombic periodic box.

| panel_member_id | target_metal | topology_class | success_status | minimized_potential_energy | final_nvt_potential_energy | final_npt_potential_energy | final_temperature_K | final_box_volume_nm3 | failure_reason |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| am1_mex_ss_only_u02 | Al | am1_monomer | success | -31031.643044 | -28716.918857 | -8369.615494 | 304.687321 | 160.327015 | - |
| am1_mex_ss_only_u02 | Dy | am1_monomer | success | -28416.369998 | -25784.243530 | -5313.145510 | 285.003095 | 165.539764 | - |
| am1_mex_ss_only_u02 | Fe | am1_monomer | success | -29693.706780 | -27532.960494 | -7101.583143 | 293.802154 | 170.758326 | - |
| am1_mex_ss_only_u02 | Nd | am1_monomer | success | -27747.392345 | -25194.120498 | -4842.620160 | 294.029430 | 173.630800 | - |
| am1_mex_ss_only_u02 | Y | am1_monomer | success | -28508.769877 | -25820.524166 | -5431.124562 | 296.584490 | 165.936874 | - |
| am1_mex_wt_reference | Al | am1_monomer | success | -33177.726937 | -30355.553310 | -10000.228924 | 302.182591 | 183.962564 | - |
| am1_mex_wt_reference | Dy | am1_monomer | success | -30341.698069 | -27268.949287 | -7153.881905 | 297.648039 | 182.377437 | - |
| am1_mex_wt_reference | Fe | am1_monomer | success | -31817.458298 | -29014.963327 | -8665.817039 | 301.450924 | 208.006205 | - |
| am1_mex_wt_reference | Nd | am1_monomer | success | -29834.152754 | -26580.475352 | -6401.000643 | 302.969445 | 191.157665 | - |
| am1_mex_wt_reference | Y | am1_monomer | success | -30550.605486 | -27423.445938 | -6906.232844 | 299.295078 | 180.347647 | - |
| hans_pocket_ss_only_u02 | Al | hans_monomer | success | -31487.191030 | -29287.683568 | -10209.751993 | 293.907784 | 203.555219 | - |
| hans_pocket_ss_only_u02 | Dy | hans_monomer | success | -28657.310278 | -26180.624779 | -7207.838924 | 292.758751 | 195.963535 | - |
| hans_pocket_ss_only_u02 | Fe | hans_monomer | success | -30320.985378 | -27762.406530 | -8571.660365 | 292.480573 | 171.853035 | - |
| hans_pocket_ss_only_u02 | Nd | hans_monomer | success | -27971.109201 | -25740.148104 | -6487.259980 | 302.691604 | 180.683506 | - |
| hans_pocket_ss_only_u02 | Y | hans_monomer | success | -28692.500977 | -26448.361384 | -7414.093052 | 300.947100 | 199.449801 | - |
| hans_interface_ss_plus_if_u04 | Al | hans_interface_multichain | success | -105661.102078 | -100929.953562 | -36929.816803 | 301.755319 | 569.909049 | - |
| hans_interface_ss_plus_if_u04 | Dy | hans_interface_multichain | success | -96396.003902 | -89710.133934 | -26060.768471 | 299.032246 | 526.230112 | - |
| hans_interface_ss_plus_if_u04 | Fe | hans_interface_multichain | success | -101595.147031 | -96152.185141 | -31419.445084 | 299.187938 | 552.980344 | - |
| hans_interface_ss_plus_if_u04 | Nd | hans_interface_multichain | success | -83793.077228 | -87874.119475 | -23968.894983 | 298.857808 | 522.266786 | - |
| hans_interface_ss_plus_if_u04 | Y | hans_interface_multichain | success | -96063.396077 | -90243.717420 | -26460.528003 | 297.463384 | 542.693315 | - |
| hans_interface_wt_reference | Al | hans_interface_multichain | success | -108812.592094 | -103210.723641 | -37014.938091 | 303.920685 | 597.733127 | - |
| hans_interface_wt_reference | Dy | hans_interface_multichain | success | -98903.262038 | -91668.162332 | -26044.198745 | 297.692296 | 510.260448 | - |
| hans_interface_wt_reference | Fe | hans_interface_multichain | success | -103898.428723 | -97916.414998 | -32194.662201 | 299.579801 | 497.049279 | - |
| hans_interface_wt_reference | Nd | hans_interface_multichain | success | -96666.459307 | -89365.488801 | -24314.715352 | 297.185712 | 527.953655 | - |
| hans_interface_wt_reference | Y | hans_interface_multichain | success | -99664.176665 | -92173.688605 | -26442.081835 | 296.472722 | 548.136356 | - |
| hans_interface_ss_plus_if_u02 | Al | hans_interface_multichain | success | -106833.523045 | -101214.485603 | -37686.590002 | 299.451130 | 579.142486 | - |
| hans_interface_ss_plus_if_u02 | Dy | hans_interface_multichain | success | -96569.560682 | -90070.144248 | -25989.396730 | 296.162858 | 520.357088 | - |
| hans_interface_ss_plus_if_u02 | Fe | hans_interface_multichain | success | -102258.950172 | -95979.361186 | -31906.918032 | 296.491628 | 565.081655 | - |
| hans_interface_ss_plus_if_u02 | Nd | hans_interface_multichain | success | -95046.751314 | -87703.319492 | -24649.529040 | 302.752989 | 532.596630 | - |
| hans_interface_ss_plus_if_u02 | Y | hans_interface_multichain | success | -96735.103467 | -90121.594478 | -26135.201938 | 300.628796 | 533.565459 | - |
