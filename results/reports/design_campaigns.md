# Design Campaigns

Phase 3A deterministic ProteinMPNN input preparation from the existing Phase 2 masks.

## Derived design sets

| design_set | count | canonical_positions | am1_positions | definition |
| --- | --- | --- | --- | --- |
| hard_fixed | 43 | 1,11,14,15,16,17,19,20,21,23,25,26,28,30,31,33,41,42,43,45,46,47,49,52,55,66,67,68,70,71,72,74,77,90,91,92,94,95,96,101,105,107,110 | 1,13,14,15,17,18,19,21,23,24,26,29,37,38,39,41,42,43,45,48,51,62,63,64,66,67,68,70,73,86,87,88,90,91,92,97,101,103 | fixed_first_shell union protected_positions |
| campaign_ss_only | 16 | 13,18,22,24,40,44,48,65,69,73,83,89,93,97,98,100 | 12,16,20,22,36,40,44,61,65,69,79,85,89,93,94,96 | mutable_second_sphere minus hard_fixed |
| campaign_ss_plus_if | 23 | 13,18,22,24,27,38,39,40,44,48,65,69,73,75,76,79,80,83,89,93,97,98,100 | 12,16,20,22,25,34,35,36,40,44,61,65,69,71,72,75,76,79,85,89,93,94,96 | union(mutable_second_sphere, mutable_interface) minus hard_fixed |
| campaign_if_only | 10 | 18,27,38,39,44,75,76,79,80,83 | 16,25,34,35,40,71,72,75,76,79 | mutable_interface minus hard_fixed |

## Backbone manifest

| backbone_id | structure_id | selected_chains | design_chains | fixed_context_chains | design_chain_residue_count | note |
| --- | --- | --- | --- | --- | --- | --- |
| am1_mex_8fns_chain_a | 8FNS | A | A | - | 105 | Phase 3A default backbone specified by the workflow: 8FNS chain A. |
| hans_pocket_8fnr_chain_a | 8FNR | A | A | - | 110 | Deterministic 8FNR representative chain selected by descending metal_site_count, descending residue_count, then ascending chain_id; this chooses chain A. |
| hans_interface_8fnr_a_b_c_d | 8FNR | A,B,C,D | A | B,C,D | 110 | Full 8FNR assembly exported with representative design chain A and the remaining chains retained as fixed interface context. |

## Campaign manifest

| campaign_id | backbone_id | design_set | designed_chains | fixed_context_chains | designable | fixed |
| --- | --- | --- | --- | --- | --- | --- |
| am1_mex_ss_only | am1_mex_8fns_chain_a | campaign_ss_only | A | - | 16 | 89 |
| hans_pocket_ss_only | hans_pocket_8fnr_chain_a | campaign_ss_only | A | - | 16 | 94 |
| hans_interface_ss_plus_if | hans_interface_8fnr_a_b_c_d | campaign_ss_plus_if | A | B,C,D | 23 | 87 |
| hans_interface_if_only | hans_interface_8fnr_a_b_c_d | campaign_if_only | A | B,C,D | 10 | 100 |

ProteinMPNN note: `fixed_positions.jsonl` uses 1-based sequence indices on the exported design chains. Use `design_campaign_positions.csv` to translate those indices back to residue numbers and canonical family positions.

## am1_mex_ss_only

- Backbone: `am1_mex_8fns_chain_a`
- Design set: `campaign_ss_only`
- Designed chains: `A`
- Fixed context chains: `-`

| chain | seq_index | residue_seq | residue_name | canonical | am1 | why_designable |
| --- | --- | --- | --- | --- | --- | --- |
| A | 6 | 34 | PHE | 13 | 12 | mutable_second_sphere because second-sphere contact observed 2x across 6MI5, 8FNS (closest 4.035 A) |
| A | 10 | 38 | LYS | 18 | 16 | mutable_second_sphere because second-sphere contact observed 7x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 4.534 A); mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.031 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 14 | 42 | ILE | 22 | 20 | mutable_second_sphere because second-sphere contact observed 8x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 3.951 A) |
| A | 16 | 44 | LEU | 24 | 22 | mutable_second_sphere because second-sphere contact observed 6x across 8DQ2, 8FNR (closest 4.204 A) |
| A | 30 | 58 | LEU | 40 | 36 | mutable_second_sphere because second-sphere contact observed 1x across 8FNS (closest 5.529 A) |
| A | 34 | 62 | LYS | 44 | 40 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 5.107 A); mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.423 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 38 | 66 | LEU | 48 | 44 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 3.895 A) |
| A | 55 | 83 | LEU | 65 | 61 | mutable_second_sphere because second-sphere contact observed 1x across 8DQ2 (closest 5.967 A) |
| A | 59 | 87 | ASN | 69 | 65 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 5.120 A) |
| A | 63 | 91 | LEU | 73 | 69 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 3.749 A) |
| A | 73 | 101 | GLU | 83 | 79 | mutable_second_sphere because second-sphere contact observed 8x across 8DQ2, 8FNR (closest 4.490 A); mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 2.459 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 79 | 107 | ALA | 89 | 85 | mutable_second_sphere because second-sphere contact observed 4x across 8FNR (closest 5.290 A) |
| A | 83 | 111 | ASN | 93 | 89 | mutable_second_sphere because second-sphere contact observed 5x across 8FNR, 8FNS (closest 4.538 A) |
| A | 87 | 115 | ILE | 97 | 93 | mutable_second_sphere because second-sphere contact observed 5x across 8FNR, 8FNS (closest 4.370 A) |
| A | 88 | 116 | ASP | 98 | 94 | mutable_second_sphere because second-sphere contact observed 5x across 8FNR, 8FNS (closest 4.728 A) |
| A | 90 | 118 | ARG | 100 | 96 | mutable_second_sphere because second-sphere contact observed 1x across 8FNS (closest 5.690 A) |

## hans_pocket_ss_only

- Backbone: `hans_pocket_8fnr_chain_a`
- Design set: `campaign_ss_only`
- Designed chains: `A`
- Fixed context chains: `-`

| chain | seq_index | residue_seq | residue_name | canonical | am1 | why_designable |
| --- | --- | --- | --- | --- | --- | --- |
| A | 7 | 30 | LEU | 13 | 12 | mutable_second_sphere because second-sphere contact observed 2x across 6MI5, 8FNS (closest 4.035 A) |
| A | 12 | 35 | LYS | 18 | 16 | mutable_second_sphere because second-sphere contact observed 7x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 4.534 A); mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.031 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 16 | 39 | ASP | 22 | 20 | mutable_second_sphere because second-sphere contact observed 8x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 3.951 A) |
| A | 18 | 41 | LEU | 24 | 22 | mutable_second_sphere because second-sphere contact observed 6x across 8DQ2, 8FNR (closest 4.204 A) |
| A | 34 | 57 | ILE | 40 | 36 | mutable_second_sphere because second-sphere contact observed 1x across 8FNS (closest 5.529 A) |
| A | 38 | 61 | GLY | 44 | 40 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 5.107 A); mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.423 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 42 | 65 | LEU | 48 | 44 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 3.895 A) |
| A | 59 | 82 | ALA | 65 | 61 | mutable_second_sphere because second-sphere contact observed 1x across 8DQ2 (closest 5.967 A) |
| A | 63 | 86 | GLY | 69 | 65 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 5.120 A) |
| A | 67 | 90 | LEU | 73 | 69 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 3.749 A) |
| A | 77 | 100 | ARG | 83 | 79 | mutable_second_sphere because second-sphere contact observed 8x across 8DQ2, 8FNR (closest 4.490 A); mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 2.459 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 83 | 106 | ALA | 89 | 85 | mutable_second_sphere because second-sphere contact observed 4x across 8FNR (closest 5.290 A) |
| A | 87 | 110 | LYS | 93 | 89 | mutable_second_sphere because second-sphere contact observed 5x across 8FNR, 8FNS (closest 4.538 A) |
| A | 91 | 114 | LEU | 97 | 93 | mutable_second_sphere because second-sphere contact observed 5x across 8FNR, 8FNS (closest 4.370 A) |
| A | 92 | 115 | THR | 98 | 94 | mutable_second_sphere because second-sphere contact observed 5x across 8FNR, 8FNS (closest 4.728 A) |
| A | 94 | 117 | ALA | 100 | 96 | mutable_second_sphere because second-sphere contact observed 1x across 8FNS (closest 5.690 A) |

## hans_interface_ss_plus_if

- Backbone: `hans_interface_8fnr_a_b_c_d`
- Design set: `campaign_ss_plus_if`
- Designed chains: `A`
- Fixed context chains: `B,C,D`

| chain | seq_index | residue_seq | residue_name | canonical | am1 | why_designable |
| --- | --- | --- | --- | --- | --- | --- |
| A | 7 | 30 | LEU | 13 | 12 | mutable_second_sphere because second-sphere contact observed 2x across 6MI5, 8FNS (closest 4.035 A) |
| A | 12 | 35 | LYS | 18 | 16 | mutable_second_sphere because second-sphere contact observed 7x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 4.534 A); mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.031 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 16 | 39 | ASP | 22 | 20 | mutable_second_sphere because second-sphere contact observed 8x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 3.951 A) |
| A | 18 | 41 | LEU | 24 | 22 | mutable_second_sphere because second-sphere contact observed 6x across 8DQ2, 8FNR (closest 4.204 A) |
| A | 21 | 44 | ALA | 27 | 25 | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 3.613 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 32 | 55 | THR | 38 | 34 | mutable_interface because Hans interchain contact observed 5x across 8DQ2, 8FNR (closest 3.031 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 33 | 56 | ALA | 39 | 35 | mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.263 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 34 | 57 | ILE | 40 | 36 | mutable_second_sphere because second-sphere contact observed 1x across 8FNS (closest 5.529 A) |
| A | 38 | 61 | GLY | 44 | 40 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 5.107 A); mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.423 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 42 | 65 | LEU | 48 | 44 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 3.895 A) |
| A | 59 | 82 | ALA | 65 | 61 | mutable_second_sphere because second-sphere contact observed 1x across 8DQ2 (closest 5.967 A) |
| A | 63 | 86 | GLY | 69 | 65 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 5.120 A) |
| A | 67 | 90 | LEU | 73 | 69 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 3.749 A) |
| A | 69 | 92 | MET | 75 | 71 | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 3.110 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 70 | 93 | ASP | 76 | 72 | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 2.459 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 73 | 96 | LEU | 79 | 75 | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 3.506 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 74 | 97 | LYS | 80 | 76 | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 4.700 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 77 | 100 | ARG | 83 | 79 | mutable_second_sphere because second-sphere contact observed 8x across 8DQ2, 8FNR (closest 4.490 A); mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 2.459 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 83 | 106 | ALA | 89 | 85 | mutable_second_sphere because second-sphere contact observed 4x across 8FNR (closest 5.290 A) |
| A | 87 | 110 | LYS | 93 | 89 | mutable_second_sphere because second-sphere contact observed 5x across 8FNR, 8FNS (closest 4.538 A) |
| A | 91 | 114 | LEU | 97 | 93 | mutable_second_sphere because second-sphere contact observed 5x across 8FNR, 8FNS (closest 4.370 A) |
| A | 92 | 115 | THR | 98 | 94 | mutable_second_sphere because second-sphere contact observed 5x across 8FNR, 8FNS (closest 4.728 A) |
| A | 94 | 117 | ALA | 100 | 96 | mutable_second_sphere because second-sphere contact observed 1x across 8FNS (closest 5.690 A) |

## hans_interface_if_only

- Backbone: `hans_interface_8fnr_a_b_c_d`
- Design set: `campaign_if_only`
- Designed chains: `A`
- Fixed context chains: `B,C,D`

| chain | seq_index | residue_seq | residue_name | canonical | am1 | why_designable |
| --- | --- | --- | --- | --- | --- | --- |
| A | 12 | 35 | LYS | 18 | 16 | mutable_second_sphere because second-sphere contact observed 7x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 4.534 A); mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.031 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 21 | 44 | ALA | 27 | 25 | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 3.613 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 32 | 55 | THR | 38 | 34 | mutable_interface because Hans interchain contact observed 5x across 8DQ2, 8FNR (closest 3.031 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 33 | 56 | ALA | 39 | 35 | mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.263 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 38 | 61 | GLY | 44 | 40 | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 5.107 A); mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.423 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 69 | 92 | MET | 75 | 71 | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 3.110 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 70 | 93 | ASP | 76 | 72 | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 2.459 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 73 | 96 | LEU | 79 | 75 | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 3.506 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 74 | 97 | LYS | 80 | 76 | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 4.700 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| A | 77 | 100 | ARG | 83 | 79 | mutable_second_sphere because second-sphere contact observed 8x across 8DQ2, 8FNR (closest 4.490 A); mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 2.459 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
