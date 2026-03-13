# Template Harmonization

Phase 2A/2B/2C deterministic template parsing, cross-template alignment, residue-role summary, and design-mask generation for four lanmodulin templates.

Sequence reference records loaded: 3

AM1 sequence note:

- Swiss-Model reports full-length LanM as 133 aa; PDB-linked structures cover residues 23-133 (6MI5) and 29-133 (8FNS).

## AM1/Mex family

| template_id | source_type | chain_ids | residue_ranges | metal_identity | metal_site_count | experimental_method |
| --- | --- | --- | --- | --- | --- | --- |
| 6MI5 | cif | X | X:23-139 | Y | 3 | SOLUTION NMR |
| 8FNS | csv_atom_table | A | A:29-133 | ND | 4 | X-RAY DIFFRACTION |

- `6MI5` uses model 1 of the NMR ensemble for deterministic counting and reports X:3 metal sites.
- `8FNS` is a single-chain crystal template spanning A:29-133 with A:4 metal sites.

## Hans family

| template_id | source_type | chain_ids | residue_ranges | metal_identity | metal_site_count | experimental_method |
| --- | --- | --- | --- | --- | --- | --- |
| 8DQ2 | csv_atom_table | A, B, C, D | A:24-133; B:24-133; C:24-133; D:24-133 | LA | 12 | X-RAY DIFFRACTION |
| 8FNR | cif | A, B, C, D | A:24-133; B:24-133; C:24-133; D:24-133 | DY | 14 | X-RAY DIFFRACTION |

- `8DQ2` and `8FNR` both resolve chains A, B, C, D over the Hans family backbone.
- Metal occupancy differs across the tetramer: `8DQ2` shows A:3, B:3, C:3, D:3, while `8FNR` shows A:4, B:3, C:3, D:4.
- `8FNR` chain D retains residues 24-133 but the observed model omits residues 34-38, reducing the chain residue count relative to `8DQ2`.

## Alignment and residue roles

- AM1 mature reference length: 111 aa from `lanmodulin_sequences.csv`.
- `6MI5` chain `X` maps residues 23-133 onto AM1 mature positions 1-111; 6 C-terminal His-tag residues remain outside the mature reference.
- `8FNS` chain `A` covers AM1 mature positions 7-111 in the representative AM1 alignment.
- Hans representative chain: `8DQ2` chain `A` aligned against `6MI5` chain `X` to define 116 canonical family positions; 5 positions are Hans-only relative to AM1 mature numbering.

| template_id | aligned_rows | outside_am1_reference | first_shell | second_sphere | solvent_contact | interchain_contact |
| --- | --- | --- | --- | --- | --- | --- |
| 6MI5 | 117 | 6 | 15 | 18 | 0 | 0 |
| 8FNS | 105 | 0 | 20 | 22 | 25 | 0 |
| 8DQ2 | 440 | 0 | 72 | 56 | 29 | 78 |
| 8FNR | 435 | 0 | 77 | 69 | 40 | 98 |

## Design mask candidates

| category | canonical_position_count |
| --- | --- |
| fixed_first_shell | 26 |
| mutable_second_sphere | 16 |
| mutable_interface | 10 |
| protected_positions | 18 |

- `fixed_first_shell` spans [14, 16, 17, 19, 21, 23, 25, 26, 28, 41, 43, 45, 47, 49, 52, 66, 68, 70, 72, 74, 77, 90, 92, 94, 96, 101].
- `mutable_second_sphere` uses the top-ranked positions below after excluding fixed and protected sites.
- `mutable_interface` applies the Hans-family interchain-contact filter within +/-3 canonical positions of metal-associated residues.

### Top mutable_second_sphere positions

| canonical_family_position | am1_reference_residue | second_sphere_obs | observed_residue_identities | rationale |
| --- | --- | --- | --- | --- |
| 44 / AM1 40 | LYS | 10 | GLY,LYS | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 5.107 A); mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.423 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| 48 / AM1 44 | LEU | 10 | LEU | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 3.895 A) |
| 69 / AM1 65 | ASN | 10 | ASN,GLY | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 5.120 A) |
| 73 / AM1 69 | LEU | 10 | LEU | mutable_second_sphere because second-sphere contact observed 10x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 3.749 A) |
| 22 / AM1 20 | ILE | 8 | ASP,ILE | mutable_second_sphere because second-sphere contact observed 8x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 3.951 A) |
| 83 / AM1 79 | GLU | 8 | ARG,GLU | mutable_second_sphere because second-sphere contact observed 8x across 8DQ2, 8FNR (closest 4.490 A); mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 2.459 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| 18 / AM1 16 | LYS | 7 | LYS | mutable_second_sphere because second-sphere contact observed 7x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 4.534 A); mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.031 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| 24 / AM1 22 | LEU | 6 | LEU | mutable_second_sphere because second-sphere contact observed 6x across 8DQ2, 8FNR (closest 4.204 A) |

### Top mutable_interface positions

| canonical_family_position | am1_reference_residue | interface_obs | observed_residue_identities | rationale |
| --- | --- | --- | --- | --- |
| 27 / AM1 25 | ALA | 8 | ALA | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 3.613 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| 75 / AM1 71 | LYS | 8 | LYS,MET | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 3.110 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| 76 / AM1 72 | LYS | 8 | ASP,LYS | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 2.459 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| 79 / AM1 75 | LEU | 8 | LEU | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 3.506 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| 80 / AM1 76 | ALA | 8 | ALA,LYS | mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 4.700 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| 83 / AM1 79 | GLU | 8 | ARG,GLU | mutable_second_sphere because second-sphere contact observed 8x across 8DQ2, 8FNR (closest 4.490 A); mutable_interface because Hans interchain contact observed 8x across 8DQ2, 8FNR (closest 2.459 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| 38 / AM1 34 | ASP | 5 | ASP,THR | mutable_interface because Hans interchain contact observed 5x across 8DQ2, 8FNR (closest 3.031 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |
| 18 / AM1 16 | LYS | 3 | LYS | mutable_second_sphere because second-sphere contact observed 7x across 6MI5, 8DQ2, 8FNR, 8FNS (closest 4.534 A); mutable_interface because Hans interchain contact observed 3x across 8DQ2, 8FNR (closest 3.031 A) and the position lies within +/-3 canonical positions of a Hans metal-associated residue |

### Positions explicitly protected

| canonical_family_position | am1_reference_residue | observed_residue_identities | protection_reasons |
| --- | --- | --- | --- |
| 1 / AM1 1 | PRO | PRO | AM1 reference residue is PRO |
| 15 / AM1 14 | PRO | ALA,PRO | AM1 reference residue is PRO |
| 20 / AM1 18 | GLY | ASN,GLY | AM1 reference residue is GLY |
| 33 / AM1 29 | GLY | GLY | AM1 reference residue is GLY |
| 42 / AM1 38 | PRO | PRO | AM1 reference residue is PRO |
| 46 / AM1 42 | GLY | GLY,THR | AM1 reference residue is GLY |
| 55 / AM1 51 | GLY | GLY | AM1 reference residue is GLY |
| 67 / AM1 63 | PRO | LYS,PRO | AM1 reference residue is PRO |
| 71 / AM1 67 | GLY | GLN,GLY | AM1 reference residue is GLY |
| 91 / AM1 87 | PRO | ALA,PRO | AM1 reference residue is PRO |
| 95 / AM1 91 | GLY | GLY | AM1 reference residue is GLY |
| 105 / AM1 101 | PRO | LYS,PRO | AM1 reference residue is PRO |
| 107 / AM1 103 | GLY | GLY | AM1 reference residue is GLY |
| 11 (gap) | - | ASP | outside AM1 mature numbering; unresolved or gap-only alignment |
| 17 (gap) | - | ASN | outside AM1 mature numbering; unresolved or gap-only alignment |
| 30 (gap) | - | ILE | outside AM1 mature numbering; unresolved or gap-only alignment |
| 31 (gap) | - | HIS | outside AM1 mature numbering; unresolved or gap-only alignment |
| 110 (gap) | - | VAL | outside AM1 mature numbering; unresolved or gap-only alignment |

- Additional noncanonical protected positions: `6MI5` chain `X` residues 134-139.
