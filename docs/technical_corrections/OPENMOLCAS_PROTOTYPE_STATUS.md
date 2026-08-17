# OpenMolcas prototype status

## Validated

- OpenMolcas v26.06 source commit `8355057f32d65706a35996b5ab07cac2962bb728`.
- H2 reference: SCF, RASSCF, and CASPT2 completed; energies matched official embedded references within 1e-8 Eh.
- Atomic Dy(III) official-test derivative: explicit RICD/CDTH 1e-6, DKH/AMFI, ANO-RCC 7s6p4d2f1g, CAS(9,7), sextet, 11 roots, RASSI/SOC with 66 states, and SINGLE_ANISO.
- Independent clean repeat: root energies, occupations, SOC energies, and principal ground-doublet g values identical at printed precision.

## Repository evidence

- `results/openmolcas_prototype/raw_outputs/` — complete successful console logs.
- `results/openmolcas_prototype/figures/` — exact-data prototype figures.
- `results/openmolcas_prototype/tables/` — extracted occupations, energies, g tensor, validation, and repeat tables.
- `results/openmolcas_prototype/provenance/` — build and input-change records.
- `results/openmolcas_prototype/orbitals/` and `anisotropy/` — selected source artifacts.
- `inputs/openmolcas_prototype/` and `scripts/openmolcas_prototype/` — exact inputs and extraction/rerun code.

## Not validated for the manuscript

The prototype is an isolated atomic Dy test with an external point-charge field. It is not the manuscript's Dy-LanM cluster, there is no project-specific Nd-LanM calculation, and there is no project Dy/Nd CASPT2 comparison or reduced-cluster sensitivity study. Prototype figures must not be called manuscript Figures 6 or 7.
