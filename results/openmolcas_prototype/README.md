# LanM OpenMolcas prototype execution package v1.0

Date executed: 16 August 2026 (America/New_York)

This is an auditable **prototype**, not a corrected manuscript-results package. It contains genuine OpenMolcas v26.06 console outputs, inputs, orbital artifacts, extracted tables, exact-data plots, provenance, and SHA-256 checksums.

## Validated executions

1. **Official H2 reference test 000** — SCF, RASSCF, and CASPT2 completed with `Happy landing`; the three printed energies agree with the official embedded reference values within 1×10⁻⁸ Eh.
2. **Dy(III) workflow prototype derived from official additional test 321** — RICD/acCD integral generation, DKH/AMFI, CAS(9,7) RASSCF with 11 sextet roots, RASSI spin–orbit coupling, and SINGLE_ANISO completed with `_RC_ALL_IS_WELL_` and `Happy landing`.
3. **Independent Dy repeat** — a second clean work directory reproduced the RASSCF energies, seven-orbital occupations, spin–orbit energies, and ground-doublet g values at every printed digit.

The Dy model is an isolated test atom with the official test's external point-charge field. It contains no LanM protein, ligand shell, waters, MD frame, or manuscript reduced cluster.

## Start here

- `PROTOTYPE_RUN_REPORT.md` — scope, settings, results, validation, and limitations.
- `TABLES_AND_FIGURES_INDEX.md` — exact inventory of extracted tables and plots.
- `raw_outputs/` — unedited console logs from the successful runs and repeat.
- `diagnostics/` — preserved unsuccessful/partial conventional-integral attempts.
- `inputs/` — original official inputs plus the explicit RICD prototype derivative.
- `orbitals/` and `anisotropy/` — generated binary/text result artifacts.
- `provenance/` — software/build/environment record and input change description.
- `scripts/` — extraction, plotting, and rerun scripts.
- `CHECKSUMS.sha256` — integrity manifest for all package files.

Read `NOT_FOR_MANUSCRIPT_SUBMISSION.txt` before reusing any result.

Source project: <https://github.com/OpenMolcas/OpenMolcas>
