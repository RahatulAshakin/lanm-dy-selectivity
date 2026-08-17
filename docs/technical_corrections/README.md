# LanM technical corrections

This directory contains the reviewer-driven technical correction reports for manuscript `JMGM-D-26-01758`. The main manuscript is intentionally unchanged.

## Updated v2.0 package

`LanM_Updated_Technical_Corrections_Report.docx` and its PDF add genuine, independently repeated OpenMolcas prototype evidence to the original correction report. The original report is retained for audit history.

The v2.0 update includes:

- all five corrected technical figures and original correction/source-audit tables;
- revised reviewer-response and quantum-evidence status tables;
- an H2 RASSCF/CASPT2 reference validation;
- an atomic Dy(III) CAS(9,7) -> RASSI/SOC -> SINGLE_ANISO prototype completed twice;
- exact prototype figures, tables, inputs, logs, provenance, and selected source artifacts.

## Corrected technical result

The repository supports a **model-prioritized Hans pocket candidate for further QM and experimental testing** and an **authentic OpenMolcas execution prototype**. It does not establish thermodynamic binding free energy, Dy selectivity, zero-waste separation, or a project-specific CASSCF/CASPT2 mechanism.

## Quantum evidence boundary

The atomic-Dy prototype is not a Dy-LanM cluster and there is no project-specific Nd-LanM production run. It must not be used to replace or authenticate manuscript Figures 6-7. See `OPENMOLCAS_PROTOTYPE_STATUS.md` and `SOURCE_DATA_GAPS.md`.

## Original technical corrections retained

- 70 design-role memberships reconciled to 66 unique residue positions.
- Sequence funnel separated into generated, unique, shortlisted, scored, and final units.
- ProteinMPNN and LigandMPNN temperatures, seeds, and counts exposed.
- Rosetta relabeled as one-structure score-only triage with no structural-realism acceptance cutoff.
- OpenMM protocol, all 45 replica-level geometry results, and generic OPC3-compatible 12-6-4 parameters reported.
- Hans chain-A pocket and A-D crystallographic-context models distinguished; Arg100 lineage labeled.
- Invalid `Delta G` transformation replaced by a dimensionless Composite Quality Index.

## Provenance

The original correction is anchored at commit `549b22918e1edf409e0c75831cb9af7069468275`; the v2.0 update was prepared from `main` commit `e425926d3b22de50beef9dda92a657767b52df66`. Prototype OpenMolcas source commit: `8355057f32d65706a35996b5ab07cac2962bb728`.
