# LanM technical corrections

This directory contains the reviewer-driven technical correction package for manuscript `JMGM-D-26-01758`. The main manuscript is intentionally unchanged.

## Corrected technical result

The repository supports a **model-prioritized Hans pocket candidate for further QM and experimental testing**. It does not currently establish thermodynamic binding free energy, Dy selectivity, zero-waste separation, or a validated CASSCF/CASPT2 mechanism.

## Corrections included

- 70 design-role memberships reconciled to 66 unique residue positions.
- Sequence funnel separated into generated, unique, shortlisted, scored, and final units.
- ProteinMPNN and LigandMPNN temperatures, seeds, and counts exposed.
- Rosetta relabeled as one-structure score-only triage with no structural-realism acceptance cutoff.
- OpenMM protocol, all 45 replica-level geometry results, and generic OPC3-compatible 12-6-4 parameters reported.
- Hans chain-A pocket and A-D crystallographic-context models distinguished; Arg100 lineage labeled.
- Invalid `Delta G` transformation replaced by a dimensionless Composite Quality Index.
- Missing ORCA/CASSCF/CASPT2 source data recorded as a submission blocker rather than reconstructed.

## Rebuild

```bash
python -m pip install -r requirements-technical-corrections.txt
python scripts/build_technical_corrections.py --output-root results/technical_corrections --report-path docs/technical_corrections/LanM_Technical_Corrections_Report.docx
```

The build is reporting-only. It does not rerun ProteinMPNN, LigandMPNN, Rosetta, OpenMM, or quantum chemistry.

## Primary provenance

- `results/tables/design_mask_candidates.csv`
- `results/proteinmpnn_smoke/**/run_command.txt`
- `results/proteinmpnn_round2_smoke/**/run_command.txt`
- `results/ligandmpnn_smoke/**/run_command.txt`
- `results/ligandmpnn_round2_smoke/**/run_command.txt`
- `results/tables/rosetta_round2_candidate_ranking.csv`
- `results/round2_openmm_screening/**/simulation_config.yaml`
- `results/round2_openmm_screening/**/screening_log.txt`
- `config/metal_parameter_values.yaml`

Structure scope is cross-checked against [RCSB 8FNR](https://www.rcsb.org/structure/8FNR) and [RCSB 8DQ2](https://www.rcsb.org/structure/8DQ2).
