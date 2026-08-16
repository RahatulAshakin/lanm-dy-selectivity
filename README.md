# lanm-dy-selectivity

## Reviewer-driven technical correction notice

This repository now separates **audited technical corrections** from legacy reporting claims. The corrected package is in `docs/technical_corrections/` and `results/technical_corrections/`.

The technically supportable result is a **model-prioritized Hans pocket candidate for further quantum and experimental testing**. The current repository does **not** establish a thermodynamic binding free energy, Dy selectivity, zero-waste separation, or a reproducible CASSCF/CASPT2 mechanism.

Important corrections:

- The previous `deltaG_pocket_proxy_kcal_mol` transformation is scientifically invalid and is not used in the correction package. It is replaced by a dimensionless Composite Quality Index.
- Rosetta was a deterministic, one-structure `score_jd2` triage with no relaxation or structural-realism cutoff; scores are compared only within topology class.
- Round-2 OpenMM used 100 ps per replica, three seeds per condition, amber19/OPC3, and a generic OPC3-compatible 12-6-4 metal model. These geometry screens are not binding free energies.
- Raw ORCA/CASSCF/CASPT2 source files are absent, so quantum figures and claims are not repository-reproducible.
- The A-D chains in 8FNR/8DQ2 are crystallographic asymmetric-unit content; RCSB reports author/PISA A2 biological assemblies.

## Build the corrected technical package

```bash
python -m pip install -r requirements-technical-corrections.txt
python scripts/build_technical_corrections.py --output-root results/technical_corrections --report-path docs/technical_corrections/LanM_Technical_Corrections_Report.docx
```

See `docs/technical_corrections/README.md` for the correction matrix and `docs/technical_corrections/SOURCE_DATA_GAPS.md` for outstanding author actions.

## Legacy workflow note

The legacy `python -m lanm.cli.build_final_project_report` path is retained for audit history but must not be used for reviewer-response submission because it emits the invalid energy-like proxy and references outputs absent from the audited `main` branch.
