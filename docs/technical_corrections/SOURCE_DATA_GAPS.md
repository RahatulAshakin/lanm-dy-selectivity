# Source-data gaps and deferred manuscript actions

## Quantum files required before claiming reproducibility

1. Raw Dy and Nd CASSCF input/output files.
2. Raw CASPT2 input/output files with state count, imaginary/IPEA shifts, and intruder-state diagnostics.
3. Exact relativistic Hamiltonian/ECP and basis sets for every element.
4. Spin-orbit calculation files or a tightly scoped justification for their absence.
5. Reduced-cluster coordinates, charge/multiplicity, caps, constraints, and atom/residue inventory.
6. Natural occupation numbers for every active orbital and every reported root.
7. Molden/GBW or equivalent orbital source files for figure regeneration.

No ORCA/CASSCF/CASPT2 files were found in the audited repository tree. Quantum Figures 6-7 are therefore not regenerated in this package.

## Manuscript-only actions intentionally not applied

- Correct the author spelling requested by Reviewer 2.
- Remove or justify references 30-36.
- Narrow title/abstract language concerning zero-waste separation and experimental validation.
- Propagate the Composite Quality Index terminology through manuscript prose and numbering.

## Additional technical extensions

- Expand the off-target panel beyond Al(III)/Fe(III) to include chemically relevant Ca(II), Mg(II), Zn(II), and other process contaminants.
- Add force-field sensitivity analysis, especially for short Al-O and Fe-O distances.
- Replace 100 ps reduced screening with converged sampling and a justified affinity/selectivity observable.
