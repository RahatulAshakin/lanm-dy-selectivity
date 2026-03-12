# Project description

Below is a tighter, execution-ready version of your project. To keep the study testable, I make ion-adsorption-clay-like leachate the primary design environment, and I use xenotime-like Dy/Y/Tb competition as an external validation set.
Gap. Lanmodulin already has a strong structural and biochemical foundation for engineering: AM1 LanM is a reviewed lanthanide-binding protein in UniProt, and experimentally determined structures exist for Y-bound LanM (6MI5), Nd-bound AM1 LanM (8FNS), and La/Dy-bound Hansschlegelia quercus LanM (8DQ2/8FNR). Hans-LanM also shows radius-sensitive dimerization and has already enabled single-stage Nd/Dy separation to >98% individual purities, while a recent computational study reproduced binding trends for wild-type LanM and a small set of EF-hand variants. What is still missing is an ore-conditioned, Dy-specific, end-to-end design pipeline that starts from these structures, explicitly optimizes against real competitor ions, and adds a quantum active-space refinement to the Dy-binding pocket after classical MD/QM screening.
Hypothesis. If the first-shell EF-hand geometry and second-sphere hydrogen-bond/dimer-interface residues of LanM are redesigned to better match the smaller, later-lanthanide coordination preferences that distinguish Dy from Nd, then at least one LanM variant will show a more favorable aqueous binding free energy for Dy than for Nd and Y, while avoiding persistent Al/Fe inner-sphere capture. This is plausible because current LanM studies show that second-sphere residues measurably affect affinity and that picometre-scale lanthanide radius differences can be amplified into large structural and separation effects in Hans-LanM.
One clear objective. Design and down-select one LanM variant that, in an ion-adsorption-clay-inspired multi-metal aqueous model, achieves

at 298 K, which corresponds to at least a 10-fold Dy-over-Nd selectivity, while showing no stable inner-sphere Al or Fe coordination during replicate MD simulations.
Detailed methodology
1) Data assembly and preprocessing
Use your uploaded local bundle as the project backbone: lanmodulin_sequences.csv, 8FNS_atoms.csv, 8DQ2_atoms.csv, lanmodulin_lanthanide_structures.csv, and lanthanide_elements_basic_properties.csv. Add the experimentally solved structures referenced in the manifest—6MI5, 8FNS, 8DQ2, and 8FNR—and standardize all structures to the mature AM1 LanM numbering from UniProt C5B164. Build residue–metal distance matrices from the atom tables, define first-shell atoms as donor atoms within 3.2 Å of the bound metal, and define second-sphere residues as residues 3.2–6.0 Å from the metal or hydrogen-bonded to first-shell residues. Use the ore compositions from your proposal to convert the competitor set into weights,

with and . Use MetalPDB for known coordination-geometry templates and MIB2 only as a fallback sanity check for designed pockets lacking confident local geometry. Your uploaded repository-links note also adds useful public rare-earth datasets for AMBER/PLUMED setup and adsorption validation.
lanthanide_repository_download_…
2) Structure-conditioned generative mutagenesis
Run a constrained multi-state design workflow. First use ProteinMPNN to preserve residues that must remain fixed while sampling allowed second-sphere and interface mutations; ProteinMPNN supports fixed-region design on a target backbone. Then use LigandMPNN, because it explicitly models nonprotein components—including metals—and performs especially well for residues contacting metals. Finally, refine each candidate with Rosetta fixbb and FastDesign so side chains are repacked and the pocket is relaxed under Dy-specific distance and angle restraints. In round 1, keep the canonical first-shell carboxylates fixed; in round 2, allow a restricted mutation set at second-sphere and dimer-interface positions. Generate about 500–1000 sequences, then retain only the top ~50 by Rosetta energy, pocket geometry, and fold preservation.
3) Classical MD and metadynamics prescreen
For each retained variant, simulate complexes with Dy, Nd, Y, Al, and Fe. Use GROMACS with AMBER19SB and OPC/OPC3 water, since current GROMACS guidance recommends OPC or OPC3 with AMBER19SB. For the metal center, use 12-6-4 LJ parameters rather than plain 12-6 LJ, and apply chelator-based tuning where needed because these models improve hydration free energies, ion–oxygen distances, coordination numbers, and LanM EF-loop binding energies relative to standard parameters. The pair potential is

Run 3 replicas per metal–variant pair, with minimization, NVT/NPT equilibration, and production trajectories. Then bias the local pocket with metadynamics using collective variables such as coordination number and average metal–oxygen distance:

Reject any design that shows long-lived Al/Fe inner-sphere binding, large pocket collapse, or Dy coordination inconsistent with the intended 8–9 donor environment.
4) Relativistic QM/QCT thermodynamics
From the best MD snapshots, extract a reduced binding-site cluster containing the metal, first-shell residues, the most important second-sphere H-bond residues, and a small number of explicit waters. Compute scalar-relativistic energies in ORCA 6.0 with DKH or X2C, then add solvation and thermal terms through a quasi-chemical/thermodynamic-cycle treatment. A practical working expression is

The formal QCT hydration term can be written as

Use as the starting hydration model for Dy-like later lanthanides and for Nd-like earlier lanthanides, then confirm with MD/metadynamics instead of forcing those values.
5) Quantum algorithms for the Dy pocket
The quantum part should be applied only to the reduced metal-pocket active space, not to the full 133-aa protein. Build the second-quantized Hamiltonian from ORCA/PySCF one- and two-electron integrals:

Reduce it with Qiskit Nature’s ActiveSpaceTransformer, then map it to qubits with JordanWignerMapper:

Use ADAPT-VQE as the primary executable quantum algorithm:

with operator selection by the largest gradient,

Then apply Quantum Subspace Expansion (QSE) to refine near-degenerate local states and improve robustness:

Use the quantum correction only as a pocket-level refinement,

and add it to the classical thermodynamic estimate,

Do not make QPE the primary project algorithm: IBM’s current documentation explicitly notes that QPE requires fault-tolerant hardware and deep circuits, whereas VQE remains the practical near-term route. In practice, run this stage first on a simulator and only send the smallest active spaces to hardware.
6) Ore-weighted ranking and decision rule
For each variant , calculate

and convert it to a selectivity factor,

Then rank each design with an ore-conditioned score,

Lower scores are better. Down-select the top 3 variants for synthesis only if they pass all three gates: Dy-selectivity target, Dy-pocket structural stability, and no persistent Al/Fe inner-sphere capture.
7) Reverse translation and experimental execution
Reverse-translate the top variants with IDT Codon Optimization or JCat. IDT screens codon usage while lowering sequence complexity and minimizing problematic secondary structures, and JCat optimizes CAI while also allowing removal of unwanted restriction sites and transcription terminators. Add a bicistronic design (BCD) upstream of the LanM ORF to reduce start-region mRNA secondary-structure blockage. Express the top constructs in E. coli, purify them, and measure apparent or competitive binding by EGTA-buffered CD/competition assays in the same spirit as prior LanM variant studies. Then validate the best construct in an immobilized format—column resin or magnetic nanoparticle format—because immobilized LanM systems have already shown selective REE capture, high desorption efficiency, reuse across cycles, and strong enrichment from complex leachates.
How to execute the project in practice
Preprocess the uploaded files into aligned Dy/Nd/La/Y templates and generate first-shell/second-sphere annotations.
Generate metal-aware variants with ProteinMPNN/LigandMPNN and relax them with Rosetta.
Run multi-metal MD/metadynamics and reject unstable or off-target-binding designs.
Compute relativistic DFT/QCT binding terms on extracted clusters.
Apply ADAPT-VQE + QSE only to the reduced Dy-pocket active space for the top candidates.
Rank by ore-weighted selectivity score and select the top 3 variants.
Codon-optimize, build BCD constructs, express, and test in competitive synthetic leachates, then in immobilized separation format.
The key practical point is that the workflow is classical-first and quantum-refined: the full protein is handled by structure design, MD, and relativistic QM, while the quantum computer is used only where it is most defensible—on the reduced Dy pocket Hamiltonian. That makes the project executable now, instead of waiting for fault-tolerant quantum hardware
