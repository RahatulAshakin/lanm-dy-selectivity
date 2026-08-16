#!/usr/bin/env python3
"""Build the LanM technical-correction submission package.

This script is deliberately reporting-only. It regenerates tables, figures,
the Word/PDF report, repository-ready documentation, and a ZIP bundle from the
audited values recorded in the repository artifacts named in each table.
It does not modify the manuscript and it does not manufacture absent QM data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import statistics
import subprocess
import sys
import textwrap
from collections import defaultdict
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


BASE_COMMIT = "549b22918e1edf409e0c75831cb9af7069468275"
REPO_URL = "https://github.com/RahatulAshakin/lanm-dy-selectivity"
RCSB_8FNR = "https://www.rcsb.org/structure/8FNR"
RCSB_8DQ2 = "https://www.rcsb.org/structure/8DQ2"

COLORS = {
    "navy": "16324F",
    "blue": "2E74B5",
    "dark_blue": "1F4D78",
    "teal": "287271",
    "cyan": "72B7B2",
    "gold": "B07C17",
    "orange": "E07A3F",
    "red": "9B1C1C",
    "green": "2F6B3B",
    "gray": "5E6670",
    "light_gray": "F2F4F7",
    "blue_gray": "E8EEF5",
    "pale_gold": "FFF6DD",
    "pale_red": "FBEAEA",
    "white": "FFFFFF",
    "black": "111111",
}

FIGURE_DPI = 240


MASK_ROWS = [
    {"role": "First-shell", "membership_count": 26, "unique_inventory_count": 66},
    {"role": "Second-sphere", "membership_count": 16, "unique_inventory_count": 66},
    {"role": "Interface", "membership_count": 10, "unique_inventory_count": 66},
    {"role": "Protected", "membership_count": 18, "unique_inventory_count": 66},
]

OVERLAP_ROWS = [
    {"canonical_position": 17, "am1_position": "", "residue": "ASN", "roles": "first-shell; protected", "membership_excess": 1},
    {"canonical_position": 18, "am1_position": 16, "residue": "LYS", "roles": "second-sphere; interface", "membership_excess": 1},
    {"canonical_position": 44, "am1_position": 40, "residue": "GLY/LYS", "roles": "second-sphere; interface", "membership_excess": 1},
    {"canonical_position": 83, "am1_position": 79, "residue": "ARG/GLU", "roles": "second-sphere; interface", "membership_excess": 1},
]

FUNNEL_ROWS = [
    {"phase": "Mask inventory", "method": "Residue-role map", "raw_generated": "", "unique_sequences": "", "retained": 66, "unit": "unique positions", "evidence": "results/tables/design_mask_candidates.csv"},
    {"phase": "Round 1", "method": "ProteinMPNN", "raw_generated": 160, "unique_sequences": 44, "retained": 12, "unit": "sequences", "evidence": "results/proteinmpnn_smoke/** and results/tables/proteinmpnn_shortlist.csv"},
    {"phase": "Round 1", "method": "LigandMPNN", "raw_generated": 24, "unique_sequences": 16, "retained": 8, "unit": "sequences", "evidence": "results/ligandmpnn_smoke/** and results/tables/ligandmpnn_shortlist.csv"},
    {"phase": "Round 1", "method": "Integrated seeds", "raw_generated": "", "unique_sequences": "", "retained": 3, "unit": "candidates", "evidence": "results/tables/integrated_candidate_ranking.csv"},
    {"phase": "Round 2", "method": "ProteinMPNN", "raw_generated": 120, "unique_sequences": 74, "retained": 6, "unit": "sequences", "evidence": "results/proteinmpnn_round2_smoke/** and results/tables/proteinmpnn_round2_shortlist.csv"},
    {"phase": "Round 2", "method": "LigandMPNN", "raw_generated": 12, "unique_sequences": 6, "retained": 4, "unit": "sequences", "evidence": "results/ligandmpnn_round2_smoke/** and results/tables/ligandmpnn_round2_shortlist.csv"},
    {"phase": "Round 2", "method": "Rosetta score-only triage", "raw_generated": 4, "unique_sequences": 4, "retained": 4, "unit": "scored candidates", "evidence": "results/tables/rosetta_round2_candidate_ranking.csv"},
    {"phase": "Round 2", "method": "OpenMM reduced panel", "raw_generated": "", "unique_sequences": "", "retained": 3, "unit": "finalists", "evidence": "results/round2_openmm_screening/**"},
]

LINEAGE_ROWS = [
    {
        "finalist": "hans_pocket_ss_only_u02_r2u01",
        "model_scope": "Hans chain-A pocket model",
        "seed_mutations_source_notation": "K12T; R77K; K87G; T92D",
        "seed_mutations_unambiguous": "Chain A: K35T; R100K; K110G; T115D",
        "round2_additional_mutations": "Chain A: A44D; T55W; M92L; K97A",
        "arg100_status": "Template Arg100 becomes Lys100 in seed",
    },
    {
        "finalist": "hans_interface_ss_plus_if_u04_r2u01",
        "model_scope": "Hans chain A designed; chains B-D retained as fixed crystallographic context",
        "seed_mutations_source_notation": "K12S; T32W; K87G; T92E",
        "seed_mutations_unambiguous": "Chain A: K35S; T55W; K110G; T115E",
        "round2_additional_mutations": "None in selected LigandMPNN design",
        "arg100_status": "Template Arg100 retained",
    },
    {
        "finalist": "am1_mex_ss_only_u02_r2u02",
        "model_scope": "AM1/Mex chain-A monomer model",
        "seed_mutations_source_notation": "K10D; L55Y; N59K; R90A",
        "seed_mutations_unambiguous": "Chain A: K38D; L83Y; N87K; R118A",
        "round2_additional_mutations": "Chain A: K57A; K93R; K94A",
        "arg100_status": "Not the Hans Arg100 comparison",
    },
]

OPENMM_REPLICATES = {
    "hans_pocket_ss_only_u02_r2u01": {
        "Dy": [2.214204, 2.215379, 2.215573],
        "Nd": [2.316352, 2.316589, 2.315862],
        "Y": [2.194464, 2.196893, 2.186353],
        "Al": [1.726393, 1.723088, 1.725170],
        "Fe": [1.926443, 1.933505, 1.923926],
    },
    "hans_interface_ss_plus_if_u04_r2u01": {
        "Dy": [2.207488, 2.203839, 2.203308],
        "Nd": [2.310945, 2.308004, 2.310995],
        "Y": [2.189468, 2.190009, 2.191884],
        "Al": [1.720404, 1.720473, 1.721261],
        "Fe": [1.930395, 1.931583, 1.931271],
    },
    "am1_mex_ss_only_u02_r2u02": {
        "Dy": [2.201751, 2.202536, 2.196603],
        "Nd": [2.300359, 2.301707, 2.304425],
        "Y": [2.188296, 2.183940, 2.185430],
        "Al": [1.714569, 1.715489, 1.714316],
        "Fe": [1.930376, 1.926456, 1.923507],
    },
}

CANDIDATE_LABELS = {
    "hans_pocket_ss_only_u02_r2u01": "Hans pocket",
    "hans_interface_ss_plus_if_u04_r2u01": "Hans interface",
    "am1_mex_ss_only_u02_r2u02": "AM1/Mex",
}

MODEL_SCORES = [
    {"candidate_id": "hans_pocket_ss_only_u02_r2u01", "ligand_confidence": 1.0000, "overall_confidence": 0.3399},
    {"candidate_id": "hans_interface_ss_plus_if_u04_r2u01", "ligand_confidence": 0.3454, "overall_confidence": 0.3303},
    {"candidate_id": "am1_mex_ss_only_u02_r2u02", "ligand_confidence": 0.2213, "overall_confidence": 0.3023},
]

METAL_PARAMETERS = [
    {"metal": "Dy", "formal_charge": 3, "rmin_half_A": 1.632, "epsilon_kcal_mol": 0.0962022, "c4_kcal_mol_A4": 183},
    {"metal": "Nd", "formal_charge": 3, "rmin_half_A": 1.712, "epsilon_kcal_mol": 0.1464093, "c4_kcal_mol_A4": 184},
    {"metal": "Y", "formal_charge": 3, "rmin_half_A": 1.626, "epsilon_kcal_mol": 0.09289608, "c4_kcal_mol_A4": 192},
    {"metal": "Al", "formal_charge": 3, "rmin_half_A": 1.361, "epsilon_kcal_mol": 0.01031847, "c4_kcal_mol_A4": 363},
    {"metal": "Fe", "formal_charge": 3, "rmin_half_A": 1.455, "epsilon_kcal_mol": 0.02662782, "c4_kcal_mol_A4": 429},
]

MODEL_SPEC_ROWS = [
    {"component": "ProteinMPNN round 1", "exact_setting": "20 sequences/target at temperatures 0.10 and 0.15; seed 37; batch size 1", "interpretation_limit": "Generative sequence sampling; not deterministic in the mathematical sense despite fixed seed."},
    {"component": "ProteinMPNN round 2", "exact_setting": "20 sequences/target at temperatures 0.10 and 0.15; seed 37; batch size 1", "interpretation_limit": "Same fixed-seed sampling settings."},
    {"component": "LigandMPNN rounds 1/2", "exact_setting": "Temperature 0.10; seed 37; 2 designs per candidate; one batch; side-chain and ligand context enabled", "interpretation_limit": "Model confidence is not a probability of thermodynamic binding."},
    {"component": "Rosetta round 2", "exact_setting": "score_jd2; nstruct=1; seed 37; fa_elec weight=0; waters retained; auto metal setup; rank within topology", "interpretation_limit": "Score-only triage. No relaxation, no structural-realism acceptance cutoff, and no cross-topology score comparison."},
    {"component": "OpenMM round 2", "exact_setting": "298 K; 1 atm; 2 fs; friction 1 ps^-1; 100 ps production; 3 replicas, seeds 101-103; 500-step reporting", "interpretation_limit": "Short reduced rescreen; not a binding free-energy calculation."},
    {"component": "OpenMM metric", "exact_setting": "Mean over frames of the minimum distance from each retained metal atom to any oxygen; 3.2 A inner-sphere cutoff", "interpretation_limit": "Occupancy=1.0 at this cutoff does not establish selectivity or experimentally persistent binding."},
    {"component": "Force field", "exact_setting": "amber19-all.xml + amber19/opc3.xml; generic published OPC3-compatible 12-6-4 metal family", "interpretation_limit": "LanM-chelator-tuned numeric parameters are absent; Al/Fe short distances require sensitivity analysis."},
]


def reviewer_rows() -> list[dict[str, str]]:
    return [
        {"reviewer_item": "R1.1 / R2-B", "issue": "Undefined and thermodynamically invalid pocket delta-G proxy", "technical_correction": "Replace with dimensionless Composite Quality Index = 0.7 ligand confidence + 0.3 overall confidence; remove -RT ln transformation.", "status": "Corrected", "deliverable": "Figure TC5; model_quality_index.csv"},
        {"reviewer_item": "R1.2", "issue": "Natural occupations and active-space suitability", "technical_correction": "Do not regenerate without raw multiroot orbital and occupation data. Require all roots, orbitals, occupations, and active-space justification.", "status": "Blocked: source data absent", "deliverable": "quantum_evidence_gap_register.csv"},
        {"reviewer_item": "R1.3", "issue": "CASPT2 shift, intruder states, number of states", "technical_correction": "Require program/version, state averaging, imaginary/IPEA shifts, state counts, convergence, and intruder-state diagnostics.", "status": "Blocked: source data absent", "deliverable": "quantum_evidence_gap_register.csv"},
        {"reviewer_item": "R1.4", "issue": "Reduced-cluster size and environmental effect", "technical_correction": "Require atom/residue inventory, total charge/multiplicity, caps/constraints, coordinates, and sensitivity comparison.", "status": "Blocked: source data absent", "deliverable": "quantum_evidence_gap_register.csv"},
        {"reviewer_item": "R1.5", "issue": "Zero-waste claim lacks experimental evidence", "technical_correction": "Technical package states that no experimental separation or waste-balance evidence is present and narrows the result to a computational candidate-prioritization claim.", "status": "Flagged; manuscript wording outside scope", "deliverable": "Technical report Sections 1 and 6"},
        {"reviewer_item": "R1.6", "issue": "Only Al and Fe competitors", "technical_correction": "State that Al(III)/Fe(III) were the tested off-target subset. Do not generalize to Ca/Mg/Zn; recommend an expanded panel.", "status": "Corrected limitation", "deliverable": "modeling_specification.csv"},
        {"reviewer_item": "R1.7", "issue": "MPNN parameters and filtering unspecified", "technical_correction": "Expose exact temperatures, seeds, generated/unique/retained counts, and selection units.", "status": "Corrected", "deliverable": "Figure TC2; sequence_funnel.csv; modeling_specification.csv"},
        {"reviewer_item": "R1.8 / R2-S4", "issue": "Rosetta realism criteria/cutoffs unspecified", "technical_correction": "Correct description to score-only triage: score_jd2, nstruct=1, fixed seed, no relax, no acceptance cutoff, within-topology ranking only.", "status": "Corrected and qualified", "deliverable": "modeling_specification.csv"},
        {"reviewer_item": "R1.9", "issue": "Hans pocket versus interface unclear", "technical_correction": "Distinguish chain-A pocket model from chain-A-designed/B-D-fixed crystallographic context; attribute score difference to model confidence, not inferred structural causality.", "status": "Corrected", "deliverable": "Figure TC3; finalist_lineage.csv"},
        {"reviewer_item": "R2-A1", "issue": "Spin-orbit coupling not documented", "technical_correction": "Require SOC/RASSI-SO/NEVPT2-SO treatment or an explicit scope-limited justification. No claim made from absent data.", "status": "Blocked: source data absent", "deliverable": "quantum_evidence_gap_register.csv"},
        {"reviewer_item": "R2-A2", "issue": "Relativistic Hamiltonian and basis sets absent", "technical_correction": "Require exact Hamiltonian/ECP, basis sets for all atoms, auxiliary bases, scalar-relativistic settings, and software version.", "status": "Blocked: source data absent", "deliverable": "quantum_evidence_gap_register.csv"},
        {"reviewer_item": "R2-A3", "issue": "Active orbitals/occupations across roots not demonstrated", "technical_correction": "Require orbital plots and natural occupation numbers for every reported root; verify seven 4f orbitals remain active.", "status": "Blocked: source data absent", "deliverable": "quantum_evidence_gap_register.csv"},
        {"reviewer_item": "R2-C1", "issue": "Metal force-field provenance and Al/Fe artifacts", "technical_correction": "Report actual generic OPC3-compatible 12-6-4 parameters; state tuned chelator values are absent; treat Al/Fe distances as model-dependent geometry, not affinity evidence.", "status": "Corrected with caution", "deliverable": "metal_parameter_values.csv; Figure TC4"},
        {"reviewer_item": "R2-Minor 1", "issue": "Christopher Stuetzle spelling", "technical_correction": "Record the requested spelling correction in the response matrix only; manuscript author line is intentionally unchanged in this package.", "status": "Deferred: manuscript outside scope", "deliverable": "reviewer_response_matrix.csv"},
        {"reviewer_item": "R2-Minor 2", "issue": "Arg100 not labeled", "technical_correction": "Identify chain-A Arg100 explicitly and show when it becomes Lys100 in the Hans pocket seed.", "status": "Corrected", "deliverable": "Figure TC3; finalist_lineage.csv"},
        {"reviewer_item": "R2-Minor 3", "issue": "References 30-36 irrelevant", "technical_correction": "Record as a manuscript bibliography action; no bibliography is modified in this technical-only package.", "status": "Deferred: manuscript outside scope", "deliverable": "reviewer_response_matrix.csv"},
        {"reviewer_item": "R2-S4 Y-distance", "issue": "Y distances converge near 2.19 Å", "technical_correction": "Report all nine replica values. Similarity is due to the common Y(III) model and near-identical local coordination; displayed equality was rounding, not duplicated data.", "status": "Corrected", "deliverable": "Figure TC4; openmm_replicates.csv"},
        {"reviewer_item": "R2-Data", "issue": "Deposit raw ORCA .in/.out", "technical_correction": "Repository audit finds no ORCA/CASSCF/CASPT2 inputs or outputs. Package lists the required filenames/metadata and marks Figs. 6-7 non-regenerable.", "status": "Blocked: author files required", "deliverable": "quantum_evidence_gap_register.csv; SOURCE_DATA_GAPS.md"},
    ]


QUANTUM_GAPS = [
    {"required_artifact": "Dy CASSCF input and full output", "minimum_content": "Coordinates, charge/multiplicity, Hamiltonian/ECP, basis sets, CAS(9,7), roots, convergence, occupations", "repository_status": "Absent", "affected_claim_or_figure": "Dy active-space claims; Figures 6-7"},
    {"required_artifact": "Nd CASSCF input and full output", "minimum_content": "Coordinates, charge/multiplicity, Hamiltonian/ECP, basis sets, CAS(3,7), roots, convergence, occupations", "repository_status": "Absent", "affected_claim_or_figure": "Nd active-space claims; Figures 6-7"},
    {"required_artifact": "Dy/Nd CASPT2 inputs and outputs", "minimum_content": "State count, imaginary and IPEA shifts, intruder-state diagnostics, energies and convergence", "repository_status": "Absent", "affected_claim_or_figure": "CASPT2 comparison"},
    {"required_artifact": "SOC calculation files", "minimum_content": "SOC method, states mixed, spin-free roots, matrices/energies, convergence", "repository_status": "Absent", "affected_claim_or_figure": "Electronic-state ordering and anisotropy claims"},
    {"required_artifact": "Reduced-cluster coordinate and definition files", "minimum_content": "Atom/residue inventory, caps, constraints, charge/multiplicity, retained waters, cluster size", "repository_status": "Absent", "affected_claim_or_figure": "Environmental-scope claims"},
    {"required_artifact": "Natural occupation tables for every root", "minimum_content": "Seven active-orbital occupations per root with orbital identities", "repository_status": "Absent", "affected_claim_or_figure": "Active-space appropriateness"},
    {"required_artifact": "Orbital image source data", "minimum_content": "Molden/gbw or equivalent files plus orbital indices and plotting settings", "repository_status": "Absent", "affected_claim_or_figure": "Figure 7"},
]

STRUCTURAL_SCOPE = [
    {"pdb_id": "8FNR", "deposited_chains": "A,B,C,D", "deposited_metal_inventory": "Dy: A=4, B=3, C=3, D=4 (14 total)", "author_biological_assembly": "A2 homodimer (C2); two author/PISA assemblies are listed", "technical_use": "Pocket model uses chain A; interface model keeps A-D crystallographic context", "source": RCSB_8FNR},
    {"pdb_id": "8DQ2", "deposited_chains": "A,B,C,D", "deposited_metal_inventory": "La: 3 per chain (12 total)", "author_biological_assembly": "A2 homodimer (C2); author/PISA; gel-filtration evidence for assembly 1", "technical_use": "Reference only; do not call A-D the biological tetramer", "source": RCSB_8DQ2},
    {"pdb_id": "8FNS", "deposited_chains": "A", "deposited_metal_inventory": "Nd: 4", "author_biological_assembly": "Monomeric reference in repository workflow", "technical_use": "AM1/Mex chain-A model", "source": "results/tables/template_chain_summary.csv"},
    {"pdb_id": "6MI5", "deposited_chains": "X", "deposited_metal_inventory": "Y: 3", "author_biological_assembly": "Monomeric reference in repository workflow", "technical_use": "AM1/Mex reference", "source": "results/tables/template_chain_summary.csv"},
]


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def openmm_replica_rows() -> list[dict]:
    rows = []
    for candidate_id, metals in OPENMM_REPLICATES.items():
        for metal, values in metals.items():
            for replica_id, value in enumerate(values, start=1):
                rows.append({
                    "candidate_id": candidate_id,
                    "candidate_label": CANDIDATE_LABELS[candidate_id],
                    "metal": metal,
                    "replica_id": replica_id,
                    "seed": 100 + replica_id,
                    "production_duration_ps": 100.0,
                    "screened_frame_count": 100,
                    "metal_atom_count": 4,
                    "mean_min_metal_oxygen_distance_A": f"{value:.6f}",
                    "inner_sphere_cutoff_A": 3.2,
                    "inner_sphere_occupancy_fraction": "1.000000",
                    "source_path": f"results/round2_openmm_screening/{candidate_id}/{metal}/replicate_{replica_id}/screening_log.txt",
                })
    return rows


def openmm_summary_rows() -> list[dict]:
    rows = []
    for candidate_id, metals in OPENMM_REPLICATES.items():
        for metal, values in metals.items():
            rows.append({
                "candidate_id": candidate_id,
                "candidate_label": CANDIDATE_LABELS[candidate_id],
                "metal": metal,
                "replicate_count": len(values),
                "mean_distance_A": f"{statistics.mean(values):.6f}",
                "sample_sd_A": f"{statistics.stdev(values):.6f}",
                "minimum_replica_mean_A": f"{min(values):.6f}",
                "maximum_replica_mean_A": f"{max(values):.6f}",
                "interpretation": "Model-dependent geometry screening; not binding free energy or selectivity",
            })
    return rows


def model_score_rows() -> list[dict]:
    rows = []
    for row in MODEL_SCORES:
        cqi = 0.7 * row["ligand_confidence"] + 0.3 * row["overall_confidence"]
        rows.append({
            "candidate_id": row["candidate_id"],
            "candidate_label": CANDIDATE_LABELS[row["candidate_id"]],
            "ligand_confidence": f"{row['ligand_confidence']:.4f}",
            "overall_confidence": f"{row['overall_confidence']:.4f}",
            "composite_quality_index": f"{cqi:.6f}",
            "rank": 0,
            "interpretation": "Dimensionless model-based prioritization; not thermodynamic free energy",
        })
    for rank, row in enumerate(sorted(rows, key=lambda x: float(x["composite_quality_index"]), reverse=True), start=1):
        row["rank"] = rank
    return sorted(rows, key=lambda x: x["rank"])


def setup_plot_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "axes.labelsize": 10.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def fig_mask(path: Path) -> None:
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.2), gridspec_kw={"width_ratios": [1.25, 1]})
    roles = [r["role"] for r in MASK_ROWS]
    counts = [r["membership_count"] for r in MASK_ROWS]
    colors = ["#287271", "#4C78A8", "#E07A3F", "#9A9EAB"]
    y = np.arange(len(roles))
    axes[0].barh(y, counts, color=colors, height=0.6)
    axes[0].set_yticks(y, labels=roles)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Role memberships")
    axes[0].set_title("A. Role membership counts")
    axes[0].set_xlim(0, 30)
    axes[0].grid(axis="x", alpha=0.2)
    for yi, count in zip(y, counts):
        axes[0].text(count + 0.5, yi, str(count), va="center", fontweight="bold")

    axes[1].axis("off")
    axes[1].set_title("B. Reconciliation", loc="left")
    blocks = [
        (0.04, 0.72, 0.40, 0.17, "70", "category memberships", "#E8EEF5"),
        (0.56, 0.72, 0.40, 0.17, "66", "unique positions", "#E7F1EC"),
        (0.30, 0.47, 0.40, 0.16, "4", "overlap excess", "#FFF6DD"),
    ]
    for x, y0, w, h, value, label, fill in blocks:
        axes[1].add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0.012,rounding_size=0.015", transform=axes[1].transAxes, facecolor=fill, edgecolor="#AEB6BF"))
        axes[1].text(x + w / 2, y0 + h * 0.62, value, ha="center", va="center", transform=axes[1].transAxes, fontsize=24, fontweight="bold", color="#16324F")
        axes[1].text(x + w / 2, y0 + h * 0.25, label, ha="center", va="center", transform=axes[1].transAxes, fontsize=9.2)
    axes[1].text(0.04, 0.37, "Overlapping canonical positions", transform=axes[1].transAxes, fontweight="bold", color="#1F4D78")
    overlap_lines = ["17  first-shell + protected", "18  second-sphere + interface", "44  second-sphere + interface", "83  second-sphere + interface"]
    for idx, line in enumerate(overlap_lines):
        axes[1].text(0.07, 0.30 - idx * 0.07, line, transform=axes[1].transAxes, fontsize=9.5)
    fig.suptitle("Corrected design-mask accounting", x=0.06, ha="left", fontsize=16, fontweight="bold", color="#16324F")
    fig.text(0.06, 0.015, "Source: results/tables/design_mask_candidates.csv. 70 is a membership sum; 66 is the unique residue inventory.", fontsize=8.7, color="#5E6670")
    fig.tight_layout(rect=[0.04, 0.05, 0.98, 0.92])
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def fig_funnel(path: Path) -> None:
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), sharey=False)
    panels = [
        (axes[0], "Round 1", [("ProteinMPNN", 160, 44, 12), ("LigandMPNN", 24, 16, 8)]),
        (axes[1], "Round 2", [("ProteinMPNN", 120, 74, 6), ("LigandMPNN", 12, 6, 4)]),
    ]
    colors = ["#4C78A8", "#72B7B2", "#E07A3F"]
    for ax, title, data in panels:
        x = np.arange(len(data))
        width = 0.22
        for j, label in enumerate(["Generated", "Unique", "Retained"]):
            vals = [row[j + 1] for row in data]
            bars = ax.bar(x + (j - 1) * width, vals, width=width, color=colors[j], label=label)
            ax.bar_label(bars, padding=2, fontsize=9, fontweight="bold")
        ax.set_xticks(x, [row[0] for row in data])
        ax.set_title(title)
        ax.set_ylabel("Sequence count")
        ax.grid(axis="y", alpha=0.18)
        ax.legend(frameon=False, ncol=3, loc="upper right")
    fig.suptitle("Corrected sequence-design funnel", x=0.06, ha="left", fontsize=16, fontweight="bold", color="#16324F")
    fig.text(0.06, 0.015, "Positions, raw generations, unique sequences, and retained candidates are reported as separate units.", fontsize=8.8, color="#5E6670")
    fig.tight_layout(rect=[0.04, 0.05, 0.98, 0.92])
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def add_box(ax, xy, width, height, title, body, fill="#F2F4F7", edge="#9AA6B2", title_color="#16324F"):
    box = FancyBboxPatch(xy, width, height, boxstyle="round,pad=0.015,rounding_size=0.018", facecolor=fill, edgecolor=edge, linewidth=1.2, transform=ax.transAxes)
    ax.add_patch(box)
    x, y = xy
    wrap_width = 52 if width >= 0.8 else 24
    wrapped_title = textwrap.fill(title, width=max(18, wrap_width))
    wrapped_body = "\n".join(textwrap.fill(line, width=wrap_width) for line in body.splitlines())
    ax.text(x + 0.02, y + height - 0.04, wrapped_title, transform=ax.transAxes, va="top", ha="left", fontsize=10.0, fontweight="bold", color=title_color)
    ax.text(x + 0.02, y + height - 0.12, wrapped_body, transform=ax.transAxes, va="top", ha="left", fontsize=8.1, color="#252A31", linespacing=1.20)


def fig_scope_arg100(path: Path) -> None:
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 7.2))
    ax = axes[0]
    ax.axis("off")
    ax.set_title("A. Structure and assembly scope", loc="left")
    add_box(ax, (0.06, 0.72), 0.88, 0.19, "8FNR / 8DQ2 deposited model", "Chains A-D are in the crystallographic asymmetric unit.\nA-D is not a biological tetramer.", fill="#E8EEF5")
    add_box(ax, (0.06, 0.43), 0.88, 0.19, "Author/PISA biological assembly", "C2 homodimer (A2). RCSB lists two A2 assemblies for each entry.", fill="#E7F1EC")
    add_box(ax, (0.06, 0.08), 0.40, 0.23, "Hans pocket model", "Chain A only\nPocket-focused design\nArg100 -> Lys100", fill="#FFF6DD", edge="#B07C17")
    add_box(ax, (0.54, 0.08), 0.40, 0.23, "Hans interface model", "Chain A designed\nB-D fixed context\nArg100 retained", fill="#FBEAEA", edge="#9B1C1C")
    for y1, y2 in [(0.71, 0.63), (0.42, 0.31)]:
        ax.add_patch(FancyArrowPatch((0.5, y1), (0.5, y2), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=12, color="#5E6670"))

    ax = axes[1]
    ax.axis("off")
    ax.set_title("B. Unambiguous chain-A mutation lineage", loc="left")
    add_box(ax, (0.04, 0.64), 0.92, 0.25, "Hans pocket finalist", "Seed: K35T; R100K; K110G; T115D\nRound 2: A44D; T55W; M92L; K97A\nArg100 is tracked as R100K.", fill="#FFF6DD", edge="#B07C17")
    add_box(ax, (0.04, 0.34), 0.92, 0.22, "Hans interface finalist", "Seed: K35S; T55W; K110G; T115E\nRound 2: no additional mutation\nTemplate Arg100 is retained.", fill="#FBEAEA", edge="#9B1C1C")
    add_box(ax, (0.04, 0.06), 0.92, 0.20, "Interpretation correction", "The higher pocket score is driven by LigandMPNN confidence (1.000 vs 0.3454), not demonstrated binding energy.", fill="#E8EEF5", edge="#4C78A8")
    fig.suptitle("Corrected Hans model scope and Arg100 labeling", x=0.055, ha="left", fontsize=16, fontweight="bold", color="#16324F")
    fig.text(0.055, 0.015, "This technical schematic clarifies model scope and lineage; it is not a substitute for a validated 3D contact-distance analysis.", fontsize=8.8, color="#5E6670")
    fig.tight_layout(rect=[0.035, 0.05, 0.98, 0.92])
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def fig_openmm(path: Path) -> None:
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(11.8, 6.2))
    metals = ["Dy", "Nd", "Y", "Al", "Fe"]
    candidates = list(OPENMM_REPLICATES)
    offsets = [-0.22, 0, 0.22]
    colors = ["#4C78A8", "#E07A3F", "#287271"]
    markers = ["o", "s", "D"]
    for ci, candidate in enumerate(candidates):
        means = []
        for mi, metal in enumerate(metals):
            vals = OPENMM_REPLICATES[candidate][metal]
            x = mi + offsets[ci]
            jitter = [-0.035, 0, 0.035]
            ax.scatter([x + j for j in jitter], vals, color=colors[ci], marker=markers[ci], s=35, alpha=0.85, zorder=3)
            means.append(statistics.mean(vals))
        ax.plot(np.arange(len(metals)) + offsets[ci], means, color=colors[ci], marker=markers[ci], linewidth=1.7, markersize=6, label=CANDIDATE_LABELS[candidate])
    ax.set_xticks(np.arange(len(metals)), [f"{m}$^{{3+}}$" for m in metals])
    ax.set_ylabel("Mean minimum metal-oxygen distance (Å)")
    ax.set_ylim(1.66, 2.36)
    ax.grid(axis="y", alpha=0.2)
    ax.axvspan(2.5, 4.5, color="#FBEAEA", alpha=0.45)
    ax.text(3.5, 2.345, "Off-target subset; short distances are\nforce-field-dependent geometry, not affinity", ha="center", va="top", fontsize=8.8, color="#9B1C1C")
    ax.text(2.0, 2.11, "Y values differ by replica;\n2.19 Å is rounded display", ha="center", fontsize=8.8, color="#1F4D78")
    ax.legend(frameon=False, ncol=3, loc="lower left")
    ax.set_title("Corrected round-2 OpenMM replica-level geometry screen", loc="left", fontsize=16, color="#16324F")
    fig.text(0.075, 0.015, "Three replicas per candidate/metal (seeds 101-103), 100 ps production each. Metric uses four retained metals and all oxygen atoms; 3.2 Å occupancy was 1.0 for every condition.", fontsize=8.6, color="#5E6670")
    fig.tight_layout(rect=[0.05, 0.055, 0.98, 0.95])
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def fig_cqi(path: Path) -> None:
    setup_plot_style()
    rows = model_score_rows()
    labels = [r["candidate_label"] for r in rows][::-1]
    ligand = [0.7 * float(r["ligand_confidence"]) for r in rows][::-1]
    overall = [0.3 * float(r["overall_confidence"]) for r in rows][::-1]
    totals = [float(r["composite_quality_index"]) for r in rows][::-1]
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(10.8, 4.8))
    ax.barh(y, ligand, color="#4C78A8", label="0.7 x ligand confidence")
    ax.barh(y, overall, left=ligand, color="#72B7B2", label="0.3 x overall confidence")
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 0.88)
    ax.set_xlabel("Composite Quality Index (dimensionless; higher = model priority)")
    ax.grid(axis="x", alpha=0.2)
    for yi, total in zip(y, totals):
        ax.text(total + 0.015, yi, f"{total:.3f}", va="center", fontweight="bold")
    ax.legend(frameon=False, ncol=2, loc="lower right")
    ax.set_title("Corrected model-quality prioritization", loc="left", fontsize=16, color="#16324F")
    fig.text(0.055, 0.025, "No RT factor and no -ln transform are used. This score is not Delta G, binding affinity, or selectivity.", fontsize=9.2, color="#9B1C1C")
    fig.tight_layout(rect=[0.05, 0.11, 0.98, 0.96])
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_cell_margins(cell, top=80, start=120, bottom=80, end=120) -> None:
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for m, v in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{m}"))
        if node is None:
            node = OxmlElement(f"w:{m}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(v))
        node.set(qn("w:type"), "dxa")


def set_table_geometry(table, widths_in: list[float]) -> None:
    if abs(sum(widths_in) - 6.5) > 1e-6:
        raise ValueError(f"Table widths must sum to 6.5 in, got {sum(widths_in)}")
    table.alignment = WD_TABLE_ALIGNMENT.LEFT
    table.autofit = False
    tbl_pr = table._tbl.tblPr
    tbl_w = tbl_pr.first_child_found_in("w:tblW")
    if tbl_w is None:
        tbl_w = OxmlElement("w:tblW")
        tbl_pr.append(tbl_w)
    tbl_w.set(qn("w:w"), "9360")
    tbl_w.set(qn("w:type"), "dxa")
    tbl_ind = tbl_pr.first_child_found_in("w:tblInd")
    if tbl_ind is None:
        tbl_ind = OxmlElement("w:tblInd")
        tbl_pr.append(tbl_ind)
    tbl_ind.set(qn("w:w"), "120")
    tbl_ind.set(qn("w:type"), "dxa")
    layout = tbl_pr.first_child_found_in("w:tblLayout")
    if layout is None:
        layout = OxmlElement("w:tblLayout")
        tbl_pr.append(layout)
    layout.set(qn("w:type"), "fixed")
    grid = table._tbl.tblGrid
    for child in list(grid):
        grid.remove(child)
    dxas = [round(w * 1440) for w in widths_in]
    for dxa in dxas:
        grid_col = OxmlElement("w:gridCol")
        grid_col.set(qn("w:w"), str(dxa))
        grid.append(grid_col)
    for row in table.rows:
        tr_pr = row._tr.get_or_add_trPr()
        cant_split = OxmlElement("w:cantSplit")
        tr_pr.append(cant_split)
        for idx, cell in enumerate(row.cells):
            cell.width = Inches(widths_in[idx])
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.TOP
            tc_pr = cell._tc.get_or_add_tcPr()
            tc_w = tc_pr.first_child_found_in("w:tcW")
            if tc_w is None:
                tc_w = OxmlElement("w:tcW")
                tc_pr.append(tc_w)
            tc_w.set(qn("w:w"), str(dxas[idx]))
            tc_w.set(qn("w:type"), "dxa")
            set_cell_margins(cell)


def font_run(run, size=None, bold=None, italic=None, color=None, name="Calibri") -> None:
    run.font.name = name
    run._element.get_or_add_rPr().get_or_add_rFonts().set(qn("w:ascii"), name)
    run._element.rPr.rFonts.set(qn("w:hAnsi"), name)
    if size is not None:
        run.font.size = Pt(size)
    if bold is not None:
        run.bold = bold
    if italic is not None:
        run.italic = italic
    if color:
        run.font.color.rgb = RGBColor.from_string(color)


def configure_doc(doc: Document) -> None:
    section = doc.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)
    section.right_margin = Inches(1)
    section.header_distance = Inches(0.492)
    section.footer_distance = Inches(0.492)
    styles = doc.styles
    normal = styles["Normal"]
    normal.font.name = "Calibri"
    normal._element.rPr.rFonts.set(qn("w:ascii"), "Calibri")
    normal._element.rPr.rFonts.set(qn("w:hAnsi"), "Calibri")
    normal.font.size = Pt(11)
    normal.paragraph_format.space_before = Pt(0)
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.10
    for name, size, color, before, after in [
        ("Heading 1", 16, COLORS["blue"], 16, 8),
        ("Heading 2", 13, COLORS["blue"], 12, 6),
        ("Heading 3", 12, COLORS["dark_blue"], 8, 4),
    ]:
        style = styles[name]
        style.font.name = "Calibri"
        style._element.rPr.rFonts.set(qn("w:ascii"), "Calibri")
        style._element.rPr.rFonts.set(qn("w:hAnsi"), "Calibri")
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = RGBColor.from_string(color)
        style.paragraph_format.space_before = Pt(before)
        style.paragraph_format.space_after = Pt(after)
        style.paragraph_format.keep_with_next = True
    if "Figure Caption" not in styles:
        styles.add_style("Figure Caption", 1)
    cap = styles["Figure Caption"]
    cap.font.name = "Calibri"
    cap.font.size = Pt(9)
    cap.font.italic = True
    cap.font.color.rgb = RGBColor.from_string(COLORS["gray"])
    cap.paragraph_format.space_before = Pt(3)
    cap.paragraph_format.space_after = Pt(8)
    cap.paragraph_format.keep_with_next = False
    if "Table Citation Text" not in styles:
        styles.add_style("Table Citation Text", 1)
    src = styles["Table Citation Text"]
    src.font.name = "Calibri"
    src.font.size = Pt(8.5)
    src.font.color.rgb = RGBColor.from_string(COLORS["gray"])
    src.paragraph_format.space_before = Pt(4)
    src.paragraph_format.space_after = Pt(4)
    configure_header_footer(section)


def configure_header_footer(section) -> None:
    header = section.header
    p = header.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    p.paragraph_format.space_after = Pt(0)
    r = p.add_run("LanM technical correction | Reproducibility package")
    font_run(r, size=8.5, color=COLORS["gray"])
    footer = section.footer
    p = footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    p.paragraph_format.space_before = Pt(0)
    r = p.add_run("Page ")
    font_run(r, size=8.5, color=COLORS["gray"])
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    p._p.append(fld)


def add_bottom_rule(paragraph, color="2E74B5", size="12") -> None:
    p_pr = paragraph._p.get_or_add_pPr()
    p_bdr = p_pr.find(qn("w:pBdr"))
    if p_bdr is None:
        p_bdr = OxmlElement("w:pBdr")
        p_pr.append(p_bdr)
    bottom = OxmlElement("w:bottom")
    bottom.set(qn("w:val"), "single")
    bottom.set(qn("w:sz"), size)
    bottom.set(qn("w:space"), "6")
    bottom.set(qn("w:color"), color)
    p_bdr.append(bottom)


def add_masthead(doc: Document) -> None:
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(12)
    p.paragraph_format.space_after = Pt(2)
    r = p.add_run("TECHNICAL CORRECTION PACKAGE")
    font_run(r, size=10.5, bold=True, color=COLORS["gold"])
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(5)
    r = p.add_run("Lanmodulin Dy-selectivity workflow")
    font_run(r, size=24, bold=True, color=COLORS["navy"])
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(14)
    r = p.add_run("Corrected figures, tables, reviewer-response map, and repository remediation")
    font_run(r, size=13, color=COLORS["gray"])
    metadata = [
        ("Manuscript", "JMGM-D-26-01758"),
        ("Scope", "Technical correction only; main manuscript intentionally unchanged"),
        ("Repository", REPO_URL),
        ("Audited base", BASE_COMMIT),
        ("Prepared", date.today().isoformat()),
    ]
    for label, value in metadata:
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(2)
        r = p.add_run(f"{label}: ")
        font_run(r, size=10.5, bold=True)
        r = p.add_run(value)
        font_run(r, size=10.5)
    rule = doc.add_paragraph()
    rule.paragraph_format.space_before = Pt(8)
    rule.paragraph_format.space_after = Pt(12)
    add_bottom_rule(rule)


def add_callout(doc, title, body, fill="E8EEF5", accent="2E74B5"):
    table = doc.add_table(rows=1, cols=1)
    table.style = "Table Grid"
    set_table_geometry(table, [6.5])
    tr_pr = table.rows[0]._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)
    cell = table.cell(0, 0)
    set_cell_shading(cell, fill)
    p = cell.paragraphs[0]
    p.paragraph_format.space_after = Pt(3)
    r = p.add_run(title)
    font_run(r, size=11, bold=True, color=accent)
    p = cell.add_paragraph()
    p.paragraph_format.space_after = Pt(0)
    r = p.add_run(body)
    font_run(r, size=10.3)
    doc.add_paragraph().paragraph_format.space_after = Pt(1)


def add_table(doc, headers, rows, widths, font_size=8.5, header_fill="F2F4F7", keep_all=False):
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    for idx, header in enumerate(headers):
        cell = table.rows[0].cells[idx]
        set_cell_shading(cell, header_fill)
        p = cell.paragraphs[0]
        p.paragraph_format.space_after = Pt(0)
        p.paragraph_format.keep_with_next = True
        r = p.add_run(str(header))
        font_run(r, size=font_size, bold=True, color=COLORS["navy"])
    for row in rows:
        cells = table.add_row().cells
        for idx, value in enumerate(row):
            p = cells[idx].paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            r = p.add_run(str(value))
            font_run(r, size=font_size)
    set_table_geometry(table, widths)
    header_pr = table.rows[0]._tr.get_or_add_trPr()
    repeat = OxmlElement("w:tblHeader")
    repeat.set(qn("w:val"), "true")
    header_pr.append(repeat)
    if keep_all:
        for row in table.rows[:-1]:
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    paragraph.paragraph_format.keep_with_next = True
    return table


def add_source(doc, text):
    p = doc.add_paragraph(style="Table Citation Text")
    p.add_run(text)


def add_figure(doc, path: Path, caption: str):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after = Pt(0)
    p.paragraph_format.keep_with_next = True
    shape = p.add_run().add_picture(str(path), width=Inches(6.35))
    shape._inline.docPr.set("title", path.stem.replace("_", " "))
    shape._inline.docPr.set("descr", caption)
    cap = doc.add_paragraph(style="Figure Caption")
    cap.add_run(caption)


def build_docx(output_path: Path, figures: dict[str, Path], tables_dir: Path) -> None:
    doc = Document()
    configure_doc(doc)
    add_masthead(doc)
    add_callout(
        doc,
        "Targeted result after correction",
        "The technically supportable result is a model-prioritized Hans pocket candidate for additional quantum and experimental testing. The present repository does not establish thermodynamic binding free energy, Dy selectivity, zero-waste separation, or a validated quantum-state mechanism.",
        fill=COLORS["pale_gold"],
        accent=COLORS["gold"],
    )
    add_table(
        doc,
        ["Area", "Status", "Submission use"],
        [
            ["Classical data correction", "Ready", "Corrected masks, funnel, lineage, model settings, OpenMM replica data, and dimensionless model score"],
            ["Quantum evidence", "Not ready", "Raw ORCA/CASSCF/CASPT2 files are absent; Figures 6-7 cannot be regenerated defensibly"],
            ["Manuscript text", "Unchanged", "Author name, title/abstract claims, citations, and prose edits are recorded but not applied"],
            ["Repository remediation", "Prepared", "New technical-correction artifacts, reproducible builder, provenance tables, and corrected README notice"],
        ],
        [1.55, 1.05, 3.9],
        font_size=9,
    )
    add_source(doc, "Status is based on repository audit at main commit " + BASE_COMMIT + ".")

    doc.add_heading("1. Scope and correction decision", level=1)
    doc.add_paragraph(
        "This package responds to the supplied reviewer comments through technical and reproducibility corrections only. It does not alter the manuscript. The package separates supported corrections from author-dependent work so that unsupported quantum values, mechanisms, or experimental claims are not introduced merely to satisfy a review point."
    )
    doc.add_heading("1.1 What can be corrected now", level=2)
    doc.paragraphs[-1].insert_paragraph_before().add_run().add_break(WD_BREAK.PAGE)
    add_table(
        doc,
        ["Correction", "Concrete update"],
        [
            ["Metric terminology", "Remove the energy-like transformation and report a dimensionless Composite Quality Index."],
            ["Design accounting", "Reconcile 70 role memberships to 66 unique positions and separate positions from generated sequences."],
            ["Sampling transparency", "Report exact MPNN temperatures, seeds, counts, and shortlist units."],
            ["Rosetta description", "Replace structural-realism wording with one-structure score-only triage and explicit limitations."],
            ["OpenMM reporting", "Expose all 45 replica/metal/candidate means, parameter values, seeds, duration, metric definition, and force-field caveats."],
            ["Structure scope", "Distinguish deposited A-D chains from author/PISA A2 assemblies and clarify chain-A versus multichain computational contexts."],
        ],
        [1.65, 4.85],
        font_size=9,
    )
    add_source(doc, "Full machine-readable tables are included in the tables/ directory.")
    doc.add_heading("1.2 What cannot be corrected without author data", level=2)
    doc.add_paragraph(
        "The repository contains no ORCA inputs or outputs and no files named for CASSCF or CASPT2. Therefore natural occupations, active-space root character, CASPT2 settings, intruder-state diagnostics, cluster sizes, relativistic treatment, and spin-orbit coupling cannot be reconstructed. Existing quantum Figures 6-7 should not be submitted as regenerated figures until the source files listed in Section 6 are deposited."
    )

    doc.add_heading("2. Reviewer-response correction matrix", level=1)
    doc.add_paragraph("The matrix below is the shareable point-by-point technical response. 'Blocked' means the correction requires source data that are not in the repository; it does not mean a value should be estimated.")
    add_table(
        doc,
        ["Item", "Reviewer issue", "Technical correction", "Status"],
        [[row["reviewer_item"], row["issue"], row["technical_correction"] + " Evidence: " + row["deliverable"] + ".", row["status"]] for row in reviewer_rows()],
        [0.75, 1.55, 2.95, 1.25],
        font_size=7.35,
        header_fill="E8EEF5",
    )
    add_source(doc, "Machine-readable version: tables/reviewer_response_matrix.csv.")

    doc.add_page_break()
    doc.add_heading("3. Corrected figures", level=1)
    add_figure(doc, figures["tc1"], "Figure TC1. Corrected design-mask accounting. Role categories overlap: 70 memberships correspond to 66 unique residue positions, with four one-count overlaps.")
    add_figure(doc, figures["tc2"], "Figure TC2. Corrected design funnel. Raw generations, unique sequences, shortlisted sequences, scored candidates, and final OpenMM finalists are reported separately.")
    add_figure(doc, figures["tc3"], "Figure TC3. Corrected Hans model scope and Arg100 labeling. The diagram distinguishes the crystallographic A-D context from A2 biological assemblies and explicitly tracks chain-A Arg100/Lys100.")
    doc.add_page_break()
    add_figure(doc, figures["tc4"], "Figure TC4. Corrected round-2 OpenMM geometry screen. Every point is a replica mean. The Al/Fe values are force-field-dependent geometry results and are not interpreted as affinity or selectivity.")
    add_figure(doc, figures["tc5"], "Figure TC5. Corrected model-quality prioritization. The Composite Quality Index is dimensionless and must not be labeled as Delta G, free energy, or thermodynamic affinity.")

    doc.add_page_break()
    doc.add_heading("4. Corrected tables and exact technical settings", level=1)
    doc.add_heading("4.1 Design-mask reconciliation", level=2)
    add_table(doc, ["Role", "Memberships"], [[r["role"], r["membership_count"]] for r in MASK_ROWS] + [["Unique inventory", 66], ["Membership sum", 70], ["Overlap excess", 4]], [4.8, 1.7], font_size=9, keep_all=True)
    add_source(doc, "Source: design_mask_candidates.csv; overlap detail is in mask_overlaps.csv.")

    doc.add_heading("4.2 MPNN and Rosetta specification", level=2)
    add_table(doc, ["Component", "Exact setting", "Interpretation limit"], [[r["component"], r["exact_setting"], r["interpretation_limit"]] for r in MODEL_SPEC_ROWS[:4]], [1.35, 2.95, 2.20], font_size=8.3, keep_all=True)
    add_source(doc, "Sources: recorded run_command.txt files; rosetta_round2_candidate_ranking.csv; repository scoring code.")

    doc.add_heading("4.3 OpenMM protocol and metric", level=2)
    add_table(doc, ["Component", "Exact setting", "Interpretation limit"], [[r["component"], r["exact_setting"].replace("3.2 A", "3.2 Å"), r["interpretation_limit"]] for r in MODEL_SPEC_ROWS[4:]], [1.35, 2.95, 2.20], font_size=8.3, keep_all=True)
    add_source(doc, "Sources: simulation_config.yaml, system_build_log.txt, and 45 screening_log.txt files under results/round2_openmm_screening/.")

    doc.add_heading("4.4 Generic OPC3-compatible 12-6-4 metal parameters", level=2)
    add_table(
        doc,
        ["Metal", "Charge", "Rmin/2 (Å)", "epsilon (kcal/mol)", "C4 (kcal/mol Å^4)"],
        [[r["metal"], r["formal_charge"], r["rmin_half_A"], r["epsilon_kcal_mol"], r["c4_kcal_mol_A4"]] for r in METAL_PARAMETERS],
        [0.75, 0.75, 1.25, 1.75, 2.0],
        font_size=8.7,
        keep_all=True,
    )
    add_source(doc, "Source: config/metal_parameter_values.yaml. The repository does not contain chelator-tuned numeric coefficients.")

    doc.add_heading("4.5 OpenMM candidate means", level=2)
    summary = openmm_summary_rows()
    rows = []
    for candidate in OPENMM_REPLICATES:
        vals = {r["metal"]: r for r in summary if r["candidate_id"] == candidate}
        rows.append([CANDIDATE_LABELS[candidate]] + [f"{float(vals[m]['mean_distance_A']):.3f} ± {float(vals[m]['sample_sd_A']):.4f}" for m in ["Dy", "Nd", "Y", "Al", "Fe"]])
    add_table(doc, ["Candidate", "Dy", "Nd", "Y", "Al", "Fe"], rows, [1.4, 1.0197, 1.0203, 1.02, 1.02, 1.02], font_size=8.3, keep_all=True)
    add_source(doc, "Values are mean ± sample SD of three replica-level means in Å. Full six-decimal data are in openmm_replicates.csv and openmm_summary.csv.")
    doc.add_paragraph(
        "The near-identical Y values are expected from the common Y(III) parameter model and similar first-shell coordination. The underlying nine replica means range from 2.183940 to 2.196893 Å; the repeated '2.19 Å' statement was only rounded display."
    )

    doc.add_heading("4.6 Dimensionless model-quality index", level=2)
    add_table(doc, ["Rank", "Candidate", "Ligand confidence", "Overall confidence", "CQI"], [[r["rank"], r["candidate_label"], r["ligand_confidence"], r["overall_confidence"], r["composite_quality_index"]] for r in model_score_rows()], [0.55, 1.85, 1.35, 1.35, 1.40], font_size=8.7, keep_all=True)
    add_source(doc, "CQI = 0.7 x ligand confidence + 0.3 x overall confidence. It is dimensionless and has no thermodynamic interpretation.")

    doc.add_heading("5. Structure and model-scope correction", level=1)
    add_table(doc, ["PDB", "Deposited model", "Biological assembly", "Technical use"], [[r["pdb_id"], r["deposited_chains"] + "; " + r["deposited_metal_inventory"], r["author_biological_assembly"], r["technical_use"]] for r in STRUCTURAL_SCOPE], [0.55, 2.05, 2.0, 1.9], font_size=8.0, keep_all=True)
    add_source(doc, f"Primary structure sources: {RCSB_8FNR} and {RCSB_8DQ2}; repository inventory: results/tables/template_chain_summary.csv.")
    doc.add_paragraph(
        "The deposited metal counts are crystallographic inventories, not numbers of equivalent functional sites. Any interface conclusion should identify the specific biological assembly or explicitly state that A-D was retained as crystallographic context."
    )

    doc.add_heading("6. Quantum evidence gap and submission gate", level=1)
    add_callout(doc, "Submission gate", "Do not present Figures 6-7, natural occupations, CASPT2 state ordering, spin-orbit conclusions, or reduced-cluster sensitivity as repository-reproduced results until the corresponding raw files are added and independently regenerated.", fill=COLORS["pale_red"], accent=COLORS["red"])
    add_table(doc, ["Required artifact", "Minimum content", "Status", "Affected output"], [[r["required_artifact"], r["minimum_content"], r["repository_status"], r["affected_claim_or_figure"]] for r in QUANTUM_GAPS], [1.45, 2.75, 0.75, 1.55], font_size=7.8, header_fill="FBEAEA")
    add_source(doc, "Repository tree audit: no paths containing ORCA, CASSCF, or CASPT2 and no .in/.inp/.out quantum files were present at the audited base commit.")

    doc.add_heading("7. Repository correction set", level=1)
    add_table(
        doc,
        ["Path", "Purpose"],
        [
            ["README.md", "Replaces the complete/reproducible and Delta-G language with a technical-correction notice and accurate limitations."],
            ["docs/technical_corrections/README.md", "Human-readable scope, reviewer mapping, and regeneration instructions."],
            ["docs/technical_corrections/SOURCE_DATA_GAPS.md", "Explicit quantum and manuscript-only actions that remain open."],
            ["docs/technical_corrections/LanM_Technical_Corrections_Report.{docx,pdf}", "Shareable Word and PDF report."],
            ["results/technical_corrections/figures/*", "Five corrected technical figures."],
            ["results/technical_corrections/tables/*", "Machine-readable corrected data and provenance."],
            ["scripts/build_technical_corrections.py", "Deterministic rebuild of tables, figures, report, and manifest."],
            ["requirements-technical-corrections.txt", "Pinned reporting dependencies."],
        ],
        [3.35, 3.15],
        font_size=8.5,
    )
    add_source(doc, "The manuscript file and its figures are intentionally not edited by this repository correction set.")

    doc.add_heading("8. Package manifest", level=1)
    doc.add_paragraph(
        "The submission bundle includes this Word report, a PDF rendering, five corrected figures, twelve corrected/source-audit CSV tables, a figure/table manifest, repository-ready documentation, a reproducible build script, pinned reporting requirements, and SHA-256 checksums. The original reviewer PDF and manuscript are not redistributed in the generated ZIP."
    )
    add_callout(doc, "Recommended share note", "This technical package resolves the reproducibility, terminology, sampling, model-scope, and classical-data presentation issues that can be corrected from the repository. Quantum-source-data and manuscript-only edits remain explicit author actions.", fill=COLORS["blue_gray"], accent=COLORS["blue"])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    doc.save(output_path)


def write_repository_docs(repo_dir: Path) -> None:
    readme = f"""# LanM technical corrections

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

Structure scope is cross-checked against [RCSB 8FNR]({RCSB_8FNR}) and [RCSB 8DQ2]({RCSB_8DQ2}).
"""
    gaps = """# Source-data gaps and deferred manuscript actions

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
"""
    notice = f"""# lanm-dy-selectivity

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
"""
    files = {
        repo_dir / "README.md": notice,
        repo_dir / "docs/technical_corrections/README.md": readme,
        repo_dir / "docs/technical_corrections/SOURCE_DATA_GAPS.md": gaps,
        repo_dir / "requirements-technical-corrections.txt": "matplotlib==3.10.5\nnumpy==2.3.2\npython-docx==1.2.0\nPillow==11.3.0\n",
    }
    for path, text in files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")


def copy_repo_artifacts(repo_dir: Path, output_root: Path, script_path: Path) -> None:
    target_fig = repo_dir / "results/technical_corrections/figures"
    target_tab = repo_dir / "results/technical_corrections/tables"
    target_doc = repo_dir / "docs/technical_corrections"
    target_script = repo_dir / "scripts/build_technical_corrections.py"
    target_fig.mkdir(parents=True, exist_ok=True)
    target_tab.mkdir(parents=True, exist_ok=True)
    target_doc.mkdir(parents=True, exist_ok=True)
    target_script.parent.mkdir(parents=True, exist_ok=True)
    for path in (output_root / "figures").glob("*"):
        shutil.copy2(path, target_fig / path.name)
    for path in (output_root / "tables").glob("*"):
        shutil.copy2(path, target_tab / path.name)
    shutil.copy2(output_root / "LanM_Technical_Corrections_Report.docx", target_doc / "LanM_Technical_Corrections_Report.docx")
    shutil.copy2(output_root / "LanM_Technical_Corrections_Report.pdf", target_doc / "LanM_Technical_Corrections_Report.pdf")
    shutil.copy2(script_path, target_script)


def render_docx(docx_path: Path, render_dir: Path) -> Path | None:
    render_dir.mkdir(parents=True, exist_ok=True)
    renderer = Path("/root/.codex/skills/builtins/documents/render_docx.py")
    cmd = [sys.executable, str(renderer), str(docx_path), "--output_dir", str(render_dir), "--emit_pdf"]
    completed = subprocess.run(cmd, check=False, capture_output=True, text=True)
    (render_dir / "render.log").write_text(completed.stdout + "\n" + completed.stderr, encoding="utf-8")
    pdf = render_dir / (docx_path.stem + ".pdf")
    return pdf if pdf.exists() else None


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(root: Path) -> None:
    rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file() and p.name not in {"CHECKSUMS.sha256", "manifest.json"} and "rendered_report" not in p.parts):
        rows.append({"path": str(path.relative_to(root)), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    (root / "manifest.json").write_text(json.dumps({"package": "LanM technical corrections", "base_commit": BASE_COMMIT, "files": rows}, indent=2) + "\n", encoding="utf-8")
    checksum_lines = [f"{row['sha256']}  {row['path']}" for row in rows]
    (root / "CHECKSUMS.sha256").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")


def write_package_readme(root: Path) -> None:
    text = """# LanM technical-correction submission package

This package contains the reviewer-driven technical corrections that can be supported by the audited repository. It does not modify or redistribute the main manuscript.

## Use these files for sharing

- `LanM_Technical_Corrections_Report.docx` — point-by-point correction report with corrected figures and summary tables.
- `LanM_Technical_Corrections_Report.pdf` — rendered share copy.
- `figures/` — five corrected technical figures at publication-ready raster resolution.
- `tables/` — machine-readable correction, provenance, and source-gap tables.
- `FIGURE_TABLE_MANIFEST.md` — quick mapping of every figure and table.
- `repository_ready/` — exact repository files prepared for the correction branch.

## Scientific boundary

The corrected result is a model-prioritized candidate for further QM and experimental testing. Missing ORCA/CASSCF/CASPT2 files prevent defensible regeneration of the manuscript's quantum figures. Those gaps are documented rather than filled with estimates.
"""
    (root / "README.md").write_text(text, encoding="utf-8")
    manifest = """# Figure and table manifest

## Corrected figures

| File | Correction |
| --- | --- |
| Figure_TC1_design_mask_inventory.png | Reconciles 70 role memberships to 66 unique positions. |
| Figure_TC2_sequence_design_funnel.png | Separates generated, unique, retained, scored, and finalist units. |
| Figure_TC3_hans_scope_and_arg100.png | Clarifies biological/crystallographic scope and labels Arg100 lineage. |
| Figure_TC4_openmm_multimetal_distances.png | Shows all replica-level means and force-field interpretation limits. |
| Figure_TC5_composite_quality_index.png | Replaces invalid Delta-G language with a dimensionless model-quality index. |

## Corrected/source-audit tables

| File | Contents |
| --- | --- |
| reviewer_response_matrix.csv | Point-by-point reviewer response and status. |
| design_mask_inventory.csv | Role memberships and unique inventory. |
| mask_overlaps.csv | Four overlapping positions that explain 70 versus 66. |
| sequence_funnel.csv | Exact funnel counts and units. |
| finalist_lineage.csv | Unambiguous chain-A mutation lineage. |
| modeling_specification.csv | MPNN, Rosetta, OpenMM, and force-field settings. |
| openmm_replicates.csv | All 45 replica/metal/candidate results. |
| openmm_summary.csv | Mean and sample SD for each condition. |
| metal_parameter_values.csv | Generic OPC3-compatible 12-6-4 values. |
| model_quality_index.csv | Corrected dimensionless prioritization. |
| structural_scope.csv | Deposited versus biological assembly and computational-use scope. |
| quantum_evidence_gap_register.csv | Missing source files and affected claims/figures. |
"""
    (root / "FIGURE_TABLE_MANIFEST.md").write_text(manifest, encoding="utf-8")


def build(output_root: Path, script_path: Path) -> dict[str, Path]:
    if output_root.exists():
        shutil.rmtree(output_root)
    figures_dir = output_root / "figures"
    tables_dir = output_root / "tables"
    figures_dir.mkdir(parents=True)
    tables_dir.mkdir(parents=True)
    table_map = {
        "reviewer_response_matrix.csv": reviewer_rows(),
        "design_mask_inventory.csv": MASK_ROWS,
        "mask_overlaps.csv": OVERLAP_ROWS,
        "sequence_funnel.csv": FUNNEL_ROWS,
        "finalist_lineage.csv": LINEAGE_ROWS,
        "modeling_specification.csv": MODEL_SPEC_ROWS,
        "openmm_replicates.csv": openmm_replica_rows(),
        "openmm_summary.csv": openmm_summary_rows(),
        "metal_parameter_values.csv": METAL_PARAMETERS,
        "model_quality_index.csv": model_score_rows(),
        "structural_scope.csv": STRUCTURAL_SCOPE,
        "quantum_evidence_gap_register.csv": QUANTUM_GAPS,
    }
    for name, rows in table_map.items():
        write_csv(tables_dir / name, rows)
    figures = {
        "tc1": figures_dir / "Figure_TC1_design_mask_inventory.png",
        "tc2": figures_dir / "Figure_TC2_sequence_design_funnel.png",
        "tc3": figures_dir / "Figure_TC3_hans_scope_and_arg100.png",
        "tc4": figures_dir / "Figure_TC4_openmm_multimetal_distances.png",
        "tc5": figures_dir / "Figure_TC5_composite_quality_index.png",
    }
    fig_mask(figures["tc1"])
    fig_funnel(figures["tc2"])
    fig_scope_arg100(figures["tc3"])
    fig_openmm(figures["tc4"])
    fig_cqi(figures["tc5"])
    report = output_root / "LanM_Technical_Corrections_Report.docx"
    build_docx(report, figures, tables_dir)
    render_dir = output_root.parent / (output_root.name + "_rendered_QA")
    pdf = render_docx(report, render_dir)
    if pdf:
        shutil.copy2(pdf, output_root / "LanM_Technical_Corrections_Report.pdf")
    if render_dir.exists():
        shutil.rmtree(render_dir)
    write_package_readme(output_root)
    repo_dir = output_root / "repository_ready"
    write_repository_docs(repo_dir)
    copy_repo_artifacts(repo_dir, output_root, script_path)
    write_manifest(output_root)
    return {"report": report, "pdf": output_root / "LanM_Technical_Corrections_Report.pdf", "repo_dir": repo_dir, **figures}


def make_zip(output_root: Path, zip_path: Path) -> None:
    if zip_path.exists():
        zip_path.unlink()
    base = zip_path.with_suffix("")
    archive = shutil.make_archive(str(base), "zip", root_dir=output_root.parent, base_dir=output_root.name)
    if Path(archive) != zip_path:
        shutil.move(archive, zip_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--zip-path", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    script_path = Path(__file__).resolve()
    outputs = build(args.output_root.resolve(), script_path)
    if args.zip_path:
        make_zip(args.output_root.resolve(), args.zip_path.resolve())
    print(json.dumps({k: str(v) for k, v in outputs.items()}, indent=2))


if __name__ == "__main__":
    main()
