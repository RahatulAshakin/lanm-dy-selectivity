"""Build the submission-ready reporting package from completed workflow outputs."""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from statistics import mean

from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_bytes, atomic_write_text, write_csv_rows
from lanm.paths import REPO_ROOT

FINALIST_CANDIDATE_IDS = (
    "hans_pocket_ss_only_u02_r2u01",
    "hans_interface_ss_plus_if_u04_r2u01",
    "am1_mex_ss_only_u02_r2u02",
)
TEMPLATE_ORDER = ("6MI5", "8FNS", "8DQ2", "8FNR")
METAL_ORDER = ("Dy", "Nd", "Y", "Al", "Fe")
ROLE_ORDER = ("first_shell", "second_sphere", "solvent_contact", "interchain_contact")
ROLE_COLORS = {
    "first_shell": "#8c2f39",
    "second_sphere": "#2f6f8f",
    "solvent_contact": "#7a8c99",
    "interchain_contact": "#4f7f3b",
}
TEMPLATE_FAMILY_LABELS = {
    "6MI5": "AM1/Mex",
    "8FNS": "AM1/Mex",
    "8DQ2": "Hans",
    "8FNR": "Hans",
}
TEMPLATE_TOPOLOGY_LABELS = {
    "6MI5": "am1_monomer_reference",
    "8FNS": "am1_monomer_reference",
    "8DQ2": "hans_multichain_reference",
    "8FNR": "hans_multichain_reference",
}
TOPOLOGY_ORDER = {
    "am1_monomer": 0,
    "hans_monomer": 1,
    "hans_interface_multichain": 2,
}
COMMON_INTERPRETATION_PARAGRAPH = (
    "Rare-earth retention was preserved across tested candidates. Persistent Al/Fe inner-sphere "
    "capture remained unresolved. No candidate advanced to full QM/QCT in the validated pipeline. "
    "The Main Quantum Pocket ranking is a deterministic next-step prioritization based on completed "
    "workflow outputs."
)


@dataclass(frozen=True, slots=True)
class FinalProjectReportInputs:
    template_chain_summary_path: Path
    residue_role_map_path: Path
    design_mask_candidates_path: Path
    proteinmpnn_shortlist_path: Path
    ligandmpnn_shortlist_path: Path
    integrated_candidate_ranking_path: Path
    md_validation_panel_path: Path
    proteinmpnn_round2_shortlist_path: Path
    ligandmpnn_round2_smoke_summary_path: Path
    ligandmpnn_round2_shortlist_path: Path
    rosetta_round2_candidate_ranking_path: Path
    round2_md_rescreen_panel_path: Path
    round2_openmm_screening_summary_path: Path
    round2_openmm_screening_panel_status_path: Path


@dataclass(frozen=True, slots=True)
class FinalProjectReportOutputs:
    final_project_report_path: Path
    project_summary_report_path: Path
    project_detailed_report_path: Path
    paper_short_path: Path
    paper_long_path: Path
    readme_path: Path
    final_table_1_path: Path
    final_table_2_path: Path
    final_table_3_path: Path
    final_figure_1_path: Path
    final_figure_2_path: Path
    final_figure_3_path: Path
    final_figure_4_path: Path


@dataclass(frozen=True, slots=True)
class FinalProjectBuildResult:
    main_quantum_pocket_candidate_id: str
    pocket_rows: tuple[dict[str, str], ...]
    output_paths: tuple[Path, ...]


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _load_csv_rows(path: Path) -> tuple[dict[str, str], ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return tuple({str(key): "" if value is None else str(value) for key, value in row.items()} for row in reader)


def _parse_int(value: str) -> int:
    return int(str(value).strip())


def _parse_float(value: str) -> float:
    return float(str(value).strip())


def _format_float(value: float, digits: int = 6) -> str:
    return f"{value:.{digits}f}"


def _format_metric(value: float | None, digits: int = 3) -> str:
    if value is None:
        return "-"
    return f"{value:.{digits}f}"


def _render_markdown_table(headers: tuple[str, ...], rows: tuple[tuple[str, ...], ...]) -> str:
    header_line = "| " + " | ".join(headers) + " |"
    divider_line = "| " + " | ".join(["---"] * len(headers)) + " |"
    body_lines = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join((header_line, divider_line, *body_lines))


def _template_sort_key(template_id: str) -> tuple[int, str]:
    try:
        return (TEMPLATE_ORDER.index(template_id), template_id)
    except ValueError:
        return (len(TEMPLATE_ORDER), template_id)


def _topology_sort_key(topology_class: str) -> tuple[int, str]:
    return (TOPOLOGY_ORDER.get(topology_class, 99), topology_class)


def _infer_topology_class(candidate_id: str, campaign_id: str, backbone_id: str) -> str:
    joined = "|".join((candidate_id, campaign_id, backbone_id)).lower()
    if "hans_interface" in joined:
        return "hans_interface_multichain"
    if "hans_pocket" in joined:
        return "hans_monomer"
    if "am1_mex" in joined:
        return "am1_monomer"
    return "unassigned_topology"


def _candidate_display_label(candidate_id: str) -> str:
    mapping = {
        "hans_pocket_ss_only_u02_r2u01": "Hans pocket r2u01",
        "hans_interface_ss_plus_if_u04_r2u01": "Hans interface r2u01",
        "am1_mex_ss_only_u02_r2u02": "AM1/Mex r2u02",
    }
    return mapping.get(candidate_id, candidate_id)


def _preferred_mutation_string(
    candidate_id: str,
    *,
    ligand_round2_by_id: dict[str, dict[str, str]],
    protein_round2_by_id: dict[str, dict[str, str]],
    ligand_round1_by_id: dict[str, dict[str, str]],
    protein_round1_by_id: dict[str, dict[str, str]],
) -> str:
    for source in (ligand_round2_by_id, protein_round2_by_id, ligand_round1_by_id, protein_round1_by_id):
        row = source.get(candidate_id)
        if row is not None:
            mutation_string = str(row.get("mutation_string", "")).strip()
            return mutation_string or "none"
    return "none"


def _mean_distance_by_candidate(
    summary_rows_by_key: dict[tuple[str, str], dict[str, str]],
    candidate_id: str,
    metals: tuple[str, ...],
) -> float:
    values = [
        _parse_float(summary_rows_by_key[(candidate_id, metal)]["median_mean_min_metal_oxygen_distance_A"])
        for metal in metals
    ]
    return float(mean(values))


def build_final_template_summary_table(
    template_chain_rows: tuple[dict[str, str], ...],
) -> tuple[dict[str, str], ...]:
    grouped_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in template_chain_rows:
        grouped_rows[row["template_id"]].append(row)

    table_rows: list[dict[str, str]] = []
    for template_id in sorted(grouped_rows, key=_template_sort_key):
        rows = grouped_rows[template_id]
        chain_ids = ",".join(sorted({str(row["chain_id"]).strip() for row in rows}))
        metal_identities = ",".join(sorted({str(row["metal_identity"]).strip() for row in rows}))
        metal_site_count = sum(_parse_int(row["metal_site_count"]) for row in rows)
        table_rows.append(
            {
                "template_id": template_id,
                "family": TEMPLATE_FAMILY_LABELS.get(template_id, "Unassigned"),
                "topology_class": TEMPLATE_TOPOLOGY_LABELS.get(template_id, "unassigned_template_topology"),
                "source_type": str(rows[0]["source_type"]).strip(),
                "chain_ids": chain_ids,
                "metal_identity": metal_identities,
                "metal_site_count": str(metal_site_count),
            }
        )
    return tuple(table_rows)


def build_candidate_progression_table(
    *,
    protein_round1_rows: tuple[dict[str, str], ...],
    ligand_round1_rows: tuple[dict[str, str], ...],
    integrated_rows: tuple[dict[str, str], ...],
    md_panel_rows: tuple[dict[str, str], ...],
    protein_round2_rows: tuple[dict[str, str], ...],
    ligand_round2_rows: tuple[dict[str, str], ...],
    finalist_rows: tuple[dict[str, str], ...],
) -> tuple[dict[str, str], ...]:
    ligand_round1_by_id = {row["candidate_id"]: row for row in ligand_round1_rows}
    ligand_round2_by_id = {row["candidate_id"]: row for row in ligand_round2_rows}
    protein_round1_by_id = {row["candidate_id"]: row for row in protein_round1_rows}
    protein_round2_by_id = {row["candidate_id"]: row for row in protein_round2_rows}
    integrated_ids = {row["candidate_id"] for row in integrated_rows}
    md_designed_ids = {
        row["candidate_id"]
        for row in md_panel_rows
        if str(row.get("panel_member_type", "")).strip() == "designed_candidate"
    }
    round2_seed_ids = {row["seed_candidate_id"] for row in protein_round2_rows}
    finalist_ids = {row["candidate_id"] for row in finalist_rows}
    ligand_round1_ids = set(ligand_round1_by_id)
    ligand_round2_ids = set(ligand_round2_by_id)

    output_rows: list[dict[str, str]] = []
    for row in sorted(protein_round1_rows, key=lambda item: _parse_int(item["shortlist_rank"])):
        candidate_id = row["candidate_id"]
        if candidate_id in round2_seed_ids:
            status = "retained"
            note = "Retained through round-1 prioritization and promoted as a round-2 seed scaffold."
        elif candidate_id in md_designed_ids:
            status = "rejected"
            note = "Retained into the round-1 MD panel but not promoted into the three-seed round-2 redesign set."
        elif candidate_id in integrated_ids:
            status = "rejected"
            note = "Entered the integrated round-1 ranking but was not kept in the MD or round-2 seed set."
        elif candidate_id in ligand_round1_ids:
            status = "rejected"
            note = "Produced a ligand-supported round-1 representative but was not retained in the integrated ranking."
        else:
            status = "rejected"
            note = "Did not survive beyond the round-1 ProteinMPNN shortlist."
        output_rows.append(
            {
                "candidate_id": candidate_id,
                "phase_origin": "round1_shortlist",
                "backbone_id": row["backbone_id"],
                "topology_class": _infer_topology_class(candidate_id, row["campaign_id"], row["backbone_id"]),
                "key_mutations": _preferred_mutation_string(
                    candidate_id,
                    ligand_round2_by_id=ligand_round2_by_id,
                    protein_round2_by_id=protein_round2_by_id,
                    ligand_round1_by_id=ligand_round1_by_id,
                    protein_round1_by_id=protein_round1_by_id,
                ),
                "retained_or_rejected": status,
                "rejection_reason_or_progress_note": note,
            }
        )

    for row in sorted(protein_round2_rows, key=lambda item: _parse_int(item["shortlist_rank"])):
        candidate_id = row["candidate_id"]
        if candidate_id in finalist_ids:
            status = "retained"
            note = "Retained through round-2 LigandMPNN/Rosetta triage and the reduced OpenMM rescreen finalist panel."
        elif candidate_id in ligand_round2_ids:
            status = "rejected"
            note = "Reached the round-2 LigandMPNN shortlist but was not retained in the three-finalist reduced rescreen panel."
        else:
            status = "rejected"
            note = "Generated in the round-2 ProteinMPNN redesign branch but did not survive the reduced LigandMPNN shortlist."
        output_rows.append(
            {
                "candidate_id": candidate_id,
                "phase_origin": "round2_redesign",
                "backbone_id": row["backbone_id"],
                "topology_class": row["topology_class"],
                "key_mutations": _preferred_mutation_string(
                    candidate_id,
                    ligand_round2_by_id=ligand_round2_by_id,
                    protein_round2_by_id=protein_round2_by_id,
                    ligand_round1_by_id=ligand_round1_by_id,
                    protein_round1_by_id=protein_round1_by_id,
                ),
                "retained_or_rejected": status,
                "rejection_reason_or_progress_note": note,
            }
        )

    return tuple(output_rows)


def _validate_finalist_sets(
    finalist_rows: tuple[dict[str, str], ...],
    panel_status_rows: tuple[dict[str, str], ...],
) -> None:
    finalist_ids = {row["candidate_id"] for row in finalist_rows}
    expected_ids = set(FINALIST_CANDIDATE_IDS)
    if finalist_ids != expected_ids:
        raise ValueError(
            "Round-2 finalist set mismatch for final reporting package: "
            f"expected {FINALIST_CANDIDATE_IDS}, found {tuple(sorted(finalist_ids))}"
        )

    status_by_candidate = {row["candidate_id"]: row for row in panel_status_rows}
    for candidate_id in FINALIST_CANDIDATE_IDS:
        status_row = status_by_candidate.get(candidate_id)
        if status_row is None:
            raise ValueError(f"Missing round-2 panel status row for finalist {candidate_id}")
        rare_earth_statuses = tuple(status_row[f"{metal.lower()}_status"] for metal in ("dy", "nd", "y"))
        if rare_earth_statuses != ("stable_bound", "stable_bound", "stable_bound"):
            raise ValueError(
                f"Finalist {candidate_id} no longer preserves Dy/Nd/Y retention in the reduced rescreen: "
                f"{rare_earth_statuses}"
            )
        if status_row["al_status"] != "persistent_capture_flag" or status_row["fe_status"] != "persistent_capture_flag":
            raise ValueError(
                f"Finalist {candidate_id} no longer matches the unresolved Al/Fe interpretation: "
                f"Al={status_row['al_status']}, Fe={status_row['fe_status']}"
            )


def build_specific_pocket_thermodynamics_table(
    finalist_rows: tuple[dict[str, str], ...],
    panel_status_rows: tuple[dict[str, str], ...],
) -> tuple[dict[str, str], ...]:
    _validate_finalist_sets(finalist_rows, panel_status_rows)
    finalist_rows_by_id = {row["candidate_id"]: row for row in finalist_rows}
    panel_status_by_id = {row["candidate_id"]: row for row in panel_status_rows}

    unsorted_rows: list[dict[str, str] | tuple[str, float, float, float, dict[str, str]]] = []
    for candidate_id in FINALIST_CANDIDATE_IDS:
        finalist_row = finalist_rows_by_id[candidate_id]
        status_row = panel_status_by_id[candidate_id]
        ligand_confidence = _parse_float(finalist_row["representative_ligand_confidence"])
        overall_confidence = _parse_float(finalist_row["representative_overall_confidence"])
        pocket_confidence = 0.7 * ligand_confidence + 0.3 * overall_confidence
        delta_g = -0.592 * math.log(max(pocket_confidence, 1e-6))
        unsorted_rows.append((candidate_id, ligand_confidence, overall_confidence, pocket_confidence, delta_g, finalist_row, status_row))

    ranked_rows: list[dict[str, str]] = []
    for rank, payload in enumerate(
        sorted(
            unsorted_rows,
            key=lambda item: (item[4], item[0]),
        ),
        start=1,
    ):
        candidate_id, ligand_confidence, overall_confidence, pocket_confidence, delta_g, finalist_row, status_row = payload
        if rank == 1:
            qm_note = (
                "Main Quantum Pocket; deterministic next-step QM/QCT priority based on completed workflow outputs. "
                "No full QM/QCT completed."
            )
        else:
            qm_note = (
                f"Pocket-priority rank {rank}; deterministic next-step QM/QCT follow-up after the Main Quantum Pocket. "
                "No full QM/QCT completed."
            )
        ranked_rows.append(
            {
                "pocket_rank": str(rank),
                "candidate_id": candidate_id,
                "backbone_id": finalist_row["backbone_id"],
                "topology_class": finalist_row["topology_class"],
                "representative_ligand_confidence": _format_float(ligand_confidence, 4),
                "representative_overall_confidence": _format_float(overall_confidence, 4),
                "pocket_confidence": _format_float(pocket_confidence, 6),
                "deltaG_pocket_proxy_kcal_mol": _format_float(delta_g, 6),
                "rare_earth_retention_status": (
                    f"Dy={status_row['dy_status']}; Nd={status_row['nd_status']}; Y={status_row['y_status']}"
                ),
                "Al_status": status_row["al_status"],
                "Fe_status": status_row["fe_status"],
                "qm_priority_note": qm_note,
            }
        )
    return tuple(ranked_rows)


def _save_matplotlib_figure(fig: object, output_path: Path) -> None:
    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=200, metadata={"Software": "lanm-dy-selectivity"})
    atomic_write_bytes(output_path, buffer.getvalue())


def _import_pyplot() -> object:
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    return plt


def render_final_figure_1_template_harmonization(
    *,
    table_1_rows: tuple[dict[str, str], ...],
    residue_role_rows: tuple[dict[str, str], ...],
    output_path: Path,
) -> None:
    plt = _import_pyplot()
    role_counts_by_template: dict[str, Counter[str]] = defaultdict(Counter)
    for row in residue_role_rows:
        role_counts_by_template[row["template_id"]][row["role"]] += 1

    template_ids = [row["template_id"] for row in table_1_rows]
    x_values = list(range(len(template_ids)))
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)

    running_bottoms = [0] * len(template_ids)
    for role in ROLE_ORDER:
        values = [role_counts_by_template[template_id].get(role, 0) for template_id in template_ids]
        axes[0].bar(
            x_values,
            values,
            bottom=running_bottoms,
            color=ROLE_COLORS[role],
            label=role,
            width=0.7,
        )
        running_bottoms = [bottom + value for bottom, value in zip(running_bottoms, values)]
    axes[0].set_title("Residue-role evidence across harmonized templates")
    axes[0].set_ylabel("Observation count")
    axes[0].set_xticks(
        x_values,
        labels=[f"{row['template_id']}\n{row['family']}" for row in table_1_rows],
    )
    axes[0].grid(axis="y", alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")

    chain_counts = [len(row["chain_ids"].split(",")) if row["chain_ids"] else 0 for row in table_1_rows]
    metal_site_counts = [_parse_int(row["metal_site_count"]) for row in table_1_rows]
    width = 0.34
    axes[1].bar(
        [value - width / 2 for value in x_values],
        chain_counts,
        width=width,
        color="#3d4f6a",
        label="chain count",
    )
    axes[1].bar(
        [value + width / 2 for value in x_values],
        metal_site_counts,
        width=width,
        color="#d95f02",
        label="metal-site count",
    )
    axes[1].set_title("Template coverage and metal-site inventory")
    axes[1].set_ylabel("Count")
    axes[1].set_xticks(
        x_values,
        labels=[f"{row['template_id']}\n{row['chain_ids']}" for row in table_1_rows],
    )
    axes[1].grid(axis="y", alpha=0.2)
    axes[1].legend(frameon=False)

    fig.suptitle("Final Figure 1. Template harmonization and residue-role summary", fontsize=14, fontweight="bold")
    _save_matplotlib_figure(fig, output_path)
    plt.close(fig)


def render_final_figure_2_design_funnel(
    *,
    design_mask_rows: tuple[dict[str, str], ...],
    protein_round1_rows: tuple[dict[str, str], ...],
    integrated_rows: tuple[dict[str, str], ...],
    ligand_round2_rows: tuple[dict[str, str], ...],
    finalist_rows: tuple[dict[str, str], ...],
    output_path: Path,
) -> None:
    plt = _import_pyplot()
    mutable_position_count = sum(
        1
        for row in design_mask_rows
        if str(row["mutable_second_sphere"]).strip() == "True" or str(row["mutable_interface"]).strip() == "True"
    )
    funnel_steps = (
        ("Mask positions evaluated", len(design_mask_rows), "#7a8c99"),
        ("Mutable positions retained", mutable_position_count, "#4f7f3b"),
        ("Round-1 ProteinMPNN shortlist", len(protein_round1_rows), "#2f6f8f"),
        ("Round-1 integrated candidates", len(integrated_rows), "#8c2f39"),
        ("Round-2 LigandMPNN shortlist", len(ligand_round2_rows), "#c97f2c"),
        ("Round-2 finalists", len(finalist_rows), "#5b3f8c"),
    )

    fig, axis = plt.subplots(figsize=(12, 6.5), constrained_layout=True)
    max_count = max(step[1] for step in funnel_steps)
    y_values = list(range(len(funnel_steps)))
    for y_value, (label, count, color) in zip(y_values, funnel_steps):
        axis.barh(y_value, count, left=-(count / 2), height=0.78, color=color, edgecolor="none")
        axis.text(0, y_value, str(count), ha="center", va="center", color="white", fontweight="bold")
        axis.text((max_count / 2) + 2, y_value, label, va="center", fontsize=10)
    axis.set_xlim(-(max_count / 2) - 4, (max_count / 2) + 26)
    axis.set_ylim(-0.75, len(funnel_steps) - 0.25)
    axis.invert_yaxis()
    axis.axis("off")
    axis.set_title("Final Figure 2. Design funnel from masks to round-2 shortlist", fontsize=14, fontweight="bold")
    axis.text(
        0.02,
        0.03,
        (
            "Round-2 local redesign expanded three seed scaffolds between the integrated round-1 stage and the "
            "reduced round-2 shortlist; this reporting figure tracks the deterministic narrowing checkpoints."
        ),
        transform=axis.transAxes,
        fontsize=9,
    )
    _save_matplotlib_figure(fig, output_path)
    plt.close(fig)


def render_final_figure_3_multimetal_validation(
    *,
    pocket_rows: tuple[dict[str, str], ...],
    round2_summary_rows: tuple[dict[str, str], ...],
    output_path: Path,
) -> None:
    plt = _import_pyplot()
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    summary_by_key = {(row["candidate_id"], row["target_metal"]): row for row in round2_summary_rows}
    candidate_ids = [row["candidate_id"] for row in pocket_rows]
    status_codes: list[list[int]] = []
    distance_text: list[list[str]] = []
    for candidate_id in candidate_ids:
        row_codes: list[int] = []
        row_text: list[str] = []
        for metal in METAL_ORDER:
            summary_row = summary_by_key[(candidate_id, metal)]
            status = summary_row["screening_status"]
            row_codes.append(0 if status == "stable_bound" else 1 if status == "persistent_capture_flag" else 2)
            row_text.append(f"{'SB' if status == 'stable_bound' else 'PC'}\n{_parse_float(summary_row['median_mean_min_metal_oxygen_distance_A']):.2f}A")
        status_codes.append(row_codes)
        distance_text.append(row_text)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
    cmap = ListedColormap(["#6baed6", "#f28e2b", "#bdbdbd"])
    image = axes[0].imshow(status_codes, cmap=cmap, aspect="auto", vmin=0, vmax=2)
    del image
    axes[0].set_xticks(range(len(METAL_ORDER)), labels=METAL_ORDER)
    axes[0].set_yticks(range(len(candidate_ids)), labels=[_candidate_display_label(candidate_id) for candidate_id in candidate_ids])
    axes[0].set_title("Reduced OpenMM rescreen status and median metal-O distance")
    for row_index, row_values in enumerate(distance_text):
        for column_index, text in enumerate(row_values):
            axes[0].text(column_index, row_index, text, ha="center", va="center", fontsize=8)
    axes[0].legend(
        handles=(
            Patch(color="#6baed6", label="stable_bound"),
            Patch(color="#f28e2b", label="persistent_capture_flag"),
        ),
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.08),
        ncols=2,
    )

    rare_earth_means = [
        _mean_distance_by_candidate(summary_by_key, candidate_id, ("Dy", "Nd", "Y"))
        for candidate_id in candidate_ids
    ]
    competitor_means = [
        _mean_distance_by_candidate(summary_by_key, candidate_id, ("Al", "Fe"))
        for candidate_id in candidate_ids
    ]
    x_values = list(range(len(candidate_ids)))
    width = 0.34
    axes[1].bar(
        [value - width / 2 for value in x_values],
        rare_earth_means,
        width=width,
        color="#1f77b4",
        label="Dy/Nd/Y mean distance",
    )
    axes[1].bar(
        [value + width / 2 for value in x_values],
        competitor_means,
        width=width,
        color="#d95f02",
        label="Al/Fe mean distance",
    )
    axes[1].set_xticks(x_values, labels=[_candidate_display_label(candidate_id) for candidate_id in candidate_ids], rotation=15)
    axes[1].set_ylabel("Median mean minimum metal-oxygen distance (A)")
    axes[1].set_title("Rare-earth retention stays wider than persistent Al/Fe capture")
    axes[1].grid(axis="y", alpha=0.2)
    axes[1].legend(frameon=False)

    fig.suptitle("Final Figure 3. Multi-metal validation comparison", fontsize=14, fontweight="bold")
    _save_matplotlib_figure(fig, output_path)
    plt.close(fig)


def render_final_figure_4_quantum_pocket_summary(
    *,
    pocket_rows: tuple[dict[str, str], ...],
    output_path: Path,
) -> None:
    plt = _import_pyplot()
    candidate_ids = [row["candidate_id"] for row in pocket_rows]
    delta_g_values = [_parse_float(row["deltaG_pocket_proxy_kcal_mol"]) for row in pocket_rows]
    overall_confidences = [_parse_float(row["representative_overall_confidence"]) for row in pocket_rows]
    ligand_confidences = [_parse_float(row["representative_ligand_confidence"]) for row in pocket_rows]
    pocket_confidences = [_parse_float(row["pocket_confidence"]) for row in pocket_rows]
    colors = ["#8c2f39" if index == 0 else "#9fb6c9" for index in range(len(candidate_ids))]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
    y_values = list(range(len(candidate_ids)))
    axes[0].barh(y_values, delta_g_values, color=colors, edgecolor="none")
    axes[0].set_yticks(y_values, labels=[_candidate_display_label(candidate_id) for candidate_id in candidate_ids])
    axes[0].invert_yaxis()
    axes[0].set_xlabel("deltaG_pocket_proxy_kcal_mol (lower is better)")
    axes[0].set_title("Deterministic thermodynamic pocket proxy ranking")
    axes[0].grid(axis="x", alpha=0.2)
    for y_value, delta_g, pocket_confidence in zip(y_values, delta_g_values, pocket_confidences):
        axes[0].text(delta_g + 0.02, y_value, f"PocketConfidence={pocket_confidence:.3f}", va="center", fontsize=9)

    sizes = [400 * value + 50 for value in pocket_confidences]
    axes[1].scatter(overall_confidences, ligand_confidences, s=sizes, c=colors, alpha=0.9)
    for candidate_id, overall_confidence, ligand_confidence in zip(candidate_ids, overall_confidences, ligand_confidences):
        axes[1].text(overall_confidence + 0.01, ligand_confidence + 0.01, _candidate_display_label(candidate_id), fontsize=9)
    axes[1].set_xlim(0.0, 1.05)
    axes[1].set_ylim(0.0, 1.05)
    axes[1].set_xlabel("Representative overall confidence")
    axes[1].set_ylabel("Representative ligand confidence")
    axes[1].set_title("Pocket-confidence components")
    axes[1].grid(alpha=0.2)
    axes[1].text(
        0.03,
        0.04,
        (
            "PocketConfidence = 0.7 x ligand + 0.3 x overall\n"
            f"Main Quantum Pocket = {_candidate_display_label(candidate_ids[0])}"
        ),
        transform=axes[1].transAxes,
        fontsize=9,
    )

    fig.suptitle("Final Figure 4. Main quantum pocket and specific pocket ranking summary", fontsize=14, fontweight="bold")
    _save_matplotlib_figure(fig, output_path)
    plt.close(fig)


def _final_findings_lines(
    *,
    main_pocket_row: dict[str, str],
) -> tuple[str, ...]:
    return (
        "Rare-earth retention was preserved across tested candidates: Dy, Nd, and Y remained stable_bound for all three finalists in the reduced rescreen.",
        "Persistent Al/Fe inner-sphere capture remained unresolved across all three finalists; no candidate solved the Al/Fe problem in the validated pipeline.",
        "No candidate advanced to full QM/QCT in the validated pipeline.",
        (
            "The Main Quantum Pocket ranking is a deterministic next-step prioritization based on completed workflow outputs. "
            f"Under that proxy, `{main_pocket_row['candidate_id']}` is the Main Quantum Pocket with "
            f"deltaG_pocket_proxy_kcal_mol={main_pocket_row['deltaG_pocket_proxy_kcal_mol']}."
        ),
    )


def _render_readme(
    *,
    inputs: FinalProjectReportInputs,
    outputs: FinalProjectReportOutputs,
    table_1_rows: tuple[dict[str, str], ...],
    pocket_rows: tuple[dict[str, str], ...],
) -> str:
    main_pocket_row = pocket_rows[0]
    workflow_outputs = (
        inputs.template_chain_summary_path,
        inputs.design_mask_candidates_path,
        inputs.proteinmpnn_shortlist_path,
        inputs.ligandmpnn_shortlist_path,
        inputs.integrated_candidate_ranking_path,
        inputs.proteinmpnn_round2_shortlist_path,
        inputs.ligandmpnn_round2_smoke_summary_path,
        inputs.ligandmpnn_round2_shortlist_path,
        inputs.rosetta_round2_candidate_ranking_path,
        inputs.round2_md_rescreen_panel_path,
        inputs.round2_openmm_screening_summary_path,
        inputs.round2_openmm_screening_panel_status_path,
    )
    final_outputs = (
        outputs.final_project_report_path,
        outputs.project_summary_report_path,
        outputs.project_detailed_report_path,
        outputs.paper_short_path,
        outputs.paper_long_path,
        outputs.final_table_1_path,
        outputs.final_table_2_path,
        outputs.final_table_3_path,
        outputs.final_figure_1_path,
        outputs.final_figure_2_path,
        outputs.final_figure_3_path,
        outputs.final_figure_4_path,
    )
    lines = [
        "# lanm-dy-selectivity",
        "",
        "## Project Overview",
        "",
        (
            "This repository now packages a complete, reproducible, goal-oriented computational workflow for "
            "lanmodulin Dy selectivity analysis through the round-2 reduced OpenMM rescreen. The final reporting "
            "phase is deterministic and reporting-only: it consumes completed workflow outputs and rewrites the "
            "submission-ready tables, figures, reports, papers, and README without rerunning ProteinMPNN, "
            "LigandMPNN, Rosetta, MD, metadynamics, QM, or quantum steps."
        ),
        "",
        "## Environments Used",
        "",
        "- Repository reporting environment: `python>=3.10`, `pytest`, `matplotlib`, CSV/Markdown generation utilities.",
        "- ProteinMPNN generation environment used earlier in the completed workflow.",
        "- LigandMPNN redesign environment used earlier in the completed workflow.",
        "- Rosetta scoring environment used earlier in the completed workflow.",
        "- OpenMM build/equilibration/screening environment used earlier in the completed workflow.",
        "",
        "## Workflow Phases Completed",
        "",
        "- Template harmonization, residue-role mapping, and design-mask generation across the AM1/Mex and Hans families.",
        "- ProteinMPNN and LigandMPNN shortlist generation for round 1, followed by integrated ranking and MD panel selection.",
        "- First-round OpenMM screening across rare-earth and competitor metals.",
        "- Round-2 local redesign, LigandMPNN/Rosetta rescreening, and reduced OpenMM finalist rescreen.",
        "- Final deterministic reporting and package assembly.",
        "",
        "## Exact Major Outputs",
        "",
        "Completed workflow inputs consumed by the final packaging command:",
        *[f"- `{_display_path(path)}`" for path in workflow_outputs],
        "",
        "Submission-ready outputs generated by the final packaging command:",
        *[f"- `{_display_path(path)}`" for path in final_outputs],
        "",
        "## Final Findings",
        "",
        *[f"- {line}" for line in _final_findings_lines(main_pocket_row=main_pocket_row)],
        "",
        f"- Main Quantum Pocket: `{main_pocket_row['candidate_id']}`.",
        "- Dominant remaining bottleneck: persistent Al/Fe inner-sphere capture during the reduced rescreen.",
        "",
        "## Template Families Covered",
        "",
        _render_markdown_table(
            ("template_id", "family", "topology_class", "chain_ids", "metal_identity", "metal_site_count"),
            tuple(
                (
                    row["template_id"],
                    row["family"],
                    row["topology_class"],
                    row["chain_ids"],
                    row["metal_identity"],
                    row["metal_site_count"],
                )
                for row in table_1_rows
            ),
        ),
        "",
        "## Reproduce The Main Analyses",
        "",
        "```bash",
        "pytest -q",
        "python -m lanm.cli.build_final_project_report",
        "```",
        "",
        (
            "The second command is reporting-only and rebuilds the final submission package from existing workflow "
            "artifacts already present under `results/` and `docs/`."
        ),
    ]
    return "\n".join(lines) + "\n"


def _render_final_project_report(
    *,
    outputs: FinalProjectReportOutputs,
    table_1_rows: tuple[dict[str, str], ...],
    table_2_rows: tuple[dict[str, str], ...],
    pocket_rows: tuple[dict[str, str], ...],
    design_mask_rows: tuple[dict[str, str], ...],
    protein_round1_rows: tuple[dict[str, str], ...],
    integrated_rows: tuple[dict[str, str], ...],
    ligand_round2_rows: tuple[dict[str, str], ...],
    finalist_rows: tuple[dict[str, str], ...],
) -> str:
    main_pocket_row = pocket_rows[0]
    mutable_position_count = sum(
        1
        for row in design_mask_rows
        if str(row["mutable_second_sphere"]).strip() == "True" or str(row["mutable_interface"]).strip() == "True"
    )
    lines = [
        "# Final Project Report",
        "",
        (
            "This report assembles the submission-ready package for the completed lanmodulin Dy-selectivity workflow. "
            "It is a reporting-only deliverable built from finished computational artifacts."
        ),
        "",
        "## Executive Interpretation",
        "",
        COMMON_INTERPRETATION_PARAGRAPH,
        "",
        *[f"- {line}" for line in _final_findings_lines(main_pocket_row=main_pocket_row)],
        "",
        "## Completed Workflow Summary",
        "",
        f"- templates harmonized: `{len(table_1_rows)}`",
        f"- canonical mask positions evaluated: `{len(design_mask_rows)}`",
        f"- mutable positions retained for design: `{mutable_position_count}`",
        f"- round-1 ProteinMPNN shortlist candidates: `{len(protein_round1_rows)}`",
        f"- round-1 integrated candidates with cross-stage evidence: `{len(integrated_rows)}`",
        f"- round-2 LigandMPNN shortlist candidates: `{len(ligand_round2_rows)}`",
        f"- reduced round-2 finalists: `{len(finalist_rows)}`",
        "",
        "## Main Quantum Pocket",
        "",
        (
            f"`{main_pocket_row['candidate_id']}` is the Main Quantum Pocket because it has the best deterministic "
            f"thermodynamic pocket proxy rank (PocketConfidence={main_pocket_row['pocket_confidence']}, "
            f"deltaG_pocket_proxy_kcal_mol={main_pocket_row['deltaG_pocket_proxy_kcal_mol']})."
        ),
        "",
        _render_markdown_table(
            (
                "pocket_rank",
                "candidate_id",
                "topology_class",
                "pocket_confidence",
                "deltaG_pocket_proxy_kcal_mol",
                "Al_status",
                "Fe_status",
            ),
            tuple(
                (
                    row["pocket_rank"],
                    row["candidate_id"],
                    row["topology_class"],
                    row["pocket_confidence"],
                    row["deltaG_pocket_proxy_kcal_mol"],
                    row["Al_status"],
                    row["Fe_status"],
                )
                for row in pocket_rows
            ),
        ),
        "",
        "## Dominant Remaining Bottleneck",
        "",
        (
            "The dominant unresolved bottleneck is persistent Al/Fe inner-sphere capture. The reduced round-2 "
            "OpenMM rescreen preserved rare-earth retention for all finalists, but it did not eliminate Al/Fe "
            "capture for any finalist."
        ),
        "",
        "## Submission Package Inventory",
        "",
        "### Figures",
        "",
        f"- `{_display_path(outputs.final_figure_1_path)}`: template harmonization and residue-role summary",
        f"- `{_display_path(outputs.final_figure_2_path)}`: design funnel from masks to round-2 shortlist",
        f"- `{_display_path(outputs.final_figure_3_path)}`: multi-metal validation comparison",
        f"- `{_display_path(outputs.final_figure_4_path)}`: main quantum pocket and thermodynamic proxy ranking",
        "",
        "### Tables",
        "",
        f"- `{_display_path(outputs.final_table_1_path)}`: template family and topology summary",
        f"- `{_display_path(outputs.final_table_2_path)}`: candidate progression across round 1 and round 2",
        f"- `{_display_path(outputs.final_table_3_path)}`: deterministic specific-pocket thermodynamic proxy ranking",
        "",
        "### Narrative Outputs",
        "",
        f"- `{_display_path(outputs.project_summary_report_path)}`",
        f"- `{_display_path(outputs.project_detailed_report_path)}`",
        f"- `{_display_path(outputs.paper_short_path)}`",
        f"- `{_display_path(outputs.paper_long_path)}`",
        "",
        "## Reproduction",
        "",
        "```bash",
        "pytest -q",
        "python -m lanm.cli.build_final_project_report",
        "```",
    ]
    del table_2_rows
    return "\n".join(lines) + "\n"


def _render_project_summary_report(
    *,
    pocket_rows: tuple[dict[str, str], ...],
) -> str:
    main_pocket_row = pocket_rows[0]
    lines = [
        "# Project Summary Report",
        "",
        (
            "The completed workflow identified a best current quantum-priority pocket and a clear remaining "
            "bottleneck without overclaiming a solved selectivity outcome."
        ),
        "",
        "## Bottom Line",
        "",
        COMMON_INTERPRETATION_PARAGRAPH,
        "",
        f"The Main Quantum Pocket is `{main_pocket_row['candidate_id']}`.",
        "",
        "## Pocket Ranking",
        "",
        _render_markdown_table(
            (
                "rank",
                "candidate_id",
                "pocket_confidence",
                "deltaG_pocket_proxy_kcal_mol",
                "qm_priority_note",
            ),
            tuple(
                (
                    row["pocket_rank"],
                    row["candidate_id"],
                    row["pocket_confidence"],
                    row["deltaG_pocket_proxy_kcal_mol"],
                    row["qm_priority_note"],
                )
                for row in pocket_rows
            ),
        ),
        "",
        "## Interpretation",
        "",
        "- Rare-earth retention was preserved across all three finalists during the reduced rescreen.",
        "- Persistent Al/Fe inner-sphere capture remained unresolved across all three finalists.",
        "- No finalist advanced to full QM/QCT in the validated pipeline.",
        "- The pocket ranking is a deterministic next-step prioritization, not a completed quantum result.",
    ]
    return "\n".join(lines) + "\n"


def _render_project_detailed_report(
    *,
    table_1_rows: tuple[dict[str, str], ...],
    table_2_rows: tuple[dict[str, str], ...],
    pocket_rows: tuple[dict[str, str], ...],
    round2_summary_rows: tuple[dict[str, str], ...],
) -> str:
    main_pocket_row = pocket_rows[0]
    summary_by_key = {(row["candidate_id"], row["target_metal"]): row for row in round2_summary_rows}
    lines = [
        "# Project Detailed Report",
        "",
        "## Scope",
        "",
        (
            "This detailed report documents the completed end-to-end computational workflow through the round-2 "
            "reduced OpenMM rescreen and explains how the final submission package was assembled deterministically."
        ),
        "",
        "## Interpretation Guardrails",
        "",
        COMMON_INTERPRETATION_PARAGRAPH,
        "",
        "## Template Basis",
        "",
        _render_markdown_table(
            ("template_id", "family", "topology_class", "source_type", "chain_ids", "metal_identity", "metal_site_count"),
            tuple(
                (
                    row["template_id"],
                    row["family"],
                    row["topology_class"],
                    row["source_type"],
                    row["chain_ids"],
                    row["metal_identity"],
                    row["metal_site_count"],
                )
                for row in table_1_rows
            ),
        ),
        "",
        "## Candidate Progression",
        "",
        _render_markdown_table(
            (
                "candidate_id",
                "phase_origin",
                "backbone_id",
                "topology_class",
                "key_mutations",
                "retained_or_rejected",
            ),
            tuple(
                (
                    row["candidate_id"],
                    row["phase_origin"],
                    row["backbone_id"],
                    row["topology_class"],
                    row["key_mutations"],
                    row["retained_or_rejected"],
                )
                for row in table_2_rows
            ),
        ),
        "",
        "## Reduced Multi-Metal Rescreen",
        "",
    ]
    for pocket_row in pocket_rows:
        candidate_id = pocket_row["candidate_id"]
        rare_earth_distances = [
            _parse_float(summary_by_key[(candidate_id, metal)]["median_mean_min_metal_oxygen_distance_A"])
            for metal in ("Dy", "Nd", "Y")
        ]
        competitor_distances = [
            _parse_float(summary_by_key[(candidate_id, metal)]["median_mean_min_metal_oxygen_distance_A"])
            for metal in ("Al", "Fe")
        ]
        lines.append(
            (
                f"- `{candidate_id}`: Dy/Nd/Y mean minimum metal-oxygen distance={mean(rare_earth_distances):.3f} A, "
                f"Al/Fe mean minimum metal-oxygen distance={mean(competitor_distances):.3f} A, "
                f"proxy deltaG={pocket_row['deltaG_pocket_proxy_kcal_mol']}."
            )
        )
    lines.extend(
        [
            "",
            "## Main Quantum Pocket Ranking Method",
            "",
            "For each finalist, the reporting package computes:",
            "",
            "- `PocketConfidence = 0.7 * representative_ligand_confidence + 0.3 * representative_overall_confidence`",
            "- `deltaG_pocket_proxy_kcal_mol = -0.592 * ln(max(PocketConfidence, 1e-6))`",
            "- finalists ranked ascending by `deltaG_pocket_proxy_kcal_mol`",
            "",
            (
                f"Using this deterministic proxy, `{main_pocket_row['candidate_id']}` is the Main Quantum Pocket. "
                "This is a prioritization rule for the next step, not evidence that full QM/QCT was completed."
            ),
            "",
            "## Remaining Bottleneck and Next-Step Framing",
            "",
            (
                "The workflow successfully narrowed the search space to a best current quantum-priority pocket while "
                "also showing that persistent Al/Fe inner-sphere capture remains the dominant unresolved bottleneck. "
                "Accordingly, any future QM/QCT effort should be framed as targeted follow-up on the Main Quantum "
                "Pocket after the unresolved competitor-metal capture problem is addressed."
            ),
        ]
    )
    return "\n".join(lines) + "\n"


def _render_short_paper(
    *,
    pocket_rows: tuple[dict[str, str], ...],
) -> str:
    main_pocket_row = pocket_rows[0]
    lines = [
        "# Short Paper: Deterministic Final Report for Lanmodulin Dy-Selectivity Workflow",
        "",
        "## Abstract",
        "",
        (
            "We assembled a submission-ready reporting package for a completed, reproducible computational workflow "
            "targeting dysprosium selectivity in lanmodulin scaffolds. The workflow combined template harmonization, "
            "design-mask construction, ProteinMPNN and LigandMPNN triage, Rosetta score-based down-selection, and "
            "first-round plus round-2 reduced OpenMM multi-metal rescreening. Rare-earth retention was preserved "
            "across tested candidates. Persistent Al/Fe inner-sphere capture remained unresolved. No candidate "
            "advanced to full QM/QCT in the validated pipeline. The Main Quantum Pocket ranking is a deterministic "
            "next-step prioritization based on completed workflow outputs. Using a pocket thermodynamic proxy "
            "defined as `PocketConfidence = 0.7 x ligand_confidence + 0.3 x overall_confidence` and "
            "`deltaG_pocket_proxy_kcal_mol = -0.592 x ln(max(PocketConfidence, 1e-6))`, the best current "
            f"quantum-priority pocket is `{main_pocket_row['candidate_id']}`."
        ),
        "",
        "## Results",
        "",
        "- Template harmonization unified four structural templates across AM1/Mex and Hans families.",
        "- Deterministic mask construction reduced the design space to a focused mutable set for round-1 and round-2 redesign.",
        "- Three round-2 finalists completed the reduced OpenMM rescreen and all preserved Dy/Nd/Y retention.",
        "- All three finalists retained persistent Al/Fe capture flags, so the Al/Fe problem remains unsolved in the validated pipeline.",
        "",
        "## Pocket Prioritization",
        "",
        _render_markdown_table(
            ("rank", "candidate_id", "pocket_confidence", "deltaG_pocket_proxy_kcal_mol"),
            tuple(
                (
                    row["pocket_rank"],
                    row["candidate_id"],
                    row["pocket_confidence"],
                    row["deltaG_pocket_proxy_kcal_mol"],
                )
                for row in pocket_rows
            ),
        ),
        "",
        "## Conclusion",
        "",
        (
            "The completed workflow should be framed as a successful reproducible prioritization campaign rather "
            "than as a solved selectivity campaign. It identified the best current quantum-priority pocket and "
            "made the dominant remaining bottleneck explicit: persistent Al/Fe inner-sphere capture."
        ),
    ]
    return "\n".join(lines) + "\n"


def _render_long_paper(
    *,
    table_1_rows: tuple[dict[str, str], ...],
    pocket_rows: tuple[dict[str, str], ...],
) -> str:
    main_pocket_row = pocket_rows[0]
    lines = [
        "# Long Paper: Reproducible Computational Workflow for Lanmodulin Dy Selectivity Prioritization",
        "",
        "## Abstract",
        "",
        (
            "This work consolidates a completed computational workflow for lanmodulin dysprosium selectivity into a "
            "submission-ready package. The pipeline harmonized four experimental templates, derived deterministic "
            "residue-role and design-mask annotations, executed staged ProteinMPNN and LigandMPNN redesign, applied "
            "topology-aware Rosetta down-selection, and completed first-round plus round-2 reduced OpenMM multi-metal "
            "rescreening. Rare-earth retention was preserved across tested candidates. Persistent Al/Fe inner-sphere "
            "capture remained unresolved. No candidate advanced to full QM/QCT in the validated pipeline. The Main "
            "Quantum Pocket ranking is a deterministic next-step prioritization based on completed workflow outputs. "
            f"Under that rule, `{main_pocket_row['candidate_id']}` is the best current quantum-priority pocket."
        ),
        "",
        "## Introduction",
        "",
        (
            "Lanmodulin already provides a strong structural basis for rare-earth recognition, but an end-to-end "
            "Dy-prioritization workflow still needs to show two things simultaneously: preservation of rare-earth "
            "binding across viable designs and suppression of undesirable competitor-metal capture. This project was "
            "therefore framed as a goal-oriented workflow rather than a single scoring step. The objective was not "
            "to overclaim success, but to build a reproducible computational path that could identify the best "
            "current pocket for future quantum follow-up while exposing the dominant remaining failure mode."
        ),
        "",
        "## Methods",
        "",
        (
            "The workflow began from four structural templates spanning AM1/Mex and Hans lanmodulin families. "
            "Template harmonization aligned residue numbering, summarized first-shell and second-sphere roles, and "
            "generated deterministic design masks. ProteinMPNN and LigandMPNN were then used in staged redesign "
            "rounds, with Rosetta providing topology-aware score-based triage. First-round OpenMM screening supplied "
            "multi-metal validation feedback for retained designs, and a focused round-2 redesign/rescreen branch "
            "completed the validated computational path. The reporting-only final phase used existing outputs only. "
            "For the three round-2 finalists, a deterministic pocket proxy was computed with "
            "`PocketConfidence = 0.7 x representative_ligand_confidence + 0.3 x representative_overall_confidence`, "
            "followed by `deltaG_pocket_proxy_kcal_mol = -0.592 x ln(max(PocketConfidence, 1e-6))`. Finalists were "
            "ranked by ascending proxy deltaG, and the top-ranked finalist was designated the Main Quantum Pocket."
        ),
        "",
        "## Results",
        "",
        (
            "Template harmonization covered the AM1/Mex monomer references 6MI5 and 8FNS together with the Hans "
            "multichain references 8DQ2 and 8FNR."
        ),
        "",
        _render_markdown_table(
            ("template_id", "family", "topology_class", "chain_ids", "metal_identity", "metal_site_count"),
            tuple(
                (
                    row["template_id"],
                    row["family"],
                    row["topology_class"],
                    row["chain_ids"],
                    row["metal_identity"],
                    row["metal_site_count"],
                )
                for row in table_1_rows
            ),
        ),
        "",
        (
            "Across the reduced round-2 rescreen, all three finalists preserved rare-earth retention, but all three "
            "also retained persistent Al/Fe capture flags. The deterministic pocket proxy therefore serves as a "
            "next-step prioritization among unfinished candidates rather than as proof that the selectivity problem "
            "has been solved."
        ),
        "",
        _render_markdown_table(
            ("rank", "candidate_id", "pocket_confidence", "deltaG_pocket_proxy_kcal_mol", "qm_priority_note"),
            tuple(
                (
                    row["pocket_rank"],
                    row["candidate_id"],
                    row["pocket_confidence"],
                    row["deltaG_pocket_proxy_kcal_mol"],
                    row["qm_priority_note"],
                )
                for row in pocket_rows
            ),
        ),
        "",
        "## Discussion",
        "",
        (
            "The completed workflow achieved the key process objective: it reduced a broad template-and-design space "
            "to a small finalist set, then identified the best current pocket for future high-cost follow-up. "
            "However, the biological and thermodynamic objective remains incomplete. Persistent Al/Fe inner-sphere "
            "capture remained unresolved across the finalist set, so no design can be presented as a solved "
            "Dy-selective outcome. Likewise, no finalist advanced to full QM/QCT in the validated pipeline. The "
            "reporting package therefore distinguishes clearly between completed evidence and next-step inference."
        ),
        "",
        "## Conclusion",
        "",
        (
            f"The workflow identified `{main_pocket_row['candidate_id']}` as the Main Quantum Pocket under a "
            "deterministic thermodynamic proxy derived from completed workflow outputs. Rare-earth retention was "
            "preserved across tested candidates, but persistent Al/Fe inner-sphere capture remained unresolved. "
            "Accordingly, the project should be viewed as a complete, reproducible prioritization workflow that "
            "identified both the best current quantum-priority pocket and the dominant remaining bottleneck."
        ),
    ]
    return "\n".join(lines) + "\n"


def build_final_project_package(
    *,
    inputs: FinalProjectReportInputs,
    outputs: FinalProjectReportOutputs,
) -> FinalProjectBuildResult:
    for path in (
        inputs.template_chain_summary_path,
        inputs.residue_role_map_path,
        inputs.design_mask_candidates_path,
        inputs.proteinmpnn_shortlist_path,
        inputs.ligandmpnn_shortlist_path,
        inputs.integrated_candidate_ranking_path,
        inputs.md_validation_panel_path,
        inputs.proteinmpnn_round2_shortlist_path,
        inputs.ligandmpnn_round2_smoke_summary_path,
        inputs.ligandmpnn_round2_shortlist_path,
        inputs.rosetta_round2_candidate_ranking_path,
        inputs.round2_md_rescreen_panel_path,
        inputs.round2_openmm_screening_summary_path,
        inputs.round2_openmm_screening_panel_status_path,
    ):
        require_path(path)

    template_chain_rows = _load_csv_rows(inputs.template_chain_summary_path)
    residue_role_rows = _load_csv_rows(inputs.residue_role_map_path)
    design_mask_rows = _load_csv_rows(inputs.design_mask_candidates_path)
    protein_round1_rows = _load_csv_rows(inputs.proteinmpnn_shortlist_path)
    ligand_round1_rows = _load_csv_rows(inputs.ligandmpnn_shortlist_path)
    integrated_rows = _load_csv_rows(inputs.integrated_candidate_ranking_path)
    md_panel_rows = _load_csv_rows(inputs.md_validation_panel_path)
    protein_round2_rows = _load_csv_rows(inputs.proteinmpnn_round2_shortlist_path)
    ligand_round2_rows = _load_csv_rows(inputs.ligandmpnn_round2_shortlist_path)
    finalist_rows = _load_csv_rows(inputs.round2_md_rescreen_panel_path)
    round2_summary_rows = _load_csv_rows(inputs.round2_openmm_screening_summary_path)
    panel_status_rows = _load_csv_rows(inputs.round2_openmm_screening_panel_status_path)

    table_1_rows = build_final_template_summary_table(template_chain_rows)
    table_2_rows = build_candidate_progression_table(
        protein_round1_rows=protein_round1_rows,
        ligand_round1_rows=ligand_round1_rows,
        integrated_rows=integrated_rows,
        md_panel_rows=md_panel_rows,
        protein_round2_rows=protein_round2_rows,
        ligand_round2_rows=ligand_round2_rows,
        finalist_rows=finalist_rows,
    )
    pocket_rows = build_specific_pocket_thermodynamics_table(finalist_rows, panel_status_rows)
    main_quantum_pocket_candidate_id = pocket_rows[0]["candidate_id"]

    write_csv_rows(outputs.final_table_1_path, table_1_rows)
    write_csv_rows(outputs.final_table_2_path, table_2_rows)
    write_csv_rows(outputs.final_table_3_path, pocket_rows)

    render_final_figure_1_template_harmonization(
        table_1_rows=table_1_rows,
        residue_role_rows=residue_role_rows,
        output_path=outputs.final_figure_1_path,
    )
    render_final_figure_2_design_funnel(
        design_mask_rows=design_mask_rows,
        protein_round1_rows=protein_round1_rows,
        integrated_rows=integrated_rows,
        ligand_round2_rows=ligand_round2_rows,
        finalist_rows=finalist_rows,
        output_path=outputs.final_figure_2_path,
    )
    render_final_figure_3_multimetal_validation(
        pocket_rows=pocket_rows,
        round2_summary_rows=round2_summary_rows,
        output_path=outputs.final_figure_3_path,
    )
    render_final_figure_4_quantum_pocket_summary(
        pocket_rows=pocket_rows,
        output_path=outputs.final_figure_4_path,
    )

    atomic_write_text(
        outputs.final_project_report_path,
        _render_final_project_report(
            outputs=outputs,
            table_1_rows=table_1_rows,
            table_2_rows=table_2_rows,
            pocket_rows=pocket_rows,
            design_mask_rows=design_mask_rows,
            protein_round1_rows=protein_round1_rows,
            integrated_rows=integrated_rows,
            ligand_round2_rows=ligand_round2_rows,
            finalist_rows=finalist_rows,
        ),
    )
    atomic_write_text(
        outputs.project_summary_report_path,
        _render_project_summary_report(pocket_rows=pocket_rows),
    )
    atomic_write_text(
        outputs.project_detailed_report_path,
        _render_project_detailed_report(
            table_1_rows=table_1_rows,
            table_2_rows=table_2_rows,
            pocket_rows=pocket_rows,
            round2_summary_rows=round2_summary_rows,
        ),
    )
    atomic_write_text(outputs.paper_short_path, _render_short_paper(pocket_rows=pocket_rows))
    atomic_write_text(outputs.paper_long_path, _render_long_paper(table_1_rows=table_1_rows, pocket_rows=pocket_rows))
    atomic_write_text(
        outputs.readme_path,
        _render_readme(
            inputs=inputs,
            outputs=outputs,
            table_1_rows=table_1_rows,
            pocket_rows=pocket_rows,
        ),
    )

    output_paths = (
        outputs.final_project_report_path,
        outputs.project_summary_report_path,
        outputs.project_detailed_report_path,
        outputs.paper_short_path,
        outputs.paper_long_path,
        outputs.readme_path,
        outputs.final_table_1_path,
        outputs.final_table_2_path,
        outputs.final_table_3_path,
        outputs.final_figure_1_path,
        outputs.final_figure_2_path,
        outputs.final_figure_3_path,
        outputs.final_figure_4_path,
    )
    missing_outputs = [path for path in output_paths if not path.exists()]
    if missing_outputs:
        raise RuntimeError(
            "Final project package is incomplete; missing outputs: "
            + ", ".join(_display_path(path) for path in missing_outputs)
        )

    return FinalProjectBuildResult(
        main_quantum_pocket_candidate_id=main_quantum_pocket_candidate_id,
        pocket_rows=pocket_rows,
        output_paths=output_paths,
    )
