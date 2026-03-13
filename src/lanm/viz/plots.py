"""Deterministic plotting helpers for Phase 1 and Phase 2 summaries."""

from __future__ import annotations

from io import BytesIO
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from lanm.filesystem import atomic_write_bytes
from lanm.configuration import ProjectConfig
from lanm.models import (
    CrossTemplateAlignmentRow,
    DesignMaskCandidate,
    MetalSiteSummary,
    ResidueRoleAssignment,
)


def _site_label(summary: MetalSiteSummary) -> str:
    return f"{summary.structure_id}:{summary.chain_id}{summary.residue_seq}"


def _site_color(summary: MetalSiteSummary, config: ProjectConfig) -> str:
    if summary.metal_element == config.target_metal.upper():
        return "#d95f02"
    if summary.metal_element in {item.upper() for item in config.competitors}:
        return "#1b9e77"
    return "#5c677d"


def render_site_overview(
    summaries: list[MetalSiteSummary],
    config: ProjectConfig,
    output_path: Path,
) -> None:
    """Render a compact overview chart for annotated metal sites."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True)
    if not summaries:
        for axis in axes:
            axis.axis("off")
        fig.suptitle("Lanmodulin metal-site overview")
        axes[0].text(0.5, 0.5, "No metal sites were annotated.", ha="center", va="center")
        fig.savefig(output_path, dpi=200)
        plt.close(fig)
        return
    ordered = sorted(summaries, key=lambda item: (item.structure_id, item.chain_id, item.residue_seq, item.atom_serial))
    labels = [_site_label(item) for item in ordered]
    colors = [_site_color(item, config) for item in ordered]
    donor_counts = [item.donor_atom_count for item in ordered]
    residue_counts = [item.residue_within_6a_count for item in ordered]
    axes[0].bar(labels, donor_counts, color=colors)
    axes[0].set_ylabel("Donor atoms <= 3.2 A")
    axes[0].set_title("First-shell donor counts")
    axes[0].tick_params(axis="x", rotation=45, labelsize=8)
    axes[1].bar(labels, residue_counts, color=colors)
    axes[1].set_ylabel("Residues <= 6.0 A")
    axes[1].set_title("Residues in metal environment")
    axes[1].tick_params(axis="x", rotation=45, labelsize=8)
    fig.suptitle(
        f"Lanmodulin site overview: target={config.target_metal}, "
        f"competitors={', '.join(config.competitors)}"
    )
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def render_template_harmonization_overview(
    *,
    cross_template_rows: tuple[CrossTemplateAlignmentRow, ...],
    residue_role_rows: tuple[ResidueRoleAssignment, ...],
    design_mask_rows: tuple[DesignMaskCandidate, ...],
    output_path: Path,
) -> None:
    """Render a deterministic Phase 2C overview of canonical role counts and mask tracks."""
    canonical_positions = sorted(
        {
            row.canonical_family_position
            for row in cross_template_rows
            if row.canonical_family_position is not None
        }
    )
    if not canonical_positions:
        fig, axis = plt.subplots(figsize=(10, 4), constrained_layout=True)
        axis.axis("off")
        axis.text(0.5, 0.5, "No canonical family positions available.", ha="center", va="center")
        buffer = BytesIO()
        fig.savefig(buffer, format="png", dpi=200, metadata={"Software": "lanm-dy-selectivity"})
        plt.close(fig)
        atomic_write_bytes(output_path, buffer.getvalue())
        return

    max_position = canonical_positions[-1]
    first_shell_counts = [0] * max_position
    second_sphere_counts = [0] * max_position
    hans_interface_counts = [0] * max_position
    for row in residue_role_rows:
        if row.canonical_family_position is None:
            continue
        index = row.canonical_family_position - 1
        if row.role == "first_shell":
            first_shell_counts[index] += 1
        elif row.role == "second_sphere":
            second_sphere_counts[index] += 1
        elif row.role == "interchain_contact" and row.template_id in {"8DQ2", "8FNR"}:
            hans_interface_counts[index] += 1

    x_values = list(range(1, max_position + 1))
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(15, 8.5),
        height_ratios=[2.8, 1.9, 0.9],
        constrained_layout=True,
    )
    axes[0].plot(x_values, first_shell_counts, color="#8c2f39", linewidth=1.8, label="first_shell")
    axes[0].plot(x_values, second_sphere_counts, color="#2f6f8f", linewidth=1.8, label="second_sphere")
    axes[0].plot(
        x_values,
        hans_interface_counts,
        color="#4f7f3b",
        linewidth=1.8,
        label="Hans interchain_contact",
    )
    axes[0].set_ylabel("Observation count")
    axes[0].set_title("Canonical role evidence across template harmonization")
    axes[0].legend(loc="upper right", frameon=False, ncols=3)
    axes[0].grid(axis="y", alpha=0.25)

    mask_specs = (
        ("fixed_first_shell", 4, "#8c2f39"),
        ("mutable_second_sphere", 3, "#2f6f8f"),
        ("mutable_interface", 2, "#4f7f3b"),
        ("protected_positions", 1, "#7a5c17"),
    )
    for attribute_name, y_value, color in mask_specs:
        positions = [
            candidate.canonical_family_position
            for candidate in design_mask_rows
            if getattr(candidate, attribute_name)
        ]
        if positions:
            axes[1].scatter(
                positions,
                [y_value] * len(positions),
                marker="s",
                s=42,
                color=color,
                edgecolors="none",
            )
    axes[1].set_yticks([4, 3, 2, 1], labels=[item[0] for item in mask_specs])
    axes[1].set_ylim(0.5, 4.5)
    axes[1].set_ylabel("Phase 2C masks")
    axes[1].grid(axis="x", alpha=0.15)

    fixed_count = sum(1 for row in design_mask_rows if row.fixed_first_shell)
    second_count = sum(1 for row in design_mask_rows if row.mutable_second_sphere)
    interface_count = sum(1 for row in design_mask_rows if row.mutable_interface)
    protected_count = sum(1 for row in design_mask_rows if row.protected_positions)
    his_tag_positions = [
        str(row.template_residue_seq)
        for row in cross_template_rows
        if row.template_id == "6MI5" and row.alignment_status == "outside_am1_mature_reference"
    ]
    axes[2].axis("off")
    axes[2].text(
        0.01,
        0.72,
        (
            f"fixed_first_shell={fixed_count} | mutable_second_sphere={second_count} | "
            f"mutable_interface={interface_count} | protected_positions={protected_count}"
        ),
        fontsize=10,
        fontweight="bold",
        transform=axes[2].transAxes,
    )
    axes[2].text(
        0.01,
        0.28,
        "Noncanonical His-tag residues protected outside the canonical axis: "
        + ", ".join(his_tag_positions),
        fontsize=9,
        transform=axes[2].transAxes,
    )
    axes[1].set_xlim(1, max_position)
    axes[1].set_xlabel("canonical_family_position")
    axes[0].set_xlim(1, max_position)
    fig.suptitle("Template harmonization overview", fontsize=14, fontweight="bold")

    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=200, metadata={"Software": "lanm-dy-selectivity"})
    plt.close(fig)
    atomic_write_bytes(output_path, buffer.getvalue())
