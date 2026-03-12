"""Deterministic plotting for Phase 1 summaries."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from lanm.configuration import ProjectConfig
from lanm.models import MetalSiteSummary


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
