#!/usr/bin/env python3
"""Create exact-data validation plots for the OpenMolcas prototype package."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


FOOTER = "Official OpenMolcas Dy atomic test derivative — prototype only; not a LanM structure or manuscript result"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def save_both(fig: plt.Figure, output_stem: Path) -> None:
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def occupation_plot(table_dir: Path, figure_dir: Path) -> None:
    rows = read_csv(table_dir / "dy3_natural_occupations_all_roots.csv")
    values = np.array(
        [[float(row[f"active_orbital_{index}"]) for index in range(1, 8)] for row in rows]
    )
    fig, ax = plt.subplots(figsize=(9.0, 7.1), constrained_layout=True)
    image = ax.imshow(values, cmap="viridis", vmin=1.0, vmax=1.7, aspect="auto")
    ax.set_title("Dy(III) CAS(9,7) natural occupations by root", fontsize=15, weight="bold", pad=14)
    ax.set_xlabel("Active-orbital index within symmetry 2")
    ax.set_ylabel("RASSCF root")
    ax.set_xticks(np.arange(7), labels=[str(index) for index in range(1, 8)])
    ax.set_yticks(np.arange(11), labels=[str(index) for index in range(1, 12)])
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            color = "white" if values[row, column] < 1.35 else "black"
            ax.text(column, row, f"{values[row, column]:.3f}", ha="center", va="center", fontsize=8.1, color=color)
    colorbar = fig.colorbar(image, ax=ax, shrink=0.82, pad=0.025)
    colorbar.set_label("Natural occupation")
    fig.text(0.5, -0.015, FOOTER, ha="center", fontsize=8.4, color="#555555")
    save_both(fig, figure_dir / "prototype_dy3_natural_occupations")


def doublet_plot(table_dir: Path, figure_dir: Path) -> None:
    rows = read_csv(table_dir / "dy3_kramers_doublet_energies.csv")
    indices = np.array([int(row["kramers_doublet"]) for row in rows])
    energies = np.array([float(row["relative_energy_cm-1"]) for row in rows])
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.7), gridspec_kw={"width_ratios": [1.7, 1.0]}, constrained_layout=True)

    axes[0].vlines(indices, 0, energies, color="#8fb8d8", linewidth=1.0)
    axes[0].scatter(indices, energies, color="#155a8a", s=25, zorder=3)
    axes[0].set_title("All 33 Kramers doublets", fontsize=12.5, weight="bold")
    axes[0].set_xlabel("Kramers-doublet index")
    axes[0].set_ylabel("Relative spin–orbit energy (cm$^{-1}$)")
    axes[0].set_xlim(0.2, 33.8)
    axes[0].grid(axis="y", alpha=0.25)

    low_indices = indices[:8]
    low_energies = energies[:8]
    axes[1].vlines(low_indices, 0, low_energies, color="#d8aa79", linewidth=1.5)
    axes[1].scatter(low_indices, low_energies, color="#a84b13", s=34, zorder=3)
    for index, energy in zip(low_indices, low_energies):
        axes[1].annotate(f"{energy:.1f}", (index, energy), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8)
    axes[1].set_title("Lowest eight doublets", fontsize=12.5, weight="bold")
    axes[1].set_xlabel("Kramers-doublet index")
    axes[1].set_ylabel("Relative energy (cm$^{-1}$)")
    axes[1].set_xticks(low_indices)
    axes[1].set_ylim(0, max(low_energies) * 1.16)
    axes[1].grid(axis="y", alpha=0.25)

    fig.suptitle("Dy(III) spin–orbit energy ladder", fontsize=15, weight="bold")
    fig.text(0.5, -0.02, FOOTER, ha="center", fontsize=8.4, color="#555555")
    save_both(fig, figure_dir / "prototype_dy3_kramers_doublet_energies")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    table_dir = root / "tables"
    figure_dir = root / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    occupation_plot(table_dir, figure_dir)
    doublet_plot(table_dir, figure_dir)


if __name__ == "__main__":
    main()
