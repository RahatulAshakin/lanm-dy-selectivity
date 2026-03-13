"""Phase 2A template harmonization summaries and report generation."""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from lanm.data.fetch import load_structure_manifest, require_path
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.models import AtomRecord
from lanm.paths import REPO_ROOT
from lanm.structure.atoms import read_atom_records
from lanm.structure.templates import parse_cif_atom_records, read_cif_experimental_method

_METAL_ELEMENTS = {"Y", "ND", "LA", "DY"}


@dataclass(frozen=True, slots=True)
class TemplateInput:
    template_id: str
    source_type: str
    path: Path


@dataclass(frozen=True, slots=True)
class TemplateChainSummaryRow:
    template_id: str
    source_type: str
    chain_id: str
    residue_start: str
    residue_end: str
    residue_range: str
    residue_count: int
    metal_identity: str
    metal_site_count: int
    experimental_method: str


@dataclass(frozen=True, slots=True)
class TemplateSiteSummaryRow:
    template_id: str
    source_type: str
    chain_id: str
    site_index: int
    site_label: str
    residue_name: str
    residue_seq: int
    insertion_code: str
    metal_identity: str
    experimental_method: str


@dataclass(frozen=True, slots=True)
class TemplateSummary:
    template_id: str
    source_type: str
    input_path: str
    chain_ids: tuple[str, ...]
    residue_ranges: tuple[str, ...]
    metal_identity: str
    metal_site_count: int
    experimental_method: str
    chain_rows: tuple[TemplateChainSummaryRow, ...]
    site_rows: tuple[TemplateSiteSummaryRow, ...]


@dataclass(frozen=True, slots=True)
class TemplateHarmonizationArtifacts:
    template_summaries: tuple[TemplateSummary, ...]
    chain_rows: tuple[TemplateChainSummaryRow, ...]
    site_rows: tuple[TemplateSiteSummaryRow, ...]
    report_markdown: str
    sequence_record_count: int


TEMPLATE_INPUTS = (
    TemplateInput("6MI5", "cif", Path("data/raw/public/structures/6MI5.cif")),
    TemplateInput("8FNS", "csv_atom_table", Path("data/raw/local_bundle/8fns_atoms.csv")),
    TemplateInput("8DQ2", "csv_atom_table", Path("data/raw/local_bundle/8dq2_atoms.csv")),
    TemplateInput("8FNR", "cif", Path("data/raw/public/structures/8FNR.cif")),
)

_FAMILY_COMPARISONS = (
    ("AM1/Mex family", ("6MI5", "8FNS")),
    ("Hans family", ("8DQ2", "8FNR")),
)


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def load_sequence_rows(path: Path) -> list[dict[str, str]]:
    """Load the sequence reference table and fail with the exact missing path."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return sorted(
            ({key: str(value) for key, value in row.items()} for row in reader),
            key=lambda row: (
                row.get("organism", ""),
                row.get("strain", ""),
                row.get("gene_symbol", ""),
            ),
        )


def _load_manifest_method_map(manifest_path: Path) -> dict[str, str]:
    return {entry.pdb_id: entry.method for entry in load_structure_manifest(manifest_path)}


def _format_residue_position(residue_seq: int, insertion_code: str) -> str:
    return f"{residue_seq}{insertion_code}" if insertion_code else str(residue_seq)


def _summarize_template(
    template_id: str,
    source_type: str,
    input_path: str,
    experimental_method: str,
    atoms: list[AtomRecord],
) -> TemplateSummary:
    polymer_residues: dict[str, set[tuple[int, str]]] = defaultdict(set)
    metal_atoms = sorted(
        [
            atom for atom in atoms
            if atom.record_type == "HETATM" and atom.element_upper in _METAL_ELEMENTS
        ],
        key=lambda atom: (atom.chain_id, atom.residue_seq, atom.insertion_code, atom.atom_serial),
    )
    for atom in atoms:
        if atom.record_type != "ATOM":
            continue
        polymer_residues[atom.chain_id].add((atom.residue_seq, atom.insertion_code))
    metal_identities = sorted({atom.element_upper for atom in metal_atoms})
    chain_ids = tuple(sorted(polymer_residues))
    chain_site_counts: dict[str, int] = defaultdict(int)
    site_rows: list[TemplateSiteSummaryRow] = []
    for site_index, atom in enumerate(metal_atoms, start=1):
        chain_site_counts[atom.chain_id] += 1
        site_rows.append(
            TemplateSiteSummaryRow(
                template_id=template_id,
                source_type=source_type,
                chain_id=atom.chain_id,
                site_index=site_index,
                site_label=f"{atom.chain_id}:{_format_residue_position(atom.residue_seq, atom.insertion_code)}",
                residue_name=atom.residue_name,
                residue_seq=atom.residue_seq,
                insertion_code=atom.insertion_code,
                metal_identity=atom.element_upper,
                experimental_method=experimental_method,
            )
        )
    chain_rows: list[TemplateChainSummaryRow] = []
    residue_ranges: list[str] = []
    for chain_id in chain_ids:
        ordered_residues = sorted(polymer_residues[chain_id], key=lambda item: (item[0], item[1]))
        start_seq, start_insertion = ordered_residues[0]
        end_seq, end_insertion = ordered_residues[-1]
        residue_start = _format_residue_position(start_seq, start_insertion)
        residue_end = _format_residue_position(end_seq, end_insertion)
        residue_range = f"{chain_id}:{residue_start}-{residue_end}"
        residue_ranges.append(residue_range)
        chain_rows.append(
            TemplateChainSummaryRow(
                template_id=template_id,
                source_type=source_type,
                chain_id=chain_id,
                residue_start=residue_start,
                residue_end=residue_end,
                residue_range=residue_range,
                residue_count=len(ordered_residues),
                metal_identity=", ".join(metal_identities),
                metal_site_count=chain_site_counts.get(chain_id, 0),
                experimental_method=experimental_method,
            )
        )
    return TemplateSummary(
        template_id=template_id,
        source_type=source_type,
        input_path=input_path,
        chain_ids=chain_ids,
        residue_ranges=tuple(residue_ranges),
        metal_identity=", ".join(metal_identities),
        metal_site_count=len(site_rows),
        experimental_method=experimental_method,
        chain_rows=tuple(chain_rows),
        site_rows=tuple(site_rows),
    )


def _load_template_summary(
    template_input: TemplateInput,
    manifest_methods: dict[str, str],
) -> TemplateSummary:
    path = REPO_ROOT / template_input.path
    require_path(path)
    if template_input.source_type == "cif":
        atoms = parse_cif_atom_records(path, structure_id=template_input.template_id)
        experimental_method = read_cif_experimental_method(path)
    else:
        atoms = read_atom_records(path, structure_id=template_input.template_id)
        experimental_method = manifest_methods[template_input.template_id]
    return _summarize_template(
        template_id=template_input.template_id,
        source_type=template_input.source_type,
        input_path=_display_path(path),
        experimental_method=experimental_method,
        atoms=atoms,
    )


def _render_markdown_table(headers: tuple[str, ...], rows: list[tuple[str, ...]]) -> str:
    header_line = "| " + " | ".join(headers) + " |"
    divider_line = "| " + " | ".join(["---"] * len(headers)) + " |"
    body_lines = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join([header_line, divider_line, *body_lines])


def _site_distribution(summary: TemplateSummary) -> str:
    return ", ".join(
        f"{row.chain_id}:{row.metal_site_count}"
        for row in summary.chain_rows
    )


def _am1_sequence_note(sequence_rows: list[dict[str, str]]) -> str:
    for row in sequence_rows:
        if row.get("strain", "").strip().upper() == "AM1":
            return row.get("notes", "").strip()
    return ""


def render_template_harmonization_markdown(
    template_summaries: tuple[TemplateSummary, ...],
    sequence_rows: list[dict[str, str]],
) -> str:
    summary_map = {summary.template_id: summary for summary in template_summaries}
    lines = [
        "# Template Harmonization",
        "",
        "Phase 2A deterministic template parsing summary for four lanmodulin templates.",
        "",
        f"Sequence reference records loaded: {len(sequence_rows)}",
        "",
        "No alignment, residue-role mapping, design masks, or figures are generated in this phase.",
        "",
    ]
    am1_note = _am1_sequence_note(sequence_rows)
    if am1_note:
        lines.extend(
            [
                "AM1 sequence note:",
                "",
                f"- {am1_note}",
                "",
            ]
        )
    for family_label, template_ids in _FAMILY_COMPARISONS:
        rows: list[tuple[str, ...]] = []
        for template_id in template_ids:
            summary = summary_map[template_id]
            rows.append(
                (
                    summary.template_id,
                    summary.source_type,
                    ", ".join(summary.chain_ids),
                    "; ".join(summary.residue_ranges),
                    summary.metal_identity,
                    str(summary.metal_site_count),
                    summary.experimental_method,
                )
            )
        lines.extend(
            [
                f"## {family_label}",
                "",
                _render_markdown_table(
                    (
                        "template_id",
                        "source_type",
                        "chain_ids",
                        "residue_ranges",
                        "metal_identity",
                        "metal_site_count",
                        "experimental_method",
                    ),
                    rows,
                ),
                "",
            ]
        )
        if family_label == "AM1/Mex family":
            lines.extend(
                [
                    f"- `6MI5` uses model 1 of the NMR ensemble for deterministic counting and reports {_site_distribution(summary_map['6MI5'])} metal sites.",
                    f"- `8FNS` is a single-chain crystal template spanning {'; '.join(summary_map['8FNS'].residue_ranges)} with {_site_distribution(summary_map['8FNS'])} metal sites.",
                    "",
                ]
            )
        else:
            lines.extend(
                [
                    f"- `8DQ2` and `8FNR` both resolve chains {', '.join(summary_map['8DQ2'].chain_ids)} over the Hans family backbone.",
                    f"- Metal occupancy differs across the tetramer: `8DQ2` shows {_site_distribution(summary_map['8DQ2'])}, while `8FNR` shows {_site_distribution(summary_map['8FNR'])}.",
                    "- `8FNR` chain D retains residues 24-133 but the observed model omits residues 34-38, reducing the chain residue count relative to `8DQ2`.",
                    "",
                ]
            )
    return "\n".join(lines).rstrip() + "\n"


def build_template_harmonization_artifacts(
    *,
    sequences_path: Path,
    manifest_path: Path,
    template_inputs: tuple[TemplateInput, ...] = TEMPLATE_INPUTS,
) -> TemplateHarmonizationArtifacts:
    sequence_rows = load_sequence_rows(sequences_path)
    manifest_methods = _load_manifest_method_map(manifest_path)
    template_summaries = tuple(
        _load_template_summary(template_input, manifest_methods)
        for template_input in template_inputs
    )
    chain_rows = tuple(
        row
        for summary in template_summaries
        for row in summary.chain_rows
    )
    site_rows = tuple(
        row
        for summary in template_summaries
        for row in summary.site_rows
    )
    return TemplateHarmonizationArtifacts(
        template_summaries=template_summaries,
        chain_rows=chain_rows,
        site_rows=site_rows,
        report_markdown=render_template_harmonization_markdown(template_summaries, sequence_rows),
        sequence_record_count=len(sequence_rows),
    )


def write_template_harmonization_outputs(
    *,
    artifacts: TemplateHarmonizationArtifacts,
    report_path: Path,
    chain_summary_path: Path,
    site_summary_path: Path,
) -> None:
    write_csv_rows(chain_summary_path, artifacts.chain_rows)
    write_csv_rows(site_summary_path, artifacts.site_rows)
    atomic_write_text(report_path, artifacts.report_markdown)
