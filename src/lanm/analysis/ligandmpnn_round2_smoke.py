"""Phase 7D deterministic LigandMPNN smoke generation for the round-2 shortlist."""

from __future__ import annotations

import csv
import itertools
import logging
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from lanm.analysis.ligandmpnn_smoke import (
    SMOKE_BATCH_SIZE,
    SMOKE_LIGAND_ATOM_CONTEXT,
    SMOKE_MODEL_TYPE,
    SMOKE_NUMBER_OF_BATCHES,
    SMOKE_NUMBER_OF_PACKS_PER_DESIGN,
    SMOKE_PACK_SIDE_CHAINS,
    SMOKE_PACK_WITH_LIGAND_CONTEXT,
    SMOKE_PACKED_SUFFIX,
    SMOKE_SEED,
    SMOKE_SIDE_CHAIN_CONTEXT,
    SMOKE_TEMPERATURE,
    SMOKE_VERBOSE,
    LigandMPNNParsedOutput,
    parse_ligandmpnn_fasta,
    validate_ligandmpnn_root,
)
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.models import (
    LigandMPNNRound2InputManifestRow,
    LigandMPNNRound2SequenceCatalogRow,
    LigandMPNNRound2SmokeSummaryRow,
)
from lanm.paths import REPO_ROOT
from lanm.structure.geometry import WATER_RESIDUES
from lanm.structure.pdb import read_pdb_atom_records

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class LigandMPNNRound2SmokeCandidate:
    shortlist_rank: int
    candidate_id: str
    seed_rank: int
    seed_candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    source_pdb_path: Path
    designed_chains: tuple[str, ...]
    fixed_context_chains: tuple[str, ...]
    preserved_metal_identity: str
    preserved_solvent_residue_count: int
    redesigned_residue_ids: tuple[str, ...]
    pdb_path: Path
    redesigned_residues_path: Path
    output_dir: Path

    @property
    def designed_chains_text(self) -> str:
        return ",".join(self.designed_chains)

    @property
    def fixed_context_chains_text(self) -> str:
        return ",".join(self.fixed_context_chains)

    @property
    def redesigned_residues_text(self) -> str:
        return " ".join(self.redesigned_residue_ids)

    @property
    def redesigned_residues_csv(self) -> str:
        return ",".join(self.redesigned_residue_ids)


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _repo_path(path_text: str) -> Path:
    candidate = Path(path_text).expanduser()
    if candidate.is_absolute():
        return candidate
    return REPO_ROOT / candidate


def _split_csv_list(value: str) -> tuple[str, ...]:
    if not value.strip():
        return ()
    return tuple(item.strip() for item in value.split(",") if item.strip())


def _read_redesigned_residue_ids(path: Path) -> tuple[str, ...]:
    require_path(path)
    return tuple(item.strip() for item in path.read_text(encoding="utf-8").split() if item.strip())


def _load_ligandmpnn_round2_input_manifest_rows(path: Path) -> tuple[LigandMPNNRound2InputManifestRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            LigandMPNNRound2InputManifestRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                seed_rank=int(str(row["seed_rank"]).strip()),
                seed_candidate_id=str(row["seed_candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                source_pdb_path=str(row["source_pdb_path"]).strip(),
                designed_chains=str(row["designed_chains"]).strip(),
                fixed_context_chains=str(row["fixed_context_chains"]).strip(),
                redesigned_residue_count=int(str(row["redesigned_residue_count"]).strip()),
                redesigned_residue_ids=str(row["redesigned_residue_ids"]).strip(),
                pdb_path=str(row["pdb_path"]).strip(),
                redesigned_residues_path=str(row["redesigned_residues_path"]).strip(),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.shortlist_rank))


def _derive_preserved_context_from_pdb(pdb_path: Path) -> tuple[str, int]:
    atoms = read_pdb_atom_records(pdb_path)
    metal_elements = {
        atom.element_upper
        for atom in atoms
        if atom.record_type == "HETATM" and atom.is_metal
    }
    if not metal_elements:
        raise ValueError(f"No preserved metal identity could be derived from {pdb_path}")
    if len(metal_elements) != 1:
        raise ValueError(
            f"Expected exactly one preserved metal identity in {pdb_path}, found {sorted(metal_elements)}"
        )
    preserved_metal_identity = next(iter(metal_elements))
    preserved_solvent_residue_count = len(
        {
            atom.residue_key
            for atom in atoms
            if atom.record_type == "HETATM" and atom.residue_name in WATER_RESIDUES
        }
    )
    return preserved_metal_identity, preserved_solvent_residue_count


def discover_round2_ligandmpnn_smoke_candidates(
    manifest_path: Path,
    output_root: Path,
) -> tuple[LigandMPNNRound2SmokeCandidate, ...]:
    """Discover the Phase 7D LigandMPNN round-2 smoke candidates from the Phase 7C manifest."""
    manifest_rows = _load_ligandmpnn_round2_input_manifest_rows(manifest_path)
    seen_candidate_ids: set[str] = set()
    candidates: list[LigandMPNNRound2SmokeCandidate] = []
    for row in manifest_rows:
        if row.candidate_id in seen_candidate_ids:
            raise ValueError(f"Duplicate LigandMPNN round-2 candidate_id in manifest: {row.candidate_id}")
        seen_candidate_ids.add(row.candidate_id)

        source_pdb_path = _repo_path(row.source_pdb_path)
        pdb_path = _repo_path(row.pdb_path)
        redesigned_residues_path = _repo_path(row.redesigned_residues_path)
        require_path(source_pdb_path)
        require_path(pdb_path)
        require_path(redesigned_residues_path)

        redesigned_residue_ids = _read_redesigned_residue_ids(redesigned_residues_path)
        expected_redesigned_residue_ids = _split_csv_list(row.redesigned_residue_ids)
        if len(redesigned_residue_ids) != row.redesigned_residue_count:
            raise ValueError(
                f"Redesigned residue count mismatch for {row.candidate_id}: "
                f"{len(redesigned_residue_ids)} vs manifest {row.redesigned_residue_count}"
            )
        if redesigned_residue_ids != expected_redesigned_residue_ids:
            raise ValueError(
                f"Redesigned residue ids mismatch for {row.candidate_id}: "
                f"{redesigned_residue_ids} vs manifest {expected_redesigned_residue_ids}"
            )

        preserved_metal_identity, preserved_solvent_residue_count = _derive_preserved_context_from_pdb(pdb_path)
        candidates.append(
            LigandMPNNRound2SmokeCandidate(
                shortlist_rank=row.shortlist_rank,
                candidate_id=row.candidate_id,
                seed_rank=row.seed_rank,
                seed_candidate_id=row.seed_candidate_id,
                campaign_id=row.campaign_id,
                backbone_id=row.backbone_id,
                topology_class=row.topology_class,
                source_pdb_path=source_pdb_path,
                designed_chains=_split_csv_list(row.designed_chains),
                fixed_context_chains=_split_csv_list(row.fixed_context_chains),
                preserved_metal_identity=preserved_metal_identity,
                preserved_solvent_residue_count=preserved_solvent_residue_count,
                redesigned_residue_ids=redesigned_residue_ids,
                pdb_path=pdb_path,
                redesigned_residues_path=redesigned_residues_path,
                output_dir=output_root / row.candidate_id,
            )
        )
    return tuple(candidates)


def build_round2_ligandmpnn_run_command(
    ligandmpnn_root: Path,
    candidate: LigandMPNNRound2SmokeCandidate,
    python_executable: str | None = None,
) -> tuple[str, ...]:
    """Build the deterministic Phase 7D LigandMPNN smoke command for one shortlisted candidate."""
    return (
        python_executable or sys.executable,
        str(ligandmpnn_root / "run.py"),
        "--model_type",
        SMOKE_MODEL_TYPE,
        "--seed",
        str(SMOKE_SEED),
        "--pdb_path",
        str(candidate.pdb_path),
        "--out_folder",
        str(candidate.output_dir),
        "--chains_to_design",
        candidate.designed_chains_text,
        "--redesigned_residues",
        candidate.redesigned_residues_text,
        "--ligand_mpnn_use_atom_context",
        str(SMOKE_LIGAND_ATOM_CONTEXT),
        "--ligand_mpnn_use_side_chain_context",
        str(SMOKE_SIDE_CHAIN_CONTEXT),
        "--pack_side_chains",
        str(SMOKE_PACK_SIDE_CHAINS),
        "--number_of_packs_per_design",
        str(SMOKE_NUMBER_OF_PACKS_PER_DESIGN),
        "--pack_with_ligand_context",
        str(SMOKE_PACK_WITH_LIGAND_CONTEXT),
        "--batch_size",
        str(SMOKE_BATCH_SIZE),
        "--number_of_batches",
        str(SMOKE_NUMBER_OF_BATCHES),
        "--temperature",
        f"{SMOKE_TEMPERATURE:g}",
        "--verbose",
        str(SMOKE_VERBOSE),
    )


def summarize_round2_ligandmpnn_output(
    candidate: LigandMPNNRound2SmokeCandidate,
    parsed_output: LigandMPNNParsedOutput,
    fasta_path: Path,
) -> tuple[LigandMPNNRound2SmokeSummaryRow, tuple[LigandMPNNRound2SequenceCatalogRow, ...]]:
    """Summarize LigandMPNN generated sequences for one round-2 shortlisted candidate."""
    native_record = parsed_output.native_record
    if native_record.name != candidate.candidate_id:
        raise ValueError(
            f"LigandMPNN FASTA name {native_record.name} did not match candidate {candidate.candidate_id}"
        )

    for record in parsed_output.generated_records:
        if record.name != candidate.candidate_id:
            raise ValueError(
                f"LigandMPNN generated sequence name {record.name} did not match candidate {candidate.candidate_id}"
            )
        if record.seed != native_record.seed:
            raise ValueError(
                f"LigandMPNN seed mismatch for {candidate.candidate_id}: "
                f"{record.seed} vs native {native_record.seed}"
            )

    stripped_sequences = [_strip_chain_separators(record.sequence) for record in parsed_output.generated_records]
    overall_confidences = [record.overall_confidence for record in parsed_output.generated_records]
    ligand_confidences = [record.ligand_confidence for record in parsed_output.generated_records]
    generated_sequence_count = len(parsed_output.generated_records)
    unique_sequence_count = len({record.sequence for record in parsed_output.generated_records})
    sequence_length = len(stripped_sequences[0]) if stripped_sequences else len(_strip_chain_separators(native_record.sequence))
    mean_identity, min_identity, max_identity = _pairwise_identity_summary(stripped_sequences)
    output_fasta_path = _display_path(fasta_path)

    catalog_rows = tuple(
        LigandMPNNRound2SequenceCatalogRow(
            shortlist_rank=candidate.shortlist_rank,
            candidate_id=candidate.candidate_id,
            seed_rank=candidate.seed_rank,
            seed_candidate_id=candidate.seed_candidate_id,
            campaign_id=candidate.campaign_id,
            backbone_id=candidate.backbone_id,
            topology_class=candidate.topology_class,
            designed_chains=candidate.designed_chains_text,
            fixed_context_chains=candidate.fixed_context_chains_text,
            preserved_metal_identity=candidate.preserved_metal_identity,
            redesigned_residues=candidate.redesigned_residues_csv,
            design_id=record.design_id,
            temperature=record.temperature,
            seed=record.seed,
            overall_confidence=record.overall_confidence,
            ligand_confidence=record.ligand_confidence,
            seq_recovery=record.seq_recovery,
            sequence_id=f"{candidate.candidate_id}_design_{record.design_id:02d}",
            designed_sequence=record.sequence,
            sequence_length=len(_strip_chain_separators(record.sequence)),
            output_fasta_path=output_fasta_path,
            backbone_pdb_path=_display_path(candidate.output_dir / "backbones" / f"{candidate.candidate_id}_{record.design_id}.pdb"),
            packed_pdb_path=_display_path(
                candidate.output_dir
                / "packed"
                / f"{candidate.candidate_id}{SMOKE_PACKED_SUFFIX}_{record.design_id}_1.pdb"
            ),
        )
        for record in parsed_output.generated_records
    )

    summary_row = LigandMPNNRound2SmokeSummaryRow(
        shortlist_rank=candidate.shortlist_rank,
        candidate_id=candidate.candidate_id,
        seed_rank=candidate.seed_rank,
        seed_candidate_id=candidate.seed_candidate_id,
        campaign_id=candidate.campaign_id,
        backbone_id=candidate.backbone_id,
        topology_class=candidate.topology_class,
        designed_chains=candidate.designed_chains_text,
        fixed_context_chains=candidate.fixed_context_chains_text,
        preserved_metal_identity=candidate.preserved_metal_identity,
        preserved_solvent_residue_count=candidate.preserved_solvent_residue_count,
        redesigned_residues=candidate.redesigned_residues_csv,
        generated_sequence_count=generated_sequence_count,
        unique_sequence_count=unique_sequence_count,
        unique_sequence_fraction=(unique_sequence_count / generated_sequence_count) if generated_sequence_count else 0.0,
        sequence_length=sequence_length,
        mean_overall_confidence=_mean_or_none(overall_confidences),
        min_overall_confidence=min(overall_confidences) if overall_confidences else None,
        max_overall_confidence=max(overall_confidences) if overall_confidences else None,
        mean_ligand_confidence=_mean_or_none(ligand_confidences),
        min_ligand_confidence=min(ligand_confidences) if ligand_confidences else None,
        max_ligand_confidence=max(ligand_confidences) if ligand_confidences else None,
        mean_pairwise_identity=mean_identity,
        min_pairwise_identity=min_identity,
        max_pairwise_identity=max_identity,
        temperature=native_record.temperature,
        seed=native_record.seed,
        output_fasta_path=output_fasta_path,
    )
    return summary_row, catalog_rows


def render_ligandmpnn_round2_smoke_markdown(
    ligandmpnn_root: Path,
    summary_rows: Sequence[LigandMPNNRound2SmokeSummaryRow],
) -> str:
    """Render the Phase 7D LigandMPNN round-2 smoke markdown report."""
    lines = [
        "# LigandMPNN Round-2 Smoke",
        "",
        "Deterministic Phase 7D LigandMPNN smoke pass across the 6 shortlisted round-2 candidates from Phase 7C.",
        "",
        f"- LigandMPNN root: `{ligandmpnn_root}`",
        f"- model_type: `{SMOKE_MODEL_TYPE}`",
        f"- seed: `{SMOKE_SEED}`",
        f"- batch_size: `{SMOKE_BATCH_SIZE}`",
        f"- number_of_batches: `{SMOKE_NUMBER_OF_BATCHES}`",
        f"- temperature: `{SMOKE_TEMPERATURE:g}`",
        f"- ligand atom context: `{bool(SMOKE_LIGAND_ATOM_CONTEXT)}`",
        f"- side-chain context: `{bool(SMOKE_SIDE_CHAIN_CONTEXT)}`",
        f"- pack_side_chains: `{bool(SMOKE_PACK_SIDE_CHAINS)}`",
        f"- pack_with_ligand_context: `{bool(SMOKE_PACK_WITH_LIGAND_CONTEXT)}`",
        "- No Rosetta, MD, QM, or quantum steps were run in this phase.",
        "",
        "| candidate_id | seed_scaffold | preserved_metal_identity | redesigned_residues | sequence_count | overall_confidence | ligand_confidence | unique_sequence_count |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in summary_rows:
        lines.append(
            "| "
            f"{row.candidate_id} | "
            f"{row.seed_candidate_id} | "
            f"{row.preserved_metal_identity} | "
            f"{row.redesigned_residues or 'native'} | "
            f"{row.generated_sequence_count} | "
            f"{_format_confidence_value(row.mean_overall_confidence)} | "
            f"{_format_confidence_value(row.mean_ligand_confidence)} | "
            f"{row.unique_sequence_count} |"
        )
    lines.append("")
    for row in summary_rows:
        lines.extend(
            [
                f"## {row.candidate_id}",
                "",
                f"- Seed scaffold: `{row.seed_candidate_id}`",
                f"- Backbone/topology: `{row.backbone_id}` / `{row.topology_class}`",
                f"- Preserved metal identity: `{row.preserved_metal_identity}`",
                f"- Preserved solvent residue count: `{row.preserved_solvent_residue_count}`",
                f"- Redesigned residues: `{row.redesigned_residues or '-'}`",
                f"- Sequence count: `{row.generated_sequence_count}` with `{row.unique_sequence_count}` unique sequences",
                f"- Overall confidence: `{_format_confidence_summary(row.mean_overall_confidence, row.min_overall_confidence, row.max_overall_confidence)}`",
                f"- Ligand confidence: `{_format_confidence_summary(row.mean_ligand_confidence, row.min_ligand_confidence, row.max_ligand_confidence)}`",
                "",
            ]
        )
    return "\n".join(lines)


def run_ligandmpnn_round2_smoke(
    ligandmpnn_root: Path,
    manifest_path: Path,
    output_root: Path,
    summary_path: Path,
    sequence_catalog_path: Path,
    report_path: Path,
) -> tuple[tuple[LigandMPNNRound2SmokeSummaryRow, ...], tuple[LigandMPNNRound2SequenceCatalogRow, ...]]:
    """Run the deterministic Phase 7D LigandMPNN smoke workflow and write repo artifacts."""
    root = validate_ligandmpnn_root(ligandmpnn_root)
    candidates = discover_round2_ligandmpnn_smoke_candidates(
        manifest_path=manifest_path,
        output_root=output_root,
    )
    if not candidates:
        raise ValueError(f"No LigandMPNN round-2 smoke candidates discovered from {manifest_path}")

    summary_rows: list[LigandMPNNRound2SmokeSummaryRow] = []
    catalog_rows: list[LigandMPNNRound2SequenceCatalogRow] = []
    for candidate in candidates:
        LOGGER.info("Running LigandMPNN round-2 smoke workflow for %s", candidate.candidate_id)
        _reset_output_directory(candidate.output_dir)

        run_command = build_round2_ligandmpnn_run_command(
            ligandmpnn_root=root,
            candidate=candidate,
        )
        atomic_write_text(candidate.output_dir / "run_command.txt", shlex.join(run_command) + "\n")
        _run_logged_command(
            command=run_command,
            cwd=root,
            stdout_path=candidate.output_dir / "run_stdout.txt",
            stderr_path=candidate.output_dir / "run_stderr.txt",
        )

        fasta_path = candidate.output_dir / "seqs" / f"{candidate.candidate_id}.fa"
        parsed_output = parse_ligandmpnn_fasta(fasta_path)
        summary_row, candidate_catalog_rows = summarize_round2_ligandmpnn_output(
            candidate=candidate,
            parsed_output=parsed_output,
            fasta_path=fasta_path,
        )
        summary_rows.append(summary_row)
        catalog_rows.extend(candidate_catalog_rows)

    write_csv_rows(summary_path, summary_rows)
    write_csv_rows(sequence_catalog_path, catalog_rows)
    atomic_write_text(
        report_path,
        render_ligandmpnn_round2_smoke_markdown(
            ligandmpnn_root=root,
            summary_rows=summary_rows,
        ),
    )
    return tuple(summary_rows), tuple(catalog_rows)


def _strip_chain_separators(sequence: str) -> str:
    return "".join(character for character in sequence if character.isalpha())


def _mean_or_none(values: Sequence[float]) -> float | None:
    if not values:
        return None
    return sum(values) / len(values)


def _pairwise_identity_summary(sequences: Sequence[str]) -> tuple[float | None, float | None, float | None]:
    if len(sequences) < 2:
        return None, None, None
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        return None, None, None
    sequence_length = next(iter(lengths))
    if sequence_length == 0:
        return None, None, None

    identities = [
        sum(left_char == right_char for left_char, right_char in zip(left, right)) / sequence_length
        for left, right in itertools.combinations(sequences, 2)
    ]
    return (
        sum(identities) / len(identities),
        min(identities),
        max(identities),
    )


def _format_confidence_value(value: float | None) -> str:
    if value is None:
        return "n/a"
    return f"{value:.4f}"


def _format_confidence_summary(
    mean_value: float | None,
    min_value: float | None,
    max_value: float | None,
) -> str:
    if mean_value is None or min_value is None or max_value is None:
        return "n/a"
    return f"mean {mean_value:.4f}; min {min_value:.4f}; max {max_value:.4f}"


def _reset_output_directory(output_dir: Path) -> None:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)


def _run_logged_command(
    command: Sequence[str],
    cwd: Path,
    stdout_path: Path,
    stderr_path: Path,
) -> None:
    completed = subprocess.run(
        list(command),
        cwd=cwd,
        check=False,
        capture_output=True,
        text=True,
    )
    atomic_write_text(stdout_path, completed.stdout)
    atomic_write_text(stderr_path, completed.stderr)
    if completed.returncode != 0:
        quoted_command = shlex.join(command)
        raise RuntimeError(f"Command failed with exit code {completed.returncode}: {quoted_command}")
