"""Phase 4B deterministic LigandMPNN smoke-run orchestration and summarization."""

from __future__ import annotations

import csv
import itertools
import logging
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.models import (
    LigandMPNNInputManifestRow,
    LigandMPNNSequenceCatalogRow,
    LigandMPNNSmokeSummaryRow,
)
from lanm.paths import REPO_ROOT

LOGGER = logging.getLogger(__name__)

SMOKE_MODEL_TYPE = "ligand_mpnn"
SMOKE_SEED = 37
SMOKE_BATCH_SIZE = 2
SMOKE_NUMBER_OF_BATCHES = 1
SMOKE_TEMPERATURE = 0.1
SMOKE_LIGAND_ATOM_CONTEXT = 1
SMOKE_SIDE_CHAIN_CONTEXT = 1
SMOKE_PACK_SIDE_CHAINS = 1
SMOKE_NUMBER_OF_PACKS_PER_DESIGN = 1
SMOKE_PACK_WITH_LIGAND_CONTEXT = 1
SMOKE_VERBOSE = 0
SMOKE_PACKED_SUFFIX = "_packed"
SMOKE_REQUIRED_CHECKPOINTS = (
    "model_params/ligandmpnn_v_32_010_25.pt",
    "model_params/ligandmpnn_sc_v_32_002_16.pt",
)

_NATIVE_HEADER_RE = re.compile(
    r"^>"
    r"(?P<name>[^,]+), "
    r"T=(?P<temperature>[^,]+), "
    r"seed=(?P<seed>\d+), "
    r"num_res=(?P<num_res>\d+), "
    r"num_ligand_res=(?P<num_ligand_res>\d+), "
    r"use_ligand_context=(?P<use_ligand_context>True|False), "
    r"ligand_cutoff_distance=(?P<ligand_cutoff_distance>[^,]+), "
    r"batch_size=(?P<batch_size>\d+), "
    r"number_of_batches=(?P<number_of_batches>\d+), "
    r"model_path=(?P<model_path>.+)$"
)
_GENERATED_HEADER_RE = re.compile(
    r"^>"
    r"(?P<name>[^,]+), "
    r"id=(?P<design_id>\d+), "
    r"T=(?P<temperature>[^,]+), "
    r"seed=(?P<seed>\d+), "
    r"overall_confidence=(?P<overall_confidence>[^,]+), "
    r"ligand_confidence=(?P<ligand_confidence>[^,]+), "
    r"seq_rec=(?P<seq_recovery>[^,]+)$"
)


@dataclass(frozen=True, slots=True)
class LigandMPNNSmokeCandidate:
    shortlist_rank: int
    candidate_id: str
    campaign_id: str
    backbone_id: str
    source_structure_id: str
    designed_chains: tuple[str, ...]
    preserved_metal_identity: str
    redesigned_residue_ids: tuple[str, ...]
    pdb_path: Path
    redesigned_residues_path: Path
    output_dir: Path

    @property
    def designed_chains_text(self) -> str:
        return ",".join(self.designed_chains)

    @property
    def redesigned_residues_text(self) -> str:
        return " ".join(self.redesigned_residue_ids)

    @property
    def redesigned_residues_csv(self) -> str:
        return ",".join(self.redesigned_residue_ids)


@dataclass(frozen=True, slots=True)
class LigandMPNNNativeRecord:
    name: str
    temperature: float
    seed: int
    num_res: int
    num_ligand_res: int
    use_ligand_context: bool
    ligand_cutoff_distance: float
    batch_size: int
    number_of_batches: int
    model_path: str
    sequence: str


@dataclass(frozen=True, slots=True)
class LigandMPNNGeneratedRecord:
    name: str
    design_id: int
    temperature: float
    seed: int
    overall_confidence: float
    ligand_confidence: float
    seq_recovery: float
    sequence: str


@dataclass(frozen=True, slots=True)
class LigandMPNNParsedOutput:
    native_record: LigandMPNNNativeRecord
    generated_records: tuple[LigandMPNNGeneratedRecord, ...]


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


def _load_ligandmpnn_input_manifest_rows(path: Path) -> tuple[LigandMPNNInputManifestRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            LigandMPNNInputManifestRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                source_structure_id=str(row["source_structure_id"]).strip(),
                source_kind=str(row["source_kind"]).strip(),
                source_path=str(row["source_path"]).strip(),
                preserved_chains=str(row["preserved_chains"]).strip(),
                designed_chains=str(row["designed_chains"]).strip(),
                fixed_context_chains=str(row["fixed_context_chains"]).strip(),
                preserved_metal_identity=str(row["preserved_metal_identity"]).strip(),
                preserved_metal_site_count=int(str(row["preserved_metal_site_count"]).strip()),
                preserved_solvent_residue_count=int(str(row["preserved_solvent_residue_count"]).strip()),
                redesigned_residue_count=int(str(row["redesigned_residue_count"]).strip()),
                pdb_path=str(row["pdb_path"]).strip(),
                redesigned_residues_path=str(row["redesigned_residues_path"]).strip(),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.shortlist_rank))


def validate_ligandmpnn_root(ligandmpnn_root: Path) -> Path:
    """Validate the LigandMPNN checkout path and required runtime files."""
    root = ligandmpnn_root.expanduser()
    if not root.is_absolute():
        raise ValueError("--ligandmpnn-root must be an absolute path")
    require_path(root)
    require_path(root / "run.py")

    missing_checkpoints = [
        root / checkpoint_path
        for checkpoint_path in SMOKE_REQUIRED_CHECKPOINTS
        if not (root / checkpoint_path).exists()
    ]
    if missing_checkpoints:
        checkpoint_list = ", ".join(_display_path(path) for path in missing_checkpoints)
        raise ValueError(
            "LigandMPNN checkout is missing required model checkpoints: "
            f"{checkpoint_list}. Run `bash get_model_params.sh ./model_params` inside the checkout first."
        )
    return root


def discover_ligandmpnn_smoke_candidates(
    manifest_path: Path,
    output_root: Path,
) -> tuple[LigandMPNNSmokeCandidate, ...]:
    """Discover deterministic LigandMPNN smoke candidates from the Phase 4A manifest."""
    manifest_rows = _load_ligandmpnn_input_manifest_rows(manifest_path)
    seen_candidate_ids: set[str] = set()
    candidates: list[LigandMPNNSmokeCandidate] = []
    for row in manifest_rows:
        if row.candidate_id in seen_candidate_ids:
            raise ValueError(f"Duplicate LigandMPNN candidate_id in manifest: {row.candidate_id}")
        seen_candidate_ids.add(row.candidate_id)

        pdb_path = _repo_path(row.pdb_path)
        redesigned_residues_path = _repo_path(row.redesigned_residues_path)
        require_path(pdb_path)
        require_path(redesigned_residues_path)

        redesigned_residue_ids = _read_redesigned_residue_ids(redesigned_residues_path)
        if len(redesigned_residue_ids) != row.redesigned_residue_count:
            raise ValueError(
                f"Redesigned residue count mismatch for {row.candidate_id}: "
                f"{len(redesigned_residue_ids)} vs manifest {row.redesigned_residue_count}"
            )

        candidates.append(
            LigandMPNNSmokeCandidate(
                shortlist_rank=row.shortlist_rank,
                candidate_id=row.candidate_id,
                campaign_id=row.campaign_id,
                backbone_id=row.backbone_id,
                source_structure_id=row.source_structure_id,
                designed_chains=_split_csv_list(row.designed_chains),
                preserved_metal_identity=row.preserved_metal_identity,
                redesigned_residue_ids=redesigned_residue_ids,
                pdb_path=pdb_path,
                redesigned_residues_path=redesigned_residues_path,
                output_dir=output_root / row.candidate_id,
            )
        )
    return tuple(candidates)


def build_ligandmpnn_run_command(
    ligandmpnn_root: Path,
    candidate: LigandMPNNSmokeCandidate,
    python_executable: str | None = None,
) -> tuple[str, ...]:
    """Build the deterministic LigandMPNN smoke command for one shortlisted candidate."""
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


def parse_ligandmpnn_fasta(fasta_path: Path) -> LigandMPNNParsedOutput:
    """Parse a LigandMPNN FASTA file into native and generated sequence records."""
    require_path(fasta_path)
    lines = [line.strip() for line in fasta_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 2 or len(lines) % 2 != 0:
        raise ValueError(f"Unexpected LigandMPNN FASTA layout in {fasta_path}")

    native_record: LigandMPNNNativeRecord | None = None
    generated_records: list[LigandMPNNGeneratedRecord] = []
    for index in range(0, len(lines), 2):
        header = lines[index]
        sequence = lines[index + 1]
        if native_record is None:
            native_record = _parse_native_record(header=header, sequence=sequence)
            continue
        generated_records.append(_parse_generated_record(header=header, sequence=sequence))
    if native_record is None:
        raise ValueError(f"No native record found in {fasta_path}")
    return LigandMPNNParsedOutput(
        native_record=native_record,
        generated_records=tuple(generated_records),
    )


def summarize_ligandmpnn_output(
    candidate: LigandMPNNSmokeCandidate,
    parsed_output: LigandMPNNParsedOutput,
    fasta_path: Path,
) -> tuple[LigandMPNNSmokeSummaryRow, tuple[LigandMPNNSequenceCatalogRow, ...]]:
    """Summarize LigandMPNN generated sequences for one shortlisted candidate."""
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
    unique_sequence_count = len({record.sequence for record in parsed_output.generated_records})
    generated_sequence_count = len(parsed_output.generated_records)
    sequence_length = len(stripped_sequences[0]) if stripped_sequences else len(_strip_chain_separators(native_record.sequence))
    mean_identity, min_identity, max_identity = _pairwise_identity_summary(stripped_sequences)
    output_fasta_path = _display_path(fasta_path)

    catalog_rows = tuple(
        LigandMPNNSequenceCatalogRow(
            shortlist_rank=candidate.shortlist_rank,
            candidate_id=candidate.candidate_id,
            campaign_id=candidate.campaign_id,
            backbone_id=candidate.backbone_id,
            source_structure_id=candidate.source_structure_id,
            designed_chains=candidate.designed_chains_text,
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

    summary_row = LigandMPNNSmokeSummaryRow(
        shortlist_rank=candidate.shortlist_rank,
        candidate_id=candidate.candidate_id,
        campaign_id=candidate.campaign_id,
        backbone_id=candidate.backbone_id,
        source_structure_id=candidate.source_structure_id,
        designed_chains=candidate.designed_chains_text,
        preserved_metal_identity=candidate.preserved_metal_identity,
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


def render_ligandmpnn_smoke_markdown(
    ligandmpnn_root: Path,
    summary_rows: Sequence[LigandMPNNSmokeSummaryRow],
) -> str:
    """Render the Phase 4B LigandMPNN smoke-run markdown report."""
    lines = [
        "# LigandMPNN Smoke Run",
        "",
        "Deterministic Phase 4B LigandMPNN smoke pass across the 12 shortlisted candidates from Phase 4A.",
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
        "",
        "| candidate_id | source_backbone | preserved_metal_identity | redesigned_residues | sequence_count | confidence_diversity_summary |",
        "| --- | --- | --- | --- | ---: | --- |",
    ]
    for row in summary_rows:
        lines.append(
            "| "
            f"{row.candidate_id} | "
            f"{row.backbone_id} ({row.source_structure_id}) | "
            f"{row.preserved_metal_identity} | "
            f"{row.redesigned_residues or 'native'} | "
            f"{row.generated_sequence_count} | "
            f"{_format_confidence_diversity_summary(row)} |"
        )
    lines.append("")
    return "\n".join(lines)


def run_ligandmpnn_smoke(
    ligandmpnn_root: Path,
    manifest_path: Path,
    output_root: Path,
    summary_path: Path,
    sequence_catalog_path: Path,
    report_path: Path,
) -> tuple[tuple[LigandMPNNSmokeSummaryRow, ...], tuple[LigandMPNNSequenceCatalogRow, ...]]:
    """Run the deterministic LigandMPNN smoke workflow and write repo artifacts."""
    root = validate_ligandmpnn_root(ligandmpnn_root)
    candidates = discover_ligandmpnn_smoke_candidates(
        manifest_path=manifest_path,
        output_root=output_root,
    )
    if not candidates:
        raise ValueError(f"No LigandMPNN smoke candidates discovered from {manifest_path}")

    summary_rows: list[LigandMPNNSmokeSummaryRow] = []
    catalog_rows: list[LigandMPNNSequenceCatalogRow] = []
    for candidate in candidates:
        LOGGER.info("Running LigandMPNN smoke workflow for %s", candidate.candidate_id)
        _reset_output_directory(candidate.output_dir)

        run_command = build_ligandmpnn_run_command(
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
        summary_row, catalog = summarize_ligandmpnn_output(
            candidate=candidate,
            parsed_output=parsed_output,
            fasta_path=fasta_path,
        )
        summary_rows.append(summary_row)
        catalog_rows.extend(catalog)

    write_csv_rows(summary_path, summary_rows)
    write_csv_rows(sequence_catalog_path, catalog_rows)
    atomic_write_text(
        report_path,
        render_ligandmpnn_smoke_markdown(
            ligandmpnn_root=root,
            summary_rows=summary_rows,
        ),
    )
    return tuple(summary_rows), tuple(catalog_rows)


def _parse_native_record(header: str, sequence: str) -> LigandMPNNNativeRecord:
    match = _NATIVE_HEADER_RE.match(header)
    if match is None:
        raise ValueError(f"Unexpected LigandMPNN native header: {header}")
    return LigandMPNNNativeRecord(
        name=match.group("name"),
        temperature=float(match.group("temperature")),
        seed=int(match.group("seed")),
        num_res=int(match.group("num_res")),
        num_ligand_res=int(match.group("num_ligand_res")),
        use_ligand_context=(match.group("use_ligand_context") == "True"),
        ligand_cutoff_distance=float(match.group("ligand_cutoff_distance")),
        batch_size=int(match.group("batch_size")),
        number_of_batches=int(match.group("number_of_batches")),
        model_path=match.group("model_path"),
        sequence=sequence,
    )


def _parse_generated_record(header: str, sequence: str) -> LigandMPNNGeneratedRecord:
    match = _GENERATED_HEADER_RE.match(header)
    if match is None:
        raise ValueError(f"Unexpected LigandMPNN generated header: {header}")
    return LigandMPNNGeneratedRecord(
        name=match.group("name"),
        design_id=int(match.group("design_id")),
        temperature=float(match.group("temperature")),
        seed=int(match.group("seed")),
        overall_confidence=float(match.group("overall_confidence")),
        ligand_confidence=float(match.group("ligand_confidence")),
        seq_recovery=float(match.group("seq_recovery")),
        sequence=sequence,
    )


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


def _format_confidence_diversity_summary(row: LigandMPNNSmokeSummaryRow) -> str:
    confidence_summary = (
        f"overall mean {row.mean_overall_confidence:.3f} "
        f"(min {row.min_overall_confidence:.3f}, max {row.max_overall_confidence:.3f}); "
        f"ligand mean {row.mean_ligand_confidence:.3f} "
        f"(min {row.min_ligand_confidence:.3f}, max {row.max_ligand_confidence:.3f})"
        if row.mean_overall_confidence is not None and row.mean_ligand_confidence is not None
        else "confidence n/a"
    )
    diversity_summary = f"{row.unique_sequence_count}/{row.generated_sequence_count} unique"
    if row.mean_pairwise_identity is None:
        return f"{confidence_summary}; {diversity_summary}; pairwise identity n/a"
    return (
        f"{confidence_summary}; {diversity_summary}; "
        f"pairwise identity mean {row.mean_pairwise_identity:.3f} "
        f"(min {row.min_pairwise_identity:.3f}, max {row.max_pairwise_identity:.3f})"
    )


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
