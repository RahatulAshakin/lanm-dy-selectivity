"""Phase 3B deterministic ProteinMPNN smoke-run orchestration and summarization."""

from __future__ import annotations

import ast
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
from lanm.models import ProteinMPNNSequenceCatalogRow, ProteinMPNNSmokeSummaryRow
from lanm.paths import REPO_ROOT

LOGGER = logging.getLogger(__name__)

SMOKE_NUM_SEQ_PER_TARGET = 20
SMOKE_SAMPLING_TEMP = "0.1 0.15"
SMOKE_SEED = 37
SMOKE_BATCH_SIZE = 1

_NATIVE_HEADER_RE = re.compile(
    r"^>"
    r"(?P<name>[^,]+), "
    r"score=(?P<score>[^,]+), "
    r"global_score=(?P<global_score>[^,]+), "
    r"fixed_chains=(?P<fixed_chains>\[[^\]]*\]), "
    r"designed_chains=(?P<designed_chains>\[[^\]]*\]), "
    r"(?P<model_field>[^=]+)=(?P<model_name>[^,]+), "
    r"git_hash=(?P<git_hash>[^,]+), "
    r"seed=(?P<seed>\d+)$"
)
_GENERATED_HEADER_RE = re.compile(
    r"^>T=(?P<temperature>[^,]+), "
    r"sample=(?P<sample_number>\d+), "
    r"score=(?P<score>[^,]+), "
    r"global_score=(?P<global_score>[^,]+), "
    r"seq_recovery=(?P<seq_recovery>[^,]+)$"
)


@dataclass(frozen=True, slots=True)
class ProteinMPNNSmokeCampaign:
    campaign_id: str
    backbone_id: str
    backbone_structure_id: str
    design_set_name: str
    designed_chains: tuple[str, ...]
    fixed_context_chains: tuple[str, ...]
    pdb_path: Path
    chain_assignment_path: Path
    fixed_positions_path: Path
    output_dir: Path

    @property
    def designed_chains_text(self) -> str:
        return ",".join(self.designed_chains)

    @property
    def fixed_context_chains_text(self) -> str:
        return ",".join(self.fixed_context_chains)


@dataclass(frozen=True, slots=True)
class ProteinMPNNNativeRecord:
    name: str
    score: float
    global_score: float
    fixed_chains: tuple[str, ...]
    designed_chains: tuple[str, ...]
    model_name: str
    git_hash: str
    seed: int
    sequence: str


@dataclass(frozen=True, slots=True)
class ProteinMPNNGeneratedRecord:
    temperature: float
    sample_number: int
    score: float
    global_score: float
    seq_recovery: float
    sequence: str


@dataclass(frozen=True, slots=True)
class ProteinMPNNParsedOutput:
    native_record: ProteinMPNNNativeRecord
    generated_records: tuple[ProteinMPNNGeneratedRecord, ...]


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


def _split_chain_list(chain_text: str) -> tuple[str, ...]:
    if not chain_text.strip():
        return ()
    return tuple(part.strip() for part in chain_text.split(",") if part.strip())


def _parse_chain_literal(payload: str) -> tuple[str, ...]:
    parsed = ast.literal_eval(payload)
    if not isinstance(parsed, list):
        raise ValueError(f"Unexpected chain list payload: {payload}")
    return tuple(str(item).strip() for item in parsed if str(item).strip())


def _load_csv_rows(path: Path) -> list[dict[str, str]]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_proteinmpnn_root(proteinmpnn_root: Path) -> Path:
    """Validate the ProteinMPNN checkout path and required entrypoints."""
    root = proteinmpnn_root.expanduser()
    if not root.is_absolute():
        raise ValueError("--proteinmpnn-root must be an absolute path")
    require_path(root)
    require_path(root / "helper_scripts" / "parse_multiple_chains.py")
    require_path(root / "protein_mpnn_run.py")
    return root


def discover_proteinmpnn_smoke_campaigns(
    campaign_manifest_path: Path,
    backbone_manifest_path: Path,
    campaigns_dir: Path,
    output_root: Path,
) -> tuple[ProteinMPNNSmokeCampaign, ...]:
    """Discover Phase 3A campaign exports for deterministic Phase 3B smoke runs."""
    manifest_rows = _load_csv_rows(campaign_manifest_path)
    backbone_rows = _load_csv_rows(backbone_manifest_path)
    require_path(campaigns_dir)

    manifest_by_id = {
        str(row["campaign_id"]).strip(): row
        for row in manifest_rows
    }
    backbone_by_id = {
        str(row["backbone_id"]).strip(): row
        for row in backbone_rows
    }

    campaigns: list[ProteinMPNNSmokeCampaign] = []
    for directory in sorted(campaigns_dir.iterdir(), key=lambda path: path.name):
        if not directory.is_dir():
            continue
        campaign_id = directory.name
        row = manifest_by_id.get(campaign_id)
        if row is None:
            raise ValueError(f"No design_campaign_manifest.csv entry for campaign {campaign_id}")
        backbone_id = str(row["backbone_id"]).strip()
        backbone_row = backbone_by_id.get(backbone_id)
        if backbone_row is None:
            raise ValueError(f"No design_backbone_manifest.csv entry for backbone {backbone_id}")

        pdb_path = _repo_path(str(row["exported_pdb_path"]).strip())
        chain_assignment_path = _repo_path(str(row["chain_assignment_path"]).strip())
        fixed_positions_path = _repo_path(str(row["fixed_positions_path"]).strip())
        require_path(pdb_path)
        require_path(chain_assignment_path)
        require_path(fixed_positions_path)

        campaigns.append(
            ProteinMPNNSmokeCampaign(
                campaign_id=campaign_id,
                backbone_id=backbone_id,
                backbone_structure_id=str(backbone_row["structure_id"]).strip(),
                design_set_name=str(row["design_set_name"]).strip(),
                designed_chains=_split_chain_list(str(row["designed_chains"]).strip()),
                fixed_context_chains=_split_chain_list(str(row["fixed_context_chains"]).strip()),
                pdb_path=pdb_path,
                chain_assignment_path=chain_assignment_path,
                fixed_positions_path=fixed_positions_path,
                output_dir=output_root / campaign_id,
            )
        )
    return tuple(campaigns)


def build_parse_multiple_chains_command(
    proteinmpnn_root: Path,
    input_path: Path,
    output_path: Path,
    python_executable: str | None = None,
) -> tuple[str, ...]:
    """Build the deterministic parse_multiple_chains.py command."""
    return (
        python_executable or sys.executable,
        str(proteinmpnn_root / "helper_scripts" / "parse_multiple_chains.py"),
        "--input_path",
        str(input_path),
        "--output_path",
        str(output_path),
    )


def build_proteinmpnn_run_command(
    proteinmpnn_root: Path,
    parsed_jsonl_path: Path,
    chain_assignment_path: Path,
    fixed_positions_path: Path,
    output_dir: Path,
    python_executable: str | None = None,
) -> tuple[str, ...]:
    """Build the deterministic protein_mpnn_run.py smoke command."""
    return (
        python_executable or sys.executable,
        str(proteinmpnn_root / "protein_mpnn_run.py"),
        "--jsonl_path",
        str(parsed_jsonl_path),
        "--chain_id_jsonl",
        str(chain_assignment_path),
        "--fixed_positions_jsonl",
        str(fixed_positions_path),
        "--out_folder",
        str(output_dir),
        "--num_seq_per_target",
        str(SMOKE_NUM_SEQ_PER_TARGET),
        "--sampling_temp",
        SMOKE_SAMPLING_TEMP,
        "--seed",
        str(SMOKE_SEED),
        "--batch_size",
        str(SMOKE_BATCH_SIZE),
    )


def parse_proteinmpnn_fasta(fasta_path: Path) -> ProteinMPNNParsedOutput:
    """Parse a ProteinMPNN FASTA file into native and generated sequence records."""
    require_path(fasta_path)
    lines = [line.strip() for line in fasta_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 2 or len(lines) % 2 != 0:
        raise ValueError(f"Unexpected ProteinMPNN FASTA layout in {fasta_path}")

    native_record: ProteinMPNNNativeRecord | None = None
    generated_records: list[ProteinMPNNGeneratedRecord] = []
    for index in range(0, len(lines), 2):
        header = lines[index]
        sequence = lines[index + 1]
        if native_record is None:
            native_record = _parse_native_record(header=header, sequence=sequence)
            continue
        generated_records.append(_parse_generated_record(header=header, sequence=sequence))
    if native_record is None:
        raise ValueError(f"No native record found in {fasta_path}")
    return ProteinMPNNParsedOutput(
        native_record=native_record,
        generated_records=tuple(generated_records),
    )


def summarize_proteinmpnn_output(
    campaign: ProteinMPNNSmokeCampaign,
    parsed_output: ProteinMPNNParsedOutput,
    fasta_path: Path,
) -> tuple[ProteinMPNNSmokeSummaryRow, tuple[ProteinMPNNSequenceCatalogRow, ...]]:
    """Summarize ProteinMPNN generated sequences for one campaign."""
    native_record = parsed_output.native_record
    if native_record.name != campaign.campaign_id:
        raise ValueError(
            f"ProteinMPNN FASTA name {native_record.name} did not match campaign {campaign.campaign_id}"
        )
    if native_record.designed_chains != campaign.designed_chains:
        raise ValueError(
            f"Designed chains {native_record.designed_chains} did not match manifest {campaign.designed_chains}"
        )

    stripped_sequences = [record.sequence.replace("/", "") for record in parsed_output.generated_records]
    unique_sequence_count = len({record.sequence for record in parsed_output.generated_records})
    generated_sequence_count = len(parsed_output.generated_records)
    sequence_length = len(stripped_sequences[0]) if stripped_sequences else len(native_record.sequence.replace("/", ""))
    mean_identity, min_identity, max_identity = _pairwise_identity_summary(stripped_sequences)
    output_fasta_path = _display_path(fasta_path)

    catalog_rows = tuple(
        ProteinMPNNSequenceCatalogRow(
            campaign_id=campaign.campaign_id,
            backbone_id=campaign.backbone_id,
            design_set_name=campaign.design_set_name,
            designed_chains=campaign.designed_chains_text,
            temperature=record.temperature,
            sample_number=record.sample_number,
            score=record.score,
            global_score=record.global_score,
            seq_recovery=record.seq_recovery,
            sequence_id=f"{campaign.campaign_id}_T{record.temperature:g}_sample_{record.sample_number:02d}",
            designed_sequence=record.sequence,
            sequence_length=len(record.sequence.replace("/", "")),
            output_fasta_path=output_fasta_path,
        )
        for record in parsed_output.generated_records
    )

    summary_row = ProteinMPNNSmokeSummaryRow(
        campaign_id=campaign.campaign_id,
        backbone_id=campaign.backbone_id,
        backbone_structure_id=campaign.backbone_structure_id,
        design_set_name=campaign.design_set_name,
        designed_chains=campaign.designed_chains_text,
        fixed_context_chains=campaign.fixed_context_chains_text,
        generated_sequence_count=generated_sequence_count,
        unique_sequence_count=unique_sequence_count,
        unique_sequence_fraction=(unique_sequence_count / generated_sequence_count) if generated_sequence_count else 0.0,
        sequence_length=sequence_length,
        native_score=native_record.score,
        native_global_score=native_record.global_score,
        mean_pairwise_identity=mean_identity,
        min_pairwise_identity=min_identity,
        max_pairwise_identity=max_identity,
        model_name=native_record.model_name,
        git_hash=native_record.git_hash,
        seed=native_record.seed,
        output_fasta_path=output_fasta_path,
    )
    return summary_row, catalog_rows


def render_proteinmpnn_smoke_markdown(
    proteinmpnn_root: Path,
    summary_rows: Sequence[ProteinMPNNSmokeSummaryRow],
) -> str:
    """Render the Phase 3B ProteinMPNN smoke-run markdown report."""
    lines = [
        "# ProteinMPNN Smoke Run",
        "",
        "Deterministic Phase 3B ProteinMPNN smoke run across the exported Phase 3A campaigns.",
        "",
        f"- ProteinMPNN root: `{proteinmpnn_root}`",
        f"- num_seq_per_target: `{SMOKE_NUM_SEQ_PER_TARGET}`",
        f"- sampling_temp: `{SMOKE_SAMPLING_TEMP}`",
        f"- seed: `{SMOKE_SEED}`",
        f"- batch_size: `{SMOKE_BATCH_SIZE}`",
        "",
        "| campaign_id | sequences_generated | designed_chains | backbone_used | diversity_summary |",
        "| --- | ---: | --- | --- | --- |",
    ]
    for row in summary_rows:
        lines.append(
            "| "
            f"{row.campaign_id} | "
            f"{row.generated_sequence_count} | "
            f"{row.designed_chains or '-'} | "
            f"{row.backbone_id} ({row.backbone_structure_id}) | "
            f"{_format_diversity_summary(row)} |"
        )
    lines.append("")
    return "\n".join(lines)


def run_proteinmpnn_smoke(
    proteinmpnn_root: Path,
    campaign_manifest_path: Path,
    backbone_manifest_path: Path,
    campaigns_dir: Path,
    output_root: Path,
    summary_path: Path,
    sequence_catalog_path: Path,
    report_path: Path,
) -> tuple[tuple[ProteinMPNNSmokeSummaryRow, ...], tuple[ProteinMPNNSequenceCatalogRow, ...]]:
    """Run the deterministic ProteinMPNN smoke workflow and write repo artifacts."""
    root = validate_proteinmpnn_root(proteinmpnn_root)
    campaigns = discover_proteinmpnn_smoke_campaigns(
        campaign_manifest_path=campaign_manifest_path,
        backbone_manifest_path=backbone_manifest_path,
        campaigns_dir=campaigns_dir,
        output_root=output_root,
    )
    if not campaigns:
        raise ValueError(f"No campaigns discovered under {campaigns_dir}")

    summary_rows: list[ProteinMPNNSmokeSummaryRow] = []
    catalog_rows: list[ProteinMPNNSequenceCatalogRow] = []
    for campaign in campaigns:
        LOGGER.info("Running ProteinMPNN smoke workflow for %s", campaign.campaign_id)
        _reset_output_directory(campaign.output_dir)
        parsed_jsonl_path = campaign.output_dir / "parsed_pdbs.jsonl"

        parse_command = build_parse_multiple_chains_command(
            proteinmpnn_root=root,
            input_path=campaign.pdb_path.parent,
            output_path=parsed_jsonl_path,
        )
        atomic_write_text(campaign.output_dir / "parse_command.txt", shlex.join(parse_command) + "\n")
        _run_logged_command(
            command=parse_command,
            cwd=root,
            stdout_path=campaign.output_dir / "parse_stdout.txt",
            stderr_path=campaign.output_dir / "parse_stderr.txt",
        )
        require_path(parsed_jsonl_path)

        run_command = build_proteinmpnn_run_command(
            proteinmpnn_root=root,
            parsed_jsonl_path=parsed_jsonl_path,
            chain_assignment_path=campaign.chain_assignment_path,
            fixed_positions_path=campaign.fixed_positions_path,
            output_dir=campaign.output_dir,
        )
        atomic_write_text(campaign.output_dir / "run_command.txt", shlex.join(run_command) + "\n")
        _run_logged_command(
            command=run_command,
            cwd=root,
            stdout_path=campaign.output_dir / "run_stdout.txt",
            stderr_path=campaign.output_dir / "run_stderr.txt",
        )

        fasta_path = campaign.output_dir / "seqs" / f"{campaign.campaign_id}.fa"
        parsed_output = parse_proteinmpnn_fasta(fasta_path)
        summary_row, catalog = summarize_proteinmpnn_output(
            campaign=campaign,
            parsed_output=parsed_output,
            fasta_path=fasta_path,
        )
        summary_rows.append(summary_row)
        catalog_rows.extend(catalog)

    write_csv_rows(summary_path, summary_rows)
    write_csv_rows(sequence_catalog_path, catalog_rows)
    atomic_write_text(
        report_path,
        render_proteinmpnn_smoke_markdown(
            proteinmpnn_root=root,
            summary_rows=summary_rows,
        ),
    )
    return tuple(summary_rows), tuple(catalog_rows)


def _parse_native_record(header: str, sequence: str) -> ProteinMPNNNativeRecord:
    match = _NATIVE_HEADER_RE.match(header)
    if match is None:
        raise ValueError(f"Unexpected ProteinMPNN native header: {header}")
    return ProteinMPNNNativeRecord(
        name=match.group("name"),
        score=float(match.group("score")),
        global_score=float(match.group("global_score")),
        fixed_chains=_parse_chain_literal(match.group("fixed_chains")),
        designed_chains=_parse_chain_literal(match.group("designed_chains")),
        model_name=match.group("model_name"),
        git_hash=match.group("git_hash"),
        seed=int(match.group("seed")),
        sequence=sequence,
    )


def _parse_generated_record(header: str, sequence: str) -> ProteinMPNNGeneratedRecord:
    match = _GENERATED_HEADER_RE.match(header)
    if match is None:
        raise ValueError(f"Unexpected ProteinMPNN generated header: {header}")
    return ProteinMPNNGeneratedRecord(
        temperature=float(match.group("temperature")),
        sample_number=int(match.group("sample_number")),
        score=float(match.group("score")),
        global_score=float(match.group("global_score")),
        seq_recovery=float(match.group("seq_recovery")),
        sequence=sequence,
    )


def _pairwise_identity_summary(sequences: Sequence[str]) -> tuple[float | None, float | None, float | None]:
    if len(sequences) < 2:
        return None, None, None
    lengths = {len(sequence) for sequence in sequences}
    if lengths != {next(iter(lengths))}:
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


def _format_diversity_summary(row: ProteinMPNNSmokeSummaryRow) -> str:
    unique_summary = f"{row.unique_sequence_count}/{row.generated_sequence_count} unique"
    if row.mean_pairwise_identity is None:
        return f"{unique_summary}; pairwise identity n/a"
    return (
        f"{unique_summary}; pairwise identity mean {row.mean_pairwise_identity:.3f} "
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
