"""Phase 7B deterministic ProteinMPNN smoke generation for focused round-2 redesign."""

from __future__ import annotations

import csv
import itertools
import json
import logging
import shlex
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

try:
    import yaml  # type: ignore
except ModuleNotFoundError:  # pragma: no cover - optional dependency in local envs
    yaml = None

from lanm.analysis.proteinmpnn_smoke import (
    SMOKE_BATCH_SIZE,
    SMOKE_NUM_SEQ_PER_TARGET,
    SMOKE_SAMPLING_TEMP,
    SMOKE_SEED,
    build_parse_multiple_chains_command,
    build_proteinmpnn_run_command,
    parse_proteinmpnn_fasta,
    validate_proteinmpnn_root,
)
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.models import ProteinMPNNRound2SequenceCatalogRow, ProteinMPNNRound2SummaryRow
from lanm.paths import REPO_ROOT
from lanm.structure.pdb import read_pdb_atom_records
from lanm.structure.residues import collect_polymer_residues

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class Round2DesignedResidue:
    chain_id: str
    sequence_index: int
    residue_seq: int
    insertion_code: str
    residue_name: str
    seed_amino_acid: str

    @property
    def residue_id(self) -> str:
        return f"{self.chain_id}{self.residue_seq}{self.insertion_code}".strip()


@dataclass(frozen=True, slots=True)
class Round2RedesignPosition:
    chain_id: str
    sequence_index: int
    residue_id: str
    seed_amino_acid: str
    canonical_family_position: int
    am1_mature_position: int


@dataclass(frozen=True, slots=True)
class ProteinMPNNRound2Campaign:
    seed_rank: int
    candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    design_set_name: str
    designed_chains: tuple[str, ...]
    fixed_context_chains: tuple[str, ...]
    pdb_path: Path
    chain_assignment_path: Path
    fixed_positions_path: Path
    output_dir: Path
    seed_scaffold_sequence: str
    designed_residues: tuple[Round2DesignedResidue, ...]
    redesignable_positions: tuple[Round2RedesignPosition, ...]

    @property
    def designed_chains_text(self) -> str:
        return ",".join(self.designed_chains)

    @property
    def fixed_context_chains_text(self) -> str:
        return ",".join(self.fixed_context_chains)

    @property
    def redesignable_residue_ids_text(self) -> str:
        return ",".join(position.residue_id for position in self.redesignable_positions)

    @property
    def redesignable_canonical_positions_text(self) -> str:
        return ",".join(str(position.canonical_family_position) for position in self.redesignable_positions)

    @property
    def redesignable_am1_positions_text(self) -> str:
        return ",".join(str(position.am1_mature_position) for position in self.redesignable_positions)


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


def _split_designed_sequence(sequence: str, designed_chains: Sequence[str]) -> tuple[str, ...]:
    if len(designed_chains) == 1:
        return (sequence.replace("/", ""),)
    parts = tuple(part.strip() for part in sequence.split("/"))
    if len(parts) != len(designed_chains):
        raise ValueError(
            f"Expected {len(designed_chains)} designed-chain parts, found {len(parts)} in {sequence!r}"
        )
    return tuple(part.replace("/", "") for part in parts)


def _load_csv_rows(path: Path) -> list[dict[str, str]]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _load_yaml_mapping(path: Path) -> dict[str, Any]:
    require_path(path)
    if yaml is None:
        raise RuntimeError("PyYAML is required to read redesign_round2 configuration inputs")
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Unexpected YAML payload in {path}")
    return payload


def _load_single_jsonl_mapping(path: Path) -> dict[str, Any]:
    require_path(path)
    rows = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one JSONL object in {path}")
    payload = json.loads(rows[0])
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return payload


def _extract_seed_scaffold(
    pdb_path: Path,
    designed_chains: Sequence[str],
) -> tuple[str, tuple[Round2DesignedResidue, ...]]:
    atoms = read_pdb_atom_records(pdb_path)
    residues_by_chain = collect_polymer_residues(atoms)
    designed_residues: list[Round2DesignedResidue] = []
    sequence_parts: list[str] = []
    for chain_id in designed_chains:
        chain_residues = residues_by_chain.get(chain_id)
        if not chain_residues:
            raise ValueError(f"No polymer residues found for designed chain {chain_id} in {pdb_path}")
        sequence_parts.append("".join(residue.one_letter_code for residue in chain_residues))
        designed_residues.extend(
            Round2DesignedResidue(
                chain_id=chain_id,
                sequence_index=index,
                residue_seq=residue.residue_seq,
                insertion_code=residue.insertion_code,
                residue_name=residue.residue_name,
                seed_amino_acid=residue.one_letter_code,
            )
            for index, residue in enumerate(chain_residues, start=1)
        )
    if len(sequence_parts) == 1:
        return sequence_parts[0], tuple(designed_residues)
    return "/".join(sequence_parts), tuple(designed_residues)


def discover_round2_proteinmpnn_campaigns(
    redesign_round2_config_path: Path,
    round2_seed_table_path: Path,
    round2_position_table_path: Path,
    campaigns_dir: Path,
    output_root: Path,
) -> tuple[ProteinMPNNRound2Campaign, ...]:
    """Discover focused round-2 ProteinMPNN campaigns from Phase 7A artifacts."""
    config_payload = _load_yaml_mapping(redesign_round2_config_path)
    configured_seeds = config_payload.get("seeds")
    if not isinstance(configured_seeds, dict):
        raise ValueError("redesign_round2.yaml is missing a seeds mapping")
    configured_seed_ids = tuple(
        str(candidate_id).strip()
        for candidate_id in config_payload.get("selection_policy", {}).get("seed_candidate_ids", [])
        if str(candidate_id).strip()
    )

    seed_rows = sorted(
        _load_csv_rows(round2_seed_table_path),
        key=lambda row: int(str(row["seed_rank"]).strip()),
    )
    position_rows = _load_csv_rows(round2_position_table_path)
    require_path(campaigns_dir)

    discovered_directory_ids = {
        directory.name
        for directory in campaigns_dir.iterdir()
        if directory.is_dir()
    }
    expected_seed_ids = {str(row["candidate_id"]).strip() for row in seed_rows}
    missing_directories = expected_seed_ids - discovered_directory_ids
    extra_directories = discovered_directory_ids - expected_seed_ids
    if missing_directories:
        raise ValueError(
            f"Missing round-2 campaign exports under {campaigns_dir}: {','.join(sorted(missing_directories))}"
        )
    if extra_directories:
        raise ValueError(
            f"Unexpected round-2 campaign exports under {campaigns_dir}: {','.join(sorted(extra_directories))}"
        )

    row_seed_ids = tuple(str(row["candidate_id"]).strip() for row in seed_rows)
    if configured_seed_ids and configured_seed_ids != row_seed_ids:
        raise ValueError(
            "Round-2 seed ordering mismatch between redesign_round2.yaml and redesign_round2_seeds.csv"
        )

    positions_by_candidate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in position_rows:
        positions_by_candidate[str(row["candidate_id"]).strip()].append(row)

    campaigns: list[ProteinMPNNRound2Campaign] = []
    for row in seed_rows:
        candidate_id = str(row["candidate_id"]).strip()
        config_seed = configured_seeds.get(candidate_id)
        if not isinstance(config_seed, dict):
            raise ValueError(f"redesign_round2.yaml is missing seed metadata for {candidate_id}")

        designed_chains = _split_csv_list(str(row["designed_chains"]).strip())
        fixed_context_chains = _split_csv_list(str(row["fixed_context_chains"]).strip())
        config_designed_chains = tuple(str(chain_id).strip() for chain_id in config_seed.get("designed_chains", []))
        config_fixed_context_chains = tuple(
            str(chain_id).strip() for chain_id in config_seed.get("fixed_context_chains", [])
        )
        if designed_chains != config_designed_chains:
            raise ValueError(f"Designed chains mismatch for {candidate_id}: {designed_chains} vs {config_designed_chains}")
        if fixed_context_chains != config_fixed_context_chains:
            raise ValueError(
                f"Fixed-context chains mismatch for {candidate_id}: {fixed_context_chains} vs {config_fixed_context_chains}"
            )

        pdb_path = _repo_path(str(row["exported_pdb_path"]).strip())
        chain_assignment_path = _repo_path(str(row["chain_assignment_path"]).strip())
        fixed_positions_path = _repo_path(str(row["fixed_positions_path"]).strip())
        require_path(pdb_path)
        require_path(chain_assignment_path)
        require_path(fixed_positions_path)
        if _display_path(pdb_path) != str(config_seed.get("exported_pdb_path", "")).strip():
            raise ValueError(f"Exported PDB mismatch for {candidate_id}")
        if _display_path(chain_assignment_path) != str(config_seed.get("chain_assignment_path", "")).strip():
            raise ValueError(f"Chain-assignment path mismatch for {candidate_id}")
        if _display_path(fixed_positions_path) != str(config_seed.get("fixed_positions_path", "")).strip():
            raise ValueError(f"Fixed-positions path mismatch for {candidate_id}")

        seed_scaffold_sequence, designed_residues = _extract_seed_scaffold(
            pdb_path=pdb_path,
            designed_chains=designed_chains,
        )
        designed_residue_map = {
            (residue.chain_id, residue.sequence_index): residue
            for residue in designed_residues
        }

        redesign_rows = sorted(
            positions_by_candidate.get(candidate_id, []),
            key=lambda entry: (str(entry["chain_id"]).strip(), int(str(entry["sequence_index"]).strip())),
        )
        if not redesign_rows:
            raise ValueError(f"No redesign_round2_positions.csv rows found for {candidate_id}")

        redesignable_positions: list[Round2RedesignPosition] = []
        for redesign_row in redesign_rows:
            chain_id = str(redesign_row["chain_id"]).strip()
            sequence_index = int(str(redesign_row["sequence_index"]).strip())
            designed_residue = designed_residue_map.get((chain_id, sequence_index))
            if designed_residue is None:
                raise ValueError(f"Unknown redesignable position for {candidate_id}: {chain_id}{sequence_index}")
            residue_id = str(redesign_row["residue_id"]).strip()
            seed_amino_acid = str(redesign_row["seed_amino_acid"]).strip()
            if residue_id != designed_residue.residue_id:
                raise ValueError(f"Residue-id mismatch for {candidate_id} at {chain_id}{sequence_index}")
            if seed_amino_acid != designed_residue.seed_amino_acid:
                raise ValueError(f"Seed amino-acid mismatch for {candidate_id} at {residue_id}")
            redesignable_positions.append(
                Round2RedesignPosition(
                    chain_id=chain_id,
                    sequence_index=sequence_index,
                    residue_id=residue_id,
                    seed_amino_acid=seed_amino_acid,
                    canonical_family_position=int(str(redesign_row["canonical_family_position"]).strip()),
                    am1_mature_position=int(str(redesign_row["am1_mature_position"]).strip()),
                )
            )

        chain_assignment_payload = _load_single_jsonl_mapping(chain_assignment_path)
        chain_assignment_entry = chain_assignment_payload.get(candidate_id)
        if not isinstance(chain_assignment_entry, list) or len(chain_assignment_entry) != 2:
            raise ValueError(f"Unexpected chain assignment payload for {candidate_id}")
        parsed_designed_chains = tuple(str(chain_id).strip() for chain_id in chain_assignment_entry[0])
        parsed_fixed_context_chains = tuple(str(chain_id).strip() for chain_id in chain_assignment_entry[1])
        if parsed_designed_chains != designed_chains:
            raise ValueError(f"Chain-assignment designed chains mismatch for {candidate_id}")
        if parsed_fixed_context_chains != fixed_context_chains:
            raise ValueError(f"Chain-assignment fixed-context chains mismatch for {candidate_id}")

        redesignable_positions_by_chain: dict[str, set[int]] = defaultdict(set)
        for position in redesignable_positions:
            redesignable_positions_by_chain[position.chain_id].add(position.sequence_index)
        expected_fixed_positions = {
            chain_id: [
                residue.sequence_index
                for residue in designed_residues
                if residue.chain_id == chain_id and residue.sequence_index not in redesignable_positions_by_chain[chain_id]
            ]
            for chain_id in designed_chains
        }
        fixed_positions_payload = _load_single_jsonl_mapping(fixed_positions_path)
        fixed_positions_entry = fixed_positions_payload.get(candidate_id)
        if not isinstance(fixed_positions_entry, dict):
            raise ValueError(f"Unexpected fixed-positions payload for {candidate_id}")
        parsed_fixed_positions = {
            str(chain_id).strip(): [int(position) for position in positions]
            for chain_id, positions in fixed_positions_entry.items()
        }
        if parsed_fixed_positions != expected_fixed_positions:
            raise ValueError(f"Fixed positions did not match redesignable-position complement for {candidate_id}")

        redesignable_residue_ids = ",".join(position.residue_id for position in redesignable_positions)
        redesignable_canonical_positions = ",".join(
            str(position.canonical_family_position) for position in redesignable_positions
        )
        redesignable_am1_positions = ",".join(str(position.am1_mature_position) for position in redesignable_positions)
        if redesignable_residue_ids != str(row["redesignable_residue_ids"]).strip():
            raise ValueError(f"Redesignable residue ids mismatch for {candidate_id}")
        if redesignable_canonical_positions != str(row["redesignable_canonical_positions"]).strip():
            raise ValueError(f"Redesignable canonical positions mismatch for {candidate_id}")
        if redesignable_am1_positions != str(row["redesignable_am1_positions"]).strip():
            raise ValueError(f"Redesignable AM1 positions mismatch for {candidate_id}")

        campaigns.append(
            ProteinMPNNRound2Campaign(
                seed_rank=int(str(row["seed_rank"]).strip()),
                candidate_id=candidate_id,
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                design_set_name=str(config_seed.get("design_set_name", "")).strip(),
                designed_chains=designed_chains,
                fixed_context_chains=fixed_context_chains,
                pdb_path=pdb_path,
                chain_assignment_path=chain_assignment_path,
                fixed_positions_path=fixed_positions_path,
                output_dir=output_root / candidate_id,
                seed_scaffold_sequence=seed_scaffold_sequence,
                designed_residues=designed_residues,
                redesignable_positions=tuple(redesignable_positions),
            )
        )
    return tuple(campaigns)


def build_round2_parse_command(
    campaign: ProteinMPNNRound2Campaign,
    proteinmpnn_root: Path,
    python_executable: str | None = None,
) -> tuple[str, ...]:
    """Build the parse_multiple_chains.py command for one round-2 campaign."""
    return build_parse_multiple_chains_command(
        proteinmpnn_root=proteinmpnn_root,
        input_path=campaign.pdb_path.parent,
        output_path=campaign.output_dir / "parsed_pdbs.jsonl",
        python_executable=python_executable,
    )


def build_round2_run_command(
    campaign: ProteinMPNNRound2Campaign,
    proteinmpnn_root: Path,
    parsed_jsonl_path: Path,
    python_executable: str | None = None,
) -> tuple[str, ...]:
    """Build the deterministic ProteinMPNN round-2 smoke command."""
    return build_proteinmpnn_run_command(
        proteinmpnn_root=proteinmpnn_root,
        parsed_jsonl_path=parsed_jsonl_path,
        chain_assignment_path=campaign.chain_assignment_path,
        fixed_positions_path=campaign.fixed_positions_path,
        output_dir=campaign.output_dir,
        python_executable=python_executable,
    )


def summarize_round2_proteinmpnn_output(
    campaign: ProteinMPNNRound2Campaign,
    fasta_path: Path,
) -> tuple[ProteinMPNNRound2SummaryRow, tuple[ProteinMPNNRound2SequenceCatalogRow, ...]]:
    """Parse and summarize one round-2 ProteinMPNN FASTA output."""
    parsed_output = parse_proteinmpnn_fasta(fasta_path)
    native_record = parsed_output.native_record
    if native_record.name != campaign.candidate_id:
        raise ValueError(
            f"ProteinMPNN FASTA name {native_record.name} did not match round-2 campaign {campaign.candidate_id}"
        )
    if native_record.designed_chains != campaign.designed_chains:
        raise ValueError(
            f"Designed chains {native_record.designed_chains} did not match round-2 manifest {campaign.designed_chains}"
        )
    if native_record.fixed_chains != campaign.fixed_context_chains:
        raise ValueError(
            f"Fixed chains {native_record.fixed_chains} did not match round-2 manifest {campaign.fixed_context_chains}"
        )
    if native_record.sequence != campaign.seed_scaffold_sequence:
        raise ValueError(
            f"Native seed sequence for {campaign.candidate_id} did not match exported round-2 scaffold"
        )

    output_fasta_path = _display_path(fasta_path)
    stripped_sequences = [record.sequence.replace("/", "") for record in parsed_output.generated_records]
    mean_identity, min_identity, max_identity = _pairwise_identity_summary(stripped_sequences)

    catalog_rows: list[ProteinMPNNRound2SequenceCatalogRow] = []
    mutation_counter: Counter[int] = Counter()
    for record in parsed_output.generated_records:
        mutation_count, mutation_string, mutated_residue_ids = _summarize_mutations_vs_seed(
            campaign=campaign,
            designed_sequence=record.sequence,
        )
        mutation_counter[mutation_count] += 1
        catalog_rows.append(
            ProteinMPNNRound2SequenceCatalogRow(
                seed_rank=campaign.seed_rank,
                candidate_id=campaign.candidate_id,
                campaign_id=campaign.campaign_id,
                backbone_id=campaign.backbone_id,
                topology_class=campaign.topology_class,
                design_set_name=campaign.design_set_name,
                designed_chains=campaign.designed_chains_text,
                temperature=record.temperature,
                sample_number=record.sample_number,
                score=record.score,
                global_score=record.global_score,
                seq_recovery=record.seq_recovery,
                sequence_id=f"{campaign.candidate_id}_T{record.temperature:g}_sample_{record.sample_number:02d}",
                designed_sequence=record.sequence,
                sequence_length=len(record.sequence.replace("/", "")),
                mutation_count_vs_seed=mutation_count,
                mutation_string_vs_seed=mutation_string,
                mutated_residue_ids=mutated_residue_ids,
                output_fasta_path=output_fasta_path,
            )
        )

    generated_sequence_count = len(parsed_output.generated_records)
    unique_sequence_count = len({record.sequence for record in parsed_output.generated_records})
    sequence_length = len(campaign.seed_scaffold_sequence.replace("/", ""))
    mutation_count_distribution = ",".join(
        f"{mutation_count}:{mutation_counter[mutation_count]}"
        for mutation_count in sorted(mutation_counter)
    )
    summary_row = ProteinMPNNRound2SummaryRow(
        seed_rank=campaign.seed_rank,
        candidate_id=campaign.candidate_id,
        campaign_id=campaign.campaign_id,
        backbone_id=campaign.backbone_id,
        topology_class=campaign.topology_class,
        design_set_name=campaign.design_set_name,
        designed_chains=campaign.designed_chains_text,
        fixed_context_chains=campaign.fixed_context_chains_text,
        seed_scaffold_path=_display_path(campaign.pdb_path),
        seed_scaffold_sequence=campaign.seed_scaffold_sequence,
        redesignable_position_count=len(campaign.redesignable_positions),
        redesignable_residue_ids=campaign.redesignable_residue_ids_text,
        redesignable_canonical_positions=campaign.redesignable_canonical_positions_text,
        redesignable_am1_positions=campaign.redesignable_am1_positions_text,
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
        mutation_count_distribution=mutation_count_distribution,
        output_fasta_path=output_fasta_path,
    )
    return summary_row, tuple(catalog_rows)


def render_proteinmpnn_round2_smoke_markdown(
    proteinmpnn_root: Path,
    summary_rows: Sequence[ProteinMPNNRound2SummaryRow],
) -> str:
    """Render the Phase 7B ProteinMPNN round-2 smoke markdown report."""
    lines = [
        "# ProteinMPNN Round-2 Smoke",
        "",
        "Deterministic Phase 7B ProteinMPNN smoke generation across the narrow round-2 redesign campaigns.",
        "",
        f"- ProteinMPNN root: `{proteinmpnn_root}`",
        f"- num_seq_per_target: `{SMOKE_NUM_SEQ_PER_TARGET}`",
        f"- sampling_temp: `{SMOKE_SAMPLING_TEMP}`",
        f"- seed: `{SMOKE_SEED}`",
        f"- batch_size: `{SMOKE_BATCH_SIZE}`",
        "",
        "| seed_scaffold | redesignable_positions | unique_sequences | mutation_count_distribution |",
        "| --- | --- | ---: | --- |",
    ]
    for row in summary_rows:
        lines.append(
            "| "
            f"{row.candidate_id} | "
            f"{row.redesignable_canonical_positions or '-'} | "
            f"{row.unique_sequence_count}/{row.generated_sequence_count} | "
            f"{row.mutation_count_distribution or '-'} |"
        )
    lines.append("")
    for row in summary_rows:
        lines.extend(
            [
                f"## {row.candidate_id}",
                "",
                f"- Seed scaffold: `{row.candidate_id}` from `{row.seed_scaffold_path}`",
                f"- Designed-chain seed sequence: `{row.seed_scaffold_sequence}`",
                f"- Redesignable positions: residue ids `{row.redesignable_residue_ids or '-'}`; canonical `{row.redesignable_canonical_positions or '-'}`; AM1 `{row.redesignable_am1_positions or '-'}`",
                f"- Unique sequences: `{row.unique_sequence_count}` of `{row.generated_sequence_count}` generated",
                f"- Mutation-count distribution vs seed scaffold: `{row.mutation_count_distribution or '-'}`",
                "",
            ]
        )
    return "\n".join(lines)


def run_proteinmpnn_round2_smoke(
    proteinmpnn_root: Path,
    redesign_round2_config_path: Path,
    round2_seed_table_path: Path,
    round2_position_table_path: Path,
    campaigns_dir: Path,
    output_root: Path,
    summary_path: Path,
    sequence_catalog_path: Path,
    report_path: Path,
) -> tuple[tuple[ProteinMPNNRound2SummaryRow, ...], tuple[ProteinMPNNRound2SequenceCatalogRow, ...]]:
    """Run the focused round-2 ProteinMPNN smoke workflow and write summary artifacts."""
    root = validate_proteinmpnn_root(proteinmpnn_root)
    campaigns = discover_round2_proteinmpnn_campaigns(
        redesign_round2_config_path=redesign_round2_config_path,
        round2_seed_table_path=round2_seed_table_path,
        round2_position_table_path=round2_position_table_path,
        campaigns_dir=campaigns_dir,
        output_root=output_root,
    )
    if not campaigns:
        raise ValueError(f"No round-2 campaigns discovered under {campaigns_dir}")

    summary_rows: list[ProteinMPNNRound2SummaryRow] = []
    catalog_rows: list[ProteinMPNNRound2SequenceCatalogRow] = []
    for campaign in campaigns:
        LOGGER.info("Running ProteinMPNN round-2 smoke workflow for %s", campaign.candidate_id)
        _reset_output_directory(campaign.output_dir)
        parsed_jsonl_path = campaign.output_dir / "parsed_pdbs.jsonl"

        parse_command = build_round2_parse_command(
            campaign=campaign,
            proteinmpnn_root=root,
        )
        atomic_write_text(campaign.output_dir / "parse_command.txt", shlex.join(parse_command) + "\n")
        _run_logged_command(
            command=parse_command,
            cwd=root,
            stdout_path=campaign.output_dir / "parse_stdout.txt",
            stderr_path=campaign.output_dir / "parse_stderr.txt",
        )
        require_path(parsed_jsonl_path)

        run_command = build_round2_run_command(
            campaign=campaign,
            proteinmpnn_root=root,
            parsed_jsonl_path=parsed_jsonl_path,
        )
        atomic_write_text(campaign.output_dir / "run_command.txt", shlex.join(run_command) + "\n")
        _run_logged_command(
            command=run_command,
            cwd=root,
            stdout_path=campaign.output_dir / "run_stdout.txt",
            stderr_path=campaign.output_dir / "run_stderr.txt",
        )

        fasta_path = campaign.output_dir / "seqs" / f"{campaign.candidate_id}.fa"
        summary_row, campaign_catalog_rows = summarize_round2_proteinmpnn_output(
            campaign=campaign,
            fasta_path=fasta_path,
        )
        summary_rows.append(summary_row)
        catalog_rows.extend(campaign_catalog_rows)

    write_csv_rows(summary_path, summary_rows)
    write_csv_rows(sequence_catalog_path, catalog_rows)
    atomic_write_text(
        report_path,
        render_proteinmpnn_round2_smoke_markdown(
            proteinmpnn_root=root,
            summary_rows=summary_rows,
        ),
    )
    return tuple(summary_rows), tuple(catalog_rows)


def _summarize_mutations_vs_seed(
    campaign: ProteinMPNNRound2Campaign,
    designed_sequence: str,
) -> tuple[int, str, str]:
    seed_parts = _split_designed_sequence(campaign.seed_scaffold_sequence, campaign.designed_chains)
    designed_parts = _split_designed_sequence(designed_sequence, campaign.designed_chains)
    redesignable_positions = {
        (position.chain_id, position.sequence_index)
        for position in campaign.redesignable_positions
    }
    designed_residue_map = {
        (residue.chain_id, residue.sequence_index): residue
        for residue in campaign.designed_residues
    }

    mutation_tokens: list[str] = []
    mutated_residue_ids: list[str] = []
    for chain_id, seed_part, designed_part in zip(campaign.designed_chains, seed_parts, designed_parts):
        if len(seed_part) != len(designed_part):
            raise ValueError(
                f"Designed sequence length mismatch for {campaign.candidate_id} chain {chain_id}: "
                f"{len(designed_part)} vs seed {len(seed_part)}"
            )
        for sequence_index, (seed_aa, designed_aa) in enumerate(zip(seed_part, designed_part), start=1):
            if seed_aa == designed_aa:
                continue
            residue = designed_residue_map[(chain_id, sequence_index)]
            if (chain_id, sequence_index) not in redesignable_positions:
                raise ValueError(
                    f"Unexpected mutation outside round-2 redesignable positions for {campaign.candidate_id}: "
                    f"{residue.residue_id}"
                )
            mutation_tokens.append(f"{residue.residue_id}:{seed_aa}>{designed_aa}")
            mutated_residue_ids.append(residue.residue_id)

    return len(mutation_tokens), ",".join(mutation_tokens), ",".join(mutated_residue_ids)


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
        raise RuntimeError(f"Command failed with exit code {completed.returncode}: {shlex.join(command)}")
