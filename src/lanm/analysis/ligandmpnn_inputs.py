"""Phase 4A LigandMPNN input preparation for shortlisted ProteinMPNN candidates."""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from lanm.configuration import load_project_config
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows, write_yaml
from lanm.models import (
    DesignBackboneManifestRow,
    DesignCampaignManifestRow,
    DesignCampaignPositionRow,
    LigandMPNNInputManifestRow,
    LigandMPNNRedesignPositionRow,
    ProteinMPNNShortlistRow,
)
from lanm.paths import REPO_ROOT
from lanm.structure.atoms import read_atom_records
from lanm.structure.geometry import WATER_RESIDUES, collect_atoms_within_cutoff
from lanm.structure.pdb import render_selected_structure_pdb, select_preferred_atom_conformers
from lanm.structure.residues import AMINO_ACID_CODES
from lanm.structure.templates import parse_cif_atom_records

_SUPPORTED_BACKBONES: dict[str, dict[str, object]] = {
    "am1_mex_8fns_chain_a": {
        "structure_id": "8FNS",
        "preserved_chains": ("A",),
        "preserved_metal_identity": "ND",
    },
    "hans_pocket_8fnr_chain_a": {
        "structure_id": "8FNR",
        "preserved_chains": ("A",),
        "preserved_metal_identity": "DY",
    },
    "hans_interface_8fnr_a_b_c_d": {
        "structure_id": "8FNR",
        "preserved_chains": ("A", "B", "C", "D"),
        "preserved_metal_identity": "DY",
    },
}


@dataclass(frozen=True, slots=True)
class ShortlistedLigandMPNNCandidate:
    shortlist_rank: int
    candidate_id: str
    campaign_id: str
    backbone_id: str
    design_set_name: str
    designed_sequence: str
    mutation_count: int
    mutation_string: str
    canonical_family_positions_mutated: str
    am1_mature_positions_mutated: str
    source_structure_id: str
    source_kind: str
    source_path: Path
    preserved_chains: tuple[str, ...]
    designed_chains: tuple[str, ...]
    fixed_context_chains: tuple[str, ...]
    preserved_metal_identity: str


@dataclass(frozen=True, slots=True)
class LigandMPNNPDBExport:
    pdb_text: str
    preserved_metal_site_count: int
    preserved_solvent_residue_count: int


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _parse_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def _parse_optional_int(value: str) -> int | None:
    stripped = value.strip()
    return int(stripped) if stripped else None


def _parse_csv_list(value: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in value.split(",") if item.strip())


def _load_shortlist_rows(path: Path) -> tuple[ProteinMPNNShortlistRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            ProteinMPNNShortlistRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                design_set_name=str(row["design_set_name"]).strip(),
                designed_sequence=str(row["designed_sequence"]).strip(),
                representative_sequence_id=str(row["representative_sequence_id"]).strip(),
                occurrence_count=int(str(row["occurrence_count"]).strip()),
                campaign_rank=int(str(row["campaign_rank"]).strip()),
                temperature=float(str(row["temperature"]).strip()),
                best_score=float(str(row["best_score"]).strip()),
                best_global_score=float(str(row["best_global_score"]).strip()),
                best_seq_recovery=float(str(row["best_seq_recovery"]).strip()),
                mutation_count=int(str(row["mutation_count"]).strip()),
                mutation_string=str(row["mutation_string"]).strip(),
                canonical_family_positions_mutated=str(row["canonical_family_positions_mutated"]).strip(),
                am1_mature_positions_mutated=str(row["am1_mature_positions_mutated"]).strip(),
                includes_second_sphere_position=_parse_bool(str(row["includes_second_sphere_position"])),
                includes_interface_position=_parse_bool(str(row["includes_interface_position"])),
                retention_reason=str(row["retention_reason"]).strip(),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.shortlist_rank))


def _load_design_backbone_manifest_rows(path: Path) -> tuple[DesignBackboneManifestRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            DesignBackboneManifestRow(
                backbone_id=str(row["backbone_id"]).strip(),
                structure_id=str(row["structure_id"]).strip(),
                source_kind=str(row["source_kind"]).strip(),
                source_path=str(row["source_path"]).strip(),
                experimental_method=str(row["experimental_method"]).strip(),
                selected_chains=str(row["selected_chains"]).strip(),
                design_chains=str(row["design_chains"]).strip(),
                fixed_context_chains=str(row["fixed_context_chains"]).strip(),
                selected_residue_count=int(str(row["selected_residue_count"]).strip()),
                design_chain_residue_count=int(str(row["design_chain_residue_count"]).strip()),
                design_chain_metal_site_count=int(str(row["design_chain_metal_site_count"]).strip()),
                representative_note=str(row["representative_note"]).strip(),
                backbone_pdb_path=str(row["backbone_pdb_path"]).strip(),
            )
            for row in reader
        ]
    return tuple(rows)


def _load_design_campaign_manifest_rows(path: Path) -> tuple[DesignCampaignManifestRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            DesignCampaignManifestRow(
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                design_set_name=str(row["design_set_name"]).strip(),
                designed_chains=str(row["designed_chains"]).strip(),
                fixed_context_chains=str(row["fixed_context_chains"]).strip(),
                exported_pdb_path=str(row["exported_pdb_path"]).strip(),
                chain_assignment_path=str(row["chain_assignment_path"]).strip(),
                fixed_positions_path=str(row["fixed_positions_path"]).strip(),
                designable_residue_count=int(str(row["designable_residue_count"]).strip()),
                fixed_residue_count=int(str(row["fixed_residue_count"]).strip()),
                designable_canonical_positions=str(row["designable_canonical_positions"]).strip(),
                designable_am1_positions=str(row["designable_am1_positions"]).strip(),
            )
            for row in reader
        ]
    return tuple(rows)


def _ligandmpnn_residue_id(chain_id: str, residue_seq: int, insertion_code: str) -> str:
    return f"{chain_id}{residue_seq}{insertion_code}".strip()


def _split_designed_sequence(sequence: str, designed_chains: tuple[str, ...]) -> tuple[str, ...]:
    if len(designed_chains) == 1:
        return (sequence.replace("/", ""),)
    parts = tuple(part.strip() for part in sequence.split("/"))
    if len(parts) != len(designed_chains):
        raise ValueError(
            f"Expected {len(designed_chains)} designed-chain sequence parts, found {len(parts)} in {sequence!r}"
        )
    return parts


def _group_position_rows_by_chain(
    position_rows: tuple[DesignCampaignPositionRow, ...],
    designed_chains: tuple[str, ...],
) -> dict[str, tuple[DesignCampaignPositionRow, ...]]:
    rows_by_chain: dict[str, list[DesignCampaignPositionRow]] = defaultdict(list)
    for row in position_rows:
        if row.chain_id not in designed_chains:
            continue
        rows_by_chain[row.chain_id].append(row)
    grouped: dict[str, tuple[DesignCampaignPositionRow, ...]] = {}
    for chain_id in designed_chains:
        chain_rows = tuple(sorted(rows_by_chain[chain_id], key=lambda row: row.sequence_index))
        if not chain_rows:
            raise ValueError(f"No design-campaign positions found for designed chain {chain_id}")
        grouped[chain_id] = chain_rows
    return grouped


def _native_sequence(chain_rows: tuple[DesignCampaignPositionRow, ...]) -> str:
    return "".join(AMINO_ACID_CODES.get(row.residue_name, "X") for row in chain_rows)


def _load_structure_atoms(
    *,
    source_kind: str,
    source_path: Path,
    structure_id: str,
) -> tuple:
    require_path(source_path)
    if source_kind == "cif":
        return tuple(parse_cif_atom_records(source_path, structure_id=structure_id))
    if source_kind == "csv_atom_table":
        return tuple(read_atom_records(source_path, structure_id=structure_id))
    raise ValueError(f"Unsupported structure source_kind: {source_kind}")


def discover_shortlisted_ligandmpnn_candidates(
    *,
    shortlist_path: Path,
    campaign_manifest_path: Path,
    backbone_manifest_path: Path,
) -> tuple[ShortlistedLigandMPNNCandidate, ...]:
    """Resolve shortlisted candidates onto deterministic LigandMPNN structural contexts."""
    shortlist_rows = _load_shortlist_rows(shortlist_path)
    campaign_map = {
        row.campaign_id: row
        for row in _load_design_campaign_manifest_rows(campaign_manifest_path)
    }
    backbone_map = {
        row.backbone_id: row
        for row in _load_design_backbone_manifest_rows(backbone_manifest_path)
    }
    candidates: list[ShortlistedLigandMPNNCandidate] = []
    for shortlist_row in shortlist_rows:
        campaign_row = campaign_map.get(shortlist_row.campaign_id)
        if campaign_row is None:
            raise ValueError(f"No design_campaign_manifest.csv entry for {shortlist_row.campaign_id}")
        backbone_row = backbone_map.get(shortlist_row.backbone_id)
        if backbone_row is None:
            raise ValueError(f"No design_backbone_manifest.csv entry for {shortlist_row.backbone_id}")
        if campaign_row.backbone_id != shortlist_row.backbone_id:
            raise ValueError(
                f"Shortlist backbone mismatch for {shortlist_row.candidate_id}: "
                f"{shortlist_row.backbone_id} vs {campaign_row.backbone_id}"
            )
        context = _SUPPORTED_BACKBONES.get(shortlist_row.backbone_id)
        if context is None:
            raise ValueError(f"Unsupported LigandMPNN backbone: {shortlist_row.backbone_id}")
        preserved_chains = tuple(context["preserved_chains"])
        if backbone_row.structure_id != str(context["structure_id"]):
            raise ValueError(
                f"Backbone {shortlist_row.backbone_id} expected structure {context['structure_id']}, "
                f"found {backbone_row.structure_id}"
            )
        if _parse_csv_list(backbone_row.selected_chains) != preserved_chains:
            raise ValueError(
                f"Backbone {shortlist_row.backbone_id} selected chains "
                f"{_parse_csv_list(backbone_row.selected_chains)} did not match {preserved_chains}"
            )
        designed_chains = _parse_csv_list(campaign_row.designed_chains)
        if _parse_csv_list(backbone_row.design_chains) != designed_chains:
            raise ValueError(
                f"Backbone {shortlist_row.backbone_id} design chains "
                f"{_parse_csv_list(backbone_row.design_chains)} did not match campaign "
                f"{designed_chains}"
            )
        fixed_context_chains = _parse_csv_list(campaign_row.fixed_context_chains)
        candidates.append(
            ShortlistedLigandMPNNCandidate(
                shortlist_rank=shortlist_row.shortlist_rank,
                candidate_id=shortlist_row.candidate_id,
                campaign_id=shortlist_row.campaign_id,
                backbone_id=shortlist_row.backbone_id,
                design_set_name=shortlist_row.design_set_name,
                designed_sequence=shortlist_row.designed_sequence,
                mutation_count=shortlist_row.mutation_count,
                mutation_string=shortlist_row.mutation_string,
                canonical_family_positions_mutated=shortlist_row.canonical_family_positions_mutated,
                am1_mature_positions_mutated=shortlist_row.am1_mature_positions_mutated,
                source_structure_id=backbone_row.structure_id,
                source_kind=backbone_row.source_kind,
                source_path=REPO_ROOT / backbone_row.source_path,
                preserved_chains=preserved_chains,
                designed_chains=designed_chains,
                fixed_context_chains=fixed_context_chains,
                preserved_metal_identity=str(context["preserved_metal_identity"]),
            )
        )
    return tuple(candidates)


def build_candidate_redesign_rows(
    *,
    candidate: ShortlistedLigandMPNNCandidate,
    position_rows: tuple[DesignCampaignPositionRow, ...],
) -> tuple[LigandMPNNRedesignPositionRow, ...]:
    """Map actual shortlist mutations back to chain and residue identifiers."""
    position_rows_by_chain = _group_position_rows_by_chain(position_rows, candidate.designed_chains)
    designed_parts = _split_designed_sequence(candidate.designed_sequence, candidate.designed_chains)

    redesign_rows: list[LigandMPNNRedesignPositionRow] = []
    canonical_positions: list[int] = []
    am1_positions: list[int] = []
    mutation_tokens: list[str] = []
    for chain_id, designed_part in zip(candidate.designed_chains, designed_parts):
        chain_rows = position_rows_by_chain[chain_id]
        native_part = _native_sequence(chain_rows)
        if len(native_part) != len(designed_part):
            raise ValueError(
                f"Designed sequence length {len(designed_part)} for {candidate.candidate_id} did not match "
                f"native chain length {len(native_part)} on chain {chain_id}"
            )
        for native_aa, designed_aa, row in zip(native_part, designed_part, chain_rows):
            if native_aa == designed_aa:
                continue
            mutation_token = f"{native_aa}{row.sequence_index}{designed_aa}"
            mutation_tokens.append(mutation_token)
            if row.canonical_family_position is not None and row.canonical_family_position not in canonical_positions:
                canonical_positions.append(row.canonical_family_position)
            if row.am1_mature_position is not None and row.am1_mature_position not in am1_positions:
                am1_positions.append(row.am1_mature_position)
            redesign_rows.append(
                LigandMPNNRedesignPositionRow(
                    shortlist_rank=candidate.shortlist_rank,
                    candidate_id=candidate.candidate_id,
                    campaign_id=candidate.campaign_id,
                    backbone_id=candidate.backbone_id,
                    chain_id=row.chain_id,
                    sequence_index=row.sequence_index,
                    residue_seq=row.residue_seq,
                    insertion_code=row.insertion_code,
                    ligandmpnn_residue_id=_ligandmpnn_residue_id(
                        row.chain_id,
                        row.residue_seq,
                        row.insertion_code,
                    ),
                    residue_name=row.residue_name,
                    native_amino_acid=native_aa,
                    designed_amino_acid=designed_aa,
                    mutation_token=mutation_token,
                    canonical_family_position=row.canonical_family_position,
                    am1_mature_position=row.am1_mature_position,
                )
            )

    mutation_string = ",".join(mutation_tokens)
    if len(redesign_rows) != candidate.mutation_count:
        raise ValueError(
            f"Mutation count mismatch for {candidate.candidate_id}: "
            f"{len(redesign_rows)} vs shortlist {candidate.mutation_count}"
        )
    if mutation_string != candidate.mutation_string:
        raise ValueError(
            f"Mutation string mismatch for {candidate.candidate_id}: "
            f"{mutation_string!r} vs shortlist {candidate.mutation_string!r}"
        )
    if ",".join(str(position) for position in canonical_positions) != candidate.canonical_family_positions_mutated:
        raise ValueError(
            f"Canonical mutation mismatch for {candidate.candidate_id}: "
            f"{canonical_positions} vs {candidate.canonical_family_positions_mutated!r}"
        )
    if ",".join(str(position) for position in am1_positions) != candidate.am1_mature_positions_mutated:
        raise ValueError(
            f"AM1 mutation mismatch for {candidate.candidate_id}: "
            f"{am1_positions} vs {candidate.am1_mature_positions_mutated!r}"
        )
    return tuple(redesign_rows)


def build_ligandmpnn_pdb_export(
    *,
    atoms: tuple,
    preserved_chains: tuple[str, ...],
    preserved_metal_identity: str,
    nearby_solvent_cutoff_A: float,
) -> LigandMPNNPDBExport:
    """Render a PDB export that keeps the protein plus selected metals and nearby solvent."""
    selected_atoms = tuple(select_preferred_atom_conformers(atoms))
    metal_atoms = tuple(
        atom
        for atom in selected_atoms
        if atom.record_type == "HETATM"
        and atom.chain_id in preserved_chains
        and atom.element_upper == preserved_metal_identity
    )
    if not metal_atoms:
        raise ValueError(
            f"No preserved metal atoms with identity {preserved_metal_identity} found on chains {preserved_chains}"
        )
    solvent_residue_keys = {
        item.atom.residue_key
        for metal_atom in metal_atoms
        for item in collect_atoms_within_cutoff(metal_atom, selected_atoms, nearby_solvent_cutoff_A)
        if item.atom.record_type == "HETATM"
        and item.atom.residue_name in WATER_RESIDUES
        and not item.atom.is_metal
    }
    included_het_residue_keys = frozenset(
        {atom.residue_key for atom in metal_atoms}
        | solvent_residue_keys
    )
    return LigandMPNNPDBExport(
        pdb_text=render_selected_structure_pdb(
            selected_atoms,
            selected_chain_ids=preserved_chains,
            included_het_residue_keys=included_het_residue_keys,
        ),
        preserved_metal_site_count=len(metal_atoms),
        preserved_solvent_residue_count=len(solvent_residue_keys),
    )


def _render_redesigned_residues_text(
    redesign_rows: tuple[LigandMPNNRedesignPositionRow, ...],
) -> str:
    if not redesign_rows:
        return ""
    return " ".join(row.ligandmpnn_residue_id for row in redesign_rows) + "\n"


def _render_report(
    *,
    manifest_rows: tuple[LigandMPNNInputManifestRow, ...],
    redesign_rows_by_candidate: dict[str, tuple[LigandMPNNRedesignPositionRow, ...]],
    nearby_solvent_cutoff_A: float,
) -> str:
    lines = [
        "# LigandMPNN Inputs",
        "",
        "Phase 4A deterministic LigandMPNN-ready input preparation for shortlisted ProteinMPNN candidates.",
        "",
        f"- Nearby solvent retained within `{nearby_solvent_cutoff_A:.1f} A` of preserved metal HETATM sites.",
        "- Redesigned residue identifiers use LigandMPNN chain-plus-residue numbering with optional insertion code.",
        "",
        "| candidate_id | source_backbone | preserved_chains | preserved_metal_identity | redesigned_residue_identifiers |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in manifest_rows:
        redesign_ids = ",".join(
            redesign_row.ligandmpnn_residue_id
            for redesign_row in redesign_rows_by_candidate[row.candidate_id]
        ) or "native"
        lines.append(
            f"| {row.candidate_id} | {row.backbone_id} | {row.preserved_chains} | "
            f"{row.preserved_metal_identity} | {redesign_ids} |"
        )
    lines.append("")
    return "\n".join(lines)


def _build_config_payload(
    *,
    manifest_rows: tuple[LigandMPNNInputManifestRow, ...],
    redesign_rows_by_candidate: dict[str, tuple[LigandMPNNRedesignPositionRow, ...]],
    nearby_solvent_cutoff_A: float,
) -> dict[str, Any]:
    return {
        "version": 1,
        "phase": "4A",
        "source_artifacts": {
            "proteinmpnn_shortlist": "results/tables/proteinmpnn_shortlist.csv",
            "design_backbone_manifest": "results/tables/design_backbone_manifest.csv",
            "design_campaign_manifest": "results/tables/design_campaign_manifest.csv",
            "design_campaign_positions": "results/tables/design_campaign_positions.csv",
            "project_config": "config/project.yaml",
        },
        "ligandmpnn": {
            "nearby_solvent_cutoff_A": nearby_solvent_cutoff_A,
            "redesigned_residue_format": "Chain plus residue number with optional insertion code, e.g. A35 or B42D.",
        },
        "candidates": {
            row.candidate_id: {
                "shortlist_rank": row.shortlist_rank,
                "campaign_id": row.campaign_id,
                "backbone_id": row.backbone_id,
                "source_structure_id": row.source_structure_id,
                "source_kind": row.source_kind,
                "source_path": row.source_path,
                "preserved_chains": _parse_csv_list(row.preserved_chains),
                "designed_chains": _parse_csv_list(row.designed_chains),
                "fixed_context_chains": _parse_csv_list(row.fixed_context_chains),
                "preserved_metal_identity": row.preserved_metal_identity,
                "preserved_metal_site_count": row.preserved_metal_site_count,
                "preserved_solvent_residue_count": row.preserved_solvent_residue_count,
                "redesigned_residue_count": row.redesigned_residue_count,
                "redesigned_residue_ids": [
                    redesign_row.ligandmpnn_residue_id
                    for redesign_row in redesign_rows_by_candidate[row.candidate_id]
                ],
                "pdb_path": row.pdb_path,
                "redesigned_residues_path": row.redesigned_residues_path,
            }
            for row in manifest_rows
        },
    }


def prepare_ligandmpnn_inputs(
    *,
    shortlist_path: Path,
    campaign_manifest_path: Path,
    backbone_manifest_path: Path,
    design_campaign_positions_path: Path,
    manifest_path: Path,
    redesign_positions_path: Path,
    report_path: Path,
    input_root: Path,
    config_path: Path,
) -> tuple[tuple[LigandMPNNInputManifestRow, ...], tuple[LigandMPNNRedesignPositionRow, ...]]:
    """Prepare deterministic LigandMPNN inputs for the shortlisted ProteinMPNN candidates."""
    candidates = discover_shortlisted_ligandmpnn_candidates(
        shortlist_path=shortlist_path,
        campaign_manifest_path=campaign_manifest_path,
        backbone_manifest_path=backbone_manifest_path,
    )
    project_config = load_project_config()
    position_rows = tuple(
        sorted(
            (
                DesignCampaignPositionRow(
                    campaign_id=row.campaign_id,
                    design_set_name=row.design_set_name,
                    backbone_id=row.backbone_id,
                    structure_id=row.structure_id,
                    chain_id=row.chain_id,
                    sequence_index=row.sequence_index,
                    residue_seq=row.residue_seq,
                    insertion_code=row.insertion_code,
                    residue_name=row.residue_name,
                    canonical_family_position=row.canonical_family_position,
                    am1_mature_position=row.am1_mature_position,
                    fixed_first_shell=row.fixed_first_shell,
                    protected_positions=row.protected_positions,
                    mutable_second_sphere=row.mutable_second_sphere,
                    mutable_interface=row.mutable_interface,
                    hard_fixed=row.hard_fixed,
                    designable=row.designable,
                    position_state=row.position_state,
                    position_reason=row.position_reason,
                    rationale=row.rationale,
                )
                for row in _load_design_campaign_position_rows(design_campaign_positions_path)
            ),
            key=lambda row: (row.campaign_id, row.chain_id, row.sequence_index),
        )
    )
    position_rows_by_campaign: dict[str, tuple[DesignCampaignPositionRow, ...]] = defaultdict(tuple)
    pending_rows_by_campaign: dict[str, list[DesignCampaignPositionRow]] = defaultdict(list)
    for row in position_rows:
        pending_rows_by_campaign[row.campaign_id].append(row)
    for campaign_id, rows in pending_rows_by_campaign.items():
        position_rows_by_campaign[campaign_id] = tuple(rows)

    atoms_by_source: dict[tuple[str, Path, str], tuple] = {}
    manifest_rows: list[LigandMPNNInputManifestRow] = []
    redesign_rows: list[LigandMPNNRedesignPositionRow] = []
    redesign_rows_by_candidate: dict[str, tuple[LigandMPNNRedesignPositionRow, ...]] = {}

    for candidate in candidates:
        campaign_rows = position_rows_by_campaign.get(candidate.campaign_id)
        if campaign_rows is None:
            raise ValueError(f"No design_campaign_positions.csv rows for {candidate.campaign_id}")
        candidate_redesign_rows = build_candidate_redesign_rows(
            candidate=candidate,
            position_rows=campaign_rows,
        )
        redesign_rows.extend(candidate_redesign_rows)
        redesign_rows_by_candidate[candidate.candidate_id] = candidate_redesign_rows

        source_key = (candidate.source_kind, candidate.source_path, candidate.source_structure_id)
        if source_key not in atoms_by_source:
            atoms_by_source[source_key] = _load_structure_atoms(
                source_kind=candidate.source_kind,
                source_path=candidate.source_path,
                structure_id=candidate.source_structure_id,
            )
        pdb_export = build_ligandmpnn_pdb_export(
            atoms=atoms_by_source[source_key],
            preserved_chains=candidate.preserved_chains,
            preserved_metal_identity=candidate.preserved_metal_identity,
            nearby_solvent_cutoff_A=project_config.second_sphere_cutoff_A,
        )

        candidate_dir = input_root / candidate.candidate_id
        pdb_path = candidate_dir / f"{candidate.candidate_id}.pdb"
        redesigned_residues_path = candidate_dir / "redesigned_residues.txt"
        atomic_write_text(pdb_path, pdb_export.pdb_text)
        atomic_write_text(
            redesigned_residues_path,
            _render_redesigned_residues_text(candidate_redesign_rows),
        )
        manifest_rows.append(
            LigandMPNNInputManifestRow(
                shortlist_rank=candidate.shortlist_rank,
                candidate_id=candidate.candidate_id,
                campaign_id=candidate.campaign_id,
                backbone_id=candidate.backbone_id,
                source_structure_id=candidate.source_structure_id,
                source_kind=candidate.source_kind,
                source_path=_display_path(candidate.source_path),
                preserved_chains=",".join(candidate.preserved_chains),
                designed_chains=",".join(candidate.designed_chains),
                fixed_context_chains=",".join(candidate.fixed_context_chains),
                preserved_metal_identity=candidate.preserved_metal_identity,
                preserved_metal_site_count=pdb_export.preserved_metal_site_count,
                preserved_solvent_residue_count=pdb_export.preserved_solvent_residue_count,
                redesigned_residue_count=len(candidate_redesign_rows),
                pdb_path=_display_path(pdb_path),
                redesigned_residues_path=_display_path(redesigned_residues_path),
            )
        )

    manifest_row_tuple = tuple(sorted(manifest_rows, key=lambda row: row.shortlist_rank))
    redesign_row_tuple = tuple(
        sorted(
            redesign_rows,
            key=lambda row: (
                row.shortlist_rank,
                row.chain_id,
                row.sequence_index,
                row.residue_seq,
                row.insertion_code,
            ),
        )
    )
    write_csv_rows(manifest_path, manifest_row_tuple)
    write_csv_rows(redesign_positions_path, redesign_row_tuple)
    atomic_write_text(
        report_path,
        _render_report(
            manifest_rows=manifest_row_tuple,
            redesign_rows_by_candidate=redesign_rows_by_candidate,
            nearby_solvent_cutoff_A=project_config.second_sphere_cutoff_A,
        ),
    )
    write_yaml(
        config_path,
        _build_config_payload(
            manifest_rows=manifest_row_tuple,
            redesign_rows_by_candidate=redesign_rows_by_candidate,
            nearby_solvent_cutoff_A=project_config.second_sphere_cutoff_A,
        ),
    )
    return manifest_row_tuple, redesign_row_tuple


def _load_design_campaign_position_rows(path: Path) -> tuple[DesignCampaignPositionRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            DesignCampaignPositionRow(
                campaign_id=str(row["campaign_id"]).strip(),
                design_set_name=str(row["design_set_name"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                structure_id=str(row["structure_id"]).strip(),
                chain_id=str(row["chain_id"]).strip(),
                sequence_index=int(str(row["sequence_index"]).strip()),
                residue_seq=int(str(row["residue_seq"]).strip()),
                insertion_code=str(row["insertion_code"]).strip(),
                residue_name=str(row["residue_name"]).strip(),
                canonical_family_position=_parse_optional_int(str(row["canonical_family_position"])),
                am1_mature_position=_parse_optional_int(str(row["am1_mature_position"])),
                fixed_first_shell=_parse_bool(str(row["fixed_first_shell"])),
                protected_positions=_parse_bool(str(row["protected_positions"])),
                mutable_second_sphere=_parse_bool(str(row["mutable_second_sphere"])),
                mutable_interface=_parse_bool(str(row["mutable_interface"])),
                hard_fixed=_parse_bool(str(row["hard_fixed"])),
                designable=_parse_bool(str(row["designable"])),
                position_state=str(row["position_state"]).strip(),
                position_reason=str(row["position_reason"]).strip(),
                rationale=str(row["rationale"]).strip(),
            )
            for row in reader
        ]
    return tuple(rows)
