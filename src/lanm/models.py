"""Typed data containers shared across Phase 1 workflows."""

from __future__ import annotations

from dataclasses import dataclass

METAL_ELEMENTS = {
    "LI", "NA", "K", "RB", "CS", "FR",
    "BE", "MG", "CA", "SR", "BA", "RA",
    "SC", "Y", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY",
    "HO", "ER", "TM", "YB", "LU", "AC", "TH", "PA", "U",
    "TI", "ZR", "HF",
    "V", "NB", "TA",
    "CR", "MO", "W",
    "MN", "TC", "RE",
    "FE", "RU", "OS",
    "CO", "RH", "IR",
    "NI", "PD", "PT",
    "CU", "AG", "AU",
    "ZN", "CD", "HG",
    "AL", "GA", "IN", "TL",
    "SN", "PB", "BI",
}
DONOR_ELEMENTS = {"O", "N", "S", "SE"}


@dataclass(frozen=True, slots=True)
class AtomRecord:
    structure_id: str
    record_type: str
    atom_serial: int
    atom_name: str
    alt_loc: str
    residue_name: str
    chain_id: str
    residue_seq: int
    insertion_code: str
    x: float
    y: float
    z: float
    occupancy: float | None
    b_factor: float | None
    element: str
    charge: str

    @property
    def element_upper(self) -> str:
        return self.element.upper()

    @property
    def is_metal(self) -> bool:
        return self.element_upper in METAL_ELEMENTS

    @property
    def residue_key(self) -> tuple[str, str, int, str]:
        return (self.chain_id, self.residue_name, self.residue_seq, self.insertion_code)

    @property
    def site_id(self) -> str:
        insertion = self.insertion_code or "-"
        return f"{self.structure_id}_{self.chain_id}_{self.residue_name}_{self.residue_seq}_{insertion}"


@dataclass(frozen=True, slots=True)
class FileSummary:
    relative_path: str
    file_type: str
    file_size_bytes: int
    row_count: int
    column_names: str
    chain_ids: str
    residue_ranges: str
    hetatm_count: int
    missing_fields: str


@dataclass(frozen=True, slots=True)
class NormalizationResult:
    normalized_paths: tuple[str, ...]
    skipped_sidecars: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class BriefExtractionResult:
    status: str
    output_path: str
    source_path: str
    note: str


@dataclass(frozen=True, slots=True)
class StructureManifestEntry:
    pdb_id: str
    title: str
    bound_metal: str
    organism: str
    method: str
    released: str
    rcsb_url: str
    cif_download_url: str

    @property
    def bound_metal_symbol(self) -> str:
        return self.bound_metal.split("(", 1)[0].strip().upper()


@dataclass(frozen=True, slots=True)
class FetchRecord:
    pdb_id: str
    bound_metal: str
    source_url: str
    local_path: str
    status: str
    size_bytes: int
    http_status: int
    message: str


@dataclass(frozen=True, slots=True)
class MetalSiteSummary:
    structure_id: str
    metal_site_id: str
    metal_element: str
    chain_id: str
    residue_name: str
    residue_seq: int
    atom_serial: int
    charge: str
    x: float
    y: float
    z: float
    donor_atom_count: int
    donor_residue_count: int
    residue_within_6a_count: int
    nearest_donor_distance_A: float | None
    farthest_donor_distance_A: float | None


@dataclass(frozen=True, slots=True)
class ShellAnnotation:
    structure_id: str
    metal_site_id: str
    metal_element: str
    shell_type: str
    chain_id: str
    residue_name: str
    residue_seq: int
    atom_name: str
    atom_serial: int
    element: str
    record_type: str
    distance_A: float
    is_water: bool


@dataclass(frozen=True, slots=True)
class TemplateChainSummary:
    structure_id: str
    source_kind: str
    input_path: str
    experimental_method: str
    chain_id: str
    residue_start: int
    residue_end: int
    residue_count: int
    sequence_length: int
    bound_metal_element: str
    bound_metal_site_count: int
    is_representative_chain: bool


@dataclass(frozen=True, slots=True)
class TemplateSiteSummaryRow:
    structure_id: str
    source_kind: str
    input_path: str
    experimental_method: str
    chain_id: str
    site_index: int
    site_label: str
    metal_site_id: str
    metal_element: str
    donor_atom_count: int
    donor_residue_count: int
    first_shell_polymer_count: int
    second_sphere_polymer_count: int
    solvent_count: int
    first_shell_residues: str
    second_sphere_residues: str
    solvent_residues: str


@dataclass(frozen=True, slots=True)
class CrossTemplateAlignmentRow:
    template_id: str
    chain_id: str
    template_residue_seq: int
    template_residue_name: str
    canonical_family_position: int | None
    am1_mature_position: int | None
    alignment_status: str


@dataclass(frozen=True, slots=True)
class ResidueRoleAssignment:
    template_id: str
    chain_id: str
    site_index: int | None
    site_label: str
    role: str
    residue_id: str
    residue_name: str
    residue_seq: int
    distance_A: float | None
    canonical_family_position: int | None
    am1_mature_position: int | None
    alignment_status: str
    note: str


@dataclass(frozen=True, slots=True)
class DesignMaskCandidate:
    canonical_family_position: int
    am1_mature_position: int | None
    am1_reference_residue: str
    observed_residue_identities: str
    template_coverage_count: int
    first_shell_observation_count: int
    second_sphere_observation_count: int
    hans_interface_observation_count: int
    interface_neighborhood: bool
    fixed_first_shell: bool
    mutable_second_sphere: bool
    mutable_interface: bool
    protected_positions: bool
    protection_reasons: str
    rationale: str


@dataclass(frozen=True, slots=True)
class DesignBackboneManifestRow:
    backbone_id: str
    structure_id: str
    source_kind: str
    source_path: str
    experimental_method: str
    selected_chains: str
    design_chains: str
    fixed_context_chains: str
    selected_residue_count: int
    design_chain_residue_count: int
    design_chain_metal_site_count: int
    representative_note: str
    backbone_pdb_path: str


@dataclass(frozen=True, slots=True)
class DesignCampaignManifestRow:
    campaign_id: str
    backbone_id: str
    design_set_name: str
    designed_chains: str
    fixed_context_chains: str
    exported_pdb_path: str
    chain_assignment_path: str
    fixed_positions_path: str
    designable_residue_count: int
    fixed_residue_count: int
    designable_canonical_positions: str
    designable_am1_positions: str


@dataclass(frozen=True, slots=True)
class DesignCampaignPositionRow:
    campaign_id: str
    design_set_name: str
    backbone_id: str
    structure_id: str
    chain_id: str
    sequence_index: int
    residue_seq: int
    insertion_code: str
    residue_name: str
    canonical_family_position: int | None
    am1_mature_position: int | None
    fixed_first_shell: bool
    protected_positions: bool
    mutable_second_sphere: bool
    mutable_interface: bool
    hard_fixed: bool
    designable: bool
    position_state: str
    position_reason: str
    rationale: str
