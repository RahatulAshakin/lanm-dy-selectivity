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
