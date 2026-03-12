from lanm.models import AtomRecord
from lanm.structure.geometry import annotate_metal_site, collect_atoms_within_cutoff, euclidean_distance


def _atom(
    atom_serial: int,
    atom_name: str,
    residue_name: str,
    residue_seq: int,
    element: str,
    x: float,
    y: float,
    z: float,
    record_type: str = "ATOM",
    chain_id: str = "A",
    charge: str = "",
) -> AtomRecord:
    return AtomRecord(
        structure_id="TEST",
        record_type=record_type,
        atom_serial=atom_serial,
        atom_name=atom_name,
        alt_loc="",
        residue_name=residue_name,
        chain_id=chain_id,
        residue_seq=residue_seq,
        insertion_code="",
        x=x,
        y=y,
        z=z,
        occupancy=1.0,
        b_factor=10.0,
        element=element,
        charge=charge,
    )


def test_collect_atoms_within_cutoff_sorts_by_distance_then_identity() -> None:
    metal = _atom(1, "ND", "ND", 201, "ND", 0.0, 0.0, 0.0, record_type="HETATM", charge="3+")
    near_b = _atom(3, "O", "ASP", 31, "O", 2.0, 0.0, 0.0)
    near_a = _atom(2, "OD1", "ASP", 30, "O", 2.0, 0.0, 0.0)

    matches = collect_atoms_within_cutoff(metal, [metal, near_b, near_a], 3.2)

    assert [item.atom.atom_serial for item in matches] == [2, 3]
    assert euclidean_distance(metal, near_a) == 2.0


def test_annotate_metal_site_builds_first_and_second_shells() -> None:
    metal = _atom(1, "ND", "ND", 201, "ND", 0.0, 0.0, 0.0, record_type="HETATM", charge="3+")
    donor_1 = _atom(2, "OD1", "ASP", 30, "O", 2.1, 0.0, 0.0)
    donor_2 = _atom(3, "OD2", "ASP", 30, "O", 2.8, 0.0, 0.0)
    water = _atom(4, "O", "HOH", 401, "O", 3.0, 0.0, 0.0, record_type="HETATM")
    second_shell = _atom(5, "CB", "ASN", 55, "C", 5.3, 0.0, 0.0)
    outside = _atom(6, "CG", "TYR", 80, "C", 7.0, 0.0, 0.0)

    summary, annotations = annotate_metal_site(
        metal_atom=metal,
        atoms=[metal, donor_1, donor_2, water, second_shell, outside],
        first_shell_cutoff_A=3.2,
        second_sphere_cutoff_A=6.0,
    )

    assert summary.donor_atom_count == 3
    assert summary.donor_residue_count == 2
    assert summary.residue_within_6a_count == 3
    assert summary.nearest_donor_distance_A == 2.1
    assert summary.farthest_donor_distance_A == 3.0
    assert [row.shell_type for row in annotations].count("first_shell") == 3
    assert [row.shell_type for row in annotations].count("second_sphere") == 1
    assert any(row.is_water for row in annotations if row.shell_type == "first_shell")
