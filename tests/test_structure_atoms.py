from pathlib import Path

from lanm.structure.atoms import parse_atom_row, read_atom_records


def test_parse_atom_row_normalizes_types() -> None:
    row = {
        "record_type": "HETATM",
        "atom_serial": "42",
        "atom_name": " ND ",
        "alt_loc": "",
        "residue_name": "nd",
        "chain_id": "A",
        "residue_seq": "201",
        "insertion_code": "",
        "x": "1.0",
        "y": "2.5",
        "z": "3.5",
        "occupancy": "0.80",
        "b_factor": "",
        "element": "nd",
        "charge": "3+",
    }

    atom = parse_atom_row(row, "8FNS")

    assert atom.structure_id == "8FNS"
    assert atom.atom_serial == 42
    assert atom.residue_name == "ND"
    assert atom.element == "ND"
    assert atom.occupancy == 0.8
    assert atom.b_factor is None


def test_read_atom_records_derives_structure_id_from_filename(tmp_path: Path) -> None:
    path = tmp_path / "8dq2_atoms.csv"
    path.write_text(
        "record_type,atom_serial,atom_name,alt_loc,residue_name,chain_id,residue_seq,"
        "insertion_code,x,y,z,occupancy,b_factor,element,charge\n"
        "ATOM,1,N,,ALA,A,1,,0.0,0.0,0.0,1.00,10.0,N,\n",
        encoding="utf-8",
    )

    atoms = read_atom_records(path)

    assert len(atoms) == 1
    assert atoms[0].structure_id == "8DQ2"
    assert atoms[0].residue_key == ("A", "ALA", 1, "")

