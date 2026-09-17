"""Tests for structure/site I/O, in particular modified-residue (HETATM) handling."""
from foldtrace.io import load_structure


def _atom(serial, name, resname, resseq, x, y, z, het=False, elem=None):
    rec = "HETATM" if het else "ATOM  "
    nm = (" " + name).ljust(4) if len(name) < 4 else name
    elem = elem or name[0]
    return (f"{rec}{serial:>5} {nm} {resname:>3} A{resseq:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {elem:>2}")


def _write_pdb(path):
    """Two standard residues, a modified catalytic residue as HETATM (KCX, carbamyl-Lys,
    with a full N/CA/C backbone), plus a water and a nickel ion as bare HETATMs."""
    lines, serial = [], 1

    def residue(resname, resseq, base, het=False):
        nonlocal serial
        for name, dx in (("N", 0.0), ("CA", 1.5), ("C", 3.0)):
            lines.append(_atom(serial, name, resname, resseq, base + dx, 0.0, 0.0, het=het))
            serial += 1

    residue("GLY", 1, 0.0)
    residue("ALA", 2, 4.5)
    residue("KCX", 217, 9.0, het=True)                        # modified Lys, full backbone
    lines.append(_atom(serial, "O", "HOH", 500, 20, 20, 20, het=True, elem="O")); serial += 1
    lines.append(_atom(serial, "NI", "NI", 600, 25, 25, 25, het=True, elem="NI"))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def test_modified_residue_hetatm_is_kept_water_and_metal_dropped(tmp_path):
    p = tmp_path / "mod.pdb"
    _write_pdb(p)
    st = load_structure(str(p))
    # KCX (carbamylated Lys217) is a modified *polymer* residue: it must appear in the
    # CA trace so the reference has no spurious gap at the catalytic position.
    assert 217 in st.resnums
    assert "KCX" in st.resnames
    # the standard residues are still present, and water + nickel are not.
    assert st.resnums == [1, 2, 217]
    assert st.seq == "GAX"          # KCX has no standard one-letter code -> X
    assert "HOH" not in st.resnames and "NI" not in st.resnames
