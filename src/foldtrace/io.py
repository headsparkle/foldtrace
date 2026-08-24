"""Structure and site-file I/O for foldtrace."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
from Bio.Data.IUPACData import protein_letters_3to1 as _THREE_TO_ONE


@dataclass
class Structure:
    """CA-level view of one chain: coordinates, one-letter sequence, residue numbers, names."""
    ca: np.ndarray          # (N, 3) float
    seq: str                # length N, one-letter (X for non-standard)
    resnums: list[int]      # length N, author residue numbers
    resnames: list[str]     # length N, three-letter names


@dataclass
class Site:
    """One chemistry-defining position in the reference, plus the residues that count as retained."""
    label: str
    resnum: int
    expected: frozenset[str]   # one-letter codes that mean "retained"
    role: str
    key: bool                  # does this residue drive the overall verdict?


def _one(resname: str) -> str:
    return _THREE_TO_ONE.get(resname.capitalize(), "X")


def load_structure(path: str, chain: str | None = None) -> Structure:
    """Load the first model of a PDB/mmCIF file as a single-chain CA Structure.

    If ``chain`` is given, only that chain is read; otherwise the first chain is used.
    Only standard-polymer residues with a CA atom are kept.
    """
    parser = MMCIFParser(QUIET=True) if path.lower().endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
    model = next(iter(parser.get_structure("s", path)))
    ca, seq, resnums, resnames = [], [], [], []
    picked = None
    for ch in model:
        if chain is not None and ch.id != chain:
            continue
        if chain is None and picked is None:
            picked = ch.id
        if chain is None and ch.id != picked:
            continue
        for res in ch:
            if res.id[0] != " " or "CA" not in res:
                continue
            ca.append(res["CA"].coord)
            resnums.append(res.id[1])
            resnames.append(res.resname)
            seq.append(_one(res.resname))
    if not ca:
        raise ValueError(f"No CA atoms found in {path}" + (f" chain {chain}" if chain else ""))
    return Structure(np.asarray(ca, dtype=float), "".join(seq), resnums, resnames)


def parse_sites(path: str) -> list[Site]:
    """Parse a catalytic-sites TSV.

    Columns (tab-separated, header required): label, ref_resnum, expected, role, key
      - expected: one or more one-letter codes, comma-separated (e.g. ``E`` or ``Y,F``).
                  Any listed residue counts as "retained"; multiple codes express a
                  conservative/altered-but-functional set.
      - role:     free text (e.g. ``catalytic``, ``pocket``).
      - key:      1/true if this residue determines the overall verdict, else 0/false.
    Lines beginning with ``#`` are ignored.
    """
    sites: list[Site] = []
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = [p.strip().lower() for p in parts]
                required = {"label", "ref_resnum", "expected", "role", "key"}
                missing = required - set(header)
                if missing:
                    raise ValueError(f"sites file missing columns: {sorted(missing)}")
                continue
            row = dict(zip(header, parts))
            key = row["key"].strip().lower() in ("1", "true", "yes", "y")
            expected = frozenset(c.strip().upper() for c in row["expected"].split(",") if c.strip())
            sites.append(Site(row["label"].strip(), int(row["ref_resnum"]),
                              expected, row["role"].strip(), key))
    if not sites:
        raise ValueError(f"no sites parsed from {path}")
    return sites
