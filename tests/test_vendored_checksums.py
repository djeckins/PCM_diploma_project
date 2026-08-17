"""Defect class: silent substitution of the method authors' vendored files."""

import hashlib
import json

from pcm2.runtime import VENDORED_DIR

EXPECTED = {
    "heuristic_grid.json":
        "bb29c21e1339e71e7a877d5095a4d72c4ac3ab06ce04c863fbcfb407e567c4f0",
    "wimley_white_1996.json":
        "d898884281bf666074ea8c2a5629896b31b9bea6bed4d4389a736ef355542b5c",
}


def test_checksums_match_provenance():
    for name, want in EXPECTED.items():
        got = hashlib.sha256((VENDORED_DIR / name).read_bytes()).hexdigest()
        assert got == want, f"{name}: file substituted ({got})"
    prov = (VENDORED_DIR / "PROVENANCE.md").read_text()
    for want in EXPECTED.values():
        assert want in prov, "PROVENANCE.md is missing a checksum"


def test_grid_structural_invariants():
    rows = json.loads((VENDORED_DIR / "heuristic_grid.json").read_text())
    hs = {r["hydrophobicity"] for r in rows}
    rs = {r["radius"] for r in rows}
    assert len(rows) == len(hs) * len(rs), "grid is not regular"
    assert min(rs) == 0.1 and max(rs) == 0.6, "radii in nm, the paper's domain"
    assert min(hs) == -0.45 and abs(max(hs) - 0.30) < 1e-9


def test_ww_scale_matches_paper_normalisation():
    doc = json.loads((VENDORED_DIR / "wimley_white_1996.json").read_text())
    vals = {e["resname"]: e["hydrophobicity"] for e in doc["hydrophobicity"]}
    assert len(vals) == 20
    assert vals["ASP"] == -1.0  # scale minimum, as in PNAS Fig. 4
    assert abs(max(vals.values()) - 0.336) < 0.01
