"""Defect class: direction canonization applied to the wrong columns — a structural
quantity inverted together with a field-like one, or canonization failing to be the
identity for a positive direction."""

import numpy as np

from pcm2.features import BLOCK_MODULES, build_ctx  # noqa: F401
from pcm2.features._common import ColSpec


def _flip(cols: dict, schema: list[ColSpec], direction: int) -> dict:
    out = {k: v.copy() for k, v in cols.items()}
    if direction is not None and direction < 0:
        for c in schema:
            if c.sign_flip:
                out[c.name] = -out[c.name]
    return out


SCHEMA = [
    ColSpec("ww_dipole_cos_mean", "water_wire", "cos", "e", "m", sign_flip=True),
    ColSpec("geo_z_constriction_A", "geometry", "A", "e", "m", sign_flip=False),
    ColSpec("ns_S0_occ", "named_sites", "count", "e", "m", sign_flip=False),
]


def test_positive_direction_is_bitwise_identity():
    cols = {c.name: np.array([0.3, -0.2, np.nan]) for c in SCHEMA}
    flipped = _flip(cols, SCHEMA, +1)
    for k in cols:
        assert np.array_equal(cols[k], flipped[k], equal_nan=True)


def test_negative_direction_flips_only_sign_bearing():
    cols = {c.name: np.array([0.3, -0.2]) for c in SCHEMA}
    flipped = _flip(cols, SCHEMA, -1)
    assert np.array_equal(flipped["ww_dipole_cos_mean"], [-0.3, 0.2])
    assert np.array_equal(flipped["geo_z_constriction_A"], [0.3, -0.2])
    assert np.array_equal(flipped["ns_S0_occ"], [0.3, -0.2])


def test_ordered_blocks_and_positions_never_flip_in_real_schema():
    """In the real schema no axial bin, no site and no structural coordinate is
    marked sign_flip: the protein is not flipped."""
    from pcm2.features import geometry, hydration, named_sites, occupancy
    from pcm2.features._common import SchemaCtx
    ctx = SchemaCtx(blocks=[], n_bins=4, n_sites=3, sites_basis="density_peaks",
                    filter_present=True, n_subunits=4, partner_cutoff_A=12.0,
                    window_ps=1000.0)
    for mod in (geometry, hydration, named_sites, occupancy):
        for c in mod.schema(ctx):
            if any(t in c.name for t in ("_bin", "_S", "site", "z_constr", "innermost_z")):
                assert not c.sign_flip, \
                    f"{c.name} is flipped — sign_flip marks direction-bearing quantities only"
