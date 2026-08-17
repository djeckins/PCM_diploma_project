"""Defect classes: a literal column count; a column without a block; an indicator
absent from the schema; two implementations of one quantity."""

import numpy as np

from pcm2.features import BLOCK_MODULES
from pcm2.features._common import NS_UNION_SITES, SchemaCtx


def _ctx(**kw):
    base = dict(blocks=list(BLOCK_MODULES), n_bins=8, n_sites=5,
                sites_basis="motif_rings", filter_present=True, n_subunits=4,
                partner_cutoff_A=12.0, window_ps=1000.0)
    base.update(kw)
    return SchemaCtx(**base)


def _full_schema(ctx):
    cols = []
    for b in ctx.blocks:
        cols.extend(BLOCK_MODULES[b].schema(ctx))
    return cols


def test_schema_is_a_function_of_bins_and_sites():
    n8 = len(_full_schema(_ctx(n_bins=8)))
    n5 = len(_full_schema(_ctx(n_bins=5)))
    # Five column families per bin: R, nsearch, water count, and the
    # system-relative forms (water density, R relative to the constriction).
    assert n8 - n5 == 3 * 5
    # Density-peak sites add columns; motif-ring sites always live in named_sites.
    d2 = len(_full_schema(_ctx(sites_basis="density_peaks", filter_present=False,
                               n_sites=2)))
    d4 = len(_full_schema(_ctx(sites_basis="density_peaks", filter_present=False,
                               n_sites=4)))
    assert d4 - d2 == 2


def test_no_duplicate_names_across_blocks():
    cols = _full_schema(_ctx())
    names = [c.name for c in cols]
    assert len(names) == len(set(names)), "two implementations of one quantity"


def test_every_column_has_block_and_indicators_exist():
    cols = _full_schema(_ctx())
    names = {c.name for c in cols}
    for c in cols:
        assert c.block, f"{c.name}: a column without a block breaks the assembly"
        if c.indicator is not None:
            assert c.indicator in names, f"{c.name}: indicator {c.indicator} does not exist"


def test_named_sites_are_union_wide_even_without_filter():
    cols = _full_schema(_ctx(filter_present=False, sites_basis=None, n_sites=None))
    ns = [c for c in cols if c.name.startswith("ns_S")]
    assert len(ns) == NS_UNION_SITES, \
        "union schema: site columns are present and flagged as missing, not removed"


def test_no_literal_column_count_in_sources():
    from pathlib import Path
    src = Path(__file__).resolve().parents[1] / "src" / "pcm2"
    import re
    pat = re.compile(r"(n_columns|len\(schema\))\s*==\s*\d{2,}")
    offenders = [str(p) for p in src.rglob("*.py") if pat.search(p.read_text())]
    assert offenders == [], f"literal column count: {offenders}"


def test_refused_estimator_yields_nan_not_substitute():
    """Where an estimator refused, the cell is missing — never another estimator's result."""
    from pcm2.features import geometry
    from pcm2.features._common import FrameData
    ctx = _ctx()
    ctx.params["bin_edges"] = np.array([-2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0])
    fd = FrameData()
    fd.prof_z_offsets = np.arange(-2.0, 14.1, 1.0)
    fd.prof_R = np.full(len(fd.prof_z_offsets), 3.0)
    fd.prof_search = np.zeros(len(fd.prof_z_offsets), dtype=bool)  # every slice refused
    fd.prof_boundary = np.zeros(len(fd.prof_z_offsets), dtype=bool)
    cols = {c.name: np.full(1, np.nan) for c in geometry.schema(ctx)}
    geometry.compute(ctx, fd, 0, cols)
    assert np.isnan(cols["geo_r_constriction_A"][0])
    assert cols["geo_n_search"][0] == 0.0  # the indicator explains the missing value
    for k in range(ctx.n_bins):
        assert np.isnan(cols[f"geo_r_bin{k}_A"][0])
