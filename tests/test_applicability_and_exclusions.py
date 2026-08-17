"""Defect classes covered: an inapplicable value gets imputed; a constant is
confused with an inapplicable one; the baseline score leaks into the feature
matrix; an affine identity is not re-verified on the data."""

import numpy as np
import pandas as pd

from pcm2.datasets import find_affine_duplicates
from pcm2.features._common import ColSpec


def _spec(name, block, **kw):
    return ColSpec(name, block, "u", "e", "m", **kw)


def test_baseline_score_never_reaches_model():
    from pcm2.features import BASELINE_COLS
    for c in BASELINE_COLS:
        assert not c.to_model, "the baseline score is a competing predictor"


def test_affine_duplicate_found_and_reverified_on_data():
    x = np.linspace(0, 10, 50)
    df = pd.DataFrame({"a": x, "b": 2.0 * x - 3.0, "c": np.sin(x)})
    pairs = find_affine_duplicates(df)
    assert ("a", "b") in pairs
    assert all(p[1] != "c" for p in pairs)
    # The identity no longer holds on the data → the pair must disappear (re-verification).
    df2 = df.copy()
    df2.loc[10, "b"] += 0.5
    assert find_affine_duplicates(df2) == []


def test_structural_and_constant_are_distinct_verdicts():
    """Post-hoc statistics cannot distinguish 'inapplicable' from 'constant' — the
    structure can; check that the verdicts stay separate and the counts add up."""
    from pcm2.features import resolve_applicability

    class FakeCfg(dict):
        def __getitem__(self, k):
            return {"system.arch_profile.n_sites": 0,
                    "system.arch_profile.filter_present": False,
                    "system.arch_profile.n_subunits": 2}[k]

    schema = [_spec("ns_S0_occ", "named_sites"),
              _spec("sym_nn_dist_cv", "symmetry"),
              _spec("sym_com_z_spread_A", "symmetry"),
              _spec("hyd_water_pore_n", "hydration"),
              _spec("geo_search_frac", "geometry")]
    table = pd.DataFrame({
        "ns_S0_occ": [np.nan, np.nan],
        "sym_nn_dist_cv": [np.nan, np.nan],
        "sym_com_z_spread_A": [1.0, 2.0],
        "hyd_water_pore_n": [5.0, 5.0],   # measured and constant
        "geo_search_frac": [0.7, 0.9],
    })
    v = resolve_applicability(FakeCfg(), table, schema)
    assert v["ns_S0_occ"] == "inapplicable_structural"      # no filter: a fact about the channel
    assert v["sym_nn_dist_cv"] == "inapplicable_structural"  # dimer: no circular neighbors
    assert v["sym_com_z_spread_A"] == "active"
    assert v["hyd_water_pore_n"] == "constant"               # out of the model, yet a measurement
    assert v["geo_search_frac"] == "active"
    counts = {s: sum(1 for x in v.values() if x == s)
              for s in ("active", "inapplicable_structural", "constant")}
    assert sum(counts.values()) == len(schema), "verdict counts must add up to the schema size"
