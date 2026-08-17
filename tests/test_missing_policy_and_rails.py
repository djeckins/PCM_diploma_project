"""Defect classes: indicators wiped by the rail; matrix width depends on the
fold; NaN·0=NaN; infinity slipped through as a missing value."""

import numpy as np
import pytest

from pcm2.models.nanpolicy import MissingEncoder, passthrough_nan_check
from pcm2.relevance import ActivityRail


def test_encoder_values_plus_indicators():
    X = np.array([[1.0, np.nan], [3.0, 5.0]])
    enc = MissingEncoder().fit(X)
    Z = enc.transform(X)
    assert Z.shape == (2, 4)
    assert Z[0, 1] == 5.0  # median of the training rows
    assert Z[:, 2].tolist() == [0.0, 0.0] and Z[:, 3].tolist() == [1.0, 0.0]


def test_all_missing_column_stays_as_constant_plus_flag():
    X = np.array([[1.0, np.nan], [2.0, np.nan]])
    Z = MissingEncoder().fit(X).transform(X)
    assert Z.shape[1] == 4  # matrix width does not depend on the fold
    assert np.all(Z[:, 3] == 1.0)


def test_infinity_is_not_a_missing_value():
    X = np.array([[np.inf, 1.0]])
    with pytest.raises(ValueError, match="nfinity"):
        MissingEncoder().fit(X)
    with pytest.raises(ValueError, match="nfinity"):
        passthrough_nan_check(X)


def test_rail_zeroes_by_assignment_not_multiplication():
    # NaN·0=NaN defect: a column with missing values zeroed by multiplication
    # becomes a new missingness indicator under a foreign name.
    X = np.array([[1.0, np.nan], [1.0, 2.0]])
    rail = ActivityRail(n_value_cols=2, allow_nan=True).fit(X)
    Z = rail.transform(X)
    assert not np.any(np.isnan(Z[:, 0])), "column must be zeroed by assignment"
    assert np.all(Z[:, 0] == 0.0)


def test_rail_exempts_indicator_block():
    # Matrix [value, indicator]: the indicator is constant (the whole fold has no
    # missing values), yet the rail must not touch it — exempt by construction.
    X = np.array([[1.0, 0.0], [2.0, 0.0]])
    rail = ActivityRail(n_value_cols=1, allow_nan=False).fit(X)
    Z = rail.transform(X)
    assert np.array_equal(Z[:, 1], X[:, 1])
    # Test self-check: a rail applied to the full matrix would wipe the indicator.
    rail_bad = ActivityRail(n_value_cols=2, allow_nan=False).fit(X)
    assert rail_bad.transform(X)[:, 1].tolist() == [0.0, 0.0]


def test_rail_keeps_constant_with_nan_pattern_for_native_branch():
    # A constant-with-missing column carries the missingness pattern — the very
    # mechanism of the native branch; it must not be zeroed.
    X = np.array([[3.0, 1.0], [np.nan, 2.0], [3.0, 3.0]])
    rail = ActivityRail(n_value_cols=None, allow_nan=True).fit(X)
    Z = rail.transform(X)
    assert np.isnan(Z[1, 0]) and Z[0, 0] == 3.0
