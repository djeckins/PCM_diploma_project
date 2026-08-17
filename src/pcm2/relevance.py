"""Rails: near-duplicate reporting and activity screening. Trainable pipeline steps.

Both zero out BY ASSIGNMENT rather than by deletion or multiplication: deletion
changes the matrix width and column positions (breaking the indicator block),
and multiplying by zero yields NaN·0=NaN, a new missingness indicator under a
foreign name.

Near-duplicates are not zeroed. Authors who established the instability of
feature importances under collinearity recommend changing how importances are
aggregated (up to column blocks); the rail only reports the pairs.

The activity rail zeroes a column whose variance is zero on the fold's training
rows — a SINGLE criterion. A "few distinct values" criterion is forbidden: a
binary indicator has exactly two values, so such a rule would wipe out the
entire indicator block at once. The indicator block is exempt from the check.
"""

from __future__ import annotations

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


class NearDuplicateReporter(BaseEstimator, TransformerMixin):
    """Reports column pairs with |corr| above the threshold; the transform is
    identity, and the threshold is reported without acting on the matrix."""

    def __init__(self, threshold: float = 0.98, feature_names: list[str] | None = None):
        # threshold on |Pearson r|; the run takes it from evaluate.corr_threshold.
        # feature_names are the columns of the model matrix BEFORE the encoder
        # doubled the width, which is what _name needs to rebuild an indicator name.
        self.threshold = threshold
        self.feature_names = feature_names

    def _name(self, i: int):
        """The reporter sits AFTER the encoder: the right half of the matrix
        holds missingness indicators, whose names derive from the value names."""
        if not self.feature_names:
            return i
        n = len(self.feature_names)
        if i < n:
            return self.feature_names[i]
        if i < 2 * n:
            return f"missing[{self.feature_names[i - n]}]"
        return i

    def fit(self, X, y=None):
        """Collect the near-duplicate pairs of the fold's training rows into report_.

        Each entry is (name_a, name_b, r) with the signed correlation, so a pair
        joined by a minus sign is recognizable in the artifact. The training step
        copies report_ into the fold record under dup_pairs.
        """
        X = np.asarray(X, float)
        with np.errstate(invalid="ignore"):
            # A column that does not vary on these rows has no correlation with
            # anything: r would be 0/0. Such columns are left out of the correlation
            # matrix, and idx below maps what remains back to the original positions.
            std = np.nanstd(X, axis=0)
            ok = std > 0
            # In the fitted pipeline the reporter sits after the encoder, so no
            # missing values reach it; the fill is for a bare use of the transformer.
            # Zero is the mean of a standardized column, so a substituted row pulls
            # the correlation towards zero — the reporter then misses a pair rather
            # than inventing one.
            Xc = np.where(np.isnan(X), 0.0, X)
            # Fewer than two varying columns leaves no off-diagonal to look at, and
            # np.corrcoef of a single column returns a scalar rather than a matrix.
            C = np.corrcoef(Xc[:, ok], rowvar=False) if ok.sum() >= 2 else np.empty((0, 0))
        pairs = []
        idx = np.flatnonzero(ok)
        for a in range(C.shape[0]):
            for b in range(a + 1, C.shape[1]):
                if abs(C[a, b]) >= self.threshold:
                    pairs.append((self._name(int(idx[a])), self._name(int(idx[b])),
                                  float(C[a, b])))
        self.report_ = pairs
        return self

    def transform(self, X):
        """The matrix passes through untouched — this step only records a report.

        The step stays inside the pipeline so that the pairs are recomputed on each
        fold's own training rows; a report made once on the whole sample would be a
        statement about a matrix no fold ever saw.
        """
        return X


class ActivityRail(BaseEstimator, TransformerMixin):
    """Zeroes (by assignment) columns with zero variance on the fold's training rows.

    n_value_cols: width of the value block; everything to its right is the
    indicator block, which is exempt from the check. For the native-NaN branch
    (n_value_cols=None, allow_nan=True) only columns constant WITHOUT missing
    values, or missing entirely, are zeroed: a constant-with-missing column
    carries a missingness pattern, which is the information that branch uses.
    """

    def __init__(self, n_value_cols: int | None = None, allow_nan: bool = False):
        # allow_nan selects the criterion, not a tolerance: False is the encoded
        # branch, where missing values have already been imputed away, True the
        # native branch, where a NaN is still a NaN and means "nothing to measure".
        self.n_value_cols = n_value_cols
        self.allow_nan = allow_nan

    def fit(self, X, y=None):
        """Record which columns are inactive on these rows, in zero_mask_.

        Fitted on the training rows of one fold only. A column may be inactive in
        one fold and informative in another — that is a property of the fold, and
        deciding it on the whole sample would let the test rows choose the columns.
        """
        X = np.asarray(X, float)
        n = X.shape[1]
        # Everything from limit onwards is the indicator block and is never checked.
        limit = self.n_value_cols if self.n_value_cols is not None else n
        zero = np.zeros(n, dtype=bool)
        for j in range(limit):
            col = X[:, j]
            if self.allow_nan:
                obs = col[np.isfinite(col)]
                # Nothing measured anywhere in the training rows: there is no value
                # to learn from, and a tree would only be splitting on missingness
                # that never varies.
                all_nan = len(obs) == 0
                # Complete and constant: no variation of any kind. A column that is
                # constant WHERE observed but missing elsewhere is deliberately kept
                # — its missingness pattern varies, and in the native branch that
                # pattern is a descriptor in its own right.
                const_no_nan = len(obs) == len(col) and (len(obs) == 0 or np.all(obs == obs[0]))
                zero[j] = all_nan or const_no_nan
            else:
                # Encoded branch: the encoder has already imputed, so the column
                # holds no NaN and exact equality to the first row is the test for a
                # constant. Standardization has already collapsed such a column to
                # zero; the mask records the fact and keeps both branches under the
                # single zero-variance rule.
                zero[j] = np.all(col == col[0])
        self.zero_mask_ = zero
        return self

    def transform(self, X):
        """Return the matrix with the inactive columns set to zero, same width as the input.

        The copy matters: the same feature matrix is indexed again by the other arms
        and the other folds, and an in-place write would silently zero a column for
        all of them.
        """
        X = np.asarray(X, float).copy()
        X[:, self.zero_mask_] = 0.0
        return X
