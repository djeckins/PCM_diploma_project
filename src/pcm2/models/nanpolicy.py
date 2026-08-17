"""The SINGLE definition of what a missing value means to the model.

The two branches see different matrices by design:
  * boosting — native: a missing value stays NaN, the tree learns the split
    direction;
  * linear — encoded: [values ++ indicators], one indicator per column; the
    value block is filled with the median computed ONLY on the fold's training
    rows (a transformer inside the pipeline, so per-fold fitting is automatic).

The encoded branch shows the model both the imputed value and the fact of
imputation. A column missing in every row REMAINS (a constant next to an
all-ones flag), otherwise the matrix width would depend on the fold. Infinity
is not a missing value and is caught by a separate check.
"""

from __future__ import annotations

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


class MissingEncoder(BaseEstimator, TransformerMixin):
    """[values (per-fold median imputation) ++ missingness indicators].

    A transform of width 2p from an input of width p: the left half carries the
    values with the training median substituted where the measurement is absent,
    the right half carries 1.0 exactly where the input was NaN. Column order is
    preserved in both halves, so column j of the input has its indicator at
    j + p; linear.linear_contributions and NearDuplicateReporter._name both rely
    on that offset.

    The median comes from the rows this step is fitted on. Inside a pipeline
    fitted per fold, that is the training part of the fold and nothing else, so
    no test-row value reaches the imputation.
    """

    def fit(self, X, y=None):
        X = np.asarray(X, float)
        if np.any(np.isinf(X)):
            raise ValueError("infinity in the matrix: this is not a missing value")
        # Median, not mean: the columns include radii, counts and waiting times
        # whose distributions are skewed, and the median of such a column is a
        # value the column actually attains.
        med = np.nanmedian(X, axis=0)
        # A column missing in all training rows has no median: use constant 0
        # next to an all-ones flag, which carries the information.
        self.medians_ = np.where(np.isfinite(med), med, 0.0)
        self.n_features_in_ = X.shape[1]
        return self

    def transform(self, X):
        """[n, 2p] float matrix: imputed values, then the missingness flags.

        The infinity check is repeated here and not only in fit, because
        transform is also called on the test rows of the fold, which fit never
        sees.
        """
        X = np.asarray(X, float)
        if np.any(np.isinf(X)):
            raise ValueError("infinity in the matrix: this is not a missing value")
        ind = np.isnan(X).astype(float)
        vals = np.where(np.isnan(X), self.medians_[None, :], X)
        return np.hstack([vals, ind])

    def restore_column_contributions(self, contrib: np.ndarray) -> np.ndarray:
        """Contributions of each (value, indicator) pair are summed back to the source column.

        Takes contributions in the encoded space, [n, 2p], and returns [n, p] on
        the source columns. The sum is exact because the linear model is additive
        in the encoded columns: the value term and the indicator term of one
        descriptor are two parts of what that descriptor did to the logit, and
        reporting them apart would split one mechanism across two names.
        """
        n = self.n_features_in_
        return contrib[:, :n] + contrib[:, n:2 * n]


def passthrough_nan_check(X: np.ndarray) -> np.ndarray:
    """Boosting branch: a missing value stays NaN; infinity is a separate error.

    Returns the matrix as float, unchanged apart from the dtype. Infinity is
    refused rather than mapped to NaN: xgboost's `missing=np.nan` would treat
    an infinity as an ordinary large value and send it down one branch, so a
    broken estimator would be silently absorbed into the fit.
    """
    X = np.asarray(X, float)
    if np.any(np.isinf(X)):
        raise ValueError("infinity in the matrix: this is not a missing value")
    return X
