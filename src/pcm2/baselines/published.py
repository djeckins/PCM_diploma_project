"""Published criterion of Rao et al. (2019) on the vendored tables.

The method has no closed-form expression: it is a tabulated 100×100 energy
surface (hydrophobicity × radius, nm) plus the 2.6 kJ/mol (~1 RT) contour that
the authors extract when plotting. The tables are vendored byte-for-byte
(external/rao2019_heuristic, checksums in PROVENANCE.md; a test recomputes them).

The method predicts LOCAL DEWETTING, not conduction: conduction is a
consequence of dewetting, and not the only one.

The method outputs an energy and a binary verdict; it has no probabilistic
scale, so calibration measures are not computed for this arm.

Declared deviations from the paper (see PROVENANCE.md): per-frame application
(a window-averaged version is computed alongside and the difference is
reported) and the input radius comes from our own inscribed-sphere profile
estimator; inputs outside the surface domain are clipped to the domain
boundary, and the fraction of such frames is a reported quantity.
"""

from __future__ import annotations

import json
from functools import lru_cache

import numpy as np

from ..runtime import VENDORED_DIR

# Decision contour from the paper: barrier ~1 RT = 2.6 kJ/mol (Rao 2019, PNAS,
# Fig. 4).
RAO_DEWETTING_THRESHOLD_KJMOL = 2.6


@lru_cache(maxsize=1)
def _grid() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The vendored surface as (hydrophobicity axis, radius axis in nm, E in kJ/mol).

    Both axes come back sorted ascending and E is [len(h), len(r)], so E[i, j] is
    the energy at hydrophobicity h[i] and radius r[j]. Cached because the file is
    read once per process and the surface is fixed input, not state.
    """
    rows = json.loads((VENDORED_DIR / "heuristic_grid.json").read_text())
    h = np.array(sorted({r["hydrophobicity"] for r in rows}))
    r = np.array(sorted({r["radius"] for r in rows}))
    E = np.full((len(h), len(r)), np.nan)
    hi = {v: i for i, v in enumerate(h)}
    ri = {v: i for i, v in enumerate(r)}
    for row in rows:
        E[hi[row["hydrophobicity"]], ri[row["radius"]]] = row["energy"]
    # The file is a flat list of nodes, so a full rectangular grid is an
    # assumption about it, not a guarantee. A single unfilled cell means the
    # vendored file is not the published 100x100 surface, and interpolating
    # across the hole would return a number with no source.
    if np.any(np.isnan(E)):
        raise RuntimeError("vendored surface is incomplete — file substituted?")
    return h, r, E


def rao_energy(hydrophobicity: np.ndarray, radius_nm: np.ndarray) -> np.ndarray:
    """Bilinear interpolation of the authors' surface; out-of-domain inputs clip to its edge.

    Returns the water free-energy barrier in kJ/mol per input pair, NaN where
    either input is not finite. The radius is in NANOMETRES, the unit of the
    published surface; the feature layer divides its own A values by ten before
    calling here. Clipping is the declared deviation and is measured separately by
    out_of_domain_frac, so a system sitting mostly outside the tabulated range
    cannot pass unnoticed as an ordinary result.

    Interpolated, not snapped to the nearest node: this value is used as a
    continuous per-frame score. The authors' own nearest-node lookup is
    rao_nearest_energy, and that is what the structure score below uses.
    """
    h_ax, r_ax, E = _grid()
    h = np.clip(np.asarray(hydrophobicity, float), h_ax[0], h_ax[-1])
    r = np.clip(np.asarray(radius_nm, float), r_ax[0], r_ax[-1])
    out = np.full(h.shape, np.nan)
    m = np.isfinite(h) & np.isfinite(r)
    if not np.any(m):
        return out
    # Lower-left node of the enclosing cell. Both clips keep the cell (ih, ih+1)
    # inside the array for a point sitting exactly on a boundary node, where
    # searchsorted would otherwise give -1 or one index past the end.
    ih = np.clip(np.searchsorted(h_ax, h[m]) - 1, 0, len(h_ax) - 2)
    ir = np.clip(np.searchsorted(r_ax, r[m]) - 1, 0, len(r_ax) - 2)
    # Fractional position inside the cell along each axis, in [0, 1].
    th = (h[m] - h_ax[ih]) / (h_ax[ih + 1] - h_ax[ih])
    tr = (r[m] - r_ax[ir]) / (r_ax[ir + 1] - r_ax[ir])
    out[m] = (E[ih, ir] * (1 - th) * (1 - tr) + E[ih + 1, ir] * th * (1 - tr)
              + E[ih, ir + 1] * (1 - th) * tr + E[ih + 1, ir + 1] * th * tr)
    return out


# The paper's structure-level call: sigma_d > 0.55 marks a non-conductive state.
RAO_SIGMA_D_THRESHOLD = 0.55


def rao_nearest_energy(hydrophobicity: np.ndarray, radius_nm: np.ndarray) -> np.ndarray:
    """The authors' own lookup: nearest grid node, no interpolation.

    Their scripts minimize (grid_h - h)^2 + (grid_r - r)^2 over all nodes; on
    a regular product grid that separates into the nearest node per axis,
    reproducing the rule exactly (including its mixed-unit metric).

    Returns kJ/mol per input pair, NaN where either input is not finite. There is
    no clipping here and none is needed: the nearest node of a point outside the
    domain is a boundary node, which is what clipping would have produced.
    """
    h_ax, r_ax, E = _grid()
    h = np.asarray(hydrophobicity, float)
    r = np.asarray(radius_nm, float)
    out = np.full(h.shape, np.nan)
    m = np.isfinite(h) & np.isfinite(r)
    if not np.any(m):
        return out
    ih = np.argmin(np.abs(h[m][:, None] - h_ax[None, :]), axis=1)
    ir = np.argmin(np.abs(r[m][:, None] - r_ax[None, :]), axis=1)
    out[m] = E[ih, ir]
    return out


@lru_cache(maxsize=1)
def _contour_segments() -> np.ndarray:
    """The 1RT contour of the vendored surface as line segments in (h, r) space.

    Returns [n_segments, 4] rows (h_a, r_a, h_b, r_b). The contour is the decision
    boundary of the published criterion and has no closed form, so it is extracted
    from the tabulated surface the same way the authors extract it for plotting.
    """
    # Agg: no display is needed and none may be required — the contour is being
    # traced, not drawn. The backend is set before pyplot is imported, which is
    # why the imports are local to this function.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    h_ax, r_ax, E = _grid()
    fig = plt.figure()
    try:
        # E is transposed because contour(X, Y, Z) wants Z indexed [len(Y), len(X)],
        # whereas _grid returns E indexed [len(h), len(r)]. Getting this wrong
        # would trace the contour of the transposed surface, which on a
        # non-symmetric grid is a different curve.
        cs = plt.contour(h_ax, r_ax, E.T, levels=[RAO_DEWETTING_THRESHOLD_KJMOL])
        segs = []
        for poly in cs.allsegs[0]:
            for a, b in zip(poly[:-1], poly[1:]):
                segs.append((a[0], a[1], b[0], b[1]))
    finally:
        plt.close(fig)
    if not segs:
        raise RuntimeError("the 1RT contour is empty — vendored surface substituted?")
    return np.asarray(segs)


def _dist_to_contour(h: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Shortest distance to the contour, in the authors' mixed-unit metric.

    One distance per input point. The metric mixes dimensionless hydrophobicity
    with a radius in nm, so the value is not a length; it is reproduced as
    published, because the 0.55 threshold is stated in exactly these units and
    rescaling either axis would move the threshold with it.
    """
    segs = _contour_segments()
    # Every point against every segment: p is [n, 1, 2], a and b are [1, m, 2],
    # and the arithmetic below broadcasts to [n, m].
    p = np.stack([h, r], axis=1)[:, None, :]
    a = segs[None, :, 0:2]
    b = segs[None, :, 2:4]
    ab = b - a
    denom = np.sum(ab * ab, axis=2)
    # Point-to-segment, not point-to-line: t is the projection parameter along the
    # segment, clipped to [0, 1] so a point beyond an end measures to the endpoint.
    # A degenerate segment (both ends equal) would divide by zero; the denominator
    # is replaced by 1 there, and t is then 0, which is the segment's own point.
    t = np.clip(np.sum((p - a) * ab, axis=2) / np.where(denom == 0, 1.0, denom), 0.0, 1.0)
    proj = a + t[:, :, None] * ab
    d = np.linalg.norm(p - proj, axis=2)
    return d.min(axis=1)


def rao_sigma_d(hydrophobicity: np.ndarray, radius_nm: np.ndarray) -> tuple[float, int]:
    """The published structure score: Sigma_d over dewetted-side residue points.

    Each pore-lining residue contributes one (hydrophobicity, radius) point;
    points whose nearest-node energy exceeds the 1RT contour are flagged, and
    the score is the sum of their shortest distances to the contour. The paper
    calls a structure non-conductive when the score exceeds 0.55.

    Returns (score, number of flagged points). The inputs are per-residue arrays
    for ONE frame or structure, hydrophobicity on the Wimley-White axis and the
    local radius in nm; the score is a sum over points, so it grows with how many
    residues sit on the dewetted side and by how far.
    """
    h = np.asarray(hydrophobicity, float)
    r = np.asarray(radius_nm, float)
    m = np.isfinite(h) & np.isfinite(r)
    if not np.any(m):
        # Nothing measured on this frame: NaN, distinct from a measured zero.
        return float("nan"), 0
    e = rao_nearest_energy(h[m], r[m])
    flagged = e > RAO_DEWETTING_THRESHOLD_KJMOL
    if not np.any(flagged):
        # Measured, and no point past the contour: the score is genuinely 0.
        return 0.0, 0
    d = _dist_to_contour(h[m][flagged], r[m][flagged])
    return float(d.sum()), int(flagged.sum())


def out_of_domain_frac(hydrophobicity: np.ndarray, radius_nm: np.ndarray) -> float:
    """Fraction of the MEASURED input pairs that fall outside the tabulated surface.

    NaN when no pair is measured at all, so "nothing to check" is not reported as
    "nothing outside". Non-finite pairs are excluded from the denominator: they
    are a measurement gap, not an out-of-domain input. This is the accounting
    behind the clipping in rao_energy — the criterion is applied to systems the
    authors did not tabulate, and the size of that extrapolation is reported
    rather than hidden inside the clip.
    """
    h_ax, r_ax, _ = _grid()
    h = np.asarray(hydrophobicity, float)
    r = np.asarray(radius_nm, float)
    m = np.isfinite(h) & np.isfinite(r)
    if not np.any(m):
        return float("nan")
    outside = ((h[m] < h_ax[0]) | (h[m] > h_ax[-1])
               | (r[m] < r_ax[0]) | (r[m] > r_ax[-1]))
    return float(np.mean(outside))
