"""Shared types of the feature layer: column spec, contexts, constants with sources.

Rules hard-wired into the types:
  * a missing value means "nothing to measure", NEVER zero;
  * a conditional column has a companion indicator, a measured descriptor recording
    whether a subject was present; the column-to-indicator map lives in the schema;
  * all windows look strictly backward in time;
  * sign_flip marks DIRECTION-bearing quantities: their canonization is the only
    legitimate sign operation, and it is applied exactly once before writing.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np

# Number of named K-filter sites in the union schema: S0-S4.
# Source: Zhou Y., Morais-Cabral J.H., Kaufman A., MacKinnon R. (2001)
# "Chemistry of ion coordination and hydration revealed by a K+ channel–Fab
# complex at 2.0 A resolution", Nature 414:43–48, doi:10.1038/35102009.
NS_UNION_SITES = 5

from ..runtime import VENDORED_DIR  # noqa: E402 — single source of truth for the path


def load_ww_scale() -> dict[str, float]:
    """Wimley–White hydrophobicity scale, read from the vendored, checksummed
    file distributed by the Rao-2019 authors.

    Keyed by residue name as written in the topology; the value is the file's
    rescaled field, in [-1, 1] with +1 the most hydrophobic residue. That is the
    axis the Rao-2019 energy surface is tabulated on, so the whole-residue
    transfer free energies the same file carries under `original_hydrophobicity`
    are deliberately not read: substituting them would move the surface lookup
    onto a different scale without any error being raised.

    Aliases carry the same value for the same side chain: the scale is per-residue
    and achiral, and a D-amino acid bears the same chemistry (DVAL/DLEU: gramicidin;
    FVA: formyl-valine). HSD/HSE/HSP are the CHARMM protonation states of
    histidine, which the scale does not distinguish."""
    doc = json.loads((VENDORED_DIR / "wimley_white_1996.json").read_text())
    scale = {e["resname"]: float(e["hydrophobicity"]) for e in doc["hydrophobicity"]}
    for alias, canon in (("DVAL", "VAL"), ("DLEU", "LEU"), ("DALA", "ALA"),
                         ("DTRP", "TRP"), ("FVA", "VAL"), ("HSD", "HIS"),
                         ("HSE", "HIS"), ("HSP", "HIS")):
        if canon in scale:
            scale.setdefault(alias, scale[canon])
    return scale


@dataclass
class ColSpec:
    """One column of the feature table together with everything needed to read it.

    `estimator` is the prose definition of what the column measures and how, written
    once here so that the definition travels with the value instead of living in a
    document beside it; `unit` and `missing_means` are the strings the generated
    SCHEMA.md table shows, and a NaN means what `missing_means` says, never "zero".
    `indicator` names the companion column that records whether the subject of the
    measurement existed at all, and `indicator_min` the value of that companion at or
    above which a NaN here is a build failure rather than an honest gap. Instances are
    round-tripped through features/schema.json by as_dict and ColSpec(**d), so every
    field name is part of the artifact format.
    """

    name: str
    block: str
    unit: str
    estimator: str
    missing_means: str
    indicator: str | None = None
    indicator_min: float = 1.0  # threshold in the column -> indicator -> threshold map
    sign_flip: bool = False
    conditional: bool = False   # excluded from the missing-value budget by name
    to_model: bool = True       # False: stored in the table, not fed to the model

    def as_dict(self) -> dict:
        """Field-for-field dict for features/schema.json; ColSpec(**d) reads it back."""
        return dict(name=self.name, block=self.block, unit=self.unit,
                    estimator=self.estimator, missing_means=self.missing_means,
                    indicator=self.indicator, indicator_min=self.indicator_min,
                    sign_flip=self.sign_flip,
                    conditional=self.conditional, to_model=self.to_model)


@dataclass
class SchemaCtx:
    """Everything the LIST of columns depends on: blocks x sites x bins x flags.

    The sizes come from the architecture profile of the config, so the column list
    is built without opening a trajectory. partner_cutoff_A is in angstroms and
    window_ps in picoseconds; the partner cutoff is the one numeric parameter that
    does reach the column list, because it is rounded into the electrostatics
    column names and a table computed with another cutoff must not be mistaken
    for this one.

    params — numeric runtime parameters of the estimators (bin edges, cutoffs,
    event times); the LIST of columns does not depend on them.
    """
    blocks: list[str]
    n_bins: int
    n_sites: int | None
    sites_basis: str | None
    filter_present: bool
    n_subunits: int | None
    partner_cutoff_A: float
    window_ps: float
    params: dict = field(default_factory=dict)


@dataclass
class FrameData:
    """Per-frame shared measurements; each computed by exactly one estimator.

    One coordinate convention holds throughout. Axial values are offsets along the
    pore axis from the anchor of that frame (the centre of mass of the anchor
    selection), in angstroms and minimum-image, so a rigid shift of the frame or a
    change of periodic image leaves them unchanged; radial values are distances
    from the axis in the two perpendicular directions; times are in picoseconds.
    The axis keeps its laboratory direction here — only columns marked sign_flip
    are turned toward the conduction direction, once, before the table is written.
    Blocks read this object and never re-measure from coordinates, which is what
    keeps one quantity to one estimator.
    """
    time_ps: float = 0.0
    anchor: np.ndarray | None = None
    # pore profile
    prof_R: np.ndarray | None = None
    prof_search: np.ndarray | None = None
    prof_boundary: np.ndarray | None = None
    prof_z_offsets: np.ndarray | None = None
    # permeant
    ion_z_rel: np.ndarray | None = None       # axial coordinates from the anchor, all ions
    ion_rxy: np.ndarray | None = None
    ion_in_pore: np.ndarray | None = None     # mask
    innermost: int | None = None              # index of the ion closest to the constriction
    innermost_coord_n: float = np.nan
    # water
    wat_z_rel: np.ndarray | None = None       # only oxygens inside the pore cylinder
    wat_dipole_cos: np.ndarray | None = None
    # lining
    lining_resnames: list[str] = field(default_factory=list)
    lining_carbonyl_cos: float = np.nan
    lining_n_carbonyls: int = 0
    # sites
    site_centers_z: np.ndarray | None = None  # from the anchor; None if there are no sites
    # subunits
    subunit_com: np.ndarray | None = None     # [n_sub, 3] relative to the anchor
    # electrostatics (around the innermost ion)
    pot_by_partner: dict | None = None
    axial_charge_asym_e: float = np.nan
    # entry region
    mouth_counts: dict | None = None          # {"upper": n, "lower": n, ...}
    entry_dist_A: float = np.nan
    entry_wide_n: int = 0
    # Rao heuristic: per-residue points along the whole pore lining
    rao_res_resnames: list[str] = field(default_factory=list)
    rao_res_s_A: np.ndarray | None = None     # residue COG axial position from the anchor
    rao_res_R_A: np.ndarray | None = None     # local pore radius at that position


def bin_edges(bin_low_offset_A: float, bin_width_A: float, n_bins: int) -> np.ndarray:
    """Bin edges are fixed physical offsets from the anchor.
    The standard deviation of the bin width across frames is zero by construction.

    Returns n_bins+1 ascending edges in angstroms. edges[0] and edges[-1] are also
    the ends of the pore range that every block measures inside, inclusive at both
    ends; anything outside is not a gap but a different region of the system."""
    return bin_low_offset_A + bin_width_A * np.arange(n_bins + 1)


def backward_window_slice(times_ps: np.ndarray, i: int, window_ps: float) -> slice:
    """Window strictly backward in time within the SAME trajectory.

    Returns the slice of frames whose times lie in [times_ps[i] - window_ps,
    times_ps[i]], frame i included. times_ps must be ascending and must belong to a
    single replica: a slice taken across concatenated trajectories would let one
    replica's window read another replica's frames. The window is defined in time
    rather than in frames, so changing data.stride does not change what a windowed
    column measures."""
    t0 = times_ps[i] - window_ps
    j = int(np.searchsorted(times_ps, t0, side="left"))
    return slice(j, i + 1)
