"""Typed config schema; an unknown key is an error.

One config per system. Every key carries a provenance -- detected / declared /
default -- in the `origins` section of the same file, so the three states stay
distinguishable in the artifact. A typo in a key name aborts the run rather
than falling back to the default.

detected: measured from the topology or the trajectory by the autodetect step.
declared: chosen by a person, with the reason for the choice written into the
`basis` string beside it. default: the schema value, which is never written to
the file at all -- save() prunes it, so a default that later turns out to be
wrong is corrected in one place and no config quietly disagrees with the
schema. The distinction is what makes a number in a published table answerable:
either the run measured it, or someone stated why it is what it is.

The schema below is also where the reasoning behind the numeric parameters
lives; the `doc` text of a Key sits next to the constraint it explains rather
than in prose somewhere else.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Callable

import yaml

# The closed set of provenance states. Nothing else is accepted in `origins`,
# so "unknown provenance" cannot be spelled.
ORIGINS = ("detected", "declared", "default")


class ConfigError(Exception):
    """A config disagrees with the schema. Always fatal: no step runs on a
    config it could not fully validate."""
    pass


class _Required:
    """Marker for a key that has no default and must be present in the file."""

    def __repr__(self) -> str:  # pragma: no cover
        return "REQUIRED"


REQUIRED = _Required()


class Key:
    """One leaf of the schema: type predicate, default, and the note explaining it.

    check is the predicate the value must satisfy; typename is the human name of
    the type, used only in the error message; default is the schema value, or
    REQUIRED when the key must be supplied; doc records why the key exists or
    the rule by which its value is chosen.
    """

    def __init__(self, check: Callable[[Any], bool], typename: str,
                 default: Any = REQUIRED, doc: str = ""):
        self.check = check
        self.typename = typename
        self.default = default
        self.doc = doc


class RecordList:
    """List of homogeneous records (conditions, replicas).

    Records describe the dataset itself -- which conditions exist and which
    trajectory files belong to them -- so they have no default, cannot be
    filled in, and are always written to the file in full.
    """

    def __init__(self, spec: dict):
        self.spec = spec


def _is_str(v):
    return isinstance(v, str)


def _is_bool(v):
    return isinstance(v, bool)


# bool is a subclass of int in Python, so a bare isinstance check would let
# `stride: true` through as 1 and `n_bins: false` through as 0. The exclusion
# keeps a typed key typed.
def _is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def _is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


# null is a state of its own, not a missing value: a key that autodetect has not
# filled in yet reads as None, and the step that needs it refuses by name rather
# than proceeding on a default.
def _opt(f):
    return lambda v: v is None or f(v)


def _enum(*vals):
    return lambda v: v in vals


def _list_of(f):
    return lambda v: isinstance(v, list) and all(f(x) for x in v)


# model.monotone: column name -> monotonicity direction. The map lists only the
# columns that carry a constraint -- the tree arm fills 0 for every column absent
# from it -- so an entry must be a real direction, +1 or -1. Admitting 0 here
# would allow a constraint that states nothing.
def _map_str_int(v):
    return isinstance(v, dict) and all(
        _is_str(k) and _is_int(x) and x in (-1, 1) for k, x in v.items())


_PROTOCOL = {
    # The simulation protocol is part of the system identity: a different
    # protocol means different physics.
    "forcefield": Key(_opt(_is_str), "str|null", None),
    "water_model": Key(_opt(_is_str), "str|null", None),
    "charge_scaling": Key(_opt(_is_num), "float|null", None,
                          "ECC scaling of ion charges; measured from the topology"),
    "voltage_mV": Key(_opt(_is_num), "float|null", None),
    "temperature_K": Key(_opt(_is_num), "float|null", None,
                         "thermostat reference temperature; without the MD engine it is declared"),
}

# One replica record per trajectory file. `lineage` is not decoration: splits are
# made over lineages, so two replicas that share an equilibrated ancestor must
# share a lineage string or the same starting structure would appear on both
# sides of a split and the reported score would be optimistic.
_REPLICA = {
    "id": Key(_is_str, "str"),
    "topology": Key(_is_str, "str"),
    "trajectory": Key(_is_str, "str"),
    "events": Key(_opt(_is_str), "str|null", None, "precomputed crossing annotations, if any"),
    "lineage": Key(_is_str, "str", doc="independent-build identifier; the unit of data splitting"),
    "protocol": _PROTOCOL,
}

_CONDITION = {
    "id": Key(_is_str, "str"),
    "direction": Key(_opt(lambda v: v in (-1, 1)), "±1|null", None,
                     "sign of conduction along the axis; all replicas of a condition must agree"),
    "replicas": RecordList(_REPLICA),
}

# The whole schema. A nested dict is a section, a Key is a leaf, a RecordList is
# dataset data. A section may be left out of the file entirely as long as none of
# its leaves is REQUIRED; what is left out is filled from the defaults here and
# is not written back on save.
SCHEMA: dict = {
    "system": {
        "id": Key(_is_str, "str"),
        "architecture": Key(_opt(_enum("k_filter_channel", "single_file_channel")),
                            "enum|null", None),
        "permeant": Key(_opt(_is_str), "resname|null", None,
                        "defined by per-particle crossings, not by abundance"),
        "arch_profile": {
            "n_subunits": Key(_opt(_is_int), "int|null", None),
            "filter_present": Key(_opt(_is_bool), "bool|null", None),
            "filter_motif": Key(_opt(_is_str), "str|null", None),
            "filter_motif_pattern": Key(_is_str, "regex", "T[VIL]G[YF]G",
                                        "selectivity-filter sequence signature; a declared "
                                        "override admits engineered filter mutants while the "
                                        "geometric ring confirmation still gates detection"),
            "n_sites": Key(_opt(_is_int), "int|null", None),
            "sites_basis": Key(_opt(_enum("motif_rings", "density_peaks", "declared")),
                               "enum|null", None,
                               "density_peaks/declared: discretization, not physical site count"),
        },
    },
    "data": {
        "data_root": Key(_is_str, "path"),
        "conditions": RecordList(_CONDITION),
        "selections": {
            "ion": Key(_opt(_is_str), "sel|null", None),
            "ion_negative": Key(_opt(_is_str), "sel|null", None),
            "water_oxygen": Key(_opt(_is_str), "sel|null", None),
            "water_hydrogen": Key(_opt(_is_str), "sel|null", None),
            "phosphate": Key(_opt(_is_str), "sel|null", None),
            "lipid": Key(_opt(_is_str), "sel|null", None),
            "channel": Key(_opt(_is_str), "sel|null", None),
            "anchor": Key(_opt(_is_str), "sel|null", None,
                          "anchor atoms: their center is the zero of axial offsets"),
        },
        "pore": {
            "axis": Key(_opt(_enum("x", "y", "z")), "enum|null", None),
            "cylinder_radius_A": Key(_opt(_is_num), "float|null", None,
                                     "pore cylinder for water/ion column features"),
            "crossing_cylinder_radius_A": Key(_opt(_is_num), "float|null", None,
                                              "part of the task definition, not a geometry detail"),
            "probe_search_radius_A": Key(_is_num, "float", 6.0,
                                         "profile-estimator resolution: midline offset from "
                                         "the axis; diagnostic: share of slices whose optimum "
                                         "sits on the boundary"),
            "probe_grid_dr_A": Key(_is_num, "float", 0.5, "radial step of probe candidate grid"),
            "probe_grid_n_theta": Key(_is_int, "int", 24, "angular probe candidates per ring"),
            "fence_max_gap_deg": Key(_is_num, "float", 120.0,
                                     "maximum angular gap in the fence; a wider gap makes "
                                     "the slice non-searchable"),
            "lipschitz_tol_A": Key(_is_num, "float", 0.05,
                                   "Lipschitz-envelope tolerance: below it the estimate is "
                                   "provably not a radius"),
            "slab_pad_A": Key(_is_num, "float", 8.0,
                              "half-thickness of the slice atom slab; caps the measurable gap"),
            "z_step_A": Key(_is_num, "float", 1.0, "profile slice step; grid anchored to protein"),
            "z_low_offset_A": Key(_opt(_is_num), "float|null", None,
                                  "lower profile bound as an offset from the anchor"),
            "z_high_offset_A": Key(_opt(_is_num), "float|null", None),
            "bin_low_offset_A": Key(_opt(_is_num), "float|null", None,
                                    "bin edges are fixed offsets from the anchor"),
            "bin_width_A": Key(_opt(_is_num), "float|null", None),
            "n_bins": Key(_is_int, "int", 8,
                          "schema constant: same across systems, else tables are incomparable"),
        },
        "stride": Key(_is_int, "int", 1,
                      "part of the run identity; the output folder name carries it"),
        "autodetect_scan": {
            "sample_ps": Key(_is_num, "ps", 200.0,
                             "sampling interval of the detection scan, in TIME units: "
                             "bulk-ion rms displacement per sample must stay far below "
                             "half the cell height, or reservoir-to-reservoir diffusion "
                             "through the periodic seam fabricates ladder passes"),
            "window_ps": Key(_is_num, "ps", 400000.0,
                             "contiguous span scanned from the start of each trajectory; "
                             "detection needs rates and signs, not full statistics"),
        },
    },
    "events": {
        "source": Key(_enum("own", "provided", "auto"), "enum", "auto",
                      "with two sources available, auto REFUSES and prints a comparison"),
    },
    "labels": {
        "tau_ps": Key(_list_of(_is_num), "list[ps]", REQUIRED,
                      "a sweep over horizons is mandatory"),
        "primary_tau_ps": Key(_is_num, "ps", REQUIRED),
        "anchor": Key(_enum("exit", "entrance"), "enum", "exit",
                      "entrance serves as a control ablation for circularity"),
    },
    "features": {
        "blocks": Key(_list_of(_is_str), "list[str]",
                      ["geometry", "hydration", "water_wire", "occupancy", "named_sites",
                       "electrostatics", "symmetry", "fluctuations", "dynamics", "delivery"]),
        "window_ps": Key(_opt(_is_num), "ps|null", None,
                         "windows are set in time units and look strictly backward"),
        "coordination_cutoff_A": Key(_opt(_is_num), "float|null", None,
                                     "first minimum of the ion-water g(r); median over "
                                     "conditions"),
        "partner_cutoff_A": Key(_is_num, "float", 12.0,
                                "cutoff for summing partner Coulomb contributions; "
                                "the column name encodes the cutoff"),
        "site_sigma_A": Key(_opt(_is_num), "float|null", None,
                            "sigma of soft site occupancy; rule: half the site spacing"),
        "site_centers_offset_A": Key(_opt(_list_of(_is_num)), "list[A]|null", None,
                                     "density-peak site centers as offsets from the anchor; "
                                     "only for sites_basis=density_peaks"),
        "lining_band_A": Key(_is_num, "float", 4.0,
                             "half-width of the axial band around the constriction "
                             "for lining columns"),
        "mouth_radius_A": Key(_opt(_is_num), "float|null", None,
                              "permeant counting radius at the mouth (delivery block)"),
        "mouth_depth_A": Key(_is_num, "float", 5.0, "axial extent of the mouth zone"),
    },
    "model": {
        "which": Key(_list_of(_is_str), "list[str]", ["linear", "trees"]),
        "head": Key(_enum("linear", "trees"), "enum", "trees",
                    "the head model is declared IN ADVANCE, before any results are seen"),
        "seed": Key(_is_int, "int", REQUIRED,
                    "one recorded seed; fold seeds are derived from it"),
        "n_threads": Key(_is_int, "int", 1,
                         "fixed for bitwise reproducibility of summations"),
        "linear": {
            "C_grid": Key(_list_of(_is_num), "list[float]",
                          [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0],
                          "grid wide enough that the optimum lies in its interior"),
            "tol": Key(_is_num, "float", 1e-8, "strict tolerance: free with an exact solver"),
            "max_iter": Key(_is_int, "int", 5000),
            "class_weight_mode": Key(_enum("measure", "none", "balanced"), "enum", "measure",
                                     "weighting is decided by measurement, not by default"),
        },
        "trees": {
            "learning_rate": Key(_is_num, "float", 0.05,
                                 "fixed, not tuned: rate and number of rounds are one "
                                 "parameter"),
            "rounds_ceiling": Key(_is_int, "int", 2000),
            "patience": Key(_is_int, "int", 50),
            "depth_grid": Key(_list_of(_is_int), "list[int]", [2, 3, 4],
                              "grid with an interior: the choice must not land on an edge"),
            "k_events_grid": Key(_list_of(_is_num), "list[float]", [2.0, 4.0, 8.0],
                                 "leaf constraint parameterized in events, not as a raw "
                                 "threshold"),
            "colsample_grid": Key(_list_of(_is_num), "list[float]", [0.4, 0.7, 1.0]),
            "subsample_grid": Key(_list_of(_is_num), "list[float]", [0.6, 0.8, 1.0]),
        },
        "plsfma": {
            "components_grid": Key(_list_of(_is_int), "list[int]", list(range(1, 21)),
                                   "PLS latent components, scanned exhaustively as the "
                                   "authors' g_fma -optDim does (1..20); chosen by "
                                   "validation Pearson correlation (Krivobokova et al. "
                                   "2012); 1 is the hard floor of the admissible range"),
        },
        "monotone": Key(_map_str_int, "map[col->±1]", {},
                        "only signs that follow from physics; keyed by column name"),
        "calibration": {
            "family": Key(_enum("sigmoid", "isotonic"), "enum", "sigmoid",
                          "calibrator monotone BY CONSTRUCTION"),
            "min_positives_isotonic": Key(_is_int, "int", 200,
                                          "below this, fall back to the parametric map"),
        },
    },
    "evaluate": {
        "ci_alpha": Key(_is_num, "float", 0.05),
        "n_boot": Key(_is_int, "int", 1000, "paired bootstrap over trajectories"),
        "corr_threshold": Key(_is_num, "float", 0.98,
                              "strong-correlation map threshold; reported, not enforced"),
        "alarm": {
            "smooth_ps": Key(_opt(_is_num), "ps|null", None, "strictly backward smoothing"),
            "on_quantile": Key(_is_num, "float", 0.95, "alarm arming (hysteresis)"),
            "off_quantile": Key(_is_num, "float", 0.80, "alarm release (hysteresis)"),
        },
    },
    "interpret": {
        "negligible_frac": Key(_is_num, "float", 0.05,
                               "axis contribution negligible against baseline -> 'unclear'"),
    },
    "transfer": {
        "max_importance_loss_frac": Key(_is_num, "float", 0.5,
                                        "if more importance is lost, refuse to answer"),
    },
    "accept": {
        "edge_fraction_max": Key(_is_num, "float", 0.34,
                                 "fraction of folds choosing a grid edge; above is a defect"),
        "missing_budget_frac": Key(_is_num, "float", 0.5,
                                   "missingness budget AFTER masking, over active columns"),
    },
}


def _walk_validate(node: Any, spec: Any, path: str, out: dict) -> Any:
    """Validated copy of `node` against `spec`, with absent keys filled from defaults.

    `path` is the dotted position in the tree and appears in every error message,
    so a rejection names the key a person has to go and fix.
    """
    if isinstance(spec, Key):
        # An explicit null is only accepted where the type says so (the _opt
        # wrapper). Writing null under a key that has a default would otherwise
        # read as "no opinion" and silently become the default.
        if node is None and spec.default is not REQUIRED and not spec.check(None):
            raise ConfigError(f"{path}: null is not allowed")
        if not spec.check(node):
            raise ConfigError(f"{path}: expected {spec.typename}, got {node!r}")
        # A copy: the validated tree must not alias the parsed document, or a
        # mutable value -- the list of horizons, the monotone map -- could be
        # changed through the other reference after validation had passed.
        return copy.deepcopy(node)
    if isinstance(spec, RecordList):
        if not isinstance(node, list):
            raise ConfigError(f"{path}: expected a list of records")
        return [_walk_validate(item, spec.spec, f"{path}[{i}]", out)
                for i, item in enumerate(node)]
    if isinstance(spec, dict):
        if not isinstance(node, dict):
            raise ConfigError(f"{path}: expected a section")
        # The rule the whole module exists for. A misspelt key would otherwise be
        # ignored, the schema default would be used, and the run would look
        # successful while measuring something other than what was asked for.
        unknown = set(node) - set(spec)
        if unknown:
            raise ConfigError(
                f"{path}: unknown keys {sorted(unknown)} -- a typo must not fall "
                f"back silently to the default")
        result = {}
        for name, sub in spec.items():
            child_path = f"{path}.{name}" if path else name
            if name in node:
                result[name] = _walk_validate(node[name], sub, child_path, out)
            else:
                result[name] = _fill_defaults(sub, child_path)
        return result
    raise AssertionError(f"malformed schema at {path}")


def _fill_defaults(spec: Any, path: str) -> Any:
    """Subtree of defaults for a section the file does not mention."""
    if isinstance(spec, Key):
        if spec.default is REQUIRED:
            raise ConfigError(f"{path}: required key is missing")
        # Deep copy per config: several defaults are mutable module-level objects
        # (the block list, the C grid, the depth grid), and handing out the object
        # itself would let one run's edit reach the next config loaded in the same
        # process.
        return copy.deepcopy(spec.default)
    if isinstance(spec, RecordList):
        # Conditions and replicas are measurements of the dataset, not settings;
        # there is no defensible default for "which trajectories exist".
        raise ConfigError(f"{path}: required record list is missing")
    return {name: _fill_defaults(sub, f"{path}.{name}") for name, sub in spec.items()}


class Config:
    """Validated config plus provenance ledger.

    The same object is handed to every step and copied verbatim into each
    step's PROVENANCE.json, so the values recorded next to an output are the
    values that produced it. source_path is kept because the autodetect step
    writes its detected keys back into the file it was loaded from.
    """

    def __init__(self, tree: dict, origins: dict[str, dict], source_path: Path | None = None):
        self._tree = tree
        self._origins = origins
        self.source_path = source_path

    def __getitem__(self, dotted: str) -> Any:
        """Value at a dotted path, e.g. cfg["data.pore.axis"]; KeyError if there is none.

        Raising is the point: a mistyped path in the code fails at the line that
        reads it instead of returning None and being treated as "not detected".
        """
        node: Any = self._tree
        for part in dotted.split("."):
            if isinstance(node, dict):
                node = node[part]
            else:
                raise KeyError(dotted)
        return node

    def get(self, dotted: str, default: Any = None) -> Any:
        """Value at a dotted path, or `default` where a missing path is genuinely allowed."""
        try:
            return self[dotted]
        except KeyError:
            return default

    @property
    def tree(self) -> dict:
        """The whole validated tree, as written into PROVENANCE.json."""
        return self._tree

    @property
    def origins(self) -> dict[str, dict]:
        """The provenance ledger: dotted key -> {origin, basis}.

        Only keys that were detected or declared appear here. Record paths are
        spelled with the record id in brackets, e.g.
        data.conditions[200MV].direction, and are written straight into this
        mapping because set_key only reaches leaves of the nested sections.
        """
        return self._origins

    def origin_of(self, dotted: str) -> tuple[str, str]:
        """(origin, basis) for a key; absence from the ledger means the schema default."""
        ent = self._origins.get(dotted)
        if ent is None:
            return "default", "schema value, never overridden"
        return ent["origin"], ent.get("basis", "")

    def set_key(self, dotted: str, value: Any, origin: str, basis: str) -> None:
        """Set a leaf key and record where the value came from.

        A basis must be passed with every write, which is what keeps the
        artifact answerable: no value enters a config without a stated reason.
        The value is re-checked against the schema, so a detected quantity is
        held to the same type as a hand-written one. Only leaves of the nested
        sections can be set this way; record fields (a condition's direction, a
        replica's protocol) are written into the tree and into `origins`
        directly by the step that measures them.
        """
        if origin not in ORIGINS:
            raise ConfigError(f"origin {origin!r} not in {ORIGINS}")
        spec: Any = SCHEMA
        node: Any = self._tree
        parts = dotted.split(".")
        for part in parts[:-1]:
            spec = spec[part] if isinstance(spec, dict) else spec
            node = node[part]
        leaf = parts[-1]
        if not (isinstance(spec, dict) and leaf in spec and isinstance(spec[leaf], Key)):
            raise ConfigError(f"{dotted}: not a leaf key of the schema")
        if not spec[leaf].check(value):
            raise ConfigError(f"{dotted}: value {value!r} fails the type check")
        node[leaf] = value
        self._origins[dotted] = {"origin": origin, "basis": basis}

    def replicas(self):
        """(condition_id, replica_dict) in the stable order of the file.

        Steps that fan work out to a process pool collect the results and then
        re-emit them in this order, so the row order of every artifact is the
        order of the config and not the order in which workers happened to
        finish.
        """
        for cond in self["data.conditions"]:
            for rep in cond["replicas"]:
                yield cond["id"], rep


def _cross_checks(cfg: Config) -> None:
    """Agreements between keys that a per-key type check cannot see."""
    taus = cfg["labels.tau_ps"]
    # The headline horizon has to be one of the horizons actually labelled, or
    # the primary result would be reported from columns that were never built.
    if cfg["labels.primary_tau_ps"] not in taus:
        raise ConfigError("labels.primary_tau_ps must be listed in labels.tau_ps")
    # The head model is declared in advance; it must be among the arms that are
    # trained, so the declaration cannot be quietly satisfied after the fact.
    if cfg["model.head"] not in cfg["model.which"]:
        raise ConfigError("model.head must be listed in model.which")
    known_blocks = set(SCHEMA["features"]["blocks"].default)
    extra = set(cfg["features.blocks"]) - known_blocks
    if extra:
        raise ConfigError(f"features.blocks: unknown blocks {sorted(extra)}")
    # condition/replica is the trajectory key: it names the group used for
    # splitting and the rows of every table. A repeated pair would silently fuse
    # two trajectories into one unit and inflate the effective sample size.
    seen = set()
    for cond in cfg["data.conditions"]:
        for rep in cond["replicas"]:
            key = (cond["id"], rep["id"])
            if key in seen:
                raise ConfigError(f"duplicate replica {key}")
            seen.add(key)


def load(path: str | Path) -> Config:
    """Read, validate and cross-check one config file; raises ConfigError on any defect."""
    path = Path(path)
    raw = yaml.safe_load(path.read_text()) or {}
    # `origins` is the ledger, not part of the schema tree: it is taken out
    # before validation, or the unknown-key rule would reject it as a stray
    # top-level section.
    origins = raw.pop("origins", {}) or {}
    # The ledger is checked before any value is looked at. An entry without a
    # non-empty basis is refused: a value whose reason has been left blank is
    # exactly the thing this mechanism is meant to prevent.
    for dotted, ent in origins.items():
        if not isinstance(ent, dict) or ent.get("origin") not in ORIGINS:
            raise ConfigError(f"origins.{dotted}: origin must be one of {ORIGINS}")
        if not str(ent.get("basis", "")).strip():
            raise ConfigError(f"origins.{dotted}: empty basis -- a value must carry provenance")
    tree = _walk_validate(raw, SCHEMA, "", {})
    cfg = Config(tree, origins, source_path=path)
    _cross_checks(cfg)
    return cfg


def _prune_defaults(node: Any, spec: Any, origins: dict, prefix: str) -> Any:
    """Only keys with a recorded provenance are written to the file; defaults
    stay in the schema and follow it.

    The effect is that a config file reads as the list of decisions taken for
    this system, and that changing a default in the schema changes it for every
    system at once instead of leaving copies behind in five files.
    """
    if isinstance(spec, Key):
        # Required keys have no default -- they are always written.
        keep = prefix in origins or spec.default is REQUIRED
        return copy.deepcopy(node) if keep else None
    if isinstance(spec, RecordList):
        return copy.deepcopy(node)  # records (conditions/replicas) are data and carry no default
    out = {}
    for name, sub in spec.items():
        p = f"{prefix}.{name}" if prefix else name
        if name not in node:
            continue
        kept = _prune_defaults(node[name], sub, origins, p)
        if kept is not None and kept != {}:
            out[name] = kept
    return out


def save(cfg: Config, path: str | Path | None = None) -> Path:
    """Rewrite the config file in place, keeping only decided keys; returns the path."""
    path = Path(path or cfg.source_path)
    doc = _prune_defaults(cfg.tree, SCHEMA, cfg.origins, "")
    # Sorted, so the ledger does not carry the order in which keys happened to be
    # written and a diff between two runs shows only what actually changed.
    doc["origins"] = {k: dict(v) for k, v in sorted(cfg.origins.items())}
    # Written to a sibling and renamed: the autodetect step rewrites the very
    # file the run was started from, and a crash halfway through the dump would
    # otherwise leave a truncated config behind that no later step could read.
    # sort_keys=False keeps the schema's own order, which is the reading order.
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(yaml.safe_dump(doc, allow_unicode=True, sort_keys=False, width=100))
    tmp.replace(path)
    return path


def iter_leaf_paths(spec: Any = None, prefix: str = "") -> list[str]:
    """Every leaf path of the schema in declaration order.

    A record list counts as one leaf: the paths inside a record depend on how
    many conditions and replicas a given dataset has, so they are not part of
    the schema's own vocabulary.
    """
    spec = SCHEMA if spec is None else spec
    out: list[str] = []
    for name, sub in spec.items():
        p = f"{prefix}.{name}" if prefix else name
        if isinstance(sub, Key):
            out.append(p)
        elif isinstance(sub, RecordList):
            out.append(p)
        else:
            out.extend(iter_leaf_paths(sub, p))
    return out
