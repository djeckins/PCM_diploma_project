"""Assemble every system's initial config from the files on disk.

Only declared keys with recorded bases are hand-written here: anything the
system can measure by itself is left out — auto-detection fills it in. Each
declared value is a (value, basis) pair and the basis travels with it into the
config and from there into every PROVENANCE.json, so a reader can ask of any
number why it has that value and get an answer that was written before the
result was known.

    python tools/init_configs.py           # all systems
    python tools/init_configs.py kcsa      # one system by id

Writes configs/<system>.yaml and reloads each file through the strict validator,
so a typo in a key name fails here rather than in the middle of a run.

Requirements: the trajectory datasets themselves. The replica list, the file
names and the topologies are discovered by walking the data root, so this step
cannot run on the published tree — the configs it would produce are committed
instead. Point PCM2_DATA_ROOT at the directory holding the dataset folders before
rerunning it on a new machine. No trajectory is opened: only the file tree is
read.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pcm2 import config  # noqa: E402
from pcm2.inputs import build_config_doc, write_config  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
# Where the trajectory datasets live. On a new machine, point PCM2_DATA_ROOT
# at the directory holding CHARMM-ECC-078, Shared_KCl_500mV, etc., then rerun
# this script — the generated configs absorb the machine-specific paths.
DIPLOM = Path(os.environ.get("PCM2_DATA_ROOT", str(ROOT.parent)))

TAU_BASIS = ("a sweep over horizons is required so conclusions are not tied to one "
             "arbitrary τ; values are cross-checked against the median transit time and "
             "the inter-crossing interval that the events step prints next to the "
             "declared τ")
CYL_BASIS = ("task definition: separates 'crossed through the pore' from 'crossed the "
             "bilayer past the channel'; 10 A covers the vestibules of both channels "
             "and is well under half the protein-protein spacing in the membrane")
LINEAGE_BASIS = ("declared: each replica is an independent assembly. In the 500mV set "
                 "the replicas carry distinct measured currents in their file names "
                 "(independent runs); MthK replicas sit in separate run folders. If a "
                 "common equilibrated ancestor is known, merge the lineage manually")
SEED_BASIS = "a single recorded run seed; fold seeds are derived from it"
WINDOW_BASIS = ("look-back window equals the primary horizon: dynamics are measured on "
                "the same time scale on which the prediction is made")
MOUTH_BASIS = ("two hydration shells plus the mouth vestibule; diagnostics: the "
               "dlv_entry_wide_n indicator and the fraction of frames with no permeant "
               "in the zone")
SMOOTH_BASIS = "half the primary τ: smoothing stays shorter than the prediction horizon"
MONO_BASIS = ("from physics, not from correlation: a wider pore cannot make a "
              "configuration less ready to pass a hydrated ion")
STRIDE_BASIS = ("development compute budget; the stride is part of run identity, the "
                "results folder carries it in its name")
TEMP_TPR_BASIS = ("ref-t of the thermostat read from the run tpr with gmx dump "
                  "(GROMACS in the separate 'gmx' conda environment)")


def _trees_flat_basis(loss_range: str, fold_spread: str) -> str:
    """Basis text for pinning the boosting grid, quoting that system's own numbers.

    loss_range is the spread of the selection criterion within a fold across the
    whole grid and fold_spread is the spread of the same criterion between folds.
    When the first is an order of magnitude smaller than the second, the criterion
    cannot tell the grid points apart and whichever point it picks is fold noise;
    pinning the parameter and recording why is the honest response. Both numbers
    are read out of the searched run's folds.json.
    """
    return ("fixed deliberately: a 3x3x3x3 grid search showed a flat selection-criterion "
            f"curve (within-fold log-loss range {loss_range} vs between-fold spread "
            f"{fold_spread}, see grid_scores in the searched run's folds.json); the "
            "criterion does not select this parameter, an edge pick would be noise "
            "rather than an optimum")


def _trees_grid_declared(basis: str) -> dict:
    """The four boosting grids as single-element lists, all carrying one basis.

    A one-element grid is how a pinned hyperparameter is expressed: the search
    machinery still runs, records the value and reports it, so nothing pretends a
    choice was made where none was.
    """
    return {
        "model.trees.depth_grid": ([3], basis),
        "model.trees.k_events_grid": ([4.0], basis),
        "model.trees.colsample_grid": ([0.7], basis),
        "model.trees.subsample_grid": ([0.8], basis),
    }


COMMON_DECLARED = {
    "labels.tau_ps": ([1000.0, 2000.0, 4000.0], TAU_BASIS),
    "labels.primary_tau_ps": (2000.0, TAU_BASIS),
    "events.source": ("own", "a single time anchor shared by both systems; the shipped "
                             "gramicidin annotation is used as a cross-check — the "
                             "events step prints the comparison"),
    "data.pore.crossing_cylinder_radius_A": (10.0, CYL_BASIS),
    "features.window_ps": (2000.0, WINDOW_BASIS),
    "features.mouth_radius_A": (12.0, MOUTH_BASIS),
    "model.seed": (20260814, SEED_BASIS),
    "model.monotone": ({"geo_r_constriction_A": 1}, MONO_BASIS),
    "evaluate.alarm.smooth_ps": (1000.0, SMOOTH_BASIS),
}

SYSTEMS = {
    "mthk": {
        "data_root": DIPLOM / "CHARMM-ECC-078",
        "stride": (5, STRIDE_BASIS + "; 50 ps × 5 = 250 ps between table frames"),
        "declared": _trees_grid_declared(_trees_flat_basis("0.0017-0.0034", "0.061")),
        "protocol": {
            "200MV": {
                "forcefield": ("CHARMM36", "dataset name CHARMM-ECC-078"),
                "voltage_mV": (200.0, "condition folder name 200MV"),
                "temperature_K": (323.15, TEMP_TPR_BASIS),
            },
            "neg-200MV": {
                "forcefield": ("CHARMM36", "dataset name CHARMM-ECC-078"),
                "voltage_mV": (-200.0, "condition folder name neg-200MV"),
                "temperature_K": (323.15, TEMP_TPR_BASIS),
            },
        },
    },
    "gramicidin": {
        "data_root": DIPLOM / "Shared_KCl_500mV",
        "stride": (4, STRIDE_BASIS + "; 40 ps × 4 = 160 ps between table frames"),
        "declared": _trees_grid_declared(_trees_flat_basis("0.001-0.0035", "0.094")),
        "protocol": {
            "330K_1Msalt": {
                "voltage_mV": (500.0, "dataset name Shared_KCl_500mV"),
                "temperature_K": (330.0, TEMP_TPR_BASIS),
            },
            "298K_2Msalt": {
                "voltage_mV": (500.0, "dataset name Shared_KCl_500mV"),
                "temperature_K": (298.15, TEMP_TPR_BASIS),
            },
        },
    },
    "kcsa": {
        "data_root": DIPLOM / "3 белок" / "kcsa_root",
        "stride": (4, STRIDE_BASIS + "; 40 ps × 4 = 160 ps between table frames"),
        "declared": {
            # A new system searches its hyperparameters: flatness measured on the
            # other systems is their property, not a transferable constant.
            "labels.tau_ps": ([500.0, 1000.0, 2000.0],
                              TAU_BASIS + "; at 300 mV conduction is fast, so the sweep "
                                          "is shifted toward shorter horizons"),
            "labels.primary_tau_ps": (1000.0, TAU_BASIS),
            "features.window_ps": (1000.0, WINDOW_BASIS),
            "evaluate.alarm.smooth_ps": (500.0, SMOOTH_BASIS),
            "system.arch_profile.filter_motif_pattern": (
                "[TA][VIL][GA][YF]G",
                "the dataset consists of engineered selectivity-filter mutants on the "
                "E71A background: T75A turns TVGYG into AVGYG and G77A turns it into "
                "TVAYG by construction; the widened signature admits those single "
                "substitutions while the geometric ring confirmation still gates "
                "detection"),
        },
        "protocol": {
            m: {
                "forcefield": ("CHARMM36m", "directory name CHARMM36m in the input archive"),
                "voltage_mV": (300.0, "directory name 300mV in the input archive"),
                "temperature_K": (290.0, TEMP_TPR_BASIS),
            } for m in ("E71A", "G77AE71A", "T75AE71A")
        },
    },
}

# The Na arm of the paired K/Na selectivity design: same three filter variants,
# NaCl instead of KCl. A near-zero event count in the selective background is
# the expected physics, not a data defect; such units serve as test-only units.
SYSTEMS["kcsa_na"] = {
    "data_root": DIPLOM / "3 белок" / "kcsa_na_root",
    "stride": SYSTEMS["kcsa"]["stride"],
    "declared": dict(SYSTEMS["kcsa"]["declared"]) | {
        "system.permeant": ("SOD",
            "declared up front: this is the Na arm of the paired K/Na design; water "
            "would out-cross the rare Na permeant per particle for the same reason "
            "it does in the K arm, and readiness is defined for the ionic carrier"),
    },
    "protocol": SYSTEMS["kcsa"]["protocol"],
}


# Connexin-43 gap junction (dodecameric connexon pair, double bilayer): a wide
# aqueous pore, structurally unrelated to the K-channel family. The topology is
# a GROMACS .top assembled from the authors' SI (initial_structure.gro matches
# the trajectories atom for atom; the water count in the SI topol.top predates
# the final solvent trim and was corrected to the gro composition).
SYSTEMS["cx43"] = {
    "data_root": DIPLOM / "protein_center_200mV" / "cx43_root",
    "stride": (2, STRIDE_BASIS + "; 50 ps × 2 = 100 ps between table frames"),
    "declared": {
        "system.permeant": ("POT",
            "declared up front: a wide aqueous pore passes water molecules per "
            "particle orders of magnitude more often than ions; readiness is "
            "defined for the ionic charge carrier K+ driven by the applied bias"),
        "data.pore.crossing_cylinder_radius_A": (20.0,
            "the connexon vestibule opens to ~35 A diameter while the protein "
            "footprint is ~90 A across; 20 A covers the funnel mouth and stays "
            "well inside the dodecamer wall"),
        "features.mouth_radius_A": (20.0,
            "the cytoplasmic funnel of the connexon is far wider than a "
            "K-channel mouth; the delivery zone must cover it, diagnostics: "
            "the dlv_entry_wide_n indicator and the fraction of frames with "
            "no permeant in the zone"),
    },
    "protocol": {
        "200mV": {
            "forcefield": ("CHARMM36", "CHARMM-GUI FF-Converter headers in the "
                                       "SI toppar files"),
            "voltage_mV": (200.0, "folder and mdp name efield_200mV; the shipped "
                                  "mdp applies electric-field-z = -0.008 V/nm, "
                                  "i.e. cations drift toward -z"),
            "temperature_K": (303.15, "ref_t in the shipped production mdp"),
        },
    },
    "lineage_basis": ("declared: R2-R4 are separate 100 ns production runs from "
                      "one assembled system (the SI ships a single "
                      "initial_structure); they share no production frames and "
                      "diverge at production start, so they are kept as separate "
                      "lineages for splitting; the shared assembly is recorded "
                      "here so a reviewer can merge them if stricter separation "
                      "is wanted"),
}


def main(only: list[str] | None = None) -> None:
    """Write one config per system, or only the ids listed on the command line.

    Per-system declarations override the common ones, which is how KcsA gets its
    shorter horizon sweep and Cx43 its wider cylinder without those values leaking
    into the other systems.
    """
    for sysid, spec in SYSTEMS.items():
        if only and sysid not in only:
            continue
        declared = dict(COMMON_DECLARED)
        declared["data.stride"] = spec["stride"]
        declared.update(spec["declared"])
        doc = build_config_doc(sysid, spec["data_root"], declared,
                               spec.get("lineage_basis", LINEAGE_BASIS),
                               protocol_declared=spec["protocol"])
        path = ROOT / "configs" / f"{sysid}.yaml"
        write_config(doc, path)
        cfg = config.load(path)  # strict validation: a typo must fail here
        n_reps = sum(len(c["replicas"]) for c in cfg["data.conditions"])
        print(f"{sysid}: {path} — {len(cfg['data.conditions'])} conditions, {n_reps} replicas; "
              f"validation passed")


if __name__ == "__main__":
    main(sys.argv[1:] or None)
