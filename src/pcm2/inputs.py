"""File discovery on disk and config assembly from what is found.

Recognized layouts:
  A) <root>/<condition>/<replica>/*.xtc — topology (.tpr or GROMACS .top) in
     the replica directory itself, or, when there is none there, the single one
     found anywhere under the root;
  B) <root>/<condition>/NN-*.xtc with a same-stem NN-*.tpr / NN-*.top, or that
     single topology under the root, and, when present, NN-*transit*.dat.
Discovered paths are written with origin "detected"; what is not on disk
(horizons, seed, crossing-cylinder radius) is supplied by a human and written
with origin "declared".

Nothing here reads a trajectory. The point of discovering the file list instead
of typing it is that a hand-written path can point at a neighbouring replica and
still open, whereas a detected list carries a basis string saying what was found
under which root, and the count of replicas per condition can be checked against
the dataset at a glance.
"""

from __future__ import annotations

from pathlib import Path

import yaml


# Order is preference. A .tpr is the actual input of the run that wrote the
# trajectory, so it is the topology guaranteed to match the frames; the GROMACS
# text topology is the fallback for datasets shipped without one.
TOPOLOGY_EXTS = (".tpr", ".top")


def _unique_topology(root: Path) -> Path | None:
    """The one topology under `root`, or None when there are none or several.

    Datasets that share a single built system among all replicas are common, and
    for those the shared file is the right answer. Several candidates are not
    resolved by guessing: pairing a trajectory with a topology from another build
    gives a different atom order, and every selection would then address the
    wrong atoms without any file failing to open.
    """
    tops = [p for ext in TOPOLOGY_EXTS for p in root.rglob(f"*{ext}")
            if not p.name.startswith(".")]
    return tops[0] if len(tops) == 1 else None


def discover(data_root: Path) -> tuple[list[dict], str]:
    """(condition records, basis string) for everything found under one data root.

    A directory holding no trajectory is not reported as an empty condition: it
    simply does not appear, which is how analysis folders sitting beside the data
    stay out of the config. Names beginning with a dot are skipped everywhere --
    they are what an interrupted copy or a filesystem leaves behind, not data.
    Everything is walked in sorted order, so the same tree gives the same config
    on any machine.
    """
    data_root = Path(data_root)
    conditions: list[dict] = []
    notes: list[str] = []
    for cond_dir in sorted(p for p in data_root.iterdir()
                           if p.is_dir() and not p.name.startswith(".")):
        replicas: list[dict] = []
        # Trajectories lying directly in the condition folder mean the flat
        # layout; the presence of at least one decides which of the two layouts
        # this condition is read as.
        direct_xtc = sorted(p for p in cond_dir.glob("*.xtc") if not p.name.startswith("."))
        if direct_xtc:
            for xtc in direct_xtc:
                # Flat layout: the replica id is the part of the file name before
                # the first dash -- a number in these datasets -- so that
                # 01-....xtc and 01-...transit.dat belong to the same replica.
                rep_id = xtc.stem.split("-")[0]
                # Same stem, topology suffix: the pairing the run itself used.
                tpr = next((c for c in (xtc.with_suffix(e) for e in TOPOLOGY_EXTS)
                            if c.exists()), None)
                if tpr is None:
                    tpr = _unique_topology(data_root)
                # Annotation files are matched within the replica prefix, so one
                # replica's crossings can never be attached to another. Where
                # several match, the last in sorted order is kept.
                ev = None
                for cand in sorted(cond_dir.glob(f"{rep_id}-*transit*.dat")):
                    ev = cand
                # A trajectory with no topology is skipped rather than guessed
                # at, and the note goes into the basis string: the config then
                # records that a file on disk was deliberately left out.
                if tpr is None:
                    notes.append(f"{xtc}: no topology found — replica skipped")
                    continue
                replicas.append({"id": rep_id, "topology": str(tpr),
                                 "trajectory": str(xtc),
                                 "events": str(ev) if ev else None})
        else:
            # Nested layout: one directory per replica, and the directory name is
            # the replica id, which is what appears in every table afterwards.
            for rep_dir in sorted(p for p in cond_dir.iterdir()
                                  if p.is_dir() and not p.name.startswith(".")):
                xtcs = sorted(p for p in rep_dir.glob("*.xtc") if not p.name.startswith("."))
                if not xtcs:
                    continue
                tprs = sorted(p for ext in TOPOLOGY_EXTS
                              for p in rep_dir.glob(f"*{ext}")
                              if not p.name.startswith("."))
                # A topology inside the replica directory wins over the shared one:
                # it belongs to this build, not to the dataset in general.
                tpr = tprs[0] if tprs else _unique_topology(data_root)
                # As above, the last match in sorted order is kept; the *events*
                # pattern is searched after *transit*.dat, so it takes precedence.
                ev = None
                for cand in sorted(rep_dir.glob("*transit*.dat")) + sorted(rep_dir.glob("*events*")):
                    ev = cand
                if tpr is None:
                    notes.append(f"{rep_dir}: no topology found — replica skipped")
                    continue
                replicas.append({"id": rep_dir.name, "topology": str(tpr),
                                 "trajectory": str(xtcs[0]),
                                 "events": str(ev) if ev else None})
        if replicas:
            conditions.append({"id": cond_dir.name, "replicas": replicas})
    # The basis is what the config will carry as the provenance of its file list:
    # the root that was walked, the replica count per condition, and every
    # skipped file. It has to be readable on its own, because it is the only
    # record of the discovery once the machine that held the data is gone.
    basis = (f"discovered on disk under {data_root}: "
             + "; ".join(f"{c['id']}: {len(c['replicas'])} replicas" for c in conditions)
             + ("; " + "; ".join(notes) if notes else ""))
    return conditions, basis


def build_config_doc(system_id: str, data_root: str | Path,
                     declared: dict, lineage_basis: str,
                     protocol_declared: dict | None = None) -> dict:
    """Assemble the YAML config document: paths are detected, the rest is declared.

    declared: {dotted_key: (value, basis)} — what cannot be measured.
    protocol_declared: {condition_id: {protocol key: (value, basis)}}.

    Nothing is written for a key that neither branch supplies: the schema default
    then applies, and the absence of the key from the file is itself the record
    that no decision was taken about it.
    """
    conditions, disc_basis = discover(Path(data_root))
    origins: dict[str, dict] = {}
    doc: dict = {"system": {"id": system_id},
                 "data": {"data_root": str(data_root), "conditions": conditions},
                 "labels": {}, "model": {}}
    origins["data.data_root"] = {"origin": "detected", "basis": disc_basis}
    origins["data.conditions"] = {"origin": "detected", "basis": disc_basis}
    for cond in conditions:
        for rep in cond["replicas"]:
            # Each replica starts out as a lineage of its own, i.e. as an
            # independently built assembly. This is a declaration and not a
            # measurement: where replicas were branched from one equilibrated
            # structure, their lineages have to be merged by hand, or the split
            # would put the same starting configuration on both sides and the
            # reported score would be optimistic. The reason for the choice
            # travels with the value in lineage_basis.
            rep["lineage"] = f"{cond['id']}/{rep['id']}"
            # Record paths address a record by its id rather than by position, so
            # that adding a replica does not shift the provenance of the others.
            origins[f"data.conditions[{cond['id']}].replicas[{rep['id']}].lineage"] = {
                "origin": "declared", "basis": lineage_basis}
            # The protocol is declared per condition but stored per replica: the
            # replica is the unit that is opened, scored and compared, so its
            # record has to state the physics it was run under without a lookup.
            rep["protocol"] = {}
            for pk, pv in (protocol_declared or {}).get(cond["id"], {}).items():
                rep["protocol"][pk] = pv[0]
                origins[f"data.conditions[{cond['id']}].replicas[{rep['id']}].protocol.{pk}"] = {
                    "origin": "declared", "basis": pv[1]}
    for dotted, (value, basis) in declared.items():
        node = doc
        parts = dotted.split(".")
        for p in parts[:-1]:
            node = node.setdefault(p, {})
        node[parts[-1]] = value
        origins[dotted] = {"origin": "declared", "basis": basis}
    doc["origins"] = origins
    return doc


def write_config(doc: dict, path: str | Path) -> Path:
    """Write the assembled document; returns the path.

    Keys are dumped in the order they were assembled rather than alphabetically,
    which keeps the file in reading order: what the system is, then where its
    data is, then what was decided about it.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(doc, allow_unicode=True, sort_keys=False, width=100))
    return path
