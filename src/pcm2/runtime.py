"""Run identity, atomic step output, provenance.

There is no cache in the system: a step, once launched, recomputes and replaces
its directory wholesale, so a step directory always holds the full set of files.
The alternative -- reusing whichever files happen to be newer than their inputs
-- would mean a published table could mix columns computed under two different
versions of an estimator, and nothing in the directory would show it.

A run is identified by the system id together with the frame stride, because a
stride is a measurement decision and not a speed setting: it sets the spacing of
the frame grid on which labels, windows and dynamics are computed. Two strides
therefore live in two directories and are never compared inside one.

Provenance is descriptive, not a control mechanism. PROVENANCE.json records the
config, its origins ledger, the library versions and a hash of the source tree
next to the output; it is never consulted to decide whether something needs
recomputing.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import sys
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from .config import Config

# This file is src/pcm2/runtime.py, so two levels up is the repository root.
# Everything the code writes or reads by convention (runs/, external/) is located
# from here rather than from the working directory, so a command gives the same
# result from any directory. Moving this module to another depth would silently
# point the whole pipeline at the wrong tree.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
# Single source of truth for the path to the vendored Rao-2019 material (the
# dewetting surface, the hydrophobicity scale, the licence). The feature layer
# and the report step both read it from here.
VENDORED_DIR = PROJECT_ROOT / "external" / "rao2019_heuristic"

# The step names step_output will accept. A name outside this set raises rather
# than creating a directory: a typo would otherwise put a plausible-looking but
# unreferenced folder next to the real artifacts.
STEPS = ("autodetect", "events", "features", "coords", "labels", "train",
         "figures", "report", "benchmark")


def run_dir(cfg: Config) -> Path:
    """Directory holding every step of one run: runs/<system>-stride<N>/."""
    # The frame stride is part of the run identity: the directory carries it in
    # its name.
    return PROJECT_ROOT / "runs" / f"{cfg['system.id']}-stride{cfg['data.stride']}"


def code_fingerprint() -> dict[str, str]:
    """SHA-256 of every module of the package, plus one hash over the whole tree.

    This is what makes an output attributable to a source state: if a figure and
    a table were produced under different code, their tree hashes differ. The
    files are sorted before hashing so the tree hash depends on the contents and
    not on the order the filesystem lists them in. Only the package sources are
    covered -- the config and the library versions are recorded separately in the
    same document.
    """
    src = Path(__file__).resolve().parent
    files = sorted(p for p in src.rglob("*.py"))
    per_file = {str(p.relative_to(src)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in files}
    whole = hashlib.sha256("".join(f"{k}:{v}" for k, v in per_file.items()).encode()).hexdigest()
    return {"tree_sha256": whole, "files": per_file}


def library_versions() -> dict[str, str]:
    """Versions of the numerical stack, as they were at the moment of the run.

    A library that cannot be imported is recorded as "absent" instead of raising:
    a step that does not need the model stack must still be able to write its
    provenance, and the absence is itself information about the environment.
    """
    out = {"python": sys.version.split()[0]}
    for name in ("numpy", "scipy", "pandas", "MDAnalysis", "sklearn", "xgboost",
                 "numba", "matplotlib", "yaml", "pyarrow"):
        try:
            mod = __import__(name)
            out[name] = getattr(mod, "__version__", "?")
        except Exception:
            out[name] = "absent"
    return out


def write_provenance(step_dir: Path, cfg: Config, step: str, extra: dict | None = None) -> None:
    """Record of what an output was computed with; it invalidates nothing.

    Written into the step's temporary directory just before it is published, so
    a directory that exists always has its record. default=str lets a value the
    JSON encoder does not know (a Path, a numpy scalar) be recorded as text
    rather than aborting a finished computation at the last line.
    """
    doc = {
        "step": step,
        "system": cfg["system.id"],
        "stride": cfg["data.stride"],
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "platform": {"platform": platform.platform(), "machine": platform.machine()},
        "n_threads": cfg["model.n_threads"],
        "libraries": library_versions(),
        "code": code_fingerprint(),
        "config": cfg.tree,
        "origins": cfg.origins,
    }
    if extra:
        doc["extra"] = extra
    (step_dir / "PROVENANCE.json").write_text(
        json.dumps(doc, indent=1, ensure_ascii=False, default=str))


@contextmanager
def step_output(cfg: Config, step: str):
    """Write to a temporary location and replace wholesale; on failure touch nothing.

    Yields the directory the step must write into. The caller never sees the
    final path, which is the whole point: a half-written step is invisible to the
    steps downstream, and the previous output stays readable until the new one is
    complete.
    """
    if step not in STEPS:
        raise ValueError(f"unknown step {step!r}")
    root = run_dir(cfg)
    root.mkdir(parents=True, exist_ok=True)
    final = root / step
    # Clear orphans left by killed runs of this step.
    # A process killed outright cannot run its own cleanup, so the debris is
    # collected at the start of the next attempt instead. The dot prefix marks it
    # as not an artifact: a shell glob over the run folder and the listing a person
    # reads both pass it over.
    for stale in root.glob(f".{step}.tmp-*"):
        shutil.rmtree(stale, ignore_errors=True)
    # The pid keeps two attempts at the same step out of one another's directory.
    # Running the same step twice at once on one run folder is outside the design
    # -- the sweep above would take the other process's directory for debris --
    # which is why a run is never overlapped with another run of itself.
    tmp = root / f".{step}.tmp-{os.getpid()}"
    tmp.mkdir()
    try:
        yield tmp
    except BaseException:
        # BaseException, not Exception: a KeyboardInterrupt in the middle of a
        # multi-hour step must also leave no partial directory behind. The error
        # is re-raised unchanged; this only removes the debris.
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    write_provenance(tmp, cfg, step)
    if final.exists():
        shutil.rmtree(final)
    tmp.rename(final)


def derive_seed(base_seed: int, *tags: str) -> int:
    """Fold seed derived by hash from the single recorded base seed.

    One seed is declared in the config and every per-fold seed follows from it
    and the fold's tags, so the whole set of seeds is reconstructible from the
    artifact and no seed was chosen after seeing a result. Four bytes keep the
    value inside the 32-bit range the numerical libraries accept.
    """
    h = hashlib.sha256(("|".join([str(base_seed), *tags])).encode()).digest()
    return int.from_bytes(h[:4], "little")


def pin_threads(cfg: Config) -> None:
    """Fix the thread count of every numerical backend to model.n_threads."""
    # A parallel reduction adds its partial sums in whatever order the threads
    # finish, so the low-order digits of a sum depend on the thread count. With
    # the count fixed the same input gives the same number down to the last bit.
    # All five variables are set because the backends read different ones:
    # OpenMP, OpenBLAS, MKL, Apple's vecLib and numexpr.
    n = str(cfg["model.n_threads"])
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[var] = n


class StepLog:
    """Print + per-step journal: mandatory run messages also live in the artifact.

    Every line is written and flushed as it is produced rather than at the end,
    so a run stopped by a hard kill still leaves a readable journal in its
    temporary directory. A step that raises deletes its directory, and its
    diagnostics survive only on the terminal.
    """

    def __init__(self, step_dir: Path, name: str = "log.txt"):
        self.path = step_dir / name
        self._fh = open(self.path, "a", encoding="utf-8")

    def say(self, msg: str) -> None:
        """Print the message and append it to the step's journal."""
        print(msg, flush=True)
        self._fh.write(msg + "\n")
        self._fh.flush()

    def close(self) -> None:
        self._fh.close()
