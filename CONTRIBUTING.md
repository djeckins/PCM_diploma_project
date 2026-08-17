# Working on this repository

This is the code behind a single-author MSc dissertation. It is published so that
the results can be inspected, checked and reused, not as a project looking for
features. That shapes what is useful here:

* Corrections are welcome — a defect in a measurement, a claim that the
  artefacts do not support, a step that does not reproduce. Open an issue and say
  which file and which number.
* Questions about method are welcome too. `docs/METHODS.md` is the place to
  start, and every number in `docs/` can be traced to a file under `runs/`.
* Large refactors are not the point. The layout is deliberate and the test
  suite encodes why (see the invariants below). A change that makes the code
  tidier but weakens one of those guarantees is not an improvement.

Adding a new protein needs no code change at all: it is a YAML file in
`configs/`, then `autodetect` and `events`. See "Applying it to a new protein" in
the README. If you find yourself editing `src/` to add a system, something has
gone wrong.

## Setting up

The portable route, which needs nothing but a Python 3.11 environment:

```bash
python -m pip install -e ".[dev]"
python -m pytest
```

The route the author uses, via conda:

```bash
conda env create -f environment.yml
conda activate chem
pip install -e . --no-deps
python -m pytest
```

`make setup` / `make test` / `make lint` wrap the second route, but the Makefile
assumes conda lives at `~/miniforge3`. If yours is elsewhere, either edit that
line or use the commands above directly.

`environment.yml` describes the analysis workstation. It pins nothing beyond
version floors; the exact versions the published results were produced with are
tabulated in `THIRD_PARTY_NOTICES.md`.

## Running the tests

```bash
python -m pytest
```

84 tests, about three seconds, and no simulation data required. The suite
either builds its own small synthetic inputs or reads the result tables that are
committed under `runs/`. Nothing in `tests/` opens a trajectory, so this works on
a clean checkout with no datasets mounted — you can point `PCM2_DATA_ROOT` at a
path that does not exist and everything still passes. The same commands run in
CI on every push (`.github/workflows/tests.yml`).

One trap worth knowing about. `src/pcm2/runtime.py` locates the project root
relative to the imported module, so if you have another checkout of `pcm2`
installed in editable mode, `import pcm2` and `python -m pcm2.run` may quietly
resolve to *that* tree, and a pipeline run will write into it. `pytest` is safe,
because `pyproject.toml` puts `src` first for it. Check with:

```bash
python -c "import pcm2; print(pcm2.__file__)"
```

If the path is not inside this checkout, run `pip install -e .` here, or prefix
commands with `PYTHONPATH=src`.

Two lint layers, matching CI:

```bash
python -m ruff check --select E9,F63,F7,F82 src tests tools notebooks   # must pass
python -m ruff check --statistics src tests tools                       # advisory
```

The first selection catches real defects — syntax errors, undefined names,
broken assertion strings — and the whole tree passes it. The second runs the
project's full configuration and currently reports 75 stylistic findings, the
three largest groups being `zip()` without an explicit `strict=` (19), import
ordering (18) and unused loop variables (12). They are advisory rather than
build-breaking; the count is meant only to go down.

## The invariants the test suite enforces

These are not style preferences. Each one is a defect class that was actually
hit during development, turned into a test so it cannot come back. If your change
makes one of these tests fail, the test is almost certainly right.

Layer boundaries, checked by grepping the sources
(`tests/test_layer_boundaries_by_grep.py`). The reasoning is that if the model
layer *can* open a trajectory, one day it will.

* MDAnalysis may be imported only in the reading layer (`io.py`,
  `autodetect.py`, `events.py`, `labels.py`, `coords.py`) and in
  `features/__init__.py`, which needs periodic-boundary distances.
* scikit-learn and xgboost must not appear anywhere in the feature modules or
  the reading layer. Features are measured physics and must not know a model
  exists.
* `report.py` and `viz.py` must contain no `.fit(` — reporting reads artefacts,
  it never trains.

No caching; every step replaces its own directory wholesale
(`tests/test_overwrite_and_artifacts.py`). There is no cache layer in this
project on purpose, because stale intermediate state caused two contamination
incidents that were hard to find. Concretely:

* a step deletes whatever it does not produce, so an artefact cannot outlive its
  purpose;
* a step that fails leaves the previous complete output untouched, and no
  partial temporary directory behind;
* `PROVENANCE.json` is written next to the output, recording the config, a hash
  of the source tree and the library versions.

If you add a step, write it through `step_output()` in `src/pcm2/runtime.py` and
you get all of this for free.

Tests must not write into `runs/` (`tests/conftest.py`). A session fixture
snapshots the modification time and size of every file under `runs/` and compares
them at the end of the session; any change fails the whole suite. `runs/` is
published evidence, so use `tmp_path` for anything that writes.

The numbers in `docs/` are generated, never typed
(`tests/test_docs_generated_blocks.py`). Anything between
`<!-- BEGIN GENERATED: <name> -->` and `<!-- END GENERATED -->` in `docs/` is
written by `tools/genblocks.py` from the artefacts under `runs/`. Only the
opening marker carries the block name, which is `results`, `acceptance` or
`schema`; the closing marker is bare, exactly as written here. The generator
matches that pair and nothing else, so repeating the name in the closing marker
makes the block invisible to it and its numbers stop being maintained.
To update:

```bash
python tools/genblocks.py            # rewrite the blocks from the artefacts
python tools/genblocks.py --check    # what the test and CI run
```

Do not hand-edit inside those markers — the check compares them against their
source and fails on any divergence. The test also verifies itself, by corrupting
a block in a scratch copy and confirming the corruption is caught. The same rule
applies to `runs/**/PROVENANCE.json`: those record what actually happened on a
particular machine, so editing one to look tidier turns a record into a fiction.

A typo in a config must crash, not fall back to a default
(`tests/test_config_unknown_key_must_crash.py`). Unknown keys, at any nesting
level, raise `ConfigError`. A silently ignored key is a result computed with
settings nobody chose.

The feature schema is computed, not counted
(`tests/test_schema_is_computed_not_counted.py`). Column counts are derived from
the number of bins and sites. Never write a literal column count into the code or
the docs — quote what the artefacts say, and state which count you mean, since
the schema, the modelled subset and the shipped feature table legitimately differ.

The vendored third-party files are checksummed
(`tests/test_vendored_checksums.py`). The two files in
`external/rao2019_heuristic/` are byte-for-byte copies from their authors, and
their sha256 sums are recomputed on every test run and cross-checked against
`PROVENANCE.md`. They are inputs to the published comparison, so a silent
substitution would quietly change what the benchmark means. Do not edit or
reformat them; see `THIRD_PARTY_NOTICES.md`.

## What must never be committed

* Trajectories and topologies — `.xtc`, `.trr`, `.tpr`, `.gro`, `.pdb`,
  `.top`, `.itp`, and the rest. They are tens of gigabytes, and mostly not the
  author's to redistribute. `.gitignore` blocks them by extension and a CI step
  fails the build if any appear.
* `runs/*/coords/coords.parquet` — superimposed Cα coordinates, which come
  close to being a re-encoding of the input trajectory rather than a derived
  measurement. It is regenerated locally by the `coords` step.
* Caches and editor droppings — already ignored; if you find yourself adding
  a `--force`, stop and read the ignore rule first.

Derived result tables *are* committed, on purpose: they are what makes the
repository useful to someone who does not have the trajectories.

## Style

Everything in English, including comments. Comments explain *why* a thing is done
where the reason is not obvious from the code, and name the defect class where
one motivated the design; they do not restate what the next line does.

Two wording rules apply to code comments, docs and any figure caption, because
they are matters of scientific honesty rather than taste:

* Explanations are associational, not causal. The mechanism readout says
  which axis looks unlike a conducting pore. Write "enriched among non-conducting
  frames", never "causes the block".
* Never quote a mechanism share, or a lift over chance, without the quantity it
  is measured against — the horizon τ, the base rate, and the number of positive
  events. Several systems here have single-digit event counts, and a ratio
  without its denominator is not a result.

## Licensing of contributions

By contributing you agree that your code is licensed under the MIT licence
(`LICENSE`), and any figures or derived tables under CC-BY-4.0
(`LICENSE-DATA.md`), consistent with the rest of the repository.
