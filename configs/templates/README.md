# Config templates

The five files here are the portable form of the five system configs. Three
forms of the same config exist, differing only in what stands where the data
root belongs:

* `configs/*.yaml`, the configs the published runs were made with, carry the
  redacted placeholder `/path/to/data`;
* `configs/templates/*.yaml`, the files here, carry `${PCM2_DATA_ROOT}` in the
  same positions;
* a *generated* config — as `tools/init_configs.py` writes it, or as either
  form becomes once a real root is substituted in — carries real absolute
  paths. No config of that form is shipped.

```yaml
data:
  data_root: ${PCM2_DATA_ROOT}/CHARMM-ECC-078
  conditions:
  - id: 200MV
    replicas:
    - id: '1'
      topology: ${PCM2_DATA_ROOT}/CHARMM-ECC-078/200MV/1/md_1.tpr
      trajectory: ${PCM2_DATA_ROOT}/CHARMM-ECC-078/200MV/1/fixed_c.xtc
```

The placeholder stands in every path field — `data.data_root`, and each
replica's `topology`, `trajectory` and, where a crossing annotation exists,
`events` — and in the two `origins` bases that quote the directory discovery
walked, those of `data.data_root` and `data.conditions`, which is ten basis
strings across the five files. Nothing else differs: system identity,
architecture profile, conditions, replicas, per-replica protocols, horizons,
detected geometry and the rest of the `origins` block are carried over
unchanged. Substituting `/path/to/data` for `${PCM2_DATA_ROOT}` in a template
reproduces the corresponding `configs/*.yaml` byte for byte, and the reverse
substitution reproduces the template. Neither of the two shipped forms is
runnable, though: only a real path makes a config that can find its
trajectories.

## Neither shipped form is expanded automatically

`${PCM2_DATA_ROOT}` is a plain string in a YAML file, and so is
`/path/to/data`. Nothing in `src/pcm2` expands environment variables inside a
config or reads the environment while loading one, so a config in either
shipped form will load, pass strict validation and only then fail to find its
trajectories. Both therefore have to be regenerated from the files on disk
(route 1) or edited by hand (route 2) before anything will run; the
substitution is a step you perform.

`PCM2_DATA_ROOT` *is* read by `tools/init_configs.py`, which is where the
variable belongs: it tells the discovery walk which directory to look in.

## Route 1 — regenerate from the files on disk (preferred)

If you have the trajectories arranged as `docs/DATA.md` describes, do not fill
in a template at all. Let discovery do it:

```bash
export PCM2_DATA_ROOT=/absolute/path/to/your/data
python tools/init_configs.py
```

This writes fresh `configs/<system>.yaml` with real paths, then validates each
one strictly and reports the condition and replica counts it found. Detected
keys are filled in later, by `python -m pcm2.run autodetect --config
configs/<system>.yaml`, which measures them from your own copy of the
trajectories and rewrites the config in place with the measurement recorded as
its basis.

One caveat: `tools/init_configs.py` hard-codes the directory names it expects
under the data root, including `3 белок` for the two KcsA archives (a local
folder name, Russian for "protein 3"). Either reproduce those directory names,
or edit the `data_root` entries near the top of that script, or use route 2.

## Route 2 — instantiate a template

Useful when your directory layout differs from the one discovery expects, or
when you want the exact configuration these results were produced with rather
than a freshly discovered one.

```bash
export PCM2_DATA_ROOT=/absolute/path/to/your/data
for f in configs/templates/*.yaml; do
    sed "s|\${PCM2_DATA_ROOT}|$PCM2_DATA_ROOT|g" "$f" > "configs/$(basename "$f")"
done
```

The same thing without `sed`, if a path might contain characters `sed` treats
specially:

```bash
python - <<'EOF'
import os, pathlib
root = os.environ["PCM2_DATA_ROOT"]
for t in sorted(pathlib.Path("configs/templates").glob("*.yaml")):
    out = pathlib.Path("configs") / t.name
    out.write_text(t.read_text().replace("${PCM2_DATA_ROOT}", root))
    print(out)
EOF
```

Write the result to `configs/` at the top level, not to a subdirectory: that is
where both the pipeline and the notebook look for a system by its bare id, and
it is where `tools/init_configs.py` writes too. Overwriting the shipped configs
costs nothing, since these templates remain the record and the substitution is
reversible.

Then check that the paths resolve, and only afterwards run anything:

```bash
python -c "
from pcm2 import config
cfg = config.load('configs/mthk.yaml')
from pathlib import Path
for cond, rep in cfg.replicas():
    for k in ('topology', 'trajectory'):
        p = Path(rep[k])
        print(('ok   ' if p.exists() else 'MISSING'), k, p)
"
```

## Filling a template in by hand

Only the file paths need editing — `data.data_root` and, for every replica,
`topology` and `trajectory` (and `events`, where a crossing annotation exists;
`null` otherwise). Adjust `id` and `lineage` if your replica directories are
named differently; `lineage` marks which replicas come from one independent
assembly and controls how folds are split, so replicas that share an
equilibrated ancestor should share a lineage.

The `origins` block records where each key's value came from — `detected`,
`declared` or `default` — with the basis in prose. Two of those bases, for
`data.data_root` and `data.conditions`, quote the directory discovery walked,
so the placeholder appears in prose there as well; being prose, it is harmless
left as it is. If you change a value, change its basis to say why.

Leave the rest alone unless you know why you are changing it. In particular:

- Detected keys — `data.selections.*`, `data.pore.*`,
  `features.coordination_cutoff_A`, `system.arch_profile.*` — are measurements
  of *these* trajectories. They are kept in the templates so the shipped
  results stay reproducible on the same data. Re-run `autodetect` and they will
  be re-measured from your copy.
- Declared keys — the horizons `labels.tau_ps` and `labels.primary_tau_ps`,
  the crossing cylinder, the feature window, the model seed and the frozen
  hyperparameter grids — are choices with recorded reasons, and changing one
  changes what the results mean rather than only how they are computed. The
  horizon in particular is part of the question being asked, not a tuning knob.

## For a different protein

Do not start from a template of a different system: the detected keys describe
another channel's geometry, and copying them silently substitutes one protein's
pore for another's. Assemble a config from the files on disk instead
(`tools/init_configs.py`, or a minimal hand-written config with
`system.id`, `data.data_root`, the conditions with their replicas and lineages,
and the horizon keys), then run `autodetect`, which fills the rest in and
refuses rather than guesses where it cannot measure something.

## Without any trajectories

Nothing here needs filling in. The configs shipped in `configs/` carry
placeholder paths and are still enough to read every result in `runs/` and to
run `notebooks/permeation_predict.ipynb` in table mode, which uses the feature
and label tables that ship with the repository. Real paths are needed only to
measure something new. See `docs/DATA.md`.
