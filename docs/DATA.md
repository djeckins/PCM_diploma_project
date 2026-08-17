# DATA — what the pipeline reads, where it comes from, and what works without it

The molecular-dynamics trajectories analysed in this work are not part of
this repository. What is included is everything the pipeline computed *from*
them: event records, feature and label tables, model scores, figures and
reports under `runs/`. That split is deliberate — the trajectories are tens of
gigabytes and are not all the author's to redistribute, while the derived tables
are small enough to ship and are what most readers actually need.

## 1. The trajectory sets

Five configs cover four proteins, because the KcsA filter mutants were
simulated twice, once in K⁺ and once in Na⁺, and each ionic condition is its
own config; the table below has one row per config. Everything in it is read
from `configs/*.yaml`, including the thermostat temperatures. Four of the five
carry those with `origin: detected`: for MthK, gramicidin A and both KcsA arms
the basis recorded against `protocol.temperature_K` is "ref-t of the thermostat
read from the run tpr with gmx dump (GROMACS in the separate 'gmx' conda
environment)". Cx43 is the exception — its temperature is `origin: declared`,
on the basis of `ref_t` in the production `.mdp` shipped with its deposit.
GROMACS is not required to run any code in this repository; recovering those
four temperatures was its only role here, and it happened outside the pipeline
(section 2).

| system | config | contents | protocol as recorded | trajectory span, frames analysed |
|---|---|---|---|---|
| MthK | `mthk.yaml` | 2 replicas at +200 mV, 2 at −200 mV | CHARMM36, charge scaling 0.78, 323.15 K | 500 ns each; stride 5 → 250 ps between table frames, 2001 frames per replica |
| Gramicidin A | `gramicidin.yaml` | 5 replicas at 298.15 K / 2 M KCl, 5 at 330 K / 1 M KCl, 500 mV | charge scaling 0.70 (298 K) and 0.75 (330 K); force field not recorded | 200 ns each; stride 4 → 160 ps between frames, 1251 frames per replica |
| KcsA filter mutants, K⁺ | `kcsa.yaml` | E71A, G77A·E71A, T75A·E71A, one replica each | CHARMM36m, no charge scaling, 300 mV, 290 K | 1 µs (E71A, T75A·E71A) and 2 µs (G77A·E71A); stride 4 → 160 ps between frames |
| KcsA filter mutants, Na⁺ | `kcsa_na.yaml` | same three variants, one replica each, Na⁺ as permeant | CHARMM36m, no charge scaling, 300 mV, 290 K | same spans as the K⁺ arm |
| Connexin-43 | `cx43.yaml` | replicas R2 and R4 at 200 mV | CHARMM36, no charge scaling, 303.15 K | 100 ns each; stride 2 → 100 ps between frames, 1001 frames per replica |

Two facts about this table are worth stating plainly rather than leaving a
reader to infer them. First, replication is uneven: gramicidin A contributes ten
trajectory units and MthK four, while every KcsA arm and each Cx43 condition is
a single replica per condition. Second, the KcsA spans are much longer than the
others because conduction there is rare, not because the runs are better
sampled per unit time.

### Provenance and credit

The two sources below are public deposits, and both are licensed CC-BY-4.0,
which is what makes it legitimate to publish derived tables computed from them
as long as they are credited. Deposit titles, authorship and licence terms were
checked against the Zenodo records, and the two accompanying papers against
PubMed, on 2026-08-17.

KcsA filter mutants (K⁺ and Na⁺ arms). Downloaded from the Zenodo deposit
*Data from Selectivity filter mutations shift ion permeation mechanism in
potassium channels*, Mironenko, de Groot & Kopec, 2024,
[10.5281/zenodo.12623504](https://doi.org/10.5281/zenodo.12623504) (CC-BY-4.0;
archives `ExampleTrajectories.tar.gz`, about 40 GB, and
`InputFilesAndScripts.tar.gz`). The accompanying paper is A. Mironenko,
B. L. de Groot and W. Kopec, *Selectivity filter mutations shift ion permeation
mechanism in potassium channels*, PNAS Nexus 2024, 3(7), pgae272,
[10.1093/pnasnexus/pgae272](https://doi.org/10.1093/pnasnexus/pgae272). The
three variants in `configs/kcsa.yaml` — E71A and the G77A and T75A substitutions
on that background — are the mutants that paper studies, which is the reason
they enter this work as engineered non-conducting controls rather than as
training material.

Connexin-43. Downloaded from the Zenodo deposit *Molecular dynamics
simulation data 1: Structure of the connexin-43 gap junction channel in a
putative closed state*, Acosta-Gutierrez & Gervasio, 2023,
[10.5281/zenodo.8191584](https://doi.org/10.5281/zenodo.8191584) (CC-BY-4.0;
the deposit ships the production trajectories, topology and parameter files,
per-voltage archives and a single initial structure). The accompanying paper is
C. Qi, S. Acosta Gutierrez, P. Lavriha *et al.*, *Structure of the connexin-43
gap junction channel in a putative closed state*, eLife 2023, 12, RP87616,
[10.7554/eLife.87616](https://doi.org/10.7554/eLife.87616) (the record is a
reviewed preprint, so its own citation line gives the article number as RP87616
rather than e87616). The topology used here is a GROMACS `.top` assembled from
that deposit's files: the shipped `initial_structure` matches the trajectories
atom for atom, while the deposited `topol.top` water count predates the final
solvent trim and was corrected to the structure's composition. That reassembly
is recorded where it was carried out, in the comment above the Cx43 entry in
`tools/init_configs.py`. The `origins` block of `configs/cx43.yaml` records two
other Cx43 findings — the excluded byte-identical replica (next subsection) and
the shared assembly behind the R2 and R4 lineages — and holds no water count or
solvent trim of its own.

Gramicidin A and MthK. Both sets were provided by the project supervisor at
Queen Mary University of London and are used here with permission. Permission
to use is not permission to redistribute: they are not published here, not
linked, and not re-hosted from this repository, and no licence over them is
granted or implied. Enquiries about access should go to the supervising
laboratory, as `README.md` and `THIRD_PARTY_NOTICES.md` §4 also
state. Unlike the two deposits above, this entry rests on the project record
rather than on anything a reader can check here: the repository holds no
independent provenance record for these two sets — no deposit, no DOI, no
licence file of their own — so read it as the author's statement.

For gramicidin A the `*-transit-times.dat` files that came with the
trajectories are a pre-existing crossing annotation. They are used only as a
cross-check: `events.source: own` in the config, so every label in this work
derives from this project's own detector, and the events step prints the
comparison between the two sources. The `summary.json` written by that step in
`runs/gramicidin-stride4/events/` records `n_own` and `n_provided` side by side
for every replica. No annotation file is redistributed here.

### One defect found in a public deposit

The Cx43 deposit ships three production replicas, and the third is not an
independent run: `protein_center_200mV_R3.xtc` is a byte-identical copy of
`protein_center_200mV_R2.xtc` — the same file size of 3,340,449,480 bytes, the
same MD5, and the full-resolution detector finds identical crossings, ions 208
and 198 entering and leaving at identical times. Treating both as replicas would
double-count one trajectory and fake a lineage split, so R3 is excluded and two
genuine replicas remain, R2 and R4. The evidence is recorded in the `origins`
block of `configs/cx43.yaml`. Anyone obtaining the same deposit will see three
replicas and should know why only two are analysed here.

## 2. What the pipeline needs on disk

Per replica, two files:

- a topology — either a GROMACS `.tpr` or a GROMACS `.top`. Both are read
  the same way, and element identities are assigned from the exact topology
  masses, so a hand-assembled `.top` is a first-class input.
- a trajectory — an `.xtc`.

Optionally a third: a pre-existing crossing annotation, picked up as a
cross-check if present. Discovery looks for a file whose name contains `transit`
and ends in `.dat`, or, in the per-replica layout below, any file whose name
contains `events`.

Nothing else is required. There is no need for `.mdp`, `.edr`, `.cpt` or
checkpoint files, and no GROMACS installation is needed to run anything in this
repository — the reading layer is MDAnalysis, which parses the `.tpr` itself.
GROMACS was used once, outside the pipeline: `gmx dump` on each run's `.tpr`
supplied the thermostat reference temperatures that four of the five configs
now record as detected (section 1).

`tools/init_configs.py` discovers files by walking a data root, and recognises
two layouts (both implemented in `src/pcm2/inputs.py`):

Layout A — one directory per replica. Used by MthK, both KcsA arms and Cx43.

```
<data_root>/<condition>/<replica>/*.xtc
<data_root>/<condition>/<replica>/*.tpr        (or *.top)
```

Discovery looks for `*.tpr` and then `*.top` in the replica directory itself,
and nowhere else — there is no lookup in the condition directory above it. A
replica holding neither falls back to a single topology found anywhere under the
data root, and only if the whole tree holds exactly one candidate: with two or
more, `inputs.py` returns nothing rather than pair a trajectory with another
build's atom order, and the replica is skipped with a note carried into the
config's basis string. MthK relies on that fallback: replica `200MV/1` holds the
only topology anywhere under its data root, and the other three replicas take
that same `.tpr`.

Layout B — numbered files in one condition directory. Used by gramicidin A.

```
<data_root>/<condition>/NN-<anything>.xtc
<data_root>/<condition>/NN-<anything>.tpr     (or .top; the same stem as the xtc)
<data_root>/<condition>/NN-*transit*.dat      (optional annotation)
```

The replica id is the part of the file stem before the first hyphen, so `NN` is
what identifies a replica. The topology is looked for under the *whole* stem of
its trajectory, only with a topology suffix, and the same single-topology
fallback applies when no such file exists; the annotation needs to match the
`NN` prefix alone.

The layouts this project ran on, as recorded in the configs:

```
<PCM2_DATA_ROOT>/
  CHARMM-ECC-078/                  MthK
    200MV/1/{md_1.tpr, fixed_c.xtc}
    200MV/2/fixed_c.xtc            topology taken from 200MV/1
    neg-200MV/{1,2}/fixed_c.xtc
  Shared_KCl_500mV/                gramicidin A
    298K_2Msalt/{01..05}-pt7scaling-*.{tpr,xtc}, {01..05}-transit-times.dat
    330K_1Msalt/{01..05}-pt75scaling-*.{tpr,xtc}, {01..05}-transit-times.dat
  3 белок/kcsa_root/               KcsA filter mutants, K+
    {E71A,G77AE71A,T75AE71A}/1/{production.tpr, production.xtc}
  3 белок/kcsa_na_root/            KcsA filter mutants, Na+
    {E71A,G77AE71A,T75AE71A}/1/{production.tpr, production.xtc}
  protein_center_200mV/cx43_root/  connexin-43
    200mV/{R2,R4}/{topol.top, R2.xtc | R4.xtc}
```

`3 белок` is a local folder name, Russian for "protein 3". Directory names are
not meaningful to the code — only the config has to agree with what is on disk.
If you rename it, either edit the two `data_root` lines in
`tools/init_configs.py` before running discovery, or start from
`configs/templates/` and set the path by hand.

## 3. Running the pipeline against a local copy

```bash
# 1. environment
conda env create -f environment.yml && conda activate chem
pip install -e . --no-deps
python -m pytest                       # must be green before anything else

# 2. point at the data
export PCM2_DATA_ROOT=/absolute/path/to/your/data

# 3. configs: either discover them from the files on disk ...
python tools/init_configs.py           # writes configs/<system>.yaml
#    ... or fill in a template by hand (see configs/templates/README.md)

# 4. one system, all steps
python -m pcm2.run all --config configs/mthk.yaml
```

Step 3 is where the machine-specific part lives, and it is the only place. The
code carries no absolute paths at all. A *generated* config does: discovery
writes the real paths of the machine it walked. The configs shipped here are
generated ones with those paths taken back out — replaced by the placeholder
`/path/to/data` before publication (section 5) — and their portable copies in
`configs/templates/` carry `${PCM2_DATA_ROOT}` in the same position instead.

Neither placeholder is expanded for you, and this is the one thing worth
knowing before a first run. Nothing in `src/pcm2` calls `os.path.expandvars` or
reads the environment while loading a config, so `${PCM2_DATA_ROOT}` stays a
literal substring of a YAML string: a template used as it stands will load,
pass strict validation, and only then fail to find its trajectories.
`PCM2_DATA_ROOT` is read by `tools/init_configs.py` alone, to decide which
directory to walk. A shipped or template config therefore has to be regenerated
by discovery (step 3 above) or edited by hand before it will run;
`configs/templates/README.md` sets out both routes.

Running the steps one at a time, in order, is often more useful than `all`,
because the first two are the ones that must succeed before anything can be
measured:

```bash
python -m pcm2.run autodetect --config configs/<system>.yaml
python -m pcm2.run events     --config configs/<system>.yaml
python -m pcm2.run features   --config configs/<system>.yaml
python -m pcm2.run labels     --config configs/<system>.yaml
python -m pcm2.run train      --config configs/<system>.yaml
```

`autodetect` measures the atom selections, the pore axis and extent, the
coordination cutoff and the architecture profile, and writes them back into the
config with their bases — nothing can be measured before it has run. `events`
detects the crossings, which are needed both for the labels and for the delivery
descriptors. `features` writes the descriptor table, `labels` turns crossings
into per-frame targets at each horizon, and `train` fits and evaluates. The full
step order, including the ones not listed here, is in `README.md`; each step
replaces its own output directory wholesale, so a rerun cannot mix old and new
state.

Two practical notes. `coords`, the step that writes superimposed Cα
coordinates, needs the trajectories, and its table is deliberately not shipped:
`runs/*/coords/` keeps only the step's provenance record and log, while
`coords.parquet` is git-ignored, because it is closer to a re-encoding of the
input trajectory than to a derived measurement. Anything downstream of it
— `tools/plsfma_selfcheck.py`, and the `plsfma_coords` comparison arm inside
`train` — is therefore only reproducible with the trajectories in hand; without
them the arm records itself as skipped rather than failing. And a new protein
needs no labels at all: to score an unseen system with the frozen pooled model,
`autodetect` and `events` are enough (see "Applying it to a new protein" in
`README.md`).

## 4. What works with no trajectories at all

The repository is deliberately usable without any of the data above. Verified on
the shipped tree, with no data root present:

- The test suite passes. 84 tests, a few seconds, no trajectories and no
  data root. This is the quickest way for a reader to confirm the code in front
  of them behaves as described.
- Every run's results are readable as shipped. `runs/<system>-stride<N>/`
  holds the events, features and labels tables, the trained models, the figures
  and the report for all five configs, and `runs/RESULTS/<role>/` arranges the
  same material by the role each system plays. `runs/FINAL/` holds the write-up
  figures and tables. `docs/RESULTS.md` is the short version.
- The interactive notebook runs in table mode. Because the feature and label
  tables ship, `notebooks/permeation_predict.ipynb` can query any frame of any
  of the five systems, get a verdict with its mechanism breakdown, and draw a
  timeline, all from the shipped tables. Only *window mode* — measuring
  descriptors for a time range that no run covers, or for a new protein — needs
  the trajectories.

What is not possible without the data: re-measuring descriptors, adding a
protein, and the `coords` step with everything downstream of it.

## 5. A note on the paths inside the shipped artefacts

The absolute dataset paths in `configs/*.yaml`, in the `runs/**/PROVENANCE.json`
records and in a few step logs were replaced with a placeholder before
publication, so that the repository does not carry the local directory layout of
the machine the runs were made on. They now read `/path/to/data` where a dataset
path was recorded and `/path/to/pcm2` where a path inside the project tree was.
Those strings are therefore redacted rather than literal: the directory
*structure* below the data root is exactly what was used, but the part above it
is not. The portable form of every config is in `configs/templates/`.

The redaction reached one place beyond the path fields, which is worth stating
so that a reader does not take the `origins` block to be verbatim throughout. An
`origins` entry whose basis quotes the directory discovery walked was rewritten
too. Each config carries two such entries — `data.data_root` and
`data.conditions` both record the discovery string — so ten basis strings in
all. In `configs/mthk.yaml` both now read "discovered on disk under
/path/to/data/CHARMM-ECC-078: 200MV: 2 replicas; neg-200MV: 2 replicas", and the
`configs/templates/` copies carry `${PCM2_DATA_ROOT}` in the same position. No
other basis string quotes a path, and none was altered: what each entry says
about *why* a value is what it is — the origin, the replica counts, the measured
quantities and their margins — is exactly as it was written. Everything else in
the provenance records is as recorded too: code state, library versions,
platform, thread counts and timestamps.

## 6. Third-party data that *is* inside this repository

One directory ships third-party files: `external/rao2019_heuristic/`, two data
tables from the CHAP project, redistributed unmodified with checksums, used for
the published-baseline comparison arm and (in the case of the hydrophobicity
scale) by the feature layer itself. Its licence (MIT), its checksums and its
citations are in `external/rao2019_heuristic/PROVENANCE.md` and in
`THIRD_PARTY_NOTICES.md` at the repository root. The licence terms of the
repository's own code and of the derived data under `runs/` are in `LICENSE` and
`LICENSE-DATA.md`.
