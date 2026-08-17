# Licence for figures, tables, reports and derived data

Applies to: the contents of `runs/` and `docs/` — figures (`*.png`), derived tables
(`*.parquet`, `*.csv`), result and provenance records (`*.json`), and the generated
reports (`*.md`, `*.txt`).

Licence: Creative Commons Attribution 4.0 International (CC-BY-4.0).
<https://creativecommons.org/licenses/by/4.0/>
Legal code: <https://creativecommons.org/licenses/by/4.0/legalcode>

Copyright (c) 2026 Angelina Ivanova

You are free to share and adapt this material for any purpose, including
commercially, provided you give appropriate credit, link to the licence, and
indicate whether changes were made.

## Why a separate licence from the code

These artefacts are not software, and a software licence fits them poorly. CC-BY-4.0 is
the standard for scientific figures and datasets, is accepted by Zenodo and by journal
data-availability policies, and — unlike the MIT licence — expressly licenses *sui
generis* database rights, which exist for structured tables in the EU and in Russia even
where copyright does not subsist in the underlying numbers.

## How to cite

Attribution should name the author, the dissertation or accompanying publication, and
this repository. A `CITATION.cff` file, and a DOI obtained by archiving a release on
Zenodo, are the practical way to make this unambiguous.

## What this licence does not cover

* The code that produced these artefacts — see `LICENSE` (MIT).
* `external/rao2019_heuristic/*.json` — third-party data under their own MIT
  licence; see `external/rao2019_heuristic/LICENSE.CHAP.md` and `THIRD_PARTY_NOTICES.md`.
* The input molecular-dynamics trajectories — not contained in this repository, not
  redistributed by it, and not licensed here. See the "Where the data came from"
  section of `README.md`.

## A note on derived quantities

The tables and figures in `runs/` are measurements *computed from* the input
trajectories: pore-radius profiles, coordination counts, event labels, model scores and
calibration curves. They are this project's own output and are licensed above.

One artefact is deliberately excluded from the repository for this reason:
`runs/*/coords/coords.parquet` holds superimposed Cα coordinates, which are close to a
re-encoding of the input trajectory rather than a derived measurement. It is excluded by
the `runs/*/coords/*` rule in `.gitignore` — which re-includes only that step's
`PROVENANCE.json` and `log.txt` — and is regenerated locally by the `coords` step. Keep it
excluded: it is the one output whose redistribution would amount to redistributing
third-party trajectory data.
