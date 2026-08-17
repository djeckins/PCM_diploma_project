# Vendored: Rao et al. (2019) hydrophobic-gate heuristic — the authors' own files

The published criterion has no closed-form expression: the method is a tabulated
energy surface plus a hydropathy scale, and the decision boundary is the
2.6 kJ/mol (≈1 RT) contour extracted from the surface. Retyping the table by hand
would substitute an approximation for the published method, so the files are
vendored byte-for-byte.

Downloaded: 2026-08-14, directly from the authors' repository (raw.githubusercontent.com).
Re-verified against upstream `master` on 2026-08-17: both files are still byte-for-byte
identical to the source URLs below (sha256 recomputed on freshly downloaded copies).

| file | sha256 | source |
|---|---|---|
| `heuristic_grid.json` | `bb29c21e1339e71e7a877d5095a4d72c4ac3ab06ce04c863fbcfb407e567c4f0` | <https://github.com/channotation/chap/blob/master/scripts/heuristic/heuristic_grid.json> |
| `wimley_white_1996.json` | `d898884281bf666074ea8c2a5629896b31b9bea6bed4d4389a736ef355542b5c` | <https://github.com/channotation/chap/blob/master/share/data/hydrophobicity/wimley_white_1996.json> |

The files were not edited. `tests/test_vendored_checksums.py` recomputes the sums,
so a substitution fails the test suite.

## Scope of the criterion

The criterion predicts local dewetting of a hydrophobic constriction, i.e. that the
water free-energy barrier at the narrowest point exceeds ~1 RT, and it was
trained and validated against that quantity computed in simulation. Conduction is
one of several consequences of dewetting.

## Sources

* Rao S., Klesse G., Stansfeld P.J., Sansom M.S.P., Tucker S.J. (2019).
  *A heuristic derived from analysis of the ion channel structural proteome permits
  the rapid identification of hydrophobic gates.* PNAS 116(28):13989–13995.
  <https://doi.org/10.1073/pnas.1902702116> — the criterion itself, the 2.6 kJ/mol
  contour, the Wimley–White axis, the 0.45 nm hydrophobicity band.
* Klesse G., Rao S., Sansom M.S.P., Tucker S.J. (2019). *CHAP: A Versatile Tool for
  the Structural and Functional Annotation of Ion Channel Pores.* JMB 431(17):3353–3365.
  <https://doi.org/10.1016/j.jmb.2019.06.003> — definitions of the pore radius and of
  pore-facing residues, on which the heuristic's inputs rest.
* Authors' implementation: <https://github.com/channotation/chap/tree/master/scripts/heuristic>
  (`heuristic_score.py`, `heuristic_prediction.R`), <https://www.channotation.org/docs/heuristic_method/>.

## License

CHAP is distributed under the MIT License, copyright (c) 2016–2018 Gianni Klesse,
Shanlin Rao, Mark S. P. Sansom and Stephen J. Tucker. The upstream licence text is
reproduced verbatim next to the files it covers, in `LICENSE.CHAP.md`
(sha256 `c64029774ae7ecfc589ef358cec0a1b16b8384528d4b7febeb997e38aac7c605`), as the
licence requires: "The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software."

The MIT terms permit the unmodified redistribution done here. Vendored are two data
files from that repository, unmodified. They are third-party *inputs*: no `pcm2` source
file is derived from them, and they remain under their own licence rather than this
repository's.

Both files are load-bearing, and not only for the comparison arm:

* `heuristic_grid.json` is read by the Rao-2019 baseline (`pcm2/baselines/published.py`);
* `wimley_white_1996.json` is read by the feature layer
  (`pcm2/features/_common.py:load_ww_scale`, consumed in `pcm2/features/__init__.py`),
  so it is required by the main `features` pipeline step, not merely by the benchmark.

Removing either file therefore breaks the pipeline itself, which is the reason they are
vendored rather than fetched on demand.

Note on the hydrophobicity scale: `wimley_white_1996.json` is CHAP's rescaled
([-1, 1]) form of the whole-residue interfacial hydrophobicity scale of
Wimley W.C. & White S.H. (1996), *Experimentally determined hydrophobicity scale for
proteins at membrane interfaces*, Nature Structural Biology 3(10):842–848,
<https://doi.org/10.1038/nsb1096-842>. The file itself points at
<http://blanco.biomol.uci.edu/hydrophobicity_scales.html>. The primary measurement is
theirs; the normalisation and the JSON encoding are CHAP's.

## Declared deviations from the paper

1. Per-frame application. The authors' criterion is defined on time averages
   (an averaged profile); here it is applied both per frame and in a window-averaged
   version. The run reports the difference between the two versions as a measure of
   how legitimate the transfer is.
2. The input radius is taken from our own profile estimator, the same quantity
   CHAP uses (inscribed sphere), in our implementation. The estimator discrepancy
   is not measured against the CHAP binary (absent from the environment); this
   degrades the comparison and is recorded as such.
