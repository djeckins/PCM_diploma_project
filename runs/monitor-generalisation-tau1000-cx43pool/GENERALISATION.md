# Conducting-channels generalisation (tau = 1000 ps)

Training pool: gramicidin + MthK + KcsA-E71A + Cx43 (conducting channels only).

## Leave-one-protein-out (median AP over the held protein's own chance)

| held-out protein | full columns | physics-only |
|---|---|---|
| gramicidin | x1.305 | x1.1 |
| mthk | x1.35 | x1.395 |
| kcsa_family | x1.25 | x1.19 |
| cx43 | x1.0 | x1.7 |

## Negative control: monitor applied to unseen mutant arms

| unit | frames | base rate | READY fraction | mean P(ready) | AP/chance |
|---|---|---|---|---|---|
| kcsa:G77AE71A/1 | 12501 | 0.0015 | 0.000 | 0.008 | 0.75 |
| kcsa:T75AE71A/1 | 6251 | 0.0000 | 0.000 | 0.010 | undefined |
| kcsa_na:E71A/1 | 6251 | 0.0000 | 0.000 | 0.033 | undefined |
| kcsa_na:G77AE71A/1 | 12501 | 0.0034 | 0.000 | 0.016 | 2.37 |
| kcsa_na:T75AE71A/1 | 6251 | 0.0010 | 0.000 | 0.010 | 3.75 |

Columns dropped in physics-only: 38 (listed in results.json).