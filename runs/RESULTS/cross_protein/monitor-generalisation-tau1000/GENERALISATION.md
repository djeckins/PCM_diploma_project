# Conducting-channels generalisation (tau = 1000 ps)

Training pool: gramicidin + MthK + KcsA-E71A (conducting channels only).

## Leave-one-protein-out (median AP over the held protein's own chance)

| held-out protein | full columns | physics-only |
|---|---|---|
| gramicidin | x1.455 | x2.0700000000000003 |
| mthk | x1.385 | x1.3450000000000002 |
| kcsa_family | x1.22 | x1.18 |

## Negative control: monitor applied to unseen mutant arms

| unit | frames | base rate | READY fraction | mean P(ready) | AP/chance |
|---|---|---|---|---|---|
| cx43:200mV/R2 | 1001 | 0.0202 | 0.018 | 0.070 | 1.17 |
| cx43:200mV/R4 | 1001 | 0.0000 | 0.028 | 0.082 | undefined |
| kcsa:G77AE71A/1 | 12501 | 0.0015 | 0.000 | 0.006 | 1.14 |
| kcsa:T75AE71A/1 | 6251 | 0.0000 | 0.000 | 0.012 | undefined |
| kcsa_na:E71A/1 | 6251 | 0.0000 | 0.000 | 0.008 | undefined |
| kcsa_na:G77AE71A/1 | 12501 | 0.0034 | 0.000 | 0.009 | 17.66 |
| kcsa_na:T75AE71A/1 | 6251 | 0.0010 | 0.000 | 0.012 | 0.65 |

Columns dropped in physics-only: 21 (listed in results.json).