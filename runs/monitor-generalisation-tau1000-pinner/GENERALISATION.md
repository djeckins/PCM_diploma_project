# Conducting-channels generalisation (tau = 1000 ps)

Training pool: gramicidin + MthK + KcsA-E71A (conducting channels only).

## Leave-one-protein-out (median AP over the held protein's own chance)

| held-out protein | full columns | physics-only |
|---|---|---|
| gramicidin | x1.775 | x1.785 |
| mthk | x1.23 | x1.05 |
| kcsa_family | x1.18 | x1.2 |

## Negative control: monitor applied to unseen mutant arms

| unit | frames | base rate | READY fraction | mean P(ready) | AP/chance |
|---|---|---|---|---|---|
| cx43:200mV/R2 | 1001 | 0.0202 | 0.000 | 0.112 | 1.19 |
| cx43:200mV/R4 | 1001 | 0.0000 | 0.000 | 0.112 | undefined |
| kcsa:G77AE71A/1 | 12501 | 0.0015 | 0.000 | 0.136 | 0.97 |
| kcsa:T75AE71A/1 | 6251 | 0.0000 | 0.000 | 0.135 | undefined |
| kcsa_na:E71A/1 | 6251 | 0.0000 | 0.000 | 0.115 | undefined |
| kcsa_na:G77AE71A/1 | 12501 | 0.0034 | 0.000 | 0.119 | 6.65 |
| kcsa_na:T75AE71A/1 | 6251 | 0.0010 | 0.000 | 0.124 | 1.78 |

Columns dropped in physics-only: 21 (listed in results.json).