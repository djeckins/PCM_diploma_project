# Conducting-channels generalisation (tau = 4000 ps)

Training pool: gramicidin + MthK + KcsA-E71A + Cx43 (conducting channels only).

## Leave-one-protein-out (median AP over the held protein's own chance)

| held-out protein | full columns | physics-only |
|---|---|---|
| gramicidin | x1.3900000000000001 | x1.375 |
| mthk | x1.065 | x1.0550000000000002 |
| kcsa_family | x1.18 | x1.21 |
| cx43 | x2.76 | x4.8 |

## Negative control: monitor applied to unseen mutant arms

| unit | frames | base rate | READY fraction | mean P(ready) | AP/chance |
|---|---|---|---|---|---|
| kcsa:G77AE71A/1 | 12501 | 0.0060 | 0.000 | 0.072 | 1.61 |
| kcsa:T75AE71A/1 | 6251 | 0.0000 | 0.000 | 0.058 | undefined |
| kcsa_na:E71A/1 | 6251 | 0.0000 | 0.000 | 0.063 | undefined |
| kcsa_na:G77AE71A/1 | 12501 | 0.0140 | 0.000 | 0.085 | 4.42 |
| kcsa_na:T75AE71A/1 | 6251 | 0.0040 | 0.000 | 0.123 | 0.75 |

Columns dropped in physics-only: 38 (listed in results.json).