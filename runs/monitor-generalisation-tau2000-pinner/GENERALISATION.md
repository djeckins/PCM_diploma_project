# Conducting-channels generalisation (tau = 2000 ps)

Training pool: gramicidin + MthK + KcsA-E71A (conducting channels only).

## Leave-one-protein-out (median AP over the held protein's own chance)

| held-out protein | full columns | physics-only |
|---|---|---|
| gramicidin | x1.695 | x1.6 |
| mthk | x1.09 | x1.04 |
| kcsa_family | x1.2 | x1.24 |

## Negative control: monitor applied to unseen mutant arms

| unit | frames | base rate | READY fraction | mean P(ready) | AP/chance |
|---|---|---|---|---|---|
| cx43:200mV/R2 | 1001 | 0.0408 | 0.239 | 0.168 | 1.52 |
| cx43:200mV/R4 | 1001 | 0.0000 | 0.235 | 0.166 | undefined |
| kcsa:G77AE71A/1 | 12501 | 0.0030 | 0.877 | 0.262 | 0.87 |
| kcsa:T75AE71A/1 | 6251 | 0.0000 | 0.999 | 0.277 | undefined |
| kcsa_na:E71A/1 | 6251 | 0.0000 | 0.057 | 0.163 | undefined |
| kcsa_na:G77AE71A/1 | 12501 | 0.0069 | 0.096 | 0.170 | 2.11 |
| kcsa_na:T75AE71A/1 | 6251 | 0.0021 | 0.419 | 0.199 | 1.76 |

Columns dropped in physics-only: 21 (listed in results.json).