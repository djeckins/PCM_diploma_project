# Conducting-channels generalisation (tau = 4000 ps)

Training pool: gramicidin + MthK + KcsA-E71A (conducting channels only).

## Leave-one-protein-out (median AP over the held protein's own chance)

| held-out protein | full columns | physics-only |
|---|---|---|
| gramicidin | x1.4649999999999999 | x1.43 |
| mthk | x1.06 | x1.0550000000000002 |
| kcsa_family | x1.16 | x1.14 |

## Negative control: monitor applied to unseen mutant arms

| unit | frames | base rate | READY fraction | mean P(ready) | AP/chance |
|---|---|---|---|---|---|
| cx43:200mV/R2 | 1001 | 0.0749 | 0.002 | 0.154 | 3.91 |
| cx43:200mV/R4 | 1001 | 0.0000 | 0.000 | 0.131 | undefined |
| kcsa:G77AE71A/1 | 12501 | 0.0060 | 0.000 | 0.080 | 1.21 |
| kcsa:T75AE71A/1 | 6251 | 0.0000 | 0.000 | 0.061 | undefined |
| kcsa_na:E71A/1 | 6251 | 0.0000 | 0.000 | 0.103 | undefined |
| kcsa_na:G77AE71A/1 | 12501 | 0.0140 | 0.000 | 0.139 | 3.3 |
| kcsa_na:T75AE71A/1 | 6251 | 0.0040 | 0.000 | 0.153 | 0.66 |

Columns dropped in physics-only: 21 (listed in results.json).