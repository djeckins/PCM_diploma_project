# Connexin-43 (zero-shot)

Zero-shot test on an unrelated architecture: a wide aqueous gap-junction pore, double membrane, no selectivity filter; the most distant transfer in this study.

System `cx43`; complete crossings recorded: **2**.

## Own run (`in_system/`)

Trajectory-level splits; `x chance` = average precision over the base rate.

| arm | AP | base rate | ×chance | n / n+ |
|---|---|---|---|---|
| clock | — | — | — | not defined |
| linear | — | — | — | not defined |
| linear_control_structural | — | — | — | not defined |
| plsfma_coords | — | — | — | not defined |
| published_rao2019 | 0.0291 | 0.0375 | 0.78 | 1922/72 |
| published_rao2019_win | 0.0243 | 0.0375 | 0.65 | 1922/72 |

## Predicted without ever seeing these units (`zero_shot/`)

- τ = 1000 ps, negative control: `200mV/R2` READY 0.02 (mean P 0.070), `200mV/R4` READY 0.03 (mean P 0.082)
- τ = 2000 ps, negative control: `200mV/R2` READY 0.00 (mean P 0.100), `200mV/R4` READY 0.00 (mean P 0.100)
- τ = 4000 ps, negative control: `200mV/R2` READY 0.00 (mean P 0.154), `200mV/R4` READY 0.00 (mean P 0.131)

Full held-out unit reports with reality checks: `unit-report-cx43-200mV-R2/`, `unit-report-cx43-200mV-R4/`.
