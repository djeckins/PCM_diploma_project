# G77A-E71A with Na+ (zero-shot)

Zero-shot test within a known family: the rare Na+ leak of the G77A double mutant. The monitor never saw this unit, its protein, or any Na+ system in training.

System `kcsa_na`, conditions G77AE71A; complete crossings recorded: **7**.

## Own run (`in_system/`)

Trajectory-level splits; `x chance` = average precision over the base rate.

| arm | AP | base rate | ×chance | n / n+ |
|---|---|---|---|---|
| clock | 0.0011 | 0.0020 | 0.58 | 24982/49 |
| linear | 0.0012 | 0.0020 | 0.62 | 24982/49 |
| linear_control_structural | 0.0018 | 0.0020 | 0.91 | 24982/49 |
| plsfma_coords | 0.0016 | 0.0020 | 0.80 | 24982/49 |
| published_rao2019 | 0.0015 | 0.0020 | 0.74 | 24982/49 |
| published_rao2019_win | 0.0014 | 0.0020 | 0.69 | 24982/49 |
| trees | — | — | — | not defined |

## Predicted without ever seeing these units (`zero_shot/`)

- τ = 1000 ps, negative control: `G77AE71A/1` READY 0.00 (mean P 0.009)
- τ = 2000 ps, negative control: `G77AE71A/1` READY 0.00 (mean P 0.049)
- τ = 4000 ps, negative control: `G77AE71A/1` READY 0.00 (mean P 0.139)

Full held-out unit reports with reality checks: `unit-report-kcsa_na-G77AE71A-1/`.
