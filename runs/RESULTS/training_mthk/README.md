# MthK

Training protein. A wide-vestibule K+ channel at both field signs; all four units are in the pool.

System `mthk`; complete crossings recorded: **607**.

## Own run (`in_system/`)

Trajectory-level splits; `x chance` = average precision over the base rate.

| arm | AP | base rate | ×chance | n / n+ |
|---|---|---|---|---|
| clock | 0.5062 | 0.5008 | 1.01 | 7980/3996 |
| linear | 0.5790 | 0.5008 | 1.16 | 7980/3996 |
| linear_control_structural | 0.4489 | 0.5008 | 0.90 | 7980/3996 |
| plsfma_coords | 0.4742 | 0.5008 | 0.95 | 7980/3996 |
| published_rao2019 | 0.4977 | 0.5008 | 0.99 | 7980/3996 |
| published_rao2019_win | 0.5104 | 0.5008 | 1.02 | 7980/3996 |
| trees | 0.6845 | 0.5008 | 1.37 | 7980/3996 |

Event-level monitor: crossings warned: 0.068 (structure-matched random null: 0.039); empty alarm episodes: 5 of 26

## Predicted without ever seeing these units (`zero_shot/`)

- τ = 1000 ps, leave-mthk-out: `200MV/1` x1.56, `200MV/2` x1.99, `neg-200MV/1` x1.21, `neg-200MV/2` x1.18
- τ = 2000 ps, leave-mthk-out: `200MV/1` x1.18, `200MV/2` x1.36, `neg-200MV/1` x1.08, `neg-200MV/2` x1.10
- τ = 4000 ps, leave-mthk-out: `200MV/1` x1.07, `200MV/2` x1.14, `neg-200MV/1` x1.03, `neg-200MV/2` x1.05

Full held-out unit reports with reality checks: `unit-report-mthk-200MV-1/`, `unit-report-mthk-200MV-2/`, `unit-report-mthk-neg-200MV-1/`.
