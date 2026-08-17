# KcsA E71A (K+)

Training protein: the conducting K+ unit of the KcsA filter-mutant family. Its sibling mutants stay out of training and serve as negative controls.

System `kcsa`, conditions E71A; complete crossings recorded: **49**.

## Own run (`in_system/`)

Trajectory-level splits; `x chance` = average precision over the base rate.

| arm | AP | base rate | ×chance | n / n+ |
|---|---|---|---|---|
| clock | 0.0071 | 0.0130 | 0.55 | 24982/324 |
| linear | 0.0073 | 0.0130 | 0.57 | 24982/324 |
| linear_control_structural | 0.0072 | 0.0130 | 0.55 | 24982/324 |
| plsfma_coords | 0.0102 | 0.0130 | 0.79 | 24982/324 |
| published_rao2019 | 0.0239 | 0.0130 | 1.84 | 24982/324 |
| published_rao2019_win | 0.0186 | 0.0130 | 1.43 | 24982/324 |
| trees | 0.0032 | 0.0010 | 3.14 | 18738/19 |

Event-level monitor: crossings warned: 1.000 (structure-matched random null: 1.000); empty alarm episodes: 0 of 1

## Predicted without ever seeing these units (`zero_shot/`)

- τ = 1000 ps, leave-kcsa_family-out: `E71A/1` x1.22
- τ = 2000 ps, leave-kcsa_family-out: `E71A/1` x1.21
- τ = 4000 ps, leave-kcsa_family-out: `E71A/1` x1.16

Full held-out unit reports with reality checks: `unit-report-kcsa-E71A-1/`.
