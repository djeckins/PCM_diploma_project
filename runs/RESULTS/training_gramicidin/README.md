# Gramicidin A

Training protein. A narrow single-file channel; both salt/temperature conditions conduct, and all ten trajectory units are in the pool.

System `gramicidin`; complete crossings recorded: **112**.

## Own run (`in_system/`)

Trajectory-level splits; `x chance` = average precision over the base rate.

| arm | AP | base rate | ×chance | n / n+ |
|---|---|---|---|---|
| clock | 0.1682 | 0.1096 | 1.53 | 12385/1357 |
| linear | 0.2912 | 0.1096 | 2.66 | 12385/1357 |
| linear_control_structural | 0.1066 | 0.1096 | 0.97 | 12385/1357 |
| plsfma_coords | 0.1134 | 0.1096 | 1.04 | 12385/1357 |
| published_rao2019 | 0.1158 | 0.1096 | 1.06 | 12385/1357 |
| published_rao2019_win | 0.1164 | 0.1096 | 1.06 | 12385/1357 |
| trees | 0.2839 | 0.1096 | 2.59 | 12385/1357 |

Event-level monitor: crossings warned: 0.350 (structure-matched random null: 0.111); empty alarm episodes: 12 of 40

## Predicted without ever seeing these units (`zero_shot/`)

- τ = 1000 ps, leave-gramicidin-out: `298K_2Msalt/01` x1.49, `298K_2Msalt/02` x1.05, `298K_2Msalt/03` x1.83, `298K_2Msalt/04` x1.42, `298K_2Msalt/05` x1.33, `330K_1Msalt/01` x1.95, `330K_1Msalt/02` x1.27, `330K_1Msalt/03` x1.25, `330K_1Msalt/04` x2.05, `330K_1Msalt/05` x2.13
- τ = 2000 ps, leave-gramicidin-out: `298K_2Msalt/01` x1.51, `298K_2Msalt/02` x1.70, `298K_2Msalt/03` x1.37, `298K_2Msalt/04` x1.48, `298K_2Msalt/05` x1.02, `330K_1Msalt/01` x2.11, `330K_1Msalt/02` x1.34, `330K_1Msalt/03` x1.96, `330K_1Msalt/04` x2.04, `330K_1Msalt/05` x1.73
- τ = 4000 ps, leave-gramicidin-out: `298K_2Msalt/01` x1.42, `298K_2Msalt/02` x1.62, `298K_2Msalt/03` x1.44, `298K_2Msalt/04` x1.73, `298K_2Msalt/05` x1.49, `330K_1Msalt/01` x1.79, `330K_1Msalt/02` x1.12, `330K_1Msalt/03` x1.11, `330K_1Msalt/04` x2.97, `330K_1Msalt/05` x1.12

Full held-out unit reports with reality checks: `unit-report-gramicidin-298K_2Msalt-01/`, `unit-report-gramicidin-330K_1Msalt-02/`, `unit-report-gramicidin-330K_1Msalt-04/`.
