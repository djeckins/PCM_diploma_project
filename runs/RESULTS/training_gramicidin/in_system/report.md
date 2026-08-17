# Run report gramicidin (stride 4)

## Per-arm metrics (pooled, primary horizon)

| arm | AP | base rate | ×chance | n / n+ |
|---|---|---|---|---|
| clock | 0.1682 | 0.1096 | 1.53 | 12385/1357 |
| linear | 0.2912 | 0.1096 | 2.66 | 12385/1357 |
| linear_control_structural | 0.1066 | 0.1096 | 0.97 | 12385/1357 |
| plsfma_coords | 0.1134 | 0.1096 | 1.04 | 12385/1357 |
| published_rao2019 | 0.1158 | 0.1096 | 1.06 | 12385/1357 |
| published_rao2019_win | 0.1164 | 0.1096 | 1.06 | 12385/1357 |
| trees | 0.2839 | 0.1096 | 2.59 | 12385/1357 |

## Paired comparison with the head arm (AP difference, paired bootstrap over trajectories)
- vs clock: Δ=0.1042, CI [0.041776573527923466, 0.15816391119187354], groups 10 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs linear: Δ=-0.0083, CI [-0.04192727778644835, 0.021809247895696736], groups 10 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs linear_control_structural: Δ=0.1791, CI [0.14286977134094303, 0.22949426864625252], groups 10 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs plsfma_coords: Δ=0.1720, CI [0.12925754661709135, 0.22433467513164396], groups 10 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs published_rao2019: Δ=0.1717, CI [0.13838549938060393, 0.21725321718291407], groups 10 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs published_rao2019_win: Δ=0.1713, CI [0.13418934257582793, 0.22486923157154406], groups 10 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden

## Monitor (event-level headline quantity)
- crossings warned: 0.350 (structure-matched random null: 0.111)
- empty alarm episodes: 12 of 40

## Anchor ablation
- resting: AP=0.1366 at base rate 0.0453 (×3.02)
- transit: AP=0.7536 at base rate 0.7016 (×1.07)

## Acceptance (computable gates)
- edge_exception2_recorded: PASS
- no_ceiling_hits: PASS
- penalty_in_df_units: PASS
- calibration_recorded_per_fold: PASS
- events_source_declared: PASS
- vendored_checksums: PASS
- figure_catalog: PASS
- answers_with_mechanism: PASS
- resting_frames_defined: PASS
- applicability_mask_written: PASS