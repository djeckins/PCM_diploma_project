# Run report mthk (stride 5)

## Per-arm metrics (pooled, primary horizon)

| arm | AP | base rate | ×chance | n / n+ |
|---|---|---|---|---|
| clock | 0.5062 | 0.5008 | 1.01 | 7980/3996 |
| linear | 0.5790 | 0.5008 | 1.16 | 7980/3996 |
| linear_control_structural | 0.4489 | 0.5008 | 0.90 | 7980/3996 |
| plsfma_coords | 0.4742 | 0.5008 | 0.95 | 7980/3996 |
| published_rao2019 | 0.4977 | 0.5008 | 0.99 | 7980/3996 |
| published_rao2019_win | 0.5104 | 0.5008 | 1.02 | 7980/3996 |
| trees | 0.6845 | 0.5008 | 1.37 | 7980/3996 |

## Paired comparison with the head arm (AP difference, paired bootstrap over trajectories)
- vs clock: Δ=0.1723, CI [0.13467567098073285, 0.21685906548210243], groups 4 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs linear: Δ=0.1031, CI [0.06771335421372704, 0.127212751993873], groups 4 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs linear_control_structural: Δ=0.2156, CI [0.15485039899579023, 0.2598588875082932], groups 4 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs plsfma_coords: Δ=0.2030, CI [0.16481972951544777, 0.2539320902805108], groups 4 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs published_rao2019: Δ=0.1887, CI [0.16120028600855796, 0.2301248179255002], groups 4 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs published_rao2019_win: Δ=0.1761, CI [0.1215192315876733, 0.22015350805661765], groups 4 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden

## Monitor (event-level headline quantity)
- crossings warned: 0.068 (structure-matched random null: 0.039)
- empty alarm episodes: 5 of 26

## Anchor ablation
- resting: AP=0.3882 at base rate 0.2373 (×1.64)
- transit: AP=0.7636 at base rate 0.6037 (×1.26)

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