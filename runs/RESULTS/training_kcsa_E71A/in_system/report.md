# Run report kcsa (stride 4)

## Per-arm metrics (pooled, primary horizon)

| arm | AP | base rate | ×chance | n / n+ |
|---|---|---|---|---|
| clock | 0.0071 | 0.0130 | 0.55 | 24982/324 |
| linear | 0.0073 | 0.0130 | 0.57 | 24982/324 |
| linear_control_structural | 0.0072 | 0.0130 | 0.55 | 24982/324 |
| plsfma_coords | 0.0102 | 0.0130 | 0.79 | 24982/324 |
| published_rao2019 | 0.0239 | 0.0130 | 1.84 | 24982/324 |
| published_rao2019_win | 0.0186 | 0.0130 | 1.43 | 24982/324 |
| trees | 0.0032 | 0.0010 | 3.14 | 18738/19 |

## Paired comparison with the head arm (AP difference, paired bootstrap over trajectories)
- vs clock: Δ=0.0173, CI [0.002468810636285308, 0.050275887473596674], groups 2 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs linear: Δ=0.0196, CI [0.002547532024589789, 0.050539899886604664], groups 2 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs linear_control_structural: Δ=0.0179, CI [0.0023468872531576587, 0.05059600562886337], groups 2 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs plsfma_coords: Δ=0.0179, CI [0.0018105329312007867, 0.05020178692709098], groups 2 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs published_rao2019: Δ=0.0174, CI [0.0023115003960634067, 0.04558970389730101], groups 2 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden
- vs published_rao2019_win: Δ=0.0189, CI [0.00231663043355031, 0.04908401370567248], groups 2 — with few groups the interval coverage is not guaranteed; the interval is shown, not hidden

## Monitor (event-level headline quantity)
- crossings warned: 1.000 (structure-matched random null: 1.000)
- empty alarm episodes: 0 of 1

## Anchor ablation
- transit: AP=0.1188 at base rate 0.0040 (×29.58)

## Acceptance (computable gates)
- edge_exception2_recorded: PASS
- no_ceiling_hits: PASS
- penalty_in_df_units: PASS
- calibration_recorded_per_fold: PASS
- events_source_declared: PASS
- vendored_checksums: PASS
- figure_catalog: PASS
- answers_with_mechanism: PASS
- resting_frames_defined: FAIL
- applicability_mask_written: PASS