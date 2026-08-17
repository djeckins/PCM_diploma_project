# Pooled leave-one-channel-out (tau = 1000 ps)

Trained on the pooled tables of all systems minus one held-out KcsA filter mutant; evaluated by ranking within the held-out unit against its own base rate.

| held-out unit | base rate | trees | linear | clock |
|---|---|---|---|---|
| kcsa:E71A/1 | 0.049 | x1.14 | x1.52 | x1.16 |
| kcsa:G77AE71A/1 | 0.002 | x6.64 | x2.35 | x11.86 |
| kcsa:T75AE71A/1 | 0.000 | — | — | — |
| kcsa_na:E71A/1 | 0.000 | — | — | — |
| kcsa_na:G77AE71A/1 | 0.003 | x10.88 | x6.14 | x1.79 |
| kcsa_na:T75AE71A/1 | 0.001 | x0.74 | x0.61 | x1.63 |

Potential system-label columns (disjoint per-system ranges): geo_lining_n, geo_r_bin0_A, geo_r_bin1_A, geo_r_bin2_A, geo_bin5_nsearch, hyd_water_pore_n, hyd_lining_n_facing, hyd_water_bin0_n, hyd_water_bin1_n, hyd_water_bin2_n, hyd_water_bin3_n, hyd_water_bin6_n, ww_n_waters, ww_max_gap_A, sym_com_z_spread_A.

coordinate-based FMA cannot pool across proteins: Cartesian spaces are incommensurable.