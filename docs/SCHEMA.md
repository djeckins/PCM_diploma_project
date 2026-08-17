# SCHEMA — column dictionary

The schema is computed from the config; the table below is generated from the
run's `features/schema.json` and `features/applicability.json` and therefore
cannot diverge from the code. The column count is derived, not literal.

<!-- BEGIN GENERATED: schema -->
### cx43-stride2 — 104 columns (count derived)

| column | block | units | missing means | indicator | verdict |
|---|---|---|---|---|---|
| `geo_r_constriction_A` | geometry | A | no search slice inside the pore | geo_n_search | active |
| `geo_z_constriction_A` | geometry | A | no search slice inside the pore | geo_n_search | active |
| `geo_n_search` | geometry | count | never missing: the counter is always defined | — | active |
| `geo_search_frac` | geometry | frac | never missing | — | active |
| `geo_boundary_frac` | geometry | frac | never missing | — | active |
| `geo_lining_carbonyl_cos` | geometry | cos | no lining carbonyls in the band | geo_lining_n | active |
| `geo_lining_n` | geometry | count | never missing | — | active |
| `geo_r_bin0_A` | geometry | A | no search slices in the bin | geo_bin0_nsearch | active |
| `geo_r_bin1_A` | geometry | A | no search slices in the bin | geo_bin1_nsearch | active |
| `geo_r_bin2_A` | geometry | A | no search slices in the bin | geo_bin2_nsearch | active |
| `geo_r_bin3_A` | geometry | A | no search slices in the bin | geo_bin3_nsearch | active |
| `geo_r_bin4_A` | geometry | A | no search slices in the bin | geo_bin4_nsearch | active |
| `geo_r_bin5_A` | geometry | A | no search slices in the bin | geo_bin5_nsearch | active |
| `geo_r_bin6_A` | geometry | A | no search slices in the bin | geo_bin6_nsearch | active |
| `geo_r_bin7_A` | geometry | A | no search slices in the bin | geo_bin7_nsearch | active |
| `geo_bin0_nsearch` | geometry | count | never missing | — | active |
| `geo_bin1_nsearch` | geometry | count | never missing | — | active |
| `geo_bin2_nsearch` | geometry | count | never missing | — | active |
| `geo_bin3_nsearch` | geometry | count | never missing | — | active |
| `geo_bin4_nsearch` | geometry | count | never missing | — | active |
| `geo_bin5_nsearch` | geometry | count | never missing | — | active |
| `geo_bin6_nsearch` | geometry | count | never missing | — | active |
| `geo_bin7_nsearch` | geometry | count | never missing | — | active |
| `geo_r_bin0_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin0_nsearch | active |
| `geo_r_bin1_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin1_nsearch | active |
| `geo_r_bin2_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin2_nsearch | active |
| `geo_r_bin3_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin3_nsearch | active |
| `geo_r_bin4_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin4_nsearch | active |
| `geo_r_bin5_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin5_nsearch | active |
| `geo_r_bin6_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin6_nsearch | active |
| `geo_r_bin7_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin7_nsearch | active |
| `hyd_water_pore_n` | hydration | count | never missing | — | active |
| `hyd_min_bin_water_n` | hydration | count | never missing | — | active |
| `hyd_wet_frac_win` | hydration | frac | never missing | — | constant |
| `hyd_wet_logodds_win` | hydration | logit | never missing | — | active |
| `hyd_lining_hydrophobicity` | hydration | WW[-1,1] | no pore-facing residues in the band | hyd_lining_n_facing | active |
| `hyd_lining_n_facing` | hydration | count | never missing | — | active |
| `hyd_water_pore_density` | hydration | A^-3 | no searched profile volume in the bin range | geo_n_search | active |
| `hyd_lining_facing_per_sub` | hydration | count/subunit | subunit count unknown | — | active |
| `hyd_water_bin0_n` | hydration | count | never missing | — | active |
| `hyd_water_bin1_n` | hydration | count | never missing | — | active |
| `hyd_water_bin2_n` | hydration | count | never missing | — | active |
| `hyd_water_bin3_n` | hydration | count | never missing | — | active |
| `hyd_water_bin4_n` | hydration | count | never missing | — | active |
| `hyd_water_bin5_n` | hydration | count | never missing | — | active |
| `hyd_water_bin6_n` | hydration | count | never missing | — | active |
| `hyd_water_bin7_n` | hydration | count | never missing | — | active |
| `hyd_water_bin0_density` | hydration | A^-3 | no search slices in the bin | geo_bin0_nsearch | active |
| `hyd_water_bin1_density` | hydration | A^-3 | no search slices in the bin | geo_bin1_nsearch | active |
| `hyd_water_bin2_density` | hydration | A^-3 | no search slices in the bin | geo_bin2_nsearch | active |
| `hyd_water_bin3_density` | hydration | A^-3 | no search slices in the bin | geo_bin3_nsearch | active |
| `hyd_water_bin4_density` | hydration | A^-3 | no search slices in the bin | geo_bin4_nsearch | active |
| `hyd_water_bin5_density` | hydration | A^-3 | no search slices in the bin | geo_bin5_nsearch | active |
| `hyd_water_bin6_density` | hydration | A^-3 | no search slices in the bin | geo_bin6_nsearch | active |
| `hyd_water_bin7_density` | hydration | A^-3 | no search slices in the bin | geo_bin7_nsearch | active |
| `ww_n_waters` | water_wire | count | never missing | — | active |
| `ww_max_gap_A` | water_wire | A | fewer than two waters in the pore | ww_n_waters | active |
| `ww_continuous` | water_wire | 0/1 | fewer than two waters in the pore | ww_n_waters | constant |
| `ww_dipole_cos_mean` | water_wire | cos | no water in the pore | ww_n_waters | active |
| `ww_dipole_cos_std` | water_wire | cos | fewer than two waters in the pore | ww_n_waters | active |
| `ww_water_linear_density` | water_wire | A^-1 | never missing | — | active |
| `ww_max_gap_rel` | water_wire | frac | fewer than two waters in the pore | ww_n_waters | active |
| `occ_n_ions_pore` | occupancy | count | never missing | — | active |
| `occ_has_innermost` | occupancy | 0/1 | never missing | — | constant |
| `occ_innermost_z_A` | occupancy | A | no innermost ion (no ion in pore, or constriction not measured) | occ_has_innermost | active |
| `occ_innermost_dz_constr_A` | occupancy | A | no innermost ion | occ_has_innermost | active |
| `occ_innermost_coord_n` | occupancy | count | no innermost ion | occ_has_innermost | active |
| `occ_desolvation` | occupancy | frac | no innermost ion, or no ions outside the pore for the bulk reference | occ_has_innermost | active |
| `ns_filter_n_ions` | named_sites | count | no filter (structural mask) | — | inapplicable_structural |
| `ns_S0_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | inapplicable_structural |
| `ns_S1_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | inapplicable_structural |
| `ns_S2_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | inapplicable_structural |
| `ns_S3_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | inapplicable_structural |
| `ns_S4_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | inapplicable_structural |
| `ele_pot_protein_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_water_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_ions_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_lipid_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | constant |
| `ele_axial_charge_asym_e` | electrostatics | e | constriction not measured | geo_n_search | active |
| `sym_com_radial_spread_A` | symmetry | A | fewer than two subunits (structural mask) | — | active |
| `sym_com_z_spread_A` | symmetry | A | fewer than two subunits (structural mask) | — | active |
| `sym_nn_dist_cv` | symmetry | frac | fewer than three subunits (structural mask) | — | active |
| `sym_com_z_spread_rel` | symmetry | frac | fewer than two subunits (structural mask) | — | active |
| `flu_water_pore_var_win` | fluctuations | count^2 | never missing | — | active |
| `flu_ion_pore_var_win` | fluctuations | count^2 | never missing | — | active |
| `flu_r_constr_var_win` | fluctuations | A^2 | no measured radius values in the window | geo_n_search | active |
| `dyn_ion_persist` | dynamics | 0/1 | never missing | — | active |
| `dyn_ion_drift_A_win` | dynamics | A | ion not in the pore now, or was not at the window start | dyn_ion_persist | active |
| `dyn_wet_state_age_ps` | dynamics | ps | never missing | — | active |
| `dyn_occ_change_age_ps` | dynamics | ps | never missing | — | active |
| `dlv_n_upper_mouth` | delivery | count | never missing | — | active |
| `dlv_n_lower_mouth` | delivery | count | never missing | — | active |
| `dlv_n_entry_mouth` | delivery | count | never missing | — | active |
| `dlv_entry_wide_n` | delivery | count | never missing | — | active |
| `dlv_dist_entry_A` | delivery | A | no permeant in the wide entry zone | dlv_entry_wide_n | active |
| `dlv_approach_pair_n` | delivery | count | never missing | — | active |
| `dlv_approach_A_win` | delivery | A | distance not measured at one of the window endpoints | dlv_approach_pair_n | active |
| `dlv_t_since_cross_ps` | delivery | ps | before the trajectory's first crossing the measured object does not exist | dlv_crossed_before | active |
| `dlv_crossed_before` | delivery | 0/1 | never missing | — | active |
| `dlv_rate_win` | delivery | 1/ns | never missing | — | active |
| `bl_rao_E_kJmol` | baseline | kJ/mol | inputs not measured on this frame | geo_n_search | active |
| `bl_rao_E_win_kJmol` | baseline | kJ/mol | no measured inputs in the window | geo_n_search | active |
| `bl_rao_sigma_d` | baseline | contour distance sum | no pore-facing residue points measured on this frame | bl_rao_flag_n | active |
| `bl_rao_flag_n` | baseline | count | no pore-facing residue points measured on this frame | — | active |

### gramicidin-stride4 — 104 columns (count derived)

| column | block | units | missing means | indicator | verdict |
|---|---|---|---|---|---|
| `geo_r_constriction_A` | geometry | A | no search slice inside the pore | geo_n_search | active |
| `geo_z_constriction_A` | geometry | A | no search slice inside the pore | geo_n_search | active |
| `geo_n_search` | geometry | count | never missing: the counter is always defined | — | active |
| `geo_search_frac` | geometry | frac | never missing | — | active |
| `geo_boundary_frac` | geometry | frac | never missing | — | active |
| `geo_lining_carbonyl_cos` | geometry | cos | no lining carbonyls in the band | geo_lining_n | active |
| `geo_lining_n` | geometry | count | never missing | — | active |
| `geo_r_bin0_A` | geometry | A | no search slices in the bin | geo_bin0_nsearch | active |
| `geo_r_bin1_A` | geometry | A | no search slices in the bin | geo_bin1_nsearch | active |
| `geo_r_bin2_A` | geometry | A | no search slices in the bin | geo_bin2_nsearch | active |
| `geo_r_bin3_A` | geometry | A | no search slices in the bin | geo_bin3_nsearch | active |
| `geo_r_bin4_A` | geometry | A | no search slices in the bin | geo_bin4_nsearch | active |
| `geo_r_bin5_A` | geometry | A | no search slices in the bin | geo_bin5_nsearch | active |
| `geo_r_bin6_A` | geometry | A | no search slices in the bin | geo_bin6_nsearch | active |
| `geo_r_bin7_A` | geometry | A | no search slices in the bin | geo_bin7_nsearch | active |
| `geo_bin0_nsearch` | geometry | count | never missing | — | active |
| `geo_bin1_nsearch` | geometry | count | never missing | — | active |
| `geo_bin2_nsearch` | geometry | count | never missing | — | active |
| `geo_bin3_nsearch` | geometry | count | never missing | — | active |
| `geo_bin4_nsearch` | geometry | count | never missing | — | active |
| `geo_bin5_nsearch` | geometry | count | never missing | — | active |
| `geo_bin6_nsearch` | geometry | count | never missing | — | active |
| `geo_bin7_nsearch` | geometry | count | never missing | — | active |
| `geo_r_bin0_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin0_nsearch | active |
| `geo_r_bin1_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin1_nsearch | active |
| `geo_r_bin2_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin2_nsearch | active |
| `geo_r_bin3_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin3_nsearch | active |
| `geo_r_bin4_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin4_nsearch | active |
| `geo_r_bin5_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin5_nsearch | active |
| `geo_r_bin6_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin6_nsearch | active |
| `geo_r_bin7_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin7_nsearch | active |
| `hyd_water_pore_n` | hydration | count | never missing | — | active |
| `hyd_min_bin_water_n` | hydration | count | never missing | — | active |
| `hyd_wet_frac_win` | hydration | frac | never missing | — | active |
| `hyd_wet_logodds_win` | hydration | logit | never missing | — | active |
| `hyd_lining_hydrophobicity` | hydration | WW[-1,1] | no pore-facing residues in the band | hyd_lining_n_facing | active |
| `hyd_lining_n_facing` | hydration | count | never missing | — | active |
| `hyd_water_pore_density` | hydration | A^-3 | no searched profile volume in the bin range | geo_n_search | active |
| `hyd_lining_facing_per_sub` | hydration | count/subunit | subunit count unknown | — | active |
| `hyd_water_bin0_n` | hydration | count | never missing | — | active |
| `hyd_water_bin1_n` | hydration | count | never missing | — | active |
| `hyd_water_bin2_n` | hydration | count | never missing | — | active |
| `hyd_water_bin3_n` | hydration | count | never missing | — | active |
| `hyd_water_bin4_n` | hydration | count | never missing | — | active |
| `hyd_water_bin5_n` | hydration | count | never missing | — | active |
| `hyd_water_bin6_n` | hydration | count | never missing | — | active |
| `hyd_water_bin7_n` | hydration | count | never missing | — | active |
| `hyd_water_bin0_density` | hydration | A^-3 | no search slices in the bin | geo_bin0_nsearch | active |
| `hyd_water_bin1_density` | hydration | A^-3 | no search slices in the bin | geo_bin1_nsearch | active |
| `hyd_water_bin2_density` | hydration | A^-3 | no search slices in the bin | geo_bin2_nsearch | active |
| `hyd_water_bin3_density` | hydration | A^-3 | no search slices in the bin | geo_bin3_nsearch | active |
| `hyd_water_bin4_density` | hydration | A^-3 | no search slices in the bin | geo_bin4_nsearch | active |
| `hyd_water_bin5_density` | hydration | A^-3 | no search slices in the bin | geo_bin5_nsearch | active |
| `hyd_water_bin6_density` | hydration | A^-3 | no search slices in the bin | geo_bin6_nsearch | active |
| `hyd_water_bin7_density` | hydration | A^-3 | no search slices in the bin | geo_bin7_nsearch | active |
| `ww_n_waters` | water_wire | count | never missing | — | active |
| `ww_max_gap_A` | water_wire | A | fewer than two waters in the pore | ww_n_waters | active |
| `ww_continuous` | water_wire | 0/1 | fewer than two waters in the pore | ww_n_waters | active |
| `ww_dipole_cos_mean` | water_wire | cos | no water in the pore | ww_n_waters | active |
| `ww_dipole_cos_std` | water_wire | cos | fewer than two waters in the pore | ww_n_waters | active |
| `ww_water_linear_density` | water_wire | A^-1 | never missing | — | active |
| `ww_max_gap_rel` | water_wire | frac | fewer than two waters in the pore | ww_n_waters | active |
| `occ_n_ions_pore` | occupancy | count | never missing | — | active |
| `occ_has_innermost` | occupancy | 0/1 | never missing | — | active |
| `occ_innermost_z_A` | occupancy | A | no innermost ion (no ion in pore, or constriction not measured) | occ_has_innermost | active |
| `occ_innermost_dz_constr_A` | occupancy | A | no innermost ion | occ_has_innermost | active |
| `occ_innermost_coord_n` | occupancy | count | no innermost ion | occ_has_innermost | active |
| `occ_desolvation` | occupancy | frac | no innermost ion, or no ions outside the pore for the bulk reference | occ_has_innermost | active |
| `ns_filter_n_ions` | named_sites | count | no filter (structural mask) | — | inapplicable_structural |
| `ns_S0_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | inapplicable_structural |
| `ns_S1_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | inapplicable_structural |
| `ns_S2_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | inapplicable_structural |
| `ns_S3_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | inapplicable_structural |
| `ns_S4_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | inapplicable_structural |
| `ele_pot_protein_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_water_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_ions_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_lipid_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_axial_charge_asym_e` | electrostatics | e | constriction not measured | geo_n_search | active |
| `sym_com_radial_spread_A` | symmetry | A | fewer than two subunits (structural mask) | — | active |
| `sym_com_z_spread_A` | symmetry | A | fewer than two subunits (structural mask) | — | active |
| `sym_nn_dist_cv` | symmetry | frac | fewer than three subunits (structural mask) | — | inapplicable_structural |
| `sym_com_z_spread_rel` | symmetry | frac | fewer than two subunits (structural mask) | — | active |
| `flu_water_pore_var_win` | fluctuations | count^2 | never missing | — | active |
| `flu_ion_pore_var_win` | fluctuations | count^2 | never missing | — | active |
| `flu_r_constr_var_win` | fluctuations | A^2 | no measured radius values in the window | geo_n_search | active |
| `dyn_ion_persist` | dynamics | 0/1 | never missing | — | active |
| `dyn_ion_drift_A_win` | dynamics | A | ion not in the pore now, or was not at the window start | dyn_ion_persist | active |
| `dyn_wet_state_age_ps` | dynamics | ps | never missing | — | active |
| `dyn_occ_change_age_ps` | dynamics | ps | never missing | — | active |
| `dlv_n_upper_mouth` | delivery | count | never missing | — | active |
| `dlv_n_lower_mouth` | delivery | count | never missing | — | active |
| `dlv_n_entry_mouth` | delivery | count | never missing | — | active |
| `dlv_entry_wide_n` | delivery | count | never missing | — | active |
| `dlv_dist_entry_A` | delivery | A | no permeant in the wide entry zone | dlv_entry_wide_n | active |
| `dlv_approach_pair_n` | delivery | count | never missing | — | active |
| `dlv_approach_A_win` | delivery | A | distance not measured at one of the window endpoints | dlv_approach_pair_n | active |
| `dlv_t_since_cross_ps` | delivery | ps | before the trajectory's first crossing the measured object does not exist | dlv_crossed_before | active |
| `dlv_crossed_before` | delivery | 0/1 | never missing | — | active |
| `dlv_rate_win` | delivery | 1/ns | never missing | — | active |
| `bl_rao_E_kJmol` | baseline | kJ/mol | inputs not measured on this frame | geo_n_search | active |
| `bl_rao_E_win_kJmol` | baseline | kJ/mol | no measured inputs in the window | geo_n_search | active |
| `bl_rao_sigma_d` | baseline | contour distance sum | no pore-facing residue points measured on this frame | bl_rao_flag_n | active |
| `bl_rao_flag_n` | baseline | count | no pore-facing residue points measured on this frame | — | active |

### kcsa-stride4 — 104 columns (count derived)

| column | block | units | missing means | indicator | verdict |
|---|---|---|---|---|---|
| `geo_r_constriction_A` | geometry | A | no search slice inside the pore | geo_n_search | active |
| `geo_z_constriction_A` | geometry | A | no search slice inside the pore | geo_n_search | active |
| `geo_n_search` | geometry | count | never missing: the counter is always defined | — | active |
| `geo_search_frac` | geometry | frac | never missing | — | active |
| `geo_boundary_frac` | geometry | frac | never missing | — | active |
| `geo_lining_carbonyl_cos` | geometry | cos | no lining carbonyls in the band | geo_lining_n | active |
| `geo_lining_n` | geometry | count | never missing | — | active |
| `geo_r_bin0_A` | geometry | A | no search slices in the bin | geo_bin0_nsearch | active |
| `geo_r_bin1_A` | geometry | A | no search slices in the bin | geo_bin1_nsearch | active |
| `geo_r_bin2_A` | geometry | A | no search slices in the bin | geo_bin2_nsearch | active |
| `geo_r_bin3_A` | geometry | A | no search slices in the bin | geo_bin3_nsearch | active |
| `geo_r_bin4_A` | geometry | A | no search slices in the bin | geo_bin4_nsearch | active |
| `geo_r_bin5_A` | geometry | A | no search slices in the bin | geo_bin5_nsearch | active |
| `geo_r_bin6_A` | geometry | A | no search slices in the bin | geo_bin6_nsearch | active |
| `geo_r_bin7_A` | geometry | A | no search slices in the bin | geo_bin7_nsearch | active |
| `geo_bin0_nsearch` | geometry | count | never missing | — | active |
| `geo_bin1_nsearch` | geometry | count | never missing | — | active |
| `geo_bin2_nsearch` | geometry | count | never missing | — | active |
| `geo_bin3_nsearch` | geometry | count | never missing | — | active |
| `geo_bin4_nsearch` | geometry | count | never missing | — | active |
| `geo_bin5_nsearch` | geometry | count | never missing | — | active |
| `geo_bin6_nsearch` | geometry | count | never missing | — | active |
| `geo_bin7_nsearch` | geometry | count | never missing | — | active |
| `geo_r_bin0_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin0_nsearch | active |
| `geo_r_bin1_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin1_nsearch | active |
| `geo_r_bin2_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin2_nsearch | active |
| `geo_r_bin3_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin3_nsearch | active |
| `geo_r_bin4_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin4_nsearch | active |
| `geo_r_bin5_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin5_nsearch | active |
| `geo_r_bin6_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin6_nsearch | active |
| `geo_r_bin7_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin7_nsearch | active |
| `hyd_water_pore_n` | hydration | count | never missing | — | active |
| `hyd_min_bin_water_n` | hydration | count | never missing | — | active |
| `hyd_wet_frac_win` | hydration | frac | never missing | — | active |
| `hyd_wet_logodds_win` | hydration | logit | never missing | — | active |
| `hyd_lining_hydrophobicity` | hydration | WW[-1,1] | no pore-facing residues in the band | hyd_lining_n_facing | active |
| `hyd_lining_n_facing` | hydration | count | never missing | — | active |
| `hyd_water_pore_density` | hydration | A^-3 | no searched profile volume in the bin range | geo_n_search | active |
| `hyd_lining_facing_per_sub` | hydration | count/subunit | subunit count unknown | — | active |
| `hyd_water_bin0_n` | hydration | count | never missing | — | active |
| `hyd_water_bin1_n` | hydration | count | never missing | — | active |
| `hyd_water_bin2_n` | hydration | count | never missing | — | active |
| `hyd_water_bin3_n` | hydration | count | never missing | — | active |
| `hyd_water_bin4_n` | hydration | count | never missing | — | active |
| `hyd_water_bin5_n` | hydration | count | never missing | — | active |
| `hyd_water_bin6_n` | hydration | count | never missing | — | active |
| `hyd_water_bin7_n` | hydration | count | never missing | — | active |
| `hyd_water_bin0_density` | hydration | A^-3 | no search slices in the bin | geo_bin0_nsearch | active |
| `hyd_water_bin1_density` | hydration | A^-3 | no search slices in the bin | geo_bin1_nsearch | active |
| `hyd_water_bin2_density` | hydration | A^-3 | no search slices in the bin | geo_bin2_nsearch | active |
| `hyd_water_bin3_density` | hydration | A^-3 | no search slices in the bin | geo_bin3_nsearch | active |
| `hyd_water_bin4_density` | hydration | A^-3 | no search slices in the bin | geo_bin4_nsearch | active |
| `hyd_water_bin5_density` | hydration | A^-3 | no search slices in the bin | geo_bin5_nsearch | active |
| `hyd_water_bin6_density` | hydration | A^-3 | no search slices in the bin | geo_bin6_nsearch | active |
| `hyd_water_bin7_density` | hydration | A^-3 | no search slices in the bin | geo_bin7_nsearch | active |
| `ww_n_waters` | water_wire | count | never missing | — | active |
| `ww_max_gap_A` | water_wire | A | fewer than two waters in the pore | ww_n_waters | active |
| `ww_continuous` | water_wire | 0/1 | fewer than two waters in the pore | ww_n_waters | active |
| `ww_dipole_cos_mean` | water_wire | cos | no water in the pore | ww_n_waters | active |
| `ww_dipole_cos_std` | water_wire | cos | fewer than two waters in the pore | ww_n_waters | active |
| `ww_water_linear_density` | water_wire | A^-1 | never missing | — | active |
| `ww_max_gap_rel` | water_wire | frac | fewer than two waters in the pore | ww_n_waters | active |
| `occ_n_ions_pore` | occupancy | count | never missing | — | active |
| `occ_has_innermost` | occupancy | 0/1 | never missing | — | constant |
| `occ_innermost_z_A` | occupancy | A | no innermost ion (no ion in pore, or constriction not measured) | occ_has_innermost | active |
| `occ_innermost_dz_constr_A` | occupancy | A | no innermost ion | occ_has_innermost | active |
| `occ_innermost_coord_n` | occupancy | count | no innermost ion | occ_has_innermost | active |
| `occ_desolvation` | occupancy | frac | no innermost ion, or no ions outside the pore for the bulk reference | occ_has_innermost | active |
| `ns_filter_n_ions` | named_sites | count | no filter (structural mask) | — | active |
| `ns_S0_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ns_S1_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ns_S2_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ns_S3_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ns_S4_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ele_pot_protein_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_water_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_ions_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_lipid_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_axial_charge_asym_e` | electrostatics | e | constriction not measured | geo_n_search | active |
| `sym_com_radial_spread_A` | symmetry | A | fewer than two subunits (structural mask) | — | active |
| `sym_com_z_spread_A` | symmetry | A | fewer than two subunits (structural mask) | — | active |
| `sym_nn_dist_cv` | symmetry | frac | fewer than three subunits (structural mask) | — | active |
| `sym_com_z_spread_rel` | symmetry | frac | fewer than two subunits (structural mask) | — | active |
| `flu_water_pore_var_win` | fluctuations | count^2 | never missing | — | active |
| `flu_ion_pore_var_win` | fluctuations | count^2 | never missing | — | active |
| `flu_r_constr_var_win` | fluctuations | A^2 | no measured radius values in the window | geo_n_search | active |
| `dyn_ion_persist` | dynamics | 0/1 | never missing | — | active |
| `dyn_ion_drift_A_win` | dynamics | A | ion not in the pore now, or was not at the window start | dyn_ion_persist | active |
| `dyn_wet_state_age_ps` | dynamics | ps | never missing | — | active |
| `dyn_occ_change_age_ps` | dynamics | ps | never missing | — | active |
| `dlv_n_upper_mouth` | delivery | count | never missing | — | active |
| `dlv_n_lower_mouth` | delivery | count | never missing | — | active |
| `dlv_n_entry_mouth` | delivery | count | never missing | — | active |
| `dlv_entry_wide_n` | delivery | count | never missing | — | active |
| `dlv_dist_entry_A` | delivery | A | no permeant in the wide entry zone | dlv_entry_wide_n | active |
| `dlv_approach_pair_n` | delivery | count | never missing | — | active |
| `dlv_approach_A_win` | delivery | A | distance not measured at one of the window endpoints | dlv_approach_pair_n | active |
| `dlv_t_since_cross_ps` | delivery | ps | before the trajectory's first crossing the measured object does not exist | dlv_crossed_before | active |
| `dlv_crossed_before` | delivery | 0/1 | never missing | — | active |
| `dlv_rate_win` | delivery | 1/ns | never missing | — | active |
| `bl_rao_E_kJmol` | baseline | kJ/mol | inputs not measured on this frame | geo_n_search | active |
| `bl_rao_E_win_kJmol` | baseline | kJ/mol | no measured inputs in the window | geo_n_search | active |
| `bl_rao_sigma_d` | baseline | contour distance sum | no pore-facing residue points measured on this frame | bl_rao_flag_n | active |
| `bl_rao_flag_n` | baseline | count | no pore-facing residue points measured on this frame | — | active |

### kcsa_na-stride4 — 104 columns (count derived)

| column | block | units | missing means | indicator | verdict |
|---|---|---|---|---|---|
| `geo_r_constriction_A` | geometry | A | no search slice inside the pore | geo_n_search | active |
| `geo_z_constriction_A` | geometry | A | no search slice inside the pore | geo_n_search | active |
| `geo_n_search` | geometry | count | never missing: the counter is always defined | — | active |
| `geo_search_frac` | geometry | frac | never missing | — | active |
| `geo_boundary_frac` | geometry | frac | never missing | — | active |
| `geo_lining_carbonyl_cos` | geometry | cos | no lining carbonyls in the band | geo_lining_n | active |
| `geo_lining_n` | geometry | count | never missing | — | active |
| `geo_r_bin0_A` | geometry | A | no search slices in the bin | geo_bin0_nsearch | active |
| `geo_r_bin1_A` | geometry | A | no search slices in the bin | geo_bin1_nsearch | active |
| `geo_r_bin2_A` | geometry | A | no search slices in the bin | geo_bin2_nsearch | active |
| `geo_r_bin3_A` | geometry | A | no search slices in the bin | geo_bin3_nsearch | active |
| `geo_r_bin4_A` | geometry | A | no search slices in the bin | geo_bin4_nsearch | active |
| `geo_r_bin5_A` | geometry | A | no search slices in the bin | geo_bin5_nsearch | active |
| `geo_r_bin6_A` | geometry | A | no search slices in the bin | geo_bin6_nsearch | active |
| `geo_r_bin7_A` | geometry | A | no search slices in the bin | geo_bin7_nsearch | active |
| `geo_bin0_nsearch` | geometry | count | never missing | — | active |
| `geo_bin1_nsearch` | geometry | count | never missing | — | active |
| `geo_bin2_nsearch` | geometry | count | never missing | — | active |
| `geo_bin3_nsearch` | geometry | count | never missing | — | active |
| `geo_bin4_nsearch` | geometry | count | never missing | — | active |
| `geo_bin5_nsearch` | geometry | count | never missing | — | active |
| `geo_bin6_nsearch` | geometry | count | never missing | — | active |
| `geo_bin7_nsearch` | geometry | count | never missing | — | active |
| `geo_r_bin0_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin0_nsearch | active |
| `geo_r_bin1_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin1_nsearch | active |
| `geo_r_bin2_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin2_nsearch | active |
| `geo_r_bin3_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin3_nsearch | active |
| `geo_r_bin4_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin4_nsearch | active |
| `geo_r_bin5_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin5_nsearch | active |
| `geo_r_bin6_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin6_nsearch | active |
| `geo_r_bin7_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin7_nsearch | active |
| `hyd_water_pore_n` | hydration | count | never missing | — | active |
| `hyd_min_bin_water_n` | hydration | count | never missing | — | active |
| `hyd_wet_frac_win` | hydration | frac | never missing | — | active |
| `hyd_wet_logodds_win` | hydration | logit | never missing | — | active |
| `hyd_lining_hydrophobicity` | hydration | WW[-1,1] | no pore-facing residues in the band | hyd_lining_n_facing | active |
| `hyd_lining_n_facing` | hydration | count | never missing | — | active |
| `hyd_water_pore_density` | hydration | A^-3 | no searched profile volume in the bin range | geo_n_search | active |
| `hyd_lining_facing_per_sub` | hydration | count/subunit | subunit count unknown | — | active |
| `hyd_water_bin0_n` | hydration | count | never missing | — | active |
| `hyd_water_bin1_n` | hydration | count | never missing | — | active |
| `hyd_water_bin2_n` | hydration | count | never missing | — | active |
| `hyd_water_bin3_n` | hydration | count | never missing | — | active |
| `hyd_water_bin4_n` | hydration | count | never missing | — | active |
| `hyd_water_bin5_n` | hydration | count | never missing | — | active |
| `hyd_water_bin6_n` | hydration | count | never missing | — | active |
| `hyd_water_bin7_n` | hydration | count | never missing | — | active |
| `hyd_water_bin0_density` | hydration | A^-3 | no search slices in the bin | geo_bin0_nsearch | active |
| `hyd_water_bin1_density` | hydration | A^-3 | no search slices in the bin | geo_bin1_nsearch | active |
| `hyd_water_bin2_density` | hydration | A^-3 | no search slices in the bin | geo_bin2_nsearch | active |
| `hyd_water_bin3_density` | hydration | A^-3 | no search slices in the bin | geo_bin3_nsearch | active |
| `hyd_water_bin4_density` | hydration | A^-3 | no search slices in the bin | geo_bin4_nsearch | active |
| `hyd_water_bin5_density` | hydration | A^-3 | no search slices in the bin | geo_bin5_nsearch | active |
| `hyd_water_bin6_density` | hydration | A^-3 | no search slices in the bin | geo_bin6_nsearch | active |
| `hyd_water_bin7_density` | hydration | A^-3 | no search slices in the bin | geo_bin7_nsearch | active |
| `ww_n_waters` | water_wire | count | never missing | — | active |
| `ww_max_gap_A` | water_wire | A | fewer than two waters in the pore | ww_n_waters | active |
| `ww_continuous` | water_wire | 0/1 | fewer than two waters in the pore | ww_n_waters | active |
| `ww_dipole_cos_mean` | water_wire | cos | no water in the pore | ww_n_waters | active |
| `ww_dipole_cos_std` | water_wire | cos | fewer than two waters in the pore | ww_n_waters | active |
| `ww_water_linear_density` | water_wire | A^-1 | never missing | — | active |
| `ww_max_gap_rel` | water_wire | frac | fewer than two waters in the pore | ww_n_waters | active |
| `occ_n_ions_pore` | occupancy | count | never missing | — | active |
| `occ_has_innermost` | occupancy | 0/1 | never missing | — | constant |
| `occ_innermost_z_A` | occupancy | A | no innermost ion (no ion in pore, or constriction not measured) | occ_has_innermost | active |
| `occ_innermost_dz_constr_A` | occupancy | A | no innermost ion | occ_has_innermost | active |
| `occ_innermost_coord_n` | occupancy | count | no innermost ion | occ_has_innermost | active |
| `occ_desolvation` | occupancy | frac | no innermost ion, or no ions outside the pore for the bulk reference | occ_has_innermost | active |
| `ns_filter_n_ions` | named_sites | count | no filter (structural mask) | — | active |
| `ns_S0_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ns_S1_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ns_S2_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ns_S3_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ns_S4_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ele_pot_protein_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_water_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_ions_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_lipid_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_axial_charge_asym_e` | electrostatics | e | constriction not measured | geo_n_search | active |
| `sym_com_radial_spread_A` | symmetry | A | fewer than two subunits (structural mask) | — | active |
| `sym_com_z_spread_A` | symmetry | A | fewer than two subunits (structural mask) | — | active |
| `sym_nn_dist_cv` | symmetry | frac | fewer than three subunits (structural mask) | — | active |
| `sym_com_z_spread_rel` | symmetry | frac | fewer than two subunits (structural mask) | — | active |
| `flu_water_pore_var_win` | fluctuations | count^2 | never missing | — | active |
| `flu_ion_pore_var_win` | fluctuations | count^2 | never missing | — | active |
| `flu_r_constr_var_win` | fluctuations | A^2 | no measured radius values in the window | geo_n_search | active |
| `dyn_ion_persist` | dynamics | 0/1 | never missing | — | active |
| `dyn_ion_drift_A_win` | dynamics | A | ion not in the pore now, or was not at the window start | dyn_ion_persist | active |
| `dyn_wet_state_age_ps` | dynamics | ps | never missing | — | active |
| `dyn_occ_change_age_ps` | dynamics | ps | never missing | — | active |
| `dlv_n_upper_mouth` | delivery | count | never missing | — | active |
| `dlv_n_lower_mouth` | delivery | count | never missing | — | active |
| `dlv_n_entry_mouth` | delivery | count | never missing | — | active |
| `dlv_entry_wide_n` | delivery | count | never missing | — | active |
| `dlv_dist_entry_A` | delivery | A | no permeant in the wide entry zone | dlv_entry_wide_n | active |
| `dlv_approach_pair_n` | delivery | count | never missing | — | active |
| `dlv_approach_A_win` | delivery | A | distance not measured at one of the window endpoints | dlv_approach_pair_n | active |
| `dlv_t_since_cross_ps` | delivery | ps | before the trajectory's first crossing the measured object does not exist | dlv_crossed_before | active |
| `dlv_crossed_before` | delivery | 0/1 | never missing | — | active |
| `dlv_rate_win` | delivery | 1/ns | never missing | — | active |
| `bl_rao_E_kJmol` | baseline | kJ/mol | inputs not measured on this frame | geo_n_search | active |
| `bl_rao_E_win_kJmol` | baseline | kJ/mol | no measured inputs in the window | geo_n_search | active |
| `bl_rao_sigma_d` | baseline | contour distance sum | no pore-facing residue points measured on this frame | bl_rao_flag_n | active |
| `bl_rao_flag_n` | baseline | count | no pore-facing residue points measured on this frame | — | active |

### mthk-stride5 — 104 columns (count derived)

| column | block | units | missing means | indicator | verdict |
|---|---|---|---|---|---|
| `geo_r_constriction_A` | geometry | A | no search slice inside the pore | geo_n_search | active |
| `geo_z_constriction_A` | geometry | A | no search slice inside the pore | geo_n_search | active |
| `geo_n_search` | geometry | count | never missing: the counter is always defined | — | active |
| `geo_search_frac` | geometry | frac | never missing | — | active |
| `geo_boundary_frac` | geometry | frac | never missing | — | active |
| `geo_lining_carbonyl_cos` | geometry | cos | no lining carbonyls in the band | geo_lining_n | active |
| `geo_lining_n` | geometry | count | never missing | — | active |
| `geo_r_bin0_A` | geometry | A | no search slices in the bin | geo_bin0_nsearch | active |
| `geo_r_bin1_A` | geometry | A | no search slices in the bin | geo_bin1_nsearch | active |
| `geo_r_bin2_A` | geometry | A | no search slices in the bin | geo_bin2_nsearch | active |
| `geo_r_bin3_A` | geometry | A | no search slices in the bin | geo_bin3_nsearch | active |
| `geo_r_bin4_A` | geometry | A | no search slices in the bin | geo_bin4_nsearch | active |
| `geo_r_bin5_A` | geometry | A | no search slices in the bin | geo_bin5_nsearch | active |
| `geo_r_bin6_A` | geometry | A | no search slices in the bin | geo_bin6_nsearch | active |
| `geo_r_bin7_A` | geometry | A | no search slices in the bin | geo_bin7_nsearch | active |
| `geo_bin0_nsearch` | geometry | count | never missing | — | active |
| `geo_bin1_nsearch` | geometry | count | never missing | — | active |
| `geo_bin2_nsearch` | geometry | count | never missing | — | active |
| `geo_bin3_nsearch` | geometry | count | never missing | — | active |
| `geo_bin4_nsearch` | geometry | count | never missing | — | active |
| `geo_bin5_nsearch` | geometry | count | never missing | — | constant |
| `geo_bin6_nsearch` | geometry | count | never missing | — | active |
| `geo_bin7_nsearch` | geometry | count | never missing | — | active |
| `geo_r_bin0_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin0_nsearch | active |
| `geo_r_bin1_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin1_nsearch | active |
| `geo_r_bin2_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin2_nsearch | active |
| `geo_r_bin3_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin3_nsearch | active |
| `geo_r_bin4_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin4_nsearch | active |
| `geo_r_bin5_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin5_nsearch | active |
| `geo_r_bin6_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin6_nsearch | active |
| `geo_r_bin7_rel` | geometry | ratio | no search slices in the bin (or no constriction) | geo_bin7_nsearch | active |
| `hyd_water_pore_n` | hydration | count | never missing | — | active |
| `hyd_min_bin_water_n` | hydration | count | never missing | — | constant |
| `hyd_wet_frac_win` | hydration | frac | never missing | — | constant |
| `hyd_wet_logodds_win` | hydration | logit | never missing | — | active |
| `hyd_lining_hydrophobicity` | hydration | WW[-1,1] | no pore-facing residues in the band | hyd_lining_n_facing | active |
| `hyd_lining_n_facing` | hydration | count | never missing | — | active |
| `hyd_water_pore_density` | hydration | A^-3 | no searched profile volume in the bin range | geo_n_search | active |
| `hyd_lining_facing_per_sub` | hydration | count/subunit | subunit count unknown | — | active |
| `hyd_water_bin0_n` | hydration | count | never missing | — | active |
| `hyd_water_bin1_n` | hydration | count | never missing | — | active |
| `hyd_water_bin2_n` | hydration | count | never missing | — | active |
| `hyd_water_bin3_n` | hydration | count | never missing | — | active |
| `hyd_water_bin4_n` | hydration | count | never missing | — | active |
| `hyd_water_bin5_n` | hydration | count | never missing | — | constant |
| `hyd_water_bin6_n` | hydration | count | never missing | — | active |
| `hyd_water_bin7_n` | hydration | count | never missing | — | active |
| `hyd_water_bin0_density` | hydration | A^-3 | no search slices in the bin | geo_bin0_nsearch | active |
| `hyd_water_bin1_density` | hydration | A^-3 | no search slices in the bin | geo_bin1_nsearch | active |
| `hyd_water_bin2_density` | hydration | A^-3 | no search slices in the bin | geo_bin2_nsearch | active |
| `hyd_water_bin3_density` | hydration | A^-3 | no search slices in the bin | geo_bin3_nsearch | active |
| `hyd_water_bin4_density` | hydration | A^-3 | no search slices in the bin | geo_bin4_nsearch | active |
| `hyd_water_bin5_density` | hydration | A^-3 | no search slices in the bin | geo_bin5_nsearch | constant |
| `hyd_water_bin6_density` | hydration | A^-3 | no search slices in the bin | geo_bin6_nsearch | active |
| `hyd_water_bin7_density` | hydration | A^-3 | no search slices in the bin | geo_bin7_nsearch | active |
| `ww_n_waters` | water_wire | count | never missing | — | active |
| `ww_max_gap_A` | water_wire | A | fewer than two waters in the pore | ww_n_waters | active |
| `ww_continuous` | water_wire | 0/1 | fewer than two waters in the pore | ww_n_waters | constant |
| `ww_dipole_cos_mean` | water_wire | cos | no water in the pore | ww_n_waters | active |
| `ww_dipole_cos_std` | water_wire | cos | fewer than two waters in the pore | ww_n_waters | active |
| `ww_water_linear_density` | water_wire | A^-1 | never missing | — | active |
| `ww_max_gap_rel` | water_wire | frac | fewer than two waters in the pore | ww_n_waters | active |
| `occ_n_ions_pore` | occupancy | count | never missing | — | active |
| `occ_has_innermost` | occupancy | 0/1 | never missing | — | constant |
| `occ_innermost_z_A` | occupancy | A | no innermost ion (no ion in pore, or constriction not measured) | occ_has_innermost | active |
| `occ_innermost_dz_constr_A` | occupancy | A | no innermost ion | occ_has_innermost | active |
| `occ_innermost_coord_n` | occupancy | count | no innermost ion | occ_has_innermost | active |
| `occ_desolvation` | occupancy | frac | no innermost ion, or no ions outside the pore for the bulk reference | occ_has_innermost | active |
| `ns_filter_n_ions` | named_sites | count | no filter (structural mask) | — | active |
| `ns_S0_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ns_S1_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ns_S2_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ns_S3_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ns_S4_occ` | named_sites | count | no filter, or site beyond the ring count (structural mask) | — | active |
| `ele_pot_protein_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_water_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_ions_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_pot_lipid_cut12_e_per_A` | electrostatics | e/A | no innermost ion | occ_has_innermost | active |
| `ele_axial_charge_asym_e` | electrostatics | e | constriction not measured | geo_n_search | active |
| `sym_com_radial_spread_A` | symmetry | A | fewer than two subunits (structural mask) | — | active |
| `sym_com_z_spread_A` | symmetry | A | fewer than two subunits (structural mask) | — | active |
| `sym_nn_dist_cv` | symmetry | frac | fewer than three subunits (structural mask) | — | active |
| `sym_com_z_spread_rel` | symmetry | frac | fewer than two subunits (structural mask) | — | active |
| `flu_water_pore_var_win` | fluctuations | count^2 | never missing | — | active |
| `flu_ion_pore_var_win` | fluctuations | count^2 | never missing | — | active |
| `flu_r_constr_var_win` | fluctuations | A^2 | no measured radius values in the window | geo_n_search | active |
| `dyn_ion_persist` | dynamics | 0/1 | never missing | — | active |
| `dyn_ion_drift_A_win` | dynamics | A | ion not in the pore now, or was not at the window start | dyn_ion_persist | active |
| `dyn_wet_state_age_ps` | dynamics | ps | never missing | — | active |
| `dyn_occ_change_age_ps` | dynamics | ps | never missing | — | active |
| `dlv_n_upper_mouth` | delivery | count | never missing | — | active |
| `dlv_n_lower_mouth` | delivery | count | never missing | — | active |
| `dlv_n_entry_mouth` | delivery | count | never missing | — | active |
| `dlv_entry_wide_n` | delivery | count | never missing | — | active |
| `dlv_dist_entry_A` | delivery | A | no permeant in the wide entry zone | dlv_entry_wide_n | active |
| `dlv_approach_pair_n` | delivery | count | never missing | — | active |
| `dlv_approach_A_win` | delivery | A | distance not measured at one of the window endpoints | dlv_approach_pair_n | active |
| `dlv_t_since_cross_ps` | delivery | ps | before the trajectory's first crossing the measured object does not exist | dlv_crossed_before | active |
| `dlv_crossed_before` | delivery | 0/1 | never missing | — | active |
| `dlv_rate_win` | delivery | 1/ns | never missing | — | active |
| `bl_rao_E_kJmol` | baseline | kJ/mol | inputs not measured on this frame | geo_n_search | active |
| `bl_rao_E_win_kJmol` | baseline | kJ/mol | no measured inputs in the window | geo_n_search | active |
| `bl_rao_sigma_d` | baseline | contour distance sum | no pore-facing residue points measured on this frame | bl_rao_flag_n | active |
| `bl_rao_flag_n` | baseline | count | no pore-facing residue points measured on this frame | — | active |
<!-- END GENERATED -->
