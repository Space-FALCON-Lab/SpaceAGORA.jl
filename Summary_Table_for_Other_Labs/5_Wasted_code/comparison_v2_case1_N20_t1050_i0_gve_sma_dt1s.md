# ORACLE Prototype and SpaceAGORA Comparison

Fresh runs of both models using the shared case below.

Historical report retained with its original results and runtimes. The analysis CSV directories linked below are no longer present; the Feather bundles remain available. See the [saved dt1s versus dt10s comparison](dt1s_vs_dt10s.md) for analysis reconstructed from those bundles, and the [dt10s summary](comparison_v2_case1_N20_t1050_i0_gve_sma_dt10s.md) for the existing ten-second outputs.

## 1. Scenario Conditions

| Setting | Value |
| --- | --- |
| helpers | 20 |
| helper_altitude_km | 1000 |
| target_altitude_km | 1050 |
| helper_inclination_deg | 0 |
| target_inclination_deg | 0 |
| target_ecc | 0 |
| target_nu_deg | 0 |
| orbits | 50 |
| laser_range_km | 200 |
| laser_power_w | 10000 |
| magnification | 100 |
| mass_kg | 227 |
| prototype_schedule | gve_sma |
| spaceagora_schedule | gve_sma |
| spaceagora_dt_max_s | 10 |
| output_interval_s | 1 |
| target_period_s | 6371.34083 |
| duration_s | 318567.042 |

Requested orbits count initial target Keplerian periods using the prototype reference Earth radius (6378137 m) and mu (3.986004418e14 m^3/s^2). Both models use the same derived duration in seconds, not a counted-revolution stopping condition. Smoke runs cap this duration at 60 seconds. Folder names retain the derived duration in seconds.

J2 enabled; drag disabled; output every 1.0 seconds plus the final endpoint. Helpers are circular and evenly spaced in anomaly; RAAN and argument of periapsis are zero. Altitudes specify initial semimajor axis minus each model's Earth radius. SpaceAGORA beta and eta are 1. ORACLE here means Kuang's prototype.

| Model | Status | Solver | Simulated duration (s) | Model run wall time (s) | Earth radius (m) | mu (m^3/s^2) | Light speed (m/s) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ORACLE prototype | Success | Vern9 | 318567.042 | 106.418467 | 6378137 | 3.98600442e+14 | 300000000 |
| SpaceAGORA | Success | Tsit5 | 318567.042 | 121.094978 | 6378136.6 | 3.98600436e+14 | 299792458 |

Model run wall time includes setup, simulation, recording, analysis, and any requested video rendering. It excludes worker startup and imports, but includes first-call compilation; it is not a warmed solver-only benchmark.

## 2. Saved Analysis Results

### ORACLE Prototype

Feather: [1_Kuang's Prototype Code/output/feather/prototype_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s](../1_Kuang's%20Prototype%20Code/output/feather/prototype_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s)/

Analysis CSVs: [1_Kuang's Prototype Code/output/CSV/prototype_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s/analysis](../1_Kuang's%20Prototype%20Code/output/CSV/prototype_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s/analysis)/

Geometric encounters: 11; saved encounter-state rows: 53898.

| Metric | Extremum | Value | sc_a | sc_b | Time (s) |
| --- | --- | --- | --- | --- | --- |
| relative_speed_mps | min | 16.2480891 | 11 | 21 | 309687 |
| relative_speed_mps | max | 211.308821 | 12 | 21 | 275713.566 |
| range_rate_mps | min | -92.8854155 | 11 | 21 | 306585.865 |
| range_rate_mps | max | 103.790955 | 11 | 21 | 311896 |
| geometry_duration_s | min | 2591.86384 | 1 | 21 | 0 |
| geometry_duration_s | max | 5368.35988 | 11 | 21 | 306585.865 |
| laser_on_s | min | 0 | 1 | 21 | 0 |
| laser_on_s | max | 3432.29205 | 12 | 21 | 275713.566 |

| Encounter | sc_a | sc_b | Start (s) | End (s) | Geometry (s) | Laser on (s) | Start clipped | End clipped |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 21 | 0 | 2591.86384 | 2591.86384 | 0 | true | false |
| 11 | 11 | 21 | 306585.865 | 311954.225 | 5368.35988 | 3409.28367 | false | false |
| 10 | 12 | 21 | 275713.566 | 280930.349 | 5216.78295 | 3432.29205 | false | false |
| 9 | 13 | 21 | 245237.991 | 250012.322 | 4774.33088 | 2852.54814 | false | false |
| 8 | 14 | 21 | 214803.066 | 219375.067 | 4572.00062 | 2170.08803 | false | false |
| 7 | 15 | 21 | 183913.422 | 188919.208 | 5005.78664 | 1956.10772 | false | false |
| 6 | 16 | 21 | 152840.212 | 158081.144 | 5240.93272 | 2004.05728 | false | false |
| 5 | 17 | 21 | 121708.159 | 127023.072 | 5314.91254 | 2246.98783 | false | false |
| 4 | 18 | 21 | 90577.9857 | 95893.0081 | 5315.02235 | 2580.71655 | false | false |
| 3 | 19 | 21 | 59487.59 | 64756.9192 | 5269.32922 | 2762.15793 | false | false |
| 2 | 20 | 21 | 28446.4929 | 33651.2895 | 5204.79663 | 2727.62398 | false | false |

### SpaceAGORA

Feather: [2_SpaceAGORA.jl/output/feather/spaceagora_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s](../2_SpaceAGORA.jl/output/feather/spaceagora_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s)/

Analysis CSVs: [2_SpaceAGORA.jl/output/CSV/spaceagora_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s/analysis](../2_SpaceAGORA.jl/output/CSV/spaceagora_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s/analysis)/

Geometric encounters: 11; saved encounter-state rows: 53894.

| Metric | Extremum | Value | sc_a | sc_b | Time (s) |
| --- | --- | --- | --- | --- | --- |
| relative_speed_mps | min | 16.3945998 | 11 | 21 | 309697 |
| relative_speed_mps | max | 211.289426 | 12 | 21 | 275720.925 |
| range_rate_mps | min | -92.7039875 | 11 | 21 | 306592.57 |
| range_rate_mps | max | 103.738403 | 11 | 21 | 311897 |
| geometry_duration_s | min | 2591.86361 | 1 | 21 | 0 |
| geometry_duration_s | max | 5367.85861 | 11 | 21 | 306592.57 |
| laser_on_s | min | 0.000379693605 | 1 | 21 | 0 |
| laser_on_s | max | 3422.33646 | 12 | 21 | 275720.925 |

| Encounter | sc_a | sc_b | Start (s) | End (s) | Geometry (s) | Laser on (s) | Start clipped | End clipped |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 21 | 0 | 2591.86361 | 2591.86361 | 0.000379693605 | true | false |
| 11 | 11 | 21 | 306592.57 | 311960.428 | 5367.85861 | 3404.93298 | false | false |
| 10 | 12 | 21 | 275720.925 | 280935.901 | 5214.97595 | 3422.33646 | false | false |
| 9 | 13 | 21 | 245247.187 | 250017.896 | 4770.70863 | 2845.65675 | false | false |
| 8 | 14 | 21 | 214808 | 219381.399 | 4573.39919 | 2165.73238 | false | false |
| 7 | 15 | 21 | 183915.956 | 188923.412 | 5007.45532 | 1952.08676 | false | false |
| 6 | 16 | 21 | 152841.67 | 158083.19 | 5241.51984 | 1996.81823 | false | false |
| 5 | 17 | 21 | 121709.01 | 127024.16 | 5315.15039 | 2243.11775 | false | false |
| 4 | 18 | 21 | 90578.4418 | 95893.5833 | 5315.14145 | 2575.68941 | false | false |
| 3 | 19 | 21 | 59487.7651 | 64757.169 | 5269.40394 | 2762.89935 | false | false |
| 2 | 20 | 21 | 28446.4884 | 33651.3591 | 5204.87065 | 2719.57648 | false | false |

## 3. Comparison

Differences are SpaceAGORA minus prototype. Percent differences use the absolute prototype value; zero or unavailable references have no percentage.

| Metric | Extremum | Prototype | SpaceAGORA | Difference | Difference (%) |
| --- | --- | --- | --- | --- | --- |
| geometry_duration_s | max | 5368.35988 | 5367.85861 | -0.501268786 | -0.00933746613 |
| geometry_duration_s | min | 2591.86384 | 2591.86361 | -0.000232509097 | -8.97072959e-06 |
| laser_on_s | max | 3432.29205 | 3422.33646 | -9.9555955 | -0.290056771 |
| laser_on_s | min | 0 | 0.000379693605 | 0.000379693605 | unavailable |
| range_rate_mps | max | 103.790955 | 103.738403 | -0.0525516316 | -0.0506321881 |
| range_rate_mps | min | -92.8854155 | -92.7039875 | 0.181428012 | 0.195324542 |
| relative_speed_mps | max | 211.308821 | 211.289426 | -0.0193947948 | -0.0091784123 |
| relative_speed_mps | min | 16.2480891 | 16.3945998 | 0.146510728 | 0.901710514 |

Encounters are matched by unordered spacecraft pair and chronological occurrence, not by run-local IDs. Unmatched occurrences show unavailable differences; a missed encounter can shift subsequent occurrence matching.

| sc_a | sc_b | Occurrence | Entry difference (s) | Exit difference (s) | Geometry difference (s) | Laser-on difference (s) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 21 | 1 | 0 | -0.000232509097 | -0.000232509097 | 0.000379693605 |
| 11 | 21 | 1 | 6.70486261 | 6.20359382 | -0.501268786 | -4.35069772 |
| 12 | 21 | 1 | 7.35913189 | 5.55213174 | -1.80700016 | -9.9555955 |
| 13 | 21 | 1 | 9.19606862 | 5.57382531 | -3.6222433 | -6.89139384 |
| 14 | 21 | 1 | 4.93422134 | 6.33278809 | 1.39856675 | -4.35564729 |
| 15 | 21 | 1 | 2.53479236 | 4.20347152 | 1.66867916 | -4.02095796 |
| 16 | 21 | 1 | 1.45868149 | 2.0457999 | 0.58711841 | -7.23905147 |
| 17 | 21 | 1 | 0.850633746 | 1.08848994 | 0.237856195 | -3.8700831 |
| 18 | 21 | 1 | 0.456090697 | 0.575193962 | 0.119103265 | -5.02713606 |
| 19 | 21 | 1 | 0.175065155 | 0.249792069 | 0.0747269148 | 0.741412144 |
| 20 | 21 | 1 | -0.00446266908 | 0.069563516 | 0.0740261851 | -8.04749713 |

No numerical acceptance tolerance was specified: these tables quantify agreement, not a physical-equivalence verdict. Motion extrema use only encounter states and entry/exit boundaries; unrelated pairs and out-of-encounter times are excluded. Clipped encounters cover only the recorded portion. Earth/light-speed constants, solvers, LOS checks, and scheduling differ. SpaceAGORA uses accepted-step scheduling and a velocity-kick callback; prototype forces are evaluated in the ODE. The report reads saved analysis CSVs after both worker processes exit; it does not compare full trajectories or equate the two delta-v diagnostics.

