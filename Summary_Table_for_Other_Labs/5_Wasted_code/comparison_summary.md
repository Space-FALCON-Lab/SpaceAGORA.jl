# ORACLE Prototype and SpaceAGORA Comparison

Fresh runs of both models using the shared case below.

## 1. Scenario Conditions

| Setting | Value |
| --- | --- |
| helpers | 10 |
| helper_altitude_km | 1000 |
| target_altitude_km | 1050 |
| helper_inclination_deg | 0 |
| target_inclination_deg | 0 |
| target_ecc | 0 |
| target_nu_deg | 0 |
| orbits | 10 |
| laser_range_km | 200 |
| laser_power_w | 10000 |
| magnification | 100 |
| mass_kg | 227 |
| prototype_schedule | gve_sma |
| spaceagora_schedule | gve_sma |
| spaceagora_dt_max_s | 10 |
| target_period_s | 6371.34083 |
| duration_s | 63713.4083 |

Requested orbits count initial target Keplerian periods using the prototype reference Earth radius (6378137 m) and mu (3.986004418e14 m^3/s^2). Both models use the same derived duration in seconds, not a counted-revolution stopping condition. Smoke runs cap this duration at 60 seconds. Folder names retain the derived duration in seconds.

J2 enabled; drag disabled; output every 10 seconds plus the final endpoint. Helpers are circular and evenly spaced in anomaly; RAAN and argument of periapsis are zero. Altitudes specify initial semimajor axis minus each model's Earth radius. SpaceAGORA beta and eta are 1. ORACLE here means Kuang's prototype.

| Model | Status | Solver | Actual duration (s) | Earth radius (m) | mu (m^3/s^2) | Light speed (m/s) |
| --- | --- | --- | --- | --- | --- | --- |
| ORACLE prototype | Success | Vern9 | 63713.4083 | 6378137 | 3.98600442e+14 | 300000000 |
| SpaceAGORA | Success | Tsit5 | 63713.4083 | 6378136.6 | 3.98600436e+14 | 299792458 |

## 2. Saved Analysis Results

### ORACLE Prototype

Feather: [1_Kuang's Prototype Code/output/feather/prototype_N10_h1000km_t1050km_ih0_it0_e0_nu0_T63713.4083s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s](1_Kuang's%20Prototype%20Code/output/feather/prototype_N10_h1000km_t1050km_ih0_it0_e0_nu0_T63713.4083s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s)/

Analysis CSVs: [1_Kuang's Prototype Code/output/CSV/prototype_N10_h1000km_t1050km_ih0_it0_e0_nu0_T63713.4083s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis](1_Kuang's%20Prototype%20Code/output/CSV/prototype_N10_h1000km_t1050km_ih0_it0_e0_nu0_T63713.4083s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis)/

Geometric encounters: 2; saved encounter-state rows: 682.

| Metric | Extremum | Value | sc_a | sc_b | Time (s) |
| --- | --- | --- | --- | --- | --- |
| relative_speed_mps | min | 20.1890097 | 10 | 11 | 62240 |
| relative_speed_mps | max | 198.116238 | 10 | 11 | 59523.2549 |
| range_rate_mps | min | -79.3456586 | 10 | 11 | 59523.2549 |
| range_rate_mps | max | 72.6511715 | 10 | 11 | 63713.4083 |
| geometry_duration_s | min | 2591.86384 | 1 | 11 | 0 |
| geometry_duration_s | max | 4190.15336 | 10 | 11 | 59523.2549 |
| laser_on_s | min | 0 | 1 | 11 | 0 |
| laser_on_s | max | 2772.47094 | 10 | 11 | 59523.2549 |

| Encounter | sc_a | sc_b | Start (s) | End (s) | Geometry (s) | Laser on (s) | Start clipped | End clipped |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 11 | 0 | 2591.86384 | 2591.86384 | 0 | true | false |
| 2 | 10 | 11 | 59523.2549 | 63713.4083 | 4190.15336 | 2772.47094 | false | true |

### SpaceAGORA

Feather: [2_SpaceAGORA.jl/output/feather/spaceagora_N10_h1000km_t1050km_ih0_it0_e0_nu0_T63713.4083s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s](2_SpaceAGORA.jl/output/feather/spaceagora_N10_h1000km_t1050km_ih0_it0_e0_nu0_T63713.4083s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s)/

Analysis CSVs: [2_SpaceAGORA.jl/output/CSV/spaceagora_N10_h1000km_t1050km_ih0_it0_e0_nu0_T63713.4083s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis](2_SpaceAGORA.jl/output/CSV/spaceagora_N10_h1000km_t1050km_ih0_it0_e0_nu0_T63713.4083s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis)/

Geometric encounters: 2; saved encounter-state rows: 682.

| Metric | Extremum | Value | sc_a | sc_b | Time (s) |
| --- | --- | --- | --- | --- | --- |
| relative_speed_mps | min | 20.188195 | 10 | 11 | 62240 |
| relative_speed_mps | max | 198.116251 | 10 | 11 | 59523.2462 |
| range_rate_mps | min | -79.3457149 | 10 | 11 | 59523.2462 |
| range_rate_mps | max | 72.6504062 | 10 | 11 | 63713.4083 |
| geometry_duration_s | min | 2591.86361 | 1 | 11 | 0 |
| geometry_duration_s | max | 4190.16206 | 10 | 11 | 59523.2462 |
| laser_on_s | min | 0.00037460584 | 1 | 11 | 0 |
| laser_on_s | max | 2767.98118 | 10 | 11 | 59523.2462 |

| Encounter | sc_a | sc_b | Start (s) | End (s) | Geometry (s) | Laser on (s) | Start clipped | End clipped |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 11 | 0 | 2591.86361 | 2591.86361 | 0.00037460584 | true | false |
| 2 | 10 | 11 | 59523.2462 | 63713.4083 | 4190.16206 | 2767.98118 | false | true |

## 3. Comparison

Differences are SpaceAGORA minus prototype. Percent differences use the absolute prototype value; zero or unavailable references have no percentage.

| Metric | Extremum | Prototype | SpaceAGORA | Difference | Difference (%) |
| --- | --- | --- | --- | --- | --- |
| geometry_duration_s | max | 4190.15336 | 4190.16206 | 0.0086987436 | 0.000207599647 |
| geometry_duration_s | min | 2591.86384 | 2591.86361 | -0.000232500904 | -8.97041348e-06 |
| laser_on_s | max | 2772.47094 | 2767.98118 | -4.48976371 | -0.161940875 |
| laser_on_s | min | 0 | 0.00037460584 | 0.00037460584 | unavailable |
| range_rate_mps | max | 72.6511715 | 72.6504062 | -0.000765364945 | -0.00105347915 |
| range_rate_mps | min | -79.3456586 | -79.3457149 | -5.63809681e-05 | -7.10574077e-05 |
| relative_speed_mps | max | 198.116238 | 198.116251 | 1.25070671e-05 | 6.31299445e-06 |
| relative_speed_mps | min | 20.1890097 | 20.188195 | -0.00081466992 | -0.00403521486 |

Encounters are matched by unordered spacecraft pair and chronological occurrence, not by run-local IDs. Unmatched occurrences show unavailable differences; a missed encounter can shift subsequent occurrence matching.

| sc_a | sc_b | Occurrence | Entry difference (s) | Exit difference (s) | Geometry difference (s) | Laser-on difference (s) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 11 | 1 | 0 | -0.000232500904 | -0.000232500904 | 0.00037460584 |
| 10 | 11 | 1 | -0.00869874359 | 7.27595761e-12 | 0.0086987436 | -4.48976371 |

No numerical acceptance tolerance was specified: these tables quantify agreement, not a physical-equivalence verdict. Extrema are sampled over all spacecraft pairs, including out-of-range pairs. Clipped encounters cover only the recorded portion. Earth/light-speed constants, solvers, LOS checks, and scheduling differ. SpaceAGORA uses accepted-step scheduling and a velocity-kick callback; prototype forces are evaluated in the ODE. The report reads saved analysis CSVs after both worker processes exit; it does not compare full trajectories or equate the two delta-v diagnostics.

