# Test Case 6: ORACLE Prototype and SpaceAGORA Comparison

Fresh runs of both models using the shared case below.

## 1. Scenario Conditions

| Setting | Value |
| --- | --- |
| helpers | 20 |
| helper_altitude_km | 1000 |
| target_altitude_km | 1050 |
| helper_inclination_deg | 0 |
| target_inclination_deg | 5 |
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
| output_interval_s | 10 |
| target_period_s | 6371.34083 |
| duration_s | 318567.042 |

Requested orbits count initial target Keplerian periods using the prototype reference Earth radius (6378137 m) and mu (3.986004418e14 m^3/s^2). Both models use the same derived duration in seconds, not a counted-revolution stopping condition. Smoke runs cap this duration at 60 seconds. Folder names retain the derived duration in seconds.

J2 enabled; drag disabled; output every 10.0 seconds plus the final endpoint. Helpers are circular and evenly spaced in anomaly; RAAN and argument of periapsis are zero. Altitudes specify initial semimajor axis minus each model's Earth radius. SpaceAGORA beta and eta are 1. ORACLE here means Kuang's prototype.

| Model | Status | Solver | Simulated duration (s) | Model run wall time (s) | Earth radius (m) | mu (m^3/s^2) | Light speed (m/s) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ORACLE prototype | Success | Vern9 | 318567.042 | 41.4454157 | 6378137 | 3.98600442e+14 | 300000000 |
| SpaceAGORA | Success | Tsit5 | 318567.042 | 65.5264655 | 6378136.6 | 3.98600436e+14 | 299792458 |

Model run wall time includes setup, simulation, recording, analysis, and any requested video rendering. It excludes worker startup and imports, but includes first-call compilation; it is not a warmed solver-only benchmark.

## 2. Saved Analysis Results

### ORACLE Prototype

Feather: [../1_Kuang's Prototype Code/output/feather/prototype_N20_h1000km_t1050km_ih0_it5_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s](../1_Kuang's%20Prototype%20Code/output/feather/prototype_N20_h1000km_t1050km_ih0_it5_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s)/

Analysis CSVs: [../1_Kuang's Prototype Code/output/CSV/prototype_N20_h1000km_t1050km_ih0_it5_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis](../1_Kuang's%20Prototype%20Code/output/CSV/prototype_N20_h1000km_t1050km_ih0_it5_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis)/

Geometric encounters: 18; saved encounter-state rows: 843.

| Metric | Extremum | Value | sc_a | sc_b | Time (s) |
| --- | --- | --- | --- | --- | --- |
| relative_speed_mps | min | 611.746378 | 17 | 21 | 124108.41 |
| relative_speed_mps | max | 670.583882 | 17 | 21 | 127000 |
| range_rate_mps | min | -601.090628 | 17 | 21 | 123530 |
| range_rate_mps | max | 598.910037 | 17 | 21 | 124070 |
| geometry_duration_s | min | 160.186783 | 17 | 21 | 126841.625 |
| geometry_duration_s | max | 611.151213 | 17 | 21 | 123497.258 |
| laser_on_s | min | 0 | 11 | 21 | 310778.424 |
| laser_on_s | max | 589.989012 | 12 | 21 | 279050.743 |

| Encounter | sc_a | sc_b | Start (s) | End (s) | Geometry (s) | Laser on (s) | Start clipped | End clipped |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 21 | 0 | 305.070847 | 305.070847 | 9.93768208e-05 | true | false |
| 2 | 20 | 21 | 28496.7475 | 28731.09 | 234.342423 | 234.342423 | false | false |
| 3 | 20 | 21 | 31434.5917 | 32026.8935 | 592.301847 | 0 | false | false |
| 4 | 19 | 21 | 60108.507 | 60578.928 | 470.420999 | 470.420999 | false | false |
| 5 | 19 | 21 | 63193.8468 | 63728.7912 | 534.944458 | 0 | false | false |
| 6 | 18 | 21 | 91786.5742 | 92359.7418 | 573.167597 | 573.167597 | false | false |
| 7 | 18 | 21 | 94980.869 | 95402.4138 | 421.544749 | 0 | false | false |
| 8 | 17 | 21 | 123497.258 | 124108.41 | 611.151213 | 541.530512 | false | false |
| 9 | 17 | 21 | 126841.625 | 127001.812 | 160.186783 | 0 | false | false |
| 10 | 16 | 21 | 155235.183 | 155831.467 | 596.283727 | 0 | false | false |
| 11 | 15 | 21 | 184074.604 | 184246.632 | 172.027812 | 172.027812 | false | false |
| 12 | 15 | 21 | 187001.961 | 187528.13 | 526.169093 | 0 | false | false |
| 13 | 14 | 21 | 215675.214 | 216106.248 | 431.033839 | 431.033839 | false | false |
| 14 | 14 | 21 | 218810.543 | 219185.883 | 375.340428 | 0 | false | false |
| 15 | 13 | 21 | 247347.894 | 247890.116 | 542.221757 | 542.221757 | false | false |
| 16 | 12 | 21 | 279050.743 | 279640.732 | 589.989012 | 589.989012 | false | false |
| 17 | 11 | 21 | 307762.165 | 308123.466 | 361.300635 | 361.300635 | false | false |
| 18 | 11 | 21 | 310778.424 | 311364.497 | 586.0736 | 0 | false | false |

### SpaceAGORA

Feather: [../2_SpaceAGORA.jl/output/feather/spaceagora_N20_h1000km_t1050km_ih0_it5_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s](../2_SpaceAGORA.jl/output/feather/spaceagora_N20_h1000km_t1050km_ih0_it5_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s)/

Analysis CSVs: [../2_SpaceAGORA.jl/output/CSV/spaceagora_N20_h1000km_t1050km_ih0_it5_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis](../2_SpaceAGORA.jl/output/CSV/spaceagora_N20_h1000km_t1050km_ih0_it5_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis)/

Geometric encounters: 18; saved encounter-state rows: 842.

| Metric | Extremum | Value | sc_a | sc_b | Time (s) |
| --- | --- | --- | --- | --- | --- |
| relative_speed_mps | min | 611.746476 | 17 | 21 | 124108.4 |
| relative_speed_mps | max | 670.568589 | 17 | 21 | 127000 |
| range_rate_mps | min | -601.085495 | 17 | 21 | 123530 |
| range_rate_mps | max | 598.900873 | 17 | 21 | 124070 |
| geometry_duration_s | min | 160.803792 | 17 | 21 | 126841.313 |
| geometry_duration_s | max | 611.141679 | 17 | 21 | 123497.258 |
| laser_on_s | min | 0 | 11 | 21 | 310778.378 |
| laser_on_s | max | 584.290693 | 12 | 21 | 279050.81 |

| Encounter | sc_a | sc_b | Start (s) | End (s) | Geometry (s) | Laser on (s) | Start clipped | End clipped |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 21 | 0 | 305.070842 | 305.070842 | 0.000379693605 | true | false |
| 2 | 20 | 21 | 28496.7445 | 28731.0884 | 234.343883 | 225.360088 | false | false |
| 3 | 20 | 21 | 31434.5882 | 32026.8931 | 592.304892 | 0 | false | false |
| 4 | 19 | 21 | 60108.5176 | 60578.9094 | 470.391831 | 463.83641 | false | false |
| 5 | 19 | 21 | 63193.8306 | 63728.801 | 534.970409 | 0 | false | false |
| 6 | 18 | 21 | 91786.5865 | 92359.7199 | 573.133374 | 564.321173 | false | false |
| 7 | 18 | 21 | 94980.8116 | 95402.464 | 421.652395 | 0 | false | false |
| 8 | 17 | 21 | 123497.258 | 124108.4 | 611.141679 | 541.299029 | false | false |
| 9 | 17 | 21 | 126841.313 | 127002.117 | 160.803792 | 0 | false | false |
| 10 | 16 | 21 | 155235.148 | 155831.49 | 596.342308 | 0 | false | false |
| 11 | 15 | 21 | 184075.113 | 184246.115 | 171.002721 | 164.032284 | false | false |
| 12 | 15 | 21 | 187001.854 | 187528.223 | 526.368916 | 0 | false | false |
| 13 | 14 | 21 | 215675.418 | 216106.047 | 430.62889 | 421.502769 | false | false |
| 14 | 14 | 21 | 218810.247 | 219186.163 | 375.916549 | 0 | false | false |
| 15 | 13 | 21 | 247348.033 | 247889.994 | 541.960645 | 534.030225 | false | false |
| 16 | 12 | 21 | 279050.81 | 279640.7 | 589.890579 | 584.290693 | false | false |
| 17 | 11 | 21 | 307762.75 | 308122.847 | 360.097071 | 351.811289 | false | false |
| 18 | 11 | 21 | 310778.378 | 311364.596 | 586.218503 | 0 | false | false |

## 3. Comparison

Differences are SpaceAGORA minus prototype. Percent differences use the absolute prototype value; zero or unavailable references have no percentage.

| Metric | Extremum | Prototype | SpaceAGORA | Difference | Difference (%) |
| --- | --- | --- | --- | --- | --- |
| geometry_duration_s | max | 611.151213 | 611.141679 | -0.00953314855 | -0.0015598674 |
| geometry_duration_s | min | 160.186783 | 160.803792 | 0.617008973 | 0.385180951 |
| laser_on_s | max | 589.989012 | 584.290693 | -5.69831883 | -0.965834738 |
| laser_on_s | min | 0 | 0 | 0 | unavailable |
| range_rate_mps | max | 598.910037 | 598.900873 | -0.00916363506 | -0.00153005201 |
| range_rate_mps | min | -601.090628 | -601.085495 | 0.00513279607 | 0.00085391384 |
| relative_speed_mps | max | 670.583882 | 670.568589 | -0.0152930044 | -0.00228055055 |
| relative_speed_mps | min | 611.746378 | 611.746476 | 9.83717847e-05 | 1.6080485e-05 |

Encounters are matched by unordered spacecraft pair and chronological occurrence, not by run-local IDs. Unmatched occurrences show unavailable differences; a missed encounter can shift subsequent occurrence matching.

| sc_a | sc_b | Occurrence | Entry difference (s) | Exit difference (s) | Geometry difference (s) | Laser-on difference (s) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 21 | 1 | 0 | -5.79383914e-06 | -5.79383914e-06 | 0.000280316784 |
| 11 | 21 | 1 | 0.584947571 | -0.618616866 | -1.20356444 | -9.48934627 |
| 11 | 21 | 2 | -0.0460598826 | 0.0988423372 | 0.14490222 | 0 |
| 12 | 21 | 1 | 0.066888004 | -0.0315452073 | -0.0984332113 | -5.69831883 |
| 13 | 21 | 1 | 0.13953469 | -0.121577596 | -0.261112286 | -8.19153162 |
| 14 | 21 | 1 | 0.203940993 | -0.201008578 | -0.40494957 | -9.53107 |
| 14 | 21 | 2 | -0.296372789 | 0.279747929 | 0.576120718 | 0 |
| 15 | 21 | 1 | 0.509013399 | -0.516077375 | -1.02509077 | -7.9955278 |
| 15 | 21 | 2 | -0.106649042 | 0.0931742237 | 0.199823266 | 0 |
| 16 | 21 | 1 | -0.0351082976 | 0.0234720554 | 0.058580353 | 0 |
| 17 | 21 | 1 | -0.000479016613 | -0.0100121652 | -0.00953314855 | -0.231482803 |
| 17 | 21 | 2 | -0.312617543 | 0.30439143 | 0.617008973 | 0 |
| 18 | 21 | 1 | 0.0123054991 | -0.0219173924 | -0.0342228916 | -8.84642348 |
| 18 | 21 | 2 | -0.0574601692 | 0.0501863101 | 0.107646479 | 0 |
| 19 | 21 | 1 | 0.010623062 | -0.0185455469 | -0.0291686088 | -6.58458845 |
| 19 | 21 | 2 | -0.0161311811 | 0.00982045821 | 0.0259516393 | 0 |
| 20 | 21 | 1 | -0.0030610107 | -0.00160134856 | 0.00145966214 | -8.98233586 |
| 20 | 21 | 2 | -0.00350246901 | -0.000458041937 | 0.00304442707 | 0 |

No numerical acceptance tolerance was specified: these tables quantify agreement, not a physical-equivalence verdict. Motion extrema use only encounter states and entry/exit boundaries; unrelated pairs and out-of-encounter times are excluded. Clipped encounters cover only the recorded portion. Earth/light-speed constants, solvers, LOS checks, and scheduling differ. SpaceAGORA uses accepted-step scheduling and a velocity-kick callback; prototype forces are evaluated in the ODE. The report reads saved analysis CSVs after both worker processes exit; it does not compare full trajectories or equate the two delta-v diagnostics.

