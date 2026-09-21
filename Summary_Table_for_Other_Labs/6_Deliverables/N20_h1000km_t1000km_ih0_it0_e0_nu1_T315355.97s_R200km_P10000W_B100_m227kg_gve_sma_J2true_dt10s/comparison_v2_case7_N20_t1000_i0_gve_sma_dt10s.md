# Test Case 7: ORACLE Prototype and SpaceAGORA Comparison

Fresh runs of both models using the shared case below.

## 1. Scenario Conditions

| Setting | Value |
| --- | --- |
| helpers | 20 |
| helper_altitude_km | 1000 |
| target_altitude_km | 1000 |
| helper_inclination_deg | 0 |
| target_inclination_deg | 0 |
| target_ecc | 0 |
| target_nu_deg | 1 |
| orbits | 50 |
| laser_range_km | 200 |
| laser_power_w | 10000 |
| magnification | 100 |
| mass_kg | 227 |
| prototype_schedule | gve_sma |
| spaceagora_schedule | gve_sma |
| spaceagora_dt_max_s | 10 |
| output_interval_s | 10 |
| target_period_s | 6307.11941 |
| duration_s | 315355.97 |

Requested orbits count initial target Keplerian periods using the prototype reference Earth radius (6378137 m) and mu (3.986004418e14 m^3/s^2). Both models use the same derived duration in seconds, not a counted-revolution stopping condition. Smoke runs cap this duration at 60 seconds. Folder names retain the derived duration in seconds.

J2 enabled; drag disabled; output every 10.0 seconds plus the final endpoint. Helpers are circular and evenly spaced in anomaly; RAAN and argument of periapsis are zero. Altitudes specify initial semimajor axis minus each model's Earth radius. SpaceAGORA beta and eta are 1. ORACLE here means Kuang's prototype.

| Model | Status | Solver | Simulated duration (s) | Model run wall time (s) | Earth radius (m) | mu (m^3/s^2) | Light speed (m/s) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ORACLE prototype | Success | Vern9 | 315355.97 | 37.2235339 | 6378137 | 3.98600442e+14 | 300000000 |
| SpaceAGORA | Success | Tsit5 | 315355.97 | 63.3346819 | 6378136.6 | 3.98600436e+14 | 299792458 |

Model run wall time includes setup, simulation, recording, analysis, and any requested video rendering. It excludes worker startup and imports, but includes first-call compilation; it is not a warmed solver-only benchmark.

## 2. Saved Analysis Results

### ORACLE Prototype

Feather: [../1_Kuang's Prototype Code/output/feather/prototype_N20_h1000km_t1000km_ih0_it0_e0_nu1_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s](../../1_Kuang's%20Prototype%20Code/output/feather/prototype_N20_h1000km_t1000km_ih0_it0_e0_nu1_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s)/

Analysis CSVs: [../1_Kuang's Prototype Code/output/CSV/prototype_N20_h1000km_t1000km_ih0_it0_e0_nu1_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis](../../1_Kuang's%20Prototype%20Code/output/CSV/prototype_N20_h1000km_t1000km_ih0_it0_e0_nu1_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis)/

Geometric encounters: 1; saved encounter-state rows: 9640.

| Metric | Extremum | Value | sc_a | sc_b | Time (s) |
| --- | --- | --- | --- | --- | --- |
| relative_speed_mps | min | 1.58610902 | 1 | 21 | 53920 |
| relative_speed_mps | max | 199.79618 | 1 | 21 | 96389.9092 |
| range_rate_mps | min | -4.52966487 | 1 | 21 | 51660 |
| range_rate_mps | max | 5.00477482 | 1 | 21 | 96250 |
| geometry_duration_s | min | 96389.9092 | 1 | 21 | 0 |
| geometry_duration_s | max | 96389.9092 | 1 | 21 | 0 |
| laser_on_s | min | 53911.3865 | 1 | 21 | 0 |
| laser_on_s | max | 53911.3865 | 1 | 21 | 0 |

| Encounter | sc_a | sc_b | Start (s) | End (s) | Geometry (s) | Laser on (s) | Start clipped | End clipped |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 21 | 0 | 96389.9092 | 96389.9092 | 53911.3865 | true | false |

### SpaceAGORA

Feather: [../2_SpaceAGORA.jl/output/feather/spaceagora_N20_h1000km_t1000km_ih0_it0_e0_nu1_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s](../../2_SpaceAGORA.jl/output/feather/spaceagora_N20_h1000km_t1000km_ih0_it0_e0_nu1_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s)/

Analysis CSVs: [../2_SpaceAGORA.jl/output/CSV/spaceagora_N20_h1000km_t1000km_ih0_it0_e0_nu1_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis](../../2_SpaceAGORA.jl/output/CSV/spaceagora_N20_h1000km_t1000km_ih0_it0_e0_nu1_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis)/

Geometric encounters: 1; saved encounter-state rows: 9638.

| Metric | Extremum | Value | sc_a | sc_b | Time (s) |
| --- | --- | --- | --- | --- | --- |
| relative_speed_mps | min | 1.58794739 | 1 | 21 | 53910 |
| relative_speed_mps | max | 199.791419 | 1 | 21 | 96362.7046 |
| range_rate_mps | min | -4.5301047 | 1 | 21 | 51640 |
| range_rate_mps | max | 5.00658064 | 1 | 21 | 96260 |
| geometry_duration_s | min | 96362.7046 | 1 | 21 | 0 |
| geometry_duration_s | max | 96362.7046 | 1 | 21 | 0 |
| laser_on_s | min | 53904.1966 | 1 | 21 | 0 |
| laser_on_s | max | 53904.1966 | 1 | 21 | 0 |

| Encounter | sc_a | sc_b | Start (s) | End (s) | Geometry (s) | Laser on (s) | Start clipped | End clipped |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 21 | 0 | 96362.7046 | 96362.7046 | 53904.1966 | true | false |

## 3. Comparison

Differences are SpaceAGORA minus prototype. Percent differences use the absolute prototype value; zero or unavailable references have no percentage.

| Metric | Extremum | Prototype | SpaceAGORA | Difference | Difference (%) |
| --- | --- | --- | --- | --- | --- |
| geometry_duration_s | max | 96389.9092 | 96362.7046 | -27.2045775 | -0.0282234704 |
| geometry_duration_s | min | 96389.9092 | 96362.7046 | -27.2045775 | -0.0282234704 |
| laser_on_s | max | 53911.3865 | 53904.1966 | -7.18987945 | -0.0133364766 |
| laser_on_s | min | 53911.3865 | 53904.1966 | -7.18987945 | -0.0133364766 |
| range_rate_mps | max | 5.00477482 | 5.00658064 | 0.0018058192 | 0.0360819269 |
| range_rate_mps | min | -4.52966487 | -4.5301047 | -0.000439827324 | -0.00970993078 |
| relative_speed_mps | max | 199.79618 | 199.791419 | -0.00476111657 | -0.00238298679 |
| relative_speed_mps | min | 1.58610902 | 1.58794739 | 0.00183837217 | 0.115904528 |

Encounters are matched by unordered spacecraft pair and chronological occurrence, not by run-local IDs. Unmatched occurrences show unavailable differences; a missed encounter can shift subsequent occurrence matching.

| sc_a | sc_b | Occurrence | Entry difference (s) | Exit difference (s) | Geometry difference (s) | Laser-on difference (s) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 21 | 1 | 0 | -27.2045775 | -27.2045775 | -7.18987945 |

No numerical acceptance tolerance was specified: these tables quantify agreement, not a physical-equivalence verdict. Motion extrema use only encounter states and entry/exit boundaries; unrelated pairs and out-of-encounter times are excluded. Clipped encounters cover only the recorded portion. Earth/light-speed constants, solvers, LOS checks, and scheduling differ. SpaceAGORA uses accepted-step scheduling and a velocity-kick callback; prototype forces are evaluated in the ODE. The report reads saved analysis CSVs after both worker processes exit; it does not compare full trajectories or equate the two delta-v diagnostics.

