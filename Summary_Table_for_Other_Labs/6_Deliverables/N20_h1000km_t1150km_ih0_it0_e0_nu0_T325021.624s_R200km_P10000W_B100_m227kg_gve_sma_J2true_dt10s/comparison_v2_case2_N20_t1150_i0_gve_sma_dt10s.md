# Test Case 2: ORACLE Prototype and SpaceAGORA Comparison

Fresh runs of both models using the shared case below.

## 1. Scenario Conditions

| Setting | Value |
| --- | --- |
| helpers | 20 |
| helper_altitude_km | 1000 |
| target_altitude_km | 1150 |
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
| target_period_s | 6500.43248 |
| duration_s | 325021.624 |

Requested orbits count initial target Keplerian periods using the prototype reference Earth radius (6378137 m) and mu (3.986004418e14 m^3/s^2). Both models use the same derived duration in seconds, not a counted-revolution stopping condition. Smoke runs cap this duration at 60 seconds. Folder names retain the derived duration in seconds.

J2 enabled; drag disabled; output every 10 seconds plus the final endpoint. Helpers are circular and evenly spaced in anomaly; RAAN and argument of periapsis are zero. Altitudes specify initial semimajor axis minus each model's Earth radius. SpaceAGORA beta and eta are 1. ORACLE here means Kuang's prototype.

| Model | Status | Solver | Actual duration (s) | Earth radius (m) | mu (m^3/s^2) | Light speed (m/s) |
| --- | --- | --- | --- | --- | --- | --- |
| ORACLE prototype | Success | Vern9 | 325021.624 | 6378137 | 3.98600442e+14 | 300000000 |
| SpaceAGORA | Success | Tsit5 | 325021.624 | 6378136.6 | 3.98600436e+14 | 299792458 |

## 2. Saved Analysis Results

### ORACLE Prototype

Video: [Prototype video](../../1_Kuang's%20Prototype%20Code/output/videos/prototype_N20_h1000km_t1150km_ih0_it0_e0_nu0_T325021.624s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/prototype_N20_h1000km_t1150km_ih0_it0_e0_nu0_T325021.624s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s.mp4)

Feather: [1_Kuang's Prototype Code/output/feather/prototype_N20_h1000km_t1150km_ih0_it0_e0_nu0_T325021.624s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s](../../1_Kuang's%20Prototype%20Code/output/feather/prototype_N20_h1000km_t1150km_ih0_it0_e0_nu0_T325021.624s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s)/

Analysis CSVs: [1_Kuang's Prototype Code/output/CSV/prototype_N20_h1000km_t1150km_ih0_it0_e0_nu0_T325021.624s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis](../../1_Kuang's%20Prototype%20Code/output/CSV/prototype_N20_h1000km_t1150km_ih0_it0_e0_nu0_T325021.624s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis)/

Geometric encounters: 31; saved encounter-state rows: 3667.

| Metric | Extremum | Value | sc_a | sc_b | Time (s) |
| --- | --- | --- | --- | --- | --- |
| relative_speed_mps | min | 57.9065854 | 9 | 21 | 126530 |
| relative_speed_mps | max | 166.660069 | 9 | 21 | 127281.656 |
| range_rate_mps | min | -158.25261 | 13 | 21 | 294334.742 |
| range_rate_mps | max | 158.758846 | 11 | 21 | 316871.316 |
| geometry_duration_s | min | 598.328079 | 1 | 21 | 0 |
| geometry_duration_s | max | 1503.31117 | 9 | 21 | 125778.345 |
| laser_on_s | min | 0 | 1 | 21 | 0 |
| laser_on_s | max | 790.012952 | 12 | 21 | 94199.6003 |

| Encounter | sc_a | sc_b | Start (s) | End (s) | Geometry (s) | Laser on (s) | Start clipped | End clipped |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 21 | 0 | 598.328079 | 598.328079 | 0 | true | false |
| 21 | 1 | 21 | 210230.365 | 211414.822 | 1184.45728 | 593.212941 | false | false |
| 20 | 2 | 21 | 199706.992 | 200848.975 | 1141.98306 | 565.150679 | false | false |
| 19 | 3 | 21 | 189094.855 | 190378.399 | 1283.54423 | 629.979344 | false | false |
| 18 | 4 | 21 | 178671.224 | 179825.612 | 1154.38799 | 610.829298 | false | false |
| 17 | 5 | 21 | 168109.332 | 169160.269 | 1050.93675 | 497.433856 | false | false |
| 16 | 6 | 21 | 157413.52 | 158840.472 | 1426.95159 | 704.990033 | false | false |
| 15 | 7 | 21 | 147134.637 | 148176.749 | 1042.11174 | 566.423756 | false | false |
| 14 | 8 | 21 | 136459.597 | 137509.311 | 1049.71401 | 472.358645 | false | false |
| 13 | 9 | 21 | 125778.345 | 127281.656 | 1503.31117 | 776.247472 | false | false |
| 12 | 10 | 21 | 115564.051 | 116499.067 | 935.015356 | 503.646988 | false | false |
| 11 | 11 | 21 | 104783.871 | 105914.368 | 1130.49633 | 493.484858 | false | false |
| 31 | 11 | 21 | 315591.934 | 316871.316 | 1279.38229 | 721.3588 | false | false |
| 10 | 12 | 21 | 94199.6003 | 95662.774 | 1463.1737 | 790.012952 | false | false |
| 30 | 12 | 21 | 305154.509 | 306032.695 | 878.186684 | 426.208592 | false | false |
| 9 | 13 | 21 | 83930.8215 | 84832.4522 | 901.630673 | 469.558812 | false | false |
| 29 | 13 | 21 | 294334.742 | 295762.519 | 1427.77756 | 652.310702 | false | false |
| 8 | 14 | 21 | 73118.3624 | 74357.4534 | 1239.09105 | 548.443945 | false | false |
| 28 | 14 | 21 | 284048.975 | 285214.611 | 1165.63533 | 649.863798 | false | false |
| 7 | 15 | 21 | 62642.2828 | 63988.4884 | 1346.20554 | 733.150994 | false | false |
| 27 | 15 | 21 | 273492.848 | 274462.418 | 969.570594 | 457.60073 | false | false |
| 6 | 16 | 21 | 52244.279 | 53207.2151 | 962.936077 | 483.724806 | false | false |
| 26 | 16 | 21 | 262770.552 | 264181.106 | 1410.55429 | 683.918316 | false | false |
| 5 | 17 | 21 | 41495.3043 | 42784.3623 | 1289.05793 | 601.229866 | false | false |
| 25 | 17 | 21 | 252442.519 | 253558.478 | 1115.9581 | 596.702633 | false | false |
| 4 | 18 | 21 | 31050.9417 | 32288.8287 | 1237.88707 | 654.196791 | false | false |
| 24 | 18 | 21 | 241835.508 | 242936.835 | 1101.32716 | 525.750996 | false | false |
| 3 | 19 | 21 | 20536.9706 | 21630.5441 | 1093.57349 | 541.555699 | false | false |
| 23 | 19 | 21 | 231242.377 | 232530.102 | 1287.72483 | 644.514398 | false | false |
| 2 | 20 | 21 | 9922.53602 | 11151.5475 | 1229.01145 | 606.209567 | false | false |
| 22 | 20 | 21 | 220777.8 | 221937.902 | 1160.10237 | 589.481959 | false | false |

### SpaceAGORA

Feather: [2_SpaceAGORA.jl/output/feather/spaceagora_N20_h1000km_t1150km_ih0_it0_e0_nu0_T325021.624s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s](../../2_SpaceAGORA.jl/output/feather/spaceagora_N20_h1000km_t1150km_ih0_it0_e0_nu0_T325021.624s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s)/

Analysis CSVs: [2_SpaceAGORA.jl/output/CSV/spaceagora_N20_h1000km_t1150km_ih0_it0_e0_nu0_T325021.624s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis](../../2_SpaceAGORA.jl/output/CSV/spaceagora_N20_h1000km_t1150km_ih0_it0_e0_nu0_T325021.624s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/analysis)/

Geometric encounters: 31; saved encounter-state rows: 3666.

| Metric | Extremum | Value | sc_a | sc_b | Time (s) |
| --- | --- | --- | --- | --- | --- |
| relative_speed_mps | min | 57.9085285 | 9 | 21 | 126530 |
| relative_speed_mps | max | 166.663715 | 9 | 21 | 127282.775 |
| range_rate_mps | min | -158.218011 | 13 | 21 | 294340.752 |
| range_rate_mps | max | 158.718375 | 11 | 21 | 316877.732 |
| geometry_duration_s | min | 598.328032 | 1 | 21 | 0 |
| geometry_duration_s | max | 1503.27233 | 9 | 21 | 125779.503 |
| laser_on_s | min | 0.000379794081 | 1 | 21 | 0 |
| laser_on_s | max | 788.316457 | 12 | 21 | 94200.2659 |

| Encounter | sc_a | sc_b | Start (s) | End (s) | Geometry (s) | Laser on (s) | Start clipped | End clipped |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 21 | 0 | 598.328032 | 598.328032 | 0.000379794081 | true | false |
| 21 | 1 | 21 | 210233.397 | 211417.924 | 1184.52699 | 587.090119 | false | false |
| 20 | 2 | 21 | 199709.643 | 200851.78 | 1142.1368 | 561.12626 | false | false |
| 19 | 3 | 21 | 189097.286 | 190381.029 | 1283.74291 | 628.737575 | false | false |
| 18 | 4 | 21 | 178673.535 | 179827.672 | 1154.13624 | 608.067821 | false | false |
| 17 | 5 | 21 | 168111.045 | 169162.298 | 1051.25324 | 494.703255 | false | false |
| 16 | 6 | 21 | 157415.244 | 158842.319 | 1427.07573 | 704.227073 | false | false |
| 15 | 7 | 21 | 147136.181 | 148178.043 | 1041.86227 | 560.956897 | false | false |
| 14 | 8 | 21 | 136460.64 | 137510.715 | 1050.07533 | 470.923447 | false | false |
| 13 | 9 | 21 | 125779.503 | 127282.775 | 1503.27233 | 775.832311 | false | false |
| 12 | 10 | 21 | 115564.939 | 116499.842 | 934.903122 | 498.586082 | false | false |
| 11 | 11 | 21 | 104784.465 | 105915.236 | 1130.77159 | 490.904325 | false | false |
| 31 | 11 | 21 | 315600.654 | 316877.732 | 1277.07845 | 712.702118 | false | false |
| 10 | 12 | 21 | 94200.2659 | 95663.3551 | 1463.08918 | 788.316457 | false | false |
| 30 | 12 | 21 | 305160.591 | 306039.367 | 878.775874 | 427.730966 | false | false |
| 9 | 13 | 21 | 83931.2478 | 84832.8686 | 901.620758 | 465.781207 | false | false |
| 29 | 13 | 21 | 294340.752 | 295770.142 | 1429.3904 | 648.79519 | false | false |
| 8 | 14 | 21 | 73118.6509 | 74357.8743 | 1239.22343 | 546.487038 | false | false |
| 28 | 14 | 21 | 284055.823 | 285219.843 | 1164.02029 | 645.693767 | false | false |
| 7 | 15 | 21 | 62642.5587 | 63988.7324 | 1346.17366 | 726.366625 | false | false |
| 27 | 15 | 21 | 273497.696 | 274468.104 | 970.40833 | 456.550287 | false | false |
| 6 | 16 | 21 | 52244.4283 | 53207.382 | 962.953685 | 480.017738 | false | false |
| 26 | 16 | 21 | 262775.655 | 264186.739 | 1411.08416 | 683.243541 | false | false |
| 5 | 17 | 21 | 41495.3912 | 42784.4849 | 1289.09366 | 597.236472 | false | false |
| 25 | 17 | 21 | 252447.589 | 253562.843 | 1115.25395 | 592.166702 | false | false |
| 4 | 18 | 21 | 31050.9897 | 32288.8831 | 1237.89339 | 650.581203 | false | false |
| 24 | 18 | 21 | 241839.584 | 242941.515 | 1101.9314 | 522.267808 | false | false |
| 3 | 19 | 21 | 20536.9841 | 21630.5654 | 1093.58134 | 535.14716 | false | false |
| 23 | 19 | 21 | 231246.483 | 232534.299 | 1287.81537 | 641.976334 | false | false |
| 2 | 20 | 21 | 9922.53449 | 11151.5491 | 1229.01458 | 603.203256 | false | false |
| 22 | 20 | 21 | 220781.555 | 221941.611 | 1160.05593 | 588.585789 | false | false |

## 3. Comparison

Differences are SpaceAGORA minus prototype. Percent differences use the absolute prototype value; zero or unavailable references have no percentage.

| Metric | Extremum | Prototype | SpaceAGORA | Difference | Difference (%) |
| --- | --- | --- | --- | --- | --- |
| geometry_duration_s | max | 1503.31117 | 1503.27233 | -0.0388357969 | -0.00258335053 |
| geometry_duration_s | min | 598.328079 | 598.328032 | -4.65623691e-05 | -7.7820799e-06 |
| laser_on_s | max | 790.012952 | 788.316457 | -1.6964948 | -0.214742657 |
| laser_on_s | min | 0 | 0.000379794081 | 0.000379794081 | unavailable |
| range_rate_mps | max | 158.758846 | 158.718375 | -0.0404703695 | -0.0254917257 |
| range_rate_mps | min | -158.25261 | -158.218011 | 0.0345992923 | 0.0218633312 |
| relative_speed_mps | max | 166.660069 | 166.663715 | 0.00364545751 | 0.0021873611 |
| relative_speed_mps | min | 57.9065854 | 57.9085285 | 0.00194310343 | 0.00335558283 |

Encounters are matched by unordered spacecraft pair and chronological occurrence, not by run-local IDs. Unmatched occurrences show unavailable differences; a missed encounter can shift subsequent occurrence matching.

| sc_a | sc_b | Occurrence | Entry difference (s) | Exit difference (s) | Geometry difference (s) | Laser-on difference (s) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 21 | 1 | 0 | -4.65623691e-05 | -4.65623691e-05 | 0.000379794081 |
| 1 | 21 | 2 | 3.03150955 | 3.10122021 | 0.0697106558 | -6.12282177 |
| 2 | 21 | 1 | 2.65112776 | 2.80487097 | 0.153743209 | -4.02441968 |
| 3 | 21 | 1 | 2.43146367 | 2.63014077 | 0.198677099 | -1.24176885 |
| 4 | 21 | 1 | 2.31089012 | 2.05914062 | -0.2517495 | -2.76147729 |
| 5 | 21 | 1 | 1.71272244 | 2.02921486 | 0.316492421 | -2.73060084 |
| 6 | 21 | 1 | 1.72350623 | 1.8476479 | 0.124141668 | -0.762959263 |
| 7 | 21 | 1 | 1.54379904 | 1.29433203 | -0.249467013 | -5.46685942 |
| 8 | 21 | 1 | 1.04259149 | 1.40391435 | 0.361322861 | -1.43519743 |
| 9 | 21 | 1 | 1.15761958 | 1.11878378 | -0.0388357969 | -0.415160982 |
| 10 | 21 | 1 | 0.887894185 | 0.77565971 | -0.112234475 | -5.06090608 |
| 11 | 21 | 1 | 0.593584171 | 0.868837993 | 0.275253822 | -2.58053295 |
| 11 | 21 | 2 | 8.72020535 | 6.41636259 | -2.30384276 | -8.65668229 |
| 12 | 21 | 1 | 0.665663469 | 0.581149649 | -0.0845138207 | -1.6964948 |
| 12 | 21 | 2 | 6.08257196 | 6.67176168 | 0.589189713 | 1.52237375 |
| 13 | 21 | 1 | 0.426301907 | 0.416387586 | -0.00991432052 | -3.77760485 |
| 13 | 21 | 2 | 6.00993862 | 7.62277216 | 1.61283354 | -3.51551234 |
| 14 | 21 | 1 | 0.288489763 | 0.420869113 | 0.13237935 | -1.95690663 |
| 14 | 21 | 2 | 6.84756499 | 5.23253308 | -1.61503191 | -4.17003088 |
| 15 | 21 | 1 | 0.275901595 | 0.244020017 | -0.0318815779 | -6.7843699 |
| 15 | 21 | 2 | 4.84829086 | 5.68602774 | 0.837736881 | -1.05044271 |
| 16 | 21 | 1 | 0.149297948 | 0.166905702 | 0.0176077539 | -3.70706784 |
| 16 | 21 | 2 | 5.10238219 | 5.63224487 | 0.529862677 | -0.674774886 |
| 17 | 21 | 1 | 0.0868835677 | 0.122622851 | 0.035739283 | -3.99339388 |
| 17 | 21 | 2 | 5.06997366 | 4.36583082 | -0.70414284 | -4.53593089 |
| 18 | 21 | 1 | 0.048057521 | 0.0543769548 | 0.00631943381 | -3.615588 |
| 18 | 21 | 2 | 4.0753947 | 4.67963885 | 0.604244151 | -3.48318757 |
| 19 | 21 | 1 | 0.0134086218 | 0.0212639251 | 0.00785530329 | -6.40853988 |
| 19 | 21 | 2 | 4.10612664 | 4.19666154 | 0.0905348996 | -2.53806425 |
| 20 | 21 | 1 | -0.00153892451 | 0.00158814297 | 0.00312706747 | -3.00631037 |
| 20 | 21 | 2 | 3.75533426 | 3.70888784 | -0.0464464221 | -0.896170161 |

No numerical acceptance tolerance was specified: these tables quantify agreement, not a physical-equivalence verdict. Motion extrema use only encounter states and entry/exit boundaries; unrelated pairs and out-of-encounter times are excluded. Clipped encounters cover only the recorded portion. Earth/light-speed constants, solvers, LOS checks, and scheduling differ. SpaceAGORA uses accepted-step scheduling and a velocity-kick callback; prototype forces are evaluated in the ODE. The report reads saved analysis CSVs after both worker processes exit; it does not compare full trajectories or equate the two delta-v diagnostics.

