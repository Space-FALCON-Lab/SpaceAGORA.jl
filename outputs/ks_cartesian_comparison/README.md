# Local native-coordinate comparison results

Times are recorded measurements from the validation machine; see docs/ks_cartesian_comparison.md for methodology and limitations.

## Fastest passing nonlinear settings tested

Limits: 0.01 m maximum position error and 0.00001 m/s maximum velocity error, adding the measured truth-refinement difference as a safety margin (not a rigorous bound on reference error). This selects among the tested settings only; it is not a global matched-accuracy search. One warmed timing sample per configuration.

| e | Orbits | Route | Solver | Mode | Setting | Solve (s) | Total (s) | Max position (m) | RMS position (m) | Max velocity (m/s) |
| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0.01 | 1 | cartesian | DP5 | adaptive | 1e-10 | 0.0008822 | 0.001182 | 0.001924 | 0.0008215 | 1.862e-06 |
| 0.01 | 1 | ks_target_s | Tsit5 | adaptive | 1e-10 | 0.0005076 | 0.001143 | 7.547e-05 | 2.877e-05 | 7.494e-08 |
| 0.01 | 5 | cartesian | DP5 | fixed | 20 | 0.005059 | 0.006449 | 0.009278 | 0.00399 | 9.333e-06 |
| 0.01 | 5 | ks_target_s | Tsit5 | adaptive | 1e-10 | 0.001757 | 0.00503 | 0.001952 | 0.000851 | 1.96e-06 |
| 0.01 | 50 | cartesian | Tsit5 | fixed | 5 | 0.328 | 0.3494 | 0.001113 | 0.0004344 | 1.12e-06 |
| 0.01 | 50 | ks_target_s | Tsit5 | fixed | 20 | 0.04451 | 0.08981 | 0.001591 | 0.0006992 | 1.601e-06 |
| 0.01 | 100 | cartesian | DP5 | fixed | 5 | 0.7138 | 0.7575 | 0.003893 | 0.001635 | 3.917e-06 |
| 0.01 | 100 | ks_target_s | Tsit5 | fixed | 20 | 0.0924 | 0.1873 | 0.006362 | 0.00279 | 6.403e-06 |
| 0.01 | 200 | cartesian | No passing setting | | | | | | | |
| 0.01 | 200 | ks_target_s | DP5 | fixed | 5 | 0.8443 | 1.066 | 0.0006492 | 0.000241 | 6.512e-07 |
| 0.3 | 1 | cartesian | DP5 | adaptive | 1e-10 | 0.0008848 | 0.001181 | 0.007535 | 0.002623 | 6.3e-06 |
| 0.3 | 1 | ks_target_s | DP5 | fixed | 100 | 0.0004415 | 0.001162 | 0.002574 | 0.0009215 | 2.34e-06 |
| 0.3 | 5 | cartesian | DP5 | adaptive | 1e-12 | 0.00981 | 0.01131 | 0.00125 | 0.0004007 | 1.101e-06 |
| 0.3 | 5 | ks_target_s | Tsit5 | adaptive | 1e-10 | 0.001808 | 0.005115 | 0.00425 | 0.001345 | 3.752e-06 |
| 0.3 | 50 | cartesian | DP5 | fixed | 5 | 1.01 | 1.027 | 0.001498 | 0.0004352 | 1.329e-06 |
| 0.3 | 50 | ks_target_s | Tsit5 | adaptive | 1e-12 | 0.05376 | 0.09553 | 0.003773 | 0.001219 | 3.339e-06 |
| 0.3 | 100 | cartesian | DP5 | fixed | 5 | 1.704 | 1.743 | 0.005959 | 0.001739 | 5.292e-06 |
| 0.3 | 100 | ks_target_s | DP5 | fixed | 20 | 0.1427 | 0.2398 | 0.005851 | 0.001847 | 5.192e-06 |
| 0.3 | 200 | cartesian | No passing setting | | | | | | | |
| 0.3 | 200 | ks_target_s | Tsit5 | fixed | 20 | 0.3379 | 0.5441 | 0.004042 | 0.0009968 | 3.585e-06 |

## 5-orbit scheduling results

Both nonlinear routes replay the KS-selected schedule. Changes subtract each route's own no-laser trajectory.

| Objective | Element | KS linear | Cartesian linear | KS nonlinear | Cartesian nonlinear | Different intervals |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| a | a_m | 728.7011883 | 728.7011883 | 728.5822958 | 728.5823006 | 0 |
| a | e | 2.999611115e-06 | 2.999611115e-06 | 2.921368954e-06 | 2.921367437e-06 | 0 |
| a | i_deg | 6.603340187e-05 | 6.603340187e-05 | 6.561435591e-05 | 6.561436258e-05 | 0 |
| a | RAAN_deg | 0.0001687100256 | 0.0001687100256 | 0.0001687459159 | 0.000168745201 | 0 |
| a | omega_deg | -0.006920925898 | -0.006920925898 | -0.006922959127 | -0.00692295206 | 0 |
| i | a_m | 150.8409299 | 150.8409299 | 90.11841122 | 90.12660272 | 0 |
| i | e | 2.530615758e-06 | 2.530615758e-06 | 2.899587267e-06 | 2.899452891e-06 | 0 |
| i | i_deg | 0.001733471568 | 0.001733471568 | 0.001738960051 | 0.001738960012 | 0 |
| i | RAAN_deg | 7.50242199e-05 | 7.50242199e-05 | 7.200918149e-05 | 7.199644283e-05 | 0 |
| i | omega_deg | 0.00241993295 | 0.00241993295 | 0.01037283411 | 0.01037212424 | 0 |
| omega | a_m | -4.69777913 | -4.69777913 | -5.339423141 | -5.333454038 | 0 |
| omega | e | -1.673790324e-06 | -1.673790324e-06 | -1.50563084e-06 | -1.505290771e-06 | 0 |
| omega | i_deg | 4.182641058e-06 | 4.182641058e-06 | 4.119498627e-06 | 4.134254594e-06 | 0 |
| omega | RAAN_deg | 1.543012635e-05 | 1.543012635e-05 | 1.524181217e-05 | 1.526897783e-05 | 0 |
| omega | omega_deg | 0.3648361213 | 0.3648361213 | 0.3651007654 | 0.365100815 | 0 |
| combined | a_m | 724.2304561 | 724.2304561 | 723.0780297 | 723.0780613 | 0 |
| combined | e | 2.91068455e-06 | 2.91068455e-06 | 2.872501122e-06 | 2.872551531e-06 | 0 |
| combined | i_deg | 6.494782863e-05 | 6.494782863e-05 | 6.463448623e-05 | 6.463248401e-05 | 0 |
| combined | RAAN_deg | 0.0001625288739 | 0.0001625288739 | 0.0001618675854 | 0.0001618574152 | 0 |
| combined | omega_deg | 0.006079451624 | 0.006079451624 | 0.007727966592 | 0.007727931501 | 0 |

## KS versus Cartesian terminal-gain scheduling: runtime and error

Wall time of one warmed call per objective, GC in parentheses. Planning is `plan_laser_schedule` with `linearization=:ks` or `:cartesian` (GC forced beforehand); replay integrates the KS-selected schedule nonlinearly (no forced GC). Ratio is Cartesian / KS. Errors are terminal delta a: linear error is each route's prediction minus its own nonlinear replay; cross-route columns are KS minus Cartesian. Other elements are in the next table.

| Objective | KS planning (s) | Cartesian planning (s) | Ratio | KS linear Δa error (m) | Cartesian linear Δa error (m) | Linear KS − Cartesian Δa (m) | KS replay (s) | Cartesian replay (s) | Ratio | Nonlinear KS − Cartesian Δa (m) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| a | 0.2939 (0.0407) | 0.2757 (0.0304) | 0.938 | 0.1189 | 0.1189 | 4.37e-11 | 0.3999 (0.12) | 0.7832 (0.477) | 1.96 | -4.77e-06 |
| i | 0.2742 (0.0229) | 0.278 (0.0225) | 1.01 | 60.72 | 60.71 | 8.19e-12 | 0.4017 (0.0986) | 0.7671 (0.479) | 1.91 | -0.00819 |
| omega | 0.2774 (0.023) | 0.2814 (0.0258) | 1.01 | 0.6416 | 0.6357 | 4.32e-12 | 0.48 (0.112) | 0.819 (0.486) | 1.71 | -0.00597 |
| combined | 0.2742 (0.0246) | 0.2802 (0.0247) | 1.02 | 1.152 | 1.152 | 3.09e-11 | 0.4348 (0.0939) | 0.8223 (0.519) | 1.89 | -3.16e-05 |

## KS versus Cartesian terminal-gain scheduling: errors by element

Linear-prediction error is each route's terminal linear prediction minus its own nonlinear replay of the KS-selected schedule (both routes selected identical schedules, see Different intervals above). Cross-route columns are KS minus Cartesian. Max gain difference is the largest terminal-gain mismatch between the two linearizations.

| Objective | Element | KS linear error | Cartesian linear error | Linear KS − Cartesian | Nonlinear KS − Cartesian | Max gain difference |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| a | a_m | 0.1189 | 0.1189 | 4.37e-11 | -4.77e-06 | 6.88e-14 |
| a | e | 7.824e-08 | 7.824e-08 | 1.52e-18 | 1.52e-12 | 1.43e-20 |
| a | i_deg | 4.19e-07 | 4.19e-07 | -1.15e-16 | -6.67e-12 | 1.67e-19 |
| a | RAAN_deg | -3.589e-08 | -3.518e-08 | -3.73e-17 | 7.15e-10 | 8.23e-20 |
| a | omega_deg | 2.033e-06 | 2.026e-06 | -6.58e-14 | -7.07e-09 | 1.26e-16 |
| i | a_m | 60.72 | 60.71 | 8.19e-12 | -0.00819 | 6.88e-14 |
| i | e | -3.69e-07 | -3.688e-07 | 2.17e-19 | 1.34e-10 | 1.43e-20 |
| i | i_deg | -5.488e-06 | -5.488e-06 | -2.56e-17 | 3.92e-11 | 1.67e-19 |
| i | RAAN_deg | 3.015e-06 | 3.028e-06 | -5.43e-18 | 1.27e-08 | 8.23e-20 |
| i | omega_deg | -0.007953 | -0.007952 | -6.18e-15 | 7.1e-07 | 1.26e-16 |
| omega | a_m | 0.6416 | 0.6357 | 4.32e-12 | -0.00597 | 6.88e-14 |
| omega | e | -1.682e-07 | -1.685e-07 | 1.36e-19 | -3.4e-10 | 1.43e-20 |
| omega | i_deg | 6.314e-08 | 4.839e-08 | -1.09e-17 | -1.48e-08 | 1.67e-19 |
| omega | RAAN_deg | 1.883e-07 | 1.611e-07 | -2.33e-18 | -2.72e-08 | 8.23e-20 |
| omega | omega_deg | -0.0002646 | -0.0002647 | -6.61e-15 | -4.96e-08 | 1.26e-16 |
| combined | a_m | 1.152 | 1.152 | 3.09e-11 | -3.16e-05 | 6.88e-14 |
| combined | e | 3.818e-08 | 3.813e-08 | -1.52e-18 | -5.04e-11 | 1.43e-20 |
| combined | i_deg | 3.133e-07 | 3.153e-07 | -1.34e-16 | 2e-09 | 1.67e-19 |
| combined | RAAN_deg | 6.613e-07 | 6.715e-07 | -4.04e-17 | 1.02e-08 | 8.23e-20 |
| combined | omega_deg | -0.001649 | -0.001648 | -9.98e-14 | 3.51e-08 | 1.26e-16 |

## KS versus Cartesian receding-horizon MPC: runtime and error

`run_laser_mpc` over 5 orbits, horizon 1200 intervals, 240 executed per replan (5 replans). Runtime is one warmed call (GC in parentheses) and includes the per-segment nonlinear KS replays. Segment error is the largest |linear − nonlinear KS| delta a over segments. Terminal delta a replays each route's executed schedule continuously in KS and Cartesian against that coordinate's no-laser run. Ratio is Cartesian / KS.

| Objective | KS MPC (s) | Cartesian MPC (s) | Ratio | KS max segment Δa error (m) | Cartesian max segment Δa error (m) | KS schedule: Δa KS / Cartesian replay (m) | Cartesian schedule: Δa KS / Cartesian replay (m) | Different intervals |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| a | 2.526 (1.2) | 2.649 (1.18) | 1.05 | 0.002174 | 0.002174 | 728.645 / 728.645 | 728.645 / 728.645 | 0 |
| i | 2.671 (1.21) | 2.703 (1.22) | 1.01 | 0.3576 | 0.3576 | 90.7519 / 90.7602 | 90.7519 / 90.7602 | 0 |
| omega | 2.738 (1.21) | 2.752 (1.2) | 1.01 | 0.1962 | 0.1962 | -4.05965 / -4.05377 | -4.05965 / -4.05377 | 0 |
| combined | 2.743 (1.24) | 2.69 (1.2) | 0.981 | 0.06446 | 0.06446 | 723.574 / 723.574 | 723.574 / 723.574 | 0 |

## Interpretation

The scheduler uses the original frozen-interval physical-time model in RTN. Its KS route maps SpaceAGORA's analytical KS tangent to Cartesian before RTN/ZOH. It is not a native fictitious-time STM scheduler.
The nonlinear benchmark is inertial: no RTN conversion, scheduling, or Jacobian assembly enters its solve timing. KS uses target fictitious time; fixed KS steps are nominal physical-second equivalents.
Equal numerical tolerances on differently scaled states are not a matched-accuracy result. Use work-precision plots and truth_refinement.csv, not one runtime alone, to judge tradeoffs.
A separate inertial_zoh.csv per objective propagates the selected schedule with native-s KS and inertial Cartesian matrices, with no RTN conversion. KS applies a fixed-physical-time clock correction after each interval. These finite-step approximations need not coincide with the RTN ZOH or with each other; they are not used to replace the original optimizer.
Nonlinear horizons completed: 1, 5, 50, 100, 200 orbits; 240 configurations. The scheduling and native inertial linear comparisons are 5-orbit runs. Historical affine/GVE/PPTX experiments have not been re-created by this driver.
