# Local native-coordinate comparison results

Times are recorded measurements from the validation machine; see docs/ks_cartesian_comparison.md for methodology and limitations.

## Fastest passing nonlinear settings tested

Limits: 0.01 m maximum position error and 0.00001 m/s maximum velocity error, adding the measured truth-refinement difference as a safety margin (not a rigorous bound on reference error). This selects among the tested settings only; it is not a global matched-accuracy search. One warmed timing sample per configuration.

| e | Orbits | Route | Solver | Mode | Setting | Solve (s) | Total (s) | Max position (m) | RMS position (m) | Max velocity (m/s) |
| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0.01 | 20 | cartesian | DP5 | adaptive | 1e-12 | 0.04185 | 0.0503 | 0.009841 | 0.004412 | 9.912e-06 |
| 0.01 | 20 | ks_target_s | DP5 | fixed | 20 | 0.01467 | 0.03059 | 0.001621 | 0.0007138 | 1.633e-06 |
| 0.3 | 20 | cartesian | Tsit5 | adaptive | 1e-12 | 0.06394 | 0.07266 | 0.006008 | 0.001989 | 5.34e-06 |
| 0.3 | 20 | ks_target_s | DP5 | adaptive | 1e-12 | 0.01493 | 0.0327 | 0.005736 | 0.001879 | 5.085e-06 |

## 20-orbit scheduling results

Both nonlinear routes replay the KS-selected schedule. Changes subtract each route's own no-laser trajectory.

| Objective | Element | KS linear | Cartesian linear | KS nonlinear | Cartesian nonlinear | Different intervals |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| a | a_m | 3252.046577 | 3252.046577 | -688.0049149 | -687.9771369 | 0 |
| a | e | 5.792231955e-05 | 5.792231955e-05 | 1.78065142e-05 | 1.780767238e-05 | 0 |
| a | i_deg | 0.001163149556 | 0.001163149556 | 0.0002058452622 | 0.0002058570854 | 0 |
| a | RAAN_deg | 0.002996718661 | 0.002996718661 | 0.0005766448015 | 0.0005766746757 | 0 |
| a | omega_deg | -0.04402916048 | -0.04402916049 | 0.0317459706 | 0.03176250657 | 0 |
| i | a_m | 230.9194961 | 230.9194962 | -15.44076079 | -15.42919528 | 0 |
| i | e | 1.394709578e-05 | 1.394709578e-05 | 4.754745129e-06 | 4.754340062e-06 | 0 |
| i | i_deg | 0.007040463162 | 0.007040463162 | 0.006839618558 | 0.006839629311 | 0 |
| i | RAAN_deg | 0.0008344301981 | 0.0008344301981 | 0.0004800976471 | 0.0004799739345 | 0 |
| i | omega_deg | 0.009825133874 | 0.009825133873 | 0.003357976691 | 0.003358919691 | 0 |
| omega | a_m | -29.05390291 | -29.05390291 | 41.34367446 | 41.25454279 | 0 |
| omega | e | -8.302394976e-07 | -8.302394974e-07 | -7.383936168e-06 | -7.383202177e-06 | 0 |
| omega | i_deg | 6.814809743e-06 | 6.814809743e-06 | 1.049313334e-06 | 1.063084449e-06 | 0 |
| omega | RAAN_deg | 3.381498065e-05 | 3.381498065e-05 | 2.865592611e-05 | 2.870377966e-05 | 0 |
| omega | omega_deg | 1.434039317 | 1.434039317 | 1.065777538 | 1.065682089 | 0 |
| combined | a_m | 3247.329879 | 3247.32988 | -106.4313265 | -106.4333804 | 0 |
| combined | e | 5.754587559e-05 | 5.754587559e-05 | 1.134560501e-05 | 1.134687283e-05 | 0 |
| combined | i_deg | 0.001159396731 | 0.001159396731 | 0.0004099372356 | 0.0004099377588 | 0 |
| combined | RAAN_deg | 0.002982370062 | 0.002982370063 | 0.001099238002 | 0.001099179689 | 0 |
| combined | omega_deg | -0.03103927974 | -0.03103927974 | -0.03004479742 | -0.03004388506 | 0 |

## KS versus Cartesian terminal-gain scheduling: runtime and error

Wall time of one warmed call per objective, GC in parentheses. Planning is `plan_laser_schedule` with `linearization=:ks` or `:cartesian` (GC forced beforehand); replay integrates the KS-selected schedule nonlinearly (no forced GC). Ratio is Cartesian / KS. Errors are terminal delta a: linear error is each route's prediction minus its own nonlinear replay; cross-route columns are KS minus Cartesian. Other elements are in the next table.

| Objective | KS planning (s) | Cartesian planning (s) | Ratio | KS linear Δa error (m) | Cartesian linear Δa error (m) | Linear KS − Cartesian Δa (m) | KS replay (s) | Cartesian replay (s) | Ratio | Nonlinear KS − Cartesian Δa (m) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| a | 1.129 (0.182) | 1.153 (0.183) | 1.02 | 3940 | 3940 | -1.93e-07 | 2.031 (0.82) | 1.998 (0.765) | 0.984 | -0.0278 |
| i | 1.146 (0.182) | 1.231 (0.201) | 1.07 | 246.4 | 246.3 | -3.25e-08 | 2.155 (0.854) | 2.246 (0.917) | 1.04 | -0.0116 |
| omega | 1.134 (0.175) | 1.165 (0.18) | 1.03 | -70.4 | -70.31 | -1.94e-09 | 2.327 (0.9) | 2.305 (0.915) | 0.991 | 0.0891 |
| combined | 1.145 (0.181) | 1.174 (0.185) | 1.03 | 3354 | 3354 | -1.92e-07 | 2.197 (0.922) | 1.762 (0.448) | 0.802 | 0.00205 |

## KS versus Cartesian terminal-gain scheduling: errors by element

Linear-prediction error is each route's terminal linear prediction minus its own nonlinear replay of the KS-selected schedule (both routes selected identical schedules, see Different intervals above). Cross-route columns are KS minus Cartesian. Max gain difference is the largest terminal-gain mismatch between the two linearizations.

| Objective | Element | KS linear error | Cartesian linear error | Linear KS − Cartesian | Nonlinear KS − Cartesian | Max gain difference |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| a | a_m | 3940 | 3940 | -1.93e-07 | -0.0278 | 1.03e-10 |
| a | e | 4.012e-05 | 4.011e-05 | -1.64e-15 | -1.16e-09 | 1.02e-18 |
| a | i_deg | 0.0009573 | 0.0009573 | -2.23e-14 | -1.18e-08 | 1.46e-17 |
| a | RAAN_deg | 0.00242 | 0.00242 | -6.82e-14 | -2.99e-08 | 4.24e-17 |
| a | omega_deg | -0.07578 | -0.07579 | 5.29e-12 | -1.65e-05 | 2.75e-15 |
| i | a_m | 246.4 | 246.3 | -3.25e-08 | -0.0116 | 1.03e-10 |
| i | e | 9.192e-06 | 9.193e-06 | -3.3e-16 | 4.05e-10 | 1.02e-18 |
| i | i_deg | 0.0002008 | 0.0002008 | -4.31e-15 | -1.08e-08 | 1.46e-17 |
| i | RAAN_deg | 0.0003543 | 0.0003545 | -1.08e-14 | 1.24e-07 | 4.24e-17 |
| i | omega_deg | 0.006467 | 0.006466 | 7.08e-13 | -9.43e-07 | 2.75e-15 |
| omega | a_m | -70.4 | -70.31 | -1.94e-09 | 0.0891 | 1.03e-10 |
| omega | e | 6.554e-06 | 6.553e-06 | -1.83e-16 | -7.34e-10 | 1.02e-18 |
| omega | i_deg | 5.765e-06 | 5.752e-06 | 9.71e-17 | -1.38e-08 | 1.46e-17 |
| omega | RAAN_deg | 5.159e-06 | 5.111e-06 | -1.57e-15 | -4.79e-08 | 4.24e-17 |
| omega | omega_deg | 0.3683 | 0.3684 | 7.59e-14 | 9.54e-05 | 2.75e-15 |
| combined | a_m | 3354 | 3354 | -1.92e-07 | 0.00205 | 1.03e-10 |
| combined | e | 4.62e-05 | 4.62e-05 | -1.74e-15 | -1.27e-09 | 1.02e-18 |
| combined | i_deg | 0.0007495 | 0.0007495 | -2.2e-14 | -5.23e-10 | 1.46e-17 |
| combined | RAAN_deg | 0.001883 | 0.001883 | -6.84e-14 | 5.83e-08 | 4.24e-17 |
| combined | omega_deg | -0.0009945 | -0.0009954 | 5.86e-12 | -9.12e-07 | 2.75e-15 |

## KS versus Cartesian receding-horizon MPC: runtime and error

`run_laser_mpc` over 20 orbits, horizon 1200 intervals, 240 executed per replan (20 replans). Runtime is one warmed call (GC in parentheses) and includes the per-segment nonlinear KS replays. Segment error is the largest |linear − nonlinear KS| delta a over segments. Terminal delta a replays each route's executed schedule continuously in KS and Cartesian against that coordinate's no-laser run. Ratio is Cartesian / KS.

| Objective | KS MPC (s) | Cartesian MPC (s) | Ratio | KS max segment Δa error (m) | Cartesian max segment Δa error (m) | KS schedule: Δa KS / Cartesian replay (m) | Cartesian schedule: Δa KS / Cartesian replay (m) | Different intervals |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| a | 12.79 (5.19) | 12.84 (5.24) | 1 | 0.614 | 0.614 | 1987.45 / 1987.44 | 1987.45 / 1987.44 | 0 |
| i | 12.99 (4.99) | 13.61 (5.45) | 1.05 | 0.3576 | 0.3576 | -53.3492 / -53.3257 | -53.3492 / -53.3257 | 0 |
| omega | 13.06 (4.98) | 13.16 (5.01) | 1.01 | 0.4277 | 0.4277 | -19.0729 / -19.0887 | -19.0729 / -19.0887 | 0 |
| combined | 11.99 (4.43) | 12.09 (4.44) | 1.01 | 1.091 | 1.091 | 2011.59 / 2011.59 | 2011.59 / 2011.59 | 0 |

## Interpretation

The scheduler uses the original frozen-interval physical-time model in RTN. Its KS route maps SpaceAGORA's analytical KS tangent to Cartesian before RTN/ZOH. It is not a native fictitious-time STM scheduler.
The nonlinear benchmark is inertial: no RTN conversion, scheduling, or Jacobian assembly enters its solve timing. KS uses target fictitious time; fixed KS steps are nominal physical-second equivalents.
Equal numerical tolerances on differently scaled states are not a matched-accuracy result. Use work-precision plots and truth_refinement.csv, not one runtime alone, to judge tradeoffs.
A separate inertial_zoh.csv per objective propagates the selected schedule with native-s KS and inertial Cartesian matrices, with no RTN conversion. KS applies a fixed-physical-time clock correction after each interval. These finite-step approximations need not coincide with the RTN ZOH or with each other; they are not used to replace the original optimizer.
Nonlinear horizons completed: 20 orbits; 48 configurations. The scheduling and native inertial linear comparisons are 20-orbit runs. Historical affine/GVE/PPTX experiments have not been re-created by this driver.
