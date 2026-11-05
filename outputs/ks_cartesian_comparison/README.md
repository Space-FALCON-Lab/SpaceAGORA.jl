# Local native-coordinate comparison results

Times are recorded measurements from the validation machine; see docs/ks_cartesian_comparison.md for methodology and limitations.

## Fastest passing nonlinear settings tested

Limits: 0.01 m maximum position error and 0.00001 m/s maximum velocity error, adding the measured truth-refinement difference as a safety margin (not a rigorous bound on reference error). This selects among the tested settings only; it is not a global matched-accuracy search. One warmed timing sample per configuration.

| e | Orbits | Route | Solver | Mode | Setting | Solve (s) | Total (s) | Max position (m) |
| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: |
| 0.01 | 1 | cartesian | DP5 | adaptive | 1e-10 | 0.0008822 | 0.001182 | 0.001924 |
| 0.01 | 1 | ks_target_s | Tsit5 | adaptive | 1e-10 | 0.0005076 | 0.001143 | 7.547e-05 |
| 0.01 | 5 | cartesian | DP5 | fixed | 20 | 0.005059 | 0.006449 | 0.009278 |
| 0.01 | 5 | ks_target_s | Tsit5 | adaptive | 1e-10 | 0.001757 | 0.00503 | 0.001952 |
| 0.01 | 50 | cartesian | Tsit5 | fixed | 5 | 0.328 | 0.3494 | 0.001113 |
| 0.01 | 50 | ks_target_s | Tsit5 | fixed | 20 | 0.04451 | 0.08981 | 0.001591 |
| 0.01 | 100 | cartesian | DP5 | fixed | 5 | 0.7138 | 0.7575 | 0.003893 |
| 0.01 | 100 | ks_target_s | Tsit5 | fixed | 20 | 0.0924 | 0.1873 | 0.006362 |
| 0.01 | 200 | cartesian | No passing setting | | | | | |
| 0.01 | 200 | ks_target_s | DP5 | fixed | 5 | 0.8443 | 1.066 | 0.0006492 |
| 0.3 | 1 | cartesian | DP5 | adaptive | 1e-10 | 0.0008848 | 0.001181 | 0.007535 |
| 0.3 | 1 | ks_target_s | DP5 | fixed | 100 | 0.0004415 | 0.001162 | 0.002574 |
| 0.3 | 5 | cartesian | DP5 | adaptive | 1e-12 | 0.00981 | 0.01131 | 0.00125 |
| 0.3 | 5 | ks_target_s | Tsit5 | adaptive | 1e-10 | 0.001808 | 0.005115 | 0.00425 |
| 0.3 | 50 | cartesian | DP5 | fixed | 5 | 1.01 | 1.027 | 0.001498 |
| 0.3 | 50 | ks_target_s | Tsit5 | adaptive | 1e-12 | 0.05376 | 0.09553 | 0.003773 |
| 0.3 | 100 | cartesian | DP5 | fixed | 5 | 1.704 | 1.743 | 0.005959 |
| 0.3 | 100 | ks_target_s | DP5 | fixed | 20 | 0.1427 | 0.2398 | 0.005851 |
| 0.3 | 200 | cartesian | No passing setting | | | | | |
| 0.3 | 200 | ks_target_s | Tsit5 | fixed | 20 | 0.3379 | 0.5441 | 0.004042 |

## Five-orbit scheduling results

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

## Interpretation

The scheduler uses the original frozen-interval physical-time model in RTN. Its KS route maps SpaceAGORA's analytical KS tangent to Cartesian before RTN/ZOH. It is not a native fictitious-time STM scheduler.
The nonlinear benchmark is inertial: no RTN conversion, scheduling, or Jacobian assembly enters its solve timing. KS uses target fictitious time; fixed KS steps are nominal physical-second equivalents.
Equal numerical tolerances on differently scaled states are not a matched-accuracy result. Use work-precision plots and truth_refinement.csv, not one runtime alone, to judge tradeoffs.
A separate inertial_zoh.csv per objective propagates the selected schedule with native-s KS and inertial Cartesian matrices, with no RTN conversion. KS applies a fixed-physical-time clock correction after each interval. These finite-step approximations need not coincide with the RTN ZOH or with each other; they are not used to replace the original optimizer.
Nonlinear horizons completed: 1, 5, 50, 100, 200 orbits; 240 configurations. The scheduling and native inertial linear comparisons remain five-orbit runs. Historical affine/GVE/PPTX experiments have not been re-created by this driver.
