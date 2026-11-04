# Target-helper encounter and laser-on distributions

## Are these quantities collected?

**Yes, all four are collected**, with the following definitions:

| Quantity | Saved product | Definition |
|---|---|---|
| Encounter speed | `encounters.feather` | Relative ECI speed `norm(v_b - v_a)`, m/s; sampled minimum and maximum for each encounter. |
| Range rate | `encounters.feather` | Signed `dot(r_b-r_a, v_b-v_a) / norm(r_b-r_a)`, m/s; negative approaching, positive receding; sampled minimum and maximum. |
| Encounter duration | `encounters.feather` | End minus start time inside the 200 km geometric range, seconds. |
| Laser-on duration | `laser_firing.feather` | Duration of each uninterrupted selected-pair firing interval, seconds. Includes any firing outside 200 km while selection is held between updates. |
| Laser-on within an encounter | `firing_range_audit.feather` | Derived by summing saved within-range segment durations for each encounter; includes zero for encounters never fired upon. |

## Scope and interpretation

- All **270 verified cases**, each spanning **100 initial target periods**; target-helper pairs only. Helpers start at 1000 km and 0 deg. Here i and h are the **initial target** inclination and altitude, not their evolving values.
- Exact grid: i = 0, 0.5, 1 deg; h = 900, 950, 1000, 1050, 1100 km; helpers = 1, 100, 200, 300, 400, 500, 600, 700, 800; schedulers = `gve_sma`, `gve_inc`.
- Each table pools events over the other three dimensions. **Event-weighted, not equal-case-weighted or time-weighted**: cases producing more events contribute more. These are marginal descriptive comparisons, not isolated causal scheduler effects.
- Speed/range-rate extrema use **10-second samples plus exact boundary states**. They are sampled estimates, not guaranteed continuous-time extrema. Full speed/range-rate time histories for every encounter are not retained; only the longest-encounter tracks are saved. Do not interpret the extrema distributions as instantaneous-speed distributions.
- Main tables include boundary-clipped intervals as observed, without extrapolating full lifetimes. Both start- and end-clipping flags are exported; `unclipped_only` distributions are also in the CSV.
- Each geometric encounter contributes one duration, two speed extrema and two range-rate extrema. Each firing interval contributes one uninterrupted laser duration; it need not match an encounter. Accumulated in-encounter laser-on may combine multiple firing intervals. No-firing cases have no fictional zero-length firing interval; their case total is zero.
- SD is population SD; percentiles use Julia's default linearly interpolated quantile (Hyndman-Fan type 7). NA means no events, not zero. Numbers shown to six significant digits; CSVs retain full precision.
- All input settings, verified status, run IDs, and SHA-256 checksums of the three consumed Feather products were checked. This is an evaluation of saved results, not a new propagation or a rerun of trajectory physics verification.

## Campaign totals

- Target-helper encounters: **137564**, of which **1610** are boundary-clipped and **30739** have no within-encounter laser-on time.
- Uninterrupted laser-on intervals: **124675**, of which **320** are boundary-clipped.
- Cases without encounters: **0**; cases without firing: **45**.
- Summed firing across independent runs: **9.27549e+07 s**, including **120629 s** outside range. These are run-summed totals, not one mission's duration.
- Geometric pair-time can overlap among helpers and must not be interpreted as target access time or duty fraction. The case CSV's laser duty fraction is total actual target firing divided by that case's duration.

| Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Encounter speed: sampled minimum (m/s) | 137564 | 63.4005 | 33.8006 | 0.102125 | 20.8898 | 42.1481 | 57.5104 | 76.086 | 130.83 | 238.481 |
| Encounter speed: sampled maximum (m/s) | 137564 | 195.613 | 18.0595 | 3.67418 | 169.794 | 181.799 | 194.623 | 207.671 | 227.925 | 249.815 |
| Signed range rate: sampled minimum (m/s) | 137564 | -110.803 | 30.8296 | -176.038 | -151.337 | -133.742 | -116.169 | -89.7923 | -61.7155 | 128.321 |
| Signed range rate: sampled maximum (m/s) | 137564 | 110.839 | 31.1601 | -152.833 | 59.8507 | 89.3876 | 115.863 | 133.703 | 152.843 | 178.793 |
| Encounter duration (s) | 137564 | 3771.06 | 16105.3 | 2.33803 | 1630.31 | 2009.72 | 2369.26 | 4490.08 | 5604.54 | 630712 |
| Uninterrupted laser-on duration (s) | 124675 | 743.974 | 1265.49 | 0.37341 | 20.2384 | 332.627 | 536.478 | 906.862 | 1799.46 | 61346.1 |
| Total laser-on within an encounter (s) | 137564 | 673.391 | 3613.57 | 0 | 0 | 30.3171 | 460.305 | 790.747 | 1766.22 | 314471 |

## 1) Distribution by Target inclination i (deg)

### Coverage and case totals

| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 90 | 45118 | 560 | 22799 | 82 | 45 | 1297.84 |
| 0.5 | 90 | 45387 | 546 | 50552 | 122 | 0 | 323098 |
| 1.0 | 90 | 47059 | 504 | 51324 | 116 | 0 | 322866 |

### Encounter speed

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | Encounter speed: sampled minimum (m/s) | 45118 | 42.4535 | 16.3628 | 0.105574 | 13.9367 | 31.829 | 44.2914 | 53.637 | 64.7426 | 194.253 |
| 0.0 | Encounter speed: sampled maximum (m/s) | 45118 | 187.556 | 12.8619 | 3.67418 | 168.154 | 177.756 | 186.45 | 197.768 | 208.662 | 215.474 |
| 0.5 | Encounter speed: sampled minimum (m/s) | 45387 | 60.7821 | 21.8745 | 0.208303 | 26.837 | 48.3243 | 61.8516 | 73.3176 | 86.0694 | 210.327 |
| 0.5 | Encounter speed: sampled maximum (m/s) | 45387 | 191.2 | 12.4733 | 69.9326 | 171.617 | 181.724 | 190.997 | 200.533 | 211.9 | 224.069 |
| 1.0 | Encounter speed: sampled minimum (m/s) | 47059 | 86.0088 | 41.0256 | 0.102125 | 29.7058 | 55.9534 | 80.3422 | 112.803 | 143.598 | 238.481 |
| 1.0 | Encounter speed: sampled maximum (m/s) | 47059 | 207.593 | 20.4986 | 110.541 | 170.402 | 195.626 | 211.23 | 222.324 | 235.633 | 249.815 |

### Signed range rate

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | Signed range rate: sampled minimum (m/s) | 45118 | -112.308 | 30.3018 | -156.873 | -148.704 | -135.133 | -121.698 | -89.9159 | -63.5935 | 125.883 |
| 0.0 | Signed range rate: sampled maximum (m/s) | 45118 | 112.356 | 30.1784 | -129.039 | 63.4031 | 89.9324 | 121.684 | 135.126 | 148.709 | 156.895 |
| 0.5 | Signed range rate: sampled minimum (m/s) | 45387 | -109.58 | 28.564 | -157.893 | -145.473 | -131.789 | -116.337 | -88.8864 | -65.1855 | 126.492 |
| 0.5 | Signed range rate: sampled maximum (m/s) | 45387 | 109.442 | 29.1981 | -137.886 | 59.9273 | 89.2967 | 115.1 | 130.84 | 148.77 | 162.293 |
| 1.0 | Signed range rate: sampled minimum (m/s) | 47059 | -110.538 | 33.2791 | -176.038 | -159.653 | -134.901 | -112.439 | -90.8221 | -50.5914 | 128.321 |
| 1.0 | Signed range rate: sampled maximum (m/s) | 47059 | 110.734 | 33.7462 | -152.833 | 53.0695 | 88.9093 | 112.319 | 135.095 | 161.802 | 178.793 |

### Encounter duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | Encounter duration (s) | 45118 | 4046.49 | 17174.3 | 17.7828 | 1880.89 | 2196.34 | 2504.49 | 4852.37 | 6039.69 | 630712 |
| 0.5 | Encounter duration (s) | 45387 | 3905.12 | 16452.9 | 2.33803 | 1794.91 | 2109.51 | 2417.03 | 4691.16 | 5818.48 | 630712 |
| 1.0 | Encounter duration (s) | 47059 | 3377.69 | 14632.5 | 4.43339 | 1500.95 | 1801.53 | 2096.84 | 4146.32 | 5082.45 | 630712 |

### Laser-on duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | Uninterrupted laser-on duration (s) | 22799 | 996.251 | 1658.99 | 28.0003 | 351.42 | 465.748 | 669.448 | 1030.79 | 2437.73 | 53895.9 |
| 0.0 | Total laser-on within an encounter (s) | 45118 | 502.581 | 1293.82 | 0 | 0 | 0 | 313.466 | 671.335 | 1648.97 | 53895.9 |
| 0.5 | Uninterrupted laser-on duration (s) | 50552 | 699.114 | 1211.72 | 0.37341 | 19.0053 | 274.944 | 512.12 | 879.294 | 1582.8 | 56649.6 |
| 0.5 | Total laser-on within an encounter (s) | 45387 | 777.802 | 4247.88 | 0 | 0 | 331.846 | 523.171 | 857.399 | 1779.26 | 313359 |
| 1.0 | Uninterrupted laser-on duration (s) | 51324 | 676.093 | 1089.32 | 1.83015 | 17.2545 | 260.391 | 500.435 | 858.861 | 1672.29 | 61346.1 |
| 1.0 | Total laser-on within an encounter (s) | 47059 | 736.454 | 4372.7 | 0 | 0 | 279.565 | 497.972 | 818.329 | 1829.81 | 314471 |

## 2) Distribution by Initial target altitude h (km)

### Coverage and case totals

| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 900.0 | 54 | 42845 | 313 | 33146 | 65 | 9 | 312829 |
| 950.0 | 54 | 20589 | 359 | 19524 | 65 | 9 | 316651 |
| 1000.0 | 54 | 4371 | 287 | 10444 | 62 | 9 | 319789 |
| 1050.0 | 54 | 23976 | 348 | 23731 | 63 | 9 | 323124 |
| 1100.0 | 54 | 45783 | 303 | 37830 | 65 | 9 | 325778 |

### Encounter speed

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 900.0 | Encounter speed: sampled minimum (m/s) | 42845 | 69.3087 | 26.8095 | 26.5684 | 37.2655 | 50.3482 | 62.2262 | 80.4662 | 130.223 | 188.86 |
| 900.0 | Encounter speed: sampled maximum (m/s) | 42845 | 192.493 | 15.7646 | 124.796 | 171.485 | 181.006 | 189.257 | 201.142 | 223.686 | 240.774 |
| 950.0 | Encounter speed: sampled minimum (m/s) | 20589 | 43.3331 | 24.8758 | 1.98732 | 11.9519 | 25.8692 | 37.3035 | 58.0876 | 88.9783 | 227.711 |
| 950.0 | Encounter speed: sampled maximum (m/s) | 20589 | 208.431 | 15.0638 | 135.047 | 189.845 | 199.738 | 206.476 | 215.307 | 236.533 | 249.815 |
| 1000.0 | Encounter speed: sampled minimum (m/s) | 4371 | 115.928 | 92.7348 | 0.102125 | 2.5152 | 11.8731 | 160.853 | 202.15 | 230.884 | 238.481 |
| 1000.0 | Encounter speed: sampled maximum (m/s) | 4371 | 210.544 | 18.4495 | 3.67418 | 192.792 | 201.492 | 208.072 | 222.132 | 235.737 | 240.58 |
| 1050.0 | Encounter speed: sampled minimum (m/s) | 23976 | 46.6248 | 24.3619 | 7.21222 | 16.6296 | 28.7097 | 40.1967 | 61.6914 | 95.2551 | 221.682 |
| 1050.0 | Encounter speed: sampled maximum (m/s) | 23976 | 203.202 | 17.5823 | 134.836 | 176.125 | 194.77 | 202.204 | 212.282 | 232.587 | 245.524 |
| 1100.0 | Encounter speed: sampled minimum (m/s) | 45783 | 70.6662 | 25.9363 | 32.0554 | 40.3485 | 52.1933 | 63.5641 | 81.2157 | 130.297 | 210.482 |
| 1100.0 | Encounter speed: sampled maximum (m/s) | 45783 | 187.368 | 15.7648 | 110.541 | 167.063 | 175.921 | 183.756 | 195.728 | 218.785 | 234.462 |

### Signed range rate

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 900.0 | Signed range rate: sampled minimum (m/s) | 42845 | -124.102 | 24.4351 | -176.038 | -154.525 | -139.809 | -127.89 | -113.276 | -77.0776 | 119.802 |
| 900.0 | Signed range rate: sampled maximum (m/s) | 42845 | 124.16 | 24.3468 | -152.833 | 76.3727 | 112.488 | 127.647 | 139.479 | 156.13 | 178.793 |
| 950.0 | Signed range rate: sampled minimum (m/s) | 20589 | -86.9873 | 21.2495 | -130.927 | -120.409 | -101.034 | -86.2666 | -73.4582 | -59.081 | 72.9268 |
| 950.0 | Signed range rate: sampled maximum (m/s) | 20589 | 87.1197 | 21.1253 | -104.896 | 58.0719 | 73.6856 | 84.8479 | 102.073 | 120.187 | 137.351 |
| 1000.0 | Signed range rate: sampled minimum (m/s) | 4371 | -45.287 | 30.1462 | -134.084 | -117.282 | -51.4607 | -42.0489 | -27.1317 | -8.05985 | 2.90846 |
| 1000.0 | Signed range rate: sampled maximum (m/s) | 4371 | 44.9672 | 30.2802 | -42.5444 | 7.98528 | 26.7433 | 41.8398 | 50.5942 | 117.316 | 134.086 |
| 1050.0 | Signed range rate: sampled minimum (m/s) | 23976 | -91.4488 | 18.892 | -136.245 | -118.988 | -104.734 | -91.1672 | -78.8891 | -66.4289 | 74.2421 |
| 1050.0 | Signed range rate: sampled maximum (m/s) | 23976 | 91.3577 | 21.2545 | -94.061 | 61.098 | 76.4625 | 91.2714 | 105.14 | 126.44 | 141.767 |
| 1100.0 | Signed range rate: sampled minimum (m/s) | 45783 | -125.457 | 24.0644 | -174.088 | -154.628 | -140.732 | -129.89 | -116.17 | -75.0633 | 128.321 |
| 1100.0 | Signed range rate: sampled maximum (m/s) | 45783 | 125.532 | 24.2884 | -143.893 | 74.2829 | 115.164 | 129.08 | 140.873 | 156.505 | 178.558 |

### Encounter duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 900.0 | Encounter duration (s) | 42845 | 2216.12 | 358.595 | 11.4819 | 1666.63 | 1962.45 | 2210.01 | 2452.97 | 2804.77 | 3229.02 |
| 950.0 | Encounter duration (s) | 20589 | 5242.24 | 923.948 | 2.33803 | 4213.72 | 4772.86 | 5141.64 | 5558.93 | 6958.24 | 8785.21 |
| 1000.0 | Encounter duration (s) | 4371 | 25768.5 | 87198.2 | 4.43339 | 510.752 | 1345.86 | 2468.33 | 20558.3 | 65058.2 | 630712 |
| 1050.0 | Encounter duration (s) | 23976 | 4509.42 | 674.378 | 152.9 | 3393.05 | 4046.91 | 4633.06 | 5028.5 | 5343.25 | 5457.04 |
| 1100.0 | Encounter duration (s) | 45783 | 2077.8 | 320.563 | 11.995 | 1578.9 | 1858.24 | 2071.38 | 2306.21 | 2610.75 | 2825.65 |

### Laser-on duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 900.0 | Uninterrupted laser-on duration (s) | 33146 | 538.005 | 317.417 | 2.42733 | 32.6096 | 361.518 | 488.759 | 693.863 | 1154.03 | 1655.37 |
| 900.0 | Total laser-on within an encounter (s) | 42845 | 414.345 | 361.73 | 0 | 0 | 0 | 413.96 | 616.112 | 1103.01 | 1649.83 |
| 950.0 | Uninterrupted laser-on duration (s) | 19524 | 945.261 | 807.709 | 1.83015 | 12.6915 | 300.394 | 809.216 | 1362.55 | 2584.45 | 4022.63 |
| 950.0 | Total laser-on within an encounter (s) | 20589 | 894.612 | 835.81 | 0 | 0 | 26.2811 | 766.521 | 1323.48 | 2589.62 | 4018.46 |
| 1000.0 | Uninterrupted laser-on duration (s) | 10444 | 1842.21 | 3841.39 | 0.37341 | 11.7872 | 25.9359 | 1538.16 | 1569.84 | 6513.76 | 61346.1 |
| 1000.0 | Total laser-on within an encounter (s) | 4371 | 4401.29 | 19682.2 | 0 | 0 | 0 | 1018.59 | 4226.82 | 12854.2 | 314471 |
| 1050.0 | Uninterrupted laser-on duration (s) | 23731 | 789.661 | 610.815 | 1.93627 | 16.5543 | 388.702 | 707.562 | 1095.68 | 1974.63 | 3318.31 |
| 1050.0 | Total laser-on within an encounter (s) | 23976 | 781.561 | 642.781 | 0 | 0 | 298.817 | 700.128 | 1097.4 | 2042.07 | 3318.31 |
| 1100.0 | Uninterrupted laser-on duration (s) | 37830 | 488.698 | 295.78 | 2.24995 | 32.5341 | 329.212 | 440.61 | 627.832 | 1065.25 | 1597.04 |
| 1100.0 | Total laser-on within an encounter (s) | 45783 | 403.771 | 326.891 | 0 | 0 | 51.9058 | 392.018 | 572.293 | 1029.05 | 1597.04 |

## 3) Distribution by Helper number

### Coverage and case totals

| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 30 | 70 | 30 | 462 | 19 | 5 | 3194.08 |
| 100 | 30 | 3749 | 46 | 4076 | 24 | 5 | 229413 |
| 200 | 30 | 7610 | 102 | 7310 | 37 | 5 | 358994 |
| 300 | 30 | 11488 | 141 | 10819 | 40 | 5 | 448527 |
| 400 | 30 | 15119 | 170 | 13943 | 40 | 5 | 471244 |
| 500 | 30 | 19165 | 220 | 17464 | 40 | 5 | 471802 |
| 600 | 30 | 22880 | 256 | 20522 | 40 | 5 | 471925 |
| 700 | 30 | 26837 | 290 | 23588 | 40 | 5 | 472035 |
| 800 | 30 | 30646 | 355 | 26491 | 40 | 5 | 472056 |

### Encounter speed

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Encounter speed: sampled minimum (m/s) | 70 | 80.5053 | 57.9455 | 0.801932 | 3.78803 | 36.8213 | 64.2991 | 127.414 | 188.86 | 236.45 |
| 1 | Encounter speed: sampled maximum (m/s) | 70 | 192.836 | 17.7603 | 128.595 | 171.172 | 181.982 | 190.733 | 202.761 | 223.915 | 237.546 |
| 100 | Encounter speed: sampled minimum (m/s) | 3749 | 62.8732 | 31.6693 | 0.832558 | 21.8384 | 42.7809 | 57.391 | 76.4217 | 128.378 | 236.45 |
| 100 | Encounter speed: sampled maximum (m/s) | 3749 | 195.862 | 17.486 | 128.595 | 171.039 | 182.101 | 194.746 | 207.411 | 227.623 | 248.214 |
| 200 | Encounter speed: sampled minimum (m/s) | 7610 | 63.1832 | 33.7203 | 0.585175 | 21.4046 | 41.8657 | 57.1315 | 75.7496 | 130.661 | 238.09 |
| 200 | Encounter speed: sampled maximum (m/s) | 7610 | 195.912 | 17.7937 | 102.876 | 170.44 | 182.077 | 194.824 | 207.847 | 227.24 | 248.593 |
| 300 | Encounter speed: sampled minimum (m/s) | 11488 | 63.5472 | 33.4953 | 0.235534 | 20.9085 | 42.2523 | 57.4935 | 76.4212 | 131.917 | 234.307 |
| 300 | Encounter speed: sampled maximum (m/s) | 11488 | 195.722 | 17.8489 | 25.7192 | 170.153 | 182.026 | 194.778 | 207.611 | 227.378 | 249.045 |
| 400 | Encounter speed: sampled minimum (m/s) | 15119 | 61.9616 | 30.9074 | 0.208303 | 20.7951 | 41.9158 | 57.2033 | 75.3323 | 126.573 | 236.29 |
| 400 | Encounter speed: sampled maximum (m/s) | 15119 | 195.381 | 17.9949 | 12.8596 | 169.879 | 181.694 | 194.384 | 207.367 | 227.405 | 248.753 |
| 500 | Encounter speed: sampled minimum (m/s) | 19165 | 63.8325 | 34.7566 | 0.370951 | 20.708 | 42.1852 | 57.6066 | 76.3154 | 131.961 | 237.897 |
| 500 | Encounter speed: sampled maximum (m/s) | 19165 | 195.695 | 18.1113 | 36.0069 | 169.749 | 181.837 | 194.728 | 207.835 | 228.267 | 249.508 |
| 600 | Encounter speed: sampled minimum (m/s) | 22880 | 62.9147 | 32.5209 | 0.105574 | 20.9167 | 42.0601 | 57.4905 | 75.8484 | 129.674 | 236.999 |
| 600 | Encounter speed: sampled maximum (m/s) | 22880 | 195.457 | 18.0016 | 25.7192 | 169.681 | 181.789 | 194.493 | 207.417 | 227.509 | 249.439 |
| 700 | Encounter speed: sampled minimum (m/s) | 26837 | 64.0377 | 35.4885 | 0.102125 | 20.8074 | 42.2066 | 57.6189 | 76.3546 | 131.938 | 238.481 |
| 700 | Encounter speed: sampled maximum (m/s) | 26837 | 195.699 | 18.2556 | 3.67418 | 169.716 | 181.73 | 194.678 | 207.811 | 228.753 | 249.78 |
| 800 | Encounter speed: sampled minimum (m/s) | 30646 | 63.6693 | 34.2547 | 0.167557 | 20.7892 | 42.1969 | 57.6405 | 76.2565 | 131.584 | 238.129 |
| 800 | Encounter speed: sampled maximum (m/s) | 30646 | 195.576 | 18.1399 | 12.8596 | 169.651 | 181.717 | 194.624 | 207.72 | 227.803 | 249.815 |

### Signed range rate

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Signed range rate: sampled minimum (m/s) | 70 | -67.0719 | 81.943 | -154.037 | -145.252 | -128.818 | -102.577 | -46.7 | 119.047 | 119.802 |
| 1 | Signed range rate: sampled maximum (m/s) | 70 | 102.304 | 35.8552 | 0.155869 | 41.5299 | 79.9371 | 108.495 | 130.573 | 150.178 | 157.119 |
| 100 | Signed range rate: sampled minimum (m/s) | 3749 | -112.336 | 30.2209 | -174.496 | -151.524 | -133.723 | -117.424 | -90.6303 | -66.1647 | 119.802 |
| 100 | Signed range rate: sampled maximum (m/s) | 3749 | 112.216 | 29.2263 | -84.9859 | 64.5124 | 90.4483 | 116.627 | 133.307 | 152.676 | 176.183 |
| 200 | Signed range rate: sampled minimum (m/s) | 7610 | -110.78 | 31.0986 | -173.822 | -151.034 | -133.566 | -116.109 | -89.4901 | -62.764 | 119.802 |
| 200 | Signed range rate: sampled maximum (m/s) | 7610 | 110.79 | 31.1256 | -112.366 | 60.1824 | 89.0324 | 116.085 | 133.533 | 152.286 | 178.71 |
| 300 | Signed range rate: sampled minimum (m/s) | 11488 | -110.636 | 30.6075 | -174.76 | -151.355 | -133.553 | -116.003 | -89.3322 | -60.6214 | 119.802 |
| 300 | Signed range rate: sampled maximum (m/s) | 11488 | 110.668 | 31.0414 | -133.606 | 58.4906 | 89.0492 | 115.598 | 133.567 | 152.827 | 177.655 |
| 400 | Signed range rate: sampled minimum (m/s) | 15119 | -111.609 | 29.722 | -175.923 | -151.366 | -133.871 | -116.574 | -90.73 | -64.3953 | 119.802 |
| 400 | Signed range rate: sampled maximum (m/s) | 15119 | 111.661 | 30.2072 | -144.909 | 62.7178 | 90.5477 | 116.381 | 133.862 | 152.973 | 178.748 |
| 500 | Signed range rate: sampled minimum (m/s) | 19165 | -110.533 | 31.1433 | -176.038 | -151.32 | -133.7 | -116.093 | -89.5399 | -60.4909 | 125.675 |
| 500 | Signed range rate: sampled maximum (m/s) | 19165 | 110.561 | 31.5641 | -152.833 | 59.0062 | 88.9005 | 115.68 | 133.699 | 152.748 | 178.53 |
| 600 | Signed range rate: sampled minimum (m/s) | 22880 | -111.066 | 30.2914 | -176.034 | -151.5 | -133.892 | -116.139 | -90.0041 | -62.2711 | 119.802 |
| 600 | Signed range rate: sampled maximum (m/s) | 22880 | 111.073 | 30.7854 | -133.641 | 60.367 | 89.7012 | 115.858 | 133.842 | 152.936 | 178.745 |
| 700 | Signed range rate: sampled minimum (m/s) | 26837 | -110.434 | 31.3451 | -175.764 | -151.253 | -133.729 | -116.067 | -89.6018 | -60.6872 | 122.07 |
| 700 | Signed range rate: sampled maximum (m/s) | 26837 | 110.438 | 31.8352 | -140.816 | 58.617 | 89.1138 | 115.726 | 133.642 | 152.952 | 178.793 |
| 800 | Signed range rate: sampled minimum (m/s) | 30646 | -110.68 | 30.9037 | -175.87 | -151.365 | -133.757 | -116.034 | -89.7064 | -60.9994 | 128.321 |
| 800 | Signed range rate: sampled maximum (m/s) | 30646 | 110.712 | 31.3039 | -147.289 | 59.2297 | 89.113 | 115.815 | 133.726 | 152.955 | 178.741 |

### Encounter duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Encounter duration (s) | 70 | 33652 | 127942 | 170.911 | 293.061 | 1868.84 | 2209.18 | 4505.72 | 104669 | 630712 |
| 100 | Encounter duration (s) | 3749 | 3835.41 | 18331.6 | 110.84 | 1684.43 | 2029.9 | 2354.95 | 4623.12 | 5337.88 | 630712 |
| 200 | Encounter duration (s) | 7610 | 3756.5 | 15250.9 | 2.33803 | 1653.29 | 2028.02 | 2362.44 | 4554.95 | 5413.43 | 630712 |
| 300 | Encounter duration (s) | 11488 | 3808.03 | 16412.2 | 4.43339 | 1664.45 | 2023.7 | 2368.74 | 4491.9 | 5551.19 | 630712 |
| 400 | Encounter duration (s) | 15119 | 3748.23 | 15764.9 | 17.7828 | 1655.86 | 2021.18 | 2374.22 | 4495.05 | 5614.99 | 630712 |
| 500 | Encounter duration (s) | 19165 | 3734.18 | 15334.1 | 8.51115 | 1622.13 | 2004.18 | 2367.48 | 4480.26 | 5645.36 | 630712 |
| 600 | Encounter duration (s) | 22880 | 3758.22 | 15472.7 | 64.5401 | 1642.13 | 2012.57 | 2376.26 | 4482.43 | 5681.83 | 630712 |
| 700 | Encounter duration (s) | 26837 | 3743.5 | 16015.8 | 16.5328 | 1601.96 | 1998.47 | 2365.83 | 4469.56 | 5664.69 | 630712 |
| 800 | Encounter duration (s) | 30646 | 3752.74 | 15874.9 | 11.4819 | 1618.55 | 2002.34 | 2370.05 | 4469.66 | 5687.93 | 630712 |

### Laser-on duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Uninterrupted laser-on duration (s) | 462 | 1874.02 | 4494.37 | 199.642 | 1020.76 | 1547.57 | 1561.19 | 1572.5 | 1585.46 | 61340.8 |
| 1 | Total laser-on within an encounter (s) | 70 | 12366.6 | 52948.6 | 0 | 0 | 293.061 | 1045.43 | 1730.03 | 55405.7 | 314279 |
| 100 | Uninterrupted laser-on duration (s) | 4076 | 1462.96 | 2491.05 | 7.59931 | 238.749 | 915.682 | 1248.61 | 1571.76 | 2890.54 | 61342.6 |
| 100 | Total laser-on within an encounter (s) | 3749 | 1588.84 | 7689.1 | 0 | 0 | 743.048 | 1122.25 | 1810.29 | 3019.46 | 314471 |
| 200 | Uninterrupted laser-on duration (s) | 7310 | 1394.5 | 2222.08 | 7.01114 | 102.174 | 826.623 | 1132.01 | 1558.09 | 2854.94 | 61342.2 |
| 200 | Total laser-on within an encounter (s) | 7610 | 1338.58 | 5530.21 | 0 | 0 | 502.089 | 1068.23 | 1553.19 | 2830.3 | 313928 |
| 300 | Uninterrupted laser-on duration (s) | 10819 | 1151.63 | 1831.94 | 2.66242 | 34.2263 | 598.979 | 967.389 | 1369.98 | 2481.18 | 61346.1 |
| 300 | Total laser-on within an encounter (s) | 11488 | 1083.71 | 4485.01 | 0 | 0 | 204.564 | 909.991 | 1230.48 | 2435.37 | 312584 |
| 400 | Uninterrupted laser-on duration (s) | 13943 | 906.34 | 1362.84 | 2.48926 | 26.7319 | 504.24 | 745.924 | 1060.13 | 1959.15 | 47098.2 |
| 400 | Total laser-on within an encounter (s) | 15119 | 834.981 | 3781.15 | 0 | 0 | 115.286 | 693.474 | 961.16 | 1890.91 | 313352 |
| 500 | Uninterrupted laser-on duration (s) | 17464 | 724.73 | 1021.51 | 2.12354 | 20.6554 | 389.306 | 596.406 | 847.157 | 1598.71 | 40490.1 |
| 500 | Total laser-on within an encounter (s) | 19165 | 659.569 | 2653.78 | 0 | 0 | 35.6565 | 551.655 | 770.108 | 1529.16 | 233920 |
| 600 | Uninterrupted laser-on duration (s) | 20522 | 616.81 | 856.6 | 1.93627 | 18.3238 | 343.255 | 502.669 | 710.469 | 1551.59 | 38184.2 |
| 600 | Total laser-on within an encounter (s) | 22880 | 552.392 | 2706.07 | 0 | 0 | 26.4082 | 463.236 | 644.364 | 1283.24 | 294315 |
| 700 | Uninterrupted laser-on duration (s) | 23588 | 536.744 | 711.349 | 1.9945 | 15.1429 | 305.312 | 433.834 | 613.596 | 1467.71 | 23521.9 |
| 700 | Total laser-on within an encounter (s) | 26837 | 470.917 | 2756.8 | 0 | 0 | 19.3283 | 396.043 | 554.82 | 1106.32 | 313612 |
| 800 | Uninterrupted laser-on duration (s) | 26491 | 477.906 | 624.719 | 0.37341 | 12.8597 | 276.909 | 381.098 | 540.955 | 1313.96 | 23144.8 |
| 800 | Total laser-on within an encounter (s) | 30646 | 412.27 | 2355.3 | 0 | 0 | 7.65969 | 345.685 | 485.714 | 960.681 | 305558 |

## 4) Distribution by Scheduler

### Coverage and case totals

| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| gve_inc | 135 | 68072 | 787 | 55314 | 78 | 45 | 312716 |
| gve_sma | 135 | 69492 | 823 | 69361 | 242 | 0 | 622547 |

### Encounter speed

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gve_inc | Encounter speed: sampled minimum (m/s) | 68072 | 64.2634 | 34.1327 | 0.102125 | 22.6667 | 42.9594 | 57.7485 | 76.5618 | 132.171 | 238.481 |
| gve_inc | Encounter speed: sampled maximum (m/s) | 68072 | 196.25 | 18.0524 | 3.67418 | 170.584 | 182.325 | 195.216 | 208.285 | 228.518 | 249.192 |
| gve_sma | Encounter speed: sampled minimum (m/s) | 69492 | 62.5552 | 33.4504 | 0.105574 | 18.6417 | 41.4074 | 57.2678 | 75.655 | 129.356 | 237.114 |
| gve_sma | Encounter speed: sampled maximum (m/s) | 69492 | 194.988 | 18.0446 | 110.541 | 169.098 | 181.266 | 194 | 207.036 | 227.364 | 249.815 |

### Signed range rate

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gve_inc | Signed range rate: sampled minimum (m/s) | 68072 | -111.031 | 30.534 | -176.038 | -151.249 | -133.448 | -116.35 | -89.6299 | -63.9773 | 128.321 |
| gve_inc | Signed range rate: sampled maximum (m/s) | 68072 | 111.027 | 31.1492 | -152.833 | 61.2463 | 88.9894 | 116.023 | 133.468 | 152.983 | 178.793 |
| gve_sma | Signed range rate: sampled minimum (m/s) | 69492 | -110.579 | 31.1149 | -174.088 | -151.376 | -134.04 | -115.924 | -89.9793 | -58.5792 | 128.321 |
| gve_sma | Signed range rate: sampled maximum (m/s) | 69492 | 110.656 | 31.1698 | -127.918 | 58.3555 | 89.7764 | 115.732 | 133.967 | 152.721 | 178.558 |

### Encounter duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gve_inc | Encounter duration (s) | 68072 | 3850.86 | 22257.2 | 2.33803 | 1652.74 | 2028.1 | 2347.55 | 4631.98 | 5276.75 | 630712 |
| gve_sma | Encounter duration (s) | 69492 | 3692.89 | 5309.21 | 21.6281 | 1614.83 | 1990.48 | 2398.91 | 4191.04 | 6648.52 | 123078 |

### Laser-on duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gve_inc | Uninterrupted laser-on duration (s) | 55314 | 459.566 | 486.636 | 1.83015 | 11.2786 | 32.97 | 339.191 | 680.939 | 1557.31 | 1598.4 |
| gve_inc | Total laser-on within an encounter (s) | 68072 | 373.379 | 4826.38 | 0 | 0 | 0 | 32.1984 | 478.828 | 1173.63 | 314471 |
| gve_sma | Uninterrupted laser-on duration (s) | 69361 | 970.783 | 1604.31 | 0.37341 | 329.012 | 457.64 | 657.337 | 1022.51 | 2390.62 | 61346.1 |
| gve_sma | Total laser-on within an encounter (s) | 69492 | 967.272 | 1690.12 | 0 | 320.158 | 450.031 | 648.881 | 1010.32 | 2370.93 | 61346.1 |

## Reproducibility and exported data

- [Encounter events](target_helper_encounters.csv): one row per target-helper encounter, with clipping flags, sampled extrema, duration, and accumulated within-encounter laser-on.
- [Laser-on events](target_helper_laser_intervals.csv): one row per uninterrupted firing interval, with within/outside-range durations.
- [Case summary and provenance](case_summary.csv): all 270 cases, including zero-event cases, totals, source directories, and source checksums.
- [Distributions](distributions.csv): all numerical distribution tables, plus unclipped-only sensitivity results. `contributing_cases`, `clipped_events` and `zero_events` accompany each distribution.
- [Evaluation code](../../Code/evaluate_encounters.jl); [tests](../../../Tests/absolute_RN/test_encounter_evaluation.jl).

From the project root:

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/evaluate_encounters.jl
```

The pipeline writes event CSVs and case summaries, reloads them to compute distributions, then reloads the distribution/case CSVs to generate this Markdown. No simulation files are modified. Use `--report-only` to regenerate Markdown from existing CSVs; `--input PATH` and `--output PATH` override the campaign and report directories.
