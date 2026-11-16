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

- All **81 verified cases**, each spanning **100 initial target periods**; target-helper pairs only. Helpers start at 1000 km and 0 deg. Here i and h are the **initial target** inclination and altitude, not their evolving values.
- Exact grid: i = 0.0, 0.5, 1.0 deg; h = 1000.0, 1050.0, 1100.0 km; helpers = 1, 100, 200, 300, 400, 500, 600, 700, 800; schedulers = no_laser.
- Each table pools events over the other three dimensions. **Event-weighted, not equal-case-weighted or time-weighted**: cases producing more events contribute more. These are marginal descriptive comparisons, not isolated causal scheduler effects.
- Speed/range-rate extrema use **10-second samples plus exact boundary states**. They are sampled estimates, not guaranteed continuous-time extrema. Full speed/range-rate time histories for every encounter are not retained; only the longest-encounter tracks are saved. Do not interpret the extrema distributions as instantaneous-speed distributions.
- Main tables include boundary-clipped intervals as observed, without extrapolating full lifetimes. Both start- and end-clipping flags are exported; `unclipped_only` distributions are also in the CSV.
- Each geometric encounter contributes one duration, two speed extrema and two range-rate extrema. Each firing interval contributes one uninterrupted laser duration; it need not match an encounter. Accumulated in-encounter laser-on may combine multiple firing intervals. No-firing cases have no fictional zero-length firing interval; their case total is zero.
- SD is population SD; percentiles use Julia's default linearly interpolated quantile (Hyndman-Fan type 7). NA means no events, not zero. Numbers shown to six significant digits; CSVs retain full precision.
- All input settings, verified status, run IDs, and SHA-256 checksums of the three consumed Feather products were checked. This is an evaluation of saved results, not a new propagation or a rerun of trajectory physics verification.

## Campaign totals

- Target-helper encounters: **34705**, of which **438** are boundary-clipped and **34705** have no within-encounter laser-on time.
- Uninterrupted laser-on intervals: **0**, of which **0** are boundary-clipped.
- Cases without encounters: **0**; cases without firing: **81**.
- Summed firing across independent runs: **0 s**, including **0 s** outside range. These are run-summed totals, not one mission's duration.
- Geometric pair-time can overlap among helpers and must not be interpreted as target access time or duty fraction. The case CSV's laser duty fraction is total actual target firing divided by that case's duration.

| Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Encounter speed: sampled minimum (m/s) | 34705 | 66.0681 | 38.2158 | 3.66526 | 22.5922 | 42.9671 | 57.7175 | 77.3259 | 135.991 | 228.615 |
| Encounter speed: sampled maximum (m/s) | 34705 | 194.893 | 18.4893 | 3.67418 | 169.015 | 180.463 | 194.043 | 207.258 | 226.49 | 245.141 |
| Signed range rate: sampled minimum (m/s) | 34705 | -108.864 | 31.0829 | -170.279 | -149.612 | -131.945 | -114.636 | -87.5594 | -58.7499 | 128.321 |
| Signed range rate: sampled maximum (m/s) | 34705 | 108.777 | 32.1214 | -143.333 | 54.9423 | 86.9692 | 114.454 | 131.775 | 151.065 | 175.059 |
| Encounter duration (s) | 34705 | 4644.13 | 31952.1 | 5.18799 | 1608.11 | 2031.5 | 2351.23 | 4639.04 | 5305.25 | 630712 |
| Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| Total laser-on within an encounter (s) | 34705 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## 1) Distribution by Target inclination i (deg)

### Coverage and case totals

| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 27 | 11165 | 150 | 0 | 0 | 27 | 0 |
| 0.5 | 27 | 11366 | 150 | 0 | 0 | 27 | 0 |
| 1.0 | 27 | 12174 | 138 | 0 | 0 | 27 | 0 |

### Encounter speed

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | Encounter speed: sampled minimum (m/s) | 11165 | 42.432 | 15.8669 | 3.66526 | 16.0622 | 32.2119 | 44.2447 | 52.5039 | 63.9603 | 194.253 |
| 0.0 | Encounter speed: sampled maximum (m/s) | 11165 | 186.018 | 13.4774 | 3.67418 | 167.064 | 177.003 | 184.453 | 195.673 | 207.561 | 211.74 |
| 0.5 | Encounter speed: sampled minimum (m/s) | 11366 | 62.3685 | 25.4905 | 3.68064 | 28.204 | 48.9582 | 61.9106 | 73.2956 | 86.6965 | 199.696 |
| 0.5 | Encounter speed: sampled maximum (m/s) | 11366 | 189.958 | 12.8688 | 64.4094 | 170.85 | 180.575 | 189.506 | 199.477 | 211.037 | 218.993 |
| 1.0 | Encounter speed: sampled minimum (m/s) | 12174 | 91.1991 | 46.574 | 3.71385 | 31.3288 | 57.1582 | 82.9082 | 117.094 | 206.064 | 228.615 |
| 1.0 | Encounter speed: sampled maximum (m/s) | 12174 | 207.639 | 19.8976 | 121.402 | 170.244 | 197.002 | 211.074 | 222.129 | 233.025 | 245.141 |

### Signed range rate

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | Signed range rate: sampled minimum (m/s) | 11165 | -111.983 | 28.5981 | -151.221 | -146.942 | -133.114 | -121.457 | -89.1836 | -65.4163 | 125.883 |
| 0.0 | Signed range rate: sampled maximum (m/s) | 11165 | 111.975 | 28.6732 | -128.225 | 65.345 | 89.4023 | 121.421 | 133.084 | 146.886 | 151.216 |
| 0.5 | Signed range rate: sampled minimum (m/s) | 11366 | -107.996 | 28.9382 | -152.408 | -142.889 | -130.706 | -115.672 | -87.2779 | -69.3945 | 126.492 |
| 0.5 | Signed range rate: sampled maximum (m/s) | 11366 | 107.642 | 30.7126 | -133.894 | 56.5398 | 89.8622 | 113.701 | 128.93 | 148.166 | 157.197 |
| 1.0 | Signed range rate: sampled minimum (m/s) | 12174 | -106.815 | 34.7637 | -170.279 | -156.917 | -133.04 | -109.682 | -86.777 | -40.5706 | 128.321 |
| 1.0 | Signed range rate: sampled maximum (m/s) | 12174 | 106.903 | 35.9527 | -143.333 | 40.6149 | 82.9855 | 109.587 | 133.134 | 160.601 | 175.059 |

### Encounter duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | Encounter duration (s) | 11165 | 5058.69 | 33571 | 54.565 | 1941.42 | 2253.2 | 2497.68 | 5092.17 | 5396.81 | 630712 |
| 0.5 | Encounter duration (s) | 11366 | 4840.68 | 32758.3 | 5.18799 | 1827.53 | 2148.89 | 2397.01 | 4908.04 | 5281.61 | 630712 |
| 1.0 | Encounter duration (s) | 12174 | 4080.41 | 29570.9 | 24.0544 | 1389.29 | 1799.28 | 2076.41 | 4336.93 | 4866.58 | 630712 |

### Laser-on duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 0.0 | Total laser-on within an encounter (s) | 11165 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 0.5 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 0.5 | Total laser-on within an encounter (s) | 11366 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 1.0 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 1.0 | Total laser-on within an encounter (s) | 12174 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## 2) Distribution by Initial target altitude h (km)

### Coverage and case totals

| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1000.0 | 27 | 1298 | 98 | 0 | 0 | 27 | 0 |
| 1050.0 | 27 | 11167 | 178 | 0 | 0 | 27 | 0 |
| 1100.0 | 27 | 22240 | 162 | 0 | 0 | 27 | 0 |

### Encounter speed

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1000.0 | Encounter speed: sampled minimum (m/s) | 1298 | 195.601 | 34.6847 | 3.66526 | 127.426 | 198.28 | 205.4 | 210.679 | 226.683 | 228.615 |
| 1000.0 | Encounter speed: sampled maximum (m/s) | 1298 | 212.285 | 28.398 | 3.67418 | 159.277 | 206.022 | 222.052 | 224.471 | 232.122 | 233.335 |
| 1050.0 | Encounter speed: sampled minimum (m/s) | 11167 | 44.5333 | 24.5276 | 7.21222 | 14.6926 | 26.1828 | 38.4763 | 59.6149 | 91.3159 | 221.682 |
| 1050.0 | Encounter speed: sampled maximum (m/s) | 11167 | 206.029 | 15.4154 | 134.215 | 188.54 | 197.341 | 204.004 | 213.476 | 233.345 | 245.141 |
| 1100.0 | Encounter speed: sampled minimum (m/s) | 22240 | 69.321 | 26.1536 | 32.0554 | 38.7332 | 50.3388 | 62.2871 | 80.0647 | 129.254 | 210.482 |
| 1100.0 | Encounter speed: sampled maximum (m/s) | 22240 | 188.286 | 15.6046 | 121.402 | 167.928 | 177.08 | 184.771 | 196.287 | 219.302 | 233.891 |

### Signed range rate

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1000.0 | Signed range rate: sampled minimum (m/s) | 1298 | -33.9648 | 14.5512 | -124.952 | -44.9732 | -41.0295 | -40.0826 | -29.5699 | -10.0999 | 2.90846 |
| 1000.0 | Signed range rate: sampled maximum (m/s) | 1298 | 34.053 | 14.5318 | -33.9966 | 10.2197 | 29.6769 | 40.104 | 40.9121 | 44.8436 | 124.954 |
| 1050.0 | Signed range rate: sampled minimum (m/s) | 11167 | -88.2596 | 18.9261 | -122.745 | -116.195 | -101.703 | -87.6609 | -75.6029 | -64.4862 | 74.2421 |
| 1050.0 | Signed range rate: sampled maximum (m/s) | 11167 | 88.0555 | 21.9313 | -93.9402 | 58.0786 | 73.0408 | 87.1501 | 102.609 | 125.722 | 135.472 |
| 1100.0 | Signed range rate: sampled minimum (m/s) | 22240 | -123.582 | 23.6272 | -170.279 | -151.933 | -138.587 | -127.988 | -114.063 | -75.2049 | 128.321 |
| 1100.0 | Signed range rate: sampled maximum (m/s) | 22240 | 123.542 | 24.4286 | -143.333 | 74.0688 | 113.006 | 127.097 | 138.606 | 154.697 | 175.059 |

### Encounter duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1000.0 | Encounter duration (s) | 1298 | 45103 | 159844 | 5.18799 | 794.47 | 1290.81 | 1429.16 | 1711.61 | 630712 | 630712 |
| 1050.0 | Encounter duration (s) | 11167 | 4875.19 | 559.248 | 111.306 | 4160.34 | 4660.03 | 4998.26 | 5214.08 | 5393.56 | 5457.04 |
| 1100.0 | Encounter duration (s) | 22240 | 2166.79 | 316.429 | 8.46914 | 1650.57 | 1941.07 | 2180.25 | 2378.91 | 2682.36 | 2825.65 |

### Laser-on duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1000.0 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 1000.0 | Total laser-on within an encounter (s) | 1298 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 1050.0 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 1050.0 | Total laser-on within an encounter (s) | 11167 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 1100.0 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 1100.0 | Total laser-on within an encounter (s) | 22240 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## 3) Distribution by Helper number

### Coverage and case totals

| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 9 | 18 | 9 | 0 | 0 | 9 | 0 |
| 100 | 9 | 931 | 15 | 0 | 0 | 9 | 0 |
| 200 | 9 | 1864 | 28 | 0 | 0 | 9 | 0 |
| 300 | 9 | 2991 | 36 | 0 | 0 | 9 | 0 |
| 400 | 9 | 3719 | 45 | 0 | 0 | 9 | 0 |
| 500 | 9 | 4650 | 59 | 0 | 0 | 9 | 0 |
| 600 | 9 | 5780 | 69 | 0 | 0 | 9 | 0 |
| 700 | 9 | 6909 | 80 | 0 | 0 | 9 | 0 |
| 800 | 9 | 7843 | 97 | 0 | 0 | 9 | 0 |

### Encounter speed

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Encounter speed: sampled minimum (m/s) | 18 | 69.077 | 40.4387 | 24.7848 | 25.0557 | 27.0773 | 56.3956 | 110.679 | 129.138 | 133.984 |
| 1 | Encounter speed: sampled maximum (m/s) | 18 | 179.552 | 17.9824 | 128.595 | 141.492 | 175.048 | 180.616 | 190.155 | 196.976 | 206.79 |
| 100 | Encounter speed: sampled minimum (m/s) | 931 | 61.6641 | 29.0742 | 10.2763 | 21.41 | 41.9544 | 57.1546 | 75.6386 | 125.644 | 196.804 |
| 100 | Encounter speed: sampled maximum (m/s) | 931 | 193.628 | 18.0045 | 128.595 | 168.432 | 179.327 | 193.45 | 205.113 | 225.791 | 242.806 |
| 200 | Encounter speed: sampled minimum (m/s) | 1864 | 61.2937 | 28.7168 | 7.93561 | 22.4718 | 42.0704 | 56.3987 | 74.7415 | 123.959 | 196.804 |
| 200 | Encounter speed: sampled maximum (m/s) | 1864 | 193.977 | 17.8873 | 102.876 | 168.776 | 180.149 | 193.146 | 206.172 | 225.621 | 245.141 |
| 300 | Encounter speed: sampled minimum (m/s) | 2991 | 70.8571 | 45.4883 | 7.21222 | 23.0102 | 43.7637 | 58.6694 | 80.3943 | 204.642 | 221.682 |
| 300 | Encounter speed: sampled maximum (m/s) | 2991 | 195.971 | 18.9401 | 25.7192 | 169.273 | 180.942 | 194.729 | 209.752 | 225.177 | 244.505 |
| 400 | Encounter speed: sampled minimum (m/s) | 3719 | 61.0864 | 28.2531 | 7.42751 | 22.3955 | 42.1961 | 56.6257 | 74.5224 | 122.417 | 196.804 |
| 400 | Encounter speed: sampled maximum (m/s) | 3719 | 194.015 | 18.1327 | 12.8596 | 168.846 | 180.116 | 193.289 | 206.24 | 225.598 | 245.141 |
| 500 | Encounter speed: sampled minimum (m/s) | 4650 | 61.116 | 28.2403 | 7.30537 | 22.4167 | 42.19 | 56.5055 | 74.4317 | 123.137 | 201.887 |
| 500 | Encounter speed: sampled maximum (m/s) | 4650 | 194.073 | 18.1482 | 36.0069 | 168.96 | 180.127 | 193.271 | 206.224 | 225.567 | 244.918 |
| 600 | Encounter speed: sampled minimum (m/s) | 5780 | 66.1432 | 38.433 | 7.21222 | 22.7606 | 43.0084 | 57.6408 | 77.2617 | 135.714 | 221.682 |
| 600 | Encounter speed: sampled maximum (m/s) | 5780 | 195.061 | 18.6035 | 25.7192 | 169.084 | 180.463 | 194.048 | 207.955 | 225.612 | 245.141 |
| 700 | Encounter speed: sampled minimum (m/s) | 6909 | 69.8495 | 44.8284 | 3.66526 | 23.0449 | 43.3684 | 58.5491 | 79.1816 | 198.484 | 228.615 |
| 700 | Encounter speed: sampled maximum (m/s) | 6909 | 195.45 | 18.8066 | 3.67418 | 169.166 | 180.767 | 194.581 | 207.084 | 230.668 | 244.684 |
| 800 | Encounter speed: sampled minimum (m/s) | 7843 | 67.8041 | 40.0583 | 7.27398 | 22.9351 | 43.3289 | 58.302 | 78.7021 | 168.377 | 214.383 |
| 800 | Encounter speed: sampled maximum (m/s) | 7843 | 195.172 | 18.4475 | 12.8596 | 169.115 | 180.696 | 194.432 | 207.578 | 225.885 | 245.141 |

### Signed range rate

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Signed range rate: sampled minimum (m/s) | 18 | -94.4447 | 41.9338 | -154.565 | -140.308 | -128.65 | -111.189 | -69.1291 | -13.0811 | -0.155869 |
| 1 | Signed range rate: sampled maximum (m/s) | 18 | 94.3647 | 39.73 | 0.155869 | 13.0588 | 74.6333 | 108.604 | 126.801 | 134.01 | 148.489 |
| 100 | Signed range rate: sampled minimum (m/s) | 931 | -112.621 | 26.7957 | -169.42 | -150.469 | -133.299 | -117.385 | -90.6751 | -68.1533 | -0.155869 |
| 100 | Signed range rate: sampled maximum (m/s) | 931 | 111.637 | 29.9353 | -114.089 | 65.809 | 91.6954 | 117.293 | 132.022 | 150.558 | 173.597 |
| 200 | Signed range rate: sampled minimum (m/s) | 1864 | -111.248 | 28.9584 | -169.42 | -150.13 | -132.402 | -115.793 | -90.0885 | -67.7444 | 109.267 |
| 200 | Signed range rate: sampled maximum (m/s) | 1864 | 111.261 | 29.7681 | -114.089 | 63.6565 | 90.9584 | 116.169 | 132.208 | 151.124 | 173.597 |
| 300 | Signed range rate: sampled minimum (m/s) | 2991 | -106.743 | 32.3654 | -169.516 | -149.584 | -131.542 | -113.428 | -84.2017 | -40.8657 | 74.2421 |
| 300 | Signed range rate: sampled maximum (m/s) | 2991 | 106.696 | 33.254 | -114.089 | 40.8519 | 83.142 | 113.154 | 131.309 | 150.856 | 174.608 |
| 400 | Signed range rate: sampled minimum (m/s) | 3719 | -111.582 | 27.8539 | -170.095 | -149.793 | -132.507 | -115.795 | -90.4447 | -67.9504 | 109.267 |
| 400 | Signed range rate: sampled maximum (m/s) | 3719 | 111.506 | 29.0675 | -114.089 | 64.1431 | 90.759 | 116.092 | 132.436 | 151.177 | 175.05 |
| 500 | Signed range rate: sampled minimum (m/s) | 4650 | -111.522 | 28.2136 | -170.266 | -149.59 | -132.584 | -115.793 | -90.3518 | -67.9248 | 125.675 |
| 500 | Signed range rate: sampled maximum (m/s) | 4650 | 111.472 | 29.2909 | -133.894 | 63.7626 | 90.9201 | 115.822 | 132.468 | 151.299 | 175.022 |
| 600 | Signed range rate: sampled minimum (m/s) | 5780 | -109.086 | 30.4139 | -170.115 | -149.625 | -131.893 | -114.683 | -87.597 | -58.7233 | 109.267 |
| 600 | Signed range rate: sampled maximum (m/s) | 5780 | 108.995 | 31.5185 | -130.2 | 54.7333 | 86.9331 | 114.488 | 131.806 | 151.059 | 175.059 |
| 700 | Signed range rate: sampled minimum (m/s) | 6909 | -106.295 | 34.5908 | -170.279 | -149.267 | -131.564 | -113.541 | -85.3872 | -31.133 | 122.07 |
| 700 | Signed range rate: sampled maximum (m/s) | 6909 | 106.192 | 35.5823 | -143.333 | 31.0737 | 84.2166 | 113.483 | 131.438 | 150.956 | 174.845 |
| 800 | Signed range rate: sampled minimum (m/s) | 7843 | -107.93 | 31.4949 | -170.095 | -149.369 | -131.653 | -113.944 | -85.9088 | -44.8598 | 128.321 |
| 800 | Signed range rate: sampled maximum (m/s) | 7843 | 107.897 | 32.3644 | -140.471 | 44.7837 | 84.9692 | 113.765 | 131.41 | 150.934 | 175.05 |

### Encounter duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Encounter duration (s) | 18 | 107663 | 233917 | 1649.83 | 1834.75 | 2076.59 | 3158.85 | 4912.06 | 630712 | 630712 |
| 100 | Encounter duration (s) | 931 | 5096.23 | 35596.3 | 214.492 | 1732.87 | 2047.45 | 2379.44 | 4678.45 | 5315.38 | 630712 |
| 200 | Encounter duration (s) | 1864 | 5087.79 | 35577.6 | 111.306 | 1696.48 | 2067.19 | 2379.12 | 4677.09 | 5310.72 | 630712 |
| 300 | Encounter duration (s) | 2991 | 4639.61 | 32450.7 | 152.901 | 1441.11 | 1979.73 | 2329.81 | 4614.26 | 5304.87 | 630712 |
| 400 | Encounter duration (s) | 3719 | 4591.03 | 30867.6 | 111.306 | 1715.99 | 2063.73 | 2378.16 | 4676.49 | 5309.9 | 630712 |
| 500 | Encounter duration (s) | 4650 | 4691.69 | 31871.2 | 8.46914 | 1715.49 | 2060.91 | 2377.28 | 4672.77 | 5308.79 | 630712 |
| 600 | Encounter duration (s) | 5780 | 4534.68 | 30884.5 | 104.151 | 1560.46 | 2023.3 | 2350.47 | 4639.45 | 5305.11 | 630712 |
| 700 | Encounter duration (s) | 6909 | 4419.99 | 30205.2 | 5.18799 | 1577.12 | 1995.52 | 2335.76 | 4620.56 | 5299.25 | 630712 |
| 800 | Encounter duration (s) | 7843 | 4525.39 | 30887.1 | 64.2157 | 1602.34 | 2040.88 | 2343.6 | 4625.12 | 5301.91 | 630712 |

### Laser-on duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 1 | Total laser-on within an encounter (s) | 18 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 100 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 100 | Total laser-on within an encounter (s) | 931 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 200 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 200 | Total laser-on within an encounter (s) | 1864 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 300 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 300 | Total laser-on within an encounter (s) | 2991 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 400 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 400 | Total laser-on within an encounter (s) | 3719 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 500 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 500 | Total laser-on within an encounter (s) | 4650 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 600 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 600 | Total laser-on within an encounter (s) | 5780 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 700 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 700 | Total laser-on within an encounter (s) | 6909 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 800 | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| 800 | Total laser-on within an encounter (s) | 7843 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## 4) Distribution by Scheduler

### Coverage and case totals

| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| no_laser | 81 | 34705 | 438 | 0 | 0 | 81 | 0 |

### Encounter speed

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| no_laser | Encounter speed: sampled minimum (m/s) | 34705 | 66.0681 | 38.2158 | 3.66526 | 22.5922 | 42.9671 | 57.7175 | 77.3259 | 135.991 | 228.615 |
| no_laser | Encounter speed: sampled maximum (m/s) | 34705 | 194.893 | 18.4893 | 3.67418 | 169.015 | 180.463 | 194.043 | 207.258 | 226.49 | 245.141 |

### Signed range rate

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| no_laser | Signed range rate: sampled minimum (m/s) | 34705 | -108.864 | 31.0829 | -170.279 | -149.612 | -131.945 | -114.636 | -87.5594 | -58.7499 | 128.321 |
| no_laser | Signed range rate: sampled maximum (m/s) | 34705 | 108.777 | 32.1214 | -143.333 | 54.9423 | 86.9692 | 114.454 | 131.775 | 151.065 | 175.059 |

### Encounter duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| no_laser | Encounter duration (s) | 34705 | 4644.13 | 31952.1 | 5.18799 | 1608.11 | 2031.5 | 2351.23 | 4639.04 | 5305.25 | 630712 |

### Laser-on duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| no_laser | Uninterrupted laser-on duration (s) | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| no_laser | Total laser-on within an encounter (s) | 34705 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## Reproducibility and exported data

- [Encounter events](target_helper_encounters.csv): one row per target-helper encounter, with clipping flags, sampled extrema, duration, and accumulated within-encounter laser-on.
- [Laser-on events](target_helper_laser_intervals.csv): one row per uninterrupted firing interval, with within/outside-range durations.
- [Case summary and provenance](case_summary.csv): all 81 cases, including zero-event cases, totals, source directories, and source checksums.
- [Distributions](distributions.csv): all numerical distribution tables, plus unclipped-only sensitivity results. `contributing_cases`, `clipped_events` and `zero_events` accompany each distribution.
- [Evaluation code](../../Code/evaluate_encounters.jl); [tests](../../../Tests/absolute_RN/test_encounter_evaluation.jl).

From the project root:

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/evaluate_encounters.jl
```

The pipeline writes event CSVs and case summaries, reloads them to compute distributions, then reloads the distribution/case CSVs to generate this Markdown. No simulation files are modified. Use `--report-only` to regenerate Markdown from existing CSVs; `--input PATH` and `--output PATH` override the campaign and report directories.
