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

- All **162 verified cases**, each spanning **100 initial target periods**; target-helper pairs only. Helpers start at 1000 km and 0 deg. Here i and h are the **initial target** inclination and altitude, not their evolving values.
- Exact grid: i = 0.0, 0.5, 1.0 deg; h = 1000.0, 1050.0, 1100.0 km; helpers = 1, 100, 200, 300, 400, 500, 600, 700, 800; schedulers = naive_next_entering, positive_along_track.
- Each table pools events over the other three dimensions. **Event-weighted, not equal-case-weighted or time-weighted**: cases producing more events contribute more. These are marginal descriptive comparisons, not isolated causal scheduler effects.
- Speed/range-rate extrema use **10-second samples plus exact boundary states**. They are sampled estimates, not guaranteed continuous-time extrema. Full speed/range-rate time histories for every encounter are not retained; only the longest-encounter tracks are saved. Do not interpret the extrema distributions as instantaneous-speed distributions.
- Main tables include boundary-clipped intervals as observed, without extrapolating full lifetimes. Both start- and end-clipping flags are exported; `unclipped_only` distributions are also in the CSV.
- Each geometric encounter contributes one duration, two speed extrema and two range-rate extrema. Each firing interval contributes one uninterrupted laser duration; it need not match an encounter. Accumulated in-encounter laser-on may combine multiple firing intervals. No-firing cases have no fictional zero-length firing interval; their case total is zero.
- SD is population SD; percentiles use Julia's default linearly interpolated quantile (Hyndman-Fan type 7). NA means no events, not zero. Numbers shown to six significant digits; CSVs retain full precision.
- All input settings, verified status, run IDs, and SHA-256 checksums of the three consumed Feather products were checked. This is an evaluation of saved results, not a new propagation or a rerun of trajectory physics verification.

## Campaign totals

- Target-helper encounters: **73419**, of which **942** are boundary-clipped and **39037** have no within-encounter laser-on time.
- Uninterrupted laser-on intervals: **34382**, of which **269** are boundary-clipped.
- Cases without encounters: **0**; cases without firing: **0**.
- Summed firing across independent runs: **8.33805e+07 s**, including **33095 s** outside range. These are run-summed totals, not one mission's duration.
- Geometric pair-time can overlap among helpers and must not be interpreted as target access time or duty fraction. The case CSV's laser duty fraction is total actual target firing divided by that case's duration.

| Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Encounter speed: sampled minimum (m/s) | 73419 | 66.0968 | 38.8308 | 0.0262687 | 20.3318 | 42.604 | 57.9969 | 77.7472 | 137.361 | 237.867 |
| Encounter speed: sampled maximum (m/s) | 73419 | 194.356 | 18.3462 | 4.3819 | 168.435 | 180.067 | 193.502 | 206.715 | 226.769 | 244.95 |
| Signed range rate: sampled minimum (m/s) | 73419 | -108.607 | 32.5385 | -172.882 | -150.269 | -133.062 | -114.96 | -87.4949 | -46.2322 | 128.321 |
| Signed range rate: sampled maximum (m/s) | 73419 | 108.559 | 33.3002 | -137.073 | 46.5718 | 86.6893 | 114.997 | 132.853 | 152.133 | 177.348 |
| Encounter duration (s) | 73419 | 4309.57 | 22018.1 | 0.749804 | 1591.23 | 1982.52 | 2324.74 | 4335.27 | 5308.67 | 630712 |
| Uninterrupted laser-on duration (s) | 34382 | 2425.12 | 17850.4 | 7.18036 | 524.518 | 774.717 | 1195.57 | 2128.32 | 5204.4 | 630712 |
| Total laser-on within an encounter (s) | 73419 | 1135.23 | 12275.2 | 0 | 0 | 0 | 0 | 1098.5 | 3097.65 | 630712 |

## 1) Distribution by Target inclination i (deg)

### Coverage and case totals

| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 54 | 23761 | 325 | 11029 | 92 | 0 | 620340 |
| 0.5 | 54 | 23953 | 319 | 11172 | 90 | 0 | 609867 |
| 1.0 | 54 | 25705 | 298 | 12181 | 87 | 0 | 587503 |

### Encounter speed

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | Encounter speed: sampled minimum (m/s) | 23761 | 42.634 | 17.1636 | 0.0262687 | 13.8637 | 31.5476 | 44.5556 | 53.6017 | 64.7826 | 198.016 |
| 0.0 | Encounter speed: sampled maximum (m/s) | 23761 | 185.611 | 12.8197 | 4.3819 | 166.618 | 176.094 | 184.205 | 195.429 | 206.499 | 211.795 |
| 0.5 | Encounter speed: sampled minimum (m/s) | 23953 | 61.9649 | 26.0201 | 0.25182 | 25.9406 | 48.7422 | 62.1061 | 73.5596 | 87.1305 | 210.111 |
| 0.5 | Encounter speed: sampled maximum (m/s) | 23953 | 189.618 | 12.5318 | 64.3996 | 170.252 | 180.29 | 189.073 | 199.262 | 210.211 | 218.996 |
| 1.0 | Encounter speed: sampled minimum (m/s) | 25705 | 91.6355 | 47.1855 | 0.182369 | 29.9353 | 57.2704 | 83.7224 | 118.908 | 194.399 | 237.867 |
| 1.0 | Encounter speed: sampled maximum (m/s) | 25705 | 206.855 | 20.378 | 113.08 | 169.214 | 196.224 | 210.767 | 221.205 | 233.682 | 244.95 |

### Signed range rate

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | Signed range rate: sampled minimum (m/s) | 23761 | -111.391 | 31.3843 | -155.319 | -148.07 | -134.749 | -121.836 | -88.7882 | -61.4786 | 125.883 |
| 0.0 | Signed range rate: sampled maximum (m/s) | 23761 | 111.364 | 31.5168 | -133.269 | 61.1017 | 88.7485 | 121.824 | 134.712 | 148.1 | 155.377 |
| 0.5 | Signed range rate: sampled minimum (m/s) | 23953 | -108.349 | 29.6346 | -156.049 | -143.898 | -131.675 | -116.078 | -87.6913 | -65.2559 | 126.492 |
| 0.5 | Signed range rate: sampled maximum (m/s) | 23953 | 108.032 | 31.0553 | -131.302 | 55.7114 | 88.9593 | 114.577 | 129.992 | 148.811 | 160.999 |
| 1.0 | Signed range rate: sampled minimum (m/s) | 25705 | -106.274 | 35.8123 | -172.882 | -157.761 | -132.123 | -110.37 | -85.3034 | -40.8474 | 128.321 |
| 1.0 | Signed range rate: sampled maximum (m/s) | 25705 | 106.458 | 36.5848 | -137.073 | 40.9469 | 82.4697 | 109.952 | 133.104 | 160.849 | 177.348 |

### Encounter duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | Encounter duration (s) | 23761 | 4695.03 | 23316.5 | 40.6091 | 1888.83 | 2172.74 | 2460.6 | 4712.76 | 5440.96 | 630712 |
| 0.5 | Encounter duration (s) | 23953 | 4506.66 | 22929.6 | 0.749804 | 1778.01 | 2080.25 | 2364.82 | 4525.11 | 5287.22 | 630712 |
| 1.0 | Encounter duration (s) | 25705 | 3769.6 | 19792.4 | 44.0145 | 1404.83 | 1769.8 | 2048.7 | 4087.09 | 4829.1 | 630712 |

### Laser-on duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | Uninterrupted laser-on duration (s) | 11029 | 2556.69 | 18167.2 | 73.8645 | 596.738 | 829.257 | 1280.94 | 2287.6 | 5363.64 | 630712 |
| 0.0 | Total laser-on within an encounter (s) | 23761 | 1186.28 | 12442.7 | 0 | 0 | 0 | 0 | 1161.83 | 4763.87 | 630712 |
| 0.5 | Uninterrupted laser-on duration (s) | 11172 | 2502.89 | 18087.1 | 7.18036 | 582.396 | 810.635 | 1225.58 | 2193.39 | 5174.97 | 630712 |
| 0.5 | Total laser-on within an encounter (s) | 23953 | 1166.93 | 12415.3 | 0 | 0 | 0 | 0 | 1120.49 | 3116.47 | 630712 |
| 1.0 | Uninterrupted laser-on duration (s) | 12181 | 2234.67 | 17333.3 | 76.5664 | 439.193 | 710.722 | 1086.66 | 1921.47 | 4752.55 | 630712 |
| 1.0 | Total laser-on within an encounter (s) | 25705 | 1058.5 | 11984 | 0 | 0 | 0 | 0 | 1016.04 | 2764.65 | 630712 |

## 2) Distribution by Initial target altitude h (km)

### Coverage and case totals

| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1000.0 | 54 | 4606 | 284 | 1117 | 75 | 0 | 630712 |
| 1050.0 | 54 | 23481 | 350 | 10428 | 98 | 0 | 537126 |
| 1100.0 | 54 | 45332 | 308 | 22837 | 96 | 0 | 572188 |

### Encounter speed

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1000.0 | Encounter speed: sampled minimum (m/s) | 4606 | 129.243 | 90.7512 | 0.0262687 | 2.6093 | 12.3003 | 174.874 | 203.717 | 230.295 | 237.867 |
| 1000.0 | Encounter speed: sampled maximum (m/s) | 4606 | 211.091 | 17.7756 | 4.3819 | 194.832 | 202.712 | 209.022 | 221.721 | 234.94 | 239.295 |
| 1050.0 | Encounter speed: sampled minimum (m/s) | 23481 | 45.8857 | 24.0145 | 7.24447 | 16.0453 | 28.4332 | 39.5663 | 60.8053 | 93.7511 | 221.682 |
| 1050.0 | Encounter speed: sampled maximum (m/s) | 23481 | 203.963 | 16.9233 | 133.759 | 180.501 | 195.71 | 202.618 | 212.64 | 232.191 | 244.95 |
| 1100.0 | Encounter speed: sampled minimum (m/s) | 45332 | 70.1497 | 25.9496 | 32.0355 | 39.8031 | 51.6226 | 63.0355 | 80.6635 | 129.84 | 210.482 |
| 1100.0 | Encounter speed: sampled maximum (m/s) | 45332 | 187.679 | 15.6151 | 113.08 | 167.54 | 176.323 | 184.158 | 195.796 | 218.761 | 234.103 |

### Signed range rate

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1000.0 | Signed range rate: sampled minimum (m/s) | 4606 | -42.0054 | 28.1241 | -128.44 | -117.088 | -46.1964 | -40.8854 | -25.6589 | -6.58285 | 2.90846 |
| 1000.0 | Signed range rate: sampled maximum (m/s) | 4606 | 41.8105 | 28.2003 | -37.3216 | 6.61728 | 26.0152 | 40.9055 | 46.4593 | 116.962 | 128.445 |
| 1050.0 | Signed range rate: sampled minimum (m/s) | 23481 | -90.306 | 18.6282 | -128.085 | -118.065 | -103.185 | -89.7257 | -78.2686 | -65.8156 | 74.2421 |
| 1050.0 | Signed range rate: sampled maximum (m/s) | 23481 | 90.1032 | 21.1722 | -96.1694 | 60.3146 | 75.667 | 89.3769 | 104.251 | 125.646 | 139.869 |
| 1100.0 | Signed range rate: sampled minimum (m/s) | 45332 | -124.854 | 23.8009 | -172.882 | -153.661 | -139.846 | -129.276 | -115.521 | -75.6298 | 128.321 |
| 1100.0 | Signed range rate: sampled maximum (m/s) | 45332 | 124.901 | 24.2235 | -137.073 | 74.2269 | 114.54 | 128.448 | 139.886 | 155.703 | 177.348 |

### Encounter duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1000.0 | Encounter duration (s) | 4606 | 24410.6 | 85278.9 | 36.4059 | 529.842 | 1305.33 | 2052.11 | 21711.5 | 64833 | 630712 |
| 1050.0 | Encounter duration (s) | 23481 | 4614.49 | 619.924 | 62.1945 | 3655.64 | 4239.47 | 4701.15 | 5059.31 | 5372.66 | 5574.12 |
| 1100.0 | Encounter duration (s) | 45332 | 2109.24 | 312.566 | 0.749804 | 1614.8 | 1891.03 | 2110.53 | 2325.5 | 2620.13 | 2826.35 |

### Laser-on duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1000.0 | Uninterrupted laser-on duration (s) | 1117 | 27350.3 | 95524.1 | 99.8337 | 1222.02 | 6913.32 | 9713.87 | 14530.4 | 45337.2 | 630712 |
| 1000.0 | Total laser-on within an encounter (s) | 4606 | 6632.62 | 48479.7 | 0 | 0 | 0 | 0 | 0 | 16739.7 | 630712 |
| 1050.0 | Uninterrupted laser-on duration (s) | 10428 | 2552.69 | 1451.63 | 7.18036 | 1085.42 | 1508.35 | 1925.77 | 3200.32 | 5267.96 | 5576.77 |
| 1050.0 | Total laser-on within an encounter (s) | 23481 | 1133.23 | 1594.15 | 0 | 0 | 0 | 0 | 1802.04 | 5062.16 | 5571.05 |
| 1100.0 | Uninterrupted laser-on duration (s) | 22837 | 1147.73 | 645.35 | 73.8645 | 492.97 | 692.761 | 883.406 | 1583.05 | 2429.55 | 2822.76 |
| 1100.0 | Total laser-on within an encounter (s) | 45332 | 577.696 | 733.163 | 0 | 0 | 0 | 356.91 | 886.405 | 2257.51 | 2817.36 |

## 3) Distribution by Helper number

### Coverage and case totals

| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 18 | 40 | 18 | 36 | 18 | 0 | 7650.3 |
| 100 | 18 | 1931 | 29 | 1919 | 27 | 0 | 332681 |
| 200 | 18 | 4018 | 60 | 2983 | 29 | 0 | 474998 |
| 300 | 18 | 6188 | 80 | 4059 | 31 | 0 | 619748 |
| 400 | 18 | 7932 | 96 | 5207 | 32 | 0 | 630712 |
| 500 | 18 | 10232 | 129 | 5546 | 33 | 0 | 630712 |
| 600 | 18 | 12190 | 150 | 4758 | 33 | 0 | 630712 |
| 700 | 18 | 14502 | 171 | 4992 | 33 | 0 | 630712 |
| 800 | 18 | 16386 | 209 | 4882 | 33 | 0 | 630712 |

### Encounter speed

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Encounter speed: sampled minimum (m/s) | 40 | 64.7214 | 57.539 | 1.20376 | 1.54318 | 25.952 | 49.4273 | 81.3593 | 192.48 | 236.45 |
| 1 | Encounter speed: sampled maximum (m/s) | 40 | 189.283 | 19.9451 | 128.942 | 168.554 | 178.415 | 188.176 | 200.007 | 223.548 | 237.546 |
| 100 | Encounter speed: sampled minimum (m/s) | 1931 | 63.6331 | 33.0306 | 1.20376 | 23.2588 | 42.6443 | 57.8769 | 76.9789 | 128.751 | 236.45 |
| 100 | Encounter speed: sampled maximum (m/s) | 1931 | 193.733 | 17.824 | 128.942 | 168.899 | 179.679 | 193.028 | 205.206 | 226.185 | 244.54 |
| 200 | Encounter speed: sampled minimum (m/s) | 4018 | 65.1236 | 37.2338 | 0.585455 | 20.4448 | 42.2897 | 57.7451 | 77.0073 | 134.402 | 237.312 |
| 200 | Encounter speed: sampled maximum (m/s) | 4018 | 194.127 | 18.0806 | 103.41 | 168.433 | 180.127 | 193.237 | 205.759 | 226.623 | 243.37 |
| 300 | Encounter speed: sampled minimum (m/s) | 6188 | 66.9551 | 38.9634 | 0.472259 | 19.9723 | 42.9848 | 58.4416 | 78.6835 | 153.507 | 234.307 |
| 300 | Encounter speed: sampled maximum (m/s) | 6188 | 194.475 | 18.0987 | 26.3034 | 168.444 | 180.122 | 193.949 | 206.735 | 224.994 | 244.415 |
| 400 | Encounter speed: sampled minimum (m/s) | 7932 | 62.7025 | 32.9939 | 0.0262687 | 19.8469 | 42.028 | 57.2865 | 75.6612 | 128.741 | 236.291 |
| 400 | Encounter speed: sampled maximum (m/s) | 7932 | 194.048 | 18.1464 | 15.0291 | 168.662 | 179.903 | 192.936 | 206.205 | 226.624 | 244.318 |
| 500 | Encounter speed: sampled minimum (m/s) | 10232 | 66.7278 | 39.9981 | 0.684772 | 20.3091 | 42.618 | 57.8416 | 77.9927 | 140.833 | 236.688 |
| 500 | Encounter speed: sampled maximum (m/s) | 10232 | 194.414 | 18.472 | 36.5834 | 168.437 | 180.209 | 193.59 | 206.647 | 226.838 | 244.091 |
| 600 | Encounter speed: sampled minimum (m/s) | 12190 | 64.9145 | 36.2618 | 0.121126 | 20.4796 | 42.4694 | 57.8456 | 77.2 | 135.099 | 235.788 |
| 600 | Encounter speed: sampled maximum (m/s) | 12190 | 194.049 | 18.161 | 26.3034 | 168.153 | 179.942 | 193.366 | 206.152 | 225.822 | 244.822 |
| 700 | Encounter speed: sampled minimum (m/s) | 14502 | 68.6601 | 43.8129 | 0.026899 | 20.3745 | 42.9035 | 58.5337 | 79.025 | 175.727 | 237.867 |
| 700 | Encounter speed: sampled maximum (m/s) | 14502 | 194.815 | 18.7454 | 4.3819 | 168.435 | 180.183 | 193.9 | 207.341 | 229.039 | 244.95 |
| 800 | Encounter speed: sampled minimum (m/s) | 16386 | 66.1651 | 38.5892 | 0.0536914 | 20.2797 | 42.6602 | 58.0648 | 77.8632 | 138.011 | 236.291 |
| 800 | Encounter speed: sampled maximum (m/s) | 16386 | 194.388 | 18.3458 | 15.0291 | 168.441 | 180.077 | 193.429 | 207.134 | 226.096 | 244.701 |

### Signed range rate

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Signed range rate: sampled minimum (m/s) | 40 | -93.2842 | 40.2704 | -154.141 | -135.962 | -126.612 | -113.812 | -68.8672 | -6.81783 | -4.52961 |
| 1 | Signed range rate: sampled maximum (m/s) | 40 | 94.8293 | 38.8824 | 4.55179 | 6.90944 | 72.72 | 108.49 | 127.875 | 135.411 | 154.345 |
| 100 | Signed range rate: sampled minimum (m/s) | 1931 | -111.309 | 28.773 | -171.2 | -150.66 | -133.393 | -116.239 | -89.3911 | -64.906 | -4.11257 |
| 100 | Signed range rate: sampled maximum (m/s) | 1931 | 109.887 | 31.9362 | -107.791 | 59.1783 | 87.7295 | 116.408 | 132.584 | 151.101 | 175.634 |
| 200 | Signed range rate: sampled minimum (m/s) | 4018 | -109.038 | 32.3785 | -171.187 | -150.14 | -133.211 | -115.019 | -86.4622 | -55.3108 | 109.267 |
| 200 | Signed range rate: sampled maximum (m/s) | 4018 | 108.853 | 33.0466 | -128.983 | 53.4918 | 86.8678 | 115.34 | 132.566 | 151.952 | 175.754 |
| 300 | Signed range rate: sampled minimum (m/s) | 6188 | -108.274 | 32.4382 | -172.882 | -150.358 | -133.053 | -114.981 | -87.1893 | -45.3262 | 74.2421 |
| 300 | Signed range rate: sampled maximum (m/s) | 6188 | 108.309 | 33.0419 | -113.134 | 45.487 | 85.6362 | 114.71 | 132.878 | 152.356 | 176.284 |
| 400 | Signed range rate: sampled minimum (m/s) | 7932 | -110.459 | 30.3627 | -172.276 | -150.419 | -133.374 | -115.741 | -89.3201 | -61.8236 | 109.267 |
| 400 | Signed range rate: sampled maximum (m/s) | 7932 | 110.481 | 31.2678 | -137.073 | 58.2066 | 88.9405 | 116.116 | 133.179 | 152.202 | 177.229 |
| 500 | Signed range rate: sampled minimum (m/s) | 10232 | -108.068 | 32.7907 | -171.954 | -150.121 | -132.743 | -114.512 | -87.3845 | -44.6337 | 125.675 |
| 500 | Signed range rate: sampled maximum (m/s) | 10232 | 108.076 | 33.4085 | -123.585 | 44.4997 | 86.2064 | 114.619 | 132.476 | 151.86 | 176.901 |
| 600 | Signed range rate: sampled minimum (m/s) | 12190 | -109.324 | 31.5812 | -172.494 | -150.371 | -133.28 | -115.319 | -87.5679 | -52.7398 | 109.267 |
| 600 | Signed range rate: sampled maximum (m/s) | 12190 | 109.283 | 32.3746 | -133.028 | 51.5816 | 87.0639 | 115.274 | 133.229 | 152.357 | 176.737 |
| 700 | Signed range rate: sampled minimum (m/s) | 14502 | -106.9 | 35.1211 | -171.914 | -150.223 | -132.899 | -114.366 | -86.4322 | -30.122 | 122.07 |
| 700 | Signed range rate: sampled maximum (m/s) | 14502 | 106.878 | 35.7869 | -124.846 | 30.0595 | 85.5166 | 114.435 | 132.669 | 152.071 | 176.12 |
| 800 | Signed range rate: sampled minimum (m/s) | 16386 | -108.765 | 32.0972 | -172.024 | -150.213 | -133.057 | -114.839 | -87.7502 | -44.7389 | 128.321 |
| 800 | Signed range rate: sampled maximum (m/s) | 16386 | 108.78 | 32.7966 | -134.126 | 44.7873 | 87.0297 | 114.788 | 132.845 | 152.176 | 177.348 |

### Encounter duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Encounter duration (s) | 40 | 57360.8 | 165326 | 170.911 | 1613.4 | 1964.15 | 2332.58 | 4712.73 | 630712 | 630712 |
| 100 | Encounter duration (s) | 1931 | 4525.77 | 25486 | 104.418 | 1667.91 | 2014.2 | 2348.56 | 4533.77 | 5184.81 | 630712 |
| 200 | Encounter duration (s) | 4018 | 4246.7 | 18574 | 0.749804 | 1616.06 | 1997.73 | 2332.45 | 4357.37 | 5361.75 | 630712 |
| 300 | Encounter duration (s) | 6188 | 4336.45 | 22291 | 129.907 | 1609.37 | 1983.88 | 2326.41 | 4319.47 | 5372.97 | 630712 |
| 400 | Encounter duration (s) | 7932 | 4320.03 | 21759.5 | 100.801 | 1640.31 | 2001.18 | 2336.82 | 4399.55 | 5389.7 | 630712 |
| 500 | Encounter duration (s) | 10232 | 4275.92 | 21270.5 | 62.1945 | 1576.55 | 1989.37 | 2331.03 | 4306.17 | 5354.53 | 630712 |
| 600 | Encounter duration (s) | 12190 | 4286.75 | 21469.4 | 40.6091 | 1624.86 | 1987.97 | 2328.47 | 4334.62 | 5284.71 | 630712 |
| 700 | Encounter duration (s) | 14502 | 4209.97 | 21513.5 | 36.4059 | 1485.31 | 1956.15 | 2309.15 | 4314.52 | 5286.34 | 630712 |
| 800 | Encounter duration (s) | 16386 | 4280.93 | 22035.1 | 35.0998 | 1607.56 | 1980.72 | 2319.13 | 4301.94 | 5279.33 | 630712 |

### Laser-on duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Uninterrupted laser-on duration (s) | 36 | 59198.4 | 172989 | 822.578 | 876.839 | 1532.67 | 2168.9 | 4385.8 | 630712 | 630712 |
| 1 | Total laser-on within an encounter (s) | 40 | 53277 | 165070 | 0 | 0 | 1108.05 | 2000.01 | 4226.45 | 630712 | 630712 |
| 100 | Uninterrupted laser-on duration (s) | 1919 | 3544.91 | 25081.3 | 103.26 | 869.438 | 1158.24 | 2027.09 | 2611.27 | 5138.54 | 630712 |
| 100 | Total laser-on within an encounter (s) | 1931 | 3520.88 | 25004.8 | 0 | 856.468 | 1152.29 | 2017.04 | 2598.47 | 5134.38 | 630712 |
| 200 | Uninterrupted laser-on duration (s) | 2983 | 2988.22 | 20202 | 102.263 | 816.582 | 1060.08 | 1799.24 | 2387.55 | 5119.85 | 630712 |
| 200 | Total laser-on within an encounter (s) | 4018 | 2217.57 | 17455.7 | 0 | 0 | 0 | 1145.05 | 2178.26 | 4948.17 | 630712 |
| 300 | Uninterrupted laser-on duration (s) | 4059 | 2570.68 | 17344.6 | 160.188 | 793.653 | 968.032 | 1415.47 | 2137.28 | 5329.03 | 630712 |
| 300 | Total laser-on within an encounter (s) | 6188 | 1685.62 | 14100.4 | 0 | 0 | 0 | 960.33 | 1766.04 | 5173.87 | 630712 |
| 400 | Uninterrupted laser-on duration (s) | 5207 | 2083.83 | 15304.3 | 98.0739 | 644.756 | 744.727 | 1017.02 | 1844.36 | 5279.67 | 630712 |
| 400 | Total laser-on within an encounter (s) | 7932 | 1367.39 | 12439.2 | 0 | 0 | 0 | 740.345 | 1331.92 | 4747.69 | 630712 |
| 500 | Uninterrupted laser-on duration (s) | 5546 | 1970.59 | 14858.3 | 7.18036 | 515.861 | 589.684 | 693.139 | 2003.46 | 5182.38 | 630712 |
| 500 | Total laser-on within an encounter (s) | 10232 | 1067.68 | 10982.9 | 0 | 0 | 0 | 527.754 | 732.747 | 2930.77 | 630712 |
| 600 | Uninterrupted laser-on duration (s) | 4758 | 2306.95 | 16028.4 | 94.4479 | 449.61 | 891.215 | 1179.05 | 2088.68 | 5084.41 | 630712 |
| 600 | Total laser-on within an encounter (s) | 12190 | 900.111 | 10076.8 | 0 | 0 | 0 | 0 | 993.593 | 2606.93 | 630712 |
| 700 | Uninterrupted laser-on duration (s) | 4992 | 2234.17 | 15644 | 76.5664 | 426.983 | 828.05 | 1151.51 | 1986.35 | 5190.15 | 630712 |
| 700 | Total laser-on within an encounter (s) | 14502 | 768.758 | 9239.61 | 0 | 0 | 0 | 0 | 841.952 | 2367.9 | 630712 |
| 800 | Uninterrupted laser-on duration (s) | 4882 | 2292.01 | 15825.6 | 73.8645 | 614.47 | 733.316 | 1211.52 | 2147.28 | 5171.49 | 630712 |
| 800 | Total laser-on within an encounter (s) | 16386 | 682.601 | 8701.5 | 0 | 0 | 0 | 0 | 691.887 | 2416.94 | 630712 |

## 4) Distribution by Scheduler

### Coverage and case totals

| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| naive_next_entering | 81 | 35138 | 446 | 8547 | 126 | 0 | 573223 |
| positive_along_track | 81 | 38281 | 496 | 25835 | 143 | 0 | 630712 |

### Encounter speed

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| naive_next_entering | Encounter speed: sampled minimum (m/s) | 35138 | 67.7721 | 40.6716 | 0.182369 | 22.8033 | 43.1829 | 58.1712 | 78.5206 | 154.744 | 237.867 |
| naive_next_entering | Encounter speed: sampled maximum (m/s) | 35138 | 195.072 | 18.3324 | 4.3819 | 169.193 | 180.696 | 194.363 | 207.396 | 226.852 | 244.95 |
| positive_along_track | Encounter speed: sampled minimum (m/s) | 38281 | 64.559 | 36.9941 | 0.0262687 | 17.1729 | 42.1021 | 57.8741 | 77.0813 | 134.318 | 236.45 |
| positive_along_track | Encounter speed: sampled maximum (m/s) | 38281 | 193.698 | 18.3343 | 113.08 | 167.781 | 179.4 | 192.578 | 206.027 | 226.674 | 244.54 |

### Signed range rate

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| naive_next_entering | Signed range rate: sampled minimum (m/s) | 35138 | -107.526 | 32.5543 | -170.147 | -149.327 | -131.766 | -113.99 | -85.5966 | -45.4799 | 128.321 |
| naive_next_entering | Signed range rate: sampled maximum (m/s) | 35138 | 107.385 | 33.5375 | -137.073 | 45.4383 | 84.4771 | 113.941 | 131.359 | 150.864 | 174.893 |
| positive_along_track | Signed range rate: sampled minimum (m/s) | 38281 | -109.599 | 32.4924 | -172.882 | -151.384 | -134.166 | -115.841 | -89.0307 | -47.7128 | 128.321 |
| positive_along_track | Signed range rate: sampled maximum (m/s) | 38281 | 109.637 | 33.0443 | -133.269 | 48.7588 | 88.2583 | 115.882 | 134.181 | 152.959 | 177.348 |

### Encounter duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| naive_next_entering | Encounter duration (s) | 35138 | 4588.68 | 30884 | 36.4059 | 1607.98 | 2028.74 | 2350.2 | 4661.26 | 5333.43 | 630712 |
| positive_along_track | Encounter duration (s) | 38281 | 4053.38 | 7358.04 | 0.749804 | 1575.33 | 1943.03 | 2281.91 | 4178.28 | 5132.15 | 136298 |

### Laser-on duration

| Group | Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| naive_next_entering | Uninterrupted laser-on duration (s) | 8547 | 4988.93 | 35249.1 | 73.8645 | 1728.03 | 2047.81 | 2342.48 | 4591.86 | 5301.45 | 630712 |
| naive_next_entering | Total laser-on within an encounter (s) | 35138 | 1212.58 | 17515.9 | 0 | 0 | 0 | 0 | 0 | 4817.97 | 630712 |
| positive_along_track | Uninterrupted laser-on duration (s) | 25835 | 1576.93 | 3177.93 | 7.18036 | 502.936 | 716.011 | 938.138 | 1468.39 | 2764.54 | 61344.8 |
| positive_along_track | Total laser-on within an encounter (s) | 38281 | 1064.23 | 2713.19 | 0 | 0 | 0 | 722.926 | 1118.5 | 2400.95 | 61344.8 |

## Reproducibility and exported data

- [Encounter events](target_helper_encounters.csv): one row per target-helper encounter, with clipping flags, sampled extrema, duration, and accumulated within-encounter laser-on.
- [Laser-on events](target_helper_laser_intervals.csv): one row per uninterrupted firing interval, with within/outside-range durations.
- [Case summary and provenance](case_summary.csv): all 162 cases, including zero-event cases, totals, source directories, and source checksums.
- [Distributions](distributions.csv): all numerical distribution tables, plus unclipped-only sensitivity results. `contributing_cases`, `clipped_events` and `zero_events` accompany each distribution.
- [Evaluation code](../../Code/evaluate_encounters.jl); [tests](../../../Tests/absolute_RN/test_encounter_evaluation.jl).

From the project root:

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/evaluate_encounters.jl
```

The pipeline writes event CSVs and case summaries, reloads them to compute distributions, then reloads the distribution/case CSVs to generate this Markdown. No simulation files are modified. Use `--report-only` to regenerate Markdown from existing CSVs; `--input PATH` and `--output PATH` override the campaign and report directories.
