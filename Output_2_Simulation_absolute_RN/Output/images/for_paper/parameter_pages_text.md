# Page 1. Encounter speed at entry

**Ordering assumption:** the request repeated encounter duration on pages 1 and 3; page 1 is interpreted as encounter speed to cover all four quantities.

**Definition:** One relative-speed value per geometric encounter, calculated from the saved entry velocities. This is not the minimum or maximum speed during the encounter.

**Scope:** target-helper pairs; 100 initial target periods; initial target h_t and i. Each Min/Avg/Max group pools all events over helper counts 1, 100, 200, 300, 400, 500, 600, 700 and 800 at that scheduler/h_t/i setting (nine cases). Boundary-clipped events are included. Avg is the arithmetic event-weighted mean, not an equal-case mean or a within-encounter time average. Values use m/s, rounded to six significant digits; N/A means no observations.

### SMA scheduler

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 3.66526 | 194.52 | 208.761 | 25.048 | 191.269 | 211.043 | 49.4447 | 174.724 | 195.124 |
| 0.5 | 64.2365 | 195.285 | 214.51 | 68.749 | 191.38 | 218.236 | 80.8089 | 174.826 | 204.686 |
| 1.0 | 128.315 | 196.478 | 240.58 | 130.247 | 189.834 | 244.553 | 108.06 | 172.53 | 232.754 |

### INC scheduler

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 3.66526 | 100.656 | 194.253 | 25.048 | 193.253 | 211.74 | 49.4447 | 177.922 | 195.921 |
| 0.5 | 64.2365 | 193.026 | 210.329 | 68.749 | 193.355 | 219.112 | 80.8089 | 177.966 | 206.23 |
| 1.0 | 128.315 | 195.849 | 238.481 | 130.473 | 191.489 | 245.524 | 109.618 | 175.162 | 234.462 |

## Conclusion

- Across all cells, the observed values span **3.66526 to 245.524 m/s**. One cell attaining the upper bound is **gve_inc, h_t = 1050 km, i = 1.0 deg**.
- The maximum observed entry speed increases from i = 0 to 0.5 to 1.0 deg in **6 of 6 scheduler/altitude groups**. This describes the upper bounds, not an increase in every encounter's speed.
- These bounds describe **speed at encounter entry**, not closest-approach speed, encounter-average speed, or the maximum speed reached later. For encounters already in range at simulation start, entry refers to the initial saved state.
- SMA and INC min-max ranges overlap in **9 of 9** matched altitude/inclination cells with observations from both. Min/max alone describe extremes, not typical behavior, statistical uncertainty, or scheduler superiority.

*Data: [full-precision table values and counts](../../encounter_evaluation/parameter_pages_min_max.csv); [source CSV checksums](../../encounter_evaluation/parameter_pages_sources.csv). Tables and conclusions are generated from the exported table CSV.*


---

# Page 2. Signed range rate at entry

**Definition:** One signed range-rate value per geometric encounter at its saved entry state: negative means approaching, positive means receding. This is not the within-encounter range-rate envelope.

**Scope:** target-helper pairs; 100 initial target periods; initial target h_t and i. Each Min/Avg/Max group pools all events over helper counts 1, 100, 200, 300, 400, 500, 600, 700 and 800 at that scheduler/h_t/i setting (nine cases). Boundary-clipped events are included. Avg is the arithmetic event-weighted mean, not an equal-case mean or a within-encounter time average. Values use m/s, rounded to six significant digits; N/A means no observations.

### SMA scheduler

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | -46.8543 | -19.9719 | 1.01286e-12 | -114.159 | -85.8162 | 71.6458 | -156.873 | -134.126 | 125.883 |
| 0.5 | -53.4075 | -20.6696 | 18.9772 | -118.042 | -83.3841 | 72.295 | -157.723 | -129.112 | 126.492 |
| 1.0 | -76.2345 | -35.7979 | 67.5531 | -134.956 | -84.989 | 74.2421 | -174.088 | -118.103 | 128.321 |

### INC scheduler

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | -7.86465e-13 | 6.38319e-15 | 1.01286e-12 | -104.174 | -75.2498 | 71.6458 | -151.221 | -128.958 | 125.883 |
| 0.5 | -11.7884 | -7.34381 | 18.9772 | -95.0671 | -73.5734 | 72.295 | -152.151 | -124.457 | 126.492 |
| 1.0 | -43.3444 | -29.2013 | 67.5531 | -116.29 | -79.2835 | 74.2421 | -170.412 | -115.856 | 128.321 |

## Conclusion

- Across all cells, the observed values span **-174.088 to 128.321 m/s**. One cell attaining the upper bound is **gve_sma, h_t = 1100 km, i = 1.0 deg**.
- There are **324 positive entry range-rate observations**, of which **324** are start-clipped. A positive value therefore must not automatically be interpreted as a normal inward range crossing; an encounter may already be in range at t = 0.
- SMA and INC min-max ranges overlap in **9 of 9** matched altitude/inclination cells with observations from both. Min/max alone describe extremes, not typical behavior, statistical uncertainty, or scheduler superiority.

*Data: [full-precision table values and counts](../../encounter_evaluation/parameter_pages_min_max.csv); [source CSV checksums](../../encounter_evaluation/parameter_pages_sources.csv). Tables and conclusions are generated from the exported table CSV.*


---

# Page 3. Encounter duration

**Definition:** Observed time from entry to exit of the 200 km geometric encounter range, with one duration per target-helper encounter.

**Scope:** target-helper pairs; 100 initial target periods; initial target h_t and i. Each Min/Avg/Max group pools all events over helper counts 1, 100, 200, 300, 400, 500, 600, 700 and 800 at that scheduler/h_t/i setting (nine cases). Boundary-clipped events are included. Avg is the arithmetic event-weighted mean, not an equal-case mean or a within-encounter time average. Values use s, rounded to six significant digits; N/A means no observations.

### SMA scheduler

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 1100.82 | 27305.8 | 116686 | 171.715 | 4351.09 | 5256.22 | 75.4108 | 2117.15 | 2785.6 |
| 0.5 | 129.064 | 26032.6 | 123078 | 166.316 | 4259.06 | 5172.5 | 74.573 | 2048.68 | 2712.37 |
| 1.0 | 92.2622 | 10940.4 | 114237 | 152.903 | 3931.39 | 4864.82 | 41.8082 | 1803.53 | 2650.28 |

### INC scheduler

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 630712 | 630712 | 630712 | 171.713 | 5138.32 | 5457.04 | 54.565 | 2318.34 | 2825.65 |
| 0.5 | 16.5328 | 66112.5 | 630712 | 166.312 | 5011.37 | 5410.93 | 11.995 | 2242.01 | 2815.5 |
| 1.0 | 4.43339 | 16657.2 | 630712 | 152.9 | 4508.96 | 5183.88 | 31.8475 | 1950.09 | 2781.89 |

## Conclusion

- Across all cells, the observed values span **4.43339 to 630712 s**. One cell attaining the upper bound is **gve_inc, h_t = 1000 km, i = 0.0 deg**.
- For **gve_inc, h_t = 1000 km, i = 0.0 deg**, all **32** saved encounters span the complete run. Its equal Min/Max values describe the observation window, not measured complete encounter lifetimes.
- **938 of 74130 encounters are boundary-clipped**, including **81 spanning the full run**. Long observed durations can be limited by the simulation window; they are not necessarily complete encounter lifetimes.
- SMA and INC min-max ranges overlap in **8 of 9** matched altitude/inclination cells with observations from both. Min/max alone describe extremes, not typical behavior, statistical uncertainty, or scheduler superiority.

*Data: [full-precision table values and counts](../../encounter_evaluation/parameter_pages_min_max.csv); [source CSV checksums](../../encounter_evaluation/parameter_pages_sources.csv). Tables and conclusions are generated from the exported table CSV.*


---

# Page 4. Uninterrupted laser-on duration

**Definition:** Duration of each uninterrupted firing interval, including held out-of-range firing. This matches the scheduler-separated laser plot, not summed laser-on time per geometric encounter.

**Scope:** target-helper pairs; 100 initial target periods; initial target h_t and i. Each Min/Avg/Max group pools all events over helper counts 1, 100, 200, 300, 400, 500, 600, 700 and 800 at that scheduler/h_t/i setting (nine cases). Boundary-clipped events are included. Avg is the arithmetic event-weighted mean, not an equal-case mean or a within-encounter time average. Values use s, rounded to six significant digits; N/A means no observations.

### SMA scheduler

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 668.764 | 6575.8 | 53895.9 | 28.0003 | 1071.28 | 3318.31 | 100.857 | 579.009 | 1498.19 |
| 0.5 | 0.37341 | 5687.28 | 56649.6 | 31.9279 | 1073.89 | 3292.63 | 93.4322 | 577.648 | 1504.93 |
| 1.0 | 2.23981 | 2887.09 | 61346.1 | 28.7316 | 1076.81 | 3278.44 | 38.7416 | 563.649 | 1506.2 |

### INC scheduler

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | N/A | N/A | N/A | N/A | N/A | N/A | N/A | N/A | N/A |
| 0.5 | 2.37622 | 768.965 | 1589.2 | 2.85106 | 462.348 | 1598.4 | 2.74885 | 350.845 | 1597.04 |
| 1.0 | 2.01106 | 772.593 | 1591.8 | 1.93627 | 453.772 | 1593.42 | 2.24995 | 347.681 | 1591.66 |

## Laser-on duration percentage within each encounter

**Definition:** For each geometric target-helper encounter, the percentage is `100 * laser_on_within_encounter_s / duration_s`: total laser-on time inside that encounter, summed across all firing segments, divided by its observed encounter duration and multiplied by 100. Out-of-range firing is excluded. Values are percentages on a 0-100% scale (100% means firing throughout the observed encounter).

**Scope and statistics:** The same altitude, inclination, scheduler and nine-helper-count grouping applies. Min/Avg/Max are calculated from the individual encounter percentages, **including zero-firing encounters**. Avg is the arithmetic event-weighted mean of those percentages, not the ratio of summed durations or the ratio of the duration-table averages above. Boundary-clipped encounters use their observed durations. Values are rounded to one decimal place; N/A would mean no encounters. Rounded 0.0% or 100.0% values need not be exactly zero or full-duration firing.

### SMA scheduler: laser-on duration percentage

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 0.0% | 24.6% | 100.0% | 0.0% | 24.5% | 99.9% | 0.0% | 27.3% | 99.9% |
| 0.5 | 0.0% | 27.1% | 100.0% | 0.0% | 25.1% | 99.8% | 0.0% | 28.1% | 99.8% |
| 1.0 | 0.0% | 40.4% | 100.0% | 0.0% | 27.1% | 100.0% | 0.0% | 31.1% | 99.6% |

### INC scheduler: laser-on duration percentage

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| 0.5 | 0.0% | 2.6% | 60.4% | 0.0% | 13.6% | 64.4% | 0.0% | 15.2% | 88.4% |
| 1.0 | 0.0% | 1.3% | 66.3% | 0.0% | 15.0% | 65.7% | 0.0% | 17.4% | 99.7% |

**Percentage data:** Calculated directly from the [per-encounter records](../../encounter_evaluation/target_helper_encounters.csv), with event and zero-firing counts checked against the [case summary](../../encounter_evaluation/case_summary.csv): 74,130 encounters across 162 cases, including 13,440 zero-firing encounters. INC at 0.0 deg therefore has 0% laser-on duration, even though its uninterrupted-firing statistics above are N/A because no firing intervals exist.

## Conclusion

- Across all cells, the observed values span **0.37341 to 61346.1 s**. One cell attaining the upper bound is **gve_sma, h_t = 1000 km, i = 1.0 deg**.
- At **h_t = 1000 km**, the longest observed uninterrupted firing is **61346.1 s for SMA**, versus **1591.8 s for INC**. SMA reaches much longer intervals here, but these maxima do not quantify how frequently they occur.
- **No firing intervals exist in 3 cells** (gve_inc, h_t = 1000 km, i = 0.0 deg; gve_inc, h_t = 1050 km, i = 0.0 deg; gve_inc, h_t = 1100 km, i = 0.0 deg), representing **27 cases**. Their Min/Max entries are N/A, not zero-duration firing intervals.
- Longer uninterrupted firing is not, by itself, evidence of better orbital control. Compare with total duty fraction and the scheduler's intended orbital-change objective.
- SMA and INC min-max ranges overlap in **6 of 6** matched altitude/inclination cells with observations from both. Min/max alone describe extremes, not typical behavior, statistical uncertainty, or scheduler superiority.

*Uninterrupted-duration data: [full-precision table values and counts](../../encounter_evaluation/parameter_pages_min_max.csv); [source CSV checksums](../../encounter_evaluation/parameter_pages_sources.csv). The uninterrupted-duration tables and conclusions are generated from the exported table CSV; the percentage tables use the per-encounter records linked above.*
