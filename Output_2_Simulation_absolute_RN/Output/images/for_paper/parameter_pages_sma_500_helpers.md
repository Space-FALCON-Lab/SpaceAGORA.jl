# 500 helper satellites: SMA parameters and SMA/INC laser-on statistics

**Scope:** The four original parameter tables use SMA (`gve_sma`) only; the additional within-encounter laser-on duration and percentage tables in section 4 cover both SMA (`gve_sma`) and INC (`gve_inc`) separately. All tables use exactly 500 helpers and target-helper pairs over 100 initial target periods. Initial target altitudes are 1000, 1050 and 1100 km, with initial target inclinations of 0.0, 0.5 and 1.0 deg. Initial helper altitude is 1000 km and inclination is 0.0 deg. Each Min/Avg/Max group summarizes one scenario, without pooling other helper counts or inclinations.

**Statistics:** Avg is the arithmetic event-weighted mean, not a within-encounter time average. Boundary-clipped events are included; their durations are observed portions, not necessarily complete lifetimes. For encounters already in range at simulation start, entry refers to the initial saved state. Values are rounded to six significant digits, except laser-on percentages, which use one decimal place; N/A means no observations.

## 1. Encounter speed at closest approach (m/s)

One relative-speed magnitude per geometric encounter, evaluated when helper-target separation is smallest within the observed encounter (`closest_relative_speed_mps`). This is not minimum relative speed or entry speed.

Closest approach is found from the replayed solver's dense states: compare encounter boundaries and accepted-step endpoints, and refine approaching-to-receding range-rate zero crossings within each accepted step (at most 10 s). The smallest separation over the encounter determines the reported speed. Boundary-clipped encounters use the closest point within the observed interval. Replay trajectories and scheduler histories match the original runs bit-for-bit; encounter boundary states are also checked.

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 0.38817 | 11.22 | 162.771 | 11.1017 | 31.1471 | 150.608 | 33.7862 | 54.0887 | 156.29 |
| 0.5 | 64.125 | 79.6525 | 207.788 | 26.0824 | 62.0826 | 163.645 | 37.8677 | 71.6085 | 168.847 |
| 1.0 | 128.24 | 188.775 | 237.552 | 112.499 | 135.961 | 197.66 | 46.3077 | 121.672 | 201.887 |

**Closest-approach data:** [Per-encounter times, separations and speeds](../../encounter_evaluation/sma500_closest_approach/encounters.csv); [full-precision summary and counts](../../encounter_evaluation/sma500_closest_approach/summary.csv). Signed range rate in parameter 2 remains measured at entry.

## 2. Signed range rate at entry (m/s)

One signed range-rate value per geometric encounter at its saved entry state (`entry_range_rate_mps`): negative means approaching, positive means receding. This is not the within-encounter range-rate envelope.

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | -46.8543 | -20.4893 | 5.42283e-13 | -114.008 | -86.1815 | 70.514 | -156.537 | -134.236 | 123.176 |
| 0.5 | -53.4075 | -20.8987 | 1.71835 | -117.035 | -83.7136 | 71.2304 | -156.993 | -129.179 | 123.801 |
| 1.0 | -72.3316 | -35.6619 | 6.86607 | -133.614 | -85.0415 | 73.3791 | -173.902 | -118.01 | 125.675 |

## 3. Encounter duration (s)

Observed time from entry to exit of the 200 km geometric encounter range, with one duration per target-helper encounter (`duration_s` in the encounter CSV).

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 1522.28 | 26444.3 | 102569 | 348.808 | 4323.58 | 5177.18 | 145.853 | 2108.71 | 2711.27 |
| 0.5 | 446.19 | 24889.5 | 104771 | 546.478 | 4235.4 | 5072.13 | 151.735 | 2041.73 | 2627.28 |
| 1.0 | 92.2622 | 10790.4 | 105134 | 460.351 | 3915.99 | 4677.29 | 134.862 | 1799.08 | 2538.98 |

## 4. Uninterrupted laser-on duration (s)

Duration of each uninterrupted firing interval, including held out-of-range firing (`duration_s` in the laser-interval CSV); not summed laser-on time per geometric encounter.

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 1515.88 | 6370.83 | 40490.1 | 343.459 | 1058.36 | 1708.86 | 144.34 | 585.603 | 790.16 |
| 0.5 | 4.93912 | 5682.09 | 38169.2 | 412.095 | 1067.23 | 1887.83 | 350.066 | 588.28 | 858.705 |
| 1.0 | 7.21679 | 2815.68 | 23157.5 | 527.241 | 1090.98 | 2320 | 127.827 | 595.716 | 1219.16 |

**Sources:** [Encounter records](../../encounter_evaluation/target_helper_encounters.csv) for parameters 2-3; [laser-interval records](../../encounter_evaluation/target_helper_laser_intervals.csv) for parameter 4. Records are filtered to `scheduler = gve_sma`, `helper_count = 500`, and the target altitude/inclination shown in each cell. Each scenario's event counts are checked against the [case summary](../../encounter_evaluation/case_summary.csv).

### Total laser-on duration inside each encounter (s)

**Definition:** One value per geometric target-helper encounter: `laser_on_within_encounter_s`, the sum of all firing segments inside its observed encounter interval. Out-of-range firing is excluded. This is total accumulated firing time, not the length of a single uninterrupted firing interval; separated firing segments can therefore produce a total longer than any individual interval.

**Scope and statistics:** SMA and INC are shown separately, using **only 500-helper cases**. Min/Avg/Max summarize all encounters in each altitude/inclination/scheduler scenario, **including zero-firing encounters**. Avg is the arithmetic mean of the per-encounter totals. Boundary-clipped encounters include only firing inside their observed durations. Values are in seconds, rounded to six significant digits; 0 means no in-range firing, while N/A would mean no encounters.

| Scheduler | Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SMA | 0.0 | 0 | 6370.83 | 40490.1 | 0 | 1053.11 | 1708.86 | 0 | 584.009 | 790.16 |
| SMA | 0.5 | 0 | 6183.29 | 39094.9 | 0 | 1061.89 | 1887.83 | 0 | 586.671 | 858.705 |
| SMA | 1.0 | 0 | 2946 | 47663.1 | 0 | 1085.41 | 2320 | 0 | 594.066 | 1219.16 |
| INC | 0.0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| INC | 0.5 | 0 | 5081.11 | 191186 | 0 | 628.348 | 1601.99 | 0 | 317.084 | 792.569 |
| INC | 1.0 | 0 | 1674.77 | 233920 | 0 | 627.246 | 1585.32 | 0 | 317.042 | 792.886 |

**Duration data:** Calculated from the [per-encounter records](../../encounter_evaluation/target_helper_encounters.csv) for the same 18 cases used in the percentage tables below. Event counts, zero-firing counts and summed in-range firing times are checked against the [case summary](../../encounter_evaluation/case_summary.csv).

### Laser-on duration percentage within each encounter

**Definition:** For each geometric target-helper encounter, the percentage is `100 * laser_on_within_encounter_s / duration_s`: total laser-on time inside that encounter, summed across all firing segments, divided by its observed encounter duration and multiplied by 100. Out-of-range firing is excluded. Values are percentages on a 0-100% scale (100% means firing throughout the observed encounter).

**Scope and statistics:** Both schedulers are shown separately, using **only 500-helper cases**. Each altitude/inclination/scheduler group summarizes one scenario. Min/Avg/Max are calculated from the individual encounter percentages, **including zero-firing encounters**. Avg is the arithmetic event-weighted mean of those percentages, not the ratio of summed durations or the ratio of the duration-table averages above. Boundary-clipped encounters use their observed durations. Values are rounded to one decimal place; N/A would mean no encounters. Rounded 0.0% or 100.0% values need not be exactly zero or full-duration firing.

#### SMA scheduler: laser-on duration percentage

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 0.0% | 25.3% | 99.6% | 0.0% | 24.6% | 98.5% | 0.0% | 27.9% | 99.0% |
| 0.5 | 0.0% | 28.9% | 99.9% | 0.0% | 25.3% | 99.8% | 0.0% | 28.9% | 99.0% |
| 1.0 | 0.0% | 41.3% | 100.0% | 0.0% | 27.8% | 98.9% | 0.0% | 33.1% | 94.8% |

#### INC scheduler: laser-on duration percentage

| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| 0.5 | 0.0% | 1.5% | 30.3% | 0.0% | 12.5% | 35.2% | 0.0% | 14.2% | 38.8% |
| 1.0 | 0.0% | 1.2% | 37.1% | 0.0% | 14.0% | 35.7% | 0.0% | 16.4% | 46.4% |

**Percentage data:** Calculated directly from the [per-encounter records](../../encounter_evaluation/target_helper_encounters.csv), filtered to `helper_count = 500` and the scheduler, altitude and inclination shown in each cell. Event and zero-firing counts are checked against the [case summary](../../encounter_evaluation/case_summary.csv): 10,376 encounters across 18 cases, including 1,873 zero-firing encounters. INC at 0.0 deg has no firing, so its encounter percentages are all 0.0%.
