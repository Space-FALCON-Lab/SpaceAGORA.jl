# All-scenario encounter and laser-on summary

## Scope and definitions

- Covers **all 162 verified 100-orbit cases** at initial target altitudes 1000, 1050 and 1100 km: all three inclinations, all nine helper counts, and both schedulers.
- For each altitude, speed, signed range rate and encounter duration pool every target-helper encounter across both schedulers. Laser-on duration is the duration of each uninterrupted firing interval, tabulated separately for `gve_sma` and `gve_inc`.
- **Avg is event-weighted**: each encounter/firing interval is one observation. Scenarios producing more events contribute more; these are descriptive aggregates, not equal-weight averages of scenarios and not independent-replicate estimates.
- Speed and signed range rate are measured **at encounter entry** from saved start-boundary states. If already in range at simulation start, that initial state is used. Range rate is signed: negative means approaching; positive means receding.
- Encounter duration is observed time inside 200 km. Firing duration is uninterrupted selected-link firing and may include held firing outside that range. Boundary-clipped events remain included; maxima that equal a whole run's duration are observation-window bounds, not complete encounter lifetimes.
- Min, Avg and Max are across the observations in each column. `N/A` indicates no observations, not zero. Units are shown in each column header. The first combined table uses full round-trip Float64 values; the duplicate below is presentation-rounded. Per-altitude tables use six significant figures. CSVs retain full precision.

## Full-precision combined summary: all three target altitudes

This table pools every event at **1000, 1050, and 1100 km together**. Entry speed, entry range rate and encounter duration include both schedulers. Laser-on intervals remain separated by scheduler.

| Statistic | Encounter speed at entry (m/s) | Signed range rate at entry (m/s) | Encounter duration (s) | Laser-on time (gve_sma) (s) | Laser-on time (gve_inc) (s) |
|:---| ---: | ---: | ---: | ---: | ---: |
| Max | 245.52440607433067 | 128.32054709527063 | 630711.893783881 | 61346.08065121457 | 1598.3992165351083 |
| Avg | 181.8704668636592 | -105.00506369295432 | 4261.160985276496 | 1035.793961057239 | 481.41640923495714 |
| Min | 3.6652605947767847 | -174.08798961327517 | 4.433393856626935 | 0.37340983064862476 | 1.9362670129630715 |

### Presentation-rounded combined table

This repeats the combined statistics above, rounded to the requested precision: **5 significant figures** for encounter speed, signed range rate (assumed), and laser-on durations; **6 significant figures** for encounter duration.

| Statistic | Encounter speed at entry (m/s) | Signed range rate at entry (m/s) | Encounter duration (s) | Laser-on time (gve_sma) (s) | Laser-on time (gve_inc) (s) |
|:---| ---: | ---: | ---: | ---: | ---: |
| Max | 245.52 | 128.32 | 630712 | 61346 | 1598.4 |
| Avg | 181.87 | -105.01 | 4261.16 | 1035.8 | 481.42 |
| Min | 3.6653 | -174.09 | 4.43339 | 0.37341 | 1.9363 |

**Combined coverage:** 74130 Encounter speed at entry events; 74130 Signed range rate at entry events; 74130 Encounter duration events; 39328 Laser-on duration (gve_sma) events; 32677 Laser-on duration (gve_inc) events.
**Combined boundary flags:** 938 clipped Encounter speed at entry events (81 span the full run); 938 clipped Signed range rate at entry events (81 span the full run); 938 clipped Encounter duration events (81 span the full run); 144 clipped Laser-on duration (gve_sma) events (0 span the full run); 46 clipped Laser-on duration (gve_inc) events (0 span the full run).

## Target altitude: 1000 km

| Statistic | Encounter speed at entry (m/s) | Signed range rate at entry (m/s) | Encounter duration (s) | Laser-on time (gve_sma) (s) | Laser-on time (gve_inc) (s) |
|:---| ---: | ---: | ---: | ---: | ---: |
| Max | 240.58 | 67.5531 | 630712 | 61346.1 | 1591.8 |
| Avg | 194.875 | -26.9432 | 25768.5 | 4468.75 | 770.776 |
| Min | 3.66526 | -76.2345 | 4.43339 | 0.37341 | 2.01106 |

**Coverage:** 4371 Encounter speed at entry events; 4371 Signed range rate at entry events; 4371 Encounter duration events; 3026 Laser-on duration (gve_sma) events; 7418 Laser-on duration (gve_inc) events.
**Boundary flags:** 287 clipped Encounter speed at entry events (81 span the full run); 287 clipped Signed range rate at entry events (81 span the full run); 287 clipped Encounter duration events (81 span the full run); 48 clipped Laser-on duration (gve_sma) events (0 span the full run); 14 clipped Laser-on duration (gve_inc) events (0 span the full run).

## Target altitude: 1050 km

| Statistic | Encounter speed at entry (m/s) | Signed range rate at entry (m/s) | Encounter duration (s) | Laser-on time (gve_sma) (s) | Laser-on time (gve_inc) (s) |
|:---| ---: | ---: | ---: | ---: | ---: |
| Max | 245.524 | 74.2421 | 5457.04 | 3318.31 | 1598.4 |
| Avg | 191.701 | -80.6937 | 4509.42 | 1073.97 | 458.041 |
| Min | 25.048 | -134.956 | 152.9 | 28.0003 | 1.93627 |

**Coverage:** 23976 Encounter speed at entry events; 23976 Signed range rate at entry events; 23976 Encounter duration events; 12777 Laser-on duration (gve_sma) events; 10954 Laser-on duration (gve_inc) events.
**Boundary flags:** 348 clipped Encounter speed at entry events (0 span the full run); 348 clipped Signed range rate at entry events (0 span the full run); 348 clipped Encounter duration events (0 span the full run); 49 clipped Laser-on duration (gve_sma) events (0 span the full run); 14 clipped Laser-on duration (gve_inc) events (0 span the full run).

## Target altitude: 1100 km

| Statistic | Encounter speed at entry (m/s) | Signed range rate at entry (m/s) | Encounter duration (s) | Laser-on time (gve_sma) (s) | Laser-on time (gve_inc) (s) |
|:---| ---: | ---: | ---: | ---: | ---: |
| Max | 234.462 | 128.321 | 2825.65 | 1506.2 | 1597.04 |
| Avg | 175.481 | -125.189 | 2077.8 | 573.483 | 349.266 |
| Min | 49.4447 | -174.088 | 11.995 | 38.7416 | 2.24995 |

**Coverage:** 45783 Encounter speed at entry events; 45783 Signed range rate at entry events; 45783 Encounter duration events; 23525 Laser-on duration (gve_sma) events; 14305 Laser-on duration (gve_inc) events.
**Boundary flags:** 303 clipped Encounter speed at entry events (0 span the full run); 303 clipped Signed range rate at entry events (0 span the full run); 303 clipped Encounter duration events (0 span the full run); 47 clipped Laser-on duration (gve_sma) events (0 span the full run); 18 clipped Laser-on duration (gve_inc) events (0 span the full run).

## Interpretation

This is a compact description of the simulated orbital/encounter regimes, not a scheduler-effect estimate. The encounter columns mix scheduler-conditioned trajectories and weight them by event count. The scheduler split in laser-on time describes the resulting firing intervals, not control success. For claims about scheduler performance, compare matched cases and report per-case orbital-change/control outcomes.

Sources were the verified evaluator CSVs in `.`. [Full-precision statistics and counts](selected_altitude_statistics.csv); [combined three-altitude statistics](selected_altitude_combined_statistics.csv); [source file checksums](selected_altitude_sources.csv).

