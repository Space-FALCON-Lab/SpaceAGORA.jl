# Combined selected-altitude summaries

The following two tables pool all scenarios at target altitudes 1000, 1050, and 1100 km. The first retains full round-trip Float64 precision; the second is presentation-rounded.

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
