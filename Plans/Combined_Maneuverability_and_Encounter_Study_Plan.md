# ORACLE: one simulation, two types of analysis

## 1. Goal and current status

Use [the new integration](<../9_SpaceAGORA.jl-main_full integration/>) to simulate
each grid case once, then reuse its saved data for orbital-element/RTN analysis
and the encounter report's analysis methods. No separate seven-case campaign
is needed.

The missing-force issue is fixed. Force, geometry, event, and mechanical-balance
checks passed in the [test harness](../Tests/STUDY_PREFLIGHT_CHECKS.md).
The [production runner](../Simulation/Code/README.md) implements the study's
configuration and recording pipeline. The
[production validation](../Tests/PRODUCTION_PIPELINE_VALIDATION.md) passed,
including three 100-orbit pilots and the 301-spacecraft benchmark.

### Separate positive/negative R/N campaign

The [absolute R/N runner](../Simulation_absolute_RN/Code/README.md) adds
solver-integrated positive and negative radial and normal laser delta-v, while
preserving the existing signed RTN, signed T and work diagnostics. For each
component, net = positive + negative and accumulated absolute = positive - negative.

Its separate grid has **270 fresh cases**: target altitudes 900, 950, 1000, 1050
and 1100 km; inclinations 0, 0.5 and 1 degree; helper counts 1, 100, 200, 300,
400, 500, 600, 700 and 800; and both schedulers. Helpers start at 1000 km.
Each case still uses 100 initial target periods, a 10-second step cap and 1001
trajectory rows. Row 801 is the 80-orbit endpoint; no separate 80-orbit runs
are needed.

Production results and controls are isolated in
[Simulation_absolute_RN/Output](../Simulation_absolute_RN/Output). Tests, pilots,
source snapshots and validation reports are under
[Tests/absolute_RN](../Tests/absolute_RN), not inside the new simulation folder.
The original results and code remain in place. The old 700-helper batch remains
stopped with four verified cases; this new campaign does not resume it.
The sections below describe the original study and its earlier extensions.

## 2. Experiment setup

| Setting | Value |
|---|---|
| Helpers | Initially circular at 1000 km, inclination 0 deg |
| Target | One spacecraft, initially circular |
| Target altitudes | 850, 900, 950, 1000, 1050, 1100, 1150 km |
| Target inclinations | 0, 0.5, 1.0 deg |
| Helper counts | 1, 50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800 at each target altitude |
| Initial phase | Target `nu = 1.0 deg`; helper k at `360 * (k - 1) / N` deg, k=1,...,N |
| Initial orientation | Same RAAN; argument of periapsis set to zero as a circular-orbit convention |
| Schedulers | `gve_sma` to raise SMA; `gve_inc` to raise inclination |
| Duration | `100 * initial_target_period` for each altitude: about 611,576-650,043 s (7.08-7.52 days) |
| Mass | 227 kg per spacecraft |
| Laser | P=10 kW, B=100, selection range 200 km; BP/c gives about 3.34 mN |
| Links | One terminal per spacecraft (`n_terminal=1`); only target-helper candidates |
| Scheduler timing | At initialization and every accepted integration step; `dt_max_orbit=10.0` |
| State/diagnostic saving | 1001 evenly spaced rows including exact start/end states |
| Encounter sampling | Every 10 s plus exact boundaries, independent of trajectory saving |
| Gravity | `InverseSquaredJ2GravityModel()` instead of central-gravity-only; equatorial reference radius; no drag |

Altitudes mean initial semimajor axis minus Earth's equatorial radius.
All orbits evolve under the applied forces; helpers are not held at constant
altitude. For exactly coplanar 0 deg cases, no inclination-raising firing is expected.

**Completed: 462 runs** at all seven target altitudes, with helpers at 1000 km
and eleven helper counts through 800 (excluding 700). The separate `helpers700`
campaign adds **42 runs**: 7 target altitudes x 3 inclinations x 2 schedulers.
It keeps helpers at 1000 km and does not rerun any completed cases.
Total planned coverage is **504 laser-on runs**.
New results use `N700/` within the existing
`Simulation/Output/paper_plot_mode/h1000km_t.../` layout. Their separate log is
`Simulation/Output/simulation_progress_helpers700.log`.
The mistaken helper-altitude batch is stopped; its retained data are excluded
from this study.
Optional laser-off references use identical initial conditions, duration, and
gravity/J2 settings; none are required for this extension.

### Accepted model assumptions

- Idealized constant laser force, with no pointing/acquisition losses.
  Hardware feasibility is outside scope.
- A selected link is held between scheduler updates. Delayed activation,
  shutoff, and helper reselection are accepted; a short encounter may receive
  no firing. No strict range cutoff is required.
- Earth occultation is nonlimiting within the stated altitude/range domain,
  as checked in the preflight report.

## 3. Saved outputs

Use one directory and unique run ID per case. Save numerical tables as Feather
and descriptive settings as TOML. CSV copies are optional; scenario progress and
elapsed wall time are recorded in the shared progress log.

| File | Contents |
|---|---|
| `trajectory.feather` | Time; every spacecraft's ECI position/velocity and existing cumulative ECI laser delta-v; the target RTN and spacecraft work columns defined below, at the same saved times |
| `encounters.feather` | Encounter ID, pair IDs, times, duration, boundary states/clipping flags, and sampled speed/range-rate/separation extrema with their times |
| `laser_firing.feather` | Firing ID, pair IDs, actual force-on/off times, duration, boundary-clipping flags, and power/force settings |
| `firing_range_audit.feather` | Post-processed within/outside-range segments of actual firing, as defined below |
| `longest_encounter.feather` | 10 s states plus boundaries for the longest all-pair and target-helper encounters; stored once if they are the same encounter |
| `run_info.toml` | Run ID, full setup, initial states, epoch, units, Earth/laser constants, solver settings, code version, checksums, and completion status |

There are **no separate duplicate RTN or mandatory scheduler-history files**.
The firing file is sufficient for actual firing history; raw selected-link
history is optional for debugging. Omit `scN_laser_delta_sma` from study output;
total SMA change is calculated from position/velocity afterward.

### RTN and work columns in the trajectory

The cumulative target R/T/N diagnostics are the following three solver-integrated
values, initialized to zero:

```
Rhat = r / norm(r)
Nhat = cross(r, v) / norm(cross(r, v))
That = cross(Nhat, Rhat)

d(target_dv_R)/dt = dot(a_laser, Rhat)
d(target_dv_T)/dt = dot(a_laser, That)
d(target_dv_N)/dt = dot(a_laser, Nhat)
```

Save them as `target_dv_R_mps`, `target_dv_T_mps`, and `target_dv_N_mps`.
Use the target's current state and actual applied laser acceleration, including
delayed shutoff but excluding gravity. These are measurements, not extra velocity
kicks. Existing `scN_laser_dv_1:3` columns are ECI integrals; rotating their final
vector cannot recover the RTN integrals because the frame changes over time.

Also accumulate positive and negative target along-track contributions separately
and each spacecraft's optical work, `integral(dot(F_laser, v) dt)`, in joules.
These share the trajectory's time column; they are not additional simulation runs.

### Encounter and firing rules

- Record all spacecraft-pair encounters, even without firing; also report the
  target-helper subset.
- Detect geometric entry/exit and save boundary states without rounding to the
  regular output times or changing scheduler behavior.
- Record nonzero force application, not merely helper selection. Preserve full
  firing intervals even when they extend outside a geometric encounter.
- Keep timestamps and IDs consistent across files; use typed empty tables when
  no encounters or firing occur.

### Auxiliary range audit

After simulation, split each firing interval at every recorded range crossing.
Each segment contains run/firing/pair IDs, nullable encounter ID, start/end time,
duration, and `within_range`. Segments must cover the entire firing interval
without gaps or overlaps, including firing wholly outside range.

Example: firing from 100-210 s with range exit at 203 s gives an inside segment
100-203 s and an outside segment 203-210 s.
This audits the actual trajectory; it does not change firing/RTN totals, add
missed opportunities, or predict a different scheduler's trajectory.

## 4. Post-processing

### Orbital elements

Calculate SMA, eccentricity, inclination, and RAAN from saved position/velocity,
using the simulation's Earth constants. Report **final minus initial**, including
gravity/J2 effects. Use the optional laser-off reference only when reporting a
matched laser-on minus laser-off comparison.

Handle angular wrapping and undefined classical elements: use orbit-plane and
eccentricity vectors where RAAN or periapsis is singular.

### Maneuverability and encounters

Use the saved cumulative RTN/work columns directly for net effects and
along-track cancellation; do not integrate them a second time.
Calculate firing duty cycle and helper usage from the firing intervals.

For every grid run, report encounter durations, relative speed, range rate
(negative approaching, positive separating), extrema, and the longest
encounter's ECI trajectory with boundary states.
Distinguish uninterrupted firing duration from total firing within an encounter:
the latter is the overlap of geometric and firing intervals.
Report outside-range firing separately from the range audit.

Use the saved encounter extrema and longest tracks, not the sparse 1001-row
trajectory, for encounter detail. Extrema are sampled estimates at 10 s plus boundaries.
If no relevant event exists, mark extrema/longest-encounter results unavailable.
Compare results by target altitude, inclination, helper count, and scheduler.
The existing plotting commands remain scoped to the completed 462 cases,
excluding the new 700-helper cases.

## 5. Production validation

The [production report](../Tests/PRODUCTION_PIPELINE_VALIDATION.md) records the
passed RTN/work, fresh-file, recorder-noninterference, mechanical-balance, and
encounter-product checks. Compact saving preserved the three pilots' final
states and firing schedules exactly. The outdated test include and serial
routing expectation were repaired.

The [preflight report](../Tests/STUDY_PREFLIGHT_CHECKS.md) preserves the passed
checks and measured sensitivity, not error bounds. The 10 s cap is settled;
further step-size comparisons are not a prerequisite for this study.
Saving at 10, 5, and 1 s left schedules, event
records, and shared-time values unchanged in both selected cases.
Sub-millisecond geometric recording consistency does not imply subsecond
physical switching accuracy. Old saved runs are comparison data, not exact
reference solutions; some lack complete settings and use a historical J2 radius.

## 6. Memory and execution

- Run cases sequentially in separate Julia processes; do not retain earlier
  trajectories in a grid-wide collection.
- Disable animation and full solver-step history. Keep only the fields above;
  make plots and CSV exports afterward.
- At 1001 rows, 301-spacecraft position/velocity occupies **14.46 MB raw**.
  Other diagnostics and longest-encounter tracks add storage; this is not peak RAM.
- The runner writes compressed Feather in chunks. A temporary dense-state spool
  (about 0.92 GB at 301 spacecraft) supplies encounter products and is deleted
  after they are saved. Budget for it during each run, not for the full grid at once.
- Use the SSH-independent launcher and check
  [simulation_progress_helpers700.log](../Simulation/Output/simulation_progress_helpers700.log)
  for the new campaign's scenario start/finish times, elapsed wall time,
  verification, and pauses. Earlier campaign logs remain separate.
- All-pair detection at this size involves 45,150 pairs, versus 300 target-helper
  pairs. Any restricted analysis must be labelled explicitly.

The batch reuses verified full-duration pilots, runs the remaining cases, and
verifies each saved output before moving on. Produce both analyses from these
files; no second simulation campaign is required.
