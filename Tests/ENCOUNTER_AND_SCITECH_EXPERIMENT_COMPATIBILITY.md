# Combining the Encounter Scenarios and SciTech Experiments

## Summary

The experiments in the supplied *Typical Satellite Encounter Scenarios for
ORACLE* report and the SciTech presentation can use the **same simulation
framework**. Their analyses complement each other, but their parameter sets
and durations are not identical.

This note records the feasibility assessment, not a claim that the seven
encounter cases or the full presentation sweep have already been reproduced.
The presentation setup is documented in
[ORACLE_SCITECH_PRESENTATION_EXPERIMENT_SETUP.md](ORACLE_SCITECH_PRESENTATION_EXPERIMENT_SETUP.md).

## 1. How the seven encounter cases relate to the presentation

Table 3 on page 2 of the encounter report specifies the following cases.
Helpers are always on initially circular, equatorial orbits at **1000 km**,
and are evenly distributed. The target orbits are also initially circular.

| Case | Target altitude (km) | Target inclination (deg) | Initial target anomaly (deg) | Helpers | Scheduler |
|---|---:|---:|---:|---:|---|
| 1 | 1050 | 0 | 0 | 20 | SMA raising |
| 2 | 1150 | 0 | 0 | 20 | SMA raising |
| 3 | 1050 | 0 | 0 | 100 | SMA raising |
| 4 | 1000 | 5 | 0 | 20 | SMA raising |
| 5 | 1000 | 5 | 0 | 20 | Inclination raising |
| 6 | 1050 | 5 | 0 | 20 | SMA raising |
| 7 | 1000 | 0 | 1 | 20 | SMA raising |

These are representative cases within, or extending, the geometry study
described in the presentation:

- **Cases 1-3:** altitude separation and helper density.
- **Cases 4-6:** inclined encounters and the effect of scheduler choice.
  Their **5-degree inclination extends beyond the presentation's stated
  0-1.5-degree range**.
- **Case 7:** same-orbit interaction with an initial phase offset.

The encounter report uses approximately **50 orbits**, whereas the
presentation uses **80 orbits**. Both share:

| Physical parameter | Common value |
|---|---:|
| Maximum link range | 200 km |
| Input laser power | 10 kW |
| OCL amplification factor | 100 |
| Geometric transmission factor | 1 |
| Spacecraft mass | 227 kg |

Shared physical parameters do not make the trajectories identical. Initial
phasing, orbit geometry, duration, gravity settings, and numerical treatment
also matter.

## 2. Measure both sets of outputs

For each geometry and scheduler, the combined experiment should record:

| Presentation analysis | Encounter-report analysis |
|---|---|
| Cumulative radial, transverse, and normal velocity increments | Every feasible encounter's start, end, and duration |
| Changes in SMA, eccentricity, inclination, and RAAN | Laser-on intervals associated with each encounter |
| Dependence on helper population and geometry | Relative-speed and range-rate extrema |
| Comparison of SMA and inclination objectives | ECI trajectories for the longest encounter, sampled every 10 seconds |

This would help explain **why** a geometry produces a particular maneuver
outcome: more encounters, longer encounters, more useful force directions,
or better timing of activation.

### Encounters and laser activations are different

The encounter report defines an encounter by **range and line of sight**.
Laser firing additionally requires selection by the scheduler.

For example, Case 4 on page 6 reports **101 encounters but only 20 laser
firings**. Those are reported PDF results, not newly reproduced values.
Logging selected links alone would miss encounters during which no laser
fires.

Useful event records should therefore distinguish:

1. Pair availability intervals, independent of selection.
2. Selected laser-on intervals and the encounter to which each belongs.
3. Encounters with no laser activation.
4. Events truncated by the simulation start or end.

Relative speed and signed range rate are also distinct:

$$
v_{\mathrm{rel}}=\lVert\mathbf{v}_i-\mathbf{v}_j\rVert,
\qquad
\dot{\rho}=
\frac{(\mathbf{r}_i-\mathbf{r}_j)\cdot(\mathbf{v}_i-\mathbf{v}_j)}
{\lVert\mathbf{r}_i-\mathbf{r}_j\rVert}.
$$

The latter is undefined at zero separation, requiring explicit treatment
rather than an arbitrary numeric replacement.

## 3. Changes needed in the current comparison runners

v9 now supports both required schedulers, `:gve_sma` and `:gve_inc`, so the
core scoring capability is present. However, the current
[v9 comparison runner](run_9_full.jl) is not yet a complete runner for these
two studies.

The remaining work includes:

- **Configurable target altitude and initial phase**, rather than fixed
  1050 km altitude and zero-degree phase.
- **Encounter detection independent of scheduler selection**, including
  the report's range and line-of-sight definition. v9's current graph
  availability checks do not themselves include occultation.
- **Relative-speed and range-rate calculations**, event extrema, and
  longest-encounter trajectory exports.
- **A combined output set**, including the RTN and orbital-element
  diagnostics as well as encounter records.
- **Consistent environment and numerical settings**, with explicit
  gravity flags, duration conventions, and scheduler timing.

Our recent comparison runs disabled $J_2$. The presentation includes $J_2$
in its governing model and explicitly enables it in its verification plot.
The intended runtime flags for a reproduction must be established rather
than inherited silently from the recent runs.

## 4. Initialization and event-timing issues

### Potential initial overlap in Cases 4 and 5

These cases put the target and helpers at the same orbital radius, with
the target initially at zero true anomaly and five-degree inclination.
**If** a helper begins at the same ascending-node position, they initially
coincide despite their different inclinations.

The complete orientation and helper phasing must therefore be established.
An inclination difference alone does not guarantee nonzero spatial
separation. Any deliberate change to the PDF setup should be documented,
not silently introduced.

Case 7 explicitly uses a one-degree anomaly offset to avoid overlap.
That does not establish the unstated helper phase origin for every other
case.

### Ten-second output is not ten-second event accuracy

The encounter report requests ten-second sampling for the longest-encounter
trajectory. This is an **output cadence**, not necessarily the appropriate
resolution for determining event start and end times.

Accurate duration comparisons require finer event localization or a
demonstrated time-resolution error bound. In particular:

- Encounter boundaries and scheduler switch times are separate events.
- A scheduler updating at accepted solver steps need not switch exactly
  at a geometric range crossing.
- Interpolated trajectory samples should not be mistaken for an exact
  history of discrete scheduler decisions.

## 5. Recommended experiment organization

Use **one configurable experiment runner with two named case sets**:

1. **Encounter cases:** the seven documented geometries, approximately
   50 orbits, with their stated scheduler choices.
2. **Presentation sweep:** the documented altitude, inclination, and
   helper-count combinations, 80 orbits, with both scheduler objectives.

Apply the same encounter and maneuverability diagnostics to both.
Keep their case identifiers, initial conditions, duration conventions,
and output directories separate so comparisons remain traceable.

For an exact reproduction, obtain or resolve:

- Original initial states and helper phasing.
- Complete discrete sweep manifest.
- Environment flags and constants.
- Solver, tolerances, and scheduler update cadence.
- Event detection and boundary conventions.
- Angle wrapping and orbital-element post-processing.
- Producing code revision and original data.

**Conclusion:** combining the two studies is feasible and scientifically
useful. The shared framework can explain encounter opportunities and their
resulting maneuver effects together, but exact agreement with the PDFs'
numerical results is not guaranteed without the missing reproduction details.

## Related notes

- [Presentation experiment setup](ORACLE_SCITECH_PRESENTATION_EXPERIMENT_SETUP.md)
- [Inclination scheduler implementations](INCLINATION_SCHEDULER_COMPARISON.md)
- [Recent 50-helper SMA comparison](SMA_SCHEDULER_COMPARISON_50_HELPERS_50_ORBITS.md)
- [Recent 50-helper inclination comparison](INCLINATION_COMPARISON_50_HELPERS_50_ORBITS.md)
