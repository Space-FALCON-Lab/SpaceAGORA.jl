# V5 Laser Interlinks: Implementation and Validation

Date: 2026-09-28. Runtime: Julia 1.12.1 on Linux.

Follow-up, 2026-10-03: this copy now includes array-backed graph updates and
exact iterative tree matching. See [the implementation plan](plan.md) and
[the 4,000-satellite validation results](scalability-validation.md). The original
143-assertion physics suite still passes; the expanded suite passes 407 assertions.

## Run

From the workspace root:

```sh
SPACEAGORA_RHS_CALIBRATE=off julia --project='8_SpaceAGORA.jl-main_v5 copy' --startup-file=no '8_SpaceAGORA.jl-main_v5 copy/test/runtests.jl'
julia --project='8_SpaceAGORA.jl-main_v5 copy' --startup-file=no '8_SpaceAGORA.jl-main_v5 copy/II_examples/gve_sma_interlinks.jl'
```

The first command runs the regressions and the one-hour comparison. The second
runs just the two example trajectories. No new dependencies, SPICE kernels,
GRAM assets, browser, or external services are required.

## Implementation Choices

- `SpacecraftModel` defaults to `n_terminal=1`, `battery_energy_index=100.0`, and
  `tempurature_index=100.0`. Zero terminals are permitted; indices must lie in
  `[0, 100]`. The battery energy index is normalized, not energy in joules;
  no battery depletion dynamics are implemented. The plan's
  `tempurature_index` spelling is retained.
- `InterLinkModel.linkgraph` has typed keys `((satellite, terminal), (satellite,
  terminal))`. Satellite numbers are **1-based spacecraft vector indices**, not
  persistent `SpacecraftModel.id` values. IDs are not reassigned.
- Register each candidate with `register_candidate!(model, spacecraft,
  first_endpoint, second_endpoint; parameters=InterLinkParameters(...))`.
  Endpoints are canonicalized as one undirected physical connection. Reverse
  registration is rejected, and registration leaves all flags false and score
  zero. Distinct terminals can connect a spacecraft to multiple partners.
- `InterLinkParameters` contains only `P` (watts), `B`, and `range` (meters).
  Defaults are 10,000 W, 100, and 200,000 m. Mutable scheduling state is separate.
- Eligibility checks range, terminal bounds, active spacecraft, battery energy,
  and forbidden satellite pairs. `battery_energy_threshold=50.0` by default.
  `tempurature_threshold=0.0` leaves temperature unrestricted by default.
  Optional `eligibility(key, spacecraft, state)::Bool` handles explicitly
  configured additional restrictions and should be read-only. Previous terminal
  occupancy does not affect availability.
- `SchedulingPolicyModel(:gve_sma; target_idx=1)` is a separate simulation-level
  setting. It **maximizes the chosen target's instantaneous semimajor-axis
  increase**, in m/s, not the sum of both endpoints. `target_idx` is a spacecraft
  vector index and defaults to 1. Connections not involving the target score zero.
  It is not a deorbit objective. Both endpoints still receive their physical
  forces. Other policy symbols are currently rejected.
- The score uses the energy-equivalent Gauss equation at the target:
  `da/dt = (2 * a^2 / mu) * dot(velocity, force) / mass`.
  It uses the configured planet's gravitational parameter and the current stage
  mass, without multiplying by a guessed future step duration.
- The optical law is `F = B * P / c`, directed away from the partner, matching
  the PDF Case 1 convention and the legacy run with `eta=1`, `beta=1`.
  Each endpoint force is evaluated separately from the same
  stage snapshot and divided by that endpoint's own mass. Equal and opposite
  forces here follow from this chosen symmetric law, not from an assumption
  that all interlinks must behave that way.
- `InterLinkModel.active_link_penalty` defaults to zero and has score units
  (m/s for `gve_sma`). It is charged for every selected physical connection.
  Exact component matching finds the maximum penalized total on general,
  including non-bipartite, graphs. Trees use iterative DFS and linear-time
  dynamic programming; cyclic components use memoized terminal-subset matching.
  Nonpositive net benefits are omitted; sorted traversal makes ties deterministic.
- Cyclic matching still has exponential worst-case cost. Large dense terminal
  graphs would need a scalable exact matching solver. Large trees, including
  one-terminal target stars, no longer use exponential subset matching.
- Attach the graph and policy using `SimulationConfiguration(interlink_model=...,
  scheduling_policy_model=..., ...)`. Scheduling runs at initialization and
  after every accepted step. The callback refreshes the derivative cache.
  Selection remains fixed during trial stages; geometry is always current.
- Continuous force is added in the existing coupled RHS. There are no velocity
  kicks and no scheduling or diagnostic accumulation inside RHS evaluations.
  `laser_dv` (integrated inertial acceleration, m/s) and `laser_delta_sma`
  (integrated laser-induced semimajor-axis rate, m) are additional solver state
  components and standard saved fields when interlinks are enabled.
  `model.history` records initial/accepted-step selections, not output times.
- The example defaults to `dt_max_orbit=10.0`; the comparison uses 5.0 at identical
  `reltol_orbit=1e-9`, `abstol_orbit=1e-11`. Existing non-interlink defaults are
  unchanged. Range/eligibility crossings are acted on at accepted endpoints, not
  localized by exact events.
- Radio can be represented but its force physics is not implemented. Split,
  multirate, and gravity-backbone solver modes explicitly reject interlinks.
  Validation covers the full-state `:tsit5` solver. Cross-spacecraft coupling
  disables the engine's otherwise block-diagonal Jacobian assumption.

## Verification Results

**PASS: 143/143 assertions; command exit code 0.** No skipped or failing tests.
The v5 snapshot had no pre-existing test suite; one focused suite was added.
Package loading and a 20-second engine run were verified before the full run.
The user's pre-existing deletions and reconstruction edits were preserved.

Coverage includes initialization, invalid registration, range, battery energy,
temperature and forbidden-pair restrictions, multiple terminals, repeated penalties, empty
selection, a greedy-failure example, exhaustive subset comparisons on eight
small general graphs, mass/power dependence, evaluation-order independence,
finite-difference orbital-rate checks, and changing stage geometry. Target-only
scoring is checked against the target's integrated derivative, including a
target at the second endpoint, negative target benefit, and invalid indices.

An analytic two-second switched-force case verifies continuous position changes
and derivative-cache refresh. Deliberately rejected solver trials verify that
history contains exactly one initial entry plus one entry per accepted step,
and diagnostics agree with a small-step reference. The normal no-interlink
engine path is also exercised. CSV diagnostics are read back and checked against
the final integrated state.

### One-Hour gve_sma Case

Three spacecraft at 1,000/1,050/1,050 km altitude, masses 227/227/454 kg,
with phases 0/+0.018/-0.018 rad. Two candidate connections compete for satellite
1's single terminal. Power is 10 kW, magnification 100, and range 200 km.
The reference run uses the same initial conditions with gravity only.

| Measurement | 10 s cap | 5 s cap |
| --- | ---: | ---: |
| Accepted steps | 370 | 728 |
| Maximum / median accepted interval | 10 / 10 s | 5 / 5 s |
| Selection switches | 2 | 2 |
| Satellite 1 laser delta SMA | +60.386128 m | +60.239849 m |
| Satellite 2 laser delta SMA | -37.302580 m | -37.303874 m |
| Satellite 3 laser delta SMA | -11.643647 m | -11.569613 m |
| Summed laser delta SMA | +11.439901 m | +11.366362 m |

The same switching sequence occurred in both runs. Switch times were
820.986254 and 1800.986254 s at the 10 s cap, versus 815.837709 and
1795.837709 s at the 5 s cap.

| Maximum difference between runs | Measured | Acceptance limit |
| --- | ---: | ---: |
| Per-satellite integrated laser delta SMA | 0.146279 m | 1 m |
| Integrated laser delta-v vector | 0.000138101 m/s | 0.001 m/s |
| Final position | 0.713949 m | 5 m |
| Corresponding switching time | 5.148544 s | 10 s |

In both runs, each satellite's final semimajor-axis difference from the
gravity-only reference agrees with its integrated laser diagnostic within
0.01 m. Both caps were reached, so this comparison genuinely tests different
scheduling resolutions. These tolerances validate this case, not all geometries
or future physical models.

Outputs are under `III_output/gve_sma_interlinks/dt_10.0s/` and
`III_output/gve_sma_interlinks/dt_5.0s/`. Each includes the normal trajectory
outputs with laser diagnostics and `accepted_step_schedule.csv` containing the
actual initial/accepted-step schedule times and active connections.

## PDF Case 1: Ten-Orbit Comparison

The requested duration is **10 initial target orbital periods**, not the PDF's
approximately 50 periods. Both implementations ran successfully for
63,713.403604987 s (17.698167668 hours). Each laser-on run has its own laser-off
reference with identical initial conditions and gravity.

### Settings and Scheduler Change

- Target: spacecraft vector index 1, altitude 1,050 km, initially circular and
  equatorial, true anomaly 0 degrees.
- Helpers: indices 2 through 21, altitude 1,000 km, initially circular and
  equatorial, evenly spaced at 0, 18, ..., 342 degrees.
- Mass: 227 kg each. Power: 10,000 W. Magnification: 100. Range: 200 km.
- Both use J2 gravity, no atmosphere or drag, the Tsit5 solver, a 10 s maximum
  step, and orbit relative/absolute tolerances of 1e-12. J2 matches the legacy
  ORACLE runner. The epoch is 2026-01-01.
- v5 now uses `SchedulingPolicyModel(:gve_sma; target_idx=1)`. Its score includes
  only the chosen target's instantaneous semimajor-axis rate. Both endpoints
  still receive physical forces. The legacy project also optimizes target 1.
- The legacy comparison loads `2_SpaceAGORA.jl` in a separate Julia process;
  no source files in that project were changed.

**Both runs now use `F = B * P / c`: 3.335640952 mN per endpoint**, matching
the PDF with beta and eta equal to 1. The v5 endpoint-force formula was updated
and its ten-orbit run regenerated; the already completed legacy run uses this
same force and was retained. v5 applies continuous stage forces; the legacy
runner applies endpoint velocity kicks, then updates its schedule.

### Results

Semimajor axis is the osculating two-body value computed from each saved
position and velocity using the configured planet's gravitational parameter.
It is sampled every 10 s plus the exact final time, giving 6,373 unique rows.

| Measurement | v5, target-only scheduler | Legacy, `2_SpaceAGORA.jl` |
| --- | ---: | ---: |
| Target initial SMA | 7,428,136.600000 m | 7,428,136.600000 m |
| Target final SMA | 7,428,262.949023 m | 7,428,262.833850 m |
| Target final-minus-initial SMA | +126.349023 m | +126.233850 m |
| Laser-off final-minus-initial SMA | +0.060099 m | +0.060099 m |
| Target gain relative to laser-off run | +126.288924 m | +126.173751 m |
| Accepted laser-on steps | 7,641 | 7,968 |
| Rejected laser-on steps | 0 | 0 |
| Solver return code | Success | Success |

The final v5-minus-legacy SMA difference is **0.115173 m**, or **0.091238%**
of the legacy target's final-minus-initial increase. The maximum absolute
difference over the saved history is 0.247768 m, and the history RMSE is
0.061835 m. The baseline-corrected gain ratio is 1.000912811.

These are close, but not identical, results at the tested settings. Continuous
versus impulsive integration and accepted-step scheduling remain numerical
differences. For example, the second substantial firing ends at 62,255.691 s
in v5 versus 62,251.591 s in the legacy run. Matching the force and objective
does not make the methods or accepted-step times identical. No ten-orbit
timestep-refinement study was performed, so this is not a claim of fully
converged or bitwise-equivalent trajectories.

The laser-off histories agree within 1.2853e-7 m, confirming agreement of the
shared baseline to well below a millimeter. Subtracting the laser-off reference
helps distinguish laser effects from J2 variation; it is not identical to v5's
integrated laser-only rate, because the laser changes the trajectory on which
J2 acts.

**Validation passed:** 11 assertions in each engine's run and 15 cross-run and
CSV readback assertions, in addition to the 143 scheduler/engine regressions.
Checks include the intended package identity, exact final time, common unique
sample grid, finite states, saved final SMA matching the solver state, matching
initial target SMA, gravitational parameter and force magnitude, baseline
agreement, and one schedule entry per initial/accepted step.

The supplied maximum-duration encounter CSV begins after 306,592 s, outside this
10-orbit run. It was left untouched and was not used as a whole-run reference;
these results do not validate the PDF's 50-orbit encounter statistics.

### Reproduce and Inspect

From the workspace root:

```sh
julia --project=8_SpaceAGORA.jl-main_v5 --startup-file=no 8_SpaceAGORA.jl-main_v5/test/compare_pdf_case1.jl v5 10
julia --project=2_SpaceAGORA.jl --startup-file=no 8_SpaceAGORA.jl-main_v5/test/compare_pdf_case1.jl legacy 10
```

The runner's default duration is also 10 orbits. An optional third argument
sets the maximum step in seconds. Sampling explicitly includes the endpoint,
so its `SavingCallback` uses `save_end=false` to avoid duplicate final rows.

Results are in `III_output/pdf_case1_comparison/`:

- [Side-by-side target history](III_output/pdf_case1_comparison/comparison_10orbits.csv).
  It joins the two target histories by their identical timestamps and includes
  each target's SMA, change from initial SMA, gain versus laser-off, and the
  v5-minus-legacy difference.
- [v5 target history](III_output/pdf_case1_comparison/v5_10.0orbits_dt10.0s_target.csv)
  and [legacy target history](III_output/pdf_case1_comparison/legacy_10.0orbits_dt10.0s_target.csv).
- Each run also saves `*_semimajor_axis.csv` for all 21 satellites,
  `*_summary.csv` for metadata and endpoint results, and the laser-on run's
  `*_schedule.csv` for actual accepted-step scheduling times.

The two commands regenerate the per-engine records. The side-by-side CSV was
derived from those records and verified by CSV readback; it is a comparison
artifact rather than an additional propagation.