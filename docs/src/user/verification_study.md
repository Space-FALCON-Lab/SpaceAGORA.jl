# Verification Study

The telemetry verification study propagates repository-owned mission
scenarios and grades them against flight and reference data, with explicit
pass/fail thresholds. Use it to validate a local setup or compare outputs,
rather than a one-off example run.

Once the prerequisites below are in place, the shortest run is:

```text
julia --project=. benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true
```

What to read next:

- [Assets & Modes](../assets.md)
- [CLI](../cli.md)
- [Recipes](recipes.md)

## Prerequisites

The study loads the `GRAMSuite` package from the `data/GRAMSuite.jl` submodule
before it reads any scenario, so the submodule must be initialized even when
you select only `--scenarios=odyssey`. On a fresh clone without it the study
stops with "Package GRAMSuite not found in current path". Every scenario also
builds its planet from the SPICE kernels in
`data/GRAMSuite.jl/GRAM Suite 2.0/SPICE`, which the submodule stores in Git
LFS, so its LFS files must be pulled as well. The VEx and Earth GMAT scenarios
run on native GRAM (VenusGRAM and EarthGRAM) and therefore also need the
native GRAM library built from the licensed GRAM distribution; the Odyssey
scenario flies a measured density table instead.
[GRAMSuite Setup](gramsuite_setup.md) covers each step.

The truth files the study grades (Odyssey, VEx and GMAT Earth) and the
scenario definitions in `test/telemetry_benchmark_manifest.toml` are part of
the repository. Without `--scenarios`, the study runs every scenario in the
manifest.

## Script entrypoint

Run the study script directly:

```text
julia --project=. benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true
```

The profile is `quick`, the short validation path, or `full`. With
`--enforce=true` the study still writes its summary and error tables, then
exits with an error if any row misses its threshold; leave it off when you
only want the tables. The script also accepts `--scenarios=a,b`,
`--out-summary=...`, `--out-errors=...` and `--plots=true`. By default it
writes both tables to `output/` and makes no plots.

## CLI entrypoint

The packaged CLI runs the same study:

```text
julia --project=. src/interfaces/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

It writes `telemetry_orbit_accuracy_summary.csv` and
`telemetry_orbit_accuracy_errors.csv` to the output directory
(`output/telemetry` by default). Plots are off by default, so the command
works even where plotting packages or scripts are unavailable; add
`--plots=1` to make them. `smoke` is accepted as a profile name and runs
`quick`.

The CLI starts the study script in the examples environment, so prepare that
first with `julia --startup-file=no examples/setup.jl`. The prerequisites above
apply unchanged: without the `data/GRAMSuite.jl` submodule the command stops
with "Package GRAMSuite not found in current path".

## What this workflow is good for

Use the verification study when you need a known, repository-owned validation
path, thresholds that fail the run when they are missed, and a fixed output
directory to inspect afterwards.

## Interpreting finite-window orbital decay

A slope fitted to the osculating semi-major axis over a finite window picks up
orbital-period variations and their modulation as well as drag. J2, for
example, can produce a nonzero fitted slope in a drag-free run. A trend in
this two-body diagnostic alone therefore does not show energy loss through
drag.

Before reading such a slope as drag decay, compare it with a drag-free
reference propagation of the same arc. Set `drag_enabled = false` in the
reference scenario's manifest table and keep its initial state, gravity field
and other force settings. Evaluate both series over the same time window and
comparison grid, including the same sampling gaps.
`SpaceAGORA.TelemetryVerification.zero_referenced_decay(t_s, sma_m, t_ref_s,
sma_ref_m; period_s)` applies the same harmonic estimator to both series and
returns their slope difference as `decay_m_per_day`, together with
`raw_m_per_day` and `reference_m_per_day`. The subtraction removes the trend
the estimator finds in the reference, but it does not guarantee that model
mismatch or other perturbations cancel; check how the result changes with the
fit window and the estimator settings before attributing the difference to
drag.

The same module provides `visviva_sma(r_m, v_mps, mu)`, which computes the
osculating semi-major axis from position and velocity magnitudes, and
`flight_density_table`, which converts windowed, reference-subtracted decay
rates into a `time_s` / `rho_kgm3` table for the `tabulated_time` density
source. The conversion uses a near-circular drag relation with the supplied
mass and effective drag area (`cd_area_m2`), so the density it returns is only
as good as those inputs.

## Scheduled state anchors

A scenario can re-anchor the simulated state to an external trajectory at
scheduled times. Each anchor replaces the satellite's inertial position and
velocity with the supplied planet-centered J2000 state (m and m/s, the same
convention as `initial_state_j2000_m`). Between anchors the engine propagates
on its own, and at each anchor it prints a `state_anchor` line with the
position and velocity drift it removed: the open-loop error of the segment
just ended. An anchor can also reset the satellite's orbit counter
(`orbit_count`). The Odyssey scenario derives that value from each burn's
mission orbit number and the maneuver offset, so the orbit-keyed burn replay
stays aligned with the anchored trajectory even when the propagated one had
drifted across an apoapsis crossing.

The Odyssey scenario anchors to the NAV reconstruction after each of the
mission's trim burns. The anchors are listed under `[scenarios.state_anchors]`
in `test/telemetry_benchmark_manifest.toml` and generated by
`scripts/telemetry/extract_odyssey_burn_anchors.jl`. Replayed as one
open-loop arc, the 340-orbit aerobraking campaign is chaotic: through the loop
from periapsis timing to periapsis longitude to density, a last-bit change in
the build or the runner's CPU grows into a large apoapsis error, and the
pass/fail outcome becomes a property of the bit pattern rather than of the
model. With anchors at the burns, the graded arc is a sequence of open-loop
segments of 4 to 33 orbits, short enough that the result depends on the
model, while the flight comparison, the measured-density atmosphere and the
thresholds stay the same. The summary reports the number of anchors in
`state_anchor_count`.

The nightly run integrates the Odyssey scenario with plain Tsit5
(`SPACEAGORA_TELEMETRY_SOLVER_MODE` in `.github/workflows/nightly-stress.yml`)
instead of the automatic stiff sequence, whose Rosenbrock fallback engages in
the final orbits after the walk-out raise and slows the run without changing
the graded metrics. The study keeps no per-step solver storage: every
trajectory sample comes from the results CSV, and the returned solution is
read only for its return code and solver trace. The automatic stiff sequence
is the exception, because the engine detects a solver switch from the
per-step algorithm record. A `SPACEAGORA_SOLVER_SAVE_EVERYSTEP` or
`SPACEAGORA_SOLVER_SAVE_ON` value already set in the environment is passed
through unchanged.

Any run can use the same mechanism. `StateAnchor(elapsed_s, spacecraft_index,
state; orbit_count=nothing)` describes one anchor, and
`get_state_anchor_callback(anchors)` builds a discrete callback that you pass
through `run_simulation(args; extra_callbacks=(cb,))`. It is supported on the
first-order solver paths; under
`SPACEAGORA_SOLVER_MODE=gravity_backbone_split` the callback throws when it
initializes.
