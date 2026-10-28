# Candidate forecast verification

28 October 2026. Experimental one-spacecraft interface, based on core
`a6fd9fe666c9e9677c54b70d55b4230ae2c9f005`. This verifies restart and process
isolation for the supported configuration; it does not qualify an MPC campaign.

## Reproducible regression checks

Prepare the usual test environment with `julia --startup-file=no test/setup.jl`.
Then run:

```text
julia --startup-file=no --project=test test/unit/simulation/control_forecast_tests.jl
julia --startup-file=no --project=test test/unit/simulation/checkpoint_periodic_tick_tests.jl
julia --startup-file=no --project=test test/unit/simulation/checkpoint_activity_tests.jl
```

The forecast test uses a force command that accumulates at every tick. Losing or
duplicating the pending tick therefore changes the answer. It checks copied
source state, repeated trajectories, unchanged parent configuration, a nonzero
clock origin against independent piecewise integration, rejected runtime
references and shared storage, and cleanup after caller failure, missing package
and interruption. Worker-side assertions also check kernel load order with
duplicates, furnishing/cache reset and rejection of changed or missing files
before the existing pool is touched. Unsupported nested captures and workers
are rejected before creating another worker.

The integrated test passed locally on Julia 1.12.7 (26 parent assertions, plus
the assertions executed inside the worker). The ordinary checkpoint suites
passed all 334 assertions. In the analytic experiment, the forecast differed
from uninterrupted propagation by 1.1e-11 m and 1.9e-12 m/s. The nonzero-origin
case differed from independent integration by 3.7e-12 m and 9.1e-13 m/s. These
are numerical agreement checks, not physical accuracy estimates.

The public API surface, documentation contract, extensibility documentation and
docstring-coverage gates passed. Two complete documentation-build attempts
passed cross-reference/document checks but stopped before rendering when GitHub
returned HTTP 429 for existing repository links. A successful full docs build
is still required before merge; those link checks were not disabled.

## Additional local Earth experiment

An optional MPC fixture used the original campaign's degree/order-50 gravity,
SPICE, third bodies, solar pressure, polynomial atmosphere, free-molecular
aerodynamics and fixed area plan with its runtime limiter. It propagated for
120 s, captured a pre-command snapshot at 10 s, and compared a dedicated
worker's 110 s continuation with uninterrupted propagation.

The final differences were 1.68e-9 m and 3.35e-12 m/s. The source snapshot was
unchanged; a repeated forecast after adding another kernel to the worker's
pool gave exactly the same state and table. Six source kernel files were
hash-checked. The worker was removed afterward.

This experiment placed the preserved incoming Cartesian state at the scenario's
original start epoch solely to test engine equivalence. It is not another
campaign-performance measurement, does not cover the complete pass/coast to
apoapsis, and is not part of core CI. Its scripts, environments, source snapshot,
kernel inventory, logs and comparison are retained in the research qualification
packet as `forecast_full_physics_probe.jl`, `ForecastQualificationControls/`,
`forecast-full-physics-v3/` and `forecast-full-physics-v4-worker.log`. Earlier
driver failures were retained: Julia world-age access and an incorrect test
expectation for SPICE's missing-variable sentinel.

## Review and remaining scope

Independent review identified and closed concrete issues in kernel cache reset,
remote references, shared array storage, concurrent use of the process-wide
SPICE pool, and recursive worker creation. The final source review found no
remaining blocker in the stated narrow path. The dedicated-worker requirement
still applies: the forecast lock does not make unrelated simulations safe to
run concurrently in that process.

Full-pass event equivalence, campaign handover, cold checkpoint/resume across
that handover, broader state/planet cases, prediction cost and mission accuracy
requirements remain separate work. Default controllers and their optimizer
equations are unchanged. No physical-apoapsis search is installed automatically.
