# EDG package boundary

`SpaceAGORAEDG.jl` is the optional, experimental implementation of the typed
energy-depletion guidance used by SpaceAGORA. The core facade is
`SpaceAGORA.EDGServices`, contract `v"1.0.0"`. This separation follows the
[optional adapter standard](adapter_contract.md). It changes where numerical
kernels live; it does not establish manuscript equivalence or make the
controller flight-qualified. The [user guide](../user/edg_guidance.md) explains
the scientific policy and its current limits.

## Ownership

| Responsibility | Owner |
|---|---|
| Public configuration, state, guidance/control model and angle-command types | SpaceAGORA; existing type identities and constructors are retained. |
| Heat-rate and heat-load decisions, structural restriction, targeting, candidate planning, passage heat, guidance and control decisions | Eight source files in SpaceAGORAEDG. |
| Observation copies, atmosphere, epoch, frames, gravity and vehicle queries | SpaceAGORA, behind the versioned facade. |
| Integrator, heat integrals, event roots, solver masks, callback order, output and panel actuation | SpaceAGORA. |
| Fresh-run initialization, per-passage bookkeeping and shared guidance/control state | Core lifecycle hooks, operating on the run's own EDG state. |
| Older E-EDG/T-EDG bridge routes | Retained in core; shared typed numerical calls forward to the optional provider. |

The private core `typed_edg/algorithms.jl` is a compatibility dispatcher. Its
explicit `COMPATIBILITY_FUNCTIONS` inventory forwards existing helper names to
the companion; it contains no second implementation of the eight kernels.
Core source gates reject a numerical implementation added back to that owner.

## Supported calls and values

`IMPLEMENTED_FUNCTIONS` names `guidance_decision!` and `control_decision!`.
Core invokes them with the existing immutable configuration, the run-owned
mutable EDG state, a copied observation and a borrowed prediction context.
`CONSUMED_SERVICES` lists the exact types and query operations used by the
companion. The retained helper inventory supports compatibility routes and
regression fixtures; it is not a general-purpose simulation API.

`control_decision!` returns the core-owned `EDGControlDecision`: an angle command
in radians, an `apply` flag, a switch/cache action and diagnostics. The flag
distinguishes an invalid-index no-op from a valid zero-angle command. A cache
hit or a returned command does not certify scientific feasibility. Core applies
the command to the configured non-root panel links.

Observations contain inertial position in metres, inertial velocity in m/s,
mass in kg and copied cumulative per-link heat in J/cm². Spacecraft indices and
controlled-link indices are 1-based. Position and velocity are immutable
`SVector` values; heat and link collections are tuples. Geometry queries return
only the masses, dimensions, reference areas, reflection coefficients and angles
consumed by the kernels. Mutable physical geometry is never passed to the provider.

| Service | Meaning |
|---|---|
| `planet_constants`, `thermal_parameters`, `entry_altitude_m`, `spacecraft_geometry` | Immutable value projections, with entry altitude in metres. |
| `_edg_environment_state` | Current geodetic altitude, atmosphere and air-relative motion for a copied spacecraft observation. |
| `_edg_targeting_prediction_environment` | The same physical conventions at an explicit inertial state and elapsed simulation time. |
| `_edg_sample_prediction_atmosphere` | Altitude-only prediction at zero latitude/longitude; it is not interchangeable with a full trajectory query. |
| `_edg_ephemeris_time`, `_edg_planet_frame_lpi`, `frame_state` | Configured epoch and frame services; `l_pi` maps inertial components to planet-fixed components. |
| `gravity_acceleration` | The configured central/J2 and N-body acceleration, preserving the existing predictor's arithmetic. |
| `density_scale_height`, `_edg_targeting_control_clock` | Existing density-model fallback and callback-grid semantics. |

Environment values use metres, kg/m³, kelvin, m/s and pascals; molecular speed
ratio is dimensionless. Heat rate is W/cm² and specific energy is J/kg. Frame
queries use the configured epoch plus elapsed time. Atmosphere queries use
elapsed time. Winds are east/north/up; airspeed uses spacecraft velocity
relative to the rotating planet minus atmospheric velocity, including the
rotation transport term.

## Prediction capability and state lifetime

`AbstractPredictionContext` is an opaque, synchronous, call-scoped capability.
Only the facade query functions are supported. Ordinary property access is
refused, and the query closure checks expiry before accessing the core. Core
revokes each lease in a `finally` block, including when a kernel throws. Retaining
a handle does not extend its lifetime. Julia reflection into private fields or
modules is outside this interface; the facade is an ownership boundary, not a
security sandbox.

The capability does not expose `ODEParams`, `shared_buffers`, the live state,
mutable geometry, an atmosphere reseeding operation or solver controls. Existing
atmosphere caches and RNG behavior remain core-owned and use the configured
query path. Packaging introduces no new prediction instance, seed or cache
policy. Any future change to those semantics needs its own scientific review.

EDG may update its own run-owned plan, counters and diagnostic state. Under the
default `run_simulation` isolation, core deep-copies the configuration while
preserving the shared state between paired guidance and control. Separate runs
receive separate state. `isolate_state=false` means caller-owned state with fresh
initialization; it does not resume a prior passage. Checkpoint writing and resume
remain refused because they do not persist the complete EDG state.

Core owns normal callback periods, lifecycle scheduling and physical integration.
The geodetic EDG entry/exit observer owns passage heat baselines; the solver's
spherical atmosphere phase boundary remains a separate concept. Neither the
package move nor a diagnostic snapshot changes the physical cumulative heat
integral, callback ordering, root settings, limits or 3 Hz scientific cases.

## Installation and errors

The [installation guide](../user/installation_environment.md) describes the
optional `edg` target and the exact pair recorded in package metadata. Load
`SpaceAGORAEDG` before running EDG. Non-EDG simulations can still load and run
with core alone. Selecting an EDG model without its provider raises
`EDGUnavailableError` during preflight, before output or stepping; it never
silently selects another guidance law.

The companion checks that the facade exists and requires exactly contract
`v"1.0.0"` before including its algorithms. The extension repeats the version
check during initialization and dispatch. An incompatible contract raises
`CompatibilityError`; install a supported pair and start a fresh process.
Core retains existing configuration, state pairing and panel-mapping validation.
Prediction failures, unreachable targets and observed constraint violations keep
their existing diagnostics and scientific meaning.

## Evidence and limits

The extracted baseline is the reviewed PR289/current-main integration recorded
in the companion's `provenance/extraction.json`. That record preserves original
file hashes and identifies the service substitutions. Numerical ownership tests,
the existing equation, targeting, heat, event and lifecycle regressions,
core-only failure tests, both package load orders and a short analytic-atmosphere
run exercise the software boundary. Native reference campaigns and their exact
settings remain separate evidence; an installation test does not qualify them.

The initial package is experimental. Exact journal reproduction, continuous-family
optimality, flight use, new checkpoint support and new concurrency guarantees
are not established by this extraction. Existing observed violations remain
violations; no threshold or scientific setting is relaxed to pass packaging tests.
