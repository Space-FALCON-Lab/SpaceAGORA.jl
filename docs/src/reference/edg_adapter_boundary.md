# EDG package boundary

`SpaceAGORAEDG.jl` is the optional, experimental package that computes
SpaceAGORA's typed energy-depletion guidance (EDG). SpaceAGORA reaches it through
the service module `SpaceAGORA.EDGServices`, contract `v"3.0.0"`: a module inside
SpaceAGORA that lists the two decisions the companion implements and the types
and queries SpaceAGORA supplies to it. The split follows the
[optional adapter standard](adapter_contract.md).

Moving the numerical kernels into their own package changed where the code
lives. Treat EDG as experimental: exact reproduction of the published results has
not been shown, and EDG is not approved for flight use. The
[user guide](../user/edg_guidance.md) explains the scientific policy and its
current limits.

## Ownership

| Responsibility | Owner |
|---|---|
| Public configuration, state, guidance/control model and angle-command types | SpaceAGORA; existing type identities and constructors are retained. |
| Heat-rate and heat-load decisions, structural restriction, targeting, candidate planning, passage heat, guidance and control decisions | SpaceAGORAEDG numerical implementation. |
| Observation copies, atmosphere, epoch, frames, gravity and vehicle queries | SpaceAGORA, through the versioned service module. |
| Integrator, heat integrals, event roots, solver masks, callback order, output and panel actuation | SpaceAGORA. |
| Fresh-run initialization, per-passage bookkeeping and shared guidance/control state | SpaceAGORA's lifecycle hooks, working on the run's own EDG state. |
| Older E-EDG/T-EDG strategy routes | Removed from core; SpaceAGORAEDG is the only EDG implementation. |

## Callers

Calls cross this boundary in two directions:

- **SpaceAGORA calls the companion.** SpaceAGORA's EDG guidance and control call
  the two decisions through `ext/SpaceAGORAEDGExt.jl`, which forwards each call
  by an allowed name. SpaceAGORA keeps the callbacks, the pairing of guidance and
  control, the panel commands and the per-run state.
- **The companion calls SpaceAGORA.** The companion's decisions import the names
  they use from `EDGServices`, one by one, in `src/SpaceAGORAEDG.jl`. During a
  decision, SpaceAGORA lends them a prediction context (see
  [below](#Prediction-capability-and-state-lifetime)).

The names are listed in `src/models/gnc/companions/edg/services.jl`. EDG forwards
calls by name; an external propagator, by contrast, adds methods to SpaceAGORA's
own functions.

SpaceAGORA reaches the companion only through the two decisions below. Some of the
companion's tests and research scripts still use private SpaceAGORA paths. These
paths are not supported;
maintainers find them listed, with the checks for a boundary change, in
`docs/contributing/companion_boundaries.md` in the
[contributor documentation](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/tree/main/docs/contributing).

## Supported calls and values

`IMPLEMENTED_FUNCTIONS` names the two decision services the companion
implements:

- `guidance_decision!(state, config, u, p, t, spacecraft_index)`;
- `control_decision!(state, config, controlled_panel_links, u, p, t, spacecraft_index)`.

Core invokes them with the run-owned mutable EDG state first, which the
companion updates, then the immutable configuration, the controlled panel links
(control only), a copied observation `u`, a borrowed prediction context `p`, the
elapsed time `t` in seconds and the spacecraft's position in the run,
`spacecraft_index`, from 1. `CONSUMED_SERVICES` lists the exact types and query
operations supplied to the companion by this contract. The retained helper list
serves older routes and regression fixtures; do not use it as a general
simulation API.

`control_decision!` returns the core-owned `EDGControlDecision`: an angle command
in radians, an `apply` flag, a switch/cache action and diagnostics. The flag
distinguishes an invalid-index no-op from a valid zero-angle command. A returned
command, or a cache hit, says only that the decision ran; judge its feasibility
from the run's results. SpaceAGORA applies the command to the configured non-root
panel links.

Observations contain inertial position in metres, inertial velocity in m/s,
mass in kg and copied cumulative per-link heat in J/cm². Spacecraft indices and
controlled-link indices are 1-based. Position and velocity are immutable
`SVector` values; heat and link collections are tuples. Geometry queries return
only the masses, dimensions, reference areas, reflection coefficients and angles
consumed by the kernels. Mutable physical geometry is never passed to the companion.

| Service | Meaning |
|---|---|
| `planet_constants`, `thermal_parameters`, `entry_altitude_m`, `spacecraft_geometry(p, spacecraft_index)` | Immutable value projections, with entry altitude in metres. |
| `_edg_environment_state(p, u, t, spacecraft_index)` | Current geodetic altitude, atmosphere and air-relative motion for a copied spacecraft observation. |
| `_edg_targeting_prediction_environment` | The same physical conventions at an explicit inertial state and elapsed simulation time. |
| `_edg_ephemeris_time`, `frame_state` | Configured epoch and frame services. |
| `gravity_acceleration` | The configured central/J2 and N-body acceleration, preserving the existing predictor's arithmetic. |
| `density_scale_height` | Existing density-model fallback. |
| `_edg_targeting_control_clock(p)` | The period and origin (s) of the EDG control ticks, on which targeting predicts its switch times. |

Every query takes the prediction context first. The control clock is the one the
engine reports to the EDG state (`gnc_clock_started!`) when the run starts,
before the models' initializers, and again at each solve's start: the context of a
companion call carries the clock of the EDG state passed to that call, and core no
longer searches the run's configuration for the EDG controller. A call without an
EDG state, or in a run without the EDG control model, has no clock, and the query
raises an `ArgumentError`.

Environment values use metres, kg/m³, kelvin, m/s and pascals; molecular speed
ratio is dimensionless. Heat rate is W/cm² and specific energy is J/kg. Frame
queries use the configured epoch plus elapsed time. Atmosphere queries use
elapsed time. Winds are east/north/up; airspeed uses spacecraft velocity
relative to the rotating planet minus atmospheric velocity, including the
rotation transport term.

## Prediction capability and state lifetime

During a decision, SpaceAGORA lends the companion a prediction context, an
`AbstractPredictionContext`. Use it only with the service module's query
functions: it is opaque, and ordinary property access is refused. The context is
a lease, valid only during that one synchronous call. SpaceAGORA revokes it in a
`finally` block when the call ends, also when a kernel throws, and each query
checks that the lease is still valid before it reads SpaceAGORA's data. Keeping a
handle does not extend its life. Julia reflection can still reach private fields
or modules; that is outside this interface, which marks ownership and is not a
security sandbox.

The context does not expose `ODEParams`, `shared_buffers`, the live state,
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
with SpaceAGORA alone. Selecting an EDG model without the companion raises
`EDGUnavailableError` during preflight, before output or stepping; it never
silently selects another guidance law.

The companion checks that the service module exists and requires exactly contract
`v"3.0.0"` before including its algorithms. The extension repeats the version
check during initialization and dispatch. An incompatible contract raises
`CompatibilityError`; install a supported pair and start a fresh process.
Core retains existing configuration, state pairing and panel-mapping validation.
Prediction failures, unreachable targets and observed constraint violations keep
their existing diagnostics and scientific meaning.

## The EDG guidance family

The guidance decision is also an algorithm of the
[GNC algorithm contract](rpo_planner_contract.md#The-general-algorithm-contract), so the generic binding
`AlgorithmGuidance` can run it in place of
`AerobrakingEnergyDepletionGuidanceModel`:

```julia
config = AerobrakingEnergyDepletionConfig(guidance_modes = (:max_energy_depletion,))
state = AerobrakingEnergyDepletionState(num_sats = 1)
guidance = AlgorithmGuidance(algorithm = EDGGuidanceAlgorithm(), spacecraft_index = 1,
    settings = EDGGuidanceSettings(config = config, state = state), replan_interval_s = 1 / 3)
control = AerobrakingEnergyDepletionControlModel(config, state)
```

| Part | The EDG guidance family (`EDGGuidanceFamily`) |
|---|---|
| Algorithm | `EDGGuidanceAlgorithm()`: the companion's `guidance_decision!`. Capabilities (`EDGGuidanceCapabilities`): the `:max_energy_depletion` and `:targeting` guidance modes. It draws no random numbers (`consumes_rng` is `false`) and needs the companion to be loaded. |
| Request | `EDGGuidanceRequest`: copies of the observation, the configuration and the EDG state (with its control clock), and an opaque handle from which each companion call borrows its prediction context. |
| Product | `EDGGuidanceDecision`: the spacecraft's guidance fields of the EDG state. The companion decides on a copy of the state; a proposal fails if it changes anything else. |
| Validation | `validate_edg_guidance_decision`: the request's identity and spacecraft, a known mode that the configuration allows (`:inactive` and `:safe_low_drag` always), and a bracketing count that has not decreased. |
| Install target | The EDG state of the settings, which the EDG control model reads. |
| On rejection | `:fail`, as when the companion raises. |
| Contract version | `family_contract_version` is this service module's contract version. |

Requests come at the binding's guidance ticks (`replan_interval_s`, with one at
the run's start), and the binding ticks before the same tick's control. The
binding runs one spacecraft; a run with several uses the EDG guidance model.
The binding keeps a record of every request, which costs memory over long
campaigns.

## Changes from contract 2.0.0

Contract 3.0.0 removes two queries the companion did not call:
`_edg_planet_frame_lpi` (the inertial-to-planet-fixed rotation) and
`_edg_sample_prediction_atmosphere` (an altitude-only atmosphere sample at zero
latitude and longitude). The companion's frame and atmosphere values come from
`_edg_targeting_prediction_environment`, `_edg_environment_state` and
`frame_state`. SpaceAGORA keeps both calculations for its own prediction code;
they are no longer part of the service module. Nothing else changes: the
decisions, the other queries, their arguments and their values are those of
contract 2.0.0.

SpaceAGORAEDG implements contract 3.0.0 from version 0.5.0. A companion built for
2.0.0 is refused at load, and SpaceAGORAEDG 0.5.0 refuses a SpaceAGORA that
provides 2.0.0.

## Changes from contract 1.0.0

Contract 2.0.0 follows SpaceAGORA's [extension conventions](extension_conventions.md):

- `guidance_decision!` and `control_decision!` take the EDG state first, then the
  configuration (1.0.0: configuration first).
- `_edg_environment_state(p, u, t, spacecraft_index)` takes the context first
  (1.0.0: `(u, p, t, i)`).
- `_edg_targeting_control_clock(p)` takes the context alone and returns the clock
  of the call's EDG state (1.0.0: `(config, p)`, which searched the configuration
  for the control model holding `config`). The EDG state gained
  `control_clock_period_s`.
- Spacecraft indices are named `spacecraft_index`.
- The other queries, `require_version`, `require_provider` and the errors are
  unchanged, and so are the numerical kernels. The values of the queries do not
  change.

SpaceAGORAEDG implements contract 2.0.0 from version 0.2.0. A companion built for
1.0.0 is refused at load.

## Evidence and limits

Treat this package as experimental. Its interface and installation checks show
that the package loads, answers calls and refuses a wrong contract version. To use
a native reference campaign or its scientific settings, run that campaign itself.

Moving the code into a package did not test exact reproduction of the journal
results, optimality over the continuous family, flight use, new checkpoint support
or new concurrency guarantees; do not claim them on the strength of the move.
Observed violations remain violations, and no threshold or scientific setting was
relaxed to make the packaging tests pass.
