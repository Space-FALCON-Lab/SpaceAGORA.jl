# Energy Depletion Guidance for aerobraking

Energy Depletion Guidance (EDG) controls the drag produced during an
aerobraking passage by changing the angle of selected spacecraft panels. A
larger exposed area generally removes orbital energy more quickly, but it also
changes aerodynamic loading and heating. EDG chooses panel commands within
that tradeoff: it can seek the greatest energy depletion allowed by its
configured constraints, or seek an orbital target that lies within its
predicted reachable range.

The original EDG formulation was developed for atmospheric passes with
articulating solar panels by [Falcone and Putnam (2022)](#References). The
generalized formulation of [Simha, Yu, and Falcone (2026)](#References) extends
the treatment to multiple planetary bodies, structural constraints, and
targeting. These papers establish the scientific context. This guide describes
the implementation currently shipped in SpaceAGORA, including the places
where a software setting, a prediction, and an accepted physical result must
be distinguished.

The numerical implementation lives in the optional, experimental
`SpaceAGORAEDG.jl` package. SpaceAGORA retains the public configuration/model
constructors, physical dynamics, lifecycle and actuation. Install the supported
pair with the [optional EDG setup](installation_environment.md), then load both:

```julia
using SpaceAGORA, SpaceAGORAEDG
```

The [EDG package boundary](../maintainer/edg_adapter_boundary.md) describes its
versioned queries and state ownership. Selecting an EDG model without loading
the provider gives a setup error before the simulation starts. The companion's
`examples/analytic_atmosphere.jl` is a short installed-package demonstration;
it is not a reproduction of either paper or a targeting-accuracy benchmark.

## What the guidance and controller do

The implementation separates two decisions. Guidance selects the mode and,
when targeting is enabled, evaluates a reachable energy bracket. Control
uses that decision to obtain a panel angle, applies the active heat-rate and
structural restrictions, and updates the selected panel geometry at each
control callback. Aerodynamic and thermal models then evaluate the resulting
geometry during propagation.

The two public model types are `AerobrakingEnergyDepletionGuidanceModel` and
`AerobrakingEnergyDepletionControlModel`. They share an immutable
`AerobrakingEnergyDepletionConfig` and the **same mutable**
`AerobrakingEnergyDepletionState`. Sharing the state is essential: a controller
with a separate, initially identical state does not receive the guidance
model's decisions. Simulation preflight checks this pairing and requires one
state entry per spacecraft.

In `:max_energy_depletion` mode, the controller begins from the configured
maximum angle. The selected submodes can reduce that command to satisfy the
internal heat-rate or structural model. Heat-load control can also schedule
a low-drag interval within the passage. The name describes the objective of
this policy; it is not a claim that every propagated case is a globally
optimal trajectory.

In `:targeting` mode, guidance evaluates executable one-switch candidates on
the controller's actual command grid. Each candidate propagates its own
atmospheric passage, evaluates the remaining heat budget and minimum-angle
admissibility, and coasts its exit state to apoapsis for the target-energy
mapping. The admissible family is a contiguous, ordered prefix of that grid.
This sampled certificate is not a continuous-family optimality proof.

Guidance owns the complete plan and its replanning clock. Before the accepted
switch it replans every two seconds by default, on a due guidance callback.
When targeting guidance and control are due together, guidance runs first;
control consumes the resulting plan without solving again or advancing that
clock. After the switch, the accepted plan remains cached for the passage.
Switch times are elapsed simulation seconds, not durations from the last update.

An upper-energy target beyond the sampled family selects minimum drag. A target
requiring greater depletion selects maximum depletion when that mode is enabled;
otherwise it selects `:safe_low_drag`. The latter name means the minimum angle,
not proof of constraint satisfaction. A failed candidate certificate remains
visible even if a fallback command is applied. A failed sampled bracket does not
prove that every physical trajectory is infeasible.

## Energy, heating, and loading

The energy used for orbital diagnostics is the planet-relative two-body
specific orbital energy,

```math
\varepsilon = \frac{\lVert\mathbf{v}_{\mathrm{i}}\rVert^2}{2}
              - \frac{\mu}{\lVert\mathbf{r}_{\mathrm{i}}\rVert},
```

where position and velocity are expressed in the inertial frame and the
result is in J/kg. For a bound Keplerian ellipse,
``\varepsilon=-\mu/(r_a+r_p)``. Consequently,
`target_apoapsis_radius_m` is an apoapsis **radius from the planet's centre**,
in metres, rather than an altitude. The target mapping uses each candidate's predicted periapsis and
its configured-gravity coast energy change. A fixed periapsis or a shared
correction for all candidates must not be assumed. With additional perturbations enabled, this two-body energy remains
a useful diagnostic, but its change should not automatically be attributed
entirely to atmospheric drag.

Atmospheric loading uses air-relative velocity in the rotating planet frame.
The environment query accounts for local wind when wind is enabled, so this
speed generally differs from the inertial speed in the energy expression. Dynamic pressure
is ``q_{\mathrm{dyn}}=\rho V_{\mathrm{air}}^2/2``, in Pa.

Heat rate and heat load are different quantities. EDG's heat-rate limit is
in W/cm², while its heat-load limit is in J/cm². Heat load is the time
integral of heat rate on a link; it is not the instantaneous heat rate at
the current state. The conversion from either square-centimetre unit to
its square-metre counterpart introduces a factor of ``10^4``.

The controller uses an internal Maxwellian heat-rate estimate to choose an
angle. That estimate takes the wall temperature to equal the sampled
atmospheric temperature. The propagated thermal model evaluates the actual
configured thermal settings and executed link geometry. Its saved heat
rates can therefore differ from `last_heat_rate_w_cm2`, which records the
controller's estimate. Assess a thermal constraint using the propagated
panel quantities as well as the guidance diagnostics.

The planning budget is passage heat: the maximum of the individual controlled
panel increments since the geodetic entry event. The controller copies each
panel's entry heat integral and subtracts that panel's own baseline. It never
subtracts the largest baseline from the largest current integral and never
resets physical heat. Starts inside the atmosphere use a zero baseline.
`last_heat_load_j_cm2` remains cumulative telemetry;
`last_pass_heat_load_j_cm2` records the passage quantity used by planning.

Heat-load control follows the reviewed two-switch plan, second-switch updates
and predictive heat-security rule. A cumulative threshold crossing does not
introduce a separate unreviewed hard-clipping law. Heat-budget status reports
exhaustion explicitly.

The structural setting also needs a precise interpretation. The controller
uses the journal plate-coefficient proxy: controlled-panel reference areas are
summed at the candidate angle, while all remaining reference area is grouped as
bus area at `π/2`. Both groups use the root link's reflection coefficient. Let
the resulting drag area at a candidate angle be
``A_D(\alpha)``. The pressure-valued setting defines the modeled force bound

```math
D(\alpha) = q_{\mathrm{dyn}} A_D(\alpha)
\;\leq\;
D_{\mathrm{limit}} = p_{\mathrm{limit}} A_D(\alpha_{\max}).
```

This structural proxy is separate from the trajectory predictor's per-link
core aerodynamic-coefficient evaluation. Thus
`structural_load_limit_pa` is a drag-based load proxy normalized by a
reference drag area. It is neither a finite-element stress limit nor a
direct upper bound on ambient dynamic pressure. Its connection to a
spacecraft's hardware limits must be established for that spacecraft.

If even the minimum angle exceeds an internal heat-rate or structural
bound, the corresponding control calculation returns the minimum angle.
The simulation is not terminated merely by that fallback. Commands are
also held between control callbacks, so checking a constraint only at those
callbacks is insufficient to establish its maximum over the passage.

## Configure a paired EDG model

The following constructs a paired model for a spacecraft with controllable
panels at link indices 2 and 3. It deliberately follows the simple panel-cap
exercise used in the Odyssey tutorial: the infinite heat-rate threshold
leaves the cap unmodified. It is a configuration example, not a thermal
protection policy.

```julia
using SpaceAGORA

edg_config = AerobrakingEnergyDepletionConfig(
    guidance_modes = (:max_energy_depletion,),
    max_energy_submodes = (:heat_rate,),
    controlled_panel_links = (2, 3),
    max_alpha_rad = deg2rad(30.0),
    min_alpha_rad = 1e-4,
    heat_rate_limit_w_cm2 = Inf,
)

edg_state = AerobrakingEnergyDepletionState(num_sats = 1)
edg_guidance = AerobrakingEnergyDepletionGuidanceModel(edg_config, edg_state)
edg_control = AerobrakingEnergyDepletionControlModel(edg_config, edg_state)

guidance_model = GuidanceModel(
    guidance_effectors = (edg_guidance,),
    guidance_rates = [1 / 3],
)
control_model = ControlModel(
    control_effectors = (edg_control,),
    control_rates = [1 / 3],
)
```

Pass these `guidance_model` and `control_model` values to the
`SimulationConfiguration` constructor along with the mission, environment,
dynamics, and initial-time settings for the case. Set `num_sats` to the
number of spacecraft in that configuration. The EDG state is indexed by
spacecraft position in the configuration, while `controlled_panel_links`
indexes the links within each spacecraft. Every selected link must exist
and must be a rotatable, non-root link. On a spacecraft with articulated
dynamics, this kinematic panel control is restricted to `:fixed` links
belonging to the root body. A panel on a moving joint body is refused;
EDG does not provide a joint motor for that body.

Despite their names, `guidance_rates` and `control_rates` contain callback
**periods in seconds**. The example runs both guidance and control at 3 Hz, every
1/3 second. This is the qualified EDG comparison cadence. Choose these periods together with the integration
and output resolution needed to resolve the passage. A smaller period is
not, by itself, evidence of a more accurate prediction or an accepted
constraint margin.

The default control effector uses the same panel selection as the EDG
configuration. Keep the prediction and actuation selections aligned unless
a different arrangement has been deliberately modeled and checked. For
fixed-attitude runs, the force and thermal models must also respond to the
changed panel geometry. The Odyssey example uses
`AerodynamicCoefficientfM(fixed_attitude_incidence=:attitude)` for this
purpose; see [panel angles and heating](simulation_configuration.md).
The EDG effector directly updates panel geometry. It does not add actuator
rate limits, hinge dynamics, or a spacecraft attitude controller.

The principal settings are:

| Setting | Meaning |
| --- | --- |
| `guidance_modes` | `(:max_energy_depletion,)`, `(:targeting,)`, or both. Including maximum depletion makes it available as the targeting fallback. |
| `max_energy_submodes` | A nonempty selection from `:heat_rate`, `:structural_load`, and `:heat_load`. These select the restrictions considered by the EDG calculations. |
| `min_alpha_rad`, `max_alpha_rad` | Command bounds in radians. Both must be finite, with `0 ≤ min_alpha_rad ≤ max_alpha_rad`. |
| `heat_rate_limit_w_cm2` | Internal heat-rate threshold in W/cm². Use a physically justified positive value when heat-rate protection is intended. |
| `heat_load_limit_j_cm2` | Passage heat-load threshold in J/cm², using per-panel entry baselines. |
| `structural_load_limit_pa` | Pressure-valued parameter defining the drag-based force bound described above. |
| `target_apoapsis_radius_m` | Desired apoapsis radius in metres. Supply a finite positive radius when enabling targeting. |
| `heat_load_switch_solver` | `:closed_form` or `:tpbvp_integration` for the heat-load switching calculation. Selecting a solver does not establish physical validity for a new atmosphere or vehicle. |
| `planning_horizon_s` | Horizon used by the prediction machinery, not the mission duration. |
| `targeting_certification_samples` | Retained compatibility utility setting, at least two; the default is nine. The supported journal planner samples the executable control grid instead. |
| `targeting_energy_order_tolerance_jkg` | Retained compatibility utility tolerance in J/kg; the supported journal family uses strict sampled ordering. |
| `targeting_heat_load_tolerance_j_cm2` | Retained compatibility utility tolerance; the journal certificate compares directly with the configured limit. |
| `switch_recompute_interval_s` | Guidance-owned targeting replan interval, default 2 seconds, before the applied switch. |
| `second_switch_reevaluation` | Enable the reviewed closed-form second-switch updates, default `true`. |
| `heat_load_security_mode` | Enable the reviewed predictive remaining-heat safeguard, default `true`. |

The physical limits default to `Inf`. A selected submode with an infinite
limit does not impose a finite protection threshold. Heat-rate and
structural restrictions require positive finite limits; zero does not
request zero heating or zero loading. Configure limits from
the vehicle and mission requirements, and validate them against the
propagated quantities. The replanning interval does not replace the 3 Hz command period. Commands
remain held between updates, and excursions between updates remain violations.

For a first complete run, use the
[Odyssey surrogate aerobraking tutorial](../tutorials/odyssey_surrogate.md).
From the repository root, its independent example environment is prepared
and run with:

```sh
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate_env/setup.jl
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl
```

The setup obtains the identified atmosphere and planetary assets required
by the example. It uses a frozen Mars atmosphere grid and does not require
a native GRAM installation. The default comparison changes the panel cap
from 90° to 30° and records the resulting trajectories and heating. It does
not exercise targeting, a finite thermal threshold, or a stochastic GRAM
campaign. Those distinctions make it a useful installation and geometry
exercise without treating it as validation of every EDG mode.

## Read the result and the guidance state

Use `run_simulation(args; return_results=true)` to retain the results table
and the configuration that actually ran. Under the default
`isolate_state=true`, SpaceAGORA deep-copies the input configuration, including
the shared guidance/control state. Consequently, the caller's original
`edg_state` is not the place to read the completed run's diagnostics.

For a configuration containing the paired guidance model above:

```julia
result = run_simulation(args; return_results = true)
run_guidance = only(filter(
    model -> model isa AerobrakingEnergyDepletionGuidanceModel,
    result.configuration.guidance_model.guidance_effectors,
))
state = run_guidance.state

state.selected_mode[1]
state.last_alpha_rad[1]
state.last_heat_load_j_cm2[1]
```

`result.table` contains the configured saved quantities; `result.files`
lists files written by the run. Passing `return_solution=true` as well
retains the ODE solution in `result.solution`. The EDG state itself holds
current or last-update diagnostics, not a complete time history:

| State fields | Interpretation |
| --- | --- |
| `selected_mode`, `targeting_active`, `safe_low_drag` | Current mode and fallback status. A passage reset can return the mode to `:inactive`. |
| `target_energy_jkg`, `bracket_min_energy_jkg`, `bracket_max_energy_jkg` | Target and predicted bracket energies in J/kg. Interpret them with the mode and bracketing flags. |
| `energy_bracketing_evaluated`, `energy_bracketing_count` | Whether the current passage has been bracketed, and the cumulative number of bracketing evaluations in the run. |
| `targeting_switch_s`, `heat_load_switches_s` | Cached switch times in elapsed simulation seconds. Infinite values are sentinels, not finite scheduled events. |
| `heat_load_switch_solved`, `last_switch_solve_t` | Bookkeeping for the heat-load solve and the latest switch calculation. A recorded calculation is not an independent feasibility certificate. |
| `last_alpha_rad`, `last_alpha_heat_rate_rad`, `last_alpha_structural_rad` | Applied command and the most recently computed angle bounds, in radians. |
| `last_heat_rate_w_cm2`, `last_heat_load_j_cm2`, `last_dynamic_pressure_pa` | Internal heat-rate estimate, controlled-link cumulative heat-load input, and sampled dynamic pressure. |

For a passage analysis, save the command, mode, switch diagnostics, panel
heat rates and heat loads alongside position, velocity, density, and time.
Use `SaveField` getters to retain custom diagnostic quantities; the
[outputs guide](outputs.md) describes the results and trajectory recorder,
and the Odyssey example demonstrates recording the executed panel geometry.
Saving that geometry matters: recomputing an old
state using the panels' final orientation does not reconstruct the forces
or heat rates that were applied at that time.

A whole-spacecraft maximum may include the bus, whereas the EDG heat-load
input covers its selected links. Label the panel set and aggregation used
in each plot or reported margin. Also distinguish a predicted exit energy
from the energy measured at the propagated exit event, and an internal
heat estimate from the thermal model's output. These comparisons are what
reveal prediction error and unmet constraints.

## State between passages and repeated runs

A fresh `run_simulation` initializes EDG state, including the counters.
Reusing a configuration therefore starts a fresh EDG run. Default state
isolation preserves the guidance/control pairing in the copied
configuration. With `isolate_state=false`, EDG still initializes its state,
but the caller retains responsibility for the configuration's other
mutable objects.

Geodetic EDG entry and exit invalidate the affected spacecraft's cached
bracket, mode flags, and switch times. Repeated notification of the same
crossing does not erase a newly established plan, and ordinary callbacks
within a passage retain their cache. Bracketing counts accumulate across
passages within a run. Last-command telemetry remains until the next update;
neither physical heat loads nor panel geometry are cleared by the cache
reset. The next control callback can change the between-passage command
according to the configured mode's default.

The engine's solver-phase events use a spherical radius; EDG uses its own
geodetic heat/prediction surface. EDG observes the accepted dense trajectory
before existing callback effects and before discrete guidance/control consumers.
It does not add a continuous root that could consume a coincident user stop.
The original callbacks retain their conditions, root-location options, order and
simultaneous-event masks. This preserves the solver's existing winner semantics;
it does not promise that every coincident user callback executes.

Each accepted interval is scanned at eleven nodes, and every detected geodetic
sign change is refined with dense interpolation. Per-panel heat is copied at the
refined crossing, before a later user callback can change physical heat. The bounded
scan does not guarantee detection of arbitrarily short excursions between nodes.
`heat_boundary_entry_time_s` and `heat_boundary_exit_time_s` record those refined
crossing times; `heat_boundary_observed_at_s` records when bookkeeping consumed
the crossing. These diagnostic clocks are separate from planning clocks.

An unregistered callback that stops just before the numerical geodetic sign
change keeps its original stop behavior; a proximity band is not treated as an
EDG crossing. For a scalar callback whose condition explicitly represents the
same geodetic EI surface, wrap the condition with
`SpaceAGORA.SimulationModel.SimulationCallbacks.EDGGeodeticBoundaryCondition(condition, spacecraft_index, altitude_m)`.
The callable delegates the original condition unchanged. Registration checks the
spacecraft index and exact configured EI altitude, and requires event geometry
within frame-transform roundoff. At such a declared event, heat is captured at
the original solver event time before its effect. Geodetic normal velocity
chooses entry or exit; a tangent does not force a transition. Registration is an
explicit surface-identity contract, not a test of an arbitrary function's
mathematical equivalence. It currently supports scalar continuous callbacks;
vector and generic callbacks retain raw geometric observation. The native-free
Odyssey example registers its known outbound EI event this way.

Only an actual geodetic transition resets the plan and records passage heat.
A spherical phase transition alone, a duplicate observation or an outside
computed trajectory touching the boundary tangentially must not invent a new passage.
A numerically integrated trajectory that actually dips below the surface is a
small passage in that computed trajectory, even if an analytic ideal is tangent.
Original callback effects that change the spacecraft's side are reconciled
before the next same-time consumer. Initialization changes define the initial
side without inventing an entry or erasing a prepared initial plan. At startup
only, the established transform-roundoff band uses geodetic normal velocity to
classify a supplied EI state; it is not used to advance later generic crossings.

Maximum-depletion control can run without a paired guidance model, but
targeting requires guidance to establish its bracket. A configuration may
contain at most one EDG guidance model and one EDG controller. EDG checkpoint
writing and resume are currently refused because the checkpoint format does
not preserve its passage history. Save analysis outputs normally, but do not
plan a restart workflow around an EDG checkpoint.

## Scope of the current implementation

Successful execution establishes that a configured calculation ran. It does
not, by itself, establish target accuracy, constraint satisfaction, or
agreement with a published study. The heat-rate estimate, trajectory
prediction, atmosphere representation, held control commands, and propagated
vehicle model all contribute to the gap between a planned switch and its
achieved result.

In particular, the current prediction paths do not all sample the atmosphere
in the same way. The targeting trajectory uses position-dependent queries
and wind, whereas an altitude-only prediction helper samples at zero
latitude and longitude and omits wind. A regional atmosphere grid that is
adequate for the propagated trajectory is not necessarily adequate for that
helper. Select an atmosphere whose domain and meaning cover the prediction
path being exercised, and compare the predicted passage with the propagated
one before interpreting its margins.

The supported release scope is bounded and experimental. The reviewed native
cases cover active targeting, depletion and a two-pass lifecycle. They do not
establish arbitrary-target accuracy, near-budget performance, real-time
throughput, universal constraint satisfaction or exact journal reproduction.
At 3 Hz the depletion and two-pass cases have between-update heat-rate
excursions. These remain reported violations. If even minimum angle exceeds
the structural limit, the original minimum-drag policy remains in effect with
explicit infeasibility; no new abort, limit relaxation or orbit change is added.

The original calibration inputs, fit settings and plotted journal run identity
remain unverified. Mars epoch placement, Earth fit normalization and Titan
interface/domain differences remain provenance limits. Current runs must retain
explicit initialization, exact source/data identities, callback periods, solver
settings and predicted-versus-propagated diagnostics. The numerical shooting
route is experimental compatibility/reference work; reduced sensitivity checks
do not establish a full Cartesian moving-event adjoint.

## References

1. Falcone, G., and Putnam, Z. R., “Energy Depletion Guidance for Aerobraking
   Atmospheric Passes,” *Journal of Guidance, Control, and Dynamics*, Vol. 45,
   No. 4, 2022, pp. 651–668.
   DOI: [10.2514/1.G006171](https://doi.org/10.2514/1.G006171).
2. Simha, N., Yu, E., and Falcone, G., “Generalized Energy Depletion Guidance
   for Aerobraking at Mars, Venus, Earth, and Titan,” *AIAA SCITECH 2026 Forum*,
   Orlando, FL, 2026, AIAA Paper 2026-1378.
   DOI: [10.2514/6.2026-1378](https://doi.org/10.2514/6.2026-1378).

The [publications page](../publications.md) lists the related SpaceAGORA
studies and provides BibTeX for the SciTech paper.
