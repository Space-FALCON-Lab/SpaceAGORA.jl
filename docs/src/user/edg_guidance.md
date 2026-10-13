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

EDG is currently part of SpaceAGORA itself. There is no separate EDG package
to install. Its proposed future adapter boundary is documented in the
[EDG ownership and lifecycle plan](../maintainer/edg_adapter_boundary.md);
that plan does not change the configuration described here.

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

In `:targeting` mode, guidance first predicts a range of exit energies and
apoapsides using low-drag and maximum-depletion trajectories. A target within
that range can proceed to a switch-time calculation: the controller applies
the admissible high-drag command before the switch and the minimum angle
after it. The switch times stored in the state are elapsed simulation times
in seconds, rather than durations measured from the latest callback.

Targeting also checks a sampled sequence of candidate switches for finite
predictions, heat-load feasibility, and the expected energy ordering. It uses
a contiguous acceptable portion of that sequence; at least two acceptable
samples are required. A failed check can therefore mean that the sampled
prediction did not establish an admissible bracket. It does not prove that
no physically feasible trajectory exists. If targeting cannot proceed, the
model falls back to maximum depletion when that mode is also enabled, or to
`:safe_low_drag` when targeting is the only configured mode. The latter name
means that the minimum panel angle is commanded; it does not certify that
all physical limits can be met.

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
in metres, rather than an altitude. The current target mapping uses the
predicted periapsis variation across the reachable bracket; a fixed
periapsis should not be assumed when translating a desired apoapsis into
energy. With additional perturbations enabled, this two-body energy remains
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

**In the current implementation, the heat-load input is cumulative.** The
controller reads the largest accumulated heat load among its controlled
links. If heat-load control is enabled and that value reaches
`heat_load_limit_j_cm2`, the command is reduced to the minimum angle. Passage
resets do not subtract an entry heat-load baseline or clear the propagated
heat state. A finite limit in a multiple-passage run must therefore be
interpreted with this cumulative behaviour in mind; it is not automatically
a renewed budget for each passage.

The structural setting also needs a precise interpretation. The controller
forms a drag area from the sum of each link's reference area multiplied by
its internal free-molecular drag coefficient. Let that sum at a candidate angle be
``A_D(\alpha)``. The pressure-valued setting defines the modeled force bound

```math
D(\alpha) = q_{\mathrm{dyn}} A_D(\alpha)
\;\leq\;
D_{\mathrm{limit}} = p_{\mathrm{limit}} A_D(\alpha_{\max}).
```

Controlled panels use the candidate angle in this calculation; other links
retain their modeled angles within the configured bounds. Thus
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
    guidance_rates = [3.0],
)
control_model = ControlModel(
    control_effectors = (edg_control,),
    control_rates = [0.1],
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
**periods in seconds**. The example runs guidance every 3 seconds and
control every 0.1 seconds. Choose these periods together with the integration
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
| `heat_load_limit_j_cm2` | Cumulative heat-load threshold in J/cm² under the current implementation. |
| `structural_load_limit_pa` | Pressure-valued parameter defining the drag-based force bound described above. |
| `target_apoapsis_radius_m` | Desired apoapsis radius in metres. Supply a finite positive radius when enabling targeting. |
| `heat_load_switch_solver` | `:closed_form` or `:tpbvp_integration` for the heat-load switching calculation. Selecting a solver does not establish physical validity for a new atmosphere or vehicle. |
| `planning_horizon_s` | Horizon used by the prediction machinery, not the mission duration. |
| `targeting_certification_samples` | Number of candidate switch samples, at least two; the default is nine. |
| `targeting_energy_order_tolerance_jkg` | Tolerance for the candidate energy-ordering check, in J/kg. |
| `targeting_heat_load_tolerance_j_cm2` | Tolerance for the candidate heat-load check, in J/cm². |

The physical limits default to `Inf`. A selected submode with an infinite
limit does not impose a finite protection threshold. Heat-rate and
structural restrictions require positive finite limits; zero does not
request zero heating or zero loading. Configure limits from
the vehicle and mission requirements, and validate them against the
propagated quantities. The stored `switch_recompute_interval_s` field is
currently inactive: changing it does not introduce periodic replanning.

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

Atmospheric entry and exit invalidate the affected spacecraft's cached
bracket, mode flags, and switch times. Repeated notification of the same
crossing does not erase a newly established plan, and ordinary callbacks
within a passage retain their cache. Bracketing counts accumulate across
passages within a run. Last-command telemetry remains until the next update;
neither physical heat loads nor panel geometry are cleared by the cache
reset. The next control callback can change the between-passage command
according to the configured mode's default.

Passage classification deserves care near the entry interface. The engine's
crossing events use a spherical radius, while the EDG in-passage query uses
geodetic altitude. At an exact-boundary start, the engine classifies outward
motion as outside and inward or tangential motion as inside, within its
floating-point boundary band. A known tangential-start case that subsequently
rises can retain its inside flag and cached plan through the coast. The
[Odyssey lifecycle discussion](../tutorials/odyssey_surrogate.md#EDG-state-between-runs-and-passes)
records this limit and the associated solver-phase behaviour.

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

Experimental changes to passage heat accounting, targeting, terminal events,
and constraint treatment are being evaluated separately. The ongoing GRAM
comparison and the experimental EDG work in PR #251 are not presented here
as scientifically accepted results. This guide reports the current software
contract and published methodological background. A publication or mission
study should additionally retain the exact source revision, atmosphere and
vehicle data, callback periods, solver settings, constraint definitions,
and predicted-versus-propagated diagnostics used for its conclusions.

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
