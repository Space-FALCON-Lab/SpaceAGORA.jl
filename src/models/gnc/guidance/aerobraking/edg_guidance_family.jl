# The EDG guidance family of the GNC algorithm contract (gnc_algorithms.jl): the
# energy-depletion guidance decision of the optional SpaceAGORAEDG provider, run
# through the request driver. The provider decides on a copy of the run's EDG
# state; trusted code checks the decision and installs it in that state, which the
# EDG control model reads. `AlgorithmGuidance` runs it like any other family.
module EDGGuidanceFamilies

using ..SimulationModel
import ..GNCAlgorithms as GA
import ..SimulationLifecycle as L
import ..GNCRuntime as GR
using ..AlgorithmLifecycle: AlgorithmGuidance
using ..SimulationModel.GuidanceHooks: AerobrakingEnergyDepletionConfig, AerobrakingEnergyDepletionState
const S = SimulationModel
const GH = SimulationModel.GuidanceHooks
const ES = SimulationModel.EDGServices

export EDGGuidanceFamily, EDGGuidanceCapabilities, EDGGuidanceAlgorithm, EDGGuidanceSettings,
    EDGGuidanceRequest, EDGGuidanceDecision, validate_edg_guidance_decision

"""
    EDGGuidanceFamily()

The EDG guidance family of the GNC algorithm contract: the energy-depletion
guidance decision for one spacecraft's drag passages. Its request is an
`EDGGuidanceRequest`, which copies the observation, the EDG configuration and the
run's EDG state; its product an `EDGGuidanceDecision`, proposed in an
`AlgorithmResult`; its trusted validation `validate_edg_guidance_decision`; its
install target the run's `AerobrakingEnergyDepletionState`, which the
`AerobrakingEnergyDepletionControlModel` reads; and its rejection policy `:fail`,
as when the provider raises: a refused decision stops the run. The family is
served by the `EDGServices` facade, and its contract version is the facade's
(`family_contract_version`).

An `AlgorithmGuidance` runs `EDGGuidanceAlgorithm()` with `EDGGuidanceSettings`.
The run has one spacecraft; a run with several uses
`AerobrakingEnergyDepletionGuidanceModel`.
"""
struct EDGGuidanceFamily <: GA.AbstractAlgorithmFamily end

# An EDG algorithm predicts with the run's gravity models (see AlgorithmGuidance's method).
GR._prediction_uses_run_gravity(::EDGGuidanceFamily) = true

"""
    EDGGuidanceCapabilities(; guidance_modes=())

What an algorithm of the EDG guidance family supports: the guidance modes
(`:max_energy_depletion`, `:targeting`) it can decide. The family's default
supports none.
"""
struct EDGGuidanceCapabilities
    guidance_modes::Tuple{Vararg{Symbol}}
    EDGGuidanceCapabilities(; guidance_modes=()) = new(Tuple(Symbol(m) for m in guidance_modes))
end

GA.default_capabilities(::EDGGuidanceFamily) = EDGGuidanceCapabilities()

"""
    EDGGuidanceAlgorithm()

The energy-depletion guidance decision of the optional SpaceAGORAEDG package, as
an algorithm of the EDG guidance family. It supports the maximum-depletion and
targeting modes, needs the package to be loaded (`require_algorithm_support`), and
draws no random numbers (`consumes_rng` is `false`). Each proposal calls the
provider's `guidance_decision!` on a copy of the request's EDG state and returns
the spacecraft's decision fields; a provider that changes anything else in that
state fails the proposal (`:failed`, termination `:changed_outside_decision`).
"""
struct EDGGuidanceAlgorithm <: GA.AbstractGNCAlgorithm end

GA.algorithm_family(::EDGGuidanceAlgorithm) = EDGGuidanceFamily()
GA.algorithm_capabilities(::EDGGuidanceAlgorithm) =
    EDGGuidanceCapabilities(guidance_modes=(:max_energy_depletion, :targeting))
GA.consumes_rng(::EDGGuidanceAlgorithm) = false
function GA.require_algorithm_support(::EDGGuidanceAlgorithm)
    ES.require_provider()
    return nothing
end

"""
    EDGGuidanceSettings(; config, state)

The settings of an `AlgorithmGuidance` that runs the EDG guidance family: the
run's `AerobrakingEnergyDepletionConfig` and `AerobrakingEnergyDepletionState`,
the same objects the `AerobrakingEnergyDepletionControlModel` holds. Decisions are
installed in `state`; the binding's `link` stays unset.
"""
struct EDGGuidanceSettings
    config::AerobrakingEnergyDepletionConfig
    state::AerobrakingEnergyDepletionState
    EDGGuidanceSettings(; config::AerobrakingEnergyDepletionConfig, state::AerobrakingEnergyDepletionState) =
        new(config, state)
end

"""
    EDGGuidanceRequest

The EDG guidance family's owned request for the spacecraft at
`spacecraft_index` (its position in the run, from 1) at elapsed time `time_s`
(s), with its `request_id` and `reason`:

- `config`: the run's `AerobrakingEnergyDepletionConfig`;
- `observation`: copied inertial position (m), velocity (m/s), mass (kg) and
  cumulative per-link heat (J/cm²) of every spacecraft, as the provider receives
  them;
- `state`: a copy of the run's `AerobrakingEnergyDepletionState`, including the
  control clock the engine reported when the run started;
- `services`: an opaque handle to the run's prediction services. Copies of the
  request share it; each provider call borrows its own call-scoped prediction
  context from it.

Trusted family code builds it (`family_request`); the algorithm receives its own
copy.
"""
struct EDGGuidanceRequest{O} <: GA.AbstractAlgorithmRequest
    request_id::Int
    reason::Symbol
    spacecraft_index::Int
    time_s::Float64
    config::AerobrakingEnergyDepletionConfig
    observation::O
    state::AerobrakingEnergyDepletionState
    services::ES._RunServices
end

# The state fields the guidance decision writes, each a per-spacecraft vector.
const DECISION_FIELDS = (:selected_mode, :targeting_active, :safe_low_drag, :energy_bracketing_evaluated,
    :energy_bracketing_count, :target_energy_jkg, :target_perturbation_energy_change_jkg,
    :bracket_min_energy_jkg, :bracket_max_energy_jkg, :targeting_switch_s, :last_switch_solve_t,
    :targeting_diagnostics)

"""
    EDGGuidanceDecision

The EDG guidance decision for the spacecraft at `spacecraft_index`: the values of
that spacecraft's guidance fields of the `AerobrakingEnergyDepletionState`. They
are `selected_mode`, `targeting_active`, `safe_low_drag`,
`energy_bracketing_evaluated` and `energy_bracketing_count`, the target and
bracket specific energies `target_energy_jkg`,
`target_perturbation_energy_change_jkg`, `bracket_min_energy_jkg` and
`bracket_max_energy_jkg` (J/kg), the targeting switch time `targeting_switch_s`
and the time of the last switch solve `last_switch_solve_t` (elapsed s), and the
copied `targeting_diagnostics`.
"""
struct EDGGuidanceDecision
    spacecraft_index::Int
    selected_mode::Symbol
    targeting_active::Bool
    safe_low_drag::Bool
    energy_bracketing_evaluated::Bool
    energy_bracketing_count::Int
    target_energy_jkg::Float64
    target_perturbation_energy_change_jkg::Float64
    bracket_min_energy_jkg::Float64
    bracket_max_energy_jkg::Float64
    targeting_switch_s::Float64
    last_switch_solve_t::Float64
    targeting_diagnostics::NamedTuple
end
EDGGuidanceDecision(state::AerobrakingEnergyDepletionState, i::Int) =
    EDGGuidanceDecision(i, (getfield(state, f)[i] for f in DECISION_FIELDS)...)

const _MODES = (:inactive, :max_energy_depletion, :targeting, :safe_low_drag)

"""
    validate_edg_guidance_decision(request::EDGGuidanceRequest, result) -> ValidationResult

The EDG guidance family's trusted validation of a proposal `result` against the
lifecycle's own copy of `request`. A result for another request is refused with
`:request_mismatch`, and one that is not a `:candidate` with its status. A
candidate is accepted, with reason `:accepted` and metrics `selected_mode` and
`targeting_switch_s`, when its product is an `EDGGuidanceDecision`
(`:invalid_product`) for the request's spacecraft (`:spacecraft_mismatch`) whose
mode is `:inactive`, `:safe_low_drag` or a guidance mode of the request's
configuration (`:unknown_mode`, `:mode_not_configured`), whose targeting plan is
consistent (`:inconsistent_targeting`: `targeting_active` exactly in the
`:targeting` mode), whose switch time is a time (`:invalid_switch_time`: never
`NaN`, and finite for an active targeting plan; `Inf` means no switch), and whose
bracketing count has not decreased (`:count_decreased`). Validation installs
nothing and draws from no random stream.
"""
function validate_edg_guidance_decision(request::EDGGuidanceRequest, result)
    refuse(reason, metrics=NamedTuple()) = GA.ValidationResult(false, reason, metrics)
    result.request_id == request.request_id ||
        return refuse(:request_mismatch, (result_request_id=result.request_id, request_id=request.request_id))
    result.status === :candidate || return refuse(result.status, (termination=result.termination,))
    decision = result.product
    decision isa EDGGuidanceDecision || return refuse(:invalid_product)
    i = request.spacecraft_index
    decision.spacecraft_index == i || return refuse(:spacecraft_mismatch)
    mode = decision.selected_mode
    mode in _MODES || return refuse(:unknown_mode, (selected_mode=mode,))
    mode in (:max_energy_depletion, :targeting) && !(mode in request.config.guidance_modes) &&
        return refuse(:mode_not_configured, (selected_mode=mode,))
    decision.targeting_active == (mode === :targeting) ||
        return refuse(:inconsistent_targeting, (selected_mode=mode, targeting_active=decision.targeting_active))
    switch_s = decision.targeting_switch_s
    (isnan(switch_s) || decision.targeting_active && !isfinite(switch_s)) &&
        return refuse(:invalid_switch_time, (targeting_switch_s=switch_s,))
    decision.energy_bracketing_count >= request.state.energy_bracketing_count[i] ||
        return refuse(:count_decreased)
    return GA.ValidationResult(true, :accepted, (selected_mode=mode, targeting_switch_s=decision.targeting_switch_s))
end

# ---- The family's methods of the general contract -----------------------------

function GA.validate_capabilities(::EDGGuidanceFamily, algorithm, request::EDGGuidanceRequest; checkpoint::Bool=false)
    checkpoint && return GA.ValidationResult(false, :restart_unsupported)
    supported = GA.algorithm_capabilities(algorithm).guidance_modes
    missing = Tuple(m for m in request.config.guidance_modes if !(m in supported))
    isempty(missing) || return GA.ValidationResult(false, :unsupported_guidance_mode, (guidance_modes=missing,))
    return GA.ValidationResult(true, :supported)
end

function GA.preflight_family(::EDGGuidanceFamily, algorithm, settings, args, spacecraft_index; isolate_state::Bool=true)
    settings isa EDGGuidanceSettings || throw(ArgumentError(
        "The EDG guidance family's settings are an EDGGuidanceSettings with the run's EDG configuration and state."))
    n = length(args.dynamics_model.spacecraft)
    n == 1 || throw(ArgumentError(
        "The EDG guidance family runs one spacecraft; for a run with $(n), use AerobrakingEnergyDepletionGuidanceModel."))
    length(settings.state.selected_mode) == n || throw(ArgumentError(
        "The EDG state of the guidance settings must have one entry per spacecraft ($(n))."))
    return nothing
end

function GA.family_link(::EDGGuidanceFamily, link, settings, p, spacecraft_index)
    link === nothing || throw(ArgumentError(
        "The EDG guidance family installs into the EDG state of its settings: leave the binding's `link` unset."))
    return settings.state
end

function GA.family_request(::EDGGuidanceFamily, settings::EDGGuidanceSettings, u, p, t, spacecraft_index, reason, request_id)
    return EDGGuidanceRequest(Int(request_id), reason, Int(spacecraft_index), Float64(t), settings.config,
        ES._observation_values(u), deepcopy(settings.state), ES._RunServices(p))
end

GA.validate_proposal(::EDGGuidanceFamily, request::EDGGuidanceRequest, result, trusted) =
    validate_edg_guidance_decision(request, result)

# The decision replaces the spacecraft's guidance fields, one assignment each.
function GA.install!(state::AerobrakingEnergyDepletionState, ::EDGGuidanceFamily, decision::EDGGuidanceDecision,
        request::EDGGuidanceRequest, t; validation, diagnostics)
    i = decision.spacecraft_index
    for f in DECISION_FIELDS
        getfield(state, f)[i] = getfield(decision, f)
    end
    return nothing
end

# A decision applies from its request on, until the next request replaces it.
GA.is_current(::EDGGuidanceFamily, request::EDGGuidanceRequest, decision::EDGGuidanceDecision, time_s) =
    Float64(time_s) >= request.time_s
GA.on_rejection(::EDGGuidanceFamily) = :fail
GA.family_contract_version(::EDGGuidanceFamily) = ES.CONTRACT_VERSION
GA.algorithm_result_type(::EDGGuidanceFamily) = GA.AlgorithmResult

# Whether `after` differs from `before` anywhere but the decision fields of row `i`.
function _changed_outside(before::AerobrakingEnergyDepletionState, after::AerobrakingEnergyDepletionState, i::Int)
    for f in fieldnames(AerobrakingEnergyDepletionState)
        a, b = getfield(before, f), getfield(after, f)
        if a isa AbstractVector && f in DECISION_FIELDS
            all(k -> k == i || isequal(a[k], b[k]), eachindex(a)) || return true
        else
            isequal(a, b) || return true
        end
    end
    return false
end

function GA.propose!(state, ::EDGGuidanceAlgorithm, request::EDGGuidanceRequest, rng)
    i = request.spacecraft_index
    # The provider decides on its own copy, so the request it was given is unchanged.
    edg = deepcopy(request.state)
    S.EDGAlgorithms.guidance_decision!(edg, request.config, request.observation, request.services, request.time_s, i)
    _changed_outside(request.state, edg, i) && return GA.AlgorithmResult(request_id=request.request_id,
        status=:failed, termination=:changed_outside_decision)
    return GA.AlgorithmResult(request_id=request.request_id, status=:candidate, termination=:decided,
        product=EDGGuidanceDecision(edg, i))
end

# ---- The binding in an EDG run ---------------------------------------------

# A binding with EDG guidance settings is the run's EDG guidance, as the EDG
# guidance model is, whichever algorithm of the family it runs: it produces the EDG
# state (a link the engine initializes once per run), the EDG checks pair it with
# the EDG control model, the passage heat observer keeps its state's baselines, and
# the run has the atmosphere events EDG needs.
const _EDGBinding = AlgorithmGuidance{<:GA.AbstractGNCAlgorithm, EDGGuidanceSettings}
GR.gnc_links(g::_EDGBinding) = (GR.Produces(g.settings.state),)
GH._edg_guidance_pairing(g::_EDGBinding) = (state=g.settings.state, config=g.settings.config)
S.ControlHooks._edg_heat_state(g::_EDGBinding) = g.settings.state
L.requires_atmosphere_events(::_EDGBinding) = true

end # module EDGGuidanceFamilies
