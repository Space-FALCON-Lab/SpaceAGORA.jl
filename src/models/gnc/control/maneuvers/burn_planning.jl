# The burn-plan family of the GNC algorithm contract (gnc_algorithms.jl): an
# algorithm proposes one timed constant-thrust burn for a spacecraft, trusted
# code checks it with the burn scheduler's own size rule, and an accepted plan is
# installed in the scheduler's burn-plan buffer, which executes it on its
# [start, stop) window. `AlgorithmGuidance` runs any planner of the family.
module BurnPlanning

using LinearAlgebra: norm
using StaticArrays: SVector
using ..SimulationModel
import ..GNCAlgorithms as GA
using ..SimulationModel.CommandTypes: PropulsiveBurnPlan
using ..SimulationModel.ThrusterModels: BaseThrusterModel
const S = SimulationModel
const H = SimulationModel.ControlHooks

export AbstractBurnPlanner, BurnPlanFamily, BurnPlanRequest, BurnPlanSlot, ApoapsisBurnPlanner
export validate_burn_plan, burn_plan_is_current, burn_plan_sizing

"""
    BurnPlanFamily()

The burn-plan family of the GNC algorithm contract. Its request is a
`BurnPlanRequest`, a snapshot of one spacecraft and of its burn scheduler; its
product a `PropulsiveBurnPlan`, proposed in an `AlgorithmResult`; its trusted
validation `validate_burn_plan`, which applies the size check the burn
scheduler (`BaseThrusterModel`) applies to its own burns and checks the window;
its install target the spacecraft's slot in the scheduler's burn-plan buffer
(`BurnPlanSlot`); its currency rule `burn_plan_is_current`; and its rejection
policy `:keep_previous`: a refused proposal leaves the installed plan, and the
binding records the refusal with the validation's reason.

An `AlgorithmGuidance` runs any planner of the family (an
`AbstractBurnPlanner`, such as `ApoapsisBurnPlanner`). The configuration needs
one `BaseThrusterModel`, whose `Δv` is zero and whose window is unset for that
spacecraft, and no thruster guidance model commanding the scheduler's own
burns. The scheduler takes up an installed plan at the end of the step it was
installed in: it stops the solver at the plan's window ends and fires the
engine inside `[start, stop)`, with the plan's direction, thrust and specific
impulse. The family takes no `settings`.
"""
struct BurnPlanFamily <: GA.AbstractAlgorithmFamily end

"""
    AbstractBurnPlanner <: AbstractGNCAlgorithm

An algorithm of the burn-plan family (`BurnPlanFamily`, inherited). A planner
implements `propose!(state, planner, request::BurnPlanRequest, rng)` and
returns an `AlgorithmResult` whose product, for a `:candidate`, is the
`PropulsiveBurnPlan` it proposes. `ApoapsisBurnPlanner` is a subtype.
"""
abstract type AbstractBurnPlanner <: GA.AbstractGNCAlgorithm end
GA.algorithm_family(::AbstractBurnPlanner) = BurnPlanFamily()

_tuple3(x, name) = length(x) == 3 ? ntuple(k -> Float64(x[k]), 3) :
    throw(ArgumentError("$name must contain 3 values."))

"""
    BurnPlanRequest(; request_id, reason, spacecraft_index, time_s, position_m, velocity_mps,
                    mass_kg, available_propellant_kg=nothing, thrust_n, isp_s, mu_m3ps2,
                    planet_radius_m, entry_interface_altitude_m, installed_window_s=(NaN, NaN))

The burn-plan family's owned request snapshot, in SI units, for the spacecraft
at `spacecraft_index` (its position in the run, from 1) at elapsed time
`time_s`: its inertial position and velocity relative to the planet's centre,
its mass, the propellant it can still use (`nothing` when the vehicle declares
no limit), the thrust and specific impulse (s) its burn scheduler configures
for it, the planet's gravitational parameter and equatorial radius, the
altitude of the atmosphere's entry interface, and the window of the plan
installed when the request was made (`(NaN, NaN)` when none is). Trusted
family code builds it (`family_request`); the planner receives its own copy.
"""
struct BurnPlanRequest <: GA.AbstractAlgorithmRequest
    request_id::Int
    reason::Symbol
    spacecraft_index::Int
    time_s::Float64
    position_m::NTuple{3, Float64}
    velocity_mps::NTuple{3, Float64}
    mass_kg::Float64
    available_propellant_kg::Union{Nothing, Float64}
    thrust_n::Float64
    isp_s::Float64
    mu_m3ps2::Float64
    planet_radius_m::Float64
    entry_interface_altitude_m::Float64
    installed_window_s::NTuple{2, Float64}
    function BurnPlanRequest(; request_id::Integer, reason::Symbol, spacecraft_index::Integer, time_s::Real,
            position_m, velocity_mps, mass_kg::Real, available_propellant_kg=nothing,
            thrust_n::Real, isp_s::Real, mu_m3ps2::Real, planet_radius_m::Real,
            entry_interface_altitude_m::Real, installed_window_s=(NaN, NaN))
        request_id > 0 || throw(ArgumentError("request_id must be positive."))
        spacecraft_index > 0 || throw(ArgumentError("spacecraft_index must be positive."))
        isfinite(time_s) || throw(ArgumentError("time_s must be finite."))
        length(installed_window_s) == 2 || throw(ArgumentError("installed_window_s must contain 2 values."))
        new(Int(request_id), reason, Int(spacecraft_index), Float64(time_s),
            _tuple3(position_m, "position_m"), _tuple3(velocity_mps, "velocity_mps"), Float64(mass_kg),
            available_propellant_kg === nothing ? nothing : Float64(available_propellant_kg),
            Float64(thrust_n), Float64(isp_s), Float64(mu_m3ps2), Float64(planet_radius_m),
            Float64(entry_interface_altitude_m), (Float64(installed_window_s[1]), Float64(installed_window_s[2])))
    end
end

"""
    BurnPlanSlot(plans, spacecraft_index)

The burn-plan family's install target: the entry for one spacecraft of the burn
scheduler's plan buffer `plans`. `AlgorithmGuidance` finds it in the run
(`family_link`); the binding's `link` stays unset.
"""
struct BurnPlanSlot
    plans::Vector{PropulsiveBurnPlan}
    spacecraft_index::Int
end

"""
    burn_plan_sizing(mass_kg, delta_v_mps, direction_rad, thrust_n, isp_s; available_propellant_kg=nothing)
        -> NamedTuple

The burn scheduler's size check of a constant-thrust burn of `delta_v_mps` (m/s)
in direction `direction_rad` for a spacecraft of `mass_kg` (kg), at `thrust_n` (N)
and specific impulse `isp_s` (s): `(reason, propellant_required_kg,
commanded_impulse_n_s, burn_duration_s)`, in kg, N·s and s, with reason
`:accepted`. Otherwise the reason names the problem (`:invalid_mass`,
`:invalid_delta_v`, `:invalid_direction`, `:invalid_thrust`, `:invalid_isp`,
`:invalid_mass_fraction`, `:invalid_burn_size`) and the values are NaN, or it is
`:insufficient_propellant`, with the values, when the burn needs more than
`available_propellant_kg` (`nothing`: no limit). `validate_burn_plan` checks a
proposed `PropulsiveBurnPlan` against it, so a planner states these values and a
window of this duration.
"""
burn_plan_sizing(mass_kg::Real, delta_v_mps::Real, direction_rad::Real, thrust_n::Real, isp_s::Real;
        available_propellant_kg::Union{Nothing, Real}=nothing) =
    H._burn_plan_sizing(Float64(mass_kg), Float64(delta_v_mps), Float64(direction_rad), Float64(thrust_n),
        Float64(isp_s), available_propellant_kg === nothing ? nothing : Float64(available_propellant_kg))

# Relative agreement of a stated value with the trusted one.
_agrees(stated, trusted) = abs(stated - trusted) <= 1e-9 * abs(trusted)

"""
    validate_burn_plan(request::BurnPlanRequest, result) -> ValidationResult

The burn-plan family's trusted validation of a proposal `result` against the
lifecycle's own copy of `request`. A result for another request (a different
`request_id`, such as a stale or cached one) is refused with `:request_mismatch`,
whatever its status. A `:candidate` for this request is accepted, with reason
`:accepted` and metrics `propellant_required_kg`, `commanded_impulse_n_s` and
`burn_duration_s`, when its product is a valid `PropulsiveBurnPlan` that:

1. does not replace a burn in progress at the request's time
   (`:burn_in_progress`);
2. passes the burn scheduler's own size check at the request's mass, available
   propellant, thrust and specific impulse (`:invalid_mass`, `:invalid_delta_v`,
   `:invalid_direction`, `:invalid_thrust`, `:invalid_isp`,
   `:invalid_mass_fraction`, `:invalid_burn_size`, `:insufficient_propellant`);
3. uses the request's thrust and specific impulse (`:thrust_mismatch`,
   `:isp_mismatch`) and states the propellant and impulse of that check to a
   relative 1e-9 (`:propellant_mismatch`, `:impulse_mismatch`);
4. has a window with finite ends and `stop > start` (`:invalid_window`) that
   starts no earlier than the request (`:window_before_request`) and lasts the
   check's burn duration, to a relative 1e-9 plus the rounding of its ends
   (`:duration_mismatch`).

Any other result is refused with its status as the reason (`:infeasible`,
`:unsupported`, `:failed`), and a candidate whose product is not a valid
`PropulsiveBurnPlan` with `:invalid_product`. Validation installs nothing and
draws from no random stream.
"""
function validate_burn_plan(request::BurnPlanRequest, result)
    refuse(reason, metrics=NamedTuple()) = GA.ValidationResult(false, reason, metrics)
    result.request_id == request.request_id ||
        return refuse(:request_mismatch, (result_request_id=result.request_id, request_id=request.request_id))
    result.status === :candidate || return refuse(result.status, (termination=result.termination,))
    plan = result.product
    (plan isa PropulsiveBurnPlan && plan.valid) || return refuse(:invalid_product)
    installed_start, installed_stop = request.installed_window_s
    installed_start <= request.time_s < installed_stop &&
        return refuse(:burn_in_progress, (installed_window_s=request.installed_window_s,))
    sizing = H._burn_plan_sizing(request.mass_kg, plan.delta_v_mps, plan.direction_rad,
        request.thrust_n, request.isp_s, request.available_propellant_kg)
    sizing.reason === :accepted || return refuse(sizing.reason,
        (propellant_required_kg=sizing.propellant_required_kg, available_propellant_kg=request.available_propellant_kg))
    plan.thrust_n == request.thrust_n || return refuse(:thrust_mismatch)
    plan.isp_s == request.isp_s || return refuse(:isp_mismatch)
    _agrees(plan.propellant_required_kg, sizing.propellant_required_kg) || return refuse(:propellant_mismatch)
    _agrees(plan.commanded_impulse_n_s, sizing.commanded_impulse_n_s) || return refuse(:impulse_mismatch)
    start_s, stop_s = plan.start_burn_s, plan.stop_burn_s
    (isfinite(start_s) && isfinite(stop_s) && stop_s > start_s) || return refuse(:invalid_window)
    start_s >= request.time_s || return refuse(:window_before_request)
    duration_s = stop_s - start_s
    abs(duration_s - sizing.burn_duration_s) <= 1e-9 * sizing.burn_duration_s + 4 * eps(max(abs(start_s), abs(stop_s))) ||
        return refuse(:duration_mismatch, (window_duration_s=duration_s, burn_duration_s=sizing.burn_duration_s))
    return GA.ValidationResult(true, :accepted, (propellant_required_kg=sizing.propellant_required_kg,
        commanded_impulse_n_s=sizing.commanded_impulse_n_s, burn_duration_s=sizing.burn_duration_s))
end

"""
    burn_plan_is_current(request::BurnPlanRequest, plan::PropulsiveBurnPlan, time_s) -> Bool

Whether an installed `plan` may still be used at `time_s`: it is valid and its
window has not ended. The scheduler executes it until then; afterwards it
clears the plan at its next scheduling pass.
"""
burn_plan_is_current(request::BurnPlanRequest, plan::PropulsiveBurnPlan, time_s) =
    plan.valid && isfinite(plan.stop_burn_s) && Float64(time_s) < plan.stop_burn_s

# ---- The family's methods of the general contract -----------------------------

_burn_schedulers(args) = [c for c in args.control_model.control_effectors if c isa BaseThrusterModel]
_thruster_guidance(g) = g isa S.GuidanceModels.AerobrakingCampaignPropulsiveManeuverGuidanceModel ||
    g isa S.GuidanceModels.ApoapsisTargetPeriapsisRaiseGuidanceModel

function GA.preflight_family(::BurnPlanFamily, planner, settings, args, spacecraft_index; isolate_state::Bool=true)
    i = Int(spacecraft_index)
    settings === nothing || throw(ArgumentError(
        "The burn-plan family takes no settings: configure the planner instead."))
    schedulers = _burn_schedulers(args)
    length(schedulers) == 1 || throw(ArgumentError(
        "The burn-plan family installs into the burn scheduler's buffer: configure exactly one BaseThrusterModel " *
        "(found $(length(schedulers)))."))
    scheduler = only(schedulers)
    i <= length(scheduler.thrust) || throw(ArgumentError(
        "The BaseThrusterModel has no entry for spacecraft $(i)."))
    own_window = isfinite(scheduler.start_burn_time[i]) && isfinite(scheduler.stop_burn_time[i]) &&
        scheduler.stop_burn_time[i] > scheduler.start_burn_time[i]
    (scheduler.Δv[i] == 0.0 && !own_window) || throw(ArgumentError(
        "The BaseThrusterModel plans its own burn for spacecraft $(i): set its Δv to zero and leave its window " *
        "unset, so that only the AlgorithmGuidance of the burn-plan family plans that spacecraft's burns."))
    any(_thruster_guidance, args.guidance_model.guidance_effectors) && throw(ArgumentError(
        "A thruster guidance model commands the burn scheduler's own burns; it cannot run beside an " *
        "AlgorithmGuidance of the burn-plan family."))
    return nothing
end

function GA.family_link(::BurnPlanFamily, link, settings, p, spacecraft_index)
    link === nothing || throw(ArgumentError(
        "The burn-plan family installs into the burn scheduler's own buffer: leave the binding's `link` unset."))
    return BurnPlanSlot(p.shared_buffers.maneuver_burn_plans, Int(spacecraft_index))
end

function GA.family_request(::BurnPlanFamily, settings, u, p, t, spacecraft_index, reason, request_id)
    i = Int(spacecraft_index)
    x = u.sc[i]
    scheduler = only(_burn_schedulers(p.args))
    environment = p.args.environment_model
    mass_kg = Float64(x.mass)
    installed = p.shared_buffers.maneuver_burn_plans[i]
    return BurnPlanRequest(; request_id=request_id, reason=reason, spacecraft_index=i, time_s=t,
        position_m=x.pos, velocity_mps=x.vel, mass_kg=mass_kg,
        available_propellant_kg=H._available_propellant_kg(p, i, mass_kg),
        thrust_n=scheduler.thrust[i], isp_s=scheduler.Isp[i], mu_m3ps2=environment.planet.μ,
        planet_radius_m=environment.planet.Rp_e, entry_interface_altitude_m=environment.EI * 1000,
        installed_window_s=installed.valid ? (installed.start_burn_s, installed.stop_burn_s) : (NaN, NaN))
end

GA.validate_proposal(::BurnPlanFamily, request::BurnPlanRequest, result, trusted) = validate_burn_plan(request, result)

# The installed plan replaces the slot's entry in one assignment.
function GA.install!(slot::BurnPlanSlot, ::BurnPlanFamily, plan::PropulsiveBurnPlan, request::BurnPlanRequest, t;
        validation, diagnostics)
    slot.plans[slot.spacecraft_index] = plan
    return nothing
end

GA.is_current(::BurnPlanFamily, request::BurnPlanRequest, plan::PropulsiveBurnPlan, time_s) =
    burn_plan_is_current(request, plan, time_s)
GA.on_rejection(::BurnPlanFamily) = :keep_previous
GA.family_contract_version(::BurnPlanFamily) = v"1.0.0"
GA.algorithm_result_type(::BurnPlanFamily) = GA.AlgorithmResult
# A rejection record names the validation's reason as `rejection`.
GA.rejected_record(::BurnPlanFamily, request::BurnPlanRequest, result, validation, reason, t) =
    (event=:rejected, time_s=t, request_id=request.request_id, reason=reason, rejection=validation.reason,
     spacecraft_index=request.spacecraft_index,
     termination=result === nothing ? :not_proposed : result.termination, validation=validation)

# ---- A planner with the scheduler's own timing rule ---------------------------

"""
    ApoapsisBurnPlanner(; delta_v_mps, direction_rad=0.0)

A planner of the burn-plan family (`AbstractBurnPlanner`) with the burn
scheduler's own timing rule: a burn of `delta_v_mps` (m/s, finite and
positive), prograde for `direction_rad = 0` and retrograde for `π`, at the
request's thrust and specific impulse, centred on the next apoapsis. Like the
scheduler, it plans only at or above the atmosphere's entry interface and
before apoapsis (or on a near-circular orbit); otherwise its result is
`:infeasible`, with termination `:inside_atmosphere`, `:past_apoapsis` or
`:invalid_orbit`. It states the burn's propellant and impulse from the
scheduler's size check without the propellant limit, so the family's
validation is what refuses a burn the vehicle cannot afford. It draws no
random numbers (`consumes_rng` is `false`).
"""
struct ApoapsisBurnPlanner <: AbstractBurnPlanner
    delta_v_mps::Float64
    direction_rad::Float64
    function ApoapsisBurnPlanner(; delta_v_mps::Real, direction_rad::Real=0.0)
        (isfinite(delta_v_mps) && delta_v_mps > 0) || throw(ArgumentError("delta_v_mps must be finite and positive."))
        isfinite(direction_rad) || throw(ArgumentError("direction_rad must be finite."))
        new(Float64(delta_v_mps), Float64(direction_rad))
    end
end

GA.consumes_rng(::ApoapsisBurnPlanner) = false

function GA.propose!(state, planner::ApoapsisBurnPlanner, request::BurnPlanRequest, rng)
    no_burn(why) = GA.AlgorithmResult(request_id=request.request_id, status=:infeasible, termination=why)
    mass_kg = request.mass_kg
    (isfinite(mass_kg) && mass_kg > 0.0) || return no_burn(:invalid_mass)
    pos = SVector{3, Float64}(request.position_m)
    vel = SVector{3, Float64}(request.velocity_mps)
    a, e, _, _, _, ν = S.FrameTransforms.rvtoorbitalelement(pos, vel, (μ=request.mu_m3ps2,))
    (isfinite(a) && isfinite(e) && isfinite(ν) && a > 0.0 && 0.0 <= e < 1.0) || return no_burn(:invalid_orbit)
    # The scheduler's eligibility: out of the atmosphere and before apoapsis,
    # where a near-circular orbit's true anomaly can be arbitrary.
    altitude_m = norm(pos) - request.planet_radius_m
    ν_wrapped = S.FrameTransforms._wrap_2pi(Float64(ν))
    pre_apoapsis = e <= 1e-8 || ν_wrapped < π - 1e-12
    altitude_m >= request.entry_interface_altitude_m - 1e-6 || return no_burn(:inside_atmosphere)
    pre_apoapsis || return no_burn(:past_apoapsis)
    sizing = H._burn_plan_sizing(mass_kg, planner.delta_v_mps, planner.direction_rad,
        request.thrust_n, request.isp_s, nothing)
    sizing.reason === :accepted || return no_burn(sizing.reason)
    burn_time_s = H._constant_thrust_burn_duration_s(mass_kg, planner.delta_v_mps, request.thrust_n, request.isp_s)
    window = H._apoapsis_burn_window(request.time_s, a, e, ν_wrapped, request.mu_m3ps2, burn_time_s)
    window === nothing && return no_burn(:invalid_orbit)
    plan = PropulsiveBurnPlan(valid=true, delta_v_mps=planner.delta_v_mps, direction_rad=planner.direction_rad,
        source_orbit=-1, start_burn_s=window[1], stop_burn_s=window[2], thrust_n=request.thrust_n,
        isp_s=request.isp_s, commanded_impulse_n_s=sizing.commanded_impulse_n_s,
        propellant_required_kg=sizing.propellant_required_kg)
    return GA.AlgorithmResult(request_id=request.request_id, status=:candidate, termination=:planned,
        product=plan, diagnostics=(semimajor_axis_m=a, eccentricity=e, true_anomaly_rad=ν_wrapped))
end

end # module BurnPlanning
