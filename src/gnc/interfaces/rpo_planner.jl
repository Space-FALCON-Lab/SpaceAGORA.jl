"""
Internal RPO planner contracts. This module loads with Julia standard libraries
alone. It is not yet connected to guidance or a stable top-level package API.
"""
module RPOPlannerInterfaces

using LinearAlgebra: norm, cross
using Random: AbstractRNG

"""An algorithm producing a timed target-relative RPO reference."""
abstract type AbstractRPOPlanner end

_finite(x, name) = isfinite(x) ? x : throw(ArgumentError("$name must be finite."))
function _nonnegative(x, name)
    y = _finite(Float64(x), name)
    y >= 0 || throw(ArgumentError("$name must be nonnegative."))
    return y
end
function _positive(x, name)
    y = _nonnegative(x, name)
    y > 0 || throw(ArgumentError("$name must be positive."))
    return y
end
function _state_tuple(x, ::Val{N}, name) where {N}
    length(x) == N || throw(ArgumentError("$name must contain $N values."))
    return ntuple(i -> _finite(Float64(x[i]), name), N)
end
_nonempty(x, name) = isempty(x) ? throw(ArgumentError("$name must not be empty.")) : String(x)

"""Explicit SI clearance and optional discrete reference speed/acceleration limits."""
struct RPOPlanningConstraints
    clearance_m::Float64
    max_speed_mps::Union{Nothing,Float64}
    max_acceleration_mps2::Union{Nothing,Float64}
    function RPOPlanningConstraints(; clearance_m, max_speed_mps=nothing,
                                    max_acceleration_mps2=nothing)
        new(_nonnegative(clearance_m, "clearance_m"),
            isnothing(max_speed_mps) ? nothing : _positive(max_speed_mps, "max_speed_mps"),
            isnothing(max_acceleration_mps2) ? nothing :
                _positive(max_acceleration_mps2, "max_acceleration_mps2"))
    end
end

const _RPO_MAX_LIMIT_ROUNDOFF_RTOL = 128 * eps(Float64)

"""
Prospective software validation settings. Clearance is sampled on the reference
polyline at spacing at most `clearance_sample_ds_m`, including endpoints. This
is not continuous collision certification. The work cap rejects, never skips,
checks that would exceed `max_clearance_samples`. Speed/acceleration comparisons
allow at most `128eps(Float64)` relative roundoff, with no absolute floor. Set
`limit_roundoff_rtol=0` for exact comparisons. Physical planning margin belongs
in the adapter, not in this bounded software allowance.
"""
struct RPOValidationSettings
    endpoint_atol_m::Float64
    time_atol_s::Float64
    limit_roundoff_rtol::Float64
    clearance_sample_ds_m::Float64
    max_clearance_samples::Int
    allow_time_budget_candidate::Bool
    function RPOValidationSettings(; endpoint_atol_m=1e-9, time_atol_s=1e-10,
                                   limit_roundoff_rtol=_RPO_MAX_LIMIT_ROUNDOFF_RTOL,
                                   clearance_sample_ds_m=0.05,
                                   max_clearance_samples::Integer=100_000,
                                   allow_time_budget_candidate::Bool=false)
        max_clearance_samples >= 2 || throw(ArgumentError("At least two clearance samples are required."))
        roundoff = _nonnegative(limit_roundoff_rtol, "limit_roundoff_rtol")
        roundoff <= _RPO_MAX_LIMIT_ROUNDOFF_RTOL ||
            throw(ArgumentError("Limit roundoff tolerance cannot exceed 128eps(Float64)."))
        new(_nonnegative(endpoint_atol_m, "endpoint_atol_m"),
            _nonnegative(time_atol_s, "time_atol_s"), roundoff,
            _positive(clearance_sample_ds_m, "clearance_sample_ds_m"),
            Int(max_clearance_samples), allow_time_budget_candidate)
    end
end

"""
Owned planning snapshot in SI units. `x_rtn` is rotating-frame relative position
and velocity. Geometry must be static and target-RTN aligned in this first
contract. The caller retains a separate snapshot for validation; pass a deep
copy to untrusted planner code. No integrator or live model belongs here.
"""
struct RPOPlanningRequest{E,G}
    request_id::Int
    reason::Symbol
    chaser_id::Int
    target_id::Int
    epoch::E
    time_s::Float64
    observation_time_s::Float64
    state_source::Symbol
    frame::Symbol
    x_rtn::NTuple{6,Float64}
    target_state_ii::NTuple{6,Float64}
    goal_rtn_m::NTuple{3,Float64}
    geometry::G
    geometry_revision::String
    constraints::RPOPlanningConstraints
    reference_dt_s::Float64
    preview_horizon_steps::Int
    valid_until_s::Float64
    validation::RPOValidationSettings
    function RPOPlanningRequest(; request_id::Integer, reason::Symbol=:initial,
            chaser_id::Integer, target_id::Integer, epoch, time_s,
            observation_time_s=time_s, state_source::Symbol=:truth,
            frame::Symbol=:target_rtn, x_rtn, target_state_ii, goal_rtn_m,
            geometry, geometry_revision, constraints::RPOPlanningConstraints,
            reference_dt_s, preview_horizon_steps::Integer, valid_until_s,
            validation::RPOValidationSettings=RPOValidationSettings())
        request_id > 0 || throw(ArgumentError("request_id must be positive."))
        chaser_id > 0 && target_id > 0 && chaser_id != target_id ||
            throw(ArgumentError("Distinct positive spacecraft IDs are required."))
        reason in (:initial, :forced, :replan, :retime) || throw(ArgumentError("Unknown planning reason."))
        state_source === :truth || throw(ArgumentError("Only truth observations are supported."))
        frame === :target_rtn || throw(ArgumentError("Only target_rtn geometry is supported."))
        t = _finite(Float64(time_s), "time_s")
        obs = _finite(Float64(observation_time_s), "observation_time_s")
        dt = _positive(reference_dt_s, "reference_dt_s")
        validation.time_atol_s < dt / 2 || throw(ArgumentError("Time tolerance must be less than half a reference step."))
        abs(obs - t) <= validation.time_atol_s || throw(ArgumentError("Stale observation."))
        stop = _finite(Float64(valid_until_s), "valid_until_s")
        preview_horizon_steps > 0 || throw(ArgumentError("Preview horizon must be positive."))
        t + dt > t || throw(ArgumentError("Reference step is below simulation-time resolution."))
        preview_end = t + preview_horizon_steps * dt
        isfinite(preview_end) && stop > t && preview_end <= stop ||
            throw(ArgumentError("Validity must cover the initial control preview."))
        target = _state_tuple(target_state_ii, Val(6), "target_state_ii")
        r, v = collect(target[1:3]), collect(target[4:6])
        h = norm(cross(r, v))
        isfinite(norm(r)) && norm(r) > eps(Float64) && isfinite(h) && h > eps(Float64) ||
            throw(ArgumentError("Target state cannot define an RTN frame."))
        ep, geom = deepcopy(epoch), deepcopy(geometry)
        new{typeof(ep),typeof(geom)}(Int(request_id), reason, Int(chaser_id), Int(target_id),
            ep, t, obs, state_source, frame, _state_tuple(x_rtn, Val(6), "x_rtn"),
            target, _state_tuple(goal_rtn_m, Val(3), "goal_rtn_m"), geom,
            _nonempty(geometry_revision, "geometry_revision"), constraints, dt,
            Int(preview_horizon_steps), stop, validation)
    end
end

"""
Owned candidate reference. Construction copies arrays; validation, rather than
construction, diagnoses malformed algorithm output. Times are relative to
`origin_time_s`. The only initial terminal policy repeats the last position AND
velocity sample; it is not a stationary-hold guarantee.
"""
struct RPOReference
    t_ref_s::Vector{Float64}
    r_ref_rtn_m::Matrix{Float64}
    v_ref_rtn_mps::Matrix{Float64}
    origin_time_s::Float64
    valid_until_s::Float64
    chaser_id::Int
    target_id::Int
    frame::Symbol
    geometry_revision::String
    terminal_policy::Symbol
    function RPOReference(; t_ref_s, r_ref_rtn_m, v_ref_rtn_mps, origin_time_s,
            valid_until_s, chaser_id::Integer, target_id::Integer,
            geometry_revision, frame::Symbol=:target_rtn,
            terminal_policy::Symbol=:repeat_last_sample)
        new(copy(Vector{Float64}(t_ref_s)), copy(Matrix{Float64}(r_ref_rtn_m)),
            copy(Matrix{Float64}(v_ref_rtn_mps)), Float64(origin_time_s),
            Float64(valid_until_s), Int(chaser_id), Int(target_id), frame,
            String(geometry_revision), terminal_policy)
    end
end

"""Algorithm outcome; `:candidate` still requires independent validation."""
struct RPOPlanningResult{D}
    request_id::Int
    status::Symbol
    termination::Symbol
    reference::Union{Nothing,RPOReference}
    diagnostics::D
    function RPOPlanningResult(; request_id::Integer, status::Symbol,
            termination::Symbol, reference::Union{Nothing,RPOReference}=nothing,
            diagnostics=NamedTuple())
        status in (:candidate, :infeasible, :unsupported, :failed) ||
            throw(ArgumentError("Unknown planner result status."))
        (status === :candidate) == (reference !== nothing) ||
            throw(ArgumentError("Only candidates must contain a reference."))
        # Diagnostics can contain an explicitly tagged algorithm path. They are
        # never interpreted as geometry or used by the controller/validator.
        d = deepcopy(diagnostics)
        new{typeof(d)}(Int(request_id), status, termination, deepcopy(reference), d)
    end
end

"""Capability declaration; an unimplemented planner supports nothing by default."""
Base.@kwdef struct RPOPlannerCapabilities
    state_sources::Tuple{Vararg{Symbol}} = ()
    frames::Tuple{Vararg{Symbol}} = ()
    retiming::Bool = false
    restart::Bool = false
end

"""Declare planner capabilities without loading algorithm configuration."""
planner_capabilities(::AbstractRPOPlanner) = RPOPlannerCapabilities()
"""Create fresh run-owned planner state. Stateless planners use `nothing`."""
initialize_planner(::AbstractRPOPlanner, run_context) = nothing
"""Generate a candidate with an explicitly supplied random stream."""
function plan_rpo!(state, planner::AbstractRPOPlanner, request::RPOPlanningRequest, rng::AbstractRNG)
    return RPOPlanningResult(request_id=request.request_id, status=:unsupported,
        termination=:planning_not_implemented)
end
"""Optional retiming hook. Unsupported retiming does not silently replan."""
function retime_rpo!(state, planner::AbstractRPOPlanner, request::RPOPlanningRequest,
        active_reference::RPOReference, rng::AbstractRNG)
    return RPOPlanningResult(request_id=request.request_id, status=:unsupported,
        termination=:retiming_not_implemented)
end

"""Pure validation outcome with a machine-readable reason and measured values."""
struct RPOValidationResult{M}
    accepted::Bool
    reason::Symbol
    metrics::M
end
_reject(reason; metrics...) = RPOValidationResult(false, reason, (; metrics...))

"""Check declared capabilities. This contract does not implement run startup."""
function validate_rpo_capabilities(planner::AbstractRPOPlanner, request::RPOPlanningRequest;
                                   checkpoint::Bool=false)
    caps = planner_capabilities(planner)
    request.state_source in caps.state_sources || return _reject(:unsupported_state_source)
    request.frame in caps.frames || return _reject(:unsupported_frame)
    request.reason === :retime && !caps.retiming && return _reject(:unsupported_retiming)
    checkpoint && !caps.restart && return _reject(:unsupported_restart)
    return RPOValidationResult(true, :supported, NamedTuple())
end

"""Check current and forecast times against the candidate's declared lifetime."""
function rpo_reference_is_current(request::RPOPlanningRequest, reference::RPOReference, time_s)
    t, tol = Float64(time_s), request.validation.time_atol_s
    preview_end = t + request.preview_horizon_steps * request.reference_dt_s
    return isfinite(t) && isfinite(preview_end) &&
        isfinite(reference.origin_time_s) && isfinite(reference.valid_until_s) &&
        t >= reference.origin_time_s - tol &&
        preview_end <= min(reference.valid_until_s, request.valid_until_s) + tol
end

include("rpo_reference_validation.jl")
end
