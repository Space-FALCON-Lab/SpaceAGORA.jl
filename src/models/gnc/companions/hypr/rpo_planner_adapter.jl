# Opt-in contract adapter for existing HYPR planning and configured retiming.
module HYPRRPOPlanning
using LinearAlgebra: norm, dot
using StaticArrays: SVector
using ..RPOPlannerInterfaces
using ..SimulationModel
import ..GNCAlgorithms as GA
const P = RPOPlannerInterfaces
const S = SimulationModel
const G = S.GuidanceHooks
import ..RPOPlannerInterfaces: planner_capabilities, plan_rpo!, retime_rpo!

"""
Wrap a copied HYPR configuration. The physical request owns clearance and the
reference interval. Enabled request limits cap the configured retiming limits
using prospective headroom. Other settings and optimizer-returned adaptations
are retained. `rrt_on_replan` explicitly enables the existing warm start for
`:replan` requests. The adapter retimes an active reference (`retime_rpo!`);
restart is not advertised.
"""
struct HYPRRPOPlanner <: P.AbstractRPOPlanner
    config::S.RPOPSOConfig
    headroom::P.RPOPlanningHeadroom
    max_reference_samples::Int
    rrt_on_replan::Bool
    function HYPRRPOPlanner(config::S.RPOPSOConfig; headroom=P.RPOPlanningHeadroom(),
                            max_reference_samples::Integer=100_000, rrt_on_replan::Bool=false)
        2 <= max_reference_samples < typemax(Int) || throw(ArgumentError("Invalid reference sample budget."))
        new(deepcopy(S.rpo_pso_config(config)), headroom, Int(max_reference_samples), rrt_on_replan)
    end
end
planner_capabilities(::HYPRRPOPlanner) = P.RPOPlannerCapabilities(state_sources=(:truth,),
    frames=(:target_rtn,), retiming=true)

function _plan_hypr_rpo! end
function plan_rpo!(state::Nothing, planner::HYPRRPOPlanner, request::P.RPOPlanningRequest, rng::P.AbstractRNG)
    S.HYPRSupport.require_hypr()
    return _plan_hypr_rpo!(state, planner, request, rng)
end

# The route a retime follows: the vertices of `route` from the point nearest
# `current` to the end, with that point replaced by `current`. This is the route
# rule of the legacy model's `rpo_retime_existing_plan`; that function keeps its own
# copy, so the legacy model's runs stay bitwise unchanged until it is removed.
function retiming_route(route_rtn, current_rtn)
    current = SVector{3, Float64}(current_rtn)
    route = G._rpo_retime_route_vertices(route_rtn)
    n = size(route, 2)
    best, segment, projection = Inf, 1, route[:, 1]
    for j in 1:n-1
        a = SVector{3, Float64}(route[:, j]); d = SVector{3, Float64}(route[:, j+1]) - a
        q = a + clamp(dot(current - a, d) / max(dot(d, d), eps(Float64)), 0.0, 1.0) * d
        distance = norm(current - q)
        distance < best && ((best, segment, projection) = (distance, j, collect(q)))
    end
    tail = n == 1 ? Matrix{Float64}(route) : hcat(projection, route[:, segment+1:end])
    # The projection locates progress; it is not a waypoint to visit first.
    norm(current - tail[:, 1]) <= 1.0e-10 && return tail
    return size(tail, 2) == 1 ? hcat(current, tail) : hcat(current, tail[:, 2:end])
end

# Retime `active_reference` for `request` without searching for a new path. The
# route is the active reference's validated polyline from the point nearest the
# chaser to its end, started at the chaser's position (`retiming_route`, the rule
# the legacy replanning model uses). The HYPR retimer then times it as a polyline
# under the request's clearance, step and limits, mapped with headroom exactly as for
# a new plan, starting at the configured initial speed and ending at rest.
#
# Draws no random numbers. Refuses an active reference for other spacecraft or
# another frame (`:unsupported`), one that does not end at the request's goal
# (`:infeasible`, `:goal_not_on_reference`: a changed goal needs a new plan), and,
# like `plan_rpo!`, a lifetime or work budget the timed reference exceeds. The
# adapter applies the unchanged validator before returning a candidate, which still
# requires lifecycle-owned validation.
function retime_rpo!(state::Nothing, planner::HYPRRPOPlanner, request::P.RPOPlanningRequest,
        active_reference::P.RPOReference, rng::P.AbstractRNG)
    S.HYPRSupport.require_hypr()
    reject(status, reason; diagnostics=NamedTuple()) = P.RPOPlanningResult(
        request_id=request.request_id, status=status, termination=reason, diagnostics=diagnostics)
    request.geometry isa S.RPOReferenceGeometry || return reject(:unsupported, :unsupported_geometry)
    active = active_reference
    active.chaser_id == request.chaser_id && active.target_id == request.target_id &&
        active.frame === request.frame || return reject(:unsupported, :reference_mismatch)
    route = active.r_ref_rtn_m
    size(route, 1) == 3 && size(route, 2) >= 1 && all(isfinite, route) ||
        return reject(:failed, :invalid_active_reference)
    norm(route[:, end] - collect(request.goal_rtn_m)) <= request.validation.endpoint_atol_m ||
        return reject(:infeasible, :goal_not_on_reference)
    path = retiming_route(route, collect(request.x_rtn[1:3]))
    return timed_candidate(request, planner.config, planner.headroom, planner.max_reference_samples,
        path, :polyline) do budget, input_cfg
        (planner=:hypr, revision=:retime, planning_budget=budget, requested_config=planner.config,
         input_config=input_cfg, path_kind=:polyline, path=path, active_origin_time_s=active.origin_time_s)
    end
end

# Time `path` for `request` with the HYPR retimer and return a candidate or a refusal.
# The request's clearance, step and limits are mapped onto `cfg` with `headroom` as for
# a new HYPR plan (HYPR's plan_rpo!); `curve_type` says how `path` is read.
# `describe(budget, input_cfg)` gives the candidate's diagnostics. Like plan_rpo!, it
# applies the unchanged validator before returning a candidate, which still requires
# lifecycle-owned validation.
function timed_candidate(describe, request, cfg, headroom, max_reference_samples, path, curve_type)
    reject(status, reason; diagnostics=NamedTuple()) = P.RPOPlanningResult(
        request_id=request.request_id, status=status, termination=reason, diagnostics=diagnostics)
    clearance = request.constraints.clearance_m
    vscale = something(request.constraints.max_speed_mps,
        isfinite(cfg.retime_max_speed_mps) ? cfg.retime_max_speed_mps : norm(collect(request.x_rtn[4:6])))
    budget = P.rpo_planning_budget(request, headroom;
        position_scale_m=max(norm(collect(request.x_rtn[1:3])), norm(collect(request.goal_rtn_m))),
        velocity_scale_mps=max(vscale, cfg.retime_initial_speed_mps))
    budget.supported || return reject(:unsupported, budget.reason; diagnostics=(planning_budget=budget,))
    vmax = budget.max_speed_mps === nothing ? cfg.retime_max_speed_mps : min(cfg.retime_max_speed_mps, budget.max_speed_mps)
    amax = budget.max_acceleration_mps2 === nothing ? cfg.retime_a_max_mps2 : min(cfg.retime_a_max_mps2, budget.max_acceleration_mps2)
    cfg.retime_min_speed_mps <= vmax && cfg.retime_initial_speed_mps <= vmax ||
        return reject(:unsupported, :configured_speed_exceeds_planning_budget)
    request.constraints.max_acceleration_mps2 !== nothing && !cfg.retime_accel_limit_enable &&
        return reject(:unsupported, :acceleration_limited_retimer_required)
    input_cfg = S.rpo_pso_config(cfg; safe_distance_m=clearance, retime_dt_s=request.reference_dt_s,
        retime_max_speed_mps=vmax, retime_a_max_mps2=amax, curve_type=curve_type,
        retime_max_steps=min(cfg.retime_max_steps, max_reference_samples - 1))
    if input_cfg.retime_accel_limit_enable
        # Inspect the profile's duration before allocating its uniform arrays.
        samples, params, clearances = G.rpo_sample_path_with_params(path, input_cfg, request.geometry;
            safe_distance_m=clearance, base_ds_m=G.rpo_retime_sampling_ds_m(input_cfg, clearance),
            curve_type=curve_type)
        profile = G.rpo_retime_profile(G.RPORetimeCurve(path, curve_type), samples, params, clearances,
            request.geometry, input_cfg; safe_distance_m=clearance)
        intervals = profile.duration_s / request.reference_dt_s
        isfinite(intervals) && intervals <= max_reference_samples - 1 ||
            return reject(:failed, :reference_work_limit)
        request.time_s + ceil(intervals) * request.reference_dt_s <= request.valid_until_s ||
            return reject(:infeasible, :insufficient_reference_lifetime)
        timed = G.rpo_retimed_reference_from_profile(profile, input_cfg.retime_dt_s)
        times, positions, velocities = timed.t_s, timed.r_rtn, timed.v_rtn
    else
        times, positions, velocities = G.rpo_reference_from_path(path, request.geometry, input_cfg;
            safe_distance_m=clearance)
    end
    length(times) <= max_reference_samples || return reject(:failed, :reference_work_limit)
    all(isfinite, times) && all(isfinite, positions) && all(isfinite, velocities) && !isempty(times) ||
        return reject(:failed, :nonfinite_reference)
    request.time_s + last(times) <= request.valid_until_s ||
        return reject(:infeasible, :insufficient_reference_lifetime)
    actual_budget = P.rpo_planning_budget(request, headroom;
        position_scale_m=maximum(norm, eachcol(positions)), velocity_scale_mps=maximum(norm, eachcol(velocities)))
    actual_budget.supported ||
        return reject(:unsupported, actual_budget.reason; diagnostics=(planning_budget=actual_budget,))
    ref = P.RPOReference(t_ref_s=times, r_ref_rtn_m=positions, v_ref_rtn_mps=velocities,
        origin_time_s=request.time_s, valid_until_s=request.valid_until_s,
        chaser_id=request.chaser_id, target_id=request.target_id, frame=request.frame,
        geometry_revision=request.geometry_revision)
    diagnostics = describe(actual_budget, input_cfg)
    candidate = P.RPOPlanningResult(request_id=request.request_id, status=:candidate,
        termination=:completed, reference=ref, diagnostics=diagnostics)
    checked = P.validate_rpo_result(request, candidate;
        clearance_at=(p, g) -> S.rpo_clearance_distance_to_station(p, g))
    checked.accepted || return reject(:failed, :reference_rejected;
        diagnostics=merge(diagnostics, (validation=checked,)))
    # Candidate still requires lifecycle-owned validation before installation.
    return candidate
end
P.require_planner_support(::HYPRRPOPlanner) = S.HYPRSupport.require_hypr()
P.initialize_planner(::HYPRRPOPlanner, context) = (S.HYPRSupport.require_hypr(); nothing)

"""
    RetainedPathRPOPlanner(path_rtn, config::RPOPSOConfig; endpoint_tolerance_m=1e-6,
        headroom=RPOPlanningHeadroom(), max_reference_samples=100_000)

Replay a saved geometric path, such as a published HYPR plan, instead of searching.
`path_rtn` holds the path's 3 by N points in target RTN metres (N ≥ 2): the
control polygon for a Bézier `config.curve_type`, the vertices for `:polyline`.
`config` is the HYPR configuration the path was planned with. Both are copied.

For a request whose chaser position is within `endpoint_tolerance_m` of the path's
first point, and whose goal is within it of the last, `plan_rpo!` sets those two
points to the request's and times the path with the HYPR retimer. The request's
clearance, reference step and limits are mapped onto `config` with `headroom`, as
`HYPRRPOPlanner` maps them for a new plan. Otherwise it refuses
(`:infeasible`, `:start_not_on_path` or `:goal_not_on_path`), and it refuses or
fails for the same budget, lifetime and validator reasons as `HYPRRPOPlanner`.
In a run, the timed reference still passes the lifecycle's validation before it is
installed.

It draws no random numbers, declares no retiming or restart, and needs the HYPR
companion.
"""
struct RetainedPathRPOPlanner <: P.AbstractRPOPlanner
    path_rtn::Matrix{Float64}
    config::S.RPOPSOConfig
    endpoint_tolerance_m::Float64
    headroom::P.RPOPlanningHeadroom
    max_reference_samples::Int
    function RetainedPathRPOPlanner(path_rtn::AbstractMatrix, config::S.RPOPSOConfig;
            endpoint_tolerance_m::Real=1e-6, headroom=P.RPOPlanningHeadroom(),
            max_reference_samples::Integer=100_000)
        size(path_rtn, 1) == 3 && size(path_rtn, 2) >= 2 && all(isfinite, path_rtn) ||
            throw(ArgumentError("A retained path needs at least two finite 3D points."))
        2 <= max_reference_samples < typemax(Int) || throw(ArgumentError("Invalid reference sample budget."))
        new(copy(Matrix{Float64}(path_rtn)), deepcopy(S.rpo_pso_config(config)),
            P._nonnegative(endpoint_tolerance_m, "endpoint_tolerance_m"), headroom, Int(max_reference_samples))
    end
end
planner_capabilities(::RetainedPathRPOPlanner) = P.RPOPlannerCapabilities(state_sources=(:truth,),
    frames=(:target_rtn,))
# Timing a fixed path is deterministic.
GA.consumes_rng(::RetainedPathRPOPlanner) = false
P.require_planner_support(::RetainedPathRPOPlanner) = S.HYPRSupport.require_hypr()
P.initialize_planner(::RetainedPathRPOPlanner, context) = (S.HYPRSupport.require_hypr(); nothing)

function plan_rpo!(state::Nothing, planner::RetainedPathRPOPlanner, request::P.RPOPlanningRequest,
        rng::P.AbstractRNG)
    S.HYPRSupport.require_hypr()
    reject(status, reason) = P.RPOPlanningResult(request_id=request.request_id, status=status,
        termination=reason)
    request.reason === :retime && return reject(:unsupported, :retiming_not_implemented)
    request.geometry isa S.RPOReferenceGeometry || return reject(:unsupported, :unsupported_geometry)
    path = copy(planner.path_rtn)
    start, goal = collect(request.x_rtn[1:3]), collect(request.goal_rtn_m)
    norm(path[:, 1] - start) <= planner.endpoint_tolerance_m || return reject(:infeasible, :start_not_on_path)
    norm(path[:, end] - goal) <= planner.endpoint_tolerance_m || return reject(:infeasible, :goal_not_on_path)
    path[:, 1] .= start
    path[:, end] .= goal
    curve = planner.config.curve_type
    return timed_candidate(request, planner.config, planner.headroom, planner.max_reference_samples,
        path, curve) do budget, input_cfg
        (planner=:retained_path, planning_budget=budget, requested_config=planner.config,
         input_config=input_cfg, path_kind=curve === :bezier ? :bezier_control_polygon : :polyline,
         path=path)
    end
end
end
