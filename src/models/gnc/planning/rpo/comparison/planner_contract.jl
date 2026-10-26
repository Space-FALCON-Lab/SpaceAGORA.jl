# RPO planner contract methods of the seven comparison planners. The objects and
# their numerical routes live with the comparison harness (planner_objects.jl).
module RPOComparisonPlanning
using LinearAlgebra: norm
using StaticArrays: SVector
using ..RPOPlannerInterfaces
using ..SimulationModel
const P = RPOPlannerInterfaces
const S = SimulationModel
const G = S.GuidanceHooks
import ..RPOPlannerInterfaces: planner_capabilities, plan_rpo!

const ComparisonPlanner = G.AbstractRPOComparisonPlanner
const SwarmPlanner = Union{G.HYPRComparisonPlanner, G.UnrefinedPSOPlanner}
const SearchPlanner = Union{G.RRTConnectPlanner, G.RRTConnectBezierPlanner, G.RRTStarPlanner}

# Truth-state, target-RTN planning only: no retiming and no restart.
planner_capabilities(::ComparisonPlanner) = P.RPOPlannerCapabilities(state_sources=(:truth,), frames=(:target_rtn,))
P.require_planner_support(::ComparisonPlanner) = S.HYPRSupport.require_hypr()
P.initialize_planner(::ComparisonPlanner, context) = (S.HYPRSupport.require_hypr(); nothing)

_termination(::SwarmPlanner, raw, cap) = raw.iteration_timed_out ? :time_budget : raw.early_stopped ? :completed : :iteration_limit
_termination(::SearchPlanner, raw, cap) = :completed
_termination(::Union{G.CHOMPPlanner, G.STOMPPlanner}, raw, cap) = raw.iterations >= cap ? :iteration_limit : :completed
_iteration_cap(::Union{SwarmPlanner, SearchPlanner}, cfg) = 0
_iteration_cap(planner::Union{G.CHOMPPlanner, G.STOMPPlanner}, cfg) =
    G._rpo_comparison_iterations(planner.optimizer, cfg, Inf, planner.settings.n_iters)

"""
Plan with the comparison mode's kernels under the request's clearance, reference
step, lifetime and limits, never the comparison's fixed 0.5 m safety distance.
Request limits with headroom cap the configured retiming limits as for
`HYPRRPOPlanner`, and the returned path is timed with the configuration the
kernel returned. The planner works without a wall-clock limit. A search that finds
no path is `:failed/:path_not_found`; it never returns the direct-segment
diagnostic as a candidate. The unchanged validator runs before a candidate is
returned, which still requires lifecycle-owned validation.
"""
function plan_rpo!(state::Nothing, planner::ComparisonPlanner, request::P.RPOPlanningRequest, rng::P.AbstractRNG)
    S.HYPRSupport.require_hypr()
    reject(status, reason; diagnostics=NamedTuple()) = P.RPOPlanningResult(
        request_id=request.request_id, status=status, termination=reason, diagnostics=diagnostics)
    request.reason === :retime && return reject(:unsupported, :retiming_not_implemented)
    request.geometry isa S.RPOReferenceGeometry || return reject(:unsupported, :unsupported_geometry)
    cfg = planner.config
    clearance = request.constraints.clearance_m
    vscale = something(request.constraints.max_speed_mps,
        isfinite(cfg.retime_max_speed_mps) ? cfg.retime_max_speed_mps : norm(collect(request.x_rtn[4:6])))
    budget = P.rpo_planning_budget(request, planner.headroom;
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
        retime_max_speed_mps=vmax, retime_a_max_mps2=amax,
        retime_max_steps=min(cfg.retime_max_steps, planner.max_reference_samples - 1))
    id = G.normalize_rpo_comparison_planner_type(planner)
    raw, iterations = G._rpo_comparison_planner_path(planner, SVector{3, Float64}(request.x_rtn[1:3]),
        SVector{3, Float64}(request.goal_rtn_m), request.geometry, input_cfg;
        safe_distance_m=clearance, rng=rng, runtime_limit_s=Inf, warmstart_path=nothing)
    diagnostics = (planner=id, planning_budget=budget, requested_config=cfg, input_config=input_cfg,
        effective_config=raw.config,
        path_kind=raw.config.curve_type === :bezier ? :bezier_control_polygon : :polyline,
        path=raw.path, cost=raw.cost, cost_history=raw.cost_history, iterations=iterations)
    hasproperty(raw, :path_found) && raw.path_found === false &&
        return reject(:failed, :path_not_found; diagnostics=diagnostics)
    # Time the path with the configuration the kernel returned, never the original.
    if raw.config.retime_accel_limit_enable
        # Inspect the profile's duration before allocating its uniform arrays.
        samples, params, clearances = G.rpo_sample_path_with_params(raw.path, raw.config, request.geometry;
            safe_distance_m=clearance, base_ds_m=G.rpo_retime_sampling_ds_m(raw.config, clearance),
            curve_type=raw.config.curve_type)
        profile = G.rpo_retime_profile(G.RPORetimeCurve(raw.path, raw.config.curve_type), samples, params,
            clearances, request.geometry, raw.config; safe_distance_m=clearance)
        intervals = profile.duration_s / request.reference_dt_s
        isfinite(intervals) && intervals <= planner.max_reference_samples - 1 ||
            return reject(:failed, :reference_work_limit; diagnostics=diagnostics)
        request.time_s + ceil(intervals) * request.reference_dt_s <= request.valid_until_s ||
            return reject(:infeasible, :insufficient_reference_lifetime; diagnostics=diagnostics)
        timed = G.rpo_retimed_reference_from_profile(profile, raw.config.retime_dt_s)
        times, positions, velocities = timed.t_s, timed.r_rtn, timed.v_rtn
    else
        times, positions, velocities = G.rpo_reference_from_path(raw.path, request.geometry, raw.config;
            safe_distance_m=clearance)
    end
    length(times) <= planner.max_reference_samples ||
        return reject(:failed, :reference_work_limit; diagnostics=diagnostics)
    all(isfinite, times) && all(isfinite, positions) && all(isfinite, velocities) && !isempty(times) ||
        return reject(:failed, :nonfinite_reference; diagnostics=diagnostics)
    request.time_s + last(times) <= request.valid_until_s ||
        return reject(:infeasible, :insufficient_reference_lifetime; diagnostics=diagnostics)
    actual_budget = P.rpo_planning_budget(request, planner.headroom;
        position_scale_m=maximum(norm, eachcol(positions)), velocity_scale_mps=maximum(norm, eachcol(velocities)))
    diagnostics = merge(diagnostics, (planning_budget=actual_budget,))
    actual_budget.supported || return reject(:unsupported, actual_budget.reason; diagnostics=diagnostics)
    ref = P.RPOReference(t_ref_s=times, r_ref_rtn_m=positions, v_ref_rtn_mps=velocities,
        origin_time_s=request.time_s, valid_until_s=request.valid_until_s,
        chaser_id=request.chaser_id, target_id=request.target_id, frame=request.frame,
        geometry_revision=request.geometry_revision)
    candidate = P.RPOPlanningResult(request_id=request.request_id, status=:candidate,
        termination=_termination(planner, raw, _iteration_cap(planner, input_cfg)), reference=ref,
        diagnostics=diagnostics)
    checked = P.validate_rpo_result(request, candidate;
        clearance_at=(p, g) -> S.rpo_clearance_distance_to_station(p, g))
    checked.accepted || return reject(:failed, :reference_rejected;
        diagnostics=merge(diagnostics, (validation=checked,)))
    # Candidate still requires lifecycle-owned validation before installation.
    return candidate
end
end
