function _plan_hypr_rpo!(::Nothing, planner::HYPRRPOPlanner, request::P.RPOPlanningRequest, rng::P.AbstractRNG)
    reject(status, reason; diagnostics=NamedTuple()) = P.RPOPlanningResult(
        request_id=request.request_id, status=status, termination=reason, diagnostics=diagnostics)
    request.reason === :retime && return reject(:unsupported, :retiming_not_implemented)
    request.geometry isa S.RPOReferenceGeometry || return reject(:unsupported, :unsupported_geometry)
    cfg = planner.config
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
    input_cfg = S.rpo_pso_config(cfg; safe_distance_m=request.constraints.clearance_m,
        retime_dt_s=request.reference_dt_s, retime_max_speed_mps=vmax, retime_a_max_mps2=amax,
        retime_max_steps=min(cfg.retime_max_steps, planner.max_reference_samples - 1),
        rrt_warmstart_enable=cfg.rrt_warmstart_enable || (planner.rrt_on_replan && request.reason === :replan))
    raw = S.rpo_pso_plan_path(request.x_rtn[1:3], request.goal_rtn_m, request.geometry, input_cfg;
        safe_distance_m=request.constraints.clearance_m, rng=rng)
    # The optimizer may adapt counts, weights and sampling. Use its returned
    # configuration verbatim, never the original configuration, for retiming.
    if raw.config.retime_accel_limit_enable
        # Reuse the same profile construction and evaluator as the existing
        # wrapper, but inspect duration before allocating its uniform arrays.
        samples, params, clearances = G.rpo_sample_path_with_params(raw.path, raw.config, request.geometry;
            safe_distance_m=request.constraints.clearance_m,
            base_ds_m=G.rpo_retime_sampling_ds_m(raw.config, request.constraints.clearance_m),
            curve_type=raw.config.curve_type)
        profile = G.rpo_retime_profile(G.RPORetimeCurve(raw.path, raw.config.curve_type),
            samples, params, clearances, request.geometry, raw.config;
            safe_distance_m=request.constraints.clearance_m)
        intervals = profile.duration_s / request.reference_dt_s
        isfinite(intervals) && intervals <= planner.max_reference_samples - 1 ||
            return reject(:failed, :reference_work_limit)
        request.time_s + ceil(intervals) * request.reference_dt_s <= request.valid_until_s ||
            return reject(:infeasible, :insufficient_reference_lifetime)
        timed = G.rpo_retimed_reference_from_profile(profile, raw.config.retime_dt_s)
        times, positions, velocities = timed.t_s, timed.r_rtn, timed.v_rtn
    else
        times, positions, velocities = G.rpo_reference_from_path(raw.path, request.geometry, raw.config;
            safe_distance_m=request.constraints.clearance_m)
    end
    length(times) <= planner.max_reference_samples || return reject(:failed, :reference_work_limit)
    all(isfinite, times) && all(isfinite, positions) && all(isfinite, velocities) && !isempty(times) ||
        return reject(:failed, :nonfinite_reference)
    # Legacy retiming has no duration profile to inspect before allocation.
    # Apply the same strict lifetime budget as acceleration-limited retiming
    # before building a candidate, including shortages within validator roundoff.
    request.time_s + last(times) <= request.valid_until_s ||
        return reject(:infeasible, :insufficient_reference_lifetime)
    actual_budget = P.rpo_planning_budget(request, planner.headroom;
        position_scale_m=maximum(norm, eachcol(positions)), velocity_scale_mps=maximum(norm, eachcol(velocities)))
    actual_budget.supported || return reject(:unsupported, actual_budget.reason; diagnostics=(planning_budget=actual_budget,))
    ref = P.RPOReference(t_ref_s=times, r_ref_rtn_m=positions, v_ref_rtn_mps=velocities,
        origin_time_s=request.time_s, valid_until_s=request.valid_until_s,
        chaser_id=request.chaser_id, target_id=request.target_id, frame=request.frame,
        geometry_revision=request.geometry_revision)
    diagnostics = (planner=:hypr, planning_budget=actual_budget,
        requested_config=cfg, input_config=input_cfg, effective_config=raw.config,
        path_kind=raw.config.curve_type === :bezier ? :bezier_control_polygon : :polyline,
        path=raw.path, cost=raw.cost, cost_history=raw.cost_history, adaptive=raw.adaptive,
        warmstart=raw.warmstart, early_stopped=raw.early_stopped,
        iteration_timed_out=raw.iteration_timed_out)
    termination = raw.iteration_timed_out ? :time_budget : raw.early_stopped ? :completed : :iteration_limit
    candidate = P.RPOPlanningResult(request_id=request.request_id, status=:candidate,
        termination=termination, reference=ref, diagnostics=diagnostics)
    checked = P.validate_rpo_result(request, candidate;
        clearance_at=(p,g)->S.rpo_clearance_distance_to_station(p,g))
    checked.accepted || return reject(:failed, :reference_rejected;
        diagnostics=merge(diagnostics, (validation=checked,)))
    # Candidate still requires lifecycle-owned validation before installation.
    return candidate
end
