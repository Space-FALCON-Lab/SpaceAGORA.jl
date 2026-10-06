"""Extract position and velocity vectors for one spacecraft from the simulation state."""
function _rpo_state_pos_vel(u, idx::Int)
    return SVector{3, Float64}(u.sc[idx].pos), SVector{3, Float64}(u.sc[idx].vel)
end

"""Build a fresh RPO plan from the model goal and the current spacecraft state."""
function build_rpo_plan(model::RPOGuidanceModel, u, t::Float64)
    model.geometry === nothing && throw(ArgumentError("RPOGuidanceModel requires reference geometry."))
    r_chaser, v_chaser = _rpo_state_pos_vel(u, model.chaser_idx)
    r_target, v_target = _rpo_state_pos_vel(u, model.target_idx)
    x_rel = inertial_to_rtn_relative_state(r_chaser, v_chaser, r_target, v_target)
    return build_rpo_plan_from_start(model, x_rel[1:3], model.geometry, t)
end

"""Build an RPO plan from an explicit start state and geometry snapshot."""
function build_rpo_plan_from_start(model::RPOGuidanceModel, start_rtn, geometry, t::Float64; safe_distance_override=nothing, force_rrt_warmstart::Bool=false)
    start = SVector{3, Float64}(start_rtn)
    base_cfg = model.pso_config === nothing ? RPOPSOConfig() : model.pso_config
    force_rrt_warmstart && (base_cfg = rpo_pso_config(base_cfg; rrt_warmstart_enable=true))
    safe_distance = safe_distance_override === nothing ?
        (model.safe_distance_m > 0.0 ? model.safe_distance_m : base_cfg.safe_distance_m) :
        Float64(safe_distance_override)
    plan_result = rpo_pso_plan_path(
        start,
        model.goal_rtn,
        geometry,
        base_cfg;
        safe_distance_m=safe_distance,
    )
    t_ref, r_ref, v_ref = rpo_reference_from_path(
        plan_result.path,
        geometry,
        plan_result.config;
        safe_distance_m=safe_distance,
    )
    return RPOPlan(
        valid=true,
        t_ref_s=t_ref,
        r_ref_rtn=r_ref,
        v_ref_rtn=v_ref,
        path_rtn=plan_result.path,
        cost=plan_result.cost,
        diagnostics=(
            components=plan_result.components,
            adaptive=plan_result.adaptive,
            refinement_improved=plan_result.refinement_improved,
            early_stopped=plan_result.early_stopped,
            early_stop_iter=plan_result.early_stop_iter,
            iteration_timed_out=plan_result.iteration_timed_out,
            iteration_timeout_iter=plan_result.iteration_timeout_iter,
            iteration_timeout_phase=plan_result.iteration_timeout_phase,
            iteration_timeout_events=plan_result.iteration_timeout_events,
            warmstart=plan_result.warmstart,
            planned_at_s=t,
        ),
    )
end

"""Return the replanning configuration attached to an RPO guidance model."""
function _rpo_replanning_config(model::RPOGuidanceModel)
    if model.replanning_config === nothing
        if model.force_replan || model.replanning_phase != :tracking
            safe = model.safe_distance_m > 0.0 ? model.safe_distance_m :
                (model.pso_config === nothing ? 0.0 : model.pso_config.safe_distance_m)
            return RPOReplanningConfig(safe_distance_m=safe, retime_clearance_m=max(0.25, safe))
        end
        return nothing
    end
    model.replanning_config isa RPOReplanningConfig && return model.replanning_config
    return RPOReplanningConfig(; model.replanning_config...)
end

"""Record a replanning action in the model history."""
function _rpo_record_replanning_event!(model::RPOGuidanceModel, action::Symbol, decision, t::Real)
    push!(
        model.replanning_events,
        (
            time_s=Float64(t),
            action=action,
            reason=decision.reason,
            min_clearance_m=decision.min_clearance,
            active_spheres=length(decision.spheres),
        ),
    )
    model.last_replanning_time_s = Float64(t)
    return model
end

"""Command a fixed RTN position and zero velocity; MPC physically brakes into this hold."""
function _rpo_begin_replanning_hold!(model::RPOGuidanceModel, start, decision, t)
    hold = RPOPlan(valid=true, t_ref_s=[0.0], r_ref_rtn=reshape(collect(start), 3, 1),
        v_ref_rtn=zeros(3, 1), path_rtn=reshape(collect(start), 3, 1),
        diagnostics=(replanning_action=:hold,))
    update_rpo_plan_buffer!(model.plan_buffer, hold, t)
    model.pending_replan = (decision=decision,)
    model.replanning_phase = :braking
    model.safe_hold_count += 1
    model.force_replan = false
    _rpo_record_replanning_event!(model, :brake, decision, t)
    return true
end

"""Advance brake/hold/plan/release, replaying measured compute latency in simulation time."""
function _rpo_advance_replanning_hold!(model::RPOGuidanceModel, x, t; planner=build_rpo_plan_from_start)
    if model.replanning_phase == :hold_failed
        return model.force_replan ?
            _rpo_begin_replanning_hold!(model, x[1:3], model.pending_replan.decision, t) : false
    end
    config = _rpo_replanning_config(model)
    pending = model.pending_replan
    safe_distance = config.safe_distance_m > 0.0 ? config.safe_distance_m :
        (model.safe_distance_m > 0.0 ? model.safe_distance_m :
            (model.pso_config === nothing ? 0.0 : model.pso_config.safe_distance_m))
    settled = norm(x[1:3] - model.plan_buffer.plan.r_ref_rtn[:, end]) <= model.hold_position_tolerance_m &&
        norm(x[4:6]) <= model.hold_speed_tolerance_mps
    settled || return false
    model.replanning_phase == :planning && t < pending.available_at_s && return false
    spheres = rpo_active_replanning_spheres(config, t)
    geometry = rpo_geometry_with_replanning_spheres(model.geometry, spheres;
        sphere_surface_samples=config.sphere_surface_samples)
    decision = merge(pending.decision, (geometry=geometry, spheres=spheres))
    if model.replanning_phase == :braking
        original_goal = model.goal_rtn
        config.desired_goal_rtn !== nothing && (model.goal_rtn = config.desired_goal_rtn)
        target_goal = model.goal_rtn
        plan = nothing
        failure = nothing
        runtime_s = @elapsed try
            plan = planner(model, x[1:3], geometry, t; safe_distance_override=safe_distance,
                force_rrt_warmstart=true)
        catch err
            failure = err
        finally
            model.goal_rtn = original_goal
        end
        model.pending_replan = (decision=decision, plan=plan, failure=failure,
            goal=target_goal, runtime_s=runtime_s, available_at_s=t + runtime_s)
        model.replanning_phase = :planning
        _rpo_record_replanning_event!(model, :hold, decision, t)
        return false
    end
    plan = pending.plan
    usable = pending.failure === nothing && plan !== nothing && plan.valid && size(plan.r_ref_rtn, 2) > 0
    failure_reason = pending.failure === nothing ? :invalid_plan : :planner_exception
    if usable && config.desired_goal_rtn !== nothing &&
            norm(config.desired_goal_rtn - pending.goal) > config.goal_change_tolerance_m
        usable = false
        failure_reason = :goal_changed_during_planning
    end
    if usable
        # Recheck against the current map, including the connection from the held state.
        path = hcat(x[1:3], plan.r_ref_rtn)
        clearance_ok = all(rpo_capsule_clearance_to_station(view(path, :, j), view(path, :, j + 1), geometry) >=
            safe_distance for j in 1:(size(path, 2) - 1))
        start_ok = norm(x[1:3] - plan.r_ref_rtn[:, 1]) <= model.hold_position_tolerance_m
        goal_ok = norm(plan.r_ref_rtn[:, end] - pending.goal) <= model.hold_position_tolerance_m
        usable = clearance_ok && start_ok && goal_ok
        failure_reason = !clearance_ok ? :unsafe_replacement : !start_ok ? :hold_drift : :goal_not_reached
    end
    if !usable
        model.replanning_phase = :hold_failed
        model.replan_failure_count += 1
        _rpo_record_replanning_event!(model, :replan_failed, merge(decision, (reason=failure_reason,)), t)
        @warn "RPO replanning failed; maintaining position hold." reason=failure_reason failure=pending.failure
        return false
    end
    update_rpo_plan_buffer!(model.plan_buffer, plan, t)
    model.goal_rtn = pending.goal
    model.replanning_phase = :tracking
    model.pending_replan = nothing
    model.replan_count += 1
    _rpo_record_replanning_event!(model, :replan, decision, t)
    return true
end

"""Retime live, ignore live, or brake and hold before computing a replacement route."""
function maybe_update_rpo_replanning!(model::RPOGuidanceModel, u, t::Float64)
    config = _rpo_replanning_config(model)
    config === nothing && return false
    config.enabled || model.force_replan || model.replanning_phase != :tracking || return false
    model.plan_buffer.valid || return false

    r_chaser, v_chaser = _rpo_state_pos_vel(u, model.chaser_idx)
    r_target, v_target = _rpo_state_pos_vel(u, model.target_idx)
    x_rel = inertial_to_rtn_relative_state(r_chaser, v_chaser, r_target, v_target)
    if model.replanning_phase != :tracking
        return _rpo_advance_replanning_hold!(model, x_rel, t)
    end
    start = x_rel[1:3]
    decision = rpo_replanning_decision(model.plan_buffer.plan, start, model.geometry, config, t)
    model.force_replan && (decision = merge(decision, (action=:replan, reason=:forced)))
    signature = rpo_replanning_signature(decision.spheres)
    if decision.action == :none
        model.last_replanning_signature = signature
        model.replanning_persistence_count = 0
        return false
    end

    if signature == model.last_replanning_signature
        model.replanning_persistence_count += 1
    else
        model.last_replanning_signature = signature
        model.replanning_persistence_count = 1
    end
    if !model.force_replan
        model.replanning_persistence_count < config.hysteresis_samples && return false
        t - model.last_replanning_time_s < config.min_replan_interval_s && return false
    end

    base_cfg = model.pso_config === nothing ? RPOPSOConfig() : model.pso_config
    safe_distance = config.safe_distance_m > 0.0 ? config.safe_distance_m :
        (model.safe_distance_m > 0.0 ? model.safe_distance_m : base_cfg.safe_distance_m)

    if decision.action == :retime
        plan = rpo_retime_existing_plan(model.plan_buffer.plan, start, decision.geometry, base_cfg, safe_distance, t)
        update_rpo_plan_buffer!(model.plan_buffer, plan, t)
        model.retime_count += 1
        _rpo_record_replanning_event!(model, :retime, decision, t)
        return true
    elseif decision.action == :replan
        return _rpo_begin_replanning_hold!(model, start, decision, t)
    end
    return false
end

"""Advance guidance-side plan state for a guidance model during simulation."""
function calcGuidanceEffect!(model::RPOGuidanceModel, u, p, t::Float64, sat_idx::Int)
    sat_idx == model.chaser_idx || return nothing
    if !model.plan_buffer.valid
        plan = build_rpo_plan(model, u, t)
        update_rpo_plan_buffer!(model.plan_buffer, plan, t)
        model.force_replan = false
    else
        maybe_update_rpo_replanning!(model, u, t)
    end
    return nothing
end
