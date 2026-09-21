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
    model.replanning_config === nothing && return nothing
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

"""Update the active RPO plan when dynamic obstacles, goals, or tracking errors require replanning."""
function maybe_update_rpo_replanning!(model::RPOGuidanceModel, u, t::Float64)
    config = _rpo_replanning_config(model)
    config === nothing && return false
    config.enabled || return false
    model.plan_buffer.valid || return false

    r_chaser, v_chaser = _rpo_state_pos_vel(u, model.chaser_idx)
    r_target, v_target = _rpo_state_pos_vel(u, model.target_idx)
    x_rel = inertial_to_rtn_relative_state(r_chaser, v_chaser, r_target, v_target)
    start = x_rel[1:3]
    decision = rpo_replanning_decision(model.plan_buffer.plan, start, model.geometry, config, t)
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
    model.replanning_persistence_count < config.hysteresis_samples && return false
    t - model.last_replanning_time_s < config.min_replan_interval_s && return false

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
        original_goal = model.goal_rtn
        goal_was_changed = false
        try
            if config.desired_goal_rtn !== nothing
                model.goal_rtn = config.desired_goal_rtn
                goal_was_changed = true
            end
            plan = build_rpo_plan_from_start(
                model,
                start,
                decision.geometry,
                t;
                safe_distance_override=safe_distance,
                force_rrt_warmstart=true,
            )
            update_rpo_plan_buffer!(model.plan_buffer, plan, t)
            model.replan_count += 1
            _rpo_record_replanning_event!(model, :replan, decision, t)
            return true
        catch err
            goal_was_changed && (model.goal_rtn = original_goal)
            model.replan_failure_count += 1
            _rpo_record_replanning_event!(model, :replan_failed, decision, t)
            @warn "RPO replanning failed; keeping the previous plan." exception=(err, catch_backtrace())
        end
    end
    return false
end

"""Advance guidance-side plan state for a guidance model during simulation."""
function calcGuidanceEffect!(model::RPOGuidanceModel, u, p, t::Float64, sat_idx::Int)
    sat_idx == model.chaser_idx || return nothing
    if !model.plan_buffer.valid || model.force_replan
        plan = build_rpo_plan(model, u, t)
        update_rpo_plan_buffer!(model.plan_buffer, plan, t)
        model.force_replan = false
    else
        maybe_update_rpo_replanning!(model, u, t)
    end
    return nothing
end
