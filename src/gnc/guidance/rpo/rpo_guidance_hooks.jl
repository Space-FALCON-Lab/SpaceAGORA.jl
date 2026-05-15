function _rpo_state_pos_vel(u, idx::Int)
    return SVector{3, Float64}(u.sc[idx].pos), SVector{3, Float64}(u.sc[idx].vel)
end

function build_rpo_plan(model::RPOGuidanceModel, u, t::Float64)
    model.geometry === nothing && throw(ArgumentError("RPOGuidanceModel requires reference geometry."))
    r_chaser, v_chaser = _rpo_state_pos_vel(u, model.chaser_idx)
    r_target, v_target = _rpo_state_pos_vel(u, model.target_idx)
    x_rel = inertial_to_rtn_relative_state(r_chaser, v_chaser, r_target, v_target)
    start = x_rel[1:3]
    base_cfg = model.pso_config === nothing ? RPOPSOConfig() : model.pso_config
    safe_distance = model.safe_distance_m > 0.0 ? model.safe_distance_m : base_cfg.safe_distance_m
    plan_result = rpo_pso_plan_path(
        start,
        model.goal_rtn,
        model.geometry,
        base_cfg;
        safe_distance_m=safe_distance,
    )
    t_ref, r_ref, v_ref = rpo_reference_from_path(
        plan_result.path,
        model.geometry,
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
            planned_at_s=t,
        ),
    )
end

function calcGuidanceEffect!(model::RPOGuidanceModel, u, p, t::Float64, sat_idx::Int)
    sat_idx == model.chaser_idx || return nothing
    if !model.plan_buffer.valid || model.force_replan
        plan = build_rpo_plan(model, u, t)
        update_rpo_plan_buffer!(model.plan_buffer, plan, t)
        model.force_replan = false
    end
    return nothing
end
