"""Return robot-arm HYPR PSO weights for an iteration."""
function _robot_arm_hypr_iteration_weights(cfg::RobotArmHYPRConfig, iter::Int)
    return hypr_iteration_weights(
        cfg.schedule_enable,
        cfg.n_iters,
        iter,
        cfg.w_inertia,
        cfg.c1,
        cfg.c2,
        cfg.schedule_transition_fraction,
        cfg.schedule_w_min,
        cfg.schedule_w_end_fraction,
        cfg.schedule_c1_end_fraction,
        cfg.schedule_c2_end_fraction,
        cfg.schedule_c_min,
        cfg.schedule_c_max,
    )
end

"""Cull weak robot-arm HYPR particles and reseed them inside joint bounds."""
function _robot_arm_hypr_cull_swarm!(
    positions,
    velocities,
    pbest,
    pbest_cost,
    lo_rep,
    hi_rep,
    gbest,
    cfg::RobotArmHYPRConfig,
    iter::Int,
    rng,
)
    if !cfg.cull_enable || iter < cfg.cull_start_iter || cfg.cull_fraction_max <= 0.0 ||
            !isfinite(minimum(pbest_cost))
        return 0
    end
    n_replace = min(size(positions, 2) - 1, floor(Int, cfg.cull_fraction_max * size(positions, 2)))
    n_replace <= 0 && return 0
    worst = sortperm(pbest_cost; rev=true)[1:n_replace]
    @inbounds for pidx in worst
        for d in axes(positions, 1)
            span = hi_rep[d] - lo_rep[d]
            jitter = cfg.cull_noise_scale * span * randn(rng)
            positions[d, pidx] = clamp(gbest[d] + jitter, lo_rep[d], hi_rep[d])
            velocities[d, pidx] = 0.0
            pbest[d, pidx] = positions[d, pidx]
        end
        pbest_cost[pidx] = Inf
    end
    return n_replace
end

"""Check whether a robot-arm HYPR cost improvement is large enough to count."""
function _robot_arm_hypr_material_improvement(new_cost::Real, reference_cost::Real, cfg::RobotArmHYPRConfig)
    return hypr_material_improvement(
        new_cost,
        reference_cost,
        cfg.early_stopping_min_abs_improvement,
        cfg.early_stopping_min_rel_improvement,
    )
end

"""Return whether robot-arm HYPR components satisfy early-stopping feasibility."""
function _robot_arm_hypr_early_stopping_feasible(components, cfg::RobotArmHYPRConfig)
    !cfg.early_stopping_require_feasible && return true
    return getproperty(components, :J_obs) <= 1.0e-9
end

"""Construct a full robot-arm control-point matrix from flattened internal points."""
function _robot_arm_control_points(q_start, q_goal, pos, n_waypoints::Int)
    q0 = Float64.(collect(q_start))
    qf = Float64.(collect(q_goal))
    n = length(q0)
    points = zeros(n, n_waypoints + 2)
    points[:, 1] .= q0
    points[:, end] .= qf
    @inbounds for j in 1:n_waypoints
        offset = n * (j - 1)
        points[:, j + 1] .= @view pos[(offset + 1):(offset + n)]
    end
    return points
end

"""Create straight-line robot-arm control points between start and goal."""
function _robot_arm_seed_control_points(q_start, q_goal, n_waypoints::Int)
    q0 = Float64.(collect(q_start))
    qf = Float64.(collect(q_goal))
    n = length(q0)
    points = zeros(n, n_waypoints + 2)
    points[:, 1] .= q0
    points[:, end] .= qf
    for j in 1:n_waypoints
        α = j / (n_waypoints + 1)
        points[:, j + 1] .= (1.0 - α) .* q0 .+ α .* qf
    end
    return points
end

"""Flatten robot-arm internal control points into a PSO particle vector."""
function _robot_arm_flatten_internal_points(points)
    n = size(points, 1)
    n_waypoints = size(points, 2) - 2
    flat = zeros(n * n_waypoints)
    for j in 1:n_waypoints
        offset = n * (j - 1)
        flat[(offset + 1):(offset + n)] .= points[:, j + 1]
    end
    return flat
end

"""Sample a robot-arm HYPR path as a Bezier curve or polyline."""
function robot_arm_sample_hypr_path(points, n_samples::Int; curve_type::Symbol=:bezier)
    return hypr_sample_count_path(points, n_samples; curve_type=curve_type)
end

"""Return joint-space length of sampled robot-arm path points."""
function _robot_arm_path_length(samples)
    return hypr_path_length(samples)
end

"""Compute squared second-difference smoothness for a robot-arm path."""
function _robot_arm_path_smoothness(samples)
    size(samples, 2) < 3 && return 0.0
    total = 0.0
    @inbounds for k in 2:(size(samples, 2) - 1)
        total += sum(abs2, samples[:, k + 1] .- 2.0 .* samples[:, k] .+ samples[:, k - 1])
    end
    return total / max(size(samples, 2) - 2, 1)
end

"""Compute local retiming scale from spacecraft reaction demand."""
function _robot_arm_hypr_reaction_scale(cfg::RobotArmHYPRConfig, demand_ratio::Real, k::Int, n_seg::Int)
    cfg.retime_reaction_time_s <= 0.0 && return 1.0
    denom = max(n_seg - 1, 1)
    phase = clamp((k - 1) / denom, 0.0, 1.0)
    taper = 4.0 * phase * (1.0 - phase)
    scale = 1.0 + cfg.retime_reaction_gain * cfg.retime_reaction_time_s * max(Float64(demand_ratio), 0.0) * taper
    return clamp(scale, cfg.retime_min_scale, cfg.retime_max_scale)
end

"""Accumulate segment time scales into a robot-arm reference time grid."""
function _robot_arm_hypr_reference_times_from_scales(nominal_dt::Float64, scales::AbstractVector{<:Real})
    t_ref = zeros(Float64, length(scales) + 1)
    @inbounds for k in eachindex(scales)
        t_ref[k + 1] = t_ref[k] + nominal_dt * Float64(scales[k])
    end
    return t_ref
end

"""Compute link center-of-mass histories for a joint reference."""
function _robot_arm_hypr_link_com_history(model::ClothArmModel, base_pose::ClothArmBasePose, q_ref::Matrix{Float64})
    n_links = length(model.links)
    nt = size(q_ref, 2)
    com = Array{Float64}(undef, 3, n_links, nt)
    @inbounds for k in 1:nt
        pose = cloth_fk(model, base_pose, q_ref[:, k])
        for i in 1:n_links
            com[:, i, k] .= pose.link_com_positions[i]
        end
    end
    return com
end

"""Estimate rigid-arm base force and torque demand ratios along a trajectory."""
function _robot_arm_hypr_rigid_base_wrench_ratios(
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    q_ref::Matrix{Float64},
    t_ref::AbstractVector{<:Real},
    cfg::RobotArmHYPRConfig,
)
    if !(isfinite(cfg.retime_max_base_force_n) || isfinite(cfg.retime_max_base_torque_nm))
        return (force_ratio=0.0, torque_ratio=0.0, node_ratio=zeros(size(q_ref, 2)))
    end
    nt = size(q_ref, 2)
    node_ratio = zeros(nt)
    force_ratio = 0.0
    torque_ratio = 0.0
    nt >= 2 || return (force_ratio=force_ratio, torque_ratio=torque_ratio, node_ratio=node_ratio)
    com = _robot_arm_hypr_link_com_history(model, base_pose, q_ref)
    masses = [link.mass_kg for link in model.links]
    @inbounds for k in 1:nt
        dt_prev = k == 1 ? max(Float64(t_ref[2] - t_ref[1]), 1.0e-9) : max(Float64(t_ref[k] - t_ref[k - 1]), 1.0e-9)
        dt_next = k == nt ? max(Float64(t_ref[end] - t_ref[end - 1]), 1.0e-9) : max(Float64(t_ref[k + 1] - t_ref[k]), 1.0e-9)
        force = SVector{3, Float64}(0.0, 0.0, 0.0)
        torque = SVector{3, Float64}(0.0, 0.0, 0.0)
        for i in eachindex(masses)
            v_prev = k == 1 ?
                SVector{3, Float64}(0.0, 0.0, 0.0) :
                (SVector{3, Float64}(com[:, i, k]) - SVector{3, Float64}(com[:, i, k - 1])) / dt_prev
            v_next = k == nt ?
                SVector{3, Float64}(0.0, 0.0, 0.0) :
                (SVector{3, Float64}(com[:, i, k + 1]) - SVector{3, Float64}(com[:, i, k])) / dt_next
            acc = 2.0 * (v_next - v_prev) / (dt_prev + dt_next)
            fi = masses[i] * acc
            ri = SVector{3, Float64}(com[:, i, k]) - base_pose.position
            force += fi
            torque += cross(ri, fi)
        end
        if isfinite(cfg.retime_max_base_force_n)
            f_ratio = cfg.retime_base_wrench_model_gain * maximum(abs.(force)) / cfg.retime_max_base_force_n
            force_ratio = max(force_ratio, f_ratio)
            node_ratio[k] = max(node_ratio[k], f_ratio)
        end
        if isfinite(cfg.retime_max_base_torque_nm)
            τ_ratio = cfg.retime_base_wrench_model_gain * maximum(abs.(torque)) / cfg.retime_max_base_torque_nm
            torque_ratio = max(torque_ratio, τ_ratio)
            node_ratio[k] = max(node_ratio[k], τ_ratio)
        end
    end
    return (force_ratio=force_ratio, torque_ratio=torque_ratio, node_ratio=node_ratio)
end

"""Extract cloth-coupled arm state vectors for reaction-load analysis."""
function _robot_arm_hypr_cloth_state_for_reaction(plan::RobotArmPlan, x, parts_fn)
    n = length(plan.model.links)
    state = (
        pos=zeros(3),
        vel=zeros(3),
        q=Float64[0.0, 0.0, 0.0, 1.0],
        ω=zeros(3),
        mass=1.0,
        heat_loads=zeros(1),
        arm_r=zeros(3, n),
        arm_q=zeros(4, n),
        arm_v=zeros(3, n),
        arm_ω=zeros(3, n),
    )
    @inbounds for i in 1:n
        part = parts_fn(x, i)
        state.arm_r[:, i] .= part.r
        state.arm_q[:, i] .= part.q
        state.arm_v[:, i] .= part.v
        state.arm_ω[:, i] .= part.ω
    end
    return state
end

"""Estimate cloth-coupled base wrench ratios along a robot-arm trajectory."""
function _robot_arm_hypr_cloth_base_wrench_ratios(
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    q_ref::Matrix{Float64},
    t_ref::AbstractVector{<:Real},
    cfg::RobotArmHYPRConfig,
)
    parent = parentmodule(@__MODULE__)
    if !isdefined(parent, :ClothRobotArmDynamics) || !isdefined(parent, :ClothMultibody)
        return _robot_arm_hypr_rigid_base_wrench_ratios(model, base_pose, q_ref, t_ref, cfg)
    end
    dynmod = getfield(parent, :ClothRobotArmDynamics)
    mbmod = getfield(parent, :ClothMultibody)
    simulate_fn = getfield(dynmod, :simulate_cloth_robot_arm_plan)
    assign_rhs_fn = getfield(dynmod, :assign_coupled_cloth_robot_arm_rhs!)
    parts_fn = getfield(mbmod, :compliant_state_parts)

    target = cloth_fk(model, base_pose, q_ref[:, end]).end_effector_position
    plan = _robot_arm_plan_from_q_reference(
        model,
        base_pose,
        q_ref[:, 1],
        q_ref[:, end],
        target,
        t_ref,
        q_ref;
        planner=:hypr,
    )
    sim = simulate_fn(
        plan;
        dt_s=cfg.retime_cloth_dt_s,
        duration_s=t_ref[end],
        integrator=:implicit_midpoint,
        k_translation_n_m=cfg.retime_cloth_k_translation_n_m,
        c_translation_n_s_m=cfg.retime_cloth_c_translation_n_s_m,
        k_rotation_n_m_rad=cfg.retime_cloth_k_rotation_n_m_rad,
        c_rotation_n_m_s_rad=cfg.retime_cloth_c_rotation_n_m_s_rad,
    )

    nt = size(q_ref, 2)
    node_ratio = zeros(nt)
    force_ratio = 0.0
    torque_ratio = 0.0
    @inbounds for (sample_idx, t) in pairs(sim.trajectory.t_s)
        state = _robot_arm_hypr_cloth_state_for_reaction(plan, sim.trajectory.states[sample_idx], parts_fn)
        du = (
            pos=zeros(3),
            vel=zeros(3),
            q=zeros(4),
            ω=zeros(3),
            mass=0.0,
            heat_loads=zeros(1),
            arm_r=zeros(3, length(plan.model.links)),
            arm_q=zeros(4, length(plan.model.links)),
            arm_v=zeros(3, length(plan.model.links)),
            arm_ω=zeros(3, length(plan.model.links)),
        )
        base_force = MVector{3, Float64}(0.0, 0.0, 0.0)
        base_torque = MVector{3, Float64}(0.0, 0.0, 0.0)
        assign_rhs_fn(
            du,
            state,
            plan,
            t,
            base_force,
            base_torque;
            k_translation_n_m=cfg.retime_cloth_k_translation_n_m,
            c_translation_n_s_m=cfg.retime_cloth_c_translation_n_s_m,
            k_rotation_n_m_rad=cfg.retime_cloth_k_rotation_n_m_rad,
            c_rotation_n_m_s_rad=cfg.retime_cloth_c_rotation_n_m_s_rad,
        )

        ratio = 0.0
        if isfinite(cfg.retime_max_base_force_n)
            f_ratio = cfg.retime_base_wrench_model_gain * maximum(abs.(base_force)) / cfg.retime_max_base_force_n
            force_ratio = max(force_ratio, f_ratio)
            ratio = max(ratio, f_ratio)
        end
        if isfinite(cfg.retime_max_base_torque_nm)
            τ_ratio = cfg.retime_base_wrench_model_gain * maximum(abs.(base_torque)) / cfg.retime_max_base_torque_nm
            torque_ratio = max(torque_ratio, τ_ratio)
            ratio = max(ratio, τ_ratio)
        end
        idx = clamp(searchsortedfirst(t_ref, t), 1, nt)
        node_ratio[idx] = max(node_ratio[idx], ratio)
        if idx > 1
            node_ratio[idx - 1] = max(node_ratio[idx - 1], ratio)
        end
    end
    return (force_ratio=force_ratio, torque_ratio=torque_ratio, node_ratio=node_ratio)
end

"""Dispatch robot-arm reaction-load estimation to rigid or cloth-coupled models."""
function _robot_arm_hypr_base_wrench_ratios(
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    q_ref::Matrix{Float64},
    t_ref::AbstractVector{<:Real},
    cfg::RobotArmHYPRConfig,
)
    if cfg.retime_cloth_physics_enable
        return _robot_arm_hypr_cloth_base_wrench_ratios(model, base_pose, q_ref, t_ref, cfg)
    end
    return _robot_arm_hypr_rigid_base_wrench_ratios(model, base_pose, q_ref, t_ref, cfg)
end

"""Retiming a robot-arm joint reference using reaction-load scaling."""
function _robot_arm_hypr_retime_reference(
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    q_ref::Matrix{Float64},
    planner_config::RobotArmPlannerConfig,
    cfg::RobotArmHYPRConfig,
)
    nt = size(q_ref, 2)
    nt >= 2 || return _reference_times(planner_config.dt_s, planner_config.duration_s)
    if !cfg.retime_enable
        return _reference_times(planner_config.dt_s, planner_config.duration_s)
    end

    n_seg = nt - 1
    nominal_dt = planner_config.duration_s / n_seg
    nominal_dt > 0.0 || throw(ArgumentError("planner duration must be positive for arm HYPR retiming."))

    scales = fill(max(cfg.retime_min_scale, 1.0), n_seg)
    vmax = cfg.retime_max_joint_velocity_rad_s
    amax = cfg.retime_max_joint_acceleration_rad_s2

    @inbounds for k in 1:n_seg
        if isfinite(vmax)
            dq = maximum(abs.(q_ref[:, k + 1] .- q_ref[:, k]))
            scales[k] = max(scales[k], dq / max(vmax * nominal_dt, 1.0e-12))
        end
    end

    if isfinite(amax) && nt >= 3
        @inbounds for k in 2:(nt - 1)
            ddq = maximum(abs.(q_ref[:, k + 1] .- 2.0 .* q_ref[:, k] .+ q_ref[:, k - 1])) / nominal_dt^2
            node_scale = sqrt(ddq / max(amax, 1.0e-12))
            scales[k - 1] = max(scales[k - 1], node_scale)
            scales[k] = max(scales[k], node_scale)
        end
    end

    if !(isfinite(cfg.retime_max_base_force_n) || isfinite(cfg.retime_max_base_torque_nm))
        @inbounds for k in 1:n_seg
            demand_ratio = max(scales[k] - 1.0, 0.0)
            scales[k] = max(scales[k], _robot_arm_hypr_reaction_scale(cfg, demand_ratio, k, n_seg))
            scales[k] = clamp(scales[k], cfg.retime_min_scale, cfg.retime_max_scale)
        end
    else
        @inbounds for k in 1:n_seg
            scales[k] = clamp(scales[k], cfg.retime_min_scale, cfg.retime_max_scale)
        end
    end

    if isfinite(cfg.retime_max_base_force_n) || isfinite(cfg.retime_max_base_torque_nm)
        for _ in 1:cfg.retime_base_wrench_iters
            t_ref_iter = _robot_arm_hypr_reference_times_from_scales(nominal_dt, scales)
            ratios = _robot_arm_hypr_base_wrench_ratios(model, base_pose, q_ref, t_ref_iter, cfg)
            max_demand = maximum(ratios.node_ratio)
            max_demand <= 1.0 && break
            changed = false
            global_scale = cfg.retime_base_wrench_margin * sqrt(max_demand)
            if isfinite(global_scale) && global_scale > 1.0
                @inbounds for k in eachindex(scales)
                    new_scale = clamp(scales[k] * global_scale, cfg.retime_min_scale, cfg.retime_max_scale)
                    changed |= new_scale > scales[k] + 1.0e-9
                    scales[k] = new_scale
                end
            end
            @inbounds for k in 1:nt
                demand = ratios.node_ratio[k]
                demand <= 1.0 && continue
                s = cfg.retime_base_wrench_margin * sqrt(demand)
                if k > 1
                    new_scale = clamp(scales[k - 1] * s, cfg.retime_min_scale, cfg.retime_max_scale)
                    changed |= new_scale > scales[k - 1] + 1.0e-9
                    scales[k - 1] = new_scale
                end
                if k <= n_seg
                    new_scale = clamp(scales[k] * s, cfg.retime_min_scale, cfg.retime_max_scale)
                    changed |= new_scale > scales[k] + 1.0e-9
                    scales[k] = new_scale
                end
            end
            changed || break
        end
    end
    return _robot_arm_hypr_reference_times_from_scales(nominal_dt, scales)
end

"""Return the distance from a point to a 3D segment."""
function _robot_arm_segment_distance(p::SVector{3, Float64}, a::SVector{3, Float64}, b::SVector{3, Float64})
    ab = b - a
    denom = dot(ab, ab)
    denom <= eps(Float64) && return norm(p - a)
    t = clamp(dot(p - a, ab) / denom, 0.0, 1.0)
    return norm(p - (a + t * ab))
end
