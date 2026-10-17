"""Extract chaser position and velocity for the RPO controller."""
function _rpo_control_state_pos_vel(u, idx::Int)
    return SVector{3, Float64}(u.sc[idx].pos), SVector{3, Float64}(u.sc[idx].vel)
end

"""Apply the RPO MPC controller at the current simulation time and cache the actuator command."""
function update_control!(model::RPOMPCControlModel, u, p::ODEParams, t::Float64, sat_idx::Int)
    sat_idx == model.chaser_index || return nothing
    # The engine has already called the planner that writes this model's plan
    # buffer (`before_reference_control!` for its linked producers).
    model.plan_buffer.valid || return nothing
    model.controller === nothing && return nothing

    r_chaser, v_chaser = _rpo_control_state_pos_vel(u, model.chaser_index)
    r_target, v_target = _rpo_control_state_pos_vel(u, model.target_index)
    x_rel = inertial_to_rtn_relative_state(r_chaser, v_chaser, r_target, v_target)
    t_plan = max(0.0, t - model.plan_buffer.updated_at_s)
    x_ref = rpo_ref_preview(model.plan_buffer.plan, t_plan, model.control_dt_s, model.controller.horizon)
    a_rtn = rpo_lqmpc_control(model.controller, x_rel, x_ref)
    a_ii = rtn_accel_to_inertial(a_rtn, r_target, v_target)
    mass = Float64(u.sc[model.chaser_index].mass)
    force_ii_des = mass * a_ii

    q = SVector{4, Float64}(u.sc[model.chaser_index].q)
    R_i_to_b = rot(q)
    force_body_des = R_i_to_b * force_ii_des
    thruster_forces = rpo_allocate_six_axis_thrusters(force_body_des, model.thrusters; minimum_throttle_fraction=model.minimum_throttle_fraction)
    force_body, thruster_torque_body = rpo_thruster_wrench_body(thruster_forces, model.thrusters)
    rw_torque_body = rpo_reaction_wheel_torque_command(model, u, sat_idx)
    force_ii = R_i_to_b' * force_body
    model.held = RPOHeldActuation(
        force_ii=force_ii,
        torque_body=thruster_torque_body + rw_torque_body,
        thruster_forces_n=thruster_forces,
        rw_torque_body=rw_torque_body,
    )
    model.actuation === :mounted &&
        _rpo_fill_mounted_levels!(model, p.args.dynamics_model.spacecraft[model.chaser_index], thruster_forces)
    log = model.command_log
    if log !== nothing
        push!(log.t_s, t)
        push!(log.x_rel_rtn, SVector{6, Float64}(x_rel))
        push!(log.accel_cmd_rtn, a_rtn)
        push!(log.qp_status, model.controller.qp_results.info.status)
        push!(log.force_body_desired_n, force_body_des)
        push!(log.thruster_forces_n, thruster_forces)
        push!(log.mass_kg, mass)
        push!(log.q_chaser, q)
    end
    return nothing
end

# Each run starts from the held command the model had before its first run (a
# planner's plan-buffer link resets it first), a zero warm start, a fresh solver
# workspace, an empty command log and no mounted levels (G-8).
function SimulationLifecycle.initialize_control!(model::RPOMPCControlModel, u, p, t)
    if model.run_start === nothing
        model.run_start = deepcopy(model.held)
    else
        model.held = deepcopy(model.run_start)
    end
    ctrl = model.controller
    if ctrl isa RpoLQMPCController
        fill!(ctrl.U_prev, 0.0)
        ctrl.qp_model = _rpo_lqmpc_workspace(ctrl.H, ctrl.G, ctrl.W)
        ctrl.qp_results = OSQP.Results()
    end
    log = model.command_log
    if log !== nothing
        for name in fieldnames(RPOControlCommandLog)
            empty!(getfield(log, name))
        end
    end
    empty!(model.mounted_levels)
    return nothing
end

"""Return the held force and torque command for a control effector."""
function control_wrench(model::RPOMPCControlModel, x::StateSample, p::ODEParams, t::Float64)
    x.spacecraft_index == model.chaser_index || return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    return model.held.force_ii, model.held.torque_body
end

function actuation_mode(model::RPOMPCControlModel)
    model.actuation === :ideal || model.actuation === :mounted ||
        throw(ArgumentError("RPOMPCControlModel actuation must be :ideal or :mounted, got $(repr(model.actuation))."))
    return model.actuation
end

# One level per actuator mounted on the chaser, in mounting order: the six
# allocated forces on the first six mounted thrusters, every other actuator idle.
function _rpo_fill_mounted_levels!(model::RPOMPCControlModel, chaser, thruster_forces)
    mounted_thruster_count(chaser) >= 6 || throw(ArgumentError(
        "Mounted RPO actuation needs the six-axis jets mounted first on the chaser: add six_axis_thrusters(model.thrusters)."))
    levels = model.mounted_levels
    resize!(levels, mounted_actuator_count(chaser))
    fill!(levels, 0.0)
    @inbounds for j in 1:6
        levels[j] = thruster_forces[j]
    end
    return levels
end

"""The chaser's mounted actuator levels held since the last control tick (`:mounted` actuation); none before the first tick."""
function actuator_levels(model::RPOMPCControlModel, i::Int)
    (i == model.chaser_index && !isempty(model.mounted_levels)) || return nothing
    return model.mounted_levels
end

"""Return the held reaction-wheel-attributable torque, held constant (ZOH) between control ticks like `control_wrench`."""
function control_wheel_torque(model::RPOMPCControlModel, x::StateSample, p::ODEParams, t::Float64)
    x.spacecraft_index == model.chaser_index || return nothing
    return model.held.rw_torque_body
end

"""
    control_thruster_levels(model::RPOMPCControlModel, i)

The held six-axis thruster command as firing levels: each jet's commanded
thrust over its rating. The chaser's remaining thrusters, if the vehicle
carries more than the six the allocator drives, stay idle.
"""
function control_thruster_levels(model::RPOMPCControlModel, i::Int)
    i == model.chaser_index || return nothing
    levels = zeros(6)
    @inbounds for j in 1:6
        rating = model.thrusters.max_thrust_n[j]
        rating > 0.0 && (levels[j] = clamp(model.held.thruster_forces_n[j] / rating, 0.0, 1.0))
    end
    return levels
end

"""Return the propellant the six jets use, in kg/s (positive while they fire)."""
function control_propellant_rate(model::RPOMPCControlModel, x::StateSample, p::ODEParams, t::Float64)::Float64
    x.spacecraft_index == model.chaser_index || return 0.0
    rate = 0.0
    @inbounds for j in 1:6
        thrust = model.held.thruster_forces_n[j]
        if thrust > 0.0
            rate += propellant_mass_rate(thrust, model.thrusters.isp_s[j])
        end
    end
    return rate
end
