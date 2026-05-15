function _rpo_control_state_pos_vel(u, idx::Int)
    return SVector{3, Float64}(u.sc[idx].pos), SVector{3, Float64}(u.sc[idx].vel)
end

function calcControlEffect!(model::RPOMPCControlModel, u, p::ODEParams, t::Float64, sat_idx::Int)
    sat_idx == model.chaser_idx || return nothing
    model.plan_buffer.valid || return nothing
    model.controller === nothing && return nothing

    r_chaser, v_chaser = _rpo_control_state_pos_vel(u, model.chaser_idx)
    r_target, v_target = _rpo_control_state_pos_vel(u, model.target_idx)
    x_rel = inertial_to_rtn_relative_state(r_chaser, v_chaser, r_target, v_target)
    t_plan = max(0.0, t - model.plan_buffer.updated_at_s)
    x_ref = rpo_ref_preview(model.plan_buffer.plan, t_plan, model.control_dt_s, model.controller.horizon)
    a_rtn = rpo_lqmpc_control(model.controller, x_rel, x_ref)
    a_ii = rtn_accel_to_inertial(a_rtn, r_target, v_target)
    mass = Float64(u.sc[model.chaser_idx].mass)
    force_ii_des = mass * a_ii

    q = SVector{4, Float64}(u.sc[model.chaser_idx].q)
    R_i_to_b = rot(q)
    force_body_des = R_i_to_b * force_ii_des
    thruster_forces = rpo_allocate_six_axis_thrusters(force_body_des, model.thrusters)
    force_body, thruster_torque_body = rpo_thruster_wrench_body(thruster_forces, model.thrusters)
    rw_torque_body = rpo_reaction_wheel_torque_command(model, u, sat_idx)
    force_ii = R_i_to_b' * force_body
    model.held = RPOHeldActuation(
        force_ii=force_ii,
        torque_body=thruster_torque_body + rw_torque_body,
        thruster_forces_n=thruster_forces,
    )
    return nothing
end

function calcControlForceTorque(model::RPOMPCControlModel, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)
    i == model.chaser_idx || return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    return model.held.force_ii, model.held.torque_body
end

function calcControlMassFlowRate(model::RPOMPCControlModel, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)::Float64
    i == model.chaser_idx || return 0.0
    mdot = 0.0
    @inbounds for j in 1:6
        thrust = model.held.thruster_forces_n[j]
        if thrust > 0.0
            mdot -= thrust / (model.thrusters.isp_s[j] * 9.80665)
        end
    end
    return mdot
end
