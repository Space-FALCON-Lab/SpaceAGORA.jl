@inline function _evaluate_dynamic_effector(
    effector,
    sc_view,
    state_sample,
    p,
    sat_idx::Int,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    if _wrench_method_available(effector)
        state_sample === nothing && error("wrench-based effector evaluation requires a typed StateSample.")
        req = SimulationModel.environment_requirements(effector)
        env = sample_environment(req, effector, sc_view, p, sat_idx, t; write_buffers=false)
        return SimulationModel.wrench(effector, state_sample, env, t)
    end
    return SimulationModel.calcForceTorque(effector, sc_view, p, sat_idx)
end

@inline function _solver_partition_validated(effector)::Symbol
    partition = SimulationModel.solver_partition(effector)
    if partition === :implicit || partition === :explicit
        return partition
    end
    throw(ArgumentError(
        "solver_partition($(nameof(typeof(effector)))) must return :implicit or :explicit, got $(repr(partition))."
    ))
end

@inline _effector_in_partition(effector, partition::Symbol)::Bool =
    _solver_partition_validated(effector) === partition

@inline function _partition_needs_state_sample(
    dynamic_effectors::Tuple,
    partition::Symbol,
)::Bool
    return any(_effector_in_partition(effector, partition) && _wrench_method_available(effector) for effector in dynamic_effectors)
end

@inline function _partition_selected_count(
    dynamic_effectors::Tuple,
    partition::Symbol,
)::Int
    count = 0
    @inbounds for effector in dynamic_effectors
        count += _effector_in_partition(effector, partition) ? 1 : 0
    end
    return count
end

@inline function _accumulate_dynamic_effectors!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    sc_view,
    p,
    sat_idx::Int,
    t::Float64,
    dynamic_effectors::Tuple,
    effector_decision
)
    # Only pay the time_ns() syscall overhead when telemetry is actually needed.
    needs_timing = effector_decision.policy_applied
    effector_started_ns = needs_timing ? time_ns() : UInt64(0)
    n_effectors = length(dynamic_effectors)
    needs_state_sample = any(_wrench_method_available(effector) for effector in dynamic_effectors)
    state_sample = if needs_state_sample
        spacecraft = p.args.dynamics_model.spacecraft[sat_idx]
        build_state_sample(sc_view, spacecraft, p.args.mission_configuration.orientation_sim)
    else
        nothing
    end
    if effector_decision.use_threads
        reduced = SimulationModel.ParallelPolicy.threaded_reduce(
            n_effectors,
            effector_decision.allotment,
            () -> MVector{6, Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (local_sum, eff_idx) -> begin
                effector = dynamic_effectors[eff_idx]
                force, torque = _evaluate_dynamic_effector(effector, sc_view, state_sample, p, sat_idx, t)
                local_sum[1] += force[1]
                local_sum[2] += force[2]
                local_sum[3] += force[3]
                local_sum[4] += torque[1]
                local_sum[5] += torque[2]
                local_sum[6] += torque[3]
                return nothing
            end,
            (dest, src) -> begin
                @inbounds for i in 1:6
                    dest[i] += src[i]
                end
                return nothing
            end
        )
        # Unpack reduced MVector directly into forces/torques without an intermediate SVector.
        @inbounds begin
            forces[1] = reduced[1]
            forces[2] = reduced[2]
            forces[3] = reduced[3]
            torques[1] = reduced[4]
            torques[2] = reduced[5]
            torques[3] = reduced[6]
        end
    else
        @inbounds for effector in dynamic_effectors
            force, torque = _evaluate_dynamic_effector(effector, sc_view, state_sample, p, sat_idx, t)
            forces .+= force
            torques .+= torque
        end
    end
    if needs_timing
        elapsed_ns = Int64(time_ns() - effector_started_ns)
        if sat_idx == 1
            _update_effector_cost_model!(p.shared_buffers, n_effectors, elapsed_ns, effector_decision.allotment)
        end
        SimulationModel.ParallelPolicy.record_policy_observation!(
            :dynamic_effectors;
            mode=effector_decision.mode,
            num_items=n_effectors,
            use_threads=effector_decision.use_threads,
            elapsed_ns=elapsed_ns
        )
    end
    return nothing
end

@inline function _accumulate_dynamic_effectors_partitioned!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    sc_view,
    p,
    sat_idx::Int,
    t::Float64,
    dynamic_effectors::Tuple,
    effector_decision,
    partition::Symbol,
)
    selected_count = _partition_selected_count(dynamic_effectors, partition)
    selected_count == 0 && return nothing

    state_sample = if _partition_needs_state_sample(dynamic_effectors, partition)
        spacecraft = p.args.dynamics_model.spacecraft[sat_idx]
        build_state_sample(sc_view, spacecraft, p.args.mission_configuration.orientation_sim)
    else
        nothing
    end

    if effector_decision.use_threads && selected_count > 1
        reduced = SimulationModel.ParallelPolicy.threaded_reduce(
            length(dynamic_effectors),
            effector_decision.allotment,
            () -> MVector{6, Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (local_sum, eff_idx) -> begin
                effector = dynamic_effectors[eff_idx]
                _effector_in_partition(effector, partition) || return nothing
                force, torque = _evaluate_dynamic_effector(effector, sc_view, state_sample, p, sat_idx, t)
                local_sum[1] += force[1]
                local_sum[2] += force[2]
                local_sum[3] += force[3]
                local_sum[4] += torque[1]
                local_sum[5] += torque[2]
                local_sum[6] += torque[3]
                return nothing
            end,
            (dest, src) -> begin
                @inbounds for i in 1:6
                    dest[i] += src[i]
                end
                return nothing
            end
        )
        forces .= SVector{3, Float64}(reduced[1], reduced[2], reduced[3])
        torques .= SVector{3, Float64}(reduced[4], reduced[5], reduced[6])
        return nothing
    end

    @inbounds for effector in dynamic_effectors
        _effector_in_partition(effector, partition) || continue
        force, torque = _evaluate_dynamic_effector(effector, sc_view, state_sample, p, sat_idx, t)
        forces .+= force
        torques .+= torque
    end
    return nothing
end

@inline function _accumulate_dynamic_effectors!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    sc_view,
    p,
    sat_idx::Int,
    dynamic_effectors::Tuple,
    effector_decision
)
    t = if p === nothing
        0.0
    else
        Float64(p.shared_buffers.current_time[])
    end
    return _accumulate_dynamic_effectors!(forces, torques, sc_view, p, sat_idx, t, dynamic_effectors, effector_decision)
end

@inline function _accumulate_control_effectors!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    sc_view,
    p,
    sat_idx::Int,
    t::Float64,
    debug_control::Bool
)::Float64
    mass_rate = 0.0
    @inbounds for control_effector in p.args.control_model.control_effectors
        control_force, control_torque = SimulationModel.calcControlForceTorque(control_effector, sc_view, p, sat_idx, t)
        control_mass_rate = SimulationModel.calcControlMassFlowRate(control_effector, sc_view, p, sat_idx, t)
        if debug_control && (norm(control_force) > 0.0 || norm(control_torque) > 0.0)
            println("Applying control effect for spacecraft $sat_idx at time $t seconds:")
            println("  Control force: $control_force")
        end
        forces .+= control_force
        torques .+= control_torque
        mass_rate += isfinite(control_mass_rate) ? control_mass_rate : 0.0
    end
    return mass_rate
end

@inline function _assign_heat_rate_derivative!(du_heat::AbstractVector, heat_rates::AbstractVector)
    if length(heat_rates) == length(du_heat)
        du_heat .= heat_rates
        return nothing
    end
    du_heat .= 0.0
    n_copy = min(length(heat_rates), length(du_heat))
    @inbounds for j in 1:n_copy
        du_heat[j] = heat_rates[j]
    end
    return nothing
end

@inline function _gravity_backbone_xyz_chunk(state, sat_idx::Int)::SVector{3, Float64}
    if hasproperty(state, :sc)
        sc_state = state.sc[sat_idx]
        if hasproperty(sc_state, :pos)
            return SVector{3, Float64}(sc_state.pos)
        elseif hasproperty(sc_state, :vel)
            return SVector{3, Float64}(sc_state.vel)
        end
        return SVector{3, Float64}(sc_state)
    end

    if sat_idx <= length(state)
        sat_state = state[sat_idx]
        if hasproperty(sat_state, :pos)
            return SVector{3, Float64}(sat_state.pos)
        elseif hasproperty(sat_state, :vel)
            return SVector{3, Float64}(sat_state.vel)
        end
    end

    base = 3 * (sat_idx - 1)
    return SVector{3, Float64}(state[base + 1], state[base + 2], state[base + 3])
end

@inline function _gravity_backbone_state_sample(q_state, dq_state, p, sat_idx::Int)::StateSample
    spacecraft = p.args.dynamics_model.spacecraft[sat_idx]
    pos_ii = _gravity_backbone_xyz_chunk(q_state, sat_idx)
    vel_ii = _gravity_backbone_xyz_chunk(dq_state, sat_idx)
    mass_kg = spacecraft.dry_mass + spacecraft.prop_mass
    return StateSample(pos_ii, vel_ii, mass_kg; spacecraft=spacecraft)
end

@inline function _gravity_backbone_state_sample(u_state, p, sat_idx::Int)::StateSample
    return _gravity_backbone_state_sample(
        _gravity_backbone_position_state(u_state),
        _gravity_backbone_velocity_state(u_state),
        p,
        sat_idx,
    )
end

@inline function _gravity_backbone_core_acceleration(
    dynamic_effectors::Tuple,
    state_sample::StateSample,
    p,
    sat_idx::Int,
    t::Float64,
)::SVector{3, Float64}
    accel_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for effector in dynamic_effectors
        _gravity_backbone_structure_validated(effector) === :position_only_static_gravity || continue
        req = SimulationModel.environment_requirements(effector)
        env = sample_environment(req, effector, state_sample, p, sat_idx, t; write_buffers=false)
        accel_ii .+= SimulationModel.gravity_backbone_acceleration_ii(effector, state_sample, env, t)
    end
    return SVector{3, Float64}(accel_ii)
end

@inline function _gravity_backbone_kick_acceleration(
    dynamic_effectors::Tuple,
    state_sample::StateSample,
    p,
    sat_idx::Int,
    t::Float64,
)::SVector{3, Float64}
    accel_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for effector in dynamic_effectors
        _gravity_backbone_kick_structure_validated(effector) === :velocity_kick_explicit || continue
        req = SimulationModel.environment_requirements(effector)
        env = sample_environment(req, effector, state_sample, p, sat_idx, t; write_buffers=false)
        accel_ii .+= SimulationModel.gravity_backbone_kick_acceleration_ii(effector, state_sample, env, t)
    end
    return SVector{3, Float64}(accel_ii)
end

function _gravity_backbone_half_kick!(u_state, p, t::Float64, half_dt::Float64)
    half_dt > 0.0 || return nothing
    dynamic_effectors = p.args.dynamics_model.dynamic_effectors
    vel_state = _gravity_backbone_velocity_state(u_state)
    p.shared_buffers.current_time[] = t
    use_rhs_batch = _rhs_batch_parallel_enabled(length(vel_state.sc))

    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(vel_state.sc) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(vel_state.sc)
            if !p.is_active[i]
                continue
            end
            state_sample = _gravity_backbone_state_sample(u_state, p, i)
            accel_ii = _gravity_backbone_kick_acceleration(dynamic_effectors, state_sample, p, i, t)
            vel_state.sc[i].vel .+= half_dt .* accel_ii
        end
    else
        @inbounds for i in eachindex(vel_state.sc)
            if !p.is_active[i]
                continue
            end
            state_sample = _gravity_backbone_state_sample(u_state, p, i)
            accel_ii = _gravity_backbone_kick_acceleration(dynamic_effectors, state_sample, p, i, t)
            vel_state.sc[i].vel .+= half_dt .* accel_ii
        end
    end
    return nothing
end

function spacecraft_dynamics_gravity_backbone!(ddu, dq, q, p, t::Float64)
    q_state = q.sc
    dq_state = dq.sc
    ddu_state = ddu.sc
    dynamic_effectors = p.args.dynamics_model.dynamic_effectors
    p.shared_buffers.current_time[] = t
    use_rhs_batch = _rhs_batch_parallel_enabled(length(q_state))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(q_state) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(q_state)
            if !p.is_active[i]
                ddu_state[i].vel .= 0.0
                continue
            end
            state_sample = _gravity_backbone_state_sample(q, dq, p, i)
            ddu_state[i].vel .= _gravity_backbone_core_acceleration(dynamic_effectors, state_sample, p, i, t)
        end
    else
        @inbounds for i in eachindex(q_state)
            if !p.is_active[i]
                ddu_state[i].vel .= 0.0
                continue
            end
            state_sample = _gravity_backbone_state_sample(q, dq, p, i)
            ddu_state[i].vel .= _gravity_backbone_core_acceleration(dynamic_effectors, state_sample, p, i, t)
        end
    end
    return nothing
end

@inline function _assign_orientation_rhs!(
    du_view,
    sc_view,
    inertia_tensor::AbstractMatrix{<:Real},
    torques::AbstractVector{<:Real};
    propagate_quaternion::Bool,
    include_gyroscopic::Bool,
)
    omega_body = SimulationModel.DynamicsRotational.body_angular_velocity(sc_view.ω)
    tau_body = SimulationModel.DynamicsRotational.body_torque(torques)

    if propagate_quaternion
        du_view.q .= SimulationModel.DynamicsRotational.quaternion_derivative(omega_body, sc_view.q)
    else
        du_view.q .= 0.0
    end

    du_view.ω .= SimulationModel.DynamicsRotational.angular_acceleration(
        omega_body,
        inertia_tensor,
        tau_body;
        include_gyroscopic=include_gyroscopic,
    )
    return nothing
end

"""Return whether an effector owns a robot-arm plan for the requested spacecraft."""
@inline function _robot_arm_effector_matches(effector, sat_idx::Int)::Bool
    return hasproperty(effector, :plan) &&
        hasproperty(effector, :spacecraft_idx) &&
        getproperty(effector, :spacecraft_idx) == sat_idx &&
        getproperty(effector, :plan) isa SimulationModel.RobotArmPlan
end

"""Read a numeric effector field as `Float64`, falling back to a default value."""
@inline _effector_float(effector, name::Symbol, default::Float64)::Float64 =
    hasproperty(effector, name) ? Float64(getproperty(effector, name)) : default

"""Read an optional effector field, falling back to a default value."""
@inline _effector_value(effector, name::Symbol, default) =
    hasproperty(effector, name) ? getproperty(effector, name) : default

"""Build the cloth-coupled robot-arm RHS configuration from an active effector."""
@inline function _robot_arm_coupling_from_effector(effector, t::Float64)
    return (
        plan=getproperty(effector, :plan),
        t_s=max(0.0, t - _effector_float(effector, :updated_at_s, 0.0)),
        k_translation_n_m=_effector_value(effector, :k_translation_n_m, 5.0e3),
        c_translation_n_s_m=_effector_value(effector, :c_translation_n_s_m, 30.0),
        k_rotation_n_m_rad=_effector_value(effector, :k_rotation_n_m_rad, 15.0),
        c_rotation_n_m_s_rad=_effector_value(effector, :c_rotation_n_m_s_rad, 0.5),
        joint_actuators=hasproperty(effector, :joint_actuators) ?
            getproperty(effector, :joint_actuators) :
            SimulationModel.CompliantJointActuator[],
    )
end

"""Find the active robot-arm coupling configuration for one spacecraft."""
function _robot_arm_coupling(args, sat_idx::Int, t::Float64)
    if hasproperty(args, :control_model) && hasproperty(args.control_model, :control_effectors)
        @inbounds for effector in args.control_model.control_effectors
            _robot_arm_effector_matches(effector, sat_idx) && return _robot_arm_coupling_from_effector(effector, t)
        end
    end
    if hasproperty(args, :dynamics_model) && hasproperty(args.dynamics_model, :dynamic_effectors)
        @inbounds for effector in args.dynamics_model.dynamic_effectors
            _robot_arm_effector_matches(effector, sat_idx) && return _robot_arm_coupling_from_effector(effector, t)
        end
    end
    return nothing
end

"""Apply coupled cloth robot-arm state derivatives to one spacecraft RHS view."""
@inline function _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, sat_idx::Int, t::Float64, forces, torques)
    coupling = _robot_arm_coupling(p.args, sat_idx, t)
    coupling === nothing && return nothing
    SimulationModel.assign_coupled_cloth_robot_arm_rhs!(
        du_view,
        sc_view,
        coupling.plan,
        coupling.t_s,
        forces,
        torques;
        k_translation_n_m=coupling.k_translation_n_m,
        c_translation_n_s_m=coupling.c_translation_n_s_m,
        k_rotation_n_m_rad=coupling.k_rotation_n_m_rad,
        c_rotation_n_m_s_rad=coupling.c_rotation_n_m_s_rad,
        joint_actuators=coupling.joint_actuators,
    )
    return nothing
end

function spacecraft_dynamics!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    debug_control = p.shared_buffers.debug_control[]
    p.shared_buffers.current_time[] = t
    effector_decision = _dynamic_effector_thread_decision(p.args, p, dynamic_effectors, length(spacecraft))
    use_rhs_batch = _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=false,
                )

                SimulationModel.DynamicsTranslational.assign_full_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=false,
                )

                SimulationModel.DynamicsTranslational.assign_full_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    end
end # function spacecraft_dynamics!

function spacecraft_dynamics_slow!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    p.shared_buffers.current_time[] = t
    effector_decision = _dynamic_effector_thread_decision(p.args, p, dynamic_effectors, length(spacecraft))
    use_rhs_batch = _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=false,
                )

                SimulationModel.DynamicsTranslational.assign_slow_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=false,
                )

                SimulationModel.DynamicsTranslational.assign_slow_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    end
end # function spacecraft_dynamics_slow!

function spacecraft_dynamics_implicit_atmosphere!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    p.shared_buffers.current_time[] = t
    effector_decision = _dynamic_effector_thread_decision(p.args, p, dynamic_effectors, length(spacecraft))
    use_rhs_batch = _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors_partitioned!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision, :implicit)

                SimulationModel.DynamicsTranslational.assign_force_only_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=false,
                        include_gyroscopic=false,
                    )
                end

                du_view.heat_loads .= 0.0
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors_partitioned!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision, :implicit)

                SimulationModel.DynamicsTranslational.assign_force_only_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=false,
                        include_gyroscopic=false,
                    )
                end

                du_view.heat_loads .= 0.0
            end
        end
    end
end # function spacecraft_dynamics_implicit_atmosphere!

function spacecraft_dynamics_explicit_remainder!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    debug_control = p.shared_buffers.debug_control[]
    p.shared_buffers.current_time[] = t
    effector_decision = _dynamic_effector_thread_decision(p.args, p, dynamic_effectors, length(spacecraft))
    use_rhs_batch = _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors_partitioned!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision, :explicit)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=false,
                )

                SimulationModel.DynamicsTranslational.assign_full_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors_partitioned!(forces, torques, sc_view, p, i, t, dynamic_effectors, effector_decision, :explicit)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)
                heat_rates = SimulationModel.SimulationCallbacks._compute_stage_heat_rates!(
                    p,
                    sc_view,
                    i,
                    t;
                    use_buffered_density=false,
                )

                SimulationModel.DynamicsTranslational.assign_full_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=true,
                        include_gyroscopic=true,
                    )
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, heat_rates)
            end
        end
    end
end # function spacecraft_dynamics_explicit_remainder!

function spacecraft_dynamics_fast_control!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    spacecraft = p.args.dynamics_model.spacecraft
    debug_control = p.shared_buffers.debug_control[]
    p.shared_buffers.current_time[] = t
    use_rhs_batch = _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)

                SimulationModel.DynamicsTranslational.assign_control_only_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=false,
                        include_gyroscopic=false,
                    )
                end

                du_view.heat_loads .= 0.0
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)
                _apply_coupled_robot_arm_rhs!(du_view, sc_view, p, i, t, forces, torques)

                SimulationModel.DynamicsTranslational.assign_control_only_translational_rhs!(
                    du_view,
                    sc_view,
                    forces,
                    mass_rate,
                )

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    _assign_orientation_rhs!(
                        du_view,
                        sc_view,
                        inertia_tensor,
                        torques;
                        propagate_quaternion=false,
                        include_gyroscopic=false,
                    )
                end

                du_view.heat_loads .= 0.0
            end
        end
    end
end # function spacecraft_dynamics_fast_control!

function build_initial_conditions(args)::ComponentVector
    # Each spacecraft can have a different body count, so the state layout has
    # to be sized per spacecraft before we pack the global ComponentVector.
    sc_shapes = map(eachindex(args.dynamics_model.spacecraft)) do i
        sc = args.dynamics_model.spacecraft[i]
        n_bodies = length(sc.links)
        mass = sc.dry_mass + sc.prop_mass
        base_shape = if args.mission_configuration.orientation_sim
            (
                pos = zeros(3), 
                vel = zeros(3), 
                mass = mass, 
                heat_loads = zeros(n_bodies),
                q = Float64[0.0, 0.0, 0.0, 1.0], 
                ω = zeros(3)
            )
        else
            (
                pos = zeros(3), 
                vel = zeros(3), 
                mass = mass, 
                heat_loads = zeros(n_bodies)
            )
        end
        coupling = _robot_arm_coupling(args, i, 0.0)
        coupling === nothing && return base_shape
        return merge(base_shape, SimulationModel.coupled_cloth_robot_arm_state_shape(coupling.plan))
    end

    state = ComponentVector(sc = sc_shapes)

    for i in eachindex(args.dynamics_model.spacecraft)
        spacecraft = args.dynamics_model.spacecraft[i]
        sc_view = state.sc[i]
        ic = spacecraft.initial_condition
        r0, v0 = if ic isa CartesianInitialCondition
            collect(ic.pos), collect(ic.vel)   # assumed J2000 frame
        else
            # Orbit-element initial conditions are propagated directly in J2000.
            r_j2000, v_j2000 = orbitalelemtorv(ic, args.environment_model.planet)
            collect(r_j2000), collect(v_j2000)
        end
        sc_view.pos .= r0
        sc_view.vel .= v0
        sc_view.heat_loads .= 0.0  
        
        if args.mission_configuration.orientation_sim
            sc_view.q .= SimulationModel.project_unit_quaternion(spacecraft.initial_condition.q)
            sc_view.ω .= spacecraft.initial_condition.ang_vel
        end
        coupling = _robot_arm_coupling(args, i, 0.0)
        if coupling !== nothing
            SimulationModel.initialize_coupled_cloth_robot_arm_state!(sc_view, coupling.plan; t_s=coupling.t_s)
        end
    end

    return state
end
