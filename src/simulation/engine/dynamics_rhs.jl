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
    effector_started_ns = time_ns()
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
        forces .= SVector{3, Float64}(reduced[1], reduced[2], reduced[3])
        torques .= SVector{3, Float64}(reduced[4], reduced[5], reduced[6])
    else
        @inbounds for effector in dynamic_effectors
            force, torque = _evaluate_dynamic_effector(effector, sc_view, state_sample, p, sat_idx, t)
            forces .+= force
            torques .+= torque
        end
    end
    elapsed_ns = Int64(time_ns() - effector_started_ns)
    if effector_decision.policy_applied && sat_idx == 1
        _update_effector_cost_model!(p.shared_buffers, n_effectors, elapsed_ns, effector_decision.allotment)
    end
    if effector_decision.policy_applied
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
    # 1. Build the structure (Axis) based on each spacecraft's unique body count
    # This identifies exactly how many heat_load slots each SC needs
    sc_shapes = map(args.dynamics_model.spacecraft) do sc
        # Get the number of bodies for this specific spacecraft
        n_bodies = length(sc.links)
        mass = sc.dry_mass + sc.prop_mass
        # Create the initial state for this spacecraft with the correct size for heat_loads
        if args.mission_configuration.orientation_sim
            return (
                pos = zeros(3), 
                vel = zeros(3), 
                mass = mass, 
                heat_loads = zeros(n_bodies), # Variable size!
                q = Float64[0.0, 0.0, 0.0, 1.0], 
                ω = zeros(3)
            )
        else
            return (
                pos = zeros(3), 
                vel = zeros(3), 
                mass = mass, 
                heat_loads = zeros(n_bodies)
            )
        end
    end

    # 2. Pack everything into one ComponentVector
    # Julia allocates ONE flat array and calculates all offsets automatically
    state = ComponentVector(sc = sc_shapes) # Add more components here as needed in the future (e.g., separate orientation state if not using quaternions, etc.)

    # 3. Fill the values (The logic remains the same)
    for i in eachindex(args.dynamics_model.spacecraft)
        spacecraft = args.dynamics_model.spacecraft[i]
        sc_view = state.sc[i]
        r0, v0 = orbitalelemtorv(spacecraft.initial_condition, args.environment_model.planet)
        sc_view.pos .= r0
        sc_view.vel .= v0
        # sc_view.mass .= spacecraft.dry_mass + spacecraft.prop_mass
        # Note: heat_loads is already the correct size for this specific i!
        sc_view.heat_loads .= 0.0  
        
        if args.mission_configuration.orientation_sim
            sc_view.q .= SimulationModel.project_unit_quaternion(spacecraft.initial_condition.q)
            sc_view.ω .= spacecraft.initial_condition.ang_vel
        end
    end

    return state
end
