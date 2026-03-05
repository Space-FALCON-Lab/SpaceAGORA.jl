@inline function _accumulate_dynamic_effectors!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    sc_view,
    p,
    sat_idx::Int,
    dynamic_effectors::Tuple,
    effector_decision
)
    effector_started_ns = time_ns()
    n_effectors = length(dynamic_effectors)
    if effector_decision.use_threads
        reduced = SimulationModel.ParallelPolicy.threaded_reduce(
            n_effectors,
            effector_decision.allotment,
            () -> MVector{6, Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (local_sum, eff_idx) -> begin
                effector = dynamic_effectors[eff_idx]
                force, torque = SimulationModel.calcForceTorque(effector, sc_view, p, sat_idx)
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
            force, torque = SimulationModel.calcForceTorque(effector, sc_view, p, sat_idx)
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
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, dynamic_effectors, effector_decision)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)

                du_view.pos .= sc_view.vel
                du_view.vel .= forces / sc_view.mass
                du_view.mass = mass_rate

                if p.args.mission_configuration.orientation_sim
                    ω_body = SVector{3, Float64}(sc_view.ω)
                    inertia_tensor = spacecraft[i].inertia_tensor
                    τ_body = SVector{3, Float64}(torques)
                    du_view.q .= 0.5 * quat_mult(SVector{4, Float64}(ω_body..., 0.0), sc_view.q)
                    du_view.ω .= inertia_tensor \ (τ_body - cross(ω_body, inertia_tensor * ω_body))
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, p.shared_buffers.heat_rates[i])
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
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, dynamic_effectors, effector_decision)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)

                du_view.pos .= sc_view.vel
                du_view.vel .= forces / sc_view.mass
                du_view.mass = mass_rate

                if p.args.mission_configuration.orientation_sim
                    ω_body = SVector{3, Float64}(sc_view.ω)
                    inertia_tensor = spacecraft[i].inertia_tensor
                    τ_body = SVector{3, Float64}(torques)
                    du_view.q .= 0.5 * quat_mult(SVector{4, Float64}(ω_body..., 0.0), sc_view.q)
                    du_view.ω .= inertia_tensor \ (τ_body - cross(ω_body, inertia_tensor * ω_body))
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, p.shared_buffers.heat_rates[i])
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
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, dynamic_effectors, effector_decision)

                du_view.pos .= sc_view.vel
                du_view.vel .= forces / sc_view.mass
                du_view.mass = 0.0

                if p.args.mission_configuration.orientation_sim
                    ω_body = SVector{3, Float64}(sc_view.ω)
                    inertia_tensor = spacecraft[i].inertia_tensor
                    τ_body = SVector{3, Float64}(torques)
                    du_view.q .= 0.5 * quat_mult(SVector{4, Float64}(ω_body..., 0.0), sc_view.q)
                    du_view.ω .= inertia_tensor \ (τ_body - cross(ω_body, inertia_tensor * ω_body))
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, p.shared_buffers.heat_rates[i])
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
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, dynamic_effectors, effector_decision)

                du_view.pos .= sc_view.vel
                du_view.vel .= forces / sc_view.mass
                du_view.mass = 0.0

                if p.args.mission_configuration.orientation_sim
                    ω_body = SVector{3, Float64}(sc_view.ω)
                    inertia_tensor = spacecraft[i].inertia_tensor
                    τ_body = SVector{3, Float64}(torques)
                    du_view.q .= 0.5 * quat_mult(SVector{4, Float64}(ω_body..., 0.0), sc_view.q)
                    du_view.ω .= inertia_tensor \ (τ_body - cross(ω_body, inertia_tensor * ω_body))
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, p.shared_buffers.heat_rates[i])
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

                du_view.pos .= 0.0
                du_view.vel .= forces / sc_view.mass
                du_view.mass = mass_rate

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    τ_body = SVector{3, Float64}(torques)
                    du_view.q .= 0.0
                    du_view.ω .= inertia_tensor \ τ_body
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

                du_view.pos .= 0.0
                du_view.vel .= forces / sc_view.mass
                du_view.mass = mass_rate

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    τ_body = SVector{3, Float64}(torques)
                    du_view.q .= 0.0
                    du_view.ω .= inertia_tensor \ τ_body
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
            sc_view.q .= spacecraft.initial_condition.q
            sc_view.ω .= spacecraft.initial_condition.ang_vel
        end
    end

    return state
end
