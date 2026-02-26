include("../utils/Reference_system.jl")
using LinearAlgebra
using StaticArrays
using LoopVectorization
using ComponentArrays
using DifferentialEquations
using CSV
using DataFrames
using Polyester

const _normalize_warning_emitted = Ref(false)

@inline _typed_normalize_warning_enabled() = get(ENV, "SPACEAGORA_WARN_NORMALIZE", "1") == "1"

function _warn_legacy_normalize_flag!(args)
    if !args.simulation_settings.normalize || !_typed_normalize_warning_enabled()
        return nothing
    end
    if _normalize_warning_emitted[]
        return nothing
    end
    _normalize_warning_emitted[] = true
    @warn "SimulationSettings.normalize=true is legacy-only in typed run_simulation; propagation is always SI-native (m, s, kg). Set normalize=false to silence this warning."
    return nothing
end

function _debug_print_nan_parameter_paths!(x, path::AbstractString="p")
    if x isa Number
        if isnan(x)
            println("NaN found in parameter: $path")
        end
        return nothing
    end

    if x isa Base.RefValue{<:Number}
        xv = x[]
        if isnan(xv)
            println("NaN found in parameter: $path[]")
        end
        return nothing
    end

    if x isa AbstractArray{<:Number}
        for (idx, xv) in pairs(x)
            if isnan(xv)
                println("NaN found in parameter: $path[$idx]")
            end
        end
        return nothing
    end

    # Skip generic arrays of non-numeric types to keep debug scans bounded.
    if x isa AbstractArray
        return nothing
    end

    T = typeof(x)
    if isstructtype(T)
        for field in fieldnames(T)
            val = getfield(x, field)
            _debug_print_nan_parameter_paths!(val, string(path, ".", field))
        end
    end
    return nothing
end

function run_simulation(args; isolate_state::Bool=true, return_solution::Bool=false)
    # Isolate mutable campaign/model state by default so repeated/concurrent runs
    # do not alias shared in-memory objects.
    args = isolate_state ? deepcopy(args) : args

    # Typed pipeline is SI-native (meters, seconds, kilograms). The
    # `simulation_settings.normalize` field is retained for legacy compatibility.
    _warn_legacy_normalize_flag!(args)

    # Set up the model and initial conditions
    initial_conditions = build_initial_conditions(args)
    if args.simulation_settings.verbose
        println("Initial conditions:")
        println(initial_conditions)
    end

    # Define the ODE problem
    p = SimulationModel.ODEParams{length(args.dynamics_model.spacecraft)}(args=args) # Define the parameters for the ODE problem, including the shared buffers for the callbacks
    p.shared_buffers.debug_control[] = get(ENV, "SPACEAGORA_DEBUG_CONTROL", "0") == "1"
    p.shared_buffers.debug_initial_derivative[] = get(ENV, "SPACEAGORA_DEBUG_INITIAL_DERIVATIVE", "0") == "1"
    callbacks = SimulationModel.get_callbacks(length(args.dynamics_model.spacecraft), args.dynamics_model.dynamic_effectors, args) # Get the callbacks based on the number of satellites and the dynamic effectors being used in the simulation
    initial_time = args.initial_time
    start_epoch = from_utc(DateTime(
            initial_time.year,
            initial_time.month,
            initial_time.day,
            initial_time.hour,
            initial_time.minute,
            initial_time.second
        ))
    et_start = lock(SimulationModel.SPICE_LOCK) do
        utc2et(to_utc(start_epoch))
    end
    p.shared_buffers.et_start[] = et_start
    lock(SimulationModel.SPICE_LOCK) do
        args.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_$(args.environment_model.planet.name)", et_start)) * args.environment_model.planet.J2000_to_pci' # Initialize the planet frame at the start of the simulation (will be updated in the callback)
    end
    # println("Initial conditions:")
    # println(initial_conditions)
    # println("ODE parameters:")
    # println(p)
    # println("args.mission_configuration.mission_time: $(args.mission_configuration.mission_time)")
    prob = ODEProblem(spacecraft_dynamics!, initial_conditions, (0.0, args.mission_configuration.mission_time), p, callback=callbacks)
    if p.shared_buffers.debug_initial_derivative[]
        # 1. Manually evaluate the derivative at the start
        du_test = copy(prob.u0)
        try
            prob.f(du_test, prob.u0, prob.p, prob.tspan[1])
        catch e
            @error "The derivative function itself crashed!" exception=e
            throw(ErrorException("Initial derivative evaluation failed; aborting solve in debug mode."))
        end

        # 2. Check for NaNs and print exactly where they are
        if any(isnan, du_test)
            println("--- INITIAL NaN DETECTED ---")

            # Check global parameters in p
            _debug_print_nan_parameter_paths!(prob.p)

            # Check the state vector (u)
            # Assuming your u has a .sc field for satellites
            for (i, sat) in enumerate(du_test.sc)
                if any(isnan, sat.pos) || any(isnan, sat.vel)
                    println("NaN found in Satellite $i derivative!")
                    println("  Pos: $(sat.pos)")
                    println("  Vel: $(sat.vel)")
                end
            end

            throw(ErrorException("Initial derivative contains NaN; aborting solve in debug mode."))
        end
    end
    sol = solve(prob, Tsit5(); reltol=args.integration_tolerances.reltol_orbit, abstol=args.integration_tolerances.abstol_orbit, dtmax=args.integration_tolerances.dt_max_orbit)
    flat_sol = Array(sol) # Convert to flat array for easier processing (optional, depends on how you want to handle the results)
    indices = CartesianIndices(initial_conditions)
    # println("Full solution:")
    # println(sol.u)
    # Process and save results
    # For now, just write to csv file, but eventually want to save in a more efficient format and generate plots
    
    if args.simulation_settings.results
        results_df = DataFrame(time=sol.t)
        for i in eachindex(args.dynamics_model.spacecraft)
            results_df[!, "sc$(i)_pos_1"] = [sol.u[t].sc[i].pos[1] for t in 1:length(sol.t)]
            results_df[!, "sc$(i)_pos_2"] = [sol.u[t].sc[i].pos[2] for t in 1:length(sol.t)]
            results_df[!, "sc$(i)_pos_3"] = [sol.u[t].sc[i].pos[3] for t in 1:length(sol.t)]
            results_df[!, "sc$(i)_vel_1"] = [sol.u[t].sc[i].vel[1] for t in 1:length(sol.t)]
            results_df[!, "sc$(i)_vel_2"] = [sol.u[t].sc[i].vel[2] for t in 1:length(sol.t)]
            results_df[!, "sc$(i)_vel_3"] = [sol.u[t].sc[i].vel[3] for t in 1:length(sol.t)]
            results_df[!, "sc$(i)_mass"] = [sol.u[t].sc[i].mass for t in 1:length(sol.t)]
            if args.mission_configuration.orientation_sim
                results_df[!, "sc$(i)_q_1"] = [sol.u[t].sc[i].q[1] for t in 1:length(sol.t)]
                results_df[!, "sc$(i)_q_2"] = [sol.u[t].sc[i].q[2] for t in 1:length(sol.t)]
                results_df[!, "sc$(i)_q_3"] = [sol.u[t].sc[i].q[3] for t in 1:length(sol.t)]
                results_df[!, "sc$(i)_q_4"] = [sol.u[t].sc[i].q[4] for t in 1:length(sol.t)]
                results_df[!, "sc$(i)_ω_1"] = [sol.u[t].sc[i].ω[1] for t in 1:length(sol.t)]
                results_df[!, "sc$(i)_ω_2"] = [sol.u[t].sc[i].ω[2] for t in 1:length(sol.t)]
                results_df[!, "sc$(i)_ω_3"] = [sol.u[t].sc[i].ω[3] for t in 1:length(sol.t)]
            end
        end
        CSV.write("simulation_results.csv", results_df)
    end

    if return_solution
        return sol
    end
    return nothing
end

function spacecraft_dynamics!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    # Unpack the state vector
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    debug_control = p.shared_buffers.debug_control[]
    p.shared_buffers.current_time[] = t
    minbatch = Int(ceil(length(spacecraft) / Polyester.num_cores())) # Determine the batch size for LoopVectorization based on the number of spacecraft and available CPU cores
    # Loop over each spacecraft and compute its dynamics
    @batch minbatch=minbatch for i in eachindex(sc_state)
        if !p.is_active[i]
            sc_du[i] .= 0.0 # Set the derivatives to zero for inactive spacecraft
            continue
        end
        @views begin
            sc_view = sc_state[i]
            du_view = sc_du[i]
            # println("Computing dynamics for spacecraft $i at time $t seconds...")
            # Compute forces and torques using the dynamic effectors
            forces = MVector{3, Float64}(0.0, 0.0, 0.0)
            torques = MVector{3, Float64}(0.0, 0.0, 0.0)
            mass_rate = 0.0
            @inbounds for effector in dynamic_effectors
                force, torque = SimulationModel.calcForceTorque(effector, sc_view, p, i)
                forces .+= force
                torques .+= torque
            end

            # Compute control forces and torques using the control effectors (if any)
                @inbounds for control_effector in p.args.control_model.control_effectors
                    control_force, control_torque = SimulationModel.calcControlForceTorque(control_effector, sc_view, p, i, t)
                    control_mass_rate = SimulationModel.calcControlMassFlowRate(control_effector, sc_view, p, i, t)
                    if debug_control && (norm(control_force) > 0.0 || norm(control_torque) > 0.0)
                        println("Applying control effect for spacecraft $i at time $t seconds:")
                        println("  Control force: $control_force")
                        # println("  Control torque: $control_torque")
                    end
                # println("Control force for spacecraft $i at time $t seconds: $control_force")
                forces .+= control_force
                torques .+= control_torque
                mass_rate += isfinite(control_mass_rate) ? control_mass_rate : 0.0
            end

            # Update the derivatives of position and velocity
            du_view.pos .= sc_view.vel
            du_view.vel .= forces / sc_view.mass
            du_view.mass = mass_rate

            if p.args.mission_configuration.orientation_sim
                # Update the derivatives of orientation (quaternion) and angular velocity
                du_view.q .= 0.5 * quat_mult(SVector{4, Float64}(sc_view.ω..., 0.0), sc_view.q)
                du_view.ω .= torques / sc_view.mass # Simplified rotational dynamics (assuming unit inertia)
            end

            # Update heat loads based on current state and environment (placeholder logic)
            du_view.heat_loads .= norm(forces) * 1e-6 # Example: heat load proportional to force magnitude
        end
    end
end # function spacecraft_dynamics!

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
