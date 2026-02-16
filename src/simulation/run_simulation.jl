using LinearAlgebra
using StaticArrays
using LoopVectorization
using ComponentArrays
using DifferentialEquations
using CSV
# using .SimulationModel: SimulationConfiguration

function run_simulation(args::SimulationConfiguration)
    # Set up the model and initial conditions
    initial_conditions = build_initial_conditions(args)
    println("Initial conditions:")
    println(initial_conditions)

    # Define the ODE problem
    p = ODEParams(args=args)
    prob = ODEProblem(spacecraft_dynamics!, initial_conditions, (0.0, args.mission_configuration.mission_time), p)
    # Solve the ODE problem
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
end

function spacecraft_dynamics!(du::ComponentVector, u::ComponentVector, p::ODEParams, t::Float64)
    # Unpack the state vector
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    # Loop over each spacecraft and compute its dynamics
    for i in eachindex(sc_state)
        sc_view = sc_state[i]
        du_view = sc_du[i]

        # Compute forces and torques using the dynamic effectors
        forces = MVector{3, Float64}(0.0, 0.0, 0.0)
        torques = MVector{3, Float64}(0.0, 0.0, 0.0)
        @inbounds for effector in dynamic_effectors
           force, torque = calcForceTorque(effector, sc_view, p)
           forces .+= force
           torques .+= torque
        end

        # Update the derivatives of position and velocity
        du_view.pos .= sc_view.vel
        du_view.vel .= forces / sc_view.mass

        if p.args.mission_configuration.orientation_sim
            # Update the derivatives of orientation (quaternion) and angular velocity
            du_view.q .= 0.5 * quat_mult(SVector{4, Float64}(sc_view.ω..., 0.0), sc_view.q)
            du_view.ω .= torques / sc_view.mass # Simplified rotational dynamics (assuming unit inertia)
        end

        # Update heat loads based on current state and environment (placeholder logic)
        du_view.heat_loads .= norm(forces) * 1e-6 # Example: heat load proportional to force magnitude
    end
end # function spacecraft_dynamics!

function build_initial_conditions(args::SimulationConfiguration)::ComponentVector
    # 1. Build the structure (Axis) based on each spacecraft's unique body count
    # This identifies exactly how many heat_load slots each SC needs
    sc_shapes = map(args.dynamics_model.spacecraft) do sc
        # Get the number of bodies for this specific spacecraft
        n_bodies = length(sc.children) + 1 # +1 for the root body itself 
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