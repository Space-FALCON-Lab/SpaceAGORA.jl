include("../utils/Define_mission.jl")
include("../utils/MonteCarlo_set.jl")
include("../utils/Initial_cond_calc.jl")
include("Set_and_run.jl")
include("../utils/Mission_anim.jl")
# include("../utils/plot_data_multi.jl")

using Dash
using DashCoreComponents
using DashHtmlComponents
import .config



function run_orbitalelements(args)
    apoapsis, periapsis_alt, inclination, Ω, ω = collect(range(start=round(args[:ra_initial_a]), stop=round(args[:ra_initial_b]), step=round(args[:ra_step]))), 
                                                 collect(range(start=round(args[:hp_initial_a]), stop=round(args[:hp_initial_b]), step=round(args[:hp_step]))), 
                                                 args[:inclination], args[:Ω], args[:ω]
    if args[:initial_condition_type] == 2
        apoapsis, periapsis_alt = ic_calculation_ae(apoapsis, periapsis_alt, args)
    end
    final_apoapsis = args[:final_apoapsis]

    for periapsis_item in periapsis_alt
        for apoapsis_item in apoapsis

            if Bool(args[:print_res])
                println("Apoapsis Radius: " * string(apoapsis_item/10^3) * " km, Periapsis Altitude: " * string(periapsis_item/10^3) * " km")  
            end

            state = Dict()

            MC, count, args = MonteCarlo_setting(args)

            for mc_index in range(args[:initial_montecarlo_number], args[:montecarlo_size], step=1)
                state[:Apoapsis], state[:Periapsis], state[:Inclination], state[:Ω], state[:ω], state[:Final_sma] = apoapsis_item, Float64(periapsis_item*1e-3), inclination, Ω, ω, final_apoapsis

                args[:simulation_filename] = "Results_ctrl=" * string(args[:control_mode]) * "_ra=" * string(Int64(round(apoapsis_item/1e3))) * "_rp=" * string(Float64(periapsis_item/1e3)) * "_hl=" * string(args[:max_heat_rate])

                if args[:montecarlo] == true
                    args = MonteCarlo_setting_passage(mc_index, args)
                end

                aerobraking_campaign(args, state)
                MonteCarlo_append(MC, args, count)
            end

            if args[:montecarlo] == true
                MonteCarlo_save(args, state, MC)
            end
        end
    end
end

function run_vgamma(args)
    γ_0, v_0, inclination, Ω, ω = collect(range(start=round(args[:γ_initial_a]*100), stop=round(args[:γ_initial_a]*100), step=round(args[:γ_step]*100))), 
                                  collect(range(start=round(args[:v_initial_a]), stop=round(args[:v_initial_b]), step=round(args[:v_step]))),
                                  args[:inclination], args[:Ω], args[:ω]
    final_apoapsis = args[:final_apoapsis]

    for γ in γ_0
        γ = -γ/100

        for v in v_0
            state = Dict()
            planet = planet_data(args[:planet])
            apoapsis, periapsis_alt = ic_calculation_rptoae(planet, γ, v, args)

            if Bool(args[:print_res])
                println("Velocity: " * string(v) * " m/s, Flight-Path Angle: " * string(γ) * " deg")
            end

            MC, count, args = MonteCarlo_setting(args)

            for mc_index in range(args[:initial_montecarlo_number], args[:montecarlo_size])
                state[:Apoapsis], state[:Periapsis], state[:Inclination], state[:Ω], state[:ω], state[:Final_sma] = apoapsis, periapsis_alt * 1e-3, inclination, Ω, ω, final_apoapsis

                args[:simulation_filename] = "Results_ctrl=" * string(args[:control_mode]) * "_v=" * string(v) * "_gamma=" * string(abs(γ)) * "_" * string(args[:α]) * "deg"

                if args[:montecarlo] == true
                    args = MonteCarlo_setting_passage(mc_index, args)
                end

                aerobraking_campaign(args, state)
                MonteCarlo_append(MC, args, count)
            end

            if args[:montecarlo] == true
                MonteCarlo_save(args, state, MC)
            end
        end
    end
end

function run_orbitalelements_ae(args)
    #Used to run orbital simulation if initial conditions are in terms of sma and eccentricity

    # Generate ranges for semi-major axis (a) and eccentricity (e) based on input arguments
    a, e, inclination, Ω, ω = collect(range(start=round(args[:a_initial_a]), stop=round(args[:a_initial_b]), step=round(args[:a_step]))), 
                             collect(range(start=args[:e_initial_a], stop=args[:e_initial_b], step=args[:e_step])), 
                             args[:inclination], args[:Ω], args[:ω]

    # Extract final apoapsis value from arguments
    final_apoapsis = args[:final_apoapsis]

    # Get planet data using the planet name from arguments
    planet = planet_data(args[:planet])

    # Loop over all eccentricity values
    for e_item in e
        # Loop over all semi-major axis values
        for a_item in a
            # Calculate apoapsis and periapsis from a and e for the given planet and arguments
            apoapsis_item, periapsis_item = ic_calculation_ae(planet, a_item, e_item, args)

            # Optionally print the current apoapsis and periapsis values in km
            if Bool(args[:print_res])
                println("Apoapsis Radius: " * string(apoapsis_item/10^3) * " km, Periapsis Altitude: " * string(periapsis_item/10^3) * " km")  
            end

            # Initialize state dictionary for simulation parameters
            state = Dict()

            # Set up Monte Carlo simulation parameters and counter
            MC, count, args = MonteCarlo_setting(args)

            # Loop over Monte Carlo runs
            for mc_index in range(args[:initial_montecarlo_number], args[:montecarlo_size], step=1)
                # Assign current simulation state values
                state[:Apoapsis], state[:Periapsis], state[:Inclination], state[:Ω], state[:ω], state[:Final_sma] = apoapsis_item, Float64(periapsis_item*1e-3), inclination, Ω, ω, final_apoapsis

                # Generate a filename for simulation results based on current parameters
                args[:simulation_filename] = "Results_ctrl=" * string(args[:control_mode]) * "_ra=" * string(Int64(round(apoapsis_item/1e3))) * "_rp=" * string(Float64(periapsis_item/1e3)) * "_hl=" * string(args[:max_heat_rate])

                # If Monte Carlo mode is enabled, update arguments for this run
                if args[:montecarlo] == true
                    args = MonteCarlo_setting_passage(mc_index, args)
                end

                # Run the aerobraking campaign simulation with current arguments and state
                aerobraking_campaign(args, state)

                # Append results of this run to the Monte Carlo collection
                MonteCarlo_append(MC, args, count)
            end

            # If Monte Carlo mode is enabled, save all results for this set of runs
            if args[:montecarlo] == true
                MonteCarlo_save(args, state, MC)
            end
        end
    end
end

function run_sc_vehicles(args)
    """
    Running Space AGORA simulation designed for constellation simulation or otherwise non atmospheric spacecraft
    inputs: args struct
    outputs: Csv file with results
    """

    #extracting initial conditions from the args dict
    target_objs_dict = args[:target_objects]
    spacecraft_dict = args[:spacecraft_buses]
    target_initial_conditions = args[:target_initial_conditions]
    spacecraft_initial_conditions = args[:spacecraft_initial_conditions]

    # Get planet data using the planet name from arguments
    planet = planet_data(args[:planet])

    #storage for state values for targets and spacecraft
    sc_states = Dict()
    target_states = Dict()
    space_objects_dict = Dict()


    #initialize SpaceAGORA dashboard
    app = dash(external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"])
    app.layout = html_div([
        html_h1("SpaceAGORA Dashboard"),
        html_div("Dash.jl: Julia interface for Dash")
    ])

    # Run Dash.jl server asynchronously
    @async run_server(app, "0.0.0.0", 8050)
    

    #assign initial state values for targets
    for target in keys(target_objs_dict)
        target_state = target_initial_conditions[target]
        state = Dict()
        apoapsis_item, periapsis_item = ic_calculation_ae(planet, args[:target_initial_conditions][target][:a_initial_a], args[:target_initial_conditions][target][:e_initial_a], args)
        inclination, Ω, ω = args[:inclination], args[:Ω], args[:ω]
        final_apoapsis = args[:final_apoapsis]

        if Bool(args[:print_res])
            println("Apoapsis Radius: " * string(apoapsis_item/10^3) * " km, Periapsis Altitude: " * string(periapsis_item/10^3) * " km")  
        end

        state[:Apoapsis], state[:Periapsis], state[:Inclination], state[:Ω], state[:ω], state[:Final_sma] = apoapsis_item, Float64(periapsis_item*1e-3), inclination, Ω, ω, final_apoapsis

        target_states[target] = state
        target_objs_dict[target].sc_states = state
        # target_objs_dict[target].sc_state_history = DataFrame([state])
        space_objects_dict[target] = target_objs_dict[target]
        print(keys(space_objects_dict))
    end

    #assign initial state values for spacecraft
    for sc in keys(spacecraft_dict)
        sc_state = spacecraft_initial_conditions[sc]
        state = Dict()
        apoapsis_item, periapsis_item = ic_calculation_ae(planet, args[:spacecraft_initial_conditions][sc][:a_initial_a], args[:spacecraft_initial_conditions][sc][:e_initial_a], args)
        inclination, Ω, ω = args[:inclination], args[:Ω], args[:ω]
        final_apoapsis = args[:final_apoapsis]

        if Bool(args[:print_res])
            println("Apoapsis Radius: " * string(apoapsis_item/10^3) * " km, Periapsis Altitude: " * string(periapsis_item/10^3) * " km")  
        end

        state[:Apoapsis], state[:Periapsis], state[:Inclination], state[:Ω], state[:ω], state[:Final_sma] = apoapsis_item, Float64(periapsis_item*1e-3), inclination, Ω, ω, final_apoapsis

        sc_states[sc] = state
        spacecraft_dict[sc].sc_states = state
        # spacecraft_dict[sc].sc_state_history = [state]
        space_objects_dict[sc] = spacecraft_dict[sc]
    end

    #resorting back into args struct
    args[:spacecraft] = spacecraft_dict
    args[:target_objects] = target_objs_dict
    args[:space_objects_dict] = space_objects_dict

    # Dictionary to store all outputs of aerobraking_campaign for each spacecraft
    aerobraking_outputs = Dict{Any, Any}()

    # This will be populated inside the simulation threads below
    args[:aerobraking_outputs] = aerobraking_outputs

    # Begin spacecraft simulation tasks (true parallelism)
    tasks = Dict{Any, Task}()

    for obj_id in collect(keys(space_objects_dict))
        obj = args[:space_objects_dict][obj_id]

        # Make a private config for this iteration
        cfg = deepcopy(config.config_data())

        # Spawn each simulation as a separate task on available threads
        tasks[obj_id] = Threads.@spawn begin
            println("begin aero campaign for obj_id: ", obj_id)
            result = config.with_config(cfg) do
                state, m, name, args, temp_name = aerobraking_campaign(args, obj.sc_states, obj_id)
                (state=state, m=m, name=name, args=args, temp_name=temp_name)
            end
            println("Aero campaign complete for " * string(obj_id))
            result
        end
    end

    # Wait for all tasks to finish and collect results
    for (obj_id, t) in tasks
        aerobraking_outputs[obj_id] = fetch(t)
    end

    #writing all simulation results to csv files

    # for obj_id in keys(space_objects_dict)
    #     obj = space_objects_dict[obj_id]
    #     sc_state_history = obj.sc_state_history

    #     filename = string(obj_id, "_state_data.csv")
    #     CSV.write(filename, Tables.table(sc_state_history), writeheader=false)

    # end

    # #visualize simulation run
    # start_viz_dashboard(args,space_objects_dict)

    
    if args[:plot] == true
        # for obj_id in keys(aerobraking_outputs)
        #     output = aerobraking_outputs[obj_id]
        #     plots(output.state, output.m, output.name, output.args, output.temp_name)
        # end
        plots_multi(aerobraking_outputs)
        # plots(state, m, name, args, temp_name)
    end
    



end

function run_analysis(args)
    config.reset_all_configs!()
    config.model().body = args[:spacecraft_model]
    args = def_miss(args)

    #create thread local configs
    # thread_configs = [deepcopy(config) for _ in 1:Threads.nthreads()]
    # get_config() = thread_configs[Threads.threadid()]
    # args[:get_config] = config.get_config()

    if args[:initial_condition_type] == 1 && (Bool(args[:drag_passage]) || args[:body_shape] == "Blunted Cone")
        run_vgamma(args)
    elseif args[:initial_condition_type] == 0
        run_orbitalelements(args)
    elseif args[:initial_condition_type] == 2
        run_orbitalelements_ae(args)
    elseif args[:initial_condition_type] == 3
        #Running the swarm/constellation mission type
        run_sc_vehicles(args)
    end

    if Bool(args[:passresults])
        return config.solution
    else
        return true
    end
end