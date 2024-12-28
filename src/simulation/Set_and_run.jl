include("../physical_models/MonteCarlo_pertrubations.jl")
include("../physical_models/Planet_data.jl")
include("../physical_models/Mission.jl")
include("../utils/Save_csv.jl")
include("../utils/Plot_data.jl")
include("Aerobraking.jl")

import .config

function aerobraking_campaign(args, state)
    save_res = args[:results]

    # Descent towards Mars
    purpose = "Aerobraking around Mars"

    mission = Dict(:Purpose => purpose, 
                   :Gravity_Model => args[:gravity_model], 
                   :Density_Model => args[:density_model], 
                   :Wind => args[:wind],
                   :Aerodynamic_Model => args[:aerodynamic_model],
                   :Thermal_Model => args[:thermal_model],
                   :Control => args[:control_mode],
                   :Firings => args[:thrust_control],
                   :Shape => args[:body_shape],
                   :Monte_Carlo => args[:montecarlo])

    if args[:print_res] == true
        println("Mission is: ", mission)
    end

    ip = mission_def(mission)
    p_class = planet_data(ip.M.planet)

    if args[:gravity_model] == "Inverse Squared"
        p_class.Rp_p = p_class.Rp_e
    end

    # Vehicle - calculation notebook page1
    # Mass
    dry_mass = args[:dry_mass]
    prop_mass = args[:prop_mass]
    mass = dry_mass + prop_mass

    # Spacecraft Shape
    if args[:body_shape] == "Spacecraft"
        area_body = 7.26 # 33.38#7.26# This is recalculated for the new sc config. 11 (look notes)# m^2 2001 Mars Odyssey Aerobraking, Smith & Bell paper
        length_sp = 3.7617 # 11.4#3.7617#5.7 # m solar array length https://www.jpl.nasa.gov/news/press_kits/odysseyarrival.pdf
        height_sp = area_body / length_sp

        # Main Body
        length_ody = 2.2 # m 
        height_ody = 1.7 # m
        width_ody = 2.6 # m

    # Blunted Body Shape
    elseif args[:body_shape] == "Blunted Cone"
        δ = 70 # deg
        nose_radius = 0.6638 # m
        base_radius = 2.65/2 # m
    end

    apoapsis = state[:Apoapsis]

    state[:Periapsis] = p_class.Rp_e + state[:Periapsis]*1e3
    state[:vi] = deg2rad(180.0001)

    if args[:montecarlo] == true
        state = monte_carlo_initial_condition(state, args)
    end

    semimajoraxis_in = (state[:Apoapsis] + state[:Periapsis])/2
    eccentricity_in = (state[:Apoapsis] - state[:Periapsis]) / (state[:Apoapsis] + state[:Periapsis])
    apoapsis = state[:Apoapsis]
    periapsis = state[:Periapsis]

    # Initial Condition
    if args[:drag_passage] == true
        h_0 = 160 * 1e3
    elseif args[:body_shape] == "Blunted Cone"
        h_0 = 120 * 1e3
        args[:AE] = h_0/1e3
        args[:EI] = h_0/1e3
    end

    if Bool(args[:drag_passage]) || args[:body_shape] == "Blunted Cone"
        r = p_class.Rp_e + h_0
        
        state[:vi] = - acos(1 / eccentricity_in * (semimajoraxis_in * (1 - eccentricity_in^2) / r - 1))
        
        if args[:montecarlo] == true
            state = monte_carlo_true_anomaly(state, args)
            apoapsis = state[:Apoapsis]
            periapsis = state[:Periapsis]
        end
    end

    # Initial Model Definition
    # Body
    if args[:body_shape] == "Spacecraft"

        Mass = mass
        length_SA = length_sp
        height_SA = height_sp
        Area_SA = length_SA * height_SA
        length_SC = length_ody
        height_SC = height_ody
        Area_SC = length_ody * height_ody
        Area_tot = Area_SA + Area_SC

        b_class = config.Body(Mass, length_SA, height_SA, Area_SA, length_SC, height_SC, Area_SC, Area_tot, 0.0, 0.0, 0.0)

    elseif args[:body_shape] == "Blunted Cone"

        Mass = mass
        Delta = δ
        NoseRadius = nose_radius
        BaseRadius = base_radius
        Area_tot = pi * BaseRadius^2

        b_class = config.Body(Mass, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Area_tot, Delta, NoseRadius, BaseRadius)

    end

    function initialconditions()
        a = semimajoraxis_in
        e = eccentricity_in
        i = deg2rad(state[:Inclination])
        Ω = deg2rad(state[:Ω])
        ω = deg2rad(state[:ω])
        vi = state[:vi]
        m0 = mass
        year = args[:year]
        month = args[:month]
        day = args[:day]
        hour = args[:hours]
        min = args[:minutes]
        second = args[:secs]
        time_rot = args[:planettime]

        ic = config.Initial_condition(a, e, i, Ω, ω, vi, m0, year, month, day, hour, min, second, time_rot)

        return ic
    end

    ic_class = initialconditions()

    function aerodynamics()
        δ = deg2rad(0)
        α = deg2rad(args[:α])
        thermal_accomodation_factor = args[:thermal_accomodation_factor]
        reflection_coefficient = args[:reflection_coefficient]
        thermal_contact = 0
        heat_rate_limit = args[:max_heat_rate]
        heat_load_limit = args[:max_heat_load]

        a = config.Aerodynamics(δ, α, thermal_accomodation_factor, reflection_coefficient, thermal_contact, heat_rate_limit, heat_load_limit)

        return a        
    end

    a_class = aerodynamics()

    # Engine
    args[:phi] = deg2rad(args[:phi])
    function engine()
        ϕ = args[:phi]
        g_e = 9.81
        T = args[:thrust]
        Isp = 200

        e = config.Engines(ϕ, g_e, T, Isp)

        return e
    end

    e_class = engine()

    function model()
        body = b_class
        planet = p_class
        initialcondition = ic_class
        aerodynamics = a_class
        engine = e_class

        m = config.Model(body, planet, aerodynamics, engine, initialcondition)

        return m
    end

    m = model()

    # Initialization - Reset all the config index for new simulation
    config.cnf.count_aerobraking = 0
    config.cnf.count_overcome_hr = 0
    config.cnf.save_index_heat = 0
    config.cnf.index_propellant_mass = 1
    config.cnf.counter_random = 0

    ##########################################################
    # RUN SIMULATION
    config.cnf.heat_rate_limit = args[:max_heat_rate]
    t_el = @elapsed begin
        aerobraking(ip, m, args)
    end
    ##########################################################

    # println(" ")
    # println(config.solution)
    # println(" ")

    if Bool(args[:print_res])
        println("ρ: " * string(maximum(config.solution.physical_properties.ρ)) * " kg/m^3")
        println("heat rate: " * string(maximum(config.solution.performance.heat_rate)) * " W/cm^2")
    end

    # Save results
    if save_res == 1
        if args[:filename] == 1
            if args[:montecarlo] == true
                folder_name = args[:simulation_filename][1:findfirst!(args[:simulation_filename], "_nMC)")]
            else
                folder_name = args[:simulation_filename]
            end

            name = args[:directory_results] * folder_name * "/" * args[:simulation_filename]
            filename = name * ".csv"
        else
            name = args[:directory_results] * "/Sim" * string(args[:MarsGram_version])

            filename = name * ".csv"
        end
    #     save_csv(filename, args)
    end

    if Bool(args[:print_res])
        println("Elapsed time: " * string(t_el) * " s")
    end

    if args[:plot] == true
        plots(state, m, name, args)
    end
end