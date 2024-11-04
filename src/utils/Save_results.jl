include("../config.jl")

function save_results(time, ratio)
    initial_time = 0

    if length(solution.orientation.time) == 0
        config.prev_step_integrator = 0.0
        initial_time = config.prev_step_integrator
        config.initial_time_saved = 0.0
    else
        if config.counter_integrator == 1
            config.prev_step_integrator =  0.0
            intial_time = config.prev_step_integrator
            config.initial_time_saved = solution.orientation.time[end]
        end
    end

    n_variable_to_save = length(config.solution_intermediate[1]) - 1
    range_time = [item[1] for item in config.solution_intermediate]
    results = zeros(n_variable_to_save, 1)

    t = []
    i = 0

    index_prev = 0

    for true_time in time
        if i % ratio == 0
            index = 
            append!(t, range_time[index+index_prev] + initial_time)
            range_solution = reshape(config.solution_intermediate[index+index_prev], (n_variable_to_save + 1, 1))

            if length(t) == 1
                results = range_solution[2:end]
            else
                results = append!(results, range_solution[2:end])
            end
            
            index_prev = index
        end

        i += 1
    end

    time_0 = t[end]
    config.prev_step_integrator = time_0
    config.solution_intermediate = []

    t = [i + config.initial_time_saved for i in t]


    ## SAVE RESULTS GLOBAL
    # Orientation
    append!(solution.orientation.time, t)
    append!(solution.orientation.year, results[1])
    append!(solution.orientation.month, results[2])
    append!(solution.orientation.day, results[3])
    append!(solution.orientation.hour, results[4])
    append!(solution.orientation.minute, results[5])
    append!(solution.orientation.second, results[6])
    append!(solution.orientation.number_of_passage, results[7])
    append!(solution.orientation.pos_ii[1], results[8])
    append!(solution.orientation.pos_ii[2], results[9])
    append!(solution.orientation.pos_ii[3], results[10])
    append!(solution.orientation.vel_ii[1], results[11])
    append!(solution.orientation.vel_ii[2], results[12])
    append!(solution.orientation.vel_ii[3], results[13])
    append!(solution.orientation.pos_ii_mag, results[14])
    append!(solution.orientation.vel_ii_mag, results[15])

    append!(solution.orientation.pos_pp[1], results[16])
    append!(solution.orientation.pos_pp[2], results[17])
    append!(solution.orientation.pos_pp[3], results[18])
    append!(solution.orientation.pos_pp_mag, results[19])
    append!(solution.orientation.vel_pp[1], results[20])
    append!(solution.orientation.vel_pp[2], results[21])
    append!(solution.orientation.vel_pp[3], results[22])
    append!(solution.orientation.vel_pp_mag, results[23])

    append!(solution.orientation.oe[1], results[24])
    append!(solution.orientation.oe[2], results[25])
    append!(solution.orientation.oe[3], results[26])
    append!(solution.orientation.oe[4], results[27])
    append!(solution.orientation.oe[5], results[28])
    append!(solution.orientation.oe[6], results[29])

    append!(solution.orientation.lat, results[30])
    append!(solution.orientation.lon, results[31])
    append!(solution.orientation.alt, results[32])
    appwnd!(solution.orientation.γ_ii, results[33])
    append!(solution.orientation.γ_pp, results[34])

    append!(solution.orientation.h_ii[1], results[35])
    append!(solution.orientation.h_ii[2], results[36])
    append!(solution.orientation.h_ii[3], results[37])
    append!(solution.orientation.h_pp[1], results[38])
    append!(solution.orientation.h_pp[2], results[39])
    append!(solution.orientation.h_pp[3], results[40])
    append!(solution.orientation.h_ii_mag, results[41])
    append!(solution.orientation.h_pp_mag, results[42])

    append!(solution.orientation.uD[1], results[43])
    append!(solution.orientation.uD[2], results[44])
    append!(solution.orientation.uD[3], results[45])
    append!(solution.orientation.uE[1], results[46])
    append!(solution.orientation.uE[2], results[47])
    append!(solution.orientation.uE[3], results[48])
    append!(solution.orientation.uN[1], results[49])
    append!(solution.orientation.uN[2], results[50])
    append!(solution.orientation.uN[3], results[51])
    append!(solution.orientation.vN, results[52])
    append!(solution.orientation.vE, results[53])
    append!(solution.orientation.azi_pp, results[54])

    # Physical properties
    append!(solution.physical_properties.ρ, results[55])
    append!(solution.physical_properties.T, results[56])
    append!(solution.physical_properties.p, results[57])
    append!(solution.physical_properties.wind[1], results[58])
    append!(solution.physical_properties.wind[2], results[59])
    append!(solution.physical_properties.wind[3], results[60])
    append!(solution.physical_properties.cL, results[61])
    append!(solution.physical_properties.cD, results[62])
    append!(solution.physical_properties.α, results[63])
    append!(solution.physical_properties.S, results[64])

    # Performance
    append!(solution.performance.mass, results[65])
    append!(solution.performance.heat_rate, results[66])
    append!(solution.performance.heat_load, results[67])
    append!(solution.performance.T_r, results[68])
    append!(solution.performance.q, results[69])

    # Forces
    append!(solution.forces.gravity_ii[1], results[70])
    append!(solution.forces.gravity_ii[2], results[71])
    append!(solution.forces.gravity_ii[3], results[72])
    append!(solution.forces.drag_pp[1], results[73])
    append!(solution.forces.drag_pp[2], results[74])
    append!(solution.forces.drag_pp[3], results[75])
    append!(solution.forces.drag_ii[1], results[76])
    append!(solution.forces.drag_ii[2], results[77])
    append!(solution.forces.drag_ii[3], results[78])
    append!(solution.forces.lift_pp[1], results[79])
    append!(solution.forces.lift_pp[2], results[80])
    append!(solution.forces.lift_pp[3], results[81])
    append!(solution.forces.lift_ii[1], results[82])
    append!(solution.forces.lift_ii[2], results[83])
    append!(solution.forces.lift_ii[3], results[84])
    append!(solution.forces.force_ii[1], results[85])
    append!(solution.forces.force_ii[2], results[86])
    append!(solution.forces.force_ii[3], results[87])
    append!(solution.forces.energy, results[88])

    # Simulation
    append!(solution.simulation.MC_seed, results[89])
    append!(solution.simulation.drag_passage, results[90])

    return time_0
end

function clean_results()
    solution.orientation.time = []
    solution.orientation.year = []
    solution.orientation.month = []
    solution.orientation.day = []
    solution.orientation.hour = []
    solution.orientation.minute = []
    solution.orientation.second = []
    solution.orientation.number_of_passage = []
    solution.orientation.pos_ii[1] = []
    solution.orientation.pos_ii[2] = []
    solution.orientation.pos_ii[3] = []
    solution.orientation.vel_ii[1] = []
    solution.orientation.vel_ii[2] = []
    solution.orientation.vel_ii[3] = []
    solution.orientation.pos_ii_mag = []
    solution.orientation.vel_ii_mag = []

    solution.orientation.pos_pp[1] = []
    solution.orientation.pos_pp[2] = []
    solution.orientation.pos_pp[3] = []
    solution.orientation.pos_pp_mag = []
    solution.orientation.vel_pp[1] = []
    solution.orientation.vel_pp[2] = []
    solution.orientation.vel_pp[3] = []
    solution.orientation.vel_pp_mag = []

    solution.orientation.oe[1] = []
    solution.orientation.oe[2] = []
    solution.orientation.oe[3] = []
    solution.orientation.oe[4] = []
    solution.orientation.oe[5] = []
    solution.orientation.oe[6] = []

    solution.orientation.lat = []
    solution.orientation.lon = []
    solution.orientation.alt = []
    solution.orientation.γ_ii = []
    solution.orientation.γ_pp = []

    solution.orientation.h_ii[1] = []
    solution.orientation.h_ii[2] = []
    solution.orientation.h_ii[3] = []
    solution.orientation.h_pp[1] = []
    solution.orientation.h_pp[2] = []
    solution.orientation.h_pp[3] = []
    solution.orientation.h_ii_mag = []
    solution.orientation.h_pp_mag = []

    solution.orientation.uD[1] = []
    solution.orientation.uD[2] = []
    solution.orientation.uD[3] = []
    solution.orientation.uE[1] = []
    solution.orientation.uE[2] = []
    solution.orientation.uE[3] = []
    solution.orientation.uN[1] = []
    solution.orientation.uN[2] = []
    solution.orientation.uN[3] = []
    solution.orientation.vN = []
    solution.orientation.vE = []
    solution.orientation.azi_pp = []

    # Physical properties
    solution.physical_properties.ρ = []
    solution.physical_properties.T = []
    solution.physical_properties.p = []
    solution.physical_properties.wind[1] = []
    solution.physical_properties.wind[2] = []
    solution.physical_properties.wind[3] = []
    solution.physical_properties.cL = []
    solution.physical_properties.cD = []
    solution.physical_properties.α = []
    solution.physical_properties.S = []

    # Performance
    solution.performance.mass = []
    solution.performance.heat_rate = []
    solution.performance.heat_load = []
    solution.performance.T_r = []
    solution.performance.q = []

    # Forces
    solution.forces.gravity_ii[1] = []
    solution.forces.gravity_ii[2] = []
    solution.forces.gravity_ii[3] = []
    solution.forces.drag_pp[1] = []
    solution.forces.drag_pp[2] = []
    solution.forces.drag_pp[3] = []
    solution.forces.drag_ii[1] = []
    solution.forces.drag_ii[2] = []
    solution.forces.drag_ii[3] = []
    solution.forces.lift_pp[1] = []
    solution.forces.lift_pp[2] = []
    solution.forces.lift_pp[3] = []
    solution.forces.lift_ii[1] = []
    solution.forces.lift_ii[2] = []
    solution.forces.lift_ii[3] = []
    solution.forces.force_ii[1] = []
    solution.forces.force_ii[2] = []
    solution.forces.force_ii[3] = []
    solution.forces.energy = []

    # Simulation
    solution.simulation.MC_seed = []
    solution.simulation.drag_passage = []

    # Closed form
    solution.closed_form.t_cf = []
    solution.closed_form.h_cf = []
    solution.closed_form.γ_cf = []
    solution.closed_form.v_cf = []

    return
end