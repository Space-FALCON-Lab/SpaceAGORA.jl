import .config

function save_results(time, ratio)
    initial_time = 0

    if length(config.solution.orientation.time) == 0
        config.cnf.prev_step_integrator = 0.0
        initial_time = config.cnf.prev_step_integrator
        config.cnf.initial_time_saved = 0.0
    else
        if config.cnf.counter_integrator == 1
            config.cnf.prev_step_integrator =  0.0
            intial_time = config.cnf.prev_step_integrator
            config.cnf.initial_time_saved = config.solution.orientation.time[end]
        end
    end

    n_variable_to_save = length(config.cnf.solution_intermediate[1]) - 1
    range_time = [item[1] for item in config.cnf.solution_intermediate]
    results = Array{Float64}(undef, (n_variable_to_save, 1)) # Array{Float64}(undef, n_variable_to_save, 1) # zeros(n_variable_to_save, 1)

    # println(length(config.cnf.solution_intermediate[1]))
    # println(n_variable_to_save)

    t = []
    i = 0

    index_prev = 1

    # println(length(config.cnf.solution_intermediate))
    # println(range_time)

    # println(time)   

    for true_time in time
        if i % ratio == 0
            index = findfirst(x -> x == true_time, range_time[index_prev:end])
            # println(index)
            append!(t, range_time[index+index_prev] + initial_time)
            range_solution = reshape(config.cnf.solution_intermediate[index+index_prev], (n_variable_to_save + 1, 1))
            # println(range_solution)

            if length(t) == 1
                results = range_solution[2:end]
            else
                results = hcat(results, range_solution[2:end])
            end
            
            index_prev = index
        end

        i += 1
    end

    # println(t)

    time_0 = t[end]
    config.cnf.prev_step_integrator = time_0
    config.cnf.solution_intermediate = []

    t = [i + config.cnf.initial_time_saved for i in t]

    # println(results)


    ## SAVE RESULTS GLOBAL
    # Orientation
    # push!(config.solution.orientation.time, t)
    # push!(config.solution.orientation.year, results[1])
    # push!(config.solution.orientation.month, results[2])
    # push!(config.solution.orientation.day, results[3])
    # push!(config.solution.orientation.hour, results[4])
    # push!(config.solution.orientation.minute, results[5])
    # push!(config.solution.orientation.second, results[6])
    # push!(config.solution.orientation.number_of_passage, results[7])
    # push!(config.solution.orientation.pos_ii[1], results[8])
    # push!(config.solution.orientation.pos_ii[2], results[9])
    # push!(config.solution.orientation.pos_ii[3], results[10])
    # push!(config.solution.orientation.vel_ii[1], results[11])
    # push!(config.solution.orientation.vel_ii[2], results[12])
    # push!(config.solution.orientation.vel_ii[3], results[13])
    # push!(config.solution.orientation.pos_ii_mag, results[14])
    # push!(config.solution.orientation.vel_ii_mag, results[15])

    # push!(config.solution.orientation.pos_pp[1], results[16])
    # push!(config.solution.orientation.pos_pp[2], results[17])
    # push!(config.solution.orientation.pos_pp[3], results[18])
    # push!(config.solution.orientation.pos_pp_mag, results[19])
    # push!(config.solution.orientation.vel_pp[1], results[20])
    # push!(config.solution.orientation.vel_pp[2], results[21])
    # push!(config.solution.orientation.vel_pp[3], results[22])
    # push!(config.solution.orientation.vel_pp_mag, results[23])

    # push!(config.solution.orientation.oe[1], results[24])
    # push!(config.solution.orientation.oe[2], results[25])
    # push!(config.solution.orientation.oe[3], results[26])
    # push!(config.solution.orientation.oe[4], results[27])
    # push!(config.solution.orientation.oe[5], results[28])
    # push!(config.solution.orientation.oe[6], results[29])

    # push!(config.solution.orientation.lat, results[30])
    # push!(config.solution.orientation.lon, results[31])
    # push!(config.solution.orientation.alt, results[32])
    # push!(config.solution.orientation.γ_ii, results[33])
    # push!(config.solution.orientation.γ_pp, results[34])

    # push!(config.solution.orientation.h_ii[1], results[35])
    # push!(config.solution.orientation.h_ii[2], results[36])
    # push!(config.solution.orientation.h_ii[3], results[37])
    # push!(config.solution.orientation.h_pp[1], results[38])
    # push!(config.solution.orientation.h_pp[2], results[39])
    # push!(config.solution.orientation.h_pp[3], results[40])
    # push!(config.solution.orientation.h_ii_mag, results[41])
    # push!(config.solution.orientation.h_pp_mag, results[42])

    # push!(config.solution.orientation.uD[1], results[43])
    # push!(config.solution.orientation.uD[2], results[44])
    # push!(config.solution.orientation.uD[3], results[45])
    # push!(config.solution.orientation.uE[1], results[46])
    # push!(config.solution.orientation.uE[2], results[47])
    # push!(config.solution.orientation.uE[3], results[48])
    # push!(config.solution.orientation.uN[1], results[49])
    # push!(config.solution.orientation.uN[2], results[50])
    # push!(config.solution.orientation.uN[3], results[51])
    # push!(config.solution.orientation.vN, results[52])
    # push!(config.solution.orientation.vE, results[53])
    # push!(config.solution.orientation.azi_pp, results[54])

    # # Physical properties
    # push!(config.solution.physical_properties.ρ, results[55])
    # push!(config.solution.physical_properties.T, results[56])
    # push!(config.solution.physical_properties.p, results[57])
    # push!(config.solution.physical_properties.wind[1], results[58])
    # push!(config.solution.physical_properties.wind[2], results[59])
    # push!(config.solution.physical_properties.wind[3], results[60])
    # push!(config.solution.physical_properties.cL, results[61])
    # push!(config.solution.physical_properties.cD, results[62])
    # push!(config.solution.physical_properties.α, results[63])
    # push!(config.solution.physical_properties.S, results[64])

    # # Performance
    # push!(config.solution.performance.mass, results[65])
    # push!(config.solution.performance.heat_rate, results[66])
    # push!(config.solution.performance.heat_load, results[67])
    # push!(config.solution.performance.T_r, results[68])
    # push!(config.solution.performance.q, results[69])

    # # Forces
    # push!(config.solution.forces.gravity_ii[1], results[70])
    # push!(config.solution.forces.gravity_ii[2], results[71])
    # push!(config.solution.forces.gravity_ii[3], results[72])
    # push!(config.solution.forces.drag_pp[1], results[73])
    # push!(config.solution.forces.drag_pp[2], results[74])
    # push!(config.solution.forces.drag_pp[3], results[75])
    # push!(config.solution.forces.drag_ii[1], results[76])
    # push!(config.solution.forces.drag_ii[2], results[77])
    # push!(config.solution.forces.drag_ii[3], results[78])
    # push!(config.solution.forces.lift_pp[1], results[79])
    # push!(config.solution.forces.lift_pp[2], results[80])
    # push!(config.solution.forces.lift_pp[3], results[81])
    # push!(config.solution.forces.lift_ii[1], results[82])
    # push!(config.solution.forces.lift_ii[2], results[83])
    # push!(config.solution.forces.lift_ii[3], results[84])
    # push!(config.solution.forces.force_ii[1], results[85])
    # push!(config.solution.forces.force_ii[2], results[86])
    # push!(config.solution.forces.force_ii[3], results[87])
    # push!(config.solution.forces.energy, results[88])

    # # Simulation
    # push!(config.solution.simulation.MC_seed, results[89])
    # push!(config.solution.simulation.drag_passage, results[90])


    append!(config.solution.orientation.time, t)
    append!(config.solution.orientation.year, results[1,:])
    append!(config.solution.orientation.month, results[2,:])
    append!(config.solution.orientation.day, results[3,:])
    append!(config.solution.orientation.hour, results[4,:])
    append!(config.solution.orientation.minute, results[5,:])
    append!(config.solution.orientation.second, results[6,:])
    append!(config.solution.orientation.number_of_passage, results[7,:])
    append!(config.solution.orientation.pos_ii[1], results[8,:])
    append!(config.solution.orientation.pos_ii[2], results[9,:])
    append!(config.solution.orientation.pos_ii[3], results[10,:])
    append!(config.solution.orientation.vel_ii[1], results[11,:])
    append!(config.solution.orientation.vel_ii[2], results[12,:])
    append!(config.solution.orientation.vel_ii[3], results[13,:])
    append!(config.solution.orientation.pos_ii_mag, results[14,:])
    append!(config.solution.orientation.vel_ii_mag, results[15,:])

    append!(config.solution.orientation.pos_pp[1], results[16,:])
    append!(config.solution.orientation.pos_pp[2], results[17,:])
    append!(config.solution.orientation.pos_pp[3], results[18,:])
    append!(config.solution.orientation.pos_pp_mag, results[19,:])
    append!(config.solution.orientation.vel_pp[1], results[20,:])
    append!(config.solution.orientation.vel_pp[2], results[21,:])
    append!(config.solution.orientation.vel_pp[3], results[22,:])
    append!(config.solution.orientation.vel_pp_mag, results[23,:])

    append!(config.solution.orientation.oe[1], results[24,:])
    append!(config.solution.orientation.oe[2], results[25,:])
    append!(config.solution.orientation.oe[3], results[26,:])
    append!(config.solution.orientation.oe[4], results[27,:])
    append!(config.solution.orientation.oe[5], results[28,:])
    append!(config.solution.orientation.oe[6], results[29,:])

    append!(config.solution.orientation.lat, results[30,:])
    append!(config.solution.orientation.lon, results[31,:])
    append!(config.solution.orientation.alt, results[32,:])
    append!(config.solution.orientation.γ_ii, results[33,:])
    append!(config.solution.orientation.γ_pp, results[34,:])

    append!(config.solution.orientation.h_ii[1], results[35,:])
    append!(config.solution.orientation.h_ii[2], results[36,:])
    append!(config.solution.orientation.h_ii[3], results[37,:])
    append!(config.solution.orientation.h_pp[1], results[38,:])
    append!(config.solution.orientation.h_pp[2], results[39,:])
    append!(config.solution.orientation.h_pp[3], results[40,:])
    append!(config.solution.orientation.h_ii_mag, results[41,:])
    append!(config.solution.orientation.h_pp_mag, results[42,:])

    append!(config.solution.orientation.uD[1], results[43,:])
    append!(config.solution.orientation.uD[2], results[44,:])
    append!(config.solution.orientation.uD[3], results[45,:])
    append!(config.solution.orientation.uE[1], results[46,:])
    append!(config.solution.orientation.uE[2], results[47,:])
    append!(config.solution.orientation.uE[3], results[48,:])
    append!(config.solution.orientation.uN[1], results[49,:])
    append!(config.solution.orientation.uN[2], results[50,:])
    append!(config.solution.orientation.uN[3], results[51,:])
    append!(config.solution.orientation.vN, results[52,:])
    append!(config.solution.orientation.vE, results[53,:])
    append!(config.solution.orientation.azi_pp, results[54,:])

    # Physical properties
    append!(config.solution.physical_properties.ρ, results[55,:])
    append!(config.solution.physical_properties.T, results[56,:])
    append!(config.solution.physical_properties.p, results[57,:])
    append!(config.solution.physical_properties.wind[1], results[58,:])
    append!(config.solution.physical_properties.wind[2], results[59,:])
    append!(config.solution.physical_properties.wind[3], results[60,:])
    append!(config.solution.physical_properties.cL, results[61,:])
    append!(config.solution.physical_properties.cD, results[62,:])
    append!(config.solution.physical_properties.α, results[63,:])
    append!(config.solution.physical_properties.S, results[64,:])

    # Performance
    append!(config.solution.performance.mass, results[65,:])
    append!(config.solution.performance.heat_rate, results[66,:])
    append!(config.solution.performance.heat_load, results[67,:])
    append!(config.solution.performance.T_r, results[68,:])
    append!(config.solution.performance.q, results[69,:])

    # Forces
    append!(config.solution.forces.gravity_ii[1], results[70,:])
    append!(config.solution.forces.gravity_ii[2], results[71,:])
    append!(config.solution.forces.gravity_ii[3], results[72,:])
    append!(config.solution.forces.drag_pp[1], results[73,:])
    append!(config.solution.forces.drag_pp[2], results[74,:])
    append!(config.solution.forces.drag_pp[3], results[75,:])
    append!(config.solution.forces.drag_ii[1], results[76,:])
    append!(config.solution.forces.drag_ii[2], results[77,:])
    append!(config.solution.forces.drag_ii[3], results[78,:])
    append!(config.solution.forces.lift_pp[1], results[79,:])
    append!(config.solution.forces.lift_pp[2], results[80,:])
    append!(config.solution.forces.lift_pp[3], results[81,:])
    append!(config.solution.forces.lift_ii[1], results[82,:])
    append!(config.solution.forces.lift_ii[2], results[83,:])
    append!(config.solution.forces.lift_ii[3], results[84,:])
    append!(config.solution.forces.force_ii[1], results[85,:])
    append!(config.solution.forces.force_ii[2], results[86,:])
    append!(config.solution.forces.force_ii[3], results[87,:])
    append!(config.solution.forces.energy, results[88,:])

    # Simulation
    append!(config.solution.simulation.MC_seed, results[89,:])
    append!(config.solution.simulation.drag_passage, results[90,:])

    println(" ")
    println(results[25,:])
    # println(config.solution.orientation.alt)
    println(" ")

    return time_0
end

function clean_results()
    config.solution.orientation.time = []
    config.solution.orientation.year = []
    config.solution.orientation.month = []
    config.solution.orientation.day = []
    config.solution.orientation.hour = []
    config.solution.orientation.minute = []
    config.solution.orientation.second = []
    config.solution.orientation.number_of_passage = []
    config.solution.orientation.pos_ii[1] = []
    config.solution.orientation.pos_ii[2] = []
    config.solution.orientation.pos_ii[3] = []
    config.solution.orientation.vel_ii[1] = []
    config.solution.orientation.vel_ii[2] = []
    config.solution.orientation.vel_ii[3] = []
    config.solution.orientation.pos_ii_mag = []
    config.solution.orientation.vel_ii_mag = []

    config.solution.orientation.pos_pp[1] = []
    config.solution.orientation.pos_pp[2] = []
    config.solution.orientation.pos_pp[3] = []
    config.solution.orientation.pos_pp_mag = []
    config.solution.orientation.vel_pp[1] = []
    config.solution.orientation.vel_pp[2] = []
    config.solution.orientation.vel_pp[3] = []
    config.solution.orientation.vel_pp_mag = []

    config.solution.orientation.oe[1] = []
    config.solution.orientation.oe[2] = []
    config.solution.orientation.oe[3] = []
    config.solution.orientation.oe[4] = []
    config.solution.orientation.oe[5] = []
    config.solution.orientation.oe[6] = []

    config.solution.orientation.lat = []
    config.solution.orientation.lon = []
    config.solution.orientation.alt = []
    config.solution.orientation.γ_ii = []
    config.solution.orientation.γ_pp = []

    config.solution.orientation.h_ii[1] = []
    config.solution.orientation.h_ii[2] = []
    config.solution.orientation.h_ii[3] = []
    config.solution.orientation.h_pp[1] = []
    config.solution.orientation.h_pp[2] = []
    config.solution.orientation.h_pp[3] = []
    config.solution.orientation.h_ii_mag = []
    config.solution.orientation.h_pp_mag = []

    config.solution.orientation.uD[1] = []
    config.solution.orientation.uD[2] = []
    config.solution.orientation.uD[3] = []
    config.solution.orientation.uE[1] = []
    config.solution.orientation.uE[2] = []
    config.solution.orientation.uE[3] = []
    config.solution.orientation.uN[1] = []
    config.solution.orientation.uN[2] = []
    config.solution.orientation.uN[3] = []
    config.solution.orientation.vN = []
    config.solution.orientation.vE = []
    config.solution.orientation.azi_pp = []

    # Physical properties
    config.solution.physical_properties.ρ = []
    config.solution.physical_properties.T = []
    config.solution.physical_properties.p = []
    config.solution.physical_properties.wind[1] = []
    config.solution.physical_properties.wind[2] = []
    config.solution.physical_properties.wind[3] = []
    config.solution.physical_properties.cL = []
    config.solution.physical_properties.cD = []
    config.solution.physical_properties.α = []
    config.solution.physical_properties.S = []

    # Performance
    config.solution.performance.mass = []
    config.solution.performance.heat_rate = []
    config.solution.performance.heat_load = []
    config.solution.performance.T_r = []
    config.solution.performance.q = []

    # Forces
    config.solution.forces.gravity_ii[1] = []
    config.solution.forces.gravity_ii[2] = []
    config.solution.forces.gravity_ii[3] = []
    config.solution.forces.drag_pp[1] = []
    config.solution.forces.drag_pp[2] = []
    config.solution.forces.drag_pp[3] = []
    config.solution.forces.drag_ii[1] = []
    config.solution.forces.drag_ii[2] = []
    config.solution.forces.drag_ii[3] = []
    config.solution.forces.lift_pp[1] = []
    config.solution.forces.lift_pp[2] = []
    config.solution.forces.lift_pp[3] = []
    config.solution.forces.lift_ii[1] = []
    config.solution.forces.lift_ii[2] = []
    config.solution.forces.lift_ii[3] = []
    config.solution.forces.force_ii[1] = []
    config.solution.forces.force_ii[2] = []
    config.solution.forces.force_ii[3] = []
    config.solution.forces.energy = []

    # Simulation
    config.solution.simulation.MC_seed = []
    config.solution.simulation.drag_passage = []

    # Closed form
    config.solution.closed_form.t_cf = []
    config.solution.closed_form.h_cf = []
    config.solution.closed_form.γ_cf = []
    config.solution.closed_form.v_cf = []

    return
end