using Roots

 # import .config

function propulsion_ic_calcs(m, args, initial_state)
    """

    """

    Δv = args[:delta_v] * (-cos(args[:phi]))

    if Bool(args[:print_res])
        if round(rad2deg(args[:phi])) == 0
            println("LOWER MANEUVER!")
        elseif round(rad2deg(args[:phi])) == 180
            println("RAISE MANEUVER!")
        end
    end

    v_exausted = m.engines.Isp * m.engines.g_e

    if length(config.solution.performance.mass) == 0
        m_i = m.body.mass
    else
        m_i = config.solution.performance.mass[end]
    end

    m_f = m_i / exp(Δv / v_exausted)

    Δm = m_i - m_f

    T = m.engines.T

    Δt = Δm * v_exausted / T

    Δt_half = Δt / 2

    vi = initial_state.vi
    e = initial_state.e
    a = initial_state.a

    # check in initial condition is not in drag pass 300 km
    r = m.planet.Rp_e + args[:EI]*1e3 + 50e3
    vi_300 = acos((a * (1 - e^2) - r) / (e * r))
    E_300 = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(vi_300 / 2))
    E_apo = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(pi / 2))
    Δt_300 = sqrt(a^3 / m.planet.μ) * ((E_apo - e*sin(E_apo)) - (E_300 - e*sin(E_300)))

    if Δt_half > Δt_300
        Δt_half = Δt_300
        Δt = 2 * Δt_half
        Δm = Δt * T / v_exausted
        m_f = m_i - Δm
        Δv = v_exausted * log(m_f / m_i)
        args[:delta_v] = Δv / (-cos(args[:phi]))
        println("-- Thrust Maximum Time Exceeded - Thrust Time and Δv adjusted - NEW Δv = ", args[:delta_v])
    end

    if length(config.solution.performance.mass) == 0
        # inverse problem of Kepler
        t_0 = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(vi / 2))
        M_0 = E_0 - e * sin(E_0)
        n = 1/(sqrt(a^3 / m.planet.μ))
        M_e = n * Δt_half + M_0

        if M_e < pi
            x_0 = M_e + e/2
        else
            x_0 = M_e - e/2
        end

        f(E) = E - e * sin(E) - M_e

        vi_in = find_zero(f, x_0, Roots.Newton())

        if abs(vi_in) >= 2*pi
            temp = round(Int, vi_in/(2*pi))
            vi_in -= temp*2*pi
        end

        if v_in < 0
            vi_in += 2*pi
        end

        if vi_in > pi
            δ = vi_in - pi
            vi_in = pi - δ
        end

        initial_state.vi = vi_in

        return inital_state
    else
        t_fin = config.solution.orientation.time[end]
        t_in_new = t_fin - Δt_half

        dist_list = [abs(t_in_new - p) for p in config.solution.orientation.time]
        temp = sort(collect(Set(dist_list)))

        index_a = findfirst(item -> item == temp[1], dist_list)
        index_b = index_a - 1

        dist_a = dist_list[index_a]
        dist_b = dist_list[index_b]

        t_in_new = (config.solution.orientation.time[index_a] + (config.solution.orientation.time[index_b] - config.solution.orientation.time[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results = ones(90)
        results[1] = (config.solution.orientation.year[index_a] + (config.solution.orientation.year[index_b] - config.solution.orientation.year[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[2] = (config.solution.orientation.month[index_a] + (config.solution.orientation.month[index_b] - config.solution.orientation.month[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[3] = (config.solution.orientation.day[index_a] + (config.solution.orientation.day[index_b] - config.solution.orientation.day[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[4] = (config.solution.orientation.hour[index_a] + (config.solution.orientation.hour[index_b] - config.solution.orientation.hour[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[5] = (config.solution.orientation.minute[index_a] + (config.solution.orientation.minute[index_b] - config.solution.orientation.minute[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[6] = (config.solution.orientation.second[index_a] + (config.solution.orientation.second[index_b] - config.solution.orientation.second[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        
        results[7] = (config.solution.orientation.number_of_passage[index_a] + (config.solution.orientation.number_of_passage[index_b] - config.solution.orientation.number_of_passage[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        
        results[8] = (config.solution.orientation.pos_ii[1][index_a] + (config.solution.orientation.pos_ii[1][index_b] - config.solution.orientation.pos_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[9] = (config.solution.orientation.pos_ii[2][index_a] + (config.solution.orientation.pos_ii[2][index_b] - config.solution.orientation.pos_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[10] = (config.solution.orientation.pos_ii[3][index_a] + (config.solution.orientation.pos_ii[3][index_b] - config.solution.orientation.pos_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[11] = (config.solution.orientation.vel_ii[1][index_a] + (config.solution.orientation.vel_ii[1][index_b] - config.solution.orientation.vel_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[12] = (config.solution.orientation.vel_ii[2][index_a] + (config.solution.orientation.vel_ii[2][index_b] - config.solution.orientation.vel_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[13] = (config.solution.orientation.vel_ii[3][index_a] + (config.solution.orientation.vel_ii[3][index_b] - config.solution.orientation.vel_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[14] = (config.solution.orientation.pos_ii_mag[index_a] + (config.solution.orientation.pos_ii_mag[index_b] - config.solution.orientation.pos_ii_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[15] = (config.solution.orientation.vel_ii_mag[index_a] + (config.solution.orientation.vel_ii_mag[index_b] - config.solution.orientation.vel_ii_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[16] = (config.solution.orientation.pos_pp[1][index_a] + (config.solution.orientation.pos_pp[1][index_b] - config.solution.orientation.pos_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[17] = (config.solution.orientation.pos_pp[2][index_a] + (config.solution.orientation.pos_pp[2][index_b] - config.solution.orientation.pos_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[18] = (config.solution.orientation.pos_pp[3][index_a] + (config.solution.orientation.pos_pp[3][index_b] - config.solution.orientation.pos_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[19] = (config.solution.orientation.pos_pp_mag[index_a] + (config.solution.orientation.pos_pp_mag[index_b] - config.solution.orientation.pos_pp_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[20] = (config.solution.orientation.vel_pp[1][index_a] + (config.solution.orientation.vel_pp[1][index_b] - config.solution.orientation.vel_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[21] = (config.solution.orientation.vel_pp[2][index_a] + (config.solution.orientation.vel_pp[2][index_b] - config.solution.orientation.vel_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[22] = (config.solution.orientation.vel_pp[3][index_a] + (config.solution.orientation.vel_pp[3][index_b] - config.solution.orientation.vel_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[23] = (config.solution.orientation.vel_pp_mag[index_a] + (config.solution.orientation.vel_pp_mag[index_b] - config.solution.orientation.vel_pp_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[24] = (config.solution.orientation.oe[1][index_a] + (config.solution.orientation.oe[1][index_b] - config.solution.orientation.oe[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[25] = (config.solution.orientation.oe[2][index_a] + (config.solution.orientation.oe[2][index_b] - config.solution.orientation.oe[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[26] = (config.solution.orientation.oe[3][index_a] + (config.solution.orientation.oe[3][index_b] - config.solution.orientation.oe[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[27] = (config.solution.orientation.oe[4][index_a] + (config.solution.orientation.oe[4][index_b] - config.solution.orientation.oe[4][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[28] = (config.solution.orientation.oe[5][index_a] + (config.solution.orientation.oe[5][index_b] - config.solution.orientation.oe[5][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[29] = (config.solution.orientation.oe[6][index_a] + (config.solution.orientation.oe[6][index_b] - config.solution.orientation.oe[6][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[30] = (config.solution.orientation.lat[index_a] + (config.solution.orientation.lat[index_b] - config.solution.orientation.lat[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[31] = (config.solution.orientation.lon[index_a] + (config.solution.orientation.lon[index_b] - config.solution.orientation.lon[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[32] = (config.solution.orientation.alt[index_a] + (config.solution.orientation.alt[index_b] - config.solution.orientation.alt[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[33] = (config.solution.orientation.γ_ii[index_a] + (config.solution.orientation.γ_ii[index_b] - config.solution.orientation.γ_ii[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[34] = (config.solution.orientation.γ_pp[index_a] + (config.solution.orientation.γ_pp[index_b] - config.solution.orientation.γ_pp[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[35] = (config.solution.orientation.h_ii[1][index_a] + (config.solution.orientation.h_ii[1][index_b] - config.solution.orientation.h_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[36] = (config.solution.orientation.h_ii[2][index_a] + (config.solution.orientation.h_ii[2][index_b] - config.solution.orientation.h_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[37] = (config.solution.orientation.h_ii[3][index_a] + (config.solution.orientation.h_ii[3][index_b] - config.solution.orientation.h_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[38] = (config.solution.orientation.h_pp[1][index_a] + (config.solution.orientation.h_pp[1][index_b] - config.solution.orientation.h_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[39] = (config.solution.orientation.h_pp[2][index_a] + (config.solution.orientation.h_pp[2][index_b] - config.solution.orientation.h_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[40] = (config.solution.orientation.h_pp[3][index_a] + (config.solution.orientation.h_pp[3][index_b] - config.solution.orientation.h_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[41] = (config.solution.orientation.h_ii_mag[index_a] + (config.solution.orientation.h_ii_mag[index_b] - config.solution.orientation.h_ii_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[42] = (config.solution.orientation.h_pp_mag[index_a] + (config.solution.orientation.h_pp_mag[index_b] - config.solution.orientation.h_pp_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[43] = (config.solution.orientation.uD[1][index_a] + (config.solution.orientation.uD[1][index_b] - config.solution.orientation.uD[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[44] = (config.solution.orientation.uD[2][index_a] + (config.solution.orientation.uD[2][index_b] - config.solution.orientation.uD[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[45] = (config.solution.orientation.uD[3][index_a] + (config.solution.orientation.uD[3][index_b] - config.solution.orientation.uD[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[46] = (config.solution.orientation.uE[1][index_a] + (config.solution.orientation.uE[1][index_b] - config.solution.orientation.uE[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[47] = (config.solution.orientation.uE[2][index_a] + (config.solution.orientation.uE[2][index_b] - config.solution.orientation.uE[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[48] = (config.solution.orientation.uE[3][index_a] + (config.solution.orientation.uE[3][index_b] - config.solution.orientation.uE[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[49] = (config.solution.orientation.uN[1][index_a] + (config.solution.orientation.uN[1][index_b] - config.solution.orientation.uN[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[50] = (config.solution.orientation.uN[2][index_a] + (config.solution.orientation.uN[2][index_b] - config.solution.orientation.uN[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[51] = (config.solution.orientation.uN[3][index_a] + (config.solution.orientation.uN[3][index_b] - config.solution.orientation.uN[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[52] = (config.solution.orientation.vN[index_a] + (config.solution.orientation.vN[index_b] - config.solution.orientation.vN[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[53] = (config.solution.orientation.vE[index_a] + (config.solution.orientation.vE[index_b] - config.solution.orientation.vE[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[54] = (config.solution.orientation.azi_pp[index_a] + (config.solution.orientation.azi_pp[index_b] - config.solution.orientation.azi_pp[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[55] = (config.solution.physical_properties.ρ[index_a] + (config.solution.physical_properties.ρ[index_b] - config.solution.physical_properties.ρ[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[56] = (config.solution.physical_properties.T[index_a] + (config.solution.physical_properties.T[index_b] - config.solution.physical_properties.T[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[57] = (config.solution.physical_properties.p[index_a] + (config.solution.physical_properties.p[index_b] - config.solution.physical_properties.p[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        
        results[58] = (config.solution.physical_properties.wind[1][index_a] + (config.solution.physical_properties.wind[1][index_b] - config.solution.physical_properties.wind[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[59] = (config.solution.physical_properties.wind[2][index_a] + (config.solution.physical_properties.wind[2][index_b] - config.solution.physical_properties.wind[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[60] = (config.solution.physical_properties.wind[3][index_a] + (config.solution.physical_properties.wind[3][index_b] - config.solution.physical_properties.wind[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        
        results[61] = (config.solution.physical_properties.cL[index_a] + (config.solution.physical_properties.cL[index_b] - config.solution.physical_properties.cL[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a)) 
        results[62] = (config.solution.physical_properties.cD[index_a] + (config.solution.physical_properties.cD[index_b] - config.solution.physical_properties.cD[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[63] = (config.solution.physical_properties.α[index_a] + (config.solution.physical_properties.α[index_b] - config.solution.physical_properties.α[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[64] = (config.solution.physical_properties.S[index_a] + (config.solution.physical_properties.S[index_b] - config.solution.physical_properties.S[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[65] = (config.solution.performance.mass[index_a] + (config.solution.performance.mass[index_b] - config.solution.performance.mass[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[66] = (config.solution.performance.heat_rate[index_a] + (config.solution.performance.heat_rate[index_b] - config.solution.performance.heat_rate[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[67] = (config.solution.performance.heat_load[index_a] + (config.solution.performance.heat_load[index_b] - config.solution.performance.heat_load[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[68] = (config.solution.performance.T_r[index_a] + (config.solution.performance.T_r[index_b] - config.solution.performance.T_r[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[69] = (config.solution.performance.q[index_a] + (config.solution.performance.q[index_b] - config.solution.performance.q[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[70] = (config.solution.forces.gravity_ii[1][index_a] + (config.solution.forces.gravity_ii[1][index_b] - config.solution.forces.gravity_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[71] = (config.solution.forces.gravity_ii[2][index_a] + (config.solution.forces.gravity_ii[2][index_b] - config.solution.forces.gravity_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[72] = (config.solution.forces.gravity_ii[3][index_a] + (config.solution.forces.gravity_ii[3][index_b] - config.solution.forces.gravity_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[73] = (config.solution.forces.drag_pp[1][index_a] + (config.solution.forces.drag_pp[1][index_b] - config.solution.forces.drag_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[74] = (config.solution.forces.drag_pp[2][index_a] + (config.solution.forces.drag_pp[2][index_b] - config.solution.forces.drag_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[75] = (config.solution.forces.drag_pp[3][index_a] + (config.solution.forces.drag_pp[3][index_b] - config.solution.forces.drag_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[76] = (config.solution.forces.drag_ii[1][index_a] + (config.solution.forces.drag_ii[1][index_b] - config.solution.forces.drag_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[77] = (config.solution.forces.drag_ii[2][index_a] + (config.solution.forces.drag_ii[2][index_b] - config.solution.forces.drag_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[78] = (config.solution.forces.drag_ii[3][index_a] + (config.solution.forces.drag_ii[3][index_b] - config.solution.forces.drag_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[79] = (config.solution.forces.lift_pp[1][index_a] + (config.solution.forces.lift_pp[1][index_b] - config.solution.forces.lift_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[80] = (config.solution.forces.lift_pp[2][index_a] + (config.solution.forces.lift_pp[2][index_b] - config.solution.forces.lift_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[81] = (config.solution.forces.lift_pp[3][index_a] + (config.solution.forces.lift_pp[3][index_b] - config.solution.forces.lift_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[82] = (config.solution.forces.lift_ii[1][index_a] + (config.solution.forces.lift_ii[1][index_b] - config.solution.forces.lift_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[83] = (config.solution.forces.lift_ii[2][index_a] + (config.solution.forces.lift_ii[2][index_b] - config.solution.forces.lift_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[84] = (config.solution.forces.lift_ii[3][index_a] + (config.solution.forces.lift_ii[3][index_b] - config.solution.forces.lift_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[85] = (config.solution.forces.force_ii[1][index_a] + (config.solution.forces.force_ii[1][index_b] - config.solution.forces.force_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[86] = (config.solution.forces.force_ii[2][index_a] + (config.solution.forces.force_ii[2][index_b] - config.solution.forces.force_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[87] = (config.solution.forces.force_ii[3][index_a] + (config.solution.forces.force_ii[3][index_b] - config.solution.forces.force_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[88] = (config.solution.forces.energy[index_a] + (config.solution.forces.energy[index_b] - config.solution.forces.energy[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[89] = (config.solution.simulation.MC_seed[index_a] + (config.solution.simulation.MC_seed[index_b] - config.solution.simulation.MC_seed[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[90] = (config.solution.simulation.drag_passage[index_a] + (config.solution.simulation.drag_passage[index_b] - config.solution.simulation.drag_passage[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        deleteat!(config.solution.orientation.time, range(start=index_b+1, step=1, stop=length(config.solution.orientation.time)))
        deleteat!(config.solution.orientation.year, range(start=index_b+1,step=1,stop=length(config.solution.orientation.year)))
        deleteat!(config.solution.orientation.month, range(start=index_b+1,step=1,stop=length(config.solution.orientation.month)))
        deleteat!(config.solution.orientation.day, range(start=index_b+1,step=1,stop=length(config.solution.orientation.day)))
        deleteat!(config.solution.orientation.hour, range(start=index_b+1,step=1,stop=length(config.solution.orientation.hour)))
        deleteat!(config.solution.orientation.minute, range(start=index_b+1,step=1,stop=length(config.solution.orientation.minute)))
        deleteat!(config.solution.orientation.second, range(start=index_b+1,step=1,stop=length(config.solution.orientation.second)))
        deleteat!(config.solution.orientation.number_of_passage, range(start=index_b+1,step=1,stop=length(config.solution.orientation.number_of_passage)))
        deleteat!(config.solution.orientation.pos_ii[1], range(start=index_b+1,step=1,stop=length(config.solution.orientation.pos_ii[1])))   
        deleteat!(config.solution.orientation.pos_ii[2], range(start=index_b+1,step=1,stop=length(config.solution.orientation.pos_ii[2])))   
        deleteat!(config.solution.orientation.pos_ii[3], range(start=index_b+1,step=1,stop=length(config.solution.orientation.pos_ii[3])))   
        deleteat!(config.solution.orientation.vel_ii[1], range(start=index_b+1,step=1,stop=length(config.solution.orientation.vel_ii[1])))   
        deleteat!(config.solution.orientation.vel_ii[2], range(start=index_b+1,step=1,stop=length(config.solution.orientation.vel_ii[2])))   
        deleteat!(config.solution.orientation.vel_ii[3], range(start=index_b+1,step=1,stop=length(config.solution.orientation.vel_ii[3])))   
        deleteat!(config.solution.orientation.pos_ii_mag, range(start=index_b+1,step=1,stop=length(config.solution.orientation.pos_ii_mag)))
        deleteat!(config.solution.orientation.vel_ii_mag, range(start=index_b+1,step=1,stop=length(config.solution.orientation.vel_ii_mag)))

        deleteat!(config.solution.orientation.pos_pp[1], range(start=index_b+1,step=1,stop=length(config.solution.orientation.pos_pp[1])))   
        deleteat!(config.solution.orientation.pos_pp[2], range(start=index_b+1,step=1,stop=length(config.solution.orientation.pos_pp[2])))   
        deleteat!(config.solution.orientation.pos_pp[3], range(start=index_b+1,step=1,stop=length(config.solution.orientation.pos_pp[3])))   
        deleteat!(config.solution.orientation.pos_pp_mag, range(start=index_b+1,step=1,stop=length(config.solution.orientation.pos_pp_mag)))
        deleteat!(config.solution.orientation.vel_pp[1], range(start=index_b+1,step=1,stop=length(config.solution.orientation.vel_pp[1])))   
        deleteat!(config.solution.orientation.vel_pp[2], range(start=index_b+1,step=1,stop=length(config.solution.orientation.vel_pp[2])))   
        deleteat!(config.solution.orientation.vel_pp[3], range(start=index_b+1,step=1,stop=length(config.solution.orientation.vel_pp[3])))   
        deleteat!(config.solution.orientation.vel_pp_mag, range(start=index_b+1,step=1,stop=length(config.solution.orientation.vel_pp_mag)))

        deleteat!(config.solution.orientation.oe[1], range(start=index_b+1,step=1,stop=length(config.solution.orientation.oe[1])))   
        deleteat!(config.solution.orientation.oe[2], range(start=index_b+1,step=1,stop=length(config.solution.orientation.oe[2])))   
        deleteat!(config.solution.orientation.oe[3], range(start=index_b+1,step=1,stop=length(config.solution.orientation.oe[3])))   
        deleteat!(config.solution.orientation.oe[4], range(start=index_b+1,step=1,stop=length(config.solution.orientation.oe[4])))   
        deleteat!(config.solution.orientation.oe[5], range(start=index_b+1,step=1,stop=length(config.solution.orientation.oe[5])))   
        deleteat!(config.solution.orientation.oe[6], range(start=index_b+1,step=1,stop=length(config.solution.orientation.oe[6])))   

        deleteat!(config.solution.orientation.lat, range(start=index_b+1,step=1,stop=length(config.solution.orientation.lat)))
        deleteat!(config.solution.orientation.lon, range(start=index_b+1,step=1,stop=length(config.solution.orientation.lon)))
        deleteat!(config.solution.orientation.alt, range(start=index_b+1,step=1,stop=length(config.solution.orientation.alt)))
        deleteat!(config.solution.orientation.γ_ii, range(start=index_b+1,step=1,stop=length(config.solution.orientation.γ_ii)))
        deleteat!(config.solution.orientation.γ_pp, range(start=index_b+1,step=1,stop=length(config.solution.orientation.γ_pp)))

        deleteat!(config.solution.orientation.h_ii[1], range(start=index_b+1,step=1,stop=length(config.solution.orientation.h_ii[1])))   
        deleteat!(config.solution.orientation.h_ii[2], range(start=index_b+1,step=1,stop=length(config.solution.orientation.h_ii[2])))   
        deleteat!(config.solution.orientation.h_ii[3], range(start=index_b+1,step=1,stop=length(config.solution.orientation.h_ii[3])))   
        deleteat!(config.solution.orientation.h_pp[1], range(start=index_b+1,step=1,stop=length(config.solution.orientation.h_pp[1])))   
        deleteat!(config.solution.orientation.h_pp[2], range(start=index_b+1,step=1,stop=length(config.solution.orientation.h_pp[2])))   
        deleteat!(config.solution.orientation.h_pp[3], range(start=index_b+1,step=1,stop=length(config.solution.orientation.h_pp[3])))   
        deleteat!(config.solution.orientation.h_ii_mag, range(start=index_b+1,step=1,stop=length(config.solution.orientation.h_ii_mag)))
        deleteat!(config.solution.orientation.h_pp_mag, range(start=index_b+1,step=1,stop=length(config.solution.orientation.h_pp_mag)))

        deleteat!(config.solution.orientation.uD[1], range(start=index_b+1,step=1,stop=length(config.solution.orientation.uD[1])))   
        deleteat!(config.solution.orientation.uD[2], range(start=index_b+1,step=1,stop=length(config.solution.orientation.uD[2])))   
        deleteat!(config.solution.orientation.uD[3], range(start=index_b+1,step=1,stop=length(config.solution.orientation.uD[3])))   
        deleteat!(config.solution.orientation.uE[1], range(start=index_b+1,step=1,stop=length(config.solution.orientation.uE[1])))   
        deleteat!(config.solution.orientation.uE[2], range(start=index_b+1,step=1,stop=length(config.solution.orientation.uE[2])))   
        deleteat!(config.solution.orientation.uE[3], range(start=index_b+1,step=1,stop=length(config.solution.orientation.uE[3])))   
        deleteat!(config.solution.orientation.uN[1], range(start=index_b+1,step=1,stop=length(config.solution.orientation.uN[1])))   
        deleteat!(config.solution.orientation.uN[2], range(start=index_b+1,step=1,stop=length(config.solution.orientation.uN[2])))   
        deleteat!(config.solution.orientation.uN[3], range(start=index_b+1,step=1,stop=length(config.solution.orientation.uN[3])))   
        deleteat!(config.solution.orientation.vN, range(start=index_b+1,step=1,stop=length(config.solution.orientation.vN)))
        deleteat!(config.solution.orientation.vE, range(start=index_b+1,step=1,stop=length(config.solution.orientation.vE)))
        deleteat!(config.solution.orientation.azi_pp, range(start=index_b+1,step=1,stop=length(config.solution.orientation.azi_pp)))

        deleteat!(config.solution.physical_properties.ρ, range(start=index_b+1,step=1,stop=length(config.solution.physical_properties.ρ)))
        deleteat!(config.solution.physical_properties.T, range(start=index_b+1,step=1,stop=length(config.solution.physical_properties.T)))
        deleteat!(config.solution.physical_properties.p, range(start=index_b+1,step=1,stop=length(config.solution.physical_properties.p)))
        deleteat!(config.solution.physical_properties.wind[1], range(start=index_b+1,step=1,stop=length(config.solution.physical_properties.wind[1])))   
        deleteat!(config.solution.physical_properties.wind[2], range(start=index_b+1,step=1,stop=length(config.solution.physical_properties.wind[2])))   
        deleteat!(config.solution.physical_properties.wind[3], range(start=index_b+1,step=1,stop=length(config.solution.physical_properties.wind[3])))   
        deleteat!(config.solution.physical_properties.cL, range(start=index_b+1,step=1,stop=length(config.solution.physical_properties.cL)))
        deleteat!(config.solution.physical_properties.cD, range(start=index_b+1,step=1,stop=length(config.solution.physical_properties.cD)))
        deleteat!(config.solution.physical_properties.α, range(start=index_b+1,step=1,stop=length(config.solution.physical_properties.α)))
        deleteat!(config.solution.physical_properties.S, range(start=index_b+1,step=1,stop=length(config.solution.physical_properties.S)))

        deleteat!(config.solution.performance.mass, range(start=index_b+1,step=1,stop=length(config.solution.performance.mass)))
        deleteat!(config.solution.performance.heat_rate, range(start=index_b+1,step=1,stop=length(config.solution.performance.heat_rate)))
        deleteat!(config.solution.performance.heat_load, range(start=index_b+1,step=1,stop=length(config.solution.performance.heat_load)))
        deleteat!(config.solution.performance.T_r, range(start=index_b+1,step=1,stop=length(config.solution.performance.T_r)))
        deleteat!(config.solution.performance.q, range(start=index_b+1,step=1,stop=length(config.solution.performance.q)))

        deleteat!(config.solution.forces.gravity_ii[1], range(start=index_b+1,step=1,stop=length(config.solution.forces.gravity_ii[1])))
        deleteat!(config.solution.forces.gravity_ii[2], range(start=index_b+1,step=1,stop=length(config.solution.forces.gravity_ii[2])))
        deleteat!(config.solution.forces.gravity_ii[3], range(start=index_b+1,step=1,stop=length(config.solution.forces.gravity_ii[3])))
        deleteat!(config.solution.forces.drag_pp[1], range(start=index_b+1,step=1,stop=length(config.solution.forces.drag_pp[1])))
        deleteat!(config.solution.forces.drag_pp[2], range(start=index_b+1,step=1,stop=length(config.solution.forces.drag_pp[2])))
        deleteat!(config.solution.forces.drag_pp[3], range(start=index_b+1,step=1,stop=length(config.solution.forces.drag_pp[3])))
        deleteat!(config.solution.forces.drag_ii[1], range(start=index_b+1,step=1,stop=length(config.solution.forces.drag_ii[1])))
        deleteat!(config.solution.forces.drag_ii[2], range(start=index_b+1,step=1,stop=length(config.solution.forces.drag_ii[2])))
        deleteat!(config.solution.forces.drag_ii[3], range(start=index_b+1,step=1,stop=length(config.solution.forces.drag_ii[3])))
        deleteat!(config.solution.forces.lift_pp[1], range(start=index_b+1,step=1,stop=length(config.solution.forces.lift_pp[1])))
        deleteat!(config.solution.forces.lift_pp[2], range(start=index_b+1,step=1,stop=length(config.solution.forces.lift_pp[2])))
        deleteat!(config.solution.forces.lift_pp[3], range(start=index_b+1,step=1,stop=length(config.solution.forces.lift_pp[3])))
        deleteat!(config.solution.forces.lift_ii[1], range(start=index_b+1,step=1,stop=length(config.solution.forces.lift_ii[1])))
        deleteat!(config.solution.forces.lift_ii[2], range(start=index_b+1,step=1,stop=length(config.solution.forces.lift_ii[2]))) 
        deleteat!(config.solution.forces.lift_ii[3], range(start=index_b+1,step=1,stop=length(config.solution.forces.lift_ii[3]))) 
        deleteat!(config.solution.forces.force_ii[1], range(start=index_b+1,step=1,stop=length(config.solution.forces.force_ii[1])))
        deleteat!(config.solution.forces.force_ii[2], range(start=index_b+1,step=1,stop=length(config.solution.forces.force_ii[2])))
        deleteat!(config.solution.forces.force_ii[3], range(start=index_b+1,step=1,stop=length(config.solution.forces.force_ii[3])))   
        deleteat!(config.solution.forces.energy, range(start=index_b+1,step=1,stop=length(config.solution.forces.energy)))

        deleteat!(config.solution.simulation.MC_seed, range(start=index_b+1,step=1,stop=length(config.solution.simulation.MC_seed)))
        deleteat!(config.solution.simulation.drag_passage, range(start=index_b+1,step=1,stop=length(config.solution.simulation.drag_passage)))

        append!(config.solution.orientation.time, t_in_new)
        append!(config.solution.orientation.year, results[1])
        append!(config.solution.orientation.month, results[2])
        append!(config.solution.orientation.day, results[3])
        append!(config.solution.orientation.hour, results[4])
        append!(config.solution.orientation.minute, results[5])
        append!(config.solution.orientation.second, results[6])
        append!(config.solution.orientation.number_of_passage, results[7])
        append!(config.solution.orientation.pos_ii[1], results[8])
        append!(config.solution.orientation.pos_ii[2], results[9])
        append!(config.solution.orientation.pos_ii[3], results[10])
        append!(config.solution.orientation.vel_ii[1], results[11])
        append!(config.solution.orientation.vel_ii[2], results[12])
        append!(config.solution.orientation.vel_ii[3], results[13])
        append!(config.solution.orientation.pos_ii_mag, results[14])
        append!(config.solution.orientation.vel_ii_mag, results[15])

        append!(config.solution.orientation.pos_pp[1], results[16])
        append!(config.solution.orientation.pos_pp[2], results[17])
        append!(config.solution.orientation.pos_pp[3], results[18])
        append!(config.solution.orientation.pos_pp_mag, results[19])
        append!(config.solution.orientation.vel_pp[1], results[20])
        append!(config.solution.orientation.vel_pp[2], results[21])
        append!(config.solution.orientation.vel_pp[3], results[22])
        append!(config.solution.orientation.vel_pp_mag, results[23])

        append!(config.solution.orientation.oe[1], results[24])
        append!(config.solution.orientation.oe[2], results[25])
        append!(config.solution.orientation.oe[3], results[26])
        append!(config.solution.orientation.oe[4], results[27])
        append!(config.solution.orientation.oe[5], results[28])
        append!(config.solution.orientation.oe[6], results[29])

        append!(config.solution.orientation.lat, results[30])
        append!(config.solution.orientation.lon, results[31])
        append!(config.solution.orientation.alt, results[32])
        append!(config.solution.orientation.γ_ii, results[33])
        append!(config.solution.orientation.γ_pp, results[34])

        append!(config.solution.orientation.h_ii[1], results[35])
        append!(config.solution.orientation.h_ii[2], results[36])
        append!(config.solution.orientation.h_ii[3], results[37])
        append!(config.solution.orientation.h_pp[1], results[38])
        append!(config.solution.orientation.h_pp[2], results[39])
        append!(config.solution.orientation.h_pp[3], results[40])
        append!(config.solution.orientation.h_ii_mag, results[41])
        append!(config.solution.orientation.h_pp_mag, results[42])

        append!(config.solution.orientation.uD[1], results[43])
        append!(config.solution.orientation.uD[2], results[44])
        append!(config.solution.orientation.uD[3], results[45])
        append!(config.solution.orientation.uE[1], results[46])
        append!(config.solution.orientation.uE[2], results[47])
        append!(config.solution.orientation.uE[3], results[48])
        append!(config.solution.orientation.uN[1], results[49])
        append!(config.solution.orientation.uN[2], results[50])
        append!(config.solution.orientation.uN[3], results[51])
        append!(config.solution.orientation.vN, results[52])
        append!(config.solution.orientation.vE, results[53])
        append!(config.solution.orientation.azi_pp, results[54])

        append!(config.solution.physical_properties.ρ, results[55])
        append!(config.solution.physical_properties.T, results[56])
        append!(config.solution.physical_properties.p, results[57])
        append!(config.solution.physical_properties.wind[1], results[58])
        append!(config.solution.physical_properties.wind[2], results[59])
        append!(config.solution.physical_properties.wind[3], results[60])
        append!(config.solution.physical_properties.cL, results[61])
        append!(config.solution.physical_properties.cD, results[62])
        append!(config.solution.physical_properties.α, results[63])
        append!(config.solution.physical_properties.S, results[64])

        append!(config.solution.performance.mass, results[65])
        append!(config.solution.performance.heat_rate, results[66])
        append!(config.solution.performance.heat_load, results[67])
        append!(config.solution.performance.T_r, results[68])
        append!(config.solution.performance.q, results[69])

        append!(config.solution.forces.gravity_ii[1], results[70])
        append!(config.solution.forces.gravity_ii[2], results[71])
        append!(config.solution.forces.gravity_ii[3], results[72])
        append!(config.solution.forces.drag_pp[1], results[73])
        append!(config.solution.forces.drag_pp[2], results[74])
        append!(config.solution.forces.drag_pp[3], results[75])
        append!(config.solution.forces.drag_ii[1], results[76])
        append!(config.solution.forces.drag_ii[2], results[77])
        append!(config.solution.forces.drag_ii[3], results[78])
        append!(config.solution.forces.lift_pp[1], results[79])
        append!(config.solution.forces.lift_pp[2], results[80])
        append!(config.solution.forces.lift_pp[3], results[81])
        append!(config.solution.forces.lift_ii[1], results[82])
        append!(config.solution.forces.lift_ii[2], results[83])
        append!(config.solution.forces.lift_ii[3], results[84])
        append!(config.solution.forces.force_ii[1], results[85])
        append!(config.solution.forces.force_ii[2], results[86])
        append!(config.solution.forces.force_ii[3], results[87])
        append!(config.solution.forces.energy, results[88])

        append!(config.solution.simulation.MC_seed, results[89])
        append!(config.solution.simulation.drag_passage, results[90])

        return initial_state
    end
end