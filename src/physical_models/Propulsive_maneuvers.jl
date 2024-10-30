include("../config.jl")

using Roots

function propulsion_ic_calcs(m, args, initial_state)
    """

    """

    Δv = args.Δv * (-cos(args.phi))

    if args.print_res
        println("LOWER MANEUVER!")
    else
        println("UPPER MANEUVER!")
    end

    v_exausted = m.engine.Isp * m.engine.g_e

    if len(solution.performance.mass) == 0
        m_i = m.body.mass
    else
        m_i = config.solution.performance.mass[end]
    end

    m_f = m_i / exp(Δv / v_exausted)

    Δm = m_i - m_f

    T = m.engine.T

    Δt = Δm * v_exausted / T

    Δt_half = Δt / 2

    vi = initial_state.vi
    e = initial_state.e
    a = initial_state.a

    # check in initial condition is not in drag pass 300 km
    r = m.planet.Rp_e + 300e3
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
        args.Δv = Δv / (-cos(args.ϕ))
        println("-- Thrust Maximum Time Exceeded - Thrust Time and Δv adjusted - NEW Δv = ", args.Δv)
    end

    if length(solution.performance.mass) == 0
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
        t_fin = solution.orientation.time[end]
        t_in_new = t_fin - Δt_half

        dist_list = [abs(t_in_new - p) for p in solution.orientation.time]
        temp = sorted(set(dist_list))

        index_a = firstindex(item -> item == temp[0], dist_list)
        index_b = index_a - 1

        dist_a = dist_list[index_a]
        dist_b = dist_list[index_b]

        t_in_new = (solution.orientation.time[index_a] + (solution.orientation.time[index_b] - solution.orientation.time[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results = [0] * 90
        results[0] = (solution.orientation.year[index_a] + (solution.orientation.year[index_b] - solution.orientation.year[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[1] = (solution.orientation.month[index_a] + (solution.orientation.month[index_b] - solution.orientation.month[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[2] = (solution.orientation.day[index_a] + (solution.orientation.day[index_b] - solution.orientation.day[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[3] = (solution.orientation.hour[index_a] + (solution.orientation.hour[index_b] - solution.orientation.hour[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[4] = (solution.orientation.minute[index_a] + (solution.orientation.minute[index_b] - solution.orientation.minute[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[5] = (solution.orientation.second[index_a] + (solution.orientation.second[index_b] - solution.orientation.second[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        
        results[6] = (solution.orientation.number_of_passage[index_a] + (solution.orientation.number_of_passage[index_b] - solution.orientation.number_of_passage[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        
        results[7] = (solution.orientation.pos_ii[1][index_a] + (solution.orientation.pos_ii[1][index_b] - solution.orientation.pos_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[8] = (solution.orientation.pos_ii[2][index_a] + (solution.orientation.pos_ii[2][index_b] - solution.orientation.pos_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[9] = (solution.orientation.pos_ii[3][index_a] + (solution.orientation.pos_ii[3][index_b] - solution.orientation.pos_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[10] = (solution.orientation.vel_ii[1][index_a] + (solution.orientation.vel_ii[1][index_b] - solution.orientation.vel_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[11] = (solution.orientation.vel_ii[2][index_a] + (solution.orientation.vel_ii[2][index_b] - solution.orientation.vel_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[12] = (solution.orientation.vel_ii[3][index_a] + (solution.orientation.vel_ii[3][index_b] - solution.orientation.vel_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[13] = (solution.orientation.pos_ii_mag[index_a] + (solution.orientation.pos_ii_mag[index_b] - solution.orientation.pos_ii_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[14] = (solution.orientation.vel_ii_mag[index_a] + (solution.orientation.vel_ii_mag[index_b] - solution.orientation.vel_ii_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[15] = (solution.orientation.pos_pp[1][index_a] + (solution.orientation.pos_pp[1][index_b] - solution.orientation.pos_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[16] = (solution.orientation.pos_pp[2][index_a] + (solution.orientation.pos_pp[2][index_b] - solution.orientation.pos_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[17] = (solution.orientation.pos_pp[3][index_a] + (solution.orientation.pos_pp[3][index_b] - solution.orientation.pos_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[18] = (solution.orientation.pos_pp_mag[index_a] + (solution.orientation.pos_pp_mag[index_b] - solution.orientation.pos_pp_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[19] = (solution.orientation.vel_pp[1][index_a] + (solution.orientation.vel_pp[1][index_b] - solution.orientation.vel_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[20] = (solution.orientation.vel_pp[2][index_a] + (solution.orientation.vel_pp[2][index_b] - solution.orientation.vel_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[21] = (solution.orientation.vel_pp[3][index_a] + (solution.orientation.vel_pp[3][index_b] - solution.orientation.vel_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[22] = (solution.orientation.vel_pp_mag[index_a] + (solution.orientation.vel_pp_mag[index_b] - solution.orientation.vel_pp_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[23] = (solution.orientation.oe[1][index_a] + (solution.orientation.oe[1][index_b] - solution.orientation.oe[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[24] = (solution.orientation.oe[2][index_a] + (solution.orientation.oe[2][index_b] - solution.orientation.oe[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[25] = (solution.orientation.oe[3][index_a] + (solution.orientation.oe[3][index_b] - solution.orientation.oe[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[26] = (solution.orientation.oe[4][index_a] + (solution.orientation.oe[4][index_b] - solution.orientation.oe[4][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[27] = (solution.orientation.oe[5][index_a] + (solution.orientation.oe[5][index_b] - solution.orientation.oe[5][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[28] = (solution.orientation.oe[6][index_a] + (solution.orientation.oe[6][index_b] - solution.orientation.oe[6][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[29] = (solution.orientation.lat[index_a] + (solution.orientation.lat[index_b] - solution.orientation.lat[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[30] = (solution.orientation.lon[index_a] + (solution.orientation.lon[index_b] - solution.orientation.lon[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[31] = (solution.orientation.alt[index_a] + (solution.orientation.alt[index_b] - solution.orientation.alt[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[32] = (solution.orientation.γ_ii[index_a] + (solution.orientation.γ_ii[index_b] - solution.orientation.γ_ii[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[33] = (solution.orientation.γ_pp[index_a] + (solution.orientation.γ_pp[index_b] - solution.orientation.γ_pp[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[34] = (solution.orientation.h_ii[1][index_a] + (solution.orientation.h_ii[1][index_b] - solution.orientation.h_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[35] = (solution.orientation.h_ii[2][index_a] + (solution.orientation.h_ii[2][index_b] - solution.orientation.h_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[36] = (solution.orientation.h_ii[3][index_a] + (solution.orientation.h_ii[3][index_b] - solution.orientation.h_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[37] = (solution.orientation.h_pp[1][index_a] + (solution.orientation.h_pp[1][index_b] - solution.orientation.h_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[38] = (solution.orientation.h_pp[2][index_a] + (solution.orientation.h_pp[2][index_b] - solution.orientation.h_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[39] = (solution.orientation.h_pp[3][index_a] + (solution.orientation.h_pp[3][index_b] - solution.orientation.h_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[40] = (solution.orientation.h_ii_mag[index_a] + (solution.orientation.h_ii_mag[index_b] - solution.orientation.h_ii_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[41] = (solution.orientation.h_pp_mag[index_a] + (solution.orientation.h_pp_mag[index_b] - solution.orientation.h_pp_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[42] = (solution.orientation.uD[1][index_a] + (solution.orientation.uD[1][index_b] - solution.orientation.uD[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[43] = (solution.orientation.uD[2][index_a] + (solution.orientation.uD[2][index_b] - solution.orientation.uD[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[44] = (solution.orientation.uD[3][index_a] + (solution.orientation.uD[3][index_b] - solution.orientation.uD[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[45] = (solution.orientation.uE[1][index_a] + (solution.orientation.uE[1][index_b] - solution.orientation.uE[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[46] = (solution.orientation.uE[2][index_a] + (solution.orientation.uE[2][index_b] - solution.orientation.uE[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[47] = (solution.orientation.uE[3][index_a] + (solution.orientation.uE[3][index_b] - solution.orientation.uE[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[48] = (solution.orientation.uN[1][index_a] + (solution.orientation.uN[1][index_b] - solution.orientation.uN[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[49] = (solution.orientation.uN[2][index_a] + (solution.orientation.uN[2][index_b] - solution.orientation.uN[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[50] = (solution.orientation.uN[3][index_a] + (solution.orientation.uN[3][index_b] - solution.orientation.uN[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[51] = (solution.orientation.vN[index_a] + (solution.orientation.vN[index_b] - solution.orientation.vN[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[52] = (solution.orientation.vE[index_a] + (solution.orientation.vE[index_b] - solution.orientation.vE[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[53] = (solution.orientation.azi_pp[index_a] + (solution.orientation.azi_pp[index_b] - solution.orientation.azi_pp[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[54] = (solution.physical_properties.ρ[index_a] + (solution.physical_properties.ρ[index_b] - solution.physical_properties.ρ[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[55] = (solution.physical_properties.T[index_a] + (solution.physical_properties.T[index_b] - solution.physical_properties.T[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[56] = (solution.physical_properties.p[index_a] + (solution.physical_properties.p[index_b] - solution.physical_properties.p[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        
        results[57] = (solution.physical_properties.wind[1][index_a] + (solution.physical_properties.wind[1][index_b] - solution.physical_properties.wind[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[58] = (solution.physical_properties.wind[2][index_a] + (solution.physical_properties.wind[2][index_b] - solution.physical_properties.wind[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[59] = (solution.physical_properties.wind[3][index_a] + (solution.physical_properties.wind[3][index_b] - solution.physical_properties.wind[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        
        results[60] = (solution.physical_properties.cL[index_a] + (solution.physical_properties.cL[index_b] - solution.physical_properties.cL[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a)) 
        results[61] = (solution.physical_properties.cD[index_a] + (solution.physical_properties.cD[index_b] - solution.physical_properties.cD[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[62] = (solution.physical_properties.α[index_a] + (solution.physical_properties.α[index_b] - solution.physical_properties.α[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[63] = (solution.physical_properties.S[index_a] + (solution.physical_properties.S[index_b] - solution.physical_properties.S[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[64] = (solution.performance.mass[index_a] + (solution.performance.mass[index_b] - solution.performance.mass[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[65] = (solution.performance.heat_rate[index_a] + (solution.performance.heat_rate[index_b] - solution.performance.heat_rate[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[66] = (solution.performance.heat_load[index_a] + (solution.performance.heat_load[index_b] - solution.performance.heat_load[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[67] = (solution.performance.T_r[index_a] + (solution.performance.T_r[index_b] - solution.performance.T_r[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[68] = (solution.performance.q[index_a] + (solution.performance.q[index_b] - solution.performance.q[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[69] = (solution.forces.gravity_ii[1][index_a] + (solution.forces.gravity_ii[1][index_b] - solution.forces.gravity_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[70] = (solution.forces.gravity_ii[2][index_a] + (solution.forces.gravity_ii[2][index_b] - solution.forces.gravity_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[71] = (solution.forces.gravity_ii[3][index_a] + (solution.forces.gravity_ii[3][index_b] - solution.forces.gravity_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[72] = (solution.forces.drag_pp[1][index_a] + (solution.forces.drag_pp[1][index_b] - solution.forces.drag_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[73] = (solution.forces.drag_pp[2][index_a] + (solution.forces.drag_pp[2][index_b] - solution.forces.drag_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[74] = (solution.forces.drag_pp[3][index_a] + (solution.forces.drag_pp[3][index_b] - solution.forces.drag_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[75] = (solution.forces.drag_ii[1][index_a] + (solution.forces.drag_ii[1][index_b] - solution.forces.drag_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[76] = (solution.forces.drag_ii[2][index_a] + (solution.forces.drag_ii[2][index_b] - solution.forces.drag_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[77] = (solution.forces.drag_ii[3][index_a] + (solution.forces.drag_ii[3][index_b] - solution.forces.drag_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[78] = (solution.forces.lift_pp[1][index_a] + (solution.forces.lift_pp[1][index_b] - solution.forces.lift_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[79] = (solution.forces.lift_pp[2][index_a] + (solution.forces.lift_pp[2][index_b] - solution.forces.lift_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[80] = (solution.forces.lift_pp[3][index_a] + (solution.forces.lift_pp[3][index_b] - solution.forces.lift_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[81] = (solution.forces.lift_ii[1][index_a] + (solution.forces.lift_ii[1][index_b] - solution.forces.lift_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[82] = (solution.forces.lift_ii[2][index_a] + (solution.forces.lift_ii[2][index_b] - solution.forces.lift_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[83] = (solution.forces.lift_ii[3][index_a] + (solution.forces.lift_ii[3][index_b] - solution.forces.lift_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[84] = (solution.forces.force_ii[1][index_a] + (solution.forces.force_ii[1][index_b] - solution.forces.force_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[85] = (solution.forces.force_ii[2][index_a] + (solution.forces.force_ii[2][index_b] - solution.forces.force_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[86] = (solution.forces.force_ii[3][index_a] + (solution.forces.force_ii[3][index_b] - solution.forces.force_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[87] = (solution.forces.energy[index_a] + (solution.forces.energy[index_b] - solution.forces.energy[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[88] = (solution.simulation.MC_seed[index_a] + (solution.simulation.MC_seed[index_b] - solution.simulation.MC_seed[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[89] = (solution.simulation.drag_passage[index_a] + (solution.simulation.drag_passage[index_b] - solution.simulation.drag_passage[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        solution.orientation.time[index_b+1:end] = nothing
        solution.orientation.year[index_b+1:end] = nothing
        solution.orientation.month[index_b+1:end] = nothing
        solution.orientation.day[index_b+1:end] = nothing
        solution.orientation.hour[index_b+1:end] = nothing
        solution.orientation.minute[index_b+1:end] = nothing
        solution.orientation.second[index_b+1:end] = nothing
        solution.orientation.number_of_passage[index_b+1:end] = nothing
        solution.orientation.pos_ii[1][index_b+1:end] = nothing
        solution.orientation.pos_ii[2][index_b+1:end] = nothing
        solution.orientation.pos_ii[3][index_b+1:end] = nothing
        solution.orientation.vel_ii[1][index_b+1:end] = nothing
        solution.orientation.vel_ii[2][index_b+1:end] = nothing
        solution.orientation.vel_ii[3][index_b+1:end] = nothing
        solution.orientation.pos_ii_mag[index_b+1:end] = nothing
        solution.orientation.vel_ii_mag[index_b+1:end] = nothing

        solution.orientation.pos_pp[1][index_b+1:end] = nothing
        solution.orientation.pos_pp[2][index_b+1:end] = nothing
        solution.orientation.pos_pp[3][index_b+1:end] = nothing
        solution.orientation.pos_pp_mag[index_b+1:end] = nothing
        solution.orientation.vel_pp[1][index_b+1:end] = nothing
        solution.orientation.vel_pp[2][index_b+1:end] = nothing
        solution.orientation.vel_pp[3][index_b+1:end] = nothing
        solution.orientation.vel_pp_mag[index_b+1:end] = nothing

        solution.orientation.oe[1][index_b+1:end] = nothing
        solution.orientation.oe[2][index_b+1:end] = nothing
        solution.orientation.oe[3][index_b+1:end] = nothing
        solution.orientation.oe[4][index_b+1:end] = nothing
        solution.orientation.oe[5][index_b+1:end] = nothing
        solution.orientation.oe[6][index_b+1:end] = nothing

        solution.orientation.lat[index_b+1:end] = nothing
        solution.orientation.lon[index_b+1:end] = nothing
        solution.orientation.alt[index_b+1:end] = nothing
        solution.orientation.γ_ii[index_b+1:end] = nothing
        solution.orientation.γ_pp[index_b+1:end] = nothing

        solution.orientation.h_ii[1][index_b+1:end] = nothing
        solution.orientation.h_ii[2][index_b+1:end] = nothing
        solution.orientation.h_ii[3][index_b+1:end] = nothing
        solution.orientation.h_pp[1][index_b+1:end] = nothing
        solution.orientation.h_pp[2][index_b+1:end] = nothing
        solution.orientation.h_pp[3][index_b+1:end] = nothing
        solution.orientation.h_ii_mag[index_b+1:end] = nothing
        solution.orientation.h_pp_mag[index_b+1:end] = nothing

        solution.orientation.uD[1][index_b+1:end] = nothing
        solution.orientation.uD[2][index_b+1:end] = nothing
        solution.orientation.uD[3][index_b+1:end] = nothing
        solution.orientation.uE[1][index_b+1:end] = nothing
        solution.orientation.uE[2][index_b+1:end] = nothing
        solution.orientation.uE[3][index_b+1:end] = nothing
        solution.orientation.uN[1][index_b+1:end] = nothing
        solution.orientation.uN[2][index_b+1:end] = nothing
        solution.orientation.uN[3][index_b+1:end] = nothing
        solution.orientation.vN[index_b+1:end] = nothing
        solution.orientation.vE[index_b+1:end] = nothing
        solution.orientation.azi_pp[index_b+1:end] = nothing

        solution.physical_properties.ρ[index_b+1:end] = nothing
        solution.physical_properties.T[index_b+1:end] = nothing
        solution.physical_properties.p[index_b+1:end] = nothing
        solution.physical_properties.wind[1][index_b+1:end] = nothing
        solution.physical_properties.wind[2][index_b+1:end] = nothing
        solution.physical_properties.wind[3][index_b+1:end] = nothing
        solution.physical_properties.cL[index_b+1:end] = nothing
        solution.physical_properties.cD[index_b+1:end] = nothing
        solution.physical_properties.α[index_b+1:end] = nothing
        solution.physical_properties.S[index_b+1:end] = nothing

        solution.performance.mass[index_b+1:end] = nothing
        solution.performance.heat_rate[index_b+1:end] = nothing
        solution.performance.heat_load[index_b+1:end] = nothing
        solution.performance.T_r[index_b+1:end] = nothing
        solution.performance.q[index_b+1:end] = nothing

        solution.forces.gravity_ii[1][index_b+1:end] = nothing
        solution.forces.gravity_ii[2][index_b+1:end] = nothing
        solution.forces.gravity_ii[3][index_b+1:end] = nothing
        solution.forces.drag_pp[1][index_b+1:end] = nothing
        solution.forces.drag_pp[2][index_b+1:end] = nothing
        solution.forces.drag_pp[3][index_b+1:end] = nothing
        solution.forces.drag_ii[1][index_b+1:end] = nothing
        solution.forces.drag_ii[2][index_b+1:end] = nothing
        solution.forces.drag_ii[3][index_b+1:end] = nothing
        solution.forces.lift_pp[1][index_b+1:end] = nothing
        solution.forces.lift_pp[2][index_b+1:end] = nothing
        solution.forces.lift_pp[3][index_b+1:end] = nothing
        solution.forces.lift_ii[1][index_b+1:end] = nothing
        solution.forces.lift_ii[2][index_b+1:end] = nothing
        solution.forces.lift_ii[3][index_b+1:end] = nothing
        solution.forces.force_ii[1][index_b+1:end] = nothing
        solution.forces.force_ii[2][index_b+1:end] = nothing
        solution.forces.force_ii[3][index_b+1:end] = nothing
        solution.forces.energy[index_b+1:end] = nothing

        solution.simulation.MC_seed[index_b+1:end] = nothing
        solution.simulation.drag_passage[index_b+1:end] = nothing

        append!(solution.orientation.time, t_in_new)
        append!(solution.orientation.year, results[0])
        append!(solution.orientation.month, results[1])
        append!(solution.orientation.day, results[2])
        append!(solution.orientation.hour, results[3])
        append!(solution.orientation.minute, results[4])
        append!(solution.orientation.second, results[5])
        append!(solution.orientation.number_of_passage, results[6])
        append!(solution.orientation.pos_ii[1], results[7])
        append!(solution.orientation.pos_ii[2], results[8])
        append!(solution.orientation.pos_ii[3], results[9])
        append!(solution.orientation.vel_ii[1], results[10])
        append!(solution.orientation.vel_ii[2], results[11])
        append!(solution.orientation.vel_ii[3], results[12])
        append!(solution.orientation.pos_ii_mag, results[13])
        append!(solution.orientation.vel_ii_mag, results[14])

        append!(solution.orientation.pos_pp[1], results[15])
        append!(solution.orientation.pos_pp[2], results[16])
        append!(solution.orientation.pos_pp[3], results[17])
        append!(solution.orientation.pos_pp_mag, results[18])
        append!(solution.orientation.vel_pp[1], results[19])
        append!(solution.orientation.vel_pp[2], results[20])
        append!(solution.orientation.vel_pp[3], results[21])
        append!(solution.orientation.vel_pp_mag, results[22])

        append!(solution.orientation.oe[1], results[23])
        append!(solution.orientation.oe[2], results[24])
        append!(solution.orientation.oe[3], results[25])
        append!(solution.orientation.oe[4], results[26])
        append!(solution.orientation.oe[5], results[27])
        append!(solution.orientation.oe[6], results[28])

        append!(solution.orientation.lat, results[29])
        append!(solution.orientation.lon, results[30])
        append!(solution.orientation.alt, results[31])
        append!(solution.orientation.γ_ii, results[32])
        append!(solution.orientation.γ_pp, results[33])

        append!(solution.orientation.h_ii[1], results[34])
        append!(solution.orientation.h_ii[2], results[35])
        append!(solution.orientation.h_ii[3], results[36])
        append!(solution.orientation.h_pp[1], results[37])
        append!(solution.orientation.h_pp[2], results[38])
        append!(solution.orientation.h_pp[3], results[39])
        append!(solution.orientation.h_ii_mag, results[40])
        append!(solution.orientation.h_pp_mag, results[41])

        append!(solution.orientation.uD[1], results[42])
        append!(solution.orientation.uD[2], results[43])
        append!(solution.orientation.uD[3], results[44])
        append!(solution.orientation.uE[1], results[45])
        append!(solution.orientation.uE[2], results[46])
        append!(solution.orientation.uE[3], results[47])
        append!(solution.orientation.uN[1], results[48])
        append!(solution.orientation.uN[2], results[49])
        append!(solution.orientation.uN[3], results[50])
        append!(solution.orientation.vN, results[51])
        append!(solution.orientation.vE, results[52])
        append!(solution.orientation.azi_pp, results[53])

        append!(solution.physical_properties.ρ, results[54])
        append!(solution.physical_properties.T, results[55])
        append!(solution.physical_properties.p, results[56])
        append!(solution.physical_properties.wind[1], results[57])
        append!(solution.physical_properties.wind[2], results[58])
        append!(solution.physical_properties.wind[3], results[59])
        append!(solution.physical_properties.cL, results[60])
        append!(solution.physical_properties.cD, results[61])
        append!(solution.physical_properties.α, results[62])
        append!(solution.physical_properties.S, results[63])

        append!(solution.performance.mass, results[64])
        append!(solution.performance.heat_rate, results[65])
        append!(solution.performance.heat_load, results[66])
        append!(solution.performance.T_r, results[67])
        append!(solution.performance.q, results[68])

        append!(solution.forces.gravity_ii[1], results[69])
        append!(solution.forces.gravity_ii[2], results[70])
        append!(solution.forces.gravity_ii[3], results[71])
        append!(solution.forces.drag_pp[1], results[72])
        append!(solution.forces.drag_pp[2], results[73])
        append!(solution.forces.drag_pp[3], results[74])
        append!(solution.forces.drag_ii[1], results[75])
        append!(solution.forces.drag_ii[2], results[76])
        append!(solution.forces.drag_ii[3], results[77])
        append!(solution.forces.lift_pp[1], results[78])
        append!(solution.forces.lift_pp[2], results[79])
        append!(solution.forces.lift_pp[3], results[80])
        append!(solution.forces.lift_ii[1], results[81])
        append!(solution.forces.lift_ii[2], results[82])
        append!(solution.forces.lift_ii[3], results[83])
        append!(solution.forces.force_ii[1], results[84])
        append!(solution.forces.force_ii[2], results[85])
        append!(solution.forces.force_ii[3], results[86])
        append!(solution.forces.energy, results[87])

        append!(solution.simulation.MC_seed, results[88])
        append!(solution.simulation.drag_passage, results[89])

        return initial_state
    end
end