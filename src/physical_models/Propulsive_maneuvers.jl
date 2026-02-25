using Roots

if !isdefined(@__MODULE__, :_legacy_get_propulsive_runtime_state)
    @inline function _legacy_get_propulsive_runtime_state(args=nothing; cnf=nothing, solution=nothing, model=nothing)
        cnf_state = if cnf !== nothing
            cnf
        elseif args isa AbstractDict && haskey(args, :cnf)
            args[:cnf]
        else
            nothing
        end

        solution_state = if solution !== nothing
            solution
        elseif args isa AbstractDict && haskey(args, :solution)
            args[:solution]
        else
            nothing
        end

        model_state = if model !== nothing
            model
        elseif args isa AbstractDict && haskey(args, :model)
            args[:model]
        else
            nothing
        end

        if solution_state === nothing || model_state === nothing
            throw(ArgumentError("Legacy runtime state missing; provide `solution` and `model`."))
        end

        return (cnf=cnf_state, solution=solution_state, model=model_state)
    end
end

function propulsion_ic_calcs(m, args, initial_state)
    """

    """
    runtime = _legacy_get_propulsive_runtime_state(args; model=m)

    Δv = args[:delta_v] * (-cos(args[:phi]))

    if Bool(args[:print_res])
        if round(rad2deg(args[:phi])) == 0
            println("LOWER MANEUVER!")
        elseif round(rad2deg(args[:phi])) == 180
            println("RAISE MANEUVER!")
        end
    end

    v_exausted = m.engines.Isp * m.engines.g_e

    if length(runtime.solution.performance.mass) == 0
        m_i = m.body.dry_mass + m.body.prop_mass
    else
        m_i = runtime.solution.performance.mass[end]
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

    if length(runtime.solution.performance.mass) == 0
        # inverse problem of Kepler
        E_0 = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(vi / 2))
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

        if vi_in < 0
            vi_in += 2*pi
        end

        if vi_in > pi
            δ = vi_in - pi
            vi_in = pi - δ
        end

        initial_state.vi = vi_in

        return initial_state
    else
        t_fin = runtime.solution.orientation.time[end]
        t_in_new = t_fin - Δt_half

        dist_list = [abs(t_in_new - p) for p in runtime.solution.orientation.time]
        temp = sort(collect(Set(dist_list)))

        index_a = findfirst(item -> item == temp[1], dist_list)
        index_b = index_a - 1

        dist_a = dist_list[index_a]
        dist_b = dist_list[index_b]

        t_in_new = (runtime.solution.orientation.time[index_a] + (runtime.solution.orientation.time[index_b] - runtime.solution.orientation.time[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results = ones(runtime.solution.simulation.solution_states)
        results[1] = (runtime.solution.orientation.year[index_a] + (runtime.solution.orientation.year[index_b] - runtime.solution.orientation.year[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[2] = (runtime.solution.orientation.month[index_a] + (runtime.solution.orientation.month[index_b] - runtime.solution.orientation.month[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[3] = (runtime.solution.orientation.day[index_a] + (runtime.solution.orientation.day[index_b] - runtime.solution.orientation.day[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[4] = (runtime.solution.orientation.hour[index_a] + (runtime.solution.orientation.hour[index_b] - runtime.solution.orientation.hour[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[5] = (runtime.solution.orientation.minute[index_a] + (runtime.solution.orientation.minute[index_b] - runtime.solution.orientation.minute[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[6] = (runtime.solution.orientation.second[index_a] + (runtime.solution.orientation.second[index_b] - runtime.solution.orientation.second[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        
        results[7] = (runtime.solution.orientation.number_of_passage[index_a] + (runtime.solution.orientation.number_of_passage[index_b] - runtime.solution.orientation.number_of_passage[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        
        results[8] = (runtime.solution.orientation.pos_ii[1][index_a] + (runtime.solution.orientation.pos_ii[1][index_b] - runtime.solution.orientation.pos_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[9] = (runtime.solution.orientation.pos_ii[2][index_a] + (runtime.solution.orientation.pos_ii[2][index_b] - runtime.solution.orientation.pos_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[10] = (runtime.solution.orientation.pos_ii[3][index_a] + (runtime.solution.orientation.pos_ii[3][index_b] - runtime.solution.orientation.pos_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[11] = (runtime.solution.orientation.vel_ii[1][index_a] + (runtime.solution.orientation.vel_ii[1][index_b] - runtime.solution.orientation.vel_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[12] = (runtime.solution.orientation.vel_ii[2][index_a] + (runtime.solution.orientation.vel_ii[2][index_b] - runtime.solution.orientation.vel_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[13] = (runtime.solution.orientation.vel_ii[3][index_a] + (runtime.solution.orientation.vel_ii[3][index_b] - runtime.solution.orientation.vel_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[14] = (runtime.solution.orientation.pos_ii_mag[index_a] + (runtime.solution.orientation.pos_ii_mag[index_b] - runtime.solution.orientation.pos_ii_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[15] = (runtime.solution.orientation.vel_ii_mag[index_a] + (runtime.solution.orientation.vel_ii_mag[index_b] - runtime.solution.orientation.vel_ii_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[16] = (runtime.solution.orientation.pos_pp[1][index_a] + (runtime.solution.orientation.pos_pp[1][index_b] - runtime.solution.orientation.pos_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[17] = (runtime.solution.orientation.pos_pp[2][index_a] + (runtime.solution.orientation.pos_pp[2][index_b] - runtime.solution.orientation.pos_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[18] = (runtime.solution.orientation.pos_pp[3][index_a] + (runtime.solution.orientation.pos_pp[3][index_b] - runtime.solution.orientation.pos_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[19] = (runtime.solution.orientation.pos_pp_mag[index_a] + (runtime.solution.orientation.pos_pp_mag[index_b] - runtime.solution.orientation.pos_pp_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[20] = (runtime.solution.orientation.vel_pp[1][index_a] + (runtime.solution.orientation.vel_pp[1][index_b] - runtime.solution.orientation.vel_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[21] = (runtime.solution.orientation.vel_pp[2][index_a] + (runtime.solution.orientation.vel_pp[2][index_b] - runtime.solution.orientation.vel_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[22] = (runtime.solution.orientation.vel_pp[3][index_a] + (runtime.solution.orientation.vel_pp[3][index_b] - runtime.solution.orientation.vel_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[23] = (runtime.solution.orientation.vel_pp_mag[index_a] + (runtime.solution.orientation.vel_pp_mag[index_b] - runtime.solution.orientation.vel_pp_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[24] = (runtime.solution.orientation.oe[1][index_a] + (runtime.solution.orientation.oe[1][index_b] - runtime.solution.orientation.oe[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[25] = (runtime.solution.orientation.oe[2][index_a] + (runtime.solution.orientation.oe[2][index_b] - runtime.solution.orientation.oe[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[26] = (runtime.solution.orientation.oe[3][index_a] + (runtime.solution.orientation.oe[3][index_b] - runtime.solution.orientation.oe[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[27] = (runtime.solution.orientation.oe[4][index_a] + (runtime.solution.orientation.oe[4][index_b] - runtime.solution.orientation.oe[4][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[28] = (runtime.solution.orientation.oe[5][index_a] + (runtime.solution.orientation.oe[5][index_b] - runtime.solution.orientation.oe[5][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[29] = (runtime.solution.orientation.oe[6][index_a] + (runtime.solution.orientation.oe[6][index_b] - runtime.solution.orientation.oe[6][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[30] = (runtime.solution.orientation.lat[index_a] + (runtime.solution.orientation.lat[index_b] - runtime.solution.orientation.lat[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[31] = (runtime.solution.orientation.lon[index_a] + (runtime.solution.orientation.lon[index_b] - runtime.solution.orientation.lon[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[32] = (runtime.solution.orientation.alt[index_a] + (runtime.solution.orientation.alt[index_b] - runtime.solution.orientation.alt[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[33] = (runtime.solution.orientation.γ_ii[index_a] + (runtime.solution.orientation.γ_ii[index_b] - runtime.solution.orientation.γ_ii[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[34] = (runtime.solution.orientation.γ_pp[index_a] + (runtime.solution.orientation.γ_pp[index_b] - runtime.solution.orientation.γ_pp[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[35] = (runtime.solution.orientation.h_ii[1][index_a] + (runtime.solution.orientation.h_ii[1][index_b] - runtime.solution.orientation.h_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[36] = (runtime.solution.orientation.h_ii[2][index_a] + (runtime.solution.orientation.h_ii[2][index_b] - runtime.solution.orientation.h_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[37] = (runtime.solution.orientation.h_ii[3][index_a] + (runtime.solution.orientation.h_ii[3][index_b] - runtime.solution.orientation.h_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[38] = (runtime.solution.orientation.h_pp[1][index_a] + (runtime.solution.orientation.h_pp[1][index_b] - runtime.solution.orientation.h_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[39] = (runtime.solution.orientation.h_pp[2][index_a] + (runtime.solution.orientation.h_pp[2][index_b] - runtime.solution.orientation.h_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[40] = (runtime.solution.orientation.h_pp[3][index_a] + (runtime.solution.orientation.h_pp[3][index_b] - runtime.solution.orientation.h_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[41] = (runtime.solution.orientation.h_ii_mag[index_a] + (runtime.solution.orientation.h_ii_mag[index_b] - runtime.solution.orientation.h_ii_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[42] = (runtime.solution.orientation.h_pp_mag[index_a] + (runtime.solution.orientation.h_pp_mag[index_b] - runtime.solution.orientation.h_pp_mag[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[43] = (runtime.solution.orientation.uD[1][index_a] + (runtime.solution.orientation.uD[1][index_b] - runtime.solution.orientation.uD[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[44] = (runtime.solution.orientation.uD[2][index_a] + (runtime.solution.orientation.uD[2][index_b] - runtime.solution.orientation.uD[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[45] = (runtime.solution.orientation.uD[3][index_a] + (runtime.solution.orientation.uD[3][index_b] - runtime.solution.orientation.uD[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[46] = (runtime.solution.orientation.uE[1][index_a] + (runtime.solution.orientation.uE[1][index_b] - runtime.solution.orientation.uE[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[47] = (runtime.solution.orientation.uE[2][index_a] + (runtime.solution.orientation.uE[2][index_b] - runtime.solution.orientation.uE[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[48] = (runtime.solution.orientation.uE[3][index_a] + (runtime.solution.orientation.uE[3][index_b] - runtime.solution.orientation.uE[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[49] = (runtime.solution.orientation.uN[1][index_a] + (runtime.solution.orientation.uN[1][index_b] - runtime.solution.orientation.uN[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[50] = (runtime.solution.orientation.uN[2][index_a] + (runtime.solution.orientation.uN[2][index_b] - runtime.solution.orientation.uN[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[51] = (runtime.solution.orientation.uN[3][index_a] + (runtime.solution.orientation.uN[3][index_b] - runtime.solution.orientation.uN[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[52] = (runtime.solution.orientation.vN[index_a] + (runtime.solution.orientation.vN[index_b] - runtime.solution.orientation.vN[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[53] = (runtime.solution.orientation.vE[index_a] + (runtime.solution.orientation.vE[index_b] - runtime.solution.orientation.vE[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[54] = (runtime.solution.orientation.azi_pp[index_a] + (runtime.solution.orientation.azi_pp[index_b] - runtime.solution.orientation.azi_pp[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[55] = (runtime.solution.physical_properties.ρ[index_a] + (runtime.solution.physical_properties.ρ[index_b] - runtime.solution.physical_properties.ρ[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[56] = (runtime.solution.physical_properties.T[index_a] + (runtime.solution.physical_properties.T[index_b] - runtime.solution.physical_properties.T[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[57] = (runtime.solution.physical_properties.p[index_a] + (runtime.solution.physical_properties.p[index_b] - runtime.solution.physical_properties.p[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        
        results[58] = (runtime.solution.physical_properties.wind[1][index_a] + (runtime.solution.physical_properties.wind[1][index_b] - runtime.solution.physical_properties.wind[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[59] = (runtime.solution.physical_properties.wind[2][index_a] + (runtime.solution.physical_properties.wind[2][index_b] - runtime.solution.physical_properties.wind[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[60] = (runtime.solution.physical_properties.wind[3][index_a] + (runtime.solution.physical_properties.wind[3][index_b] - runtime.solution.physical_properties.wind[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        
        results[61] = (runtime.solution.physical_properties.cL[index_a] + (runtime.solution.physical_properties.cL[index_b] - runtime.solution.physical_properties.cL[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a)) 
        results[62] = (runtime.solution.physical_properties.cD[index_a] + (runtime.solution.physical_properties.cD[index_b] - runtime.solution.physical_properties.cD[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[63] = (runtime.solution.physical_properties.S[index_a] + (runtime.solution.physical_properties.S[index_b] - runtime.solution.physical_properties.S[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[98] = (runtime.solution.physical_properties.α_control[index_a] + (runtime.solution.physical_properties.α_control[index_b] - runtime.solution.physical_properties.α_control[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[99] = (runtime.solution.physical_properties.inertia_tensor[1][index_a] + (runtime.solution.physical_properties.inertia_tensor[1][index_b] - runtime.solution.physical_properties.inertia_tensor[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[100] = (runtime.solution.physical_properties.inertia_tensor[2][index_a] + (runtime.solution.physical_properties.inertia_tensor[2][index_b] - runtime.solution.physical_properties.inertia_tensor[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[101] = (runtime.solution.physical_properties.inertia_tensor[3][index_a] + (runtime.solution.physical_properties.inertia_tensor[3][index_b] - runtime.solution.physical_properties.inertia_tensor[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[102] = (runtime.solution.physical_properties.inertia_tensor[4][index_a] + (runtime.solution.physical_properties.inertia_tensor[4][index_b] - runtime.solution.physical_properties.inertia_tensor[4][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[103] = (runtime.solution.physical_properties.inertia_tensor[5][index_a] + (runtime.solution.physical_properties.inertia_tensor[5][index_b] - runtime.solution.physical_properties.inertia_tensor[5][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[104] = (runtime.solution.physical_properties.inertia_tensor[6][index_a] + (runtime.solution.physical_properties.inertia_tensor[6][index_b] - runtime.solution.physical_properties.inertia_tensor[6][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[105] = (runtime.solution.physical_properties.inertia_tensor[7][index_a] + (runtime.solution.physical_properties.inertia_tensor[7][index_b] - runtime.solution.physical_properties.inertia_tensor[7][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[106] = (runtime.solution.physical_properties.inertia_tensor[8][index_a] + (runtime.solution.physical_properties.inertia_tensor[8][index_b] - runtime.solution.physical_properties.inertia_tensor[8][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[107] = (runtime.solution.physical_properties.inertia_tensor[9][index_a] + (runtime.solution.physical_properties.inertia_tensor[9][index_b] - runtime.solution.physical_properties.inertia_tensor[9][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[108] = (runtime.solution.physical_properties.τ_rw[1][index_a] + (runtime.solution.physical_properties.τ_rw[1][index_b] - runtime.solution.physical_properties.τ_rw[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[109] = (runtime.solution.physical_properties.τ_rw[2][index_a] + (runtime.solution.physical_properties.τ_rw[2][index_b] - runtime.solution.physical_properties.τ_rw[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[110] = (runtime.solution.physical_properties.τ_rw[3][index_a] + (runtime.solution.physical_properties.τ_rw[3][index_b] - runtime.solution.physical_properties.τ_rw[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        n_bodies = length(runtime.model.body.links)
        for i in 1:n_bodies
            results[110 + i] = (runtime.solution.physical_properties.α[i][index_a] + (runtime.solution.physical_properties.α[i][index_b] - runtime.solution.physical_properties.α[i][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
            results[110 + n_bodies + i] = (runtime.solution.physical_properties.β[i][index_a] + (runtime.solution.physical_properties.β[i][index_b] - runtime.solution.physical_properties.β[i][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
            results[110 + 2*n_bodies + i] = (runtime.solution.performance.heat_rate[i][index_a] + (runtime.solution.performance.heat_rate[i][index_b] - runtime.solution.performance.heat_rate[i][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
            results[110 + 3*n_bodies + i] = (runtime.solution.performance.heat_load[i][index_a] + (runtime.solution.performance.heat_load[i][index_b] - runtime.solution.performance.heat_load[i][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        end
        n_reaction_wheels = runtime.model.body.n_reaction_wheels
        for i in 1:n_reaction_wheels
            results[110 + 4*n_bodies + i] = (runtime.solution.physical_properties.rw_h[i][index_a] + (runtime.solution.physical_properties.rw_h[i][index_b] - runtime.solution.physical_properties.rw_h[i][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
            results[110 + 4*n_bodies + n_reaction_wheels + i] = (runtime.solution.physical_properties.rw_τ[i][index_a] + (runtime.solution.physical_properties.rw_τ[i][index_b] - runtime.solution.physical_properties.rw_τ[i][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        end

        n_thrusters = runtime.model.body.n_thrusters
        for i in 1:n_thrusters
            results[110 + 4*n_bodies + 2*n_reaction_wheels + i] = (runtime.solution.physical_properties.thruster_h[i][index_a] + (runtime.solution.physical_properties.thruster_h[i][index_b] - runtime.solution.physical_properties.thruster_h[i][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        end
        # results[98 + 4*n_bodies + 2*n_reaction_wheels + 1] = (runtime.solution.physical_properties.τ_rw[1][index_a] + (runtime.solution.physical_properties.τ_rw[1][index_b] - runtime.solution.physical_properties.τ_rw[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        # results[98 + 4*n_bodies + 2*n_reaction_wheels + 2] = (runtime.solution.physical_properties.τ_rw[2][index_a] + (runtime.solution.physical_properties.τ_rw[2][index_b] - runtime.solution.physical_properties.τ_rw[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        # results[98 + 4*n_bodies + 2*n_reaction_wheels + 3] = (runtime.solution.physical_properties.τ_rw[3][index_a] + (runtime.solution.physical_properties.τ_rw[3][index_b] - runtime.solution.physical_properties.τ_rw[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[64] = (runtime.solution.performance.mass[index_a] + (runtime.solution.performance.mass[index_b] - runtime.solution.performance.mass[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        # results[65] = (runtime.solution.performance.heat_rate[index_a] + (runtime.solution.performance.heat_rate[index_b] - runtime.solution.performance.heat_rate[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        # results[66] = (runtime.solution.performance.heat_load[index_a] + (runtime.solution.performance.heat_load[index_b] - runtime.solution.performance.heat_load[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[65] = (runtime.solution.performance.T_r[index_a] + (runtime.solution.performance.T_r[index_b] - runtime.solution.performance.T_r[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[66] = (runtime.solution.performance.q[index_a] + (runtime.solution.performance.q[index_b] - runtime.solution.performance.q[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[67] = (runtime.solution.forces.gravity_ii[1][index_a] + (runtime.solution.forces.gravity_ii[1][index_b] - runtime.solution.forces.gravity_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[68] = (runtime.solution.forces.gravity_ii[2][index_a] + (runtime.solution.forces.gravity_ii[2][index_b] - runtime.solution.forces.gravity_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[69] = (runtime.solution.forces.gravity_ii[3][index_a] + (runtime.solution.forces.gravity_ii[3][index_b] - runtime.solution.forces.gravity_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[70] = (runtime.solution.forces.drag_pp[1][index_a] + (runtime.solution.forces.drag_pp[1][index_b] - runtime.solution.forces.drag_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[71] = (runtime.solution.forces.drag_pp[2][index_a] + (runtime.solution.forces.drag_pp[2][index_b] - runtime.solution.forces.drag_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[72] = (runtime.solution.forces.drag_pp[3][index_a] + (runtime.solution.forces.drag_pp[3][index_b] - runtime.solution.forces.drag_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[73] = (runtime.solution.forces.drag_ii[1][index_a] + (runtime.solution.forces.drag_ii[1][index_b] - runtime.solution.forces.drag_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[74] = (runtime.solution.forces.drag_ii[2][index_a] + (runtime.solution.forces.drag_ii[2][index_b] - runtime.solution.forces.drag_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[75] = (runtime.solution.forces.drag_ii[3][index_a] + (runtime.solution.forces.drag_ii[3][index_b] - runtime.solution.forces.drag_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[76] = (runtime.solution.forces.lift_pp[1][index_a] + (runtime.solution.forces.lift_pp[1][index_b] - runtime.solution.forces.lift_pp[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[77] = (runtime.solution.forces.lift_pp[2][index_a] + (runtime.solution.forces.lift_pp[2][index_b] - runtime.solution.forces.lift_pp[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[78] = (runtime.solution.forces.lift_pp[3][index_a] + (runtime.solution.forces.lift_pp[3][index_b] - runtime.solution.forces.lift_pp[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[79] = (runtime.solution.forces.lift_ii[1][index_a] + (runtime.solution.forces.lift_ii[1][index_b] - runtime.solution.forces.lift_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[80] = (runtime.solution.forces.lift_ii[2][index_a] + (runtime.solution.forces.lift_ii[2][index_b] - runtime.solution.forces.lift_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[81] = (runtime.solution.forces.lift_ii[3][index_a] + (runtime.solution.forces.lift_ii[3][index_b] - runtime.solution.forces.lift_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[82] = (runtime.solution.forces.force_ii[1][index_a] + (runtime.solution.forces.force_ii[1][index_b] - runtime.solution.forces.force_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[83] = (runtime.solution.forces.force_ii[2][index_a] + (runtime.solution.forces.force_ii[2][index_b] - runtime.solution.forces.force_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[84] = (runtime.solution.forces.force_ii[3][index_a] + (runtime.solution.forces.force_ii[3][index_b] - runtime.solution.forces.force_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[85] = (runtime.solution.forces.τ_ii[1][index_a] + (runtime.solution.forces.τ_ii[1][index_b] - runtime.solution.forces.τ_ii[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[86] = (runtime.solution.forces.τ_ii[2][index_a] + (runtime.solution.forces.τ_ii[2][index_b] - runtime.solution.forces.τ_ii[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[87] = (runtime.solution.forces.τ_ii[3][index_a] + (runtime.solution.forces.τ_ii[3][index_b] - runtime.solution.forces.τ_ii[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[88] = (runtime.solution.forces.energy[index_a] + (runtime.solution.forces.energy[index_b] - runtime.solution.forces.energy[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[89] = (runtime.solution.simulation.MC_seed[index_a] + (runtime.solution.simulation.MC_seed[index_b] - runtime.solution.simulation.MC_seed[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[90] = (runtime.solution.simulation.drag_passage[index_a] + (runtime.solution.simulation.drag_passage[index_b] - runtime.solution.simulation.drag_passage[index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        results[91] = (runtime.solution.orientation.quaternion[1][index_a] + (runtime.solution.orientation.quaternion[1][index_b] - runtime.solution.orientation.quaternion[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[92] = (runtime.solution.orientation.quaternion[2][index_a] + (runtime.solution.orientation.quaternion[2][index_b] - runtime.solution.orientation.quaternion[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[93] = (runtime.solution.orientation.quaternion[3][index_a] + (runtime.solution.orientation.quaternion[3][index_b] - runtime.solution.orientation.quaternion[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[94] = (runtime.solution.orientation.quaternion[4][index_a] + (runtime.solution.orientation.quaternion[4][index_b] - runtime.solution.orientation.quaternion[4][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[95] = (runtime.solution.orientation.ω[1][index_a] + (runtime.solution.orientation.ω[1][index_b] - runtime.solution.orientation.ω[1][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[96] = (runtime.solution.orientation.ω[2][index_a] + (runtime.solution.orientation.ω[2][index_b] - runtime.solution.orientation.ω[2][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))
        results[97] = (runtime.solution.orientation.ω[3][index_a] + (runtime.solution.orientation.ω[3][index_b] - runtime.solution.orientation.ω[3][index_a]) / (dist_b + abs(dist_a)) * abs(dist_a))

        deleteat!(runtime.solution.orientation.time, range(start=index_b+1, step=1, stop=length(runtime.solution.orientation.time)))
        deleteat!(runtime.solution.orientation.year, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.year)))
        deleteat!(runtime.solution.orientation.month, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.month)))
        deleteat!(runtime.solution.orientation.day, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.day)))
        deleteat!(runtime.solution.orientation.hour, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.hour)))
        deleteat!(runtime.solution.orientation.minute, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.minute)))
        deleteat!(runtime.solution.orientation.second, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.second)))
        deleteat!(runtime.solution.orientation.number_of_passage, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.number_of_passage)))
        deleteat!(runtime.solution.orientation.pos_ii[1], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.pos_ii[1])))   
        deleteat!(runtime.solution.orientation.pos_ii[2], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.pos_ii[2])))   
        deleteat!(runtime.solution.orientation.pos_ii[3], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.pos_ii[3])))   
        deleteat!(runtime.solution.orientation.vel_ii[1], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.vel_ii[1])))   
        deleteat!(runtime.solution.orientation.vel_ii[2], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.vel_ii[2])))   
        deleteat!(runtime.solution.orientation.vel_ii[3], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.vel_ii[3])))   
        deleteat!(runtime.solution.orientation.pos_ii_mag, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.pos_ii_mag)))
        deleteat!(runtime.solution.orientation.vel_ii_mag, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.vel_ii_mag)))
        deleteat!(runtime.solution.orientation.quaternion[1], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.quaternion[1])))   
        deleteat!(runtime.solution.orientation.quaternion[2], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.quaternion[2])))   
        deleteat!(runtime.solution.orientation.quaternion[3], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.quaternion[3])))   
        deleteat!(runtime.solution.orientation.quaternion[4], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.quaternion[4])))
        deleteat!(runtime.solution.orientation.ω[1], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.ω[1])))   
        deleteat!(runtime.solution.orientation.ω[2], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.ω[2])))   
        deleteat!(runtime.solution.orientation.ω[3], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.ω[3])))

        deleteat!(runtime.solution.orientation.pos_pp[1], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.pos_pp[1])))   
        deleteat!(runtime.solution.orientation.pos_pp[2], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.pos_pp[2])))   
        deleteat!(runtime.solution.orientation.pos_pp[3], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.pos_pp[3])))   
        deleteat!(runtime.solution.orientation.pos_pp_mag, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.pos_pp_mag)))
        deleteat!(runtime.solution.orientation.vel_pp[1], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.vel_pp[1])))   
        deleteat!(runtime.solution.orientation.vel_pp[2], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.vel_pp[2])))   
        deleteat!(runtime.solution.orientation.vel_pp[3], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.vel_pp[3])))   
        deleteat!(runtime.solution.orientation.vel_pp_mag, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.vel_pp_mag)))

        deleteat!(runtime.solution.orientation.oe[1], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.oe[1])))   
        deleteat!(runtime.solution.orientation.oe[2], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.oe[2])))   
        deleteat!(runtime.solution.orientation.oe[3], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.oe[3])))   
        deleteat!(runtime.solution.orientation.oe[4], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.oe[4])))   
        deleteat!(runtime.solution.orientation.oe[5], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.oe[5])))   
        deleteat!(runtime.solution.orientation.oe[6], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.oe[6])))   

        deleteat!(runtime.solution.orientation.lat, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.lat)))
        deleteat!(runtime.solution.orientation.lon, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.lon)))
        deleteat!(runtime.solution.orientation.alt, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.alt)))
        deleteat!(runtime.solution.orientation.γ_ii, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.γ_ii)))
        deleteat!(runtime.solution.orientation.γ_pp, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.γ_pp)))

        deleteat!(runtime.solution.orientation.h_ii[1], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.h_ii[1])))   
        deleteat!(runtime.solution.orientation.h_ii[2], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.h_ii[2])))   
        deleteat!(runtime.solution.orientation.h_ii[3], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.h_ii[3])))   
        deleteat!(runtime.solution.orientation.h_pp[1], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.h_pp[1])))   
        deleteat!(runtime.solution.orientation.h_pp[2], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.h_pp[2])))   
        deleteat!(runtime.solution.orientation.h_pp[3], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.h_pp[3])))   
        deleteat!(runtime.solution.orientation.h_ii_mag, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.h_ii_mag)))
        deleteat!(runtime.solution.orientation.h_pp_mag, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.h_pp_mag)))

        deleteat!(runtime.solution.orientation.uD[1], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.uD[1])))   
        deleteat!(runtime.solution.orientation.uD[2], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.uD[2])))   
        deleteat!(runtime.solution.orientation.uD[3], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.uD[3])))   
        deleteat!(runtime.solution.orientation.uE[1], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.uE[1])))   
        deleteat!(runtime.solution.orientation.uE[2], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.uE[2])))   
        deleteat!(runtime.solution.orientation.uE[3], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.uE[3])))   
        deleteat!(runtime.solution.orientation.uN[1], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.uN[1])))   
        deleteat!(runtime.solution.orientation.uN[2], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.uN[2])))   
        deleteat!(runtime.solution.orientation.uN[3], range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.uN[3])))   
        deleteat!(runtime.solution.orientation.vN, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.vN)))
        deleteat!(runtime.solution.orientation.vE, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.vE)))
        deleteat!(runtime.solution.orientation.azi_pp, range(start=index_b+1,step=1,stop=length(runtime.solution.orientation.azi_pp)))

        deleteat!(runtime.solution.physical_properties.ρ, range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.ρ)))
        deleteat!(runtime.solution.physical_properties.T, range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.T)))
        deleteat!(runtime.solution.physical_properties.p, range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.p)))
        deleteat!(runtime.solution.physical_properties.wind[1], range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.wind[1])))   
        deleteat!(runtime.solution.physical_properties.wind[2], range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.wind[2])))   
        deleteat!(runtime.solution.physical_properties.wind[3], range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.wind[3])))   
        deleteat!(runtime.solution.physical_properties.cL, range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.cL)))
        deleteat!(runtime.solution.physical_properties.cD, range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.cD)))
        deleteat!(runtime.solution.physical_properties.S, range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.S)))
        deleteat!(runtime.solution.physical_properties.α_control, range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.α_control)))
        
        for i in 1:n_bodies
            deleteat!(runtime.solution.physical_properties.α[i], range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.α[i])))
            deleteat!(runtime.solution.physical_properties.β[i], range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.β[i])))
        end
        for i in 1:runtime.model.body.n_reaction_wheels
            deleteat!(runtime.solution.physical_properties.rw_h[i], range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.rw_h[i])))
            deleteat!(runtime.solution.physical_properties.rw_τ[i], range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.rw_τ[i])))
        end
        
        deleteat!(runtime.solution.physical_properties.τ_rw[1], range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.τ_rw[1])))
        deleteat!(runtime.solution.physical_properties.τ_rw[2], range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.τ_rw[2])))
        deleteat!(runtime.solution.physical_properties.τ_rw[3], range(start=index_b+1,step=1,stop=length(runtime.solution.physical_properties.τ_rw[3])))

        deleteat!(runtime.solution.performance.mass, range(start=index_b+1,step=1,stop=length(runtime.solution.performance.mass)))
        deleteat!(runtime.solution.performance.heat_rate, range(start=index_b+1,step=1,stop=length(runtime.solution.performance.heat_rate)))
        deleteat!(runtime.solution.performance.heat_load, range(start=index_b+1,step=1,stop=length(runtime.solution.performance.heat_load)))
        deleteat!(runtime.solution.performance.T_r, range(start=index_b+1,step=1,stop=length(runtime.solution.performance.T_r)))
        deleteat!(runtime.solution.performance.q, range(start=index_b+1,step=1,stop=length(runtime.solution.performance.q)))

        deleteat!(runtime.solution.forces.gravity_ii[1], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.gravity_ii[1])))
        deleteat!(runtime.solution.forces.gravity_ii[2], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.gravity_ii[2])))
        deleteat!(runtime.solution.forces.gravity_ii[3], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.gravity_ii[3])))
        deleteat!(runtime.solution.forces.drag_pp[1], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.drag_pp[1])))
        deleteat!(runtime.solution.forces.drag_pp[2], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.drag_pp[2])))
        deleteat!(runtime.solution.forces.drag_pp[3], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.drag_pp[3])))
        deleteat!(runtime.solution.forces.drag_ii[1], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.drag_ii[1])))
        deleteat!(runtime.solution.forces.drag_ii[2], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.drag_ii[2])))
        deleteat!(runtime.solution.forces.drag_ii[3], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.drag_ii[3])))
        deleteat!(runtime.solution.forces.lift_pp[1], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.lift_pp[1])))
        deleteat!(runtime.solution.forces.lift_pp[2], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.lift_pp[2])))
        deleteat!(runtime.solution.forces.lift_pp[3], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.lift_pp[3])))
        deleteat!(runtime.solution.forces.lift_ii[1], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.lift_ii[1])))
        deleteat!(runtime.solution.forces.lift_ii[2], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.lift_ii[2]))) 
        deleteat!(runtime.solution.forces.lift_ii[3], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.lift_ii[3]))) 
        deleteat!(runtime.solution.forces.force_ii[1], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.force_ii[1])))
        deleteat!(runtime.solution.forces.force_ii[2], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.force_ii[2])))
        deleteat!(runtime.solution.forces.force_ii[3], range(start=index_b+1,step=1,stop=length(runtime.solution.forces.force_ii[3])))   
        deleteat!(runtime.solution.forces.energy, range(start=index_b+1,step=1,stop=length(runtime.solution.forces.energy)))

        deleteat!(runtime.solution.simulation.MC_seed, range(start=index_b+1,step=1,stop=length(runtime.solution.simulation.MC_seed)))
        deleteat!(runtime.solution.simulation.drag_passage, range(start=index_b+1,step=1,stop=length(runtime.solution.simulation.drag_passage)))

        append!(runtime.solution.orientation.time, t_in_new)
        append!(runtime.solution.orientation.year, results[1])
        append!(runtime.solution.orientation.month, results[2])
        append!(runtime.solution.orientation.day, results[3])
        append!(runtime.solution.orientation.hour, results[4])
        append!(runtime.solution.orientation.minute, results[5])
        append!(runtime.solution.orientation.second, results[6])
        append!(runtime.solution.orientation.number_of_passage, results[7])
        append!(runtime.solution.orientation.pos_ii[1], results[8])
        append!(runtime.solution.orientation.pos_ii[2], results[9])
        append!(runtime.solution.orientation.pos_ii[3], results[10])
        append!(runtime.solution.orientation.vel_ii[1], results[11])
        append!(runtime.solution.orientation.vel_ii[2], results[12])
        append!(runtime.solution.orientation.vel_ii[3], results[13])
        append!(runtime.solution.orientation.pos_ii_mag, results[14])
        append!(runtime.solution.orientation.vel_ii_mag, results[15])
        append!(runtime.solution.orientation.quaternion[1], results[88])
        append!(runtime.solution.orientation.quaternion[2], results[89])
        append!(runtime.solution.orientation.quaternion[3], results[90])
        append!(runtime.solution.orientation.quaternion[4], results[91])
        append!(runtime.solution.orientation.ω[1], results[92])
        append!(runtime.solution.orientation.ω[2], results[93])
        append!(runtime.solution.orientation.ω[3], results[94])

        append!(runtime.solution.orientation.pos_pp[1], results[16])
        append!(runtime.solution.orientation.pos_pp[2], results[17])
        append!(runtime.solution.orientation.pos_pp[3], results[18])
        append!(runtime.solution.orientation.pos_pp_mag, results[19])
        append!(runtime.solution.orientation.vel_pp[1], results[20])
        append!(runtime.solution.orientation.vel_pp[2], results[21])
        append!(runtime.solution.orientation.vel_pp[3], results[22])
        append!(runtime.solution.orientation.vel_pp_mag, results[23])

        append!(runtime.solution.orientation.oe[1], results[24])
        append!(runtime.solution.orientation.oe[2], results[25])
        append!(runtime.solution.orientation.oe[3], results[26])
        append!(runtime.solution.orientation.oe[4], results[27])
        append!(runtime.solution.orientation.oe[5], results[28])
        append!(runtime.solution.orientation.oe[6], results[29])

        append!(runtime.solution.orientation.lat, results[30])
        append!(runtime.solution.orientation.lon, results[31])
        append!(runtime.solution.orientation.alt, results[32])
        append!(runtime.solution.orientation.γ_ii, results[33])
        append!(runtime.solution.orientation.γ_pp, results[34])

        append!(runtime.solution.orientation.h_ii[1], results[35])
        append!(runtime.solution.orientation.h_ii[2], results[36])
        append!(runtime.solution.orientation.h_ii[3], results[37])
        append!(runtime.solution.orientation.h_pp[1], results[38])
        append!(runtime.solution.orientation.h_pp[2], results[39])
        append!(runtime.solution.orientation.h_pp[3], results[40])
        append!(runtime.solution.orientation.h_ii_mag, results[41])
        append!(runtime.solution.orientation.h_pp_mag, results[42])

        append!(runtime.solution.orientation.uD[1], results[43])
        append!(runtime.solution.orientation.uD[2], results[44])
        append!(runtime.solution.orientation.uD[3], results[45])
        append!(runtime.solution.orientation.uE[1], results[46])
        append!(runtime.solution.orientation.uE[2], results[47])
        append!(runtime.solution.orientation.uE[3], results[48])
        append!(runtime.solution.orientation.uN[1], results[49])
        append!(runtime.solution.orientation.uN[2], results[50])
        append!(runtime.solution.orientation.uN[3], results[51])
        append!(runtime.solution.orientation.vN, results[52])
        append!(runtime.solution.orientation.vE, results[53])
        append!(runtime.solution.orientation.azi_pp, results[54])

        append!(runtime.solution.physical_properties.ρ, results[55])
        append!(runtime.solution.physical_properties.T, results[56])
        append!(runtime.solution.physical_properties.p, results[57])
        append!(runtime.solution.physical_properties.wind[1], results[58])
        append!(runtime.solution.physical_properties.wind[2], results[59])
        append!(runtime.solution.physical_properties.wind[3], results[60])
        append!(runtime.solution.physical_properties.cL, results[61])
        append!(runtime.solution.physical_properties.cD, results[62])
        append!(runtime.solution.physical_properties.S, results[63])
        append!(runtime.solution.physical_properties.α_control, results[95])
        append!(runtime.solution.physical_properties.τ_rw[1], results[96,:]) # h_rw_mag
        append!(runtime.solution.physical_properties.τ_rw[2], results[97,:]) # τ_rw_x
        append!(runtime.solution.physical_properties.τ_rw[3], results[98,:]) # τ_rw_y
        
        # Initialize α and β if they are not already initialized
        n_bodies = length(runtime.model.body.links)
        if isempty(runtime.solution.physical_properties.α)
            for i in 1:n_bodies
                append!(runtime.solution.physical_properties.α, [[]])
                append!(runtime.solution.physical_properties.β, [[]])
                append!(runtime.solution.performance.heat_rate, [[]])
                append!(runtime.solution.performance.heat_load, [[]])
            end
        end

        # Append α and β for each link
        for i in 1:n_bodies
            append!(runtime.solution.physical_properties.α[i], results[98 + i,:]) # α
            append!(runtime.solution.physical_properties.β[i], results[98 + n_bodies + i,:]) # β
            append!(runtime.solution.performance.heat_rate[i], results[98 + 2*n_bodies + i,:]) # heat_rate
            append!(runtime.solution.performance.heat_load[i], results[98 + 3*n_bodies + i,:]) # heat_load
        end
        
        n_reaction_wheels = runtime.model.body.n_reaction_wheels
        # Initialize the reaction wheel properties if they are not already initialized
        if isempty(runtime.solution.physical_properties.rw_h)
            for i in 1:n_reaction_wheels
                append!(runtime.solution.physical_properties.rw_h, [[]])
                append!(runtime.solution.physical_properties.rw_τ, [[]])
            end
        end

        for i in 1:n_reaction_wheels
            append!(runtime.solution.physical_properties.rw_h[i], results[98 + 4*n_bodies + i,:])
            append!(runtime.solution.physical_properties.rw_τ[i], results[98 + 4*n_bodies + n_reaction_wheels + i,:]) # rw_τ
        end
 
    
        append!(runtime.solution.performance.mass, results[64])
        # append!(runtime.solution.performance.heat_rate, results[65])
        # append!(runtime.solution.performance.heat_load, results[66])
        append!(runtime.solution.performance.T_r, results[65])
        append!(runtime.solution.performance.q, results[66])

        append!(runtime.solution.forces.gravity_ii[1], results[67])
        append!(runtime.solution.forces.gravity_ii[2], results[68])
        append!(runtime.solution.forces.gravity_ii[3], results[69])
        append!(runtime.solution.forces.drag_pp[1], results[70])
        append!(runtime.solution.forces.drag_pp[2], results[71])
        append!(runtime.solution.forces.drag_pp[3], results[72])
        append!(runtime.solution.forces.drag_ii[1], results[73])
        append!(runtime.solution.forces.drag_ii[2], results[74])
        append!(runtime.solution.forces.drag_ii[3], results[75])
        append!(runtime.solution.forces.lift_pp[1], results[76])
        append!(runtime.solution.forces.lift_pp[2], results[77])
        append!(runtime.solution.forces.lift_pp[3], results[78])
        append!(runtime.solution.forces.lift_ii[1], results[79])
        append!(runtime.solution.forces.lift_ii[2], results[80])
        append!(runtime.solution.forces.lift_ii[3], results[81])
        append!(runtime.solution.forces.force_ii[1], results[82])
        append!(runtime.solution.forces.force_ii[2], results[83])
        append!(runtime.solution.forces.force_ii[3], results[84])
        append!(runtime.solution.forces.energy, results[85])

        append!(runtime.solution.simulation.MC_seed, results[86])
        append!(runtime.solution.simulation.drag_passage, results[87])

        return initial_state
    end
end
