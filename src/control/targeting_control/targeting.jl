include("sim_targeting.jl")
include("../utils/Eom_ctrl.jl")
# include("../heatload_control/Utils_timeswitch.jl")
# include("../../utils/Closed_form_solution.jl")

using Roots

function target_planning(f!, ip, m, args, param, OE, initial_time, final_time, a_tol, r_tol, method, events, in_cond)

    # OE = SVector{7, Float64}([initial_state.a, initial_state.e, initial_state.i, initial_state.Ω, initial_state.ω, initial_state.vi, initial_state.m])
    r0, v0 = orbitalelemtorv(OE, m.planet)

    # Run simulation
    prob = ODEProblem(f!, in_cond, (initial_time, final_time), param)
    sol_max_dratio = solve(prob, method, abstol=a_tol, reltol=r_tol, callback=events)

    println("m.aero.alpha before deepcopying: ", m.aerodynamics.α)

    # Minimum drag ratio Run
    ip_temp = deepcopy(ip)
    ip_temp.cm = 0

    m_temp = deepcopy(m)
    m_temp.aerodynamics.α = 0.0

    println("m.aero.alpha after deepcopying: ", m.aerodynamics.α)

    config.cnf.ascending_phase = false
    config.cnf.drag_state = true

    # param_temp = param

    # println(param_temp[1])
    # println(" ")
    # println(m_temp)

    param_temp = (m_temp, param[2], ip_temp, param[4:end]...)

    # param_temp[1] = m_temp
    # param_temp[3] = ip_temp

    # Run simulation
    prob = ODEProblem(f!, in_cond, (initial_time, final_time), param_temp)
    sol_min_dratio = solve(prob, method, abstol=a_tol, reltol=r_tol, callback=events)

    println("m.aero.alpha after deepcopying and running zero aoa case: ", m.aerodynamics.α)

    energy_target_min = norm(sol_max_dratio[4:6, end])^2/2 - m.planet.μ / norm(sol_max_dratio[1:3, end]) # Lowest possible energy with maximum drag ratio, most negative
    energy_target_max = norm(sol_min_dratio[4:6, end])^2/2 - m.planet.μ / norm(sol_min_dratio[1:3, end]) # Highest possible energy with minimum drag ratio, least negative

    println("Energy target min: ", energy_target_min)
    println("Energy target max: ", energy_target_max)

    h0 = norm(cross(r0, v0))
    r_p = h0^2 / (m.planet.μ * (1 + OE[2]))

    # println("Current periapsis: ", r_p - m.planet.Rp_e)

    target_energy = -m.planet.μ / (args[:ra_fin_orbit] + r_p) # change to current periapsis

    println("Target energy: ", target_energy)

    if target_energy < energy_target_max && target_energy > energy_target_min
        config.cnf.targeting = 1
    elseif target_energy < energy_target_min
        config.cnf.targeting = 0
    else
        println("Cannot target energy level that is larger than possible with minimum drag ratio")
    end

    # config.cnf.hf = args[:AE]*1e3

    # config.cnf.Vf = sqrt(2*(target_energy + m.planet.μ / (m.planet.Rp_e + config.cnf.hf)))

    # a_max = -m.planet.μ / (2 * target_energy) # maximum semi-major axis for the target energy level

    # r_af = 2*a_max - r_p

    # e_max = (r_af - r_p) / (r_af + r_p)

    # h_max = sqrt(m.planet.μ * a_max * (1 - e_max^2))

    # vi_EI = acos((h_max^2/(m.planet.μ * (m.planet.Rp_e + config.cnf.hf)) - 1)/e_max)     # true anomaly at entry interface

    # V_perp_f = m.planet.μ / h_max * (1 + e_max*cos(vi_EI))          # perpendicular velocity at atmospheric interface

    # V_rad_f = m.planet.μ / h_max * e_max * sin(vi_EI)               # radial velocity at atmospheric interface

    # config.cnf.γf = atan(V_rad_f / V_perp_f)                        # flight path angle at atmospheric interface

    # ip.cm = ip_cm_copy
    # m.aerodynamics.α = deg2rad(args[:α])

    return target_energy
end

# function control_solarpanels_targeting(f!, energy_f, ip, m, time_0, OE, args, gram_atmosphere)
function control_solarpanels_targeting_num_int(energy_f, param, time_0, in_cond)

    function func_targeting_num_int(t_switch)

        sol = asim_ctrl_targeting(t_switch, param, time_0, in_cond)

        m = param[1]

        energy_fin = norm(sol[4:6,end])^2/2 - m.planet.μ / norm(sol[1:3,end])

        println("t_switch: ", t_switch, " energy_fin: ", energy_fin)

        return (energy_fin - energy_f) / 1e6
    end

    t_switch = find_zero(ts -> func_targeting_num_int(ts), [0, 600], Roots.Brent(), verbose=true, rtol=1e-5)

    return t_switch 
end

function control_solarpanels_targeting_heatload(energy_f, param, OE)

    function func_targeting_heatload(v_E)
        m = param[1]
        ip = param[3]
        time_0 = param[7]
        args = param[8]
        gram_atmosphere = param[10]

        sol, _ = asim_ctrl_rf(ip, m, time_0, OE, args, v_E, 1.0, false, gram_atmosphere)

        energy_fin = norm(sol[4:6,end])^2/2 - m.planet.μ / norm(sol[1:3,end])

        println("v_E: ", v_E, " energy_fin: ", energy_fin)

        return (energy_fin - energy_f) / 1e6
    end

    v_E_fin = find_zero(v_E -> func_targeting_heatload(v_E), [1, 1000], Roots.Brent(), verbose=true, rtol=1e-8)

    return v_E_fin
end

# function control_solarpanels_targeting_closed_form(energy_f, param, initialcondition)

#     m = param[1]
#     args = param[8]

#     T = m.planet.T

#     # t_cf, _, _, _ = closed_form(args, m, OE, T, true, m.aerodynamics.α)

#     if config.cnf.count_numberofpassage != 1
#         t_prev = config.solution.orientation.time[end]
#     else
#         t_prev = m.initial_condition.time_rot # value(seconds(date_initial - from_utc(DateTime(2000, 1, 1, 12, 0, 0)))) # mission.initial_condition.time_rot
#     end

#     pos_ii_org, vel_ii_org = orbitalelemtorv(initialcondition, m.planet)
#     pos_ii = SVector{3, Float64}(pos_ii_org)
#     vel_ii = SVector{3, Float64}(vel_ii_org)

#     r0 = norm(pos_ii) # Inertial position magnitude
#     v0 = norm(vel_ii) # Inertial velocity magnitude
#     h0 = r0 - m.planet.Rp_e

#     pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, config.cnf.et)

#     LatLong = rtolatlong(pos_pp, m.planet)
#     lat = LatLong[2]
#     lon = LatLong[3]
#     h0 = LatLong[1]

#     h_ii = cross(pos_ii, vel_ii)
#     arg = median([-1, 1, norm(h_ii)/(r0*v0)])   # limit to[-1, 1]
#     γ0 = acos(arg)

#     if dot(pos_ii, vel_ii) < 0
#         γ0 = -γ0
#     end

#     initial_state_angle = initialcondition[6]
#     e = initialcondition[2]
#     a = initialcondition[1]
#     final_state_angle = -initial_state_angle
#     E_initialstate = 2 * atan(sqrt((1-e)/(1+e)) * tan(initial_state_angle/2))
#     E_finalstate = 2 * atan(sqrt((1-e)/(1+e)) * tan(final_state_angle/2))

#     # Evaluate time to reach next state
#     Δt = sqrt(a^3 / m.planet.μ) * ((E_finalstate - e*sin(E_finalstate)) - (E_initialstate - e*sin(E_initialstate)))
#     t_p = Δt/2

#     mass = initialcondition[end]

#     if h0 < args[:EI]*1e3 #if initial condition are lower than drag passage initial condition # this happens only running MC cases
#         # let's calculate pos_ii,v_ii for the point of trajectory corresponding to h = 160 km
#         h0 = args[:EI]*1e3
#         r = m.planet.Rp_e + h0
#         OE = rvtoorbitalelement(pos_ii, vel_ii, mass, m.planet)
#         a, e, i, Ω, ω, vi = OE[1], OE[2], OE[3], OE[4], OE[5], OE[6]
#         vi = 2*pi - acos(((a*(1 - e^2)/r)-1)/e)
#         E_real_finalstate = 2 * atan(sqrt((1-e)/(1+e)) * tan(-vi/2)) # eccentric anomaly
#         Δt = sqrt(a^3 / m.planet.μ) * ((E_real_finalstate - e*sin(E_real_finalstate)) - (E_initialstate - e*sin(E_initialstate)))
#         t_p = Δt/2
#     end

#     println("h0 = ", h0)

#     function func_targeting_cf(t_switch)
#         # t_switch =350.0

#         t_cf = collect(range(start=0, stop=t_switch, step=0.1))

#         aoa_cf1 = ones(length(t_cf)) * pi/2

#         t_cf1, h_cf1, γ_cf1, v_cf1 = closed_form_targeting(0, m, (v0, γ0, h0), T, t_cf, t_p, mass, aoa_cf1)

#         function cf_switching!(resid, z)
#             v0_n = z[1]
#             γ0_n = z[2]

#             aoa_cf = zeros(length(t_cf))

#             t_cf_n, h_cf_n, γ_cf_n, v_cf_n = closed_form_targeting(0, m, (v0_n, γ0_n, h0), T, t_cf, t_p, mass, aoa_cf)

#             resid[1] = v_cf_n[end] - v_cf1[end]
#             resid[2] = γ_cf_n[end] - γ_cf1[end]
#         end

#         # println(v_cf1[end], " ", γ_cf1[end], " ", t_cf1[end])

#         # v0_2, γ0_2, h0_2 = v_cf1[end], γ_cf1[end], h0 # h_cf1[end]

#         # t_cf = collect(range(start=t_cf1[end], stop=2*t_p, step=0.1))
#         # aoa_cf2 = zeros(length(t_cf))

#         # t_cf2, h_cf2, γ_cf2, v_cf2 = closed_form_targeting(t_cf1[end], m, (v0_2, γ0_2, h0_2), T, t_cf, t_p, mass, aoa_cf2)

#         # println(t_cf[end])
#         # println(h_cf2[end])
#         # println(v_cf2[end])

#         z0 = [v0, γ0]

#         sol_NL = nlsolve((resid, z) -> cf_switching!(resid, z), z0, show_trace=false)

#         # println(sol_NL)

#         v0_2, γ0_2, h0_2 = sol_NL.zero[1], sol_NL.zero[2], h0 # h_cf1[end]

#         t_cf2 = collect(range(start=0, stop=2*t_p, step=0.1))
#         aoa_cf2 = zeros(length(t_cf2))

#         t_cf2, h_cf2, γ_cf2, v_cf2 = closed_form_targeting(t_cf2[1], m, (v0_2, γ0_2, h0_2), T, t_cf2, t_p, mass, aoa_cf2)

#         energy_fin = v_cf2[end]^2/2 - m.planet.μ / (m.planet.Rp_e + h_cf2[end])

#         # println("Final energy (closed form): ", energy_fin)
#         # println(" ")

#         return energy_fin - energy_f
#     end

#     t_switch = find_zero(ts -> func_targeting_cf(ts), [1, 2*t_p - 1], Roots.Brent(), verbose=true)

#     return [t_switch, 2*t_p]
# end

function control_solarpanels_targeting_closed_form(energy_target, ip, m, position, args, t, heat_rate_control, reevaluation_mode, current_position=0)
    k_cf = 1.0

    T = m.planet.T

    t_cf, h_cf, γ_cf, v_cf = closed_form(args, m, position, T, true, m.aerodynamics.α)

    RT = T * m.planet.R
    S = (v_cf / sqrt(2 * RT))

    CL_90, CD_90 = aerodynamic_coefficient_fM(pi/2, m.body, T, S[1], m.aerodynamics, 0)
    CL_0, CD_0 = aerodynamic_coefficient_fM(0, m.body, T, S[1], m.aerodynamics, 0)

    CD_slope = (CD_90 - CD_0) / (pi/2)

    coeff = (CD_slope, CL_0, CD_0)

    approx_sol = (t_cf, h_cf, γ_cf, v_cf)

    # aoa_cf = aoa(m, k_Cf, t_cf, h_cf, γ_cf, v_cf, coeff, 1)[1]

    # delta_E_max = func_e(0, m, args, coeff, position, heat_rate_control, approx_sol, aoa_cf, energy_target)
    # delta_E_min = func_e(100, m, args, coeff, position, heat_rate_control, approx_sol, zeros(length(aoa_cf)), energy_target)

    # println("delta_E_max: ", delta_E_max)
    # println("delta_E_min: ", delta_E_min)

    # if delta_E_max * delta_E_min < 0
    nu_E_root = fzero(nu_E -> func_e(nu_E, m, args, coeff, position, heat_rate_control, approx_sol, energy_target), [1, 100], Roots.Brent())
    # elseif delta_E_max < 0.0
    #     return [0.0, 0.0]
    # elseif delta_E_min > 0.0
    #     return [0.0, t_cf[end]/2]
    # end

    t_cf, v_cf, γ_cf, h_cf = func_e(k_cf, m, args, coeff, position, heat_rate_control, approx_sol, energy_target, false, true)

    lambda_switch, lambdav = lambdas(m, aoa_cf, k_cf, t_cf, h_cf, γ_cf, v_cf, coeff, nu_E_root)[1:2]

    index_array = lambdav .< lambda_switch

    temp = t_cf .* index_array

    temp = filter(!iszero, temp)

    t_switch = [temp[1], temp[end]]

    if abs(t_switch[end] - t_cf[end]) < 5
        t_switch[1] -= t_cf[end]*0.04
    elseif abs(t_switch[end] - t_cf[end]) < 60
        t_switch[1] -= 5
    end

    t_switch[2] -= t_switch[2]*0.1

    return t_switch
end

function func_e(nu_E_var, m, args, coeff, position, heat_rate_control, approx_sol, energy_target, initial_guess=false, approx_calc=false)
        
        k_cf = 1.0
        
        t_cf, h_cf, γ_cf, v_cf = approx_sol

        E_fin = 0
        temp_E_diff = 1e10
        E_prev = 0

        count = 0

        aoa_cf = []
        
        while temp_E_diff > 1e-3 # stop if two following trajectory evaluation provides the same E
            # define angle of attack lagrangian multipliers
            T = m.planet.T  # fixed temperature
            aoa_cf, in_cond_lambda = aoa(m, k_cf, t_cf, h_cf, γ_cf, v_cf, coeff, nu_E_var, aoa_cf)  # update angle of attack profile with new k
            t_cf, h_cf, γ_cf, v_cf = closed_form(args, m, position, T, true, m.aerodynamics.α, aoa_cf)  # re-evaluate the closed form solution using previous angle of attack profile

            a = sqrt(m.planet.γ * m.planet.R * T)
            M = v_cf / a
            S = sqrt(m.planet.γ/2) * M
            ρ = density_exp(h_cf, m.planet)[1]  # density calculated through exponential density

            heat_rate = heat_rate_calc(args[:multiplicative_factor_heatload] * m.aerodynamics.thermal_accomodation_factor, ρ, T, T, m.planet.R, m.planet.γ, S, aoa_cf)

            if heat_rate_control == true # Account for max heat rate possible
                index_hr = heat_rate .> args[:max_heat_rate]
                heat_rate[index_hr] .= args[:max_heat_rate]
            end

            E_fin = norm(v_cf[end])^2/2 - m.planet.μ / (m.planet.Rp_e + h_cf[end]) 

            # update error
            temp_E_diff = E_fin - E_prev

            E_prev = E_fin
            count += 1
        end

        if initial_guess == true
            return in_cond_lambda
        elseif approx_calc == true
            return t_cf, v_cf, γ_cf, h_cf
        end

        delta_E = E_fin - energy_target

        return delta_E
    end