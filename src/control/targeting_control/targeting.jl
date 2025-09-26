include("sim_targeting.jl")
include("../utils/Eom_ctrl.jl")

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

    v_E_fin = find_zero(v_E -> func_targeting_heatload(v_E), [1, 100], Roots.Brent(), verbose=true, rtol=1e-8)

    return v_E_fin
end

function control_solarpanels_targeting_closed_form(energy_f, param, initialcondition)

    m = param[1]
    args = param[8]

    T = m.planet.T

    # t_cf, _, _, _ = closed_form(args, m, OE, T, true, m.aerodynamics.α)

    if config.cnf.count_numberofpassage != 1
        t_prev = config.solution.orientation.time[end]
    else
        t_prev = m.initial_condition.time_rot # value(seconds(date_initial - from_utc(DateTime(2000, 1, 1, 12, 0, 0)))) # mission.initial_condition.time_rot
    end

    pos_ii_org, vel_ii_org = orbitalelemtorv(initialcondition, m.planet)
    pos_ii = SVector{3, Float64}(pos_ii_org)
    vel_ii = SVector{3, Float64}(vel_ii_org)

    r0 = norm(pos_ii) # Inertial position magnitude
    v0 = norm(vel_ii) # Inertial velocity magnitude
    h0 = r0 - m.planet.Rp_e

    pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, config.cnf.et)

    LatLong = rtolatlong(pos_pp, m.planet)
    lat = LatLong[2]
    lon = LatLong[3]
    h0 = LatLong[1]

    h_ii = cross(pos_ii, vel_ii)
    arg = median([-1, 1, norm(h_ii)/(r0*v0)])   # limit to[-1, 1]
    γ0 = acos(arg)

    if dot(pos_ii, vel_ii) < 0
        γ0 = -γ0
    end

    initial_state_angle = initialcondition[6]
    e = initialcondition[2]
    a = initialcondition[1]
    final_state_angle = -initial_state_angle
    E_initialstate = 2 * atan(sqrt((1-e)/(1+e)) * tan(initial_state_angle/2))
    E_finalstate = 2 * atan(sqrt((1-e)/(1+e)) * tan(final_state_angle/2))

    # Evaluate time to reach next state
    Δt = sqrt(a^3 / m.planet.μ) * ((E_finalstate - e*sin(E_finalstate)) - (E_initialstate - e*sin(E_initialstate)))
    t_p = Δt/2

    mass = initialcondition[end]

    if h0 < args[:EI]*1e3 #if initial condition are lower than drag passage initial condition # this happens only running MC cases
        # let's calculate pos_ii,v_ii for the point of trajectory corresponding to h = 160 km
        h0 = args[:EI]*1e3
        r = m.planet.Rp_e + h0
        OE = rvtoorbitalelement(pos_ii, vel_ii, mass, m.planet)
        a, e, i, Ω, ω, vi = OE[1], OE[2], OE[3], OE[4], OE[5], OE[6]
        vi = 2*pi - acos(((a*(1 - e^2)/r)-1)/e)
        E_real_finalstate = 2 * atan(sqrt((1-e)/(1+e)) * tan(-vi/2)) # eccentric anomaly
        Δt = sqrt(a^3 / m.planet.μ) * ((E_real_finalstate - e*sin(E_real_finalstate)) - (E_initialstate - e*sin(E_initialstate)))
        t_p = Δt/2
    end

    function func_targeting_cf(t_switch)
        # t_switch =350.0

        t_cf = collect(range(start=0, stop=t_switch, step=0.1))

        aoa_cf1 = ones(length(t_cf)) * m.aerodynamics.α

        # println(aoa_cf)
        println("t_switch", t_switch)
        println(" ")

        t_cf1, h_cf1, γ_cf1, v_cf1 = closed_form_targeting(0, m, (v0, γ0, h0), T, t_cf, t_p, mass, aoa_cf1)

        println(v_cf1[end], " ", γ_cf1[end], " ", t_cf1[end])

        v0_2, γ0_2, h0_2 = v_cf1[end], γ_cf1[end], h0 # h_cf1[end]

        t_cf = collect(range(start=t_cf1[end], stop=2*t_p, step=0.1))
        aoa_cf2 = zeros(length(t_cf))

        t_cf2, h_cf2, γ_cf2, v_cf2 = closed_form_targeting(t_cf1[end], m, (v0_2, γ0_2, h0_2), T, t_cf, t_p, mass, aoa_cf2)

        # println(t_cf[end])
        println(h_cf2[end])
        println(v_cf2[end])

        energy_fin = v_cf2[end]^2/2 - m.planet.μ / (m.planet.Rp_e + h_cf2[end])

        println("Final energy (closed form): ", energy_fin)
        println(" ")

        return energy_fin - energy_f
    end

    t_switch = find_zero(ts -> func_targeting_cf(ts), [1, 2*t_p - 1], Roots.Brent(), verbose=true)

    return t_switch
end