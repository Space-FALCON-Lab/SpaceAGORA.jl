include("../utils/Reference_system.jl")
include("../utils/Ref_system_conf.jl")


function impact(t, y, m, solution, args)
    if sqrt(y[1]^2 + y[2]^2 + y[3]^2) <= m.planet.Rp_e + 35
        breaker = true
        print("ATTENTION: AMPACT!")
        solution.t_events[end] = ["true"]
    else
        breaker = false
    end

    return breaker, solution
    
end

function apoapsispoint(t, y_0, y_h, m, args)
    pos_ii = [y_n[1], y_n[2], y_n[3]]   # Inertial position
    vel_ii = [y_0[4], y_0[5], y_0[6]]   # Inertial velocity
    r_i = cartesian(pos_ii[1], pos_ii[2], pos_ii[3])
    v_i = cartesian(vel_ii[1], vel_ii[2], vel_ii[3])
    OE_n = rvtoorbitalelement(r_i, v_i, y[7], m.planet)

    pos_ii = [y_0[1], y_0[2], y_0[3]]   # Inertial position
    vel_ii = [y_0[4], y_0[5], y_0[6]]   # Inertial velocity
    r_i = cartesian(pos_ii[1], pos_ii[2], pos_ii[3])
    v_i = cartesian(vel_ii[1], vel_ii[2], vel_ii[3])
    OE_0 = rvtoorbitalelement(r_i, v_i, y[7], m.planet)

    if OE_n[6] - pi > 0 && OE_0[6] - pi < 0
        breaker = true
    else
        breaker = false
    end

    return breaker
end 

function apoapsisgreaterperiapsis(t, y, m, solution, args)
    r = y[1:3]
    v = y[4:6]

    Energy = dot(v,v) / 2 - m.planet.μ / sqrt(dot(r,r))
    a = m.planet.μ / (2 * Energy)
    h = cross(r, v)
    e = sqrt(1 + (2 * Energy * dot(h,h)/m.planet.μ^2))

    r_a = a * (1 + e)
    r_p = a * (1 - e)

    if r_a < r_p
        breaker = true
        solution.t_events[end] = ["true"]
    else
        breaker = false
    end

    return breaker, solution
end

function eventsecondstep(t_n, y_n, t_0, y_0, m, args)
    alt_y_n = norm(y_n[1:3]) - m.planet.Rp_e - 250*1e3
    alt_y_0 = norm(y_0[1:3]) - m.planet.Rp_e - 250*1e3

    if alt_y_0 < 0 && alt_y_n >= 0
        breaker = true
    else
        breaker = false
    end

    return breaker
end

function heat_check(t_n, y_n, t_0, y_0, m, args)
    alt_y_n = norm(y_n[1:3]) - m.planet.Rp_e - 120*1e3
    alt_y_0 = norm(y_0[1:3]) - m.planet.Rp_e - 120*1e3

    if alt_y_0 < 0 && alt_y_n > 0
        breaker = true
    else
        breaker = false
    end

    return breaker
end

function out_drag_passage(t_n, y_n, t_0, y_0, m, args)
    alt_y_n = norm(y_n[1:3]) - m.planet.Rp_e - args[:AE]
    alt_y_0 = norm(y_0[1:3]) - m.planet.Rp_e - args[:AE]

    if alt_y_0 < 0 && alt_y_n > 0
        breaker = true
    else
        breaker = false
    end

    return breaker
end

function stop_firing(t_n, y_n, m, args)
    mass = y_n[7]
    Δv = (m.engine.g_e * m.engine.Isp) * log(m.body.Mass/mass)
    Δv_max = Odyssey_firing_plan(Δv, args)

    if Δv >= (Δv_max + config.cnf.delta_v_man)
        config.cnf.firing_on = 0
    end

    return Δv - (Δv_max + config.cnf.delta_v_man)
end

function events(t_n, y_n, t_0, y_0, m, T_ijk, index_phase_aerobraking, args, solution)
    
    while True
        impact_breaker = false

        # First event
        breaker, solution = impact(T_n, y_n, m, solution, args) # Stop simulation if impact occurs

        if breaker == True
            impact_breaker = true
            break
        end

        # Second event
        if index_phase_aerobraking == 1
            breaker = apoapsispoint(t_n, y_0, y_n, m, args)

            stop_firing(t_n, y_n, m, args)
        elseif index_phase_aerobraking == 3
            if args[:control_mode] != 0
                breaker = heat_check(t_n, y_n, t_0, y_0, m, args)
            elseif args[:drag_passage] == true || args[:body_shape] == "Blunted Cone"
                breaker = out_drag_passage(t_n, y_n, t_0, y_0, m, args)
            else
                breaker = eventsecondstep(t_n, y_n, t_0, y_0, m, args)
            end
        elseif index_phase_aerobraking == 2
            stop_firing(t_n, y_n, m, args)
        end

        if breaker == true
            break
        end

        # Third event
        breaker, solution = apoapsisgreaterperiapsis(t_n, y_n, m, solution, args)

        break
    end

    return breaker, impact_breaker, solution
end


function RK4(f, h, t, y, m, T_ijk, index_phase_aerobraking, args, solution)
    k_1 = f(t, y, m, index_phase_aerobraking)
    k_2 = f(t + h/2, y + k_1 * h/2, m, index_phase_aerobraking)
    k_3 = f(t + h/2, y + k_2 * h/2, m, index_phase_aerobraking)
    k_4 = f(t + h, y + h * k_3, m, index_phase_aerobraking)

    summation = []
    for i in range(length(k_1))
        append!(summation, k_1[i] + 2*k_2[i] + k_3[i] + k_4[i])
    end

    y_n = [0, 0, 0, 0, 0, 0, 0, 0]
    for i in range(length(k_1))
        y_n[i] = y[i] + h*summation[i]/6
    end

    t_n = t + h

    # Check terminal events
    breaker, imapct_breaker, solution = events(t_n, y_n, t, y, m, T_ijk, index_phase_aerobraking, args, solution)

    return y_n, t_n, breaker, solution
end