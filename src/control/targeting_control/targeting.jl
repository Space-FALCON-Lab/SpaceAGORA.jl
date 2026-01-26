include("../utils/Eoms.jl")

function target_planning(ip, m, args, param, initial_state, initial_time, final_time, a_tol, r_tol, method, events, )

    # Run simulation
    prob = ODEProblem(f!, in_cond, (initial_time, final_time), param)
    sol_max_dratio = solve(prob, method, abstol=a_tol, reltol=r_tol, callback=events)

    # Minimum drag ratio Run
    ip_temp = ip
    ip_temp.cm = 0

    m_temp = m
    m_temp.aerodynamics.α = 0.0

    param_temp = param
    param_temp[1] = m_temp
    param_temp[3] = ip_temp

    # Run simulation
    prob = ODEProblem(f!, in_cond, (initial_time, final_time), param_temp)
    sol_min_dratio = solve(prob, method, abstol=a_tol, reltol=r_tol, callback=events)

    energy_target_max = norm(sol_max_dratio[4:6, end])^2 * 0.5 - m.planet.μ / norm(sol_max_dratio[1:3, end])
    energy_target_min = norm(sol_min_dratio[4:6, end])^2 * 0.5 - m.planet.μ / norm(sol_min_dratio[1:3, end])
    
    h0 = norm(cross(r0, v0))
    r_p = h0^2 / (m.planet.μ * (1 + initial_state.e))               # periapsis radius

    target_energy = -m.planet.μ / (args[:ra_fin_orbit] + (r_p + m.planet.Rp_e)) # change to current periapsis

    if target_energy < energy_target_max
        target_energy = energy_target_max
        config.cnf.targeting = 0
    elseif target_energy < energy_target_min && target_energy > energy_target_max
        config.cnf.targeting = 1
    else
        println("Cannot target energy level that is larger than possible with minimum drag ratio")
    end

    config.cnf.hf = args[:AE]*1e3

    config.cnf.Vf = sqrt(2*(target_energy + m.planet.μ / (m.planet.Rp_e + config.cnf.hf)))

    a_max = -m.planet.μ / (2 * target_energy) # maximum semi-major axis for the target energy level

    r_af = 2*a_max - r_p

    e_max = (r_af - r_p) / (r_af + r_p) 

    h_max = sqrt(m.planet.μ * a_max * (1 - e_max^2))

    vi_EI = acos((h_max^2/(m.planet.μ * (m.planet.Rp_e + config.cnf.hf)) - 1)/e_max)     # true anomaly at entry interface

    V_perp_f = m.planet.μ / h_max * (1 + e_max*cos(vi_EI))          # perpendicular velocity at atmospheric interface

    V_rad_f = m.planet.μ / h_max * e_max * sin(vi_EI)               # radial velocity at atmospheric interface

    config.cnf.γf = atan(V_rad_f / V_perp_f)                        # flight path angle at atmospheric interface
end

function control_solarpanels_targeting()
    
end