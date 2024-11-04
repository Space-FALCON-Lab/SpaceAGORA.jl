include("../utils/Reference_system.jl")
include("../utils/Ref_system_conf.jl")
include("../integrator/Integrators.jl")
include("../integrator/Events.jl")
include("../utils/Save_results.jl")
include("../utils/Odyssey_maneuver_plan.jl")
include("../config.jl")

include("../physical_models/Gravity_models.jl")
include("../physical_models/Density_models.jl")
include("../physical_models/Aerodynamic_models.jl")
include("../physical_models/Thermal_models.jl")
# Include Control models
# Inlcude Thrust Maneuvers for control models

using LinearAlgebra
using DifferentialEquations

function asim(ip, m, initial_state, numberofpassage, args)

    if ip.wm == 1
        wind_m = true
    end

    if ip.mc == 1
        MonteCarlo = true
    end

    OE = [initial_State.a, initial_state.e, initial_state.i, intital_state.Ω, initial_state.ω, initial_state.vi, intital_state.m]

    if (OE[1] > (3400 + 50 + 200)*1e3) && (args[:drag_passage] == false) && (args[:body_shape] == "Spacecraft")
        index_steps_EOM = 3
    else
        index_steps_EOM = 1
    end

    i = OE[3]
    Ω = OE[4]
    ω = OE[5]

    T_ijk = [[cos(Ω)*cos(ω) - sin(Ω)*sin(ω)*cos(i), sin(Ω)*cos(ω) + cos(Ω)*sin(ω)*cos(i), sin(ω) * sin(i)],
             [-cos(Ω)*sin(ω) - sin(Ω)*cos(ω)*cos(i), -sin(Ω)*sin(ω) + cos(Ω)*cos(ω)*cos(i), cos(ω)*sin(i)], 
             [sin(Ω)*sin(i), -cos(Ω)*sin(i), cos(i)]]

    T_ijk = [[x == -0 ? 0.0 : x for x in row] for row in T_ijk]

    [r0, v0] = orbitalelemtorv(OE, m.planet)
    mass = OE[end]

    # Clock
    date_initial = DateTime(m.initialcondition.year, m.initialcondition.month, m.initialcondition.day, m.initialcondition.hour, m.initialcondition.minute, m.initialcondition.second)

    config.count_number_of_passage = config.count_aerobraking + 1

    if config.count_number_of_passage ! = 1
        t_prev = config.solution.orientation.time[end]
    else
        t_prev = m.initialcondition.time_rot
    end

    function f(t0, in_cond, m, index_phase_aerobraking, ip)
        ## Counters
        # Counter for all along the simulation of all passages
        config.count_aerobraking = config.count_aerobraking + 1
        passage_number = config.count_aerobraking
        # Counter for one entire passage
        config.count_dori = config.count_dori + 1
        # Counter for one phase
        config.count_phase = congif.count_phase + 1

        # Clock
        time_real = date_initial + Second(t0)
        timereal = clock(Dates.year(time_real), Dates.month(time_real), Dates.day(time_real), Dates.hour(time_real), Dates.minute(time_real), Dates.second(time_real))

        # Assign state
        pos_ii = in_cond[1:3]       # Inertial position 
        pos_ii += 0
        vel_ii = in_cond[4:7]       # Inertial velocity
        vel_ii += 0
        mass = in_cond[8]           # Mass kg
        pos_ii_mag = norm(pos_ii)   # Magnitude of the inertial position
        vel_ii_mag = norm(vel_ii)   # Magnitude of the inertial velocity

        # Assign parameters
        ω_planet = m.planet.ω
        γ = m.planet.γ
        μ_fluid = m.planet.μ_fluid
        area_tot = m.body.area_tot

        # TRANSFORM THE STATE
        # Inertial to planet relative transformation
        [pos_pp, vel_pp] = r_intor_p(pos_ii, vel_ii, m.planet, t0, t_prev) # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
        pos_pp_mag = norm(pos_pp) # Magnitude of the planet relative position
        pos_pp_hat = pos_pp / pos_pp_mag # Unit vector of the planet relative position
        pos_ii_hat = pos_ii / pos_ii_mag # Unit vector of the inertial position

        vel_pp_mag = norm(vel_pp)
        vel_pp_hat = vel_pp / vel_pp_mag

        # Orbital Elements
        OE = rvtoorbitalelement(pos_ii, vel_ii, mass, m.planet)
        vi = OE[6]

        Mars_Gram_recalled_at_periapsis = false

        if vi > 0 && vi < pi/2 && config.ascending_phase == false
            config.ascending_phase = true
        elseif vi >= pi/2 && vi < pi && config.ascending_phase == true && args[:body_shape] == "Blunted Cone"
            config.ascending_phase = false
        end

        if config.ascending_phase == true && config.MarsGram_recall == false
            config.atmospheric_data = Dict()
        end

        # Angular Momentum Calculations 
        h_ii = cross(pos_ii, vel_ii)    # Inertial angular momentum vector [m ^ 2 / s]

        index = 1
        for item in h_ii
            if item == -0
                h_ii[index] = 0
                index = index + 1
            end
        end

        h_ii_mag = norm(h_ii)           # Magnitude of the inertial angular momentum [m ^ 2 / s]
        h_pp = cross(pos_pp, vel_pp)

        index = 1
        for item in h_pp
            if item == -0
                h_pp[index] = 0.0
                index = index + 1
            end
        end
        
        h_pp_mag = norm(h_pp)
        h_pp_hat = h_pp / h_pp_mag

        # Inertial flight path angle 
        arg = median([-1, 1, h_pp_mag / (pos_pp_mag * vel_pp_mag)])     # limit to[-1, 1]
        γ_ii = acos(arg)    
        if dot(pos_ii, vel_ii) < 0
            γ_ii = -γ_ii
        end

        # Relative flight path angle
        arg = median([-1, 1, h_pp_mag / (pos_pp_mag * vel_pp_mag)])     # limit to[-1, 1]
        γ_pp = acos(arg)
        if dot(pos_pp, vel_pp) < 0
            γ_pp = -γ_pp
        end

        ## Derived Quantity Calculations

        # Compute latitude and longitude
        LatLong = rtolatlong(pos_pp, m.planet)
        lat = LatLong[2]
        lon = LatLong[3]
        alt = LatLong[1]

        if aerobraking_phase == 2 || aerobraking_phase == 0
            if (pos_ii_mag - m.planet.Rp_e - args[:EI] * 1e3) <= 0 && config.drag_state == false && config.ascending_phase == false
                config.drag_state
                config.time_IEI = t0
            elseif (pos_ii_mag - m.planet.Rp_e >= args[:EI] * 1e3) && config.drag_state == true && config.ascending_phase
                config.drag_state = false
                config.time_OEI = t0
            end
        end

        if aerobraking_phase == 2 || aerobraking_phase == 0
            if args[:control_mode] == 1
                x = 120
            else
                x=140
            end

            if (config.hear_rate_prev > 0.005 || abs(pos_ii_mag - m.planet.Rp_e <= x*1e3)) && config.sensible_loads == false && config.ascending_phase == false
                config.sensible_loads = true
            elseif config.hear_rate_prev > 0.005 && config.sensible_loads == true && config.ascending_phase
                config.sensible_loads = false
            end
        end

        # Compute NED basis unit vectors
        uDuNuE = latlongtoNED(LatLong)
        uD = uDuNuE[1]
        uE = uDuNuE[3]
        uN = uDuNuE[2]

        # copmute azimuth
        vN = dot(vel_pp, uN)
        vE = dot(vel_pp, uE)
        azi_pp = atan(vE, vN)

        # Get density, pressure , temperature and winds
        config.MarsGram_justrecalled = 0
        if ip.dm == 0
            ρ, T_p, wind = density_constant(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 1
            ρ, T_p, wind = density_exp(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 2
            ρ, T_p, wind = density_no(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 3
            ρ, T_p, wind = marsgram(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        end

        # Define output.txt containing density datap = 0.0
        if args[:body_shape] == "Spacecraft"
            length_car = m.body.length_SA + m.body.length_SC
        elseif args[:body_shape] == "Blunted Cone"
            length_car = m.body.BaseRadius * 2
        end

        Re = vel_pp_mag * ρ * length_car / μ_fluid  # Reynolds number

        # Mach Number
        sound_velocity = sqrt(γ * m.planet.R * T_p)
        Mach = vel_pp_mag / sound_velocity
        S = sqrt(γ/2) * Mach    # Molecular speed ratio
        heat_load = in_cond[8]

        if config.drag_state == true
            ## Check type of fluid
            Kn = 1.26 * sqrt(γ) * Mach / (Re*1e-5)
            if index_phase_aerobraking == 2
                if (alt < 80000) && (config.index_warning_alt == 0)
                    println("WARNING: Altitude < 80 km!")
                end

                config.index_warning_alt = 1
            elseif alt > 100000
                config.index_warning_alt = 0
            end

            if Kn < 0.1 && config.index_warning_flow == 0
                if args[:print_res]
                    println("WARNING: Transitional flow passage!")
                end
                
                config.index_warning_flow = 1
            elseif Kn >= 0.1
                config.index_warning_flow = 0
            end
        end

        config.heat_load_past = heat_load
        # Heat rate and Control
        if (index_phase_aserobraking == 2 || index_phase_aerobraking == 1.75 || index_phase_aerobraking == 2.25) && config.drag_state && config.initial_position_closed_form
            # evaluates the closed form solution the first time at EI km
            if abs(pos_ii_mag - m.planet.Rp_e - args.EI * 1e3) <= 1e-2 && (args.control_mode == 2 || args.control_mode == 3) && config.time_switch_1 == 0
                if ip.cm == 3

                elseif ip.cm == 2

                elseif ip.cm == 1

                elseif ip.cm == 0

                end
            end

            if config.MarsGram_justrecalled == true && config.index_Mars_Gram_call != 1
                if ip.cm == 3
                    # config.α =
                elseif ip.cm == 2
                    # config.α =
                elseif ip.cm == 1
                    # config.α =
                elseif ip.cm == 0
                    # config.α =
                end
                
                append!(config.state_flesh1, [T_p, ρ, S])

                if ip.tm == 1
                    heat_rate = heatrate_convective_radiative(S, T_p, m, ρ, vel_pp_mag, config.α)
                elseif ip.tm == 2
                    heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, config.α)
                end
            end

            if indec_phase_aerobraking == 2
                if args[:control_in_loop]
                    config.state_flesh1 = [T_p, ρ, S]
                    if ip.cm == 3
                        # config.α =
                    elseif ip.cm == 2
                        # config.α =
                    elseif ip.cm == 1
                        # config.α =
                    elseif ip.cm == 0
                        # config.α =
                    end
                elseif args[:control_in_loop] == false && args[:integrator] == "Julia"
                    if controller.count_controller != controller.count_prev_controller && controller.stored_state == 0 && t0 != controller.prev_time
                        append!(config.state_flesh1, [T_p, ρ, S])

                        if controller.count_controller == 2
                            state = config.state_flesh1[end]
                        else
                            state = config.state_flesh1[end-1]
                            deleteat!(config.state_flesh1, 1)
                        end

                        controller.stored_state = 1
                        controller.prev_time = time_0

                        if ip.cm == 3
                            # config.α =
                        elseif ip.cm == 2
                            # config.α =
                        elseif ip.cm == 1
                            # config.α =
                        elseif ip.cm == 0
                            # config.α =
                        end
                    end
                end
            end

            # Heat Rate 
            if ip.tm == 1
                heat_rate = heatrate_convective_radiative(S, T_p, m, ρ, vel_pp_mag, config.α)
            elseif ip.tm == 2
                heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, config.α)
            end
            
            cp = m.planet.γ / (m.planet.γ - 1) * m.planet.R

            T_r = 0.0
        else
            T_r = 0.0
            heat_rate = 0.0
        end

        α = config.α
        config.hear_rate_prev = heat_rate # save current heat rate

        # Convert wind to pp(PCPF) frame
        wE = wind[1]
        wN = wind[2]
        wU = wind[3]

        wind_pp = wN * uN + wE * uE - wU * uD         # wind velocity in pp frame, m / s 
        vel_pp_rw = vel_pp + wind_pp                  # relative wind vector, m / s
        vel_pp_rw_hat = vel_pp_rw / norm(vel_pp_rw)   # relative wind unit vector 

        # Dynamic Pressure, CHANGE THE VELOCITY WITH THE WIND VELOCITY
        q = 0.5 * ρ * norm(vel_pp_rw)^2               # dynamic pressure based on wind, Pa
        
        ## Rotation Calculation
        rot_angle = norm(ω_planet) * (t0 + t_prev)    # angle of rotation, rad
        L_PI = [[cos(rot_angle), sin(rot_angle), 0.0], [-sin(rot_angle), cos(rot_angle), 0.0], [0.0, 0.0, 1.0]] # rotation matrix

        L_PI = [[x + 0.0 for x in row] for row in L_PI]

        if ip.gm == 0
            gravity_ii = gravity_const(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        elseif ip.gm == 1
            gravity_ii = gravity_invsquared(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        elseif ip.gm == 2
            gravity_ii = gravity_invsquared_J2(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        end

        bank_angle = 0.0

        lift_pp_hat = cross(h_pp_hat, vel_pp_rw_hat)

        # Vehicle Aerodynamic Forces
        # CL and CD
        if ip.am == 0
            [CL, CD] = aerodynamic_coefficient_constant(α, m.body, T_p, S, args, MonteCarlo)
        elseif ip.am == 1
            [CL, CD] = aerodynamic_coefficent_fM(α, m.body, T_p, S, args, MonteCarlo)
        elseif ip.am == 2
            [CL, CD] = aerodynamic_coefficient_no_ballistic_flight(α, m.body, T_p, S, args, MonteCarlo)
        end

        β = mass / (CD*area_tot)    # ballistic coefficient, kg / m ^ 2

        # Force Calculation
        drag_pp_hat = -vel_pp_rw_hat    # PLanet relative drag force direction

        drag_pp = q * CD * area_tot * drag_pp_hat                       # PLanet relative drag force vector
        lift_pp = q * CL * area_tot * lift_pp_hat * cos(bank_angle)     # PLanet relative lift force vector

        drag_ii = dot(L_PI', drag_pp)   # Inertial drag force vector
        lift_ii = dot(L_PI', lift_pp)   # Inertial lift force vector

        # Check if propellant mass is greater than 0 kg
        if config.index_propellant_mass == 1
            if mass - args[:dry_mass] <= 0.5
                config.index_propellant_mass = 0
                m.engine.T = 0

                if args[:print_res]
                    println("WARNING: No fuel left!")
                end
            end
        end

        # Thrust
        Δv = m.engine.g_e * m.engine.Isp * log(initial_state.m / mass)

        if ip.tc == 0
            # thrust_pp_mag =
        elseif ip.tc == 1
            # thrust_pp_mag =
        elseif ip.tc == 2
            # thrust_pp_mag =
        end

        # Rodrigues rotation formula to rotate thrust vector of angle phi around angular vector from D direction
        D_L_per_pp_hat = cross(drag_pp_hat, lift_pp_hat)
        thrust_pp_hat =  drag_pp_hat * cos(args[:phi]) + cross(D_L_per_pp_hat, drag_pp_hat) * sin(args[:phi]) + D_L_per_pp_hat * dot(D_L_per_pp_hat, drag_pp_hat) * (1 - cos(args[:phi]))
        #these two ways give the same direction
        thrust_pp = thrust_pp_mag * thrust_pp_hat
        thrust_ii = dot(L_PI', thrust_pp)

        # Total Force
        # Total inetrial external force vector on body [N]
        force_ii = draag_ii + lift_ii + gravity_ii + thrust_ii

        # index_steps_EOM
        y_dot = zeros(8)
        
        y_dot[1:3] = vel_ii
        y_dot[4:6] = force_ii / mass
        y_dot[7] = -norm(thrust_ii) / (m.engine.g_e * m.engine.Isp) # mass variation
        y_dot[8] = heat_rate
        energy = (vel_ii_mag^2)*0.5 - m.planet.μ / pos_ii_mag

        ## SAVE RESULTS
        if config.results_save
            sol = hcat([[t0], [timereal.year], [timereal.month], [timereal.day], [timereal.hour], [timereal.minute],
                        [timereal.second], [numberofpassage], pos_ii, vel_ii, [pos_ii_mag], [vel_ii_mag], pos_pp, 
                        [pos_pp_mag], vel_pp, [vel_pp_mag], [OE[1]], [OE[2]], [OE[3]], [OE[4]], [OE[5]], [OE[6]],
                        [lat], [lon], [alt], [γ_ii], [γ_pp], h_ii, h_pp, [h_ii_mag], [h_pp_mag], uD, uE, uN, [vN], [vE],
                        [azi_pp], [ρ], [T_p], [p], wind, [CL], [CD], [α], [S], [mass], [heat_rate], [heat_load], [T_r], 
                        [q], gravity_ii, drag_pp, drag_ii, lift_pp, lift_ii, force_ii, [energy], config.index_MonteCarlo, Int64(config.drag_state)])
            
            sol = reshape(sol, (1, 91))

            append!(config.solution_intermediate, sol)
        end

        return y_dot
    end

    ## EVENTS 
    function eventfirststep(t, y)
        return (sqrt(y[1]^2 + y[2]^2 + y[3]^2) - m.planet.Rp_e - 250*1e3)
    end
    eventfirststep.direction = -1
    eventfirststep.terminal = true

    function eventsecondstep(t, y)
        return (sqrt(y[1]^2 + y[2]^2 + y[3]^2) - m.planet.Rp_e - 260*1e3)
    end
    eventsecondstep.terminal = true
    eventsecondstep.direcion = 1

    function reached_EI(t, y)
        return (sqrt(y[1]^2 + y[2]^2 + y[3]^2) - m.planet.Rp_e - args[:EI]*1e3)
    end
    eventsecondstep.terminal = false
    eventsecondstep.direcion = -1

    function reached_AE(t, y)
        return (sqrt(y[1]^2 + y[2]^2 + y[3]^2) - m.planet.Rp_e - args[:AE]*1e3)
    end
    eventsecondstep.terminal = false
    eventsecondstep.direcion = 1

    function out_drag_passage(t, y)
        if abs((sqrt(y[1]^2 + y[2]^2 + y[3]^2) - m.planet.Rp_e - args[:AE]*1e3)) <= 1e-5
            if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 2
                config.α = m.aerodynamics.α
            elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 3
                config.α = 0.0
            end
        end

        return (sqrt(y[1]^2 + y[2]^2 + y[3]^2) - m.planet.Rp_e - args[:AE]*1e3)
    end
    out_drag_passage.terminal = true
    out_drag_passage.direction = 1

    function in_drag_passage(t, y)
        pos_ii = [y[1], y[2], y[3]]  # Inertial position
        vel_ii = [y[4], y[5], y[6]]  # Inertial velocity
        [pos_pp, vel_pp] = r_intor_p(pos_ii, vel_ii, m.planet, t, t_prev)

        LatLong = rtolatlong(pos_pp, m.planet)

        h0 = LatLong[1]

        if ip.gm == 2
            cond = h0 - args[:EI] * 1e3
            thr = 500
        else
            cond = sqrt(y[1]^2 + y[2]^2 + y[3]^2) - m.planet.Rp_e - args[:EI]*1e3
            thr = 1e-5
        end

        if abs(cond) <= thr
            controller.guidance_t_eval = range(t, t+1500, 1/args[:flash1_rate])

            # State definition for control 2, 3 State used by closed-form solution
            pos_ii = [y[1], y[2], y[3]]   # Inertial position
            vel_ii = [y[4], y[5], y[6]]   # Inertial velocity
            OE_closedform = rvtoorbitalelement(pos_ii, vel_ii, y[7], m.planet)
            config.inital_position_closed_form = OE_closedform
        end

        return cond
    end
    in_drag_passage.terminal = true
    in_drag_passage.direction = -1

    function in_drag_passage_nt(t, y)
        pos_ii = [y[1], y[2], y[3]]  # Inertial position
        vel_ii = [y[4], y[5], y[6]]  # Inertial velocity
        [pos_pp, vel_pp] = r_intor_p(pos_ii, vel_ii, m.planet, t, t_prev)

        LatLong = rtolatlong(pos_pp, m.planet)

        h0 = LatLong[1]

        if ip.gm == 2 && !args[:drag_passage]
            cond = (h0 - args[:EI]*1e3)
            thr = 500
        else
            cond = (sqrt(y[1]^2 + y[2]^2 + y[3]^2) - m.planet.Rp_e - args[:EI]*1e3)
            thr = 1e-3
        end

        if abs(cond) <= thr && !config.initial_position_closed_form
            controller.guidance_t_eval = range(t, t+1500, 1/args[:flash1_rate])

            # State definition for control 2, 3 State used by closed-form solution
            OE_closedform = rvtoorbitalelement(pos_ii, vel_ii, y[7], m.planet)
            config.inital_position_closed_form = OE_closedform
        end

        return (sqrt(y[1]^2 + y[2]^2 + y[3]^2) - m.planet.Rp_e - args[:EI]*1e3)
    end
    in_drag_passage_nt.terminal = false
    in_drag_passage_nt.direction = -1

    function apoapsispoint(t, y)
        pos_ii = [y[1], y[2], y[3]]   # Inertial position
        vel_ii = [y[4], y[5], y[6]]   # Inertial Velocity

        vi = rvtoorbitalelement(pos_ii, vel_ii, y[7], m.planet)[6]

        return rad2deg(vi) - 180
    end
    apoapsispoint.terminal = true
    apoapsispoint.direction = 1

    function periapsispoint(t, y)
        pos_ii = [y[1], y[2], y[3]]   # Inertial position
        vel_ii = [y[4], y[5], y[6]]   # Inertial Velocity

        vi = rvtoorbitalelement(pos_ii, vel_ii, y[7], m.planet)[6]

        return vi
    end
    periapsispoint.terminal = false
    periapsispoint.direction = 1

    function impact(t, y)
        if args[:body_shape] == "Blunted Body"
            min_alt = 1 * 10^3
        else
            min_alt = 35 * 10^3
        end

        return (sqrt(y[1]^2 + y[2]^2 + y[3]^2) - (m.planet.Rp_e + min_alt))
    end
    impact.terminal = true

    function apoapsisgreaterperiapsis(t, y)
        y = [x + 0 for x in y]
        r = y[1:3]
        v = y[4:6]
        Energy = norm(v)^2 * 0.5 - m.planet.μ / norm(r)
        a = -m.planet.μ / (2 * Energy)
        h = cross(r, v)
        h += 0
        e = sqrt(1 + (2 * Energy * dot(h,h)/ (m.planet.μ)^2))

        r_a = a * (1 + e)
        r_p = a * (1 - e)

        if r_a < r_p && args[:body_shape] == "Spacecraft"
            println("Periapsis greater than apoapsis!")
        end

        return r_a - r_p 
    end
    apoapsisgreaterperiapsis.terminal = true

    function stop_firing(t, y)
        mass = y[7]
        Δv = (m.engine.g_e * m.engine.Isp) * log(initial_state.m/mass)
        
        return Δv - args[:delta_v]
    end
    stop_firing.terminal = true

    function guidance(t, y)

    end
    guidance.terminal = false
    guidance.direction = 1

    function heat_rate_check(t, y)
        if args[:control_mode] == 1
            x = 120
        else
            x = 160
        end

        if abs(sqrt(y[1]^2 + y[2]^2 + y[3]^2) - m.planet.Rp_e - x*1e3) <= 1e-5
            controller.guidance_t_eval = range(t, t+2500, 1/args[:flash1_rate])
        end

        return (sqrt(y[1]^2 + y[2]^2 + y[3]^2) - m.planet.Rp_e - x*1e3)
    end
    heat_rate_check.terminal = true

    function heat_load_check_exit(t, y)
        if args[:control_mode] == 1
            x = 120
        else
            x = 160
        end

        return (sqrt(y[1]^2 + y[2]^2 + y[3]^2) - m.planet.Rp_e - x*1e3)
    end
    heat_load_check_exit.terminal = true
    heat_load_check_exit.direction = 1

    time_0 = 0
    if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 2
        config.α = m.aerodynamics.α
    elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 3
        config.α = 0.0
    end

    stop_simulation = false
    save_pre_index = 0
    save_post_index = 0

    config.impact = false
    config.solution_intermediate = []
    config.count_dori = 0
    config.atmospheric_data = Dict()
    config.previous_atmospheric_data = Dict()
    config.ascending_phase = false
    config.evaluate_switch_heat_load = false
    config.state_inner_boundary_atmosphere = [] # used in Density model for vi def
    config.time_IP = config.time_OP
    config.time_IEI = 0
    config.time_OEI = 0
    config.time_switch_1 = 0.0
    config.time_switch_2 = 0.0

    if args[:heat_load_sol] == 2
        config.time_switch_2 = 1000.0
    end

    config.timer_revaluation = 0
    config.closed_form_solution_off = 1         # used in closed form solution online to run the solution only once
    config.initial_position_closed_form = []
    config.heat_rate_list = []                  # checked this - if used
    config.aoa_list = []                        # checked this - if used
    controller.guidance_t_eval = []
    controller.count_prev_controller = 0
    controller.count_controller = 1
    controller.stored_state = 1
    controller.prev_time = 0
    controller.t = 0

    config.security_mode = False
    config.stop_simulation=False
    config.results_save = 1
    config.drag_state = False
    config.α_past = m.aerodynamics.α

    if norm(r0) - m.planet.Rp_e <= args[:EI]*1e3
        config.drag_state = True
        config.initial_position_closed_form = OE
    end

    config.sensible_loads = false
    config.counter_integrator = 0

    # Def initial conditions
    in_cond = [r0[1], r0[2], r0[3], v0[1], v0[2], v0[3], mass, 0.0]

    # If aerobraking maneuver allowed, add a prephase 0
    range_phase_i = 1
    if args[:thrust_control] == "Aerobraking Maneuver"
        range_phase_i = 0
    end

    # Solve Equations of Motion 
    for aerobraking_phase in range(range_phase_i, 4)
        index_phase_aerobraking = aerobraking_phase

        if (index_steps_EOM == 1 || args[:drag_passage]) && (aerobraking_phase == 2 || aerobraking_phase == 4 || aerobraking_phase == 1)
            continue
        end

        # Definition of eventsecondstep
        if aerobraking_phase == 1
            events = [stop_firing, apoapsisgreaterperiapsis, impact]
        elseif aerobraking_phase == 2
            events = [eventfirststep, apoapsisgreaterperiapsis. impact]
        elseif aerobraking_phase == 3 && args[:drag_passage]
            events = [out_drag_passage, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, impact]
        elseif aerobraking_phase == 3 && index_steps_EOM == 1 && args[:body_shape] == "Blunted Cone"
            events = [out_drag_passage, apoapsispoint, periapsispoint, impact]
        elseif aerobraking_phase == 3 && index_steps_EOM == 1 && args[:drag_passage] == false
            events = [apoapsispoint, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, impact]
        elseif aerobraking_phase == 3
            events = [eventsecondstep, periapsispoint, reached_EI, reached_AE, in_drag_passage, apoapsisgreaterperiapsis, impact]
        elseif aerobraking_phase == 4
            events = [apoapsispoint, periapsispoint, apoapsisgreaterperiapsis, impact]
        end

        if index_phase_aerobraking == 3
            save_pre_index = length(solution.orientation.time)
            simulator = args[:integrator]
        end

        # Definition dtep montecarlo_size
        if aerobraking_phase == 2 || aerobraking_phase == 4
            step = 5
            r_tol = 1e-10
            a_tol = 1e-12
            simulator = "Julia"
            method = "BDF"
            save_ratio = 5
        elseif aerobraking_phase == 1
            step = 0.1
            r_tol = 1e-9
            a_tol = 1e-11
            simulator = "Julia"
            method = "BDF"
            save_ratio = 5
        elseif aerobraking_phase == 3
            if args[:integrator] == "Julia"
                step = 0.1
                r_tol = 1e-9
                a_tol = 1e-11

                if MonteCarlo
                    method = "RK45"
                else
                    method = "BDF"
                end

                save_ratio = args[:save_rate]
            else
                step = 1/args[:trajectory_rate]
                save_ratio = round(1args[:trajectory_rate] / args[:save_rate])
            end
        end

        # Definition length simulation
        length_sim = 1e10
        i_sim = 0
        time_solution = []

        config.continue_simulation = true

        while config.continute_simulation
            index_phase_aerobraking = aerobraking_phase
            # if control mode =! 0, redefine sim setting and creates two more phases until reaching EI and out of the AE phase 2: between 120 km alt
            if aerobraking_phase == 3 && (args[:control_mode] != 0 && args[:control_in_loop] == 0 && config.drag_state == true && config.sensible_loads == true && config.ascending_phase == false)
                simulator = args[:integrator]
                events = [out_drag_passage, heat_load_check_exit, periapsispoint, apoapsisgreaterperiapsis, impact]

                if args[:integrator] == "Costumed"
                    step = 1/args[:trajectory_rate]
                else
                    method = "RK23"
                    if args[:control_in_loop]
                        step = 10/args[:trajectory_rate]
                    else
                        step = 1/(4*args[:flash1_rate])
                        r_tol = 1e-6
                        a_tol = 1e-7

                        if step >= args[:flash1_rate]
                            step = 1/(2*args[:flash1_rate])
                        end
                        insert!(events, 4, guidance)
                    end
                end

                controller.t = controller.guidance_t_eval[controller.count_controller]

            elseif aerobraking_phase == 3 && args[:control_mode] != 0 && args[:control_in_loop] == 0 && config.drag_state == false && config.ascending_phase == false
                events = [periaoapsispoint, out_drag_passage, heat_rate_check, apoapsisgreaterperiapsis, impact]
                simulator = "Julia"
                index_phase_aerobraking = 1.75
                step = 0.05
                r_tol = 1e-6
                a_tol = 1e-7
                method = "RK23"
            elseif aerobraking_phase == 3 && args[:control_mode] != 0 && args[:control_in_loop] == 0 && config.drag_state == true && config.ascending_phase == false
                events = [periapsispoint, in_drag_passage, apoapsisgreaterperiapsis, impact]
                simulator = "Julia"
                step = 0.1
                r_tol = 1e-9
                a_tol = 1e-10
                index_phase_aerobraking = 1.5
                method = "RK23"
            elseif aerobraking_phase == 3 && args[:control_mode] != 0 && args[:control_in_loop] == 0 && config.drag_state == true && config.ascending_phase == true
                events = [periapsispoint, out_drag_passage, apoapsisgreaterperiapsis, impact]
                simulator = "Julia"
                step = 0.1
                r_tol = 1e-9
                a_tol = 1e-10
                index_phase_aerobraking = 2.25
                method = "RK23"
            elseif aerobraking_phase == 3 && args[:control_mode] != 0 && args[:control_in_loop] == 0 && config.ascending_phase == true
                events = [eventsecondstep, periapsispoint, eventsecondstep, apoapsisgreaterperiapsis, impact]
                simulator = "Julia"
                index_phase_aerobraking = 2.5
                step = 0.5
                r_tol = 1e-8
                a_tol = 1e-9
                method = "RK23"
            end

            if args[:print_res]
                println("Step #", index_phase_aerobraking)
            end

            if simulator == "Julia"
                ## Julia Integrator
                # Time initialization
                initial_time, final_time = time_0, time_0 + length_sim

                # Run simulation
                eomfunction = 


                config.counter_integrator += 1
                in_cond = []

                # Save results 
                extend!(time_solution, solution.t)
                time_0 = time_solution[end]

                if aerobraking_phase == 1
                    new_periapsis(m, in_cond[1:3], in_cond[4:6], args)
                end
            elseif simulator == "Costumed"
                if args[:integrator] == "Costumed"
                    mutable struct Sol
                        t_events::Vector{Vector{Float64}}
                    end

                    sol = Sol([[],[]])
                end

                while stop_simulation == false
                    initial_time = time_0
                    ## Costumed Integrator
                    config.MarsGram_recall = 1
                    config.results_save = 0
                    y, t, stop_simulation, sol = RK4(f, step, initial_time, in_cond, m, T_ijk, index_phase_aerobraking, args, sol)

                    in_cond = [y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8]]
                    config.counter_integrator += 1
                    config.results_save = 1
                    f(t, in_cond, m, index_phase_aerobraking, ip)

                    # New initial Condition
                    time_0 = t_prev

                    # Save Results
                    append!(time_solution, t)

                    ## Guidance, Navigation and Control
                    if args[:control_mode] != 0
                        # DO LATER
                    end
                end

                i_sim += 1
            end

            # Define breaker campaign impact or apoapsis greater than periapsis
            continue_campaign = event(solution)

            if continue_campaign == false
                config.impact = true
                break
            end

            # Breaker conditions
            if simulator == "Julia"
                if length(sol.t_events[1]) != 0 || (args[:drag_passage] && index_phase_aerobraking == 2.25 && length(sol.t_events[2]) != 0)
                    config.continue_simulation = false
                    break
                end
            else
                if args[:control_mode] == 0 && stop_simulation == true
                    config.continue_simulation = false
                    break
                end
            end
        end

        # Save Results
        time_0 = save_results(time_solution, save_ratio)

        if index_phase_aerobraking == 3 || index_phase_aerobraking == 2.5 || (index_phase_aerobraking == 2.25 && args[:drag_passage])
            save_post_index = length(solution.orientation.time)
        end

        config.count_phase = 0

        # Define breaker campaign
        if continue_campaign == false
            if save_post_index == 0
                save_post_index = length(solution.orientation.time)
            end

            break
        end
    end

    # Re-set count index to 0
    config.count_dori = 0

    if args[:drag_passage] == false && pi - solution.orientation.oe[end][end] > 1e-4 && continue_campaign == true && args[:body_shape] != "Blunted Cone"
        final_conditions_notmet = true
        events = [apoapsispoint]
    else
        final_conditions_notmet = false
    end

    count_temp = 0

    while final_conditions_notmet
        in_cond = [solution.orientation.pos_ii[1][end], solution.orientation.pos_ii[2][end], solution.orientation.pos_ii[3][end], 
                   solution.orientation.vel_ii[1][end], solution.orientation.vel_ii[2][end], solution.orientation.vel_ii[3][end],
                   solution.performance.mass[end], solution.performance.heat_load[end]]
                    
        initial_time, final_time = time_0, time_0  + 100
        step = 0.005
        r_tol = 1e-12
        a_tol = 1e-13

        try
            eomfunction = 

            config.counter_integrator += 1
            time_0 = save_results()
            count_temp += 1
        catch
            break
        end

        if count_temp > 5
            break
        end

        if args[:drag_passage] == false && pi - solution.orientation.oe[end][end] > 1e-4 && continue_campaign == true
            final_conditions_notmet = true
            events = [apoapsispoint]
        else
            final_conditions_notmet = false
        end
    end

    config.save_index_heat = length(solution.orientation.time)
    config.time_OP = length(solution.orientation.time)

    append!(config.altitudeperiapsis, minimum(solution.orientation.alt[save_pre_index:save_post_index])*1e-3)
    append!(config.max_heatrate, maximum(solution.performance.heat_rate[save_pre_index:save_post_index]))
    config.delta_v_man = (m.engine.g_e * m.engine.Isp) * log(m.body.Mass / solution.performance.mass[end])

    if args[:print_res]
        try
            println()
        catch
            println("Problem in the indexes")
        end
    end

    append!(config.periapsis_list, minimum(solution.orientation.alt[save_pre_index:save_post_index])*1e-3)
    append!(config.orbit_number_list, config.count_number_of_passage + 1)
    append!(config.delta_v_list, config.delta_v_man)

    return continue_campaign
end