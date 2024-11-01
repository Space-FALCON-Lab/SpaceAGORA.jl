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
        

    end
end