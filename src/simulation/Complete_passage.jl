include("../utils/Reference_system.jl")
include("../integrator/Integrators.jl")
include("../integrator/Events.jl")
include("../utils/Save_results.jl")

include("../physical_models/Gravity_models.jl")
include("../physical_models/Density_models.jl")
include("../physical_models/Aerodynamic_models.jl")
include("../physical_models/Thermal_models.jl")
include("../physical_models/Perturbations.jl")

include("../control/Control.jl")
include("../control/Propulsive_maneuvers.jl")

using LinearAlgebra
using DifferentialEquations
using Dates
using AstroTime
using SPICE
using PythonCall
using StaticArrays
sys = pyimport("sys")

import .config
import .ref_sys

function asim(ip, m, initial_state, numberofpassage, args, gram_atmosphere=nothing, gram=nothing)
    wind_m = false
    if ip.wm == 1
        wind_m = true
    end

    MonteCarlo = false
    if ip.mc == 1
        MonteCarlo = true
    end

    OE = SVector{7, Float64}([initial_state.a, initial_state.e, initial_state.i, initial_state.Ω, initial_state.ω, initial_state.vi, initial_state.m])

    if (OE[1] > (m.planet.Rp_e*1e-3 + args[:EI])*1e3) && (args[:drag_passage] == false) && (args[:body_shape] == "Spacecraft")
        index_steps_EOM = 3
    else
        index_steps_EOM = 1
    end

    i = OE[3]
    Ω = OE[4]
    ω = OE[5]

    T_ijk = [cos(Ω)*cos(ω)-sin(Ω)*sin(ω)*cos(i)   sin(Ω)*cos(ω)+cos(Ω)*sin(ω)*cos(i)    sin(ω)*sin(i);
             -cos(Ω)*sin(ω)-sin(Ω)*cos(ω)*cos(i)  -sin(Ω)*sin(ω)+cos(Ω)*cos(ω)*cos(i)   cos(ω)*sin(i);
             sin(Ω)*sin(i)                        -cos(Ω)*sin(i)                        cos(i)]


    r0, v0 = orbitalelemtorv(OE, m.planet)

    Mass = OE[end]

    # Clock
    date_initial = from_utc(DateTime(m.initial_condition.year, 
                                    m.initial_condition.month,
                                    m.initial_condition.day, 
                                    m.initial_condition.hour, 
                                    m.initial_condition.minute, 
                                    m.initial_condition.second))

    config.cnf.count_numberofpassage = config.cnf.count_numberofpassage + 1

    if config.cnf.count_numberofpassage != 1
        t_prev = config.solution.orientation.time[end]
    else
        t_prev = m.initial_condition.time_rot # value(seconds(date_initial - from_utc(DateTime(2000, 1, 1, 12, 0, 0)))) # m.initial_condition.time_rot
    end

    function f!(y_dot, in_cond, param, t0)
        m = param[1]
        index_phase_aerobraking = param[2]
        ip = param[3]
        aerobraking_phase = param[4]
        t_prev = param[5]
        date_initial = param[6]
        time_0 = param[7]
        args = param[8]
        initial_state = param[9]
        gram_atmosphere = param[10]
        gram = param[11]

        
        ## Counters
        # Counter for all along the simulation of all passages
        config.cnf.count_aerobraking = config.cnf.count_aerobraking + 1
        passage_number = config.cnf.count_aerobraking
        # Counter for one entire passage
        config.cnf.count_dori = config.cnf.count_dori + 1
        # Counter for one phase
        config.cnf.count_phase = config.cnf.count_phase + 1

        t0 = t0 * config.cnf.TU

        # Clock
        current_epoch = date_initial + t0*seconds # Precompute the current epoch
        time_real = DateTime(current_epoch) # date_initial + Second(t0)
        timereal = ref_sys.clock(Dates.year(time_real), Dates.month(time_real), Dates.day(time_real), Dates.hour(time_real), Dates.minute(time_real), Dates.second(time_real))

        # Timing variables
        el_time = value(seconds(current_epoch - m.initial_condition.DateTimeIC)) # Elapsed time since the beginning of the simulation
        current_time =  value(seconds(current_epoch - m.initial_condition.DateTimeJ2000)) # current time in seconds since J2000
        time_real_utc = to_utc(time_real) # Current time in UTC as a DateTime object
        config.cnf.et = utc2et(time_real_utc) # Current time in Ephemeris Time
        m.planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_"*uppercase(m.planet.name), config.cnf.et))*m.planet.J2000_to_pci' # Construct a rotation matrix from J2000 (Planet-fixed frame 0.0 seconds past the J2000 epoch) to planet-fixed frame
        

        # Assign state
        pos_ii = SVector{3, Float64}(in_cond[1:3] * config.cnf.DU)                      # Inertial position 
        vel_ii = SVector{3, Float64}(in_cond[4:6] * config.cnf.DU / config.cnf.TU)      # Inertial velocity
        mass = in_cond[7] * config.cnf.MU                                          # Mass kg
        pos_ii_mag = norm(pos_ii)                                  # Magnitude of the inertial position
        vel_ii_mag = norm(vel_ii)                                  # Magnitude of the inertial velocity

        # Assign parameters
        ω_planet = m.planet.ω
        γ = m.planet.γ
        μ_fluid = m.planet.μ_fluid
        area_tot = m.body.area_tot

        # TRANSFORM THE STATE
        # Inertial to planet relative transformation
        pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, config.cnf.et) # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
        pos_pp_mag = norm(pos_pp) # Magnitude of the planet relative position
        pos_pp_hat = pos_pp / pos_pp_mag # Unit vector of the planet relative position
        pos_ii_hat = pos_ii / pos_ii_mag # Unit vector of the inertial position

        vel_pp_mag = norm(vel_pp)
        vel_pp_hat = vel_pp / vel_pp_mag

        # Orbital Elements
        # println("pos_ii: ", pos_ii, " vel_ii: ", vel_ii, " mass: ", mass)
        # # sleep(5.0)
        OE = rvtoorbitalelement(pos_ii, vel_ii, mass, m.planet)
        vi = OE[6]

        
        # m.planet.L_PI .= pxform("J2000", "IAU_"*uppercase(m.planet.name), current_time)*m.planet.J2000_to_pci' # Construct a rotation matrix from J2000 (Planet-fixed frame 0.0 seconds past the J2000 epoch) to planet-fixed frame
        # println("L_PI: ", m.planet.L_PI)
        # # sleep(5.0)
        Mars_Gram_recalled_at_periapsis = false

        if vi > 0 && vi < pi/2 && config.cnf.ascending_phase == false
            config.cnf.ascending_phase = true
        elseif vi >= pi/2 && vi <= pi && config.cnf.ascending_phase == true && args[:body_shape] == "Blunted Cone"
            config.cnf.ascending_phase = false
        end

        if config.cnf.ascending_phase == true && config.cnf.MarsGram_recall == false
            config.cnf.atmospheric_data = Dict()
        end

        # Angular Momentum Calculations 
        h_ii = cross(pos_ii, vel_ii)    # Inertial angular momentum vector [m ^ 2 / s]

        h_ii_mag = norm(h_ii)           # Magnitude of the inertial angular momentum [m ^ 2 / s]
        h_pp = cross(pos_pp, vel_pp)
        
        h_pp_mag = norm(h_pp)
        h_pp_hat = h_pp / h_pp_mag

        # Inertial flight path angle 
        arg = median([-1, 1, h_ii_mag / (pos_ii_mag * vel_ii_mag)])     # limit to[-1, 1]
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
        LatLong = rtolatlong(pos_pp, m.planet, args[:topography_model] == "Spherical Harmonics" && norm(pos_ii) - m.planet.Rp_e < args[:EI] * 1e3)
        
        lat = LatLong[2]
        lon = LatLong[3]
        alt = LatLong[1]
        
        # println(" ")
        # println(" Altitude: ", alt)
        # println(" ")

        if aerobraking_phase == 2 || aerobraking_phase == 0
            if (pos_ii_mag - m.planet.Rp_e - args[:EI] * 1e3) <= 0 && config.cnf.drag_state == false && config.cnf.ascending_phase == false
                config.cnf.drag_state = true
                config.cnf.time_IEI = t0
            elseif (pos_ii_mag - m.planet.Rp_e >= args[:EI] * 1e3) && config.cnf.drag_state == true && config.cnf.ascending_phase
                config.cnf.drag_state = false
                config.cnf.time_OEI = t0
            end
        end

        if aerobraking_phase == 2 || aerobraking_phase == 0
            if args[:control_mode] == 1
                x = 120
            else
                x = 140
            end

            if (config.cnf.heat_rate_prev > 0.005 || abs(pos_ii_mag - m.planet.Rp_e <= x*1e3)) && config.cnf.sensible_loads == false && config.cnf.ascending_phase == false
                config.cnf.sensible_loads = true
            elseif config.cnf.heat_rate_prev > 0.005 && config.cnf.sensible_loads == true && config.cnf.ascending_phase
                config.cnf.sensible_loads = false
            end
        end

        # Compute NED basis unit vectors
        uD, uN, uE = latlongtoNED(LatLong)

        # copmute azimuth
        vN = dot(vel_pp, uN)
        vE = dot(vel_pp, uE)
        azi_pp = atan(vE, vN)

        # Get density, pressure , temperature and winds
        config.cnf.Gram_justrecalled = 0
        if ip.dm == 0
            ρ, T_p, wind = density_constant(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 1
            ρ, T_p, wind = density_exp(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 2
            ρ, T_p, wind = density_no(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 3
            ρ, T_p, wind = density_gram(alt, m.planet, lat, lon, MonteCarlo, wind_m, args, el_time, gram_atmosphere, gram)
            ρ, T_p, wind = pyconvert(Float64, ρ), pyconvert(Float32, T_p), SVector{3, Float32}([pyconvert(Float32, wind[1]), pyconvert(Float32, wind[2]), pyconvert(Float32, wind[3])])
        elseif ip.dm == 4
            ρ, T_p, wind = density_nrlmsise(alt, m.planet, lat, lon, MonteCarlo, wind_m, args, time_real)
        end

        # Define output.txt containing density data
        p = 0.0
        if args[:body_shape] == "Spacecraft"
            length_car = m.body.length_SA + m.body.length_SC
        elseif args[:body_shape] == "Blunted Cone"
            length_car = m.body.base_radius * 2
        end

        Re = vel_pp_mag * ρ * length_car / μ_fluid  # Reynolds number

        # Mach Number
        sound_velocity = sqrt(γ * m.planet.R * T_p)
        Mach = vel_pp_mag / sound_velocity
        S = sqrt(γ/2) * Mach    # Molecular speed ratio
        heat_load = in_cond[8] * config.cnf.MU / config.cnf.TU^2 # * 1e4

        if config.cnf.drag_state == true
            ## Check type of fluid and check if this changes for different planets
            Kn = 1.26 * sqrt(γ) * Mach / (Re + 1e-5)
            if index_phase_aerobraking == 2
                if (alt < 80000) && (config.cnf.index_warning_alt == 0)
                    println("WARNING: Altitude < 80 km!")
                end

                config.cnf.index_warning_alt = 1
            elseif alt > 100000
                config.cnf.index_warning_alt = 0
            end

            if Kn < 0.1 && config.cnf.index_warning_flow == 0
                if Bool(args[:print_res])
                    println("WARNING: Transitional flow passage!")
                end
                
                config.cnf.index_warning_flow = 1
            elseif Kn >= 0.1
                config.cnf.index_warning_flow = 0
            end
        end

        config.cnf.heat_load_past = heat_load

        # Heat rate and Control
        if (index_phase_aerobraking == 2 || index_phase_aerobraking == 1.75 || index_phase_aerobraking == 2.25) && config.cnf.drag_state && length(config.cnf.initial_position_closed_form) != 0
            # evaluates the closed form solution the first time at EI km
            if abs(pos_ii_mag - m.planet.Rp_e - args[:EI] * 1e3) <= 1e-2 && (args[:control_mode] == 2 || args[:control_mode] == 3) && config.cnf.time_switch_1 == 0
                if ip.cm == 3
                    control_solarpanels_openloop(ip, m, args, [1,0], [T_p, ρ, S], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, 0, true, gram_atmosphere)
                elseif ip.cm == 2
                    control_solarpanels_heatload(ip, m, args, [1,0], [T_p, ρ, S], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, 0, gram_atmosphere)
                elseif ip.cm == 1
                    control_solarpanels_heatrate(ip, m, args, [1,0], [T_p, ρ, S], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form)
                elseif ip.cm == 0
                    no_control(ip, m, args, [1,0], [T_p, ρ, S], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form)
                end
            end

            # if config.cnf.Gram_justrecalled == true && config.cnf.index_Mars_Gram_call != 1  # in MC, when we reavaluate Mars Gram there is a discontinuity with the density which is created by how the density data are created. This discontinuity create really high peaks. We recalculate aoa for the new density data.
            #     if ip.cm == 3
            #         config.cnf.α = control_solarpanels_openloop(ip, m, args, [1,0], [T_p, ρ, S], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form)
            #     elseif ip.cm == 2
            #         config.cnf.α = control_solarpanels_heatload(ip, m, args, [1,0], [T_p, ρ, S], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form)
            #     elseif ip.cm == 1
            #         config.cnf.α = control_solarpanels_heatrate(ip, m, args, [1,0], [T_p, ρ, S], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form)
            #     elseif ip.cm == 0
            #         config.cnf.α = no_control(ip, m, args, [1,0], [T_p, ρ, S], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form)
            #     end
                
            #     push!(config.cnf.state_flesh1, [T_p, ρ, S])

            #     if ip.tm == 1
            #         heat_rate = heatrate_convective_radiative(S, T_p, m, ρ, vel_pp_mag, config.cnf.α)
            #     elseif ip.tm == 2
            #         heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, config.cnf.α)
            #     end
            # end

            if index_phase_aerobraking == 2
                if Bool(args[:control_in_loop])
                    config.cnf.state_flesh1 = [[T_p, ρ, S]]
                    if ip.cm == 3
                        config.cnf.α = control_solarpanels_openloop(ip, m, args, [1,1], config.cnf.state_flesh1[1], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE, true, gram_atmosphere)
                    elseif ip.cm == 2
                        config.cnf.α = control_solarpanels_heatload(ip, m, args, [1,1], config.cnf.state_flesh1[1], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE, gram_atmosphere)
                        # println("control_solarpanels_heatload: ", config.cnf.α)
                    elseif ip.cm == 1
                        config.cnf.α = control_solarpanels_heatrate(ip, m, args, [1,1], config.cnf.state_flesh1[1], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE)
                    elseif ip.cm == 0
                        config.cnf.α = no_control(ip, m, args, [1,1], config.cnf.state_flesh1[1], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE)
                    end
                elseif args[:control_in_loop] == false && args[:integrator] == "Julia"
                    if config.controller.count_controller != config.controller.count_prev_controller && config.controller.stored_state == 0 && t0 != config.controller.prev_time
                        push!(config.cnf.state_flesh1, [T_p, ρ, S]) # might have to change to push!

                        if config.controller.count_controller == 2
                            state = config.cnf.state_flesh1[end]
                        else
                            state = config.cnf.state_flesh1[end-1]
                            deleteat!(config.cnf.state_flesh1, 1)
                        end

                        config.controller.stored_state = 1
                        config.controller.prev_time = time_0

                        if ip.cm == 3
                            config.cnf.α = control_solarpanels_openloop(ip, m, args, [1,1], state, t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE, true, gram_atmosphere)
                        elseif ip.cm == 2
                            config.cnf.α = control_solarpanels_heatload(ip, m, args, [1,1], state, t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE, gram_atmosphere)
                        elseif ip.cm == 1
                            config.cnf.α = control_solarpanels_heatrate(ip, m, args, [1,1], state, t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE)
                        elseif ip.cm == 0
                            config.cnf.α = no_control(ip, m, args, [1,1], state, t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE)
                        end
                    end
                end
            end

            # Heat Rate 
            if ip.tm == 1
                heat_rate = heatrate_convective_radiative(S, T_p, m, ρ, vel_pp_mag, config.cnf.α)
            elseif ip.tm == 2
                heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, config.cnf.α)
            end
            
            cp = m.planet.γ / (m.planet.γ - 1) * m.planet.R

            T_r = 0.0
        else
            T_r = 0.0
            heat_rate = 0.0
        end

        α = config.cnf.α
        config.cnf.heat_rate_prev = heat_rate # save current heat rate

        # Convert wind to pp(PCPF) frame
        wE = wind[1] # positive to the east , m / s
        wN = wind[2] # positive to the north , m / s
        wU = wind[3] # positive up , m / s

        wind_pp = wN * uN + wE * uE - wU * uD         # wind velocity in pp frame, m / s 
        vel_pp_rw = vel_pp + wind_pp                  # relative wind vector, m / s
        vel_pp_rw_hat = vel_pp_rw / norm(vel_pp_rw)   # relative wind unit vector 

        # Dynamic Pressure, CHANGE THE VELOCITY WITH THE WIND VELOCITY
        q = 0.5 * ρ * norm(vel_pp_rw)^2               # dynamic pressure based on wind, Pa
               
        # Nominal gravity calculation
        if ip.gm == 0
            gravity_ii = mass * gravity_const(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        elseif ip.gm == 1
            gravity_ii = mass * gravity_invsquared(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        elseif ip.gm == 2
            gravity_ii = mass * (args[:gravity_harmonics] == 1 ? gravity_invsquared(pos_ii_mag, pos_ii, m.planet, mass, vel_ii) : gravity_invsquared_J2(pos_ii_mag, pos_ii, m.planet, mass, vel_ii))
        elseif ip.gm == 3
            gravity_ii = mass * gravity_GRAM(pos_ii, lat, lon, alt, m.planet, mass, vel_ii, el_time, gram_atmosphere, args, gram)
        end

        if length(args[:n_bodies]) != 0
            for k = 1:length(args[:n_bodies])  
                gravity_ii += mass * gravity_n_bodies(config.cnf.et, pos_ii, m.planet, config.cnf.n_bodies_list[k])
            end
        end

        srp_ii = zeros(3) # solar radiation pressure vector
        if args[:srp] == true
            p_srp_unscaled = 4.56e-6  # N / m ^ 2, solar radiation pressure at 1 AU
            srp_ii = mass * srp(m.planet, p_srp_unscaled, m.aerodynamics.reflection_coefficient, m.body.area_tot, m.body.mass, pos_ii, config.cnf.et)
        end

        if args[:gravity_harmonics] == 1
            gravity_ii += mass * m.planet.L_PI' * acc_gravity_pines!(pos_pp, m.planet.Clm, m.planet.Slm, args[:L], args[:M], m.planet.μ, m.planet.Rp_e, m.planet)
        end

        bank_angle = deg2rad(0.0)

        lift_pp_hat = cross(h_pp_hat, vel_pp_rw_hat)

        # Vehicle Aerodynamic Forces
        # CL and CD
        if ip.am == 0
            CL, CD = aerodynamic_coefficient_constant(α, m.body, T_p, S, m.aerodynamics, MonteCarlo)
        elseif ip.am == 1
            CL, CD = aerodynamic_coefficient_fM(α, m.body, T_p, S, m.aerodynamics, MonteCarlo)
        elseif ip.am == 2
            CL, CD = aerodynamic_coefficient_no_ballistic_flight(α, m.body, args, T_p, S, m.aerodynamics, MonteCarlo)
        end

        β = mass / (CD*area_tot)    # ballistic coefficient, kg / m ^ 2

        # Force Calculation
        drag_pp_hat = -vel_pp_rw_hat    # PLanet relative drag force direction

        drag_pp = q * CD * area_tot * drag_pp_hat                       # PLanet relative drag force vector
        lift_pp = q * CL * area_tot * lift_pp_hat * cos(bank_angle)     # PLanet relative lift force vector

        drag_ii = m.planet.L_PI' * drag_pp   # Inertial drag force vector
        lift_ii = m.planet.L_PI' * lift_pp   # Inertial lift force vector

        # Check if propellant mass is greater than 0 kg
        if config.cnf.index_propellant_mass == 1
            if mass - args[:dry_mass] <= 0.5
                config.cnf.index_propellant_mass = 0
                m.engines.T = 0

                if Bool(args[:print_res])
                    println("WARNING: No fuel left!")
                end
            end
        end

        # Thrust
        Δv = m.engines.g_e * m.engines.Isp * log(initial_state.m / mass)

        if ip.tc == 0
            thrust_pp_mag = no_maneuver(t0, m.engines.T, Δv, args, index_phase_aerobraking)
        elseif ip.tc == 1
            thrust_pp_mag = abms(t0, m.engines.T, Δv, args, index_phase_aerobraking)
        elseif ip.tc == 2
            thrust_pp_mag = deceleration_drag_passage(t0, m.engines.T, Δv, args, index_phase_aerobraking)
        end

        # Rodrigues rotation formula to rotate thrust vector of angle phi around angular vector from D direction
        D_L_per_pp_hat = cross(drag_pp_hat, lift_pp_hat)
        thrust_pp_hat =  drag_pp_hat * cos(args[:phi]) + cross(D_L_per_pp_hat, drag_pp_hat) * sin(args[:phi]) + D_L_per_pp_hat * dot(D_L_per_pp_hat, drag_pp_hat) * (1 - cos(args[:phi]))
        #these two ways give the same direction
        thrust_pp = thrust_pp_mag * thrust_pp_hat
        thrust_ii = m.planet.L_PI' * thrust_pp

        # Total Force
        # Total inetrial external force vector on body [N]
        force_ii = drag_ii + lift_ii + gravity_ii + thrust_ii + srp_ii
        
        y_dot[1:3] = vel_ii * (config.cnf.TU / config.cnf.DU)
        y_dot[4:6] = force_ii / mass * (config.cnf.TU^2 / config.cnf.DU) 
        y_dot[7] = -norm(thrust_ii) / (m.engines.g_e * m.engines.Isp) * config.cnf.TU / config.cnf.MU       # mass variation
        y_dot[8] = heat_rate * config.cnf.TU^3 / config.cnf.MU # * 1e-4
        energy = (vel_ii_mag^2)/2 - (m.planet.μ / pos_ii_mag)

        ## SAVE RESULTS
        if Bool(config.cnf.results_save)
            sol = [t0, timereal.year, timereal.month, timereal.day, timereal.hour, timereal.minute,
                        timereal.second, numberofpassage, pos_ii..., vel_ii..., pos_ii_mag, vel_ii_mag, pos_pp..., 
                        pos_pp_mag, vel_pp..., vel_pp_mag, OE[1], OE[2], OE[3], OE[4], OE[5], OE[6],
                        lat, lon, alt, γ_ii, γ_pp, h_ii..., h_pp..., h_ii_mag, h_pp_mag, uD..., uE..., uN..., vN, vE,
                        azi_pp, ρ, T_p, p, wind..., CL, CD, α, S, mass, heat_rate, heat_load, T_r, 
                        q, gravity_ii..., drag_pp..., drag_ii..., lift_pp..., lift_ii..., force_ii..., energy, config.cnf.index_MonteCarlo, Int64(config.cnf.drag_state)]

            sol = reshape(sol, (1, 91))

            push!(config.cnf.solution_intermediate, sol)
        end

        return y_dot
    end

    ## EVENTS: 
    function eventfirststep_condition(y, t, integrator)
        """
        Event function to detect the entry interface downcrossing.
        """
        m = integrator.p[1]
        norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - (args[:EI])*1e3   #  downcrossing
    end
    function eventfirststep_affect!(integrator)
        """
        Event function to terminate the integration at the entry interface downcrossing.
        """
        config.cnf.count_eventfirststep += 1
        terminate!(integrator)
    end
    eventfirststep = ContinuousCallback(eventfirststep_condition, nothing, eventfirststep_affect!)

    function eventfirststep_periapsis_condition(y, t, integrator)
        """
        Event function to detect the periapsis crossing.
        """
        m = integrator.p[1]
        pos_ii = [y[1], y[2], y[3]] * config.cnf.DU  # Inertial position
        vel_ii = [y[4], y[5], y[6]] * config.cnf.DU / config.cnf.TU  # Inertial Velocity

        vi = rvtoorbitalelement(pos_ii, vel_ii, y[7] * config.cnf.MU, m.planet)[6]

        rad2deg(vi) - 180  # downcrossing
    end

    function eventfirststep_periapsis_affect!(integrator)
        """
        Event function to terminate the integration at the periapsis crossing.
        """
        config.cnf.eventfirststep_periapsis += 1
        terminate!(integrator)
    end
    eventfirststep_periapsis = ContinuousCallback(eventfirststep_periapsis_condition, eventfirststep_periapsis_affect!)

    function eventsecondstep_condition(y, t, integrator)
        """
        Event function to detect the atmospheric exit upcrossing.
        """
        m = integrator.p[1]
        norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - (args[:AE])*1e3   # upcrossing
    end
    function eventsecondstep_affect!(integrator)
        """
        Event function to terminate the integration at the atmospheric exit upcrossing.
        """
        config.cnf.count_eventsecondstep += 1
        terminate!(integrator)
    end
    eventsecondstep = ContinuousCallback(eventsecondstep_condition, eventsecondstep_affect!, nothing)

    function reached_EI_condition(y, t, integrator)
        """
        Event function to detect the entry interface downcrossing.
        """
        m = integrator.p[1]
        args = integrator.p[8]
        norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - args[:EI]*1e3  # downcrossing
        # norm(y[1:3]) - m.planet.Rp_e - args[:EI]*1e3  # downcrossing
    end
    function reached_EI_affect!(integrator)
        """
        Event function that doesn't terminate the integration at EI.
        """
        nothing
    end
    reached_EI = ContinuousCallback(reached_EI_condition, nothing, reached_EI_affect!)

    function reached_AE_condition(y, t, integrator)
        """
        Event function to detect the atmospheric exit upcrossing.
        """
        m = integrator.p[1]
        args = integrator.p[8]
        norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - args[:AE]*1e3  # upcrossing
    end
    function reached_AE_affect!(integrator)
        """
        Event function that doesn't terminate the integration at AE.
        """
        nothing
    end
    reached_AE = ContinuousCallback(reached_AE_condition, reached_AE_affect!, nothing)

    function out_drag_passage_condition(y, t, integrator)
        m = integrator.p[1]
        args = integrator.p[8]

        if abs(norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - args[:AE]*1e3) <= 1e-5  # abs(norm(y[1:3]) - m.planet.Rp_e - args[:AE]*1e3) <= 1e-5
            if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 2
                config.cnf.α = m.aerodynamics.α
            elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 3
                config.cnf.α = 0.0
            end
        end

        norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - args[:AE]*1e3  # upcrossing
    end
    function out_drag_passage_affect!(integrator)
        config.cnf.count_out_drag_passage += 1
        terminate!(integrator)
    end
    out_drag_passage = ContinuousCallback(out_drag_passage_condition, out_drag_passage_affect!, nothing)

    function in_drag_passage_condition(y, t, integrator)
        m = integrator.p[1]
        args = integrator.p[8]
        t_prev = integrator.p[5]
        ip = integrator.p[3]
        date_initial = integrator.p[6]
        pos_ii = SVector{3, Float64}([y[1], y[2], y[3]] * config.cnf.DU)  # Inertial position
        vel_ii = SVector{3, Float64}([y[4], y[5], y[6]] * config.cnf.DU / config.cnf.TU)  # Inertial velocity
        pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, config.cnf.et)

        LatLong = rtolatlong(pos_pp, m.planet)

        h0 = LatLong[1]

        if ip.gm == 2
            cond = h0 - args[:EI] * 1e3
            thr = 500
        else
            cond = norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - args[:EI]*1e3
            thr = 1e-5
        end

        if abs(cond) <= thr
            config.controller.guidance_t_eval = collect(t*config.cnf.TU:1/args[:flash1_rate]:(t*config.cnf.TU)+1500)
        
            # State definition for control 2, 3 State used by closed-form solution
            pos_ii = [y[1], y[2], y[3]] * config.cnf.DU                     # Inertial position
            vel_ii = [y[4], y[5], y[6]] * config.cnf.DU / config.cnf.TU     # Inertial velocity
            OE_closedform = rvtoorbitalelement(pos_ii, vel_ii, y[7] * config.cnf.MU, m.planet)
            config.cnf.initial_position_closed_form = OE_closedform
        end

        cond  # downcrossing
    end
    function in_drag_passage_affect!(integrator) 
        config.cnf.count_in_drag_passage += 1
        terminate!(integrator)
    end
    in_drag_passage = ContinuousCallback(in_drag_passage_condition, nothing, in_drag_passage_affect!)

    function in_drag_passage_nt_condition(y, t, integrator)
        m = integrator.p[1]
        args = integrator.p[8]
        t_prev = integrator.p[5]
        ip = integrator.p[3]
        date_initial = integrator.p[6]
        pos_ii = SVector{3, Float64}([y[1], y[2], y[3]] * config.cnf.DU)                    # Inertial position
        vel_ii = SVector{3, Float64}([y[4], y[5], y[6]] * config.cnf.DU / config.cnf.TU)    # Inertial velocity
        pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, config.cnf.et)

        LatLong = rtolatlong(pos_pp, m.planet)#, args[:topography_model] == "Spherical Harmonics" && norm(pos_ii) < m.planet.Rp_e + args[:EI]*1e3)

        h0 = LatLong[1]

        if ip.gm == 2 && !Bool(args[:drag_passage])
            cond = (h0 - args[:EI]*1e3)
            thr = 500
        else
            cond = norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - args[:EI]*1e3
            thr = 1e-3
        end

        if abs(cond) <= thr && length(config.cnf.initial_position_closed_form) == 0
            config.controller.guidance_t_eval = collect(t*config.cnf.TU:1/args[:flash1_rate]:(t*config.cnf.TU)+1500)

            # State definition for control 2, 3 State used by closed-form solution
            OE_closedform = rvtoorbitalelement(pos_ii, vel_ii, y[7] * config.cnf.MU, m.planet)
            config.cnf.initial_position_closed_form = OE_closedform
        end

        norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - args[:EI]*1e3  # downcrossing
    end
    function in_drag_passage_nt_affect!(integrator)
        config.cnf.count_in_drag_passage_nt += 1
        nothing
    end
    in_drag_passage_nt = ContinuousCallback(in_drag_passage_nt_condition, nothing, in_drag_passage_nt_affect!)

    function apoapsispoint_condition(y, t, integrator)
        m = integrator.p[1]
        pos_ii = SVector{3, Float64}([y[1], y[2], y[3]]) * config.cnf.DU  # Inertial position
        vel_ii = SVector{3, Float64}([y[4], y[5], y[6]]) * config.cnf.DU / config.cnf.TU  # Inertial Velocity

        vi = rvtoorbitalelement(pos_ii, vel_ii, y[7] * config.cnf.MU, m.planet)[6]
        
        rad2deg(vi) - 180 # upcrossing
    end
    function apoapsispoint_affect!(integrator)
        config.cnf.count_apoapsispoint += 1
        terminate!(integrator)
    end
    apoapsispoint = ContinuousCallback(apoapsispoint_condition, apoapsispoint_affect!, nothing)

    function periapsispoint_condition(y, t, integrator)
        m = integrator.p[1]
        pos_ii = SVector{3, Float64}([y[1], y[2], y[3]]) * config.cnf.DU  # Inertial position
        vel_ii = SVector{3, Float64}([y[4], y[5], y[6]]) * config.cnf.DU / config.cnf.TU  # Inertial Velocity

        vi = rvtoorbitalelement(pos_ii, vel_ii, y[7] * config.cnf.MU, m.planet)[6]

        rad2deg(vi) - 180  # downcrossing
    end
    function periapsispoint_affect!(integrator)
        config.cnf.count_periapsispoint += 1
        r_p, _ = r_intor_p!(SVector{3, Float64}(integrator.u[1:3] * config.cnf.DU), 
                        SVector{3, Float64}(integrator.u[4:6] * config.cnf.DU / config.cnf.TU), 
                        integrator.p[1].planet, 
                        config.cnf.et)
        r, lat, lon = rtolatlong(r_p, integrator.p[1].planet, args[:topography_model] == "Spherical Harmonics")
        append!(config.cnf.altitude_periapsis, r*1e-3)
        append!(config.cnf.latitude_periapsis, rad2deg(lat))
        append!(config.cnf.longitude_periapsis, rad2deg(lon))
        nothing
    end
    periapsispoint = ContinuousCallback(periapsispoint_condition, nothing, periapsispoint_affect!)

    function impact_condition(y, t, integrator)
        m = integrator.p[1]
        args = integrator.p[8]

        if args[:body_shape] == "Blunted Body" || args[:type_of_mission] == "Entry"
            min_alt = 0 * 1e3
        else
            min_alt = 35 * 1e3
        end

        norm(y[1:3]) * config.cnf.DU - (m.planet.Rp_e + min_alt) # upcrossing and downcrossing
    end
    function impact_affect!(integrator)
        config.cnf.count_impact += 1
        terminate!(integrator)
    end
    impact = ContinuousCallback(impact_condition, impact_affect!)

    function apoapsisgreaterperiapsis_condition(y, t, integrator)
        m = integrator.p[1]
        args = integrator.p[8]

        r = y[1:3] * config.cnf.DU
        v = y[4:6] * config.cnf.DU / config.cnf.TU
        Energy = norm(v)^2 * 0.5 - m.planet.μ / norm(r)
        a = -m.planet.μ / (2 * Energy)
        h = cross(r, v)
        e = sqrt(1 + (2 * Energy * dot(h,h) / (m.planet.μ)^2))

        r_a = a * (1 + e)
        r_p = a * (1 - e)

        if r_a < r_p && args[:body_shape] == "Spacecraft"
            println("Periapsis greater than apoapsis!")
        end

        r_a - r_p  # upcrossing and downcrossing
    end
    function apoapsisgreaterperiapsis_affect!(integrator)
        config.cnf.count_apoapsisgreaterperiapsis += 1
        terminate!(integrator)
    end
    apoapsisgreaterperiapsis = ContinuousCallback(apoapsisgreaterperiapsis_condition, apoapsisgreaterperiapsis_affect!)

    function stop_firing_condition(y, t, integrator)
        m = integrator.p[1]
        args = integrator.p[8]
        initial_state = integrator.p[9]

        mass = y[7] * config.cnf.MU
        Δv = (m.engines.g_e * m.engines.Isp) * log(initial_state.m/mass)
        
        Δv - args[:delta_v]  # upcrossing and downcrossing
    end
    function stop_firing_affect!(integrator)
        config.cnf.count_stop_firing += 1
        terminate!(integrator)
    end
    stop_firing = ContinuousCallback(stop_firing_condition, stop_firing_affect!)

    function guidance_condition(y, t, integrator)
        if config.controller.stored_state == 1
            config.controller.count_prev_controller = config.controller.count_controller
        end

        if t * config.cnf.TU - config.controller.t > 1 # t - config.controller.t > 1
            println("Decrease step size of integration and tolerance")
        end

        if abs(t * config.cnf.TU - config.controller.t) <= 1e-8 # abs(t - config.controller.t) <= 1e-8 
            config.controller.stored_state = 0
            config.controller.count_controller += 1
            config.controller.t = config.controller.guidance_t_eval[config.controller.count_controller]
        end

        t * config.cnf.TU - config.controller.t  # upcrossing
    end
    function guidance_affect!(integrator)
        nothing
    end
    guidance = ContinuousCallback(guidance_condition, guidance_affect!, nothing)

    function heat_rate_check_condition(y, t, integrator)
        m = integrator.p[1]
        args = integrator.p[8]

        if args[:control_mode] == 1
            x = 120
        else
            x = 160
        end

        if abs(norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - x*1e3) <= 1e-5 # abs(norm(y[1:3]) - m.planet.Rp_e - x*1e3) <= 1e-5
            config.controller.guidance_t_eval = collect(range(start=t * config.cnf.TU, stop=(t * config.cnf.TU)+2500, step=1/args[:flash1_rate]))
        end

        norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - x*1e3  # upcrossing and downcrossing
    end
    function heat_rate_check_affect!(integrator)
        config.cnf.count_heat_rate_check += 1
        terminate!(integrator)
    end
    heat_rate_check = ContinuousCallback(heat_rate_check_condition, heat_rate_check_affect!)

    function heat_load_check_exit_condition(y, t, integrator)
        m = integrator.p[1]
        args = integrator.p[8]

        if args[:control_mode] == 1
            x = 120
        else
            x = 160
        end

        norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - x*1e3  # upcrossing
    end
    function heat_load_check_exit_affect!(integrator)
        config.cnf.count_heat_load_check_exit += 1
        terminate!(integrator)
    end
    heat_load_check_exit = ContinuousCallback(heat_load_check_exit_condition, heat_load_check_exit_affect!, nothing)

    function final_entry_altitude_reached_condition(y, t, integrator)
        m = integrator.p[1]
        args = integrator.p[8]
        min_alt = args[:final_altitude]

        r_p, _ = r_intor_p!(SVector{3, Float64}(y[1:3] * config.cnf.DU), 
                        SVector{3, Float64}(y[4:6] * config.cnf.DU / config.cnf.TU), 
                        integrator.p[1].planet, 
                        config.cnf.et)
        alt, lat, lon = rtolatlong(r_p, integrator.p[1].planet, args[:topography_model] == "Spherical Harmonics")
        alt - min_alt # upcrossing and downcrossing
    end

    function final_entry_altitude_reached_affect!(integrator)
        config.cnf.count_final_entry_altitude_reached += 1
        terminate!(integrator)
    end
    final_entry_altitude_reached = ContinuousCallback(final_entry_altitude_reached_condition, final_entry_altitude_reached_affect!)

    time_0 = 0
    if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 2
        config.cnf.α = m.aerodynamics.α
    elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 3
        config.cnf.α = 0.0
    end

    stop_simulation = false

    # Double check these 2
    save_pre_index = 1
    save_post_index = 1

    config.cnf.impact = false
    config.cnf.solution_intermediate = []
    config.cnf.count_dori = 0
    config.cnf.atmospheric_data = Dict()
    config.cnf.previous_atmospheric_data = Dict()
    config.cnf.ascending_phase = false
    config.cnf.evaluate_switch_heat_load = false
    config.cnf.state_inner_boundary_atmosphere = [] # used in Density model for vi def
    config.cnf.time_IP = config.cnf.time_OP
    config.cnf.time_IEI = 0
    config.cnf.time_OEI = 0
    config.cnf.time_switch_1 = 0.0
    config.cnf.time_switch_2 = 0.0

    if args[:heat_load_sol] == 2
        config.cnf.time_switch_2 = 1000.0
    end

    config.cnf.timer_revaluation = 0
    config.cnf.closed_form_solution_off = 1         # used in closed form solution online to run the solution only once
    config.cnf.initial_position_closed_form = []
    config.cnf.heat_rate_list = []                  # checked this - if used
    config.cnf.α_list = []                          # checked this - if used

    config.controller.guidance_t_eval = []
    config.controller.count_prev_controller = 0
    config.controller.count_controller = 1
    config.controller.stored_state = 1
    config.controller.prev_time = 0
    config.controller.t = 0

    config.cnf.security_mode = false
    config.cnf.stop_simulation = false
    config.cnf.results_save = 1
    config.cnf.drag_state = false
    config.cnf.α_past = m.aerodynamics.α

    if norm(r0) - m.planet.Rp_e <= args[:EI]*1e3
        config.cnf.drag_state = true
        config.cnf.initial_position_closed_form = OE
    end

    config.cnf.sensible_loads = false
    config.cnf.counter_integrator = 0

    continue_campaign = false

    # Def initial conditions
    in_cond = [r0[1], r0[2], r0[3], v0[1], v0[2], v0[3], Mass+1e-10, 0.0]

    # non dimensionalization
    in_cond[1:3] = in_cond[1:3] / config.cnf.DU
    in_cond[4:6] = in_cond[4:6] * config.cnf.TU / config.cnf.DU
    in_cond[7] = in_cond[7] / config.cnf.MU
    in_cond[8] = in_cond[8] * config.cnf.TU^2 / config.cnf.MU # * 1e4

    # If aerobraking maneuver allowed, add a prephase 0
    range_phase_i = 1
    if args[:thrust_control] == "Aerobraking Maneuver" && args[:keplerian] == false
        range_phase_i = 0
    end

    index_phase_aerobraking = range_phase_i # for scope reasons
    aerobraking_phase = range_phase_i       # for scope reasons
    method = Tsit5()                        # for scope reasons

    # Solve Equations of Motion 
    for aerobraking_phase in range(range_phase_i, 3)
        index_phase_aerobraking = aerobraking_phase

        if (index_steps_EOM == 1 || Bool(args[:drag_passage])) && (aerobraking_phase == 1 || aerobraking_phase == 3 || aerobraking_phase == 0)
            continue
        elseif cmp(lowercase(args[:type_of_mission]), "orbits") == 0 && args[:keplerian] == true && aerobraking_phase == 2  
            continue
        end

        # Definition of eventsecondstep
        if aerobraking_phase == 0
            events = CallbackSet(stop_firing, apoapsisgreaterperiapsis, impact)
            t_event_0 = "stop_firing"
            t_event_1 = "apoapsisgreaterperiapsis"
        elseif aerobraking_phase == 1 && args[:keplerian] == true
            events = CallbackSet(eventfirststep_periapsis, apoapsisgreaterperiapsis, impact, periapsispoint)
            t_event_0 = "eventfirststep_periapsis"
            t_event_1 = "apoapsisgreaterperiapsis"
        elseif aerobraking_phase == 1
            events = CallbackSet(eventfirststep, apoapsisgreaterperiapsis, impact)
            t_event_0 = "eventfirststep"
            t_event_1 = "apoapsisgreaterperiapsis"
        elseif aerobraking_phase == 2 && Bool(args[:drag_passage]) && args[:type_of_mission] != "Entry"
            events = CallbackSet(out_drag_passage, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, impact)
            t_event_0 = "out_drag_passage"
            t_event_1 = "periapsispoint"
        elseif aerobraking_phase == 2 && Bool(args[:drag_passage]) && args[:type_of_mission] == "Entry"
            events = CallbackSet(out_drag_passage, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, final_entry_altitude_reached)
            t_event_0 = "final_altitude_reached"
            t_event_1 = "out_drag_passage"
        elseif aerobraking_phase == 2 && index_steps_EOM == 1 && args[:body_shape] == "Blunted Cone"
            events = CallbackSet(out_drag_passage, apoapsispoint, periapsispoint, impact)
            t_event_0 = "out_drag_passage"
            t_event_1 = "apoapsispoint"
        elseif aerobraking_phase == 2 && index_steps_EOM == 1 && args[:drag_passage] == false
            events = CallbackSet(apoapsispoint, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, impact)
            t_event_0 = "apoapsispoint"
            t_event_1 = "periapsispoint"
        elseif aerobraking_phase == 2
            events = CallbackSet(eventsecondstep, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, impact)
            t_event_0 = "eventsecondstep"
            t_event_1 = "periapsispoint"
        elseif aerobraking_phase == 3
            events = CallbackSet(apoapsispoint, periapsispoint, apoapsisgreaterperiapsis, impact)
            t_event_0 = "apoapsispoint"
            t_event_1 = "periapsispoint"
        end

        if index_phase_aerobraking == 2 || (index_phase_aerobraking == 1 && args[:keplerian] == 1)
            save_pre_index = length(config.solution.orientation.time) + 1
            simulator = args[:integrator]
        end

        # Definition step size
        if aerobraking_phase == 1 || aerobraking_phase == 3
            step = 5.0
            r_tol = 1e-10
            a_tol = 1e-12
            simulator = "Julia"
            method = Tsit5()
            save_ratio = 2
        elseif aerobraking_phase == 0
            step = 1.0
            r_tol = 1e-9
            a_tol = 1e-11
            simulator = "Julia"
            method = Tsit5()
            save_ratio = 5
        elseif aerobraking_phase == 2
            if args[:integrator] == "Julia"
                step = 0.1
                r_tol = 1e-9
                a_tol = 1e-11

                if MonteCarlo
                    method = Tsit5()
                else
                    method = Tsit5()
                end

                save_ratio = args[:save_rate]
            else
                step = 1/args[:trajectory_rate]
                save_ratio = round(args[:trajectory_rate] / args[:save_rate])
            end
        end

        # Definition length simulation
        length_sim = 1e8
        i_sim = 0
        time_solution = []

        config.cnf.continue_simulation = true

        while config.cnf.continue_simulation
            index_phase_aerobraking = aerobraking_phase
            # if control mode =! 0, redefine sim setting and creates two more phases until reaching EI and out of the AE phase 2: between 120 km alt
            if aerobraking_phase == 2 && (args[:control_mode] != 0 && args[:control_in_loop] == 0 && config.cnf.drag_state == true && config.cnf.sensible_loads == true && config.cnf.ascending_phase == false)
                simulator = args[:integrator]
                events = CallbackSet(out_drag_passage, heat_load_check_exit, periapsispoint, apoapsisgreaterperiapsis, impact)
                t_event_0 = "out_drag_passage"
                t_event_1 = "heat_load_check_exit"

                if args[:integrator] == "Costumed"
                    step = 1/args[:trajectory_rate]
                else
                    method = Tsit5()
                    if Bool(args[:control_in_loop])
                        step = 10/args[:trajectory_rate]
                    else
                        step = 1/(4*args[:flash1_rate])
                        r_tol = 1e-6
                        a_tol = 1e-7

                        if step >= args[:flash1_rate]
                            step = 1/(2*args[:flash1_rate])
                        end
                        
                        events = CallbackSet(out_drag_passage, heat_load_check_exit, periapsispoint, guidance, apoapsisgreaterperiapsis, impact)
                    end
                end

                config.controller.t = config.controller.guidance_t_eval[config.controller.count_controller]
            
            #phase 1.75: between EI km alt and 120 kmconfig.controller.t
            elseif aerobraking_phase == 2 && args[:control_mode] != 0 && args[:control_in_loop] == 0 && config.cnf.drag_state == true && config.cnf.ascending_phase == false
                events = CallbackSet(periapsispoint, out_drag_passage, heat_rate_check, apoapsisgreaterperiapsis, impact)
                t_event_0 = "periapsispoint"
                t_event_1 = "out_drag_passage"
                simulator = "Julia"
                index_phase_aerobraking = 1.75
                step = 0.05
                r_tol = 1e-6
                a_tol = 1e-7
                method = Tsit5() # KenCarp58(autodiff = false) # Tsit5()
            # phase 1.5: between 250 km alt and EI km
            elseif aerobraking_phase == 2 && args[:control_mode] != 0 && args[:control_in_loop] == 0 && config.cnf.drag_state == false && config.cnf.ascending_phase == false
                events = CallbackSet(periapsispoint, in_drag_passage, apoapsisgreaterperiapsis, impact)
                t_event_0 = "periapsispoint"
                t_event_1 = "in_drag_passage"
                simulator = "Julia"
                step = 0.1
                r_tol = 1e-9
                a_tol = 1e-10
                index_phase_aerobraking = 1.5
                method = Tsit5()
            # phase 2.25: between 120 km alt and AE km
            elseif aerobraking_phase == 2 && args[:control_mode] != 0 && args[:control_in_loop] == 0 && config.cnf.drag_state == true && config.cnf.ascending_phase == true
                events = CallbackSet(periapsispoint, out_drag_passage, apoapsisgreaterperiapsis, impact)
                t_event_0 = "periapsispoint"
                t_event_1 = "out_drag_passage"
                simulator = "Julia"
                step = 0.1
                r_tol = 1e-9
                a_tol = 1e-10
                index_phase_aerobraking = 2.25
                method = Tsit5()
            # phase 2.5: between AE km alt and 250 km
            elseif aerobraking_phase == 2 && args[:control_mode] != 0 && args[:control_in_loop] == 0 && config.cnf.ascending_phase == true
                events = CallbackSet(eventsecondstep, periapsispoint, eventsecondstep, apoapsisgreaterperiapsis, impact)
                t_event_0 = "eventsecondstep"
                t_event_1 = "periapsispoint"
                simulator = "Julia"
                index_phase_aerobraking = 2.5
                step = 0.5
                r_tol = 1e-8
                a_tol = 1e-9
                method = Tsit5()
            end

            if Bool(args[:print_res])
                println("Step #", index_phase_aerobraking)
            end

            if simulator == "Julia"
                # counter for events
                config.cnf.count_eventfirststep = 0
                config.cnf.eventfirststep_periapsis = 0
                config.cnf.count_eventsecondstep = 0
                config.cnf.count_reached_EI = 0
                config.cnf.count_reached_AE = 0
                config.cnf.count_out_drag_passage = 0
                config.cnf.count_in_drag_passage = 0
                config.cnf.count_in_drag_passage_nt = 0
                config.cnf.count_apoapsispoint = 0
                config.cnf.count_periapsispoint = 0
                config.cnf.count_impact = 0
                config.cnf.count_apoapsisgreaterperiapsis = 0
                config.cnf.count_stop_firing = 0
                config.cnf.count_guidance = 0
                config.cnf.count_heat_rate_check = 0
                config.cnf.count_heat_load_check_exit = 0
                config.cnf.count_final_entry_altitude_reached = 0

                ## Julia Integrator
                # Time initialization
                initial_time, final_time = time_0 / config.cnf.TU, (time_0 + length_sim) / config.cnf.TU 

                # Parameter Definition
                param = (m, index_phase_aerobraking, ip, aerobraking_phase, t_prev, date_initial, time_0, args, initial_state, gram_atmosphere, gram)

                # Run simulation
                prob = ODEProblem(f!, in_cond, (initial_time, final_time), param)
                sol = solve(prob, method, abstol=a_tol, reltol=r_tol, callback=events)
                config.cnf.counter_integrator += 1
                in_cond = [sol[1,end], sol[2,end], sol[3, end], sol[4, end], sol[5, end], sol[6, end], sol[7, end], sol[8, end]]
              

                # Save results 
                push!(time_solution, (sol.t * config.cnf.TU)...)
                time_0 = time_solution[end]

                if aerobraking_phase == 0
                    new_periapsis(m, in_cond[1:3] * config.cnf.DU, in_cond[4:6] * config.cnf.DU / config.cnf.TU, args)
                end
            elseif simulator == "Costumed"
                # if args[:integrator] == "Costumed"
                #     # mutable struct Sol
                #     #     t_events::Vector{Vector{Float64}}
                #     # end

                #     # sol = Sol([[],[]])
                # end

                # while stop_simulation == false
                #     initial_time = time_0
                #     ## Costumed Integrator
                #     config.cnf.MarsGram_recall = 1
                #     config.cnf.results_save = 0
                #     y, t, stop_simulation, sol = RK4(f, step, initial_time, in_cond, m, T_ijk, index_phase_aerobraking, args, sol)

                #     in_cond = [y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8]]
                #     config.cnf.counter_integrator += 1
                #     config.cnf.results_save = 1
                #     f(t, in_cond, m, index_phase_aerobraking, ip)

                #     # New initial Condition
                #     time_0 = t_prev

                #     # Save Results
                #     append!(time_solution, t)

                #     ## Guidance, Navigation and Control
                #     if args[:control_mode] != 0
                #         # DO LATER
                #     end
                # end

                # i_sim += 1
            end

            # Define breaker campaign impact km or apoapsis greater than periapsis
            continue_campaign = event(config.cnf.count_impact, config.cnf.count_apoapsisgreaterperiapsis)

            if continue_campaign == false
                config.cnf.impact = true
                break
            end

            time_ev_0, time_ev_1 = 0, 0

            if t_event_0 == "stop_firing"
                time_ev_0 = config.cnf.count_stop_firing
            elseif t_event_0 == "eventfirststep"
                time_ev_0 = config.cnf.count_eventfirststep
            elseif t_event_0 == "eventfirststep_periapsis"
                time_ev_0 = config.cnf.eventfirststep_periapsis
            elseif t_event_0 == "eventsecondstep"
                time_ev_0 = config.cnf.count_eventsecondstep
            elseif t_event_0 == "out_drag_passage"
                time_ev_0 = config.cnf.count_out_drag_passage
            elseif t_event_0 == "in_drag_passage"
                time_ev_0 = config.cnf.count_in_drag_passage
            elseif t_event_0 == "in_drag_passage_nt"
                time_ev_0 = config.cnf.count_in_drag_passage_nt
            elseif t_event_0 == "apoapsispoint"
                time_ev_0 = config.cnf.count_apoapsispoint
            elseif t_event_0 == "periapsispoint"
                time_ev_0 = config.cnf.count_periapsispoint
            elseif t_event_0 == "apoapsisgreaterperiapsis"
                time_ev_0 = config.cnf.count_apoapsisgreaterperiapsis
            elseif t_event_0 == "heat_load_check_exit"
                time_ev_0 = config.cnf.count_heat_load_check_exit
            elseif t_event_0 == "final_altitude_reached"
                time_ev_0 = config.cnf.count_final_entry_altitude_reached
            end

            if t_event_1 == "stop_firing"
                time_ev_1 = config.cnf.count_stop_firing
            elseif t_event_1 == "eventfirststep"
                time_ev_1 = config.cnf.count_eventfirststep
            elseif t_event_1 == "eventfirststep_periapsis"
                time_ev_1 = config.cnf.eventfirststep_periapsis
            elseif t_event_1 == "eventsecondstep"
                time_ev_1 = config.cnf.count_eventsecondstep
            elseif t_event_1 == "out_drag_passage"
                time_ev_1 = config.cnf.count_out_drag_passage
            elseif t_event_1 == "in_drag_passage"
                time_ev_1 = config.cnf.count_in_drag_passage
            elseif t_event_1 == "in_drag_passage_nt"
                time_ev_1 = config.cnf.count_in_drag_passage_nt
            elseif t_event_1 == "apoapsispoint"
                time_ev_1 = config.cnf.count_apoapsispoint
            elseif t_event_1 == "periapsispoint"
                time_ev_1 = config.cnf.count_periapsispoint
            elseif t_event_1 == "apoapsisgreaterperiapsis"
                time_ev_1 = config.cnf.count_apoapsisgreaterperiapsis
            elseif t_event_1 == "heat_load_check_exit"
                time_ev_1 = config.cnf.count_heat_load_check_exit
            elseif t_event_1 == "final_altitude_reached"
                time_ev_1 = config.cnf.count_final_entry_altitude_reached
            end

            # Breaker conditions
            if simulator == "Julia"
                if time_ev_0 != 0 || (Bool(args[:drag_passage]) && index_phase_aerobraking == 2.25 && time_ev_1 != 0)
                    config.cnf.continue_simulation = false
                    break
                end
            else
                if args[:control_mode] == 0 && stop_simulation == true
                    config.cnf.continue_simulation = false
                    break
                end
            end
        end

        # Save Results
        time_0 = save_results(time_solution, save_ratio)

        if index_phase_aerobraking == 2 || index_phase_aerobraking == 2.5 || (index_phase_aerobraking == 2.25 && Bool(args[:drag_passage])) || (index_phase_aerobraking == 3 && Bool(args[:keplerian]))
            save_post_index = length(config.solution.orientation.time)
        end

        # Re-Set count index to 0
        config.cnf.count_phase = 0

        # Define breaker campaign
        if continue_campaign == false
            if save_post_index == 1
                save_post_index = length(config.solution.orientation.time)
            end

            break
        end
    end

    # Re-set count index to 0
    config.cnf.count_dori = 0

    # Check tolerance here on true anomaly
    if args[:drag_passage] == false && (pi - config.solution.orientation.oe[end][end] > 1e-5) && continue_campaign == true && args[:body_shape] != "Blunted Cone"
        final_conditions_notmet = true
        events = CallbackSet(apoapsispoint)
    else
        final_conditions_notmet = false
    end

    count_temp = 0

    while final_conditions_notmet
        in_cond = [config.solution.orientation.pos_ii[1][end], config.solution.orientation.pos_ii[2][end], config.solution.orientation.pos_ii[3][end], 
                   config.solution.orientation.vel_ii[1][end], config.solution.orientation.vel_ii[2][end], config.solution.orientation.vel_ii[3][end],
                   config.solution.performance.mass[end], config.solution.performance.heat_load[end]]

        # non dimensionalization
        in_cond[1:3] = in_cond[1:3] / config.cnf.DU
        in_cond[4:6] = in_cond[4:6] * config.cnf.TU / config.cnf.DU
        in_cond[7] = in_cond[7] / config.cnf.MU
        in_cond[8] = in_cond[8] * config.cnf.TU^2 / config.cnf.MU # * 1e4

                    
        initial_time, final_time = time_0 / config.cnf.TU, (time_0 + 1000) / config.cnf.TU 
        step = 0.05
        r_tol = 1e-12
        a_tol = 1e-13

        # counter for events
        config.cnf.count_eventfirststep = 0
        config.cnf.eventfirststep_periapsis = 0
        config.cnf.count_eventsecondstep = 0
        config.cnf.count_reached_EI = 0
        config.cnf.count_reached_AE = 0
        config.cnf.count_out_drag_passage = 0
        config.cnf.count_in_drag_passage = 0
        config.cnf.count_in_drag_passage_nt = 0
        config.cnf.count_apoapsispoint = 0
        config.cnf.count_periapsispoint = 0
        config.cnf.count_impact = 0
        config.cnf.count_apoapsisgreaterperiapsis = 0
        config.cnf.count_stop_firing = 0
        config.cnf.count_guidance = 0
        config.cnf.count_heat_rate_check = 0
        config.cnf.count_heat_load_check_exit = 0

        # Parameter Definition
        param = (m, index_phase_aerobraking, ip, aerobraking_phase, t_prev, date_initial, time_0, args, initial_state, gram_atmosphere, gram)

        # Run simulation
        prob = ODEProblem(f!, in_cond, (initial_time, final_time), param)
        sol = solve(prob, method, abstol=a_tol, reltol=r_tol, callback=events)

        config.cnf.counter_integrator += 1
        time_0 = save_results(sol.t * config.cnf.TU, 1.0)
        count_temp += 1

        if count_temp > 15
            break
        end

        if args[:drag_passage] == false && pi - config.solution.orientation.oe[end][end] > 1e-5 && continue_campaign == true
            final_conditions_notmet = true
            events = CallbackSet(apoapsispoint)
        else
            final_conditions_notmet = false
        end
    end

    config.cnf.save_index_heat = length(config.solution.orientation.time)
    config.cnf.time_OP = length(config.solution.orientation.time)


    append!(config.cnf.max_heatrate, maximum(config.solution.performance.heat_rate[save_pre_index:save_post_index]))
    config.cnf.Δv_man = (m.engines.g_e * m.engines.Isp) * log(m.body.mass / config.solution.performance.mass[end])


    if Bool(args[:print_res])
        # Print Actual periapsis altitude and Vacuum periapsis altitude
        if args[:type_of_mission] != "Entry"
            println("Actual periapsis altitude " * string(config.cnf.altitude_periapsis[end]) * " km - Vacuum periapsis altitude = " * string((config.solution.orientation.oe[1][end] * (1 - config.solution.orientation.oe[2][end]) - m.planet.Rp_e)*1e-3) * " km")

        # Print Ra new (Apoapsis)
            println("Ra new = " * string((config.solution.orientation.oe[1][end] * (1 + config.solution.orientation.oe[2][end]))*1e-3) * " km")
        end

        # Print Heat Rate and Heat Load
        if args[:keplerian] == false
            println("HEAT RATE IS " * string(maximum(config.solution.performance.heat_rate[save_pre_index:save_post_index-1])) * " W/cm^2")
            println("HEAT LOAD IS " * string(config.solution.performance.heat_load[save_post_index-1]) * " J/cm^2")
        end

        # Print Fuel Mass
        println("Fuel Mass is " * string(config.solution.performance.mass[end] - args[:dry_mass]) * " kg")

        # Print Total Time
        if args[:keplerian] == false
            println("Total time is " * string(config.solution.orientation.time[save_post_index-1] - config.solution.orientation.time[save_pre_index]) * " s")
        else
            println("Total time is " * string(config.solution.orientation.time[end] - config.solution.orientation.time[save_pre_index]) * " s")
        end

        # Print Delta-v and Delta-E
        println("Delta-v is " * string(config.cnf.Δv_man) * " m/s")
        println("Delta-E is " * string((config.solution.forces.energy[Int64(config.cnf.time_OP)-1] - config.solution.forces.energy[Int64(config.cnf.time_IP)]) * 1e-3) * " kJ")

        # Find periapsis latitude and longitude
        min_index = argmin(config.solution.orientation.alt[save_pre_index:save_post_index]) + save_pre_index - 1
        println("Latitude of periapsis is " * string(rad2deg(config.solution.orientation.lat[min_index])) * " deg")
        println("Longitude of periapsis is " * string(rad2deg(config.solution.orientation.lon[min_index])) * " deg")

        # If the body shape is 'Blunted Cone', print additional information
        if args[:body_shape] == "Blunted Cone"
            # Max Dynamic Pressure
            max_value = maximum(config.solution.performance.q)
            max_index = argmax(config.solution.performance.q)
            println("Max Dynamic Pressure " * string(max_value) * " N/m^2 at time " * string(config.solution.orientation.time[max_index]) * " s")

            # Max Heat Rate
            max_value = maximum(config.solution.performance.heat_rate)
            max_index = argmax(config.solution.performance.heat_rate)
            println("Max Heat Rate " * string(max_value) * " W/cm^2 at time " * string(config.solution.orientation.time[max_index]) * " s")

            # Max Heat Load
            max_value = maximum(config.solution.performance.heat_load)
            max_index = argmax(config.solution.performance.heat_load)
            println("Max Heat Load " * string(max_value) * " J/cm^2 at time " * string(config.solution.orientation.time[max_index]) * " s")
        end
    end

    append!(config.cnf.periapsis_list, minimum(config.solution.orientation.alt[save_pre_index:save_post_index])*1e-3)
    append!(config.cnf.orbit_number_list, config.cnf.count_numberofpassage + 1)
    append!(config.cnf.Δv_list, config.cnf.Δv_man)
    return continue_campaign
end