include("../../physical_models/Density_models.jl")
include("../../physical_models/Aerodynamic_models.jl")
include("../../physical_models/Gravity_models.jl")
include("../../physical_models/Thermal_models.jl")

include("../../utils/Reference_system.jl")
include("../../utils/Closed_form_solution.jl")

# include("Control.jl")
# include("heatload_control/Utils_timeswitch.jl")

using PreallocationTools
# using NonlinearSolve
using NLsolve
using LinearAlgebra
using DifferentialEquations
using Dates
using AstroTime
using SPICE
using PythonCall
sys = pyimport("sys")

# import .config
# import .ref_sys

function asim_ctrl_targeting_plot(ip, m, time_0, OE, args, hf, vf, γf, energy_f, v_E, k_cf, heat_rate_control, gram_atmosphere=nothing)
    sys.path.append(args[:directory_Gram])
    gram = pyimport("gram")

    wind_m = false
    if ip.wm == 1
        wind_m = true
    end

    MonteCarlo = false
    if ip.mc == 1
        MonteCarlo = true
    end

    version = args[:Gram_version]

    r0, v0 = orbitalelemtorv(OE, m.planet)

    r0 = SVector{3, Float64}(r0)
    v0 = SVector{3, Float64}(v0)

    # Clock
    date_initial = from_utc(DateTime(m.initial_condition.year, 
                                    m.initial_condition.month,
                                    m.initial_condition.day, 
                                    m.initial_condition.hour, 
                                    m.initial_condition.minute, 
                                    m.initial_condition.second))

    if config.cnf.count_numberofpassage != 1
        t_prev = config.solution.orientation.time[end]
    else
        t_prev = m.initial_condition.time_rot # value(seconds(date_initial - from_utc(DateTime(2000, 1, 1, 12, 0, 0)))) # m.initialcondition.time_rot
    end

    r0_pp, v0_pp = r_intor_p!(r0, v0, m.planet, config.cnf.et)

    T = m.planet.T    # fixed temperature
    RT = T * m.planet.R

    S = norm(v0_pp) / sqrt(2 * RT)

    # println("v0_pp: ", norm(v0_pp))
    # println("S at initial condition: ", S)
    # println("T at initial condition: ", T)

    CL_90, CD_90 = aerodynamic_coefficient_fM(pi/2, m.body, T, S, m.aerodynamics, 0)
    CL_0, CD_0 = aerodynamic_coefficient_fM(0, m.body, T, S, m.aerodynamics, 0)
    CD_slope = (CD_90 - CD_0) / (pi/2)

    # println("CD_0 outside integrator: ", CD_0)
    # println("CL_0 outside integrator: ", CL_0)
    # println("CD_90 outside integrator: ", CD_90)
    # println("CL_90 outside integrator: ", CL_90)
    # println("CD_slope outside integrator: ", CD_slope)
    # println("k_cf outside integrator: ", k_cf)

    function f_ctrl!(y_dot, in_cond, param, t0)

        # param = get_tmp(param, first(in_cond) * t0)

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
        k_cf = param[12]

        # println("date_initial: ", date_initial)
        # println("t0: ", t0)

        # Clock
        current_epoch = date_initial + t0*seconds # Precompute the current epoch
        time_real = DateTime(current_epoch) # date_initial + Second(t0)
        timereal = ref_sys.clock(Dates.year(time_real), Dates.month(time_real), Dates.day(time_real), Dates.hour(time_real), Dates.minute(time_real), Dates.second(time_real))

        # Timing variables
        el_time = value(seconds(current_epoch - m.initial_condition.DateTimeIC)) # Elapsed time since the beginning of the simulation
        current_time =  value(seconds(current_epoch - m.initial_condition.DateTimeJ2000)) # current time in seconds since J2000
        time_real_utc = to_utc(time_real) # Current time in UTC as a DateTime object
        et = utc2et(time_real_utc) # Current time in Ephemeris Time
        m.planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_"*uppercase(m.planet.name), et))*m.planet.J2000_to_pci' # Construct a rotation matrix from J2000 (Planet-fixed frame 0.0 seconds past the J2000 epoch) to planet-fixed frame
        
        pos_ii = SVector{3, Float64}(in_cond[1:3])       # Inertial position 
        vel_ii = SVector{3, Float64}(in_cond[4:6])       # Inertial velocity
        mass = m.body.mass                             # Mass kg
        pos_ii_mag = norm(pos_ii)   # Magnitude of the inertial position
        vel_ii_mag = norm(vel_ii)   # Magnitude of the inertial velocity
        lambdav_ii = in_cond[7]
        lambdagamma_ii = in_cond[8]
        lambdah_ii = in_cond[9]

        # Assign Parameters
        ω_planet = m.planet.ω
        γ = m.planet.γ
        area_tot = m.body.area_tot

        # TRANSFORM THE STATE
        # Inertial to planet relative transformation
        pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, et) # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
        pos_pp_mag = norm(pos_pp) # Magnitude of the planet relative position
        vel_pp_mag = norm(vel_pp)

        # Orbital Elements
        OE = rvtoorbitalelement(pos_ii, vel_ii, mass, m.planet)

        # Timing variables
        el_time = value(seconds((date_initial + t0*seconds) - from_utc(DateTime(args[:year], args[:month], args[:day], args[:hours], args[:minutes], args[:secs])))) # Elapsed time since the beginning of the simulation
        current_time =  value(seconds(date_initial + t0*seconds - TAIEpoch(2000, 1, 1, 12, 0, 0.0))) # current time in seconds since J2000
        time_real_utc = to_utc(time_real) # Current time in UTC as a DateTime object
        et = utc2et(time_real_utc) # Current time in Ephemeris Time

        # Angular Momentum Calculations
        h_ii = cross(pos_ii, vel_ii)        # Inertial angular momentum vector[m ^ 2 / s]
        h_ii_mag = norm(h_ii)               # Magnitude of the inertial angular momentum
        h_pp = cross(pos_pp, vel_pp)        # Planet relative angular momentum vector
        h_pp_mag = norm(h_pp)               # Magnitude of the planet relative angular momentum
        h_pp_hat = h_pp / h_pp_mag          # Unit vector of the planet relative angular momentum

        # Inertial flight path angle
        arg = median([-1, 1, h_ii_mag / (pos_ii_mag * vel_ii_mag)])
        γ_ii = acos(arg)

        if dot(pos_ii, vel_ii) < 0
            γ_ii = -γ_ii
        end

        # Relative flight path angle
        arg = median([-1, 1, h_pp_mag / (pos_pp_mag * vel_pp_mag)])
        γ_pp = acos(arg)

        if dot(pos_pp, vel_pp) < 0
            γ_pp = -γ_pp
        end

        # Derived Quantity Calculations

        # Compute Latitude and Longitude
        LatLong = rtolatlong(pos_pp, m.planet)
        lat = LatLong[2]
        lon = LatLong[3]
        alt = LatLong[1]

        # Compute NED basis unit vectors
        uDuNuE = latlongtoNED(LatLong)  # nd
        uD = uDuNuE[1]
        uE = uDuNuE[3]
        uN = uDuNuE[2]

        # Get density, pressure , temperature and winds
        if ip.dm == 0
            ρ, T_p, wind = density_constant(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 1
            ρ, T_p, wind = density_exp(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 2
            ρ, T_p, wind = density_no(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 3
            ρ, T_p, wind = density_gram(alt, m.planet, lat, lon, MonteCarlo, wind_m, args, el_time, gram_atmosphere, gram)
            ρ, T_p, wind = pyconvert(Any, ρ), pyconvert(Any, T_p), [pyconvert(Any, wind[1]), pyconvert(Any, wind[2]), pyconvert(Any, wind[3])]
        end

        # Mach Number
        sound_velocity = sqrt(γ * m.planet.R * T_p)
        Mach = vel_pp_mag / sound_velocity
        S = sqrt(γ/2) * Mach   # molecular speed ratio

        
        lambda_switch = (2.0 * m.body.mass * vel_ii_mag) / (area_tot * CD_slope * pi)

        if lambdav_ii < lambda_switch
            aoa = 0.0
        else
            aoa = m.aerodynamics.α
        end

        # if lambdav_ii > 0
        #     aoa = 0.0
        # else
        #     aoa = m.aerodynamics.α
        # end

        # Convert wind to pp(PCPF) frame
        wE = wind[1] # positive to the east , m / s
        wN = wind[2] # positive to the north , m / s
        wU = wind[3] # positive up , m / s

        wind_pp = wN * uN + wE * uE - wU * uD        # wind velocity in pp frame , m / s
        vel_pp_rw = vel_pp + wind_pp                 # relative wind vector , m / s
        vel_pp_rw_hat = vel_pp_rw / norm(vel_pp_rw)  # relative wind unit vector , nd

        # Dynamic pressure, CHANGE THE VELOCITY WITH THE WIND VELOCITY
        q = 0.5 * ρ * norm(vel_pp_rw)^2            # base on wind - relative velocity

        # Heat Rate
        heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, aoa)

        # Add the control for the heat rate if flash == 3
        if heat_rate_control == true && heat_rate > args[:max_heat_rate]
            state = [T_p, ρ, S]
            index_ratio = [1]
            aoa_hr = control_solarpanels_heatrate(ip, m, args, index_ratio, state)

            if args[:struct_ctrl] == 1
                α_struct = control_struct_load(ip, m, args, S, T_p, q, MonteCarlo)

                aoa = min(aoa_hr, α_struct) # limit the angle of attack to the structural load control
            end

            if aoa != aoa_hr
                heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, aoa)
            else
                aoa = aoa_hr
                heat_rate = args[:max_heat_rate]  
            end

            # heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, aoa)
        end

        # Rotation Calculation
        L_PI = pxform("J2000", "IAU_"*uppercase(m.planet.name), current_time)*m.planet.J2000_to_pci'
        # rot_angle = norm(ω_planet) * t0     # rad
        # L_PI = [cos(rot_angle)  sin(rot_angle)  0.0;
        #         -sin(rot_angle) cos(rot_angle)  0.0; 
        #         0.0             0.0             1.0]    # rotation matrix

        if ip.gm == 0
            gravity_ii = mass * gravity_const(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        elseif ip.gm == 1
            gravity_ii = mass * gravity_invsquared(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        elseif ip.gm == 2
            gravity_ii = mass * gravity_invsquared_J2(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        elseif ip.gm == 3
            gravity_ii = mass * gravity_GRAM(pos_ii, lat, lon, alt, m.planet, mass, vel_ii, el_time, gram_atmosphere, args, gram)
        end

        if length(args[:n_bodies]) != 0

            for k = 1:length(args[:n_bodies])  
                gravity_ii += mass * gravity_n_bodies(et, pos_ii, m.planet, config.cnf.n_bodies_list[k])
            end
        end

        if args[:srp] == true
            p_srp_unscaled = 4.56e-6  # N / m ^ 2, solar radiation pressure at 1 AU
            srp_ii = mass * srp(m.planet, p_srp_unscaled, m.aerodynamics.reflection_coefficient, m.body.area_tot, m.body.mass, pos_ii, et)
        end

        bank_angle = 0.0
        lift_pp_hat = cross(h_pp_hat,vel_pp_rw_hat)     # perpendicular vector to angular vector and velocity

        # Vehicle Aerodynamic Forces
        # CL and CD
        CL, CD = aerodynamic_coefficient_fM(aoa, m.body, T_p, S, m.aerodynamics, 0)

        # Force calculations
        drag_pp_hat = -vel_pp_rw_hat                    # Planet relative drag force direction

        drag_pp = q * CD * area_tot * drag_pp_hat                       # PLanet relative drag force vector
        lift_pp = q * CL * area_tot * lift_pp_hat * cos(bank_angle)     # PLanet relative lift force vector

        drag_ii = L_PI' * drag_pp                                       # Inertial drag force vector
        lift_ii = L_PI' * lift_pp                                       # Inertial lift force vector

        # Total Force
        # Total inertial external force vector on body [N]
        force_ii = drag_ii + lift_ii + gravity_ii

        g_ii = norm(gravity_ii) / mass

        # EOM
        lambdav_dot = -3 * ρ * vel_ii_mag^2 * aoa / pi + 
                      lambdav_ii * (ρ* area_tot * CD * vel_ii_mag) / mass - 
                      lambdagamma_ii * ((ρ * area_tot * CL) / (2 * mass) + g_ii /  vel_ii_mag^2 + 1 / (pos_ii_mag)) - 
                      lambdah_ii * γ_ii
        
        lambdag_dot = lambdav_ii * g_ii - lambdah_ii * vel_ii_mag
        
        lambdah_dot = ρ * vel_ii_mag^3 * aoa/ (pi * m.planet.H) - 
                      lambdav_ii * ((ρ * area_tot *CD * vel_ii_mag^2) / (2 * mass * m.planet.H) + 2 * g_ii * γ_ii/ (pos_ii_mag)) + 
                      lambdagamma_ii * (ρ * area_tot * CL * vel_ii_mag / (2 * mass * m.planet.H) - 2 * g_ii / (pos_ii_mag * vel_ii_mag) + vel_ii_mag / (pos_ii_mag)^2)

        # lambdav_dot = lambdav_ii * (ρ* area_tot * CD * vel_ii_mag) / mass - 
        #               lambdagamma_ii * ((ρ * area_tot * CL) / (2 * mass) + g_ii /  vel_ii_mag^2 + 1 / (pos_ii_mag)) - 
        #               lambdah_ii * γ_ii
        
        # lambdag_dot = lambdav_ii * g_ii - lambdah_ii * vel_ii_mag
        
        # lambdah_dot = -lambdav_ii * ((ρ * area_tot * CD * vel_ii_mag^2) / (2 * mass * m.planet.H) + 2 * g_ii * γ_ii/ (pos_ii_mag)) + 
        #               lambdagamma_ii * (ρ * area_tot * CL * vel_ii_mag / (2 * mass * m.planet.H) - 2 * g_ii / (pos_ii_mag * vel_ii_mag) + vel_ii_mag / (pos_ii_mag)^2)

        y_dot[1:3] = vel_ii
        y_dot[4:6] = force_ii / mass
        y_dot[7] = lambdav_dot
        y_dot[8] = lambdag_dot
        y_dot[9] = lambdah_dot
        y_dot[10] = heat_rate

        # m = param[1]
        # index_phase_aerobraking = param[2]
        # ip = param[3]
        # aerobraking_phase = param[4]
        # t_prev = param[5]
        # date_initial = param[6]
        # time_0 = param[7]
        # args = param[8]
        # initial_state = param[9]
        # gram_atmosphere = param[10]
        # gram = param[11]
        
        # ## Counters
        # # Counter for all along the simulation of all passages
        # config.cnf.count_aerobraking = config.cnf.count_aerobraking + 1
        # passage_number = config.cnf.count_aerobraking
        # # Counter for one entire passage
        # config.cnf.count_dori = config.cnf.count_dori + 1
        # # Counter for one phase
        # config.cnf.count_phase = config.cnf.count_phase + 1

        # t0 = t0 * config.cnf.TU

        # # Clock
        # current_epoch = date_initial + t0*seconds # Precompute the current epoch
        # time_real = DateTime(current_epoch) # date_initial + Second(t0)
        # timereal = ref_sys.clock(Dates.year(time_real), Dates.month(time_real), Dates.day(time_real), Dates.hour(time_real), Dates.minute(time_real), Dates.second(time_real))

        # # Timing variables
        # el_time = value(seconds(current_epoch - m.initial_condition.DateTimeIC)) # Elapsed time since the beginning of the simulation
        # current_time =  value(seconds(current_epoch - m.initial_condition.DateTimeJ2000)) # current time in seconds since J2000
        # time_real_utc = to_utc(time_real) # Current time in UTC as a DateTime object
        # config.cnf.et = utc2et(time_real_utc) # Current time in Ephemeris Time
        # m.planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_"*uppercase(m.planet.name), config.cnf.et))*m.planet.J2000_to_pci' # Construct a rotation matrix from J2000 (Planet-fixed frame 0.0 seconds past the J2000 epoch) to planet-fixed frame
        

        # # Assign state
        # pos_ii = SVector{3, Float64}(in_cond[1:3] * config.cnf.DU)                      # Inertial position 
        # vel_ii = SVector{3, Float64}(in_cond[4:6] * config.cnf.DU / config.cnf.TU)      # Inertial velocity
        # mass = m.body.mass                                          # Mass kg
        # pos_ii_mag = norm(pos_ii)                                  # Magnitude of the inertial position
        # vel_ii_mag = norm(vel_ii)                                  # Magnitude of the inertial velocity

        # lambdav_ii = in_cond[7]
        # lambdagamma_ii = in_cond[8]
        # lambdah_ii = in_cond[9]

        # # Assign parameters
        # ω_planet = m.planet.ω
        # γ = m.planet.γ
        # μ_fluid = m.planet.μ_fluid
        # area_tot = m.body.area_tot

        # # TRANSFORM THE STATE
        # # Inertial to planet relative transformation
        # pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, config.cnf.et) # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
        # pos_pp_mag = norm(pos_pp) # Magnitude of the planet relative position
        # pos_pp_hat = pos_pp / pos_pp_mag # Unit vector of the planet relative position
        # pos_ii_hat = pos_ii / pos_ii_mag # Unit vector of the inertial position

        # vel_pp_mag = norm(vel_pp)
        # vel_pp_hat = vel_pp / vel_pp_mag

        # # Orbital Elements
        # OE = rvtoorbitalelement(pos_ii, vel_ii, mass, m.planet)
        # vi = OE[6]

        # # if vi > 0 && vi < pi/2 && config.cnf.ascending_phase == false
        # #     config.cnf.ascending_phase = true
        # # elseif vi >= pi/2 && vi <= pi && config.cnf.ascending_phase == true && args[:body_shape] == "Blunted Cone"
        # #     config.cnf.ascending_phase = false
        # # end

        # if config.cnf.ascending_phase == true && config.cnf.MarsGram_recall == false
        #     config.cnf.atmospheric_data = Dict()
        # end

        # # Angular Momentum Calculations 
        # h_ii = cross(pos_ii, vel_ii)    # Inertial angular momentum vector [m ^ 2 / s]

        # h_ii_mag = norm(h_ii)           # Magnitude of the inertial angular momentum [m ^ 2 / s]
        # h_pp = cross(pos_pp, vel_pp)
        
        # h_pp_mag = norm(h_pp)
        # h_pp_hat = h_pp / h_pp_mag

        # # Inertial flight path angle 
        # arg = median([-1, 1, h_ii_mag / (pos_ii_mag * vel_ii_mag)])     # limit to[-1, 1]
        # γ_ii = acos(arg)    
        # if dot(pos_ii, vel_ii) < 0
        #     γ_ii = -γ_ii
        # end

        # # Relative flight path angle
        # arg = median([-1, 1, h_pp_mag / (pos_pp_mag * vel_pp_mag)])     # limit to[-1, 1]
        # γ_pp = acos(arg)
        # if dot(pos_pp, vel_pp) < 0
        #     γ_pp = -γ_pp
        # end

        # ## Derived Quantity Calculations

        # # Compute latitude and longitude
        # LatLong = rtolatlong(pos_pp, m.planet, args[:topography_model] == "Spherical Harmonics" && norm(pos_ii) - m.planet.Rp_e < args[:EI] * 1e3)
        
        # lat = LatLong[2]
        # lon = LatLong[3]
        # alt = LatLong[1]

        # # if aerobraking_phase == 2 || aerobraking_phase == 0
        # #     if (pos_ii_mag - m.planet.Rp_e - args[:EI] * 1e3) <= 0 && config.cnf.drag_state == false && config.cnf.ascending_phase == false
        # #         config.cnf.drag_state = true
        # #         config.cnf.time_IEI = t0
        # #     elseif (pos_ii_mag - m.planet.Rp_e >= args[:EI] * 1e3) && config.cnf.drag_state == true && config.cnf.ascending_phase
        # #         config.cnf.drag_state = false
        # #         config.cnf.time_OEI = t0
        # #     end
        # # end

        # if aerobraking_phase == 2 || aerobraking_phase == 0
        #     if args[:control_mode] == 1
        #         x = 120
        #     else
        #         x = 140
        #     end

        #     if (config.cnf.heat_rate_prev > 0.005 || abs(pos_ii_mag - m.planet.Rp_e <= x*1e3)) && config.cnf.sensible_loads == false && config.cnf.ascending_phase == false
        #         config.cnf.sensible_loads = true
        #     elseif config.cnf.heat_rate_prev > 0.005 && config.cnf.sensible_loads == true && config.cnf.ascending_phase
        #         config.cnf.sensible_loads = false
        #     end
        # end

        # # Compute NED basis unit vectors
        # uD, uN, uE = latlongtoNED(LatLong)

        # # copmute azimuth
        # vN = dot(vel_pp, uN)
        # vE = dot(vel_pp, uE)
        # azi_pp = atan(vE, vN)

        # # Get density, pressure , temperature and winds
        # config.cnf.Gram_justrecalled = 0
        # if ip.dm == 0
        #     ρ, T_p, wind = density_constant(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        # elseif ip.dm == 1
        #     ρ, T_p, wind = density_exp(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        # elseif ip.dm == 2
        #     ρ, T_p, wind = density_no(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        # elseif ip.dm == 3
        #     ρ, T_p, wind = density_gram(alt, m.planet, lat, lon, MonteCarlo, wind_m, args, el_time, gram_atmosphere, gram)
        #     ρ, T_p, wind = pyconvert(Float64, ρ), pyconvert(Float32, T_p), SVector{3, Float32}([pyconvert(Float32, wind[1]), pyconvert(Float32, wind[2]), pyconvert(Float32, wind[3])])
        # elseif ip.dm == 4
        #     ρ, T_p, wind = density_nrlmsise(alt, m.planet, lat, lon, MonteCarlo, wind_m, args, time_real)
        # end

        # # Define output.txt containing density data
        # p = 0.0
        # if args[:body_shape] == "Spacecraft"
        #     length_car = m.body.length_SA + m.body.length_SC
        # elseif args[:body_shape] == "Blunted Cone"
        #     length_car = m.body.base_radius * 2
        # end

        # Re = vel_pp_mag * ρ * length_car / μ_fluid  # Reynolds number

        # # Mach Number
        # sound_velocity = sqrt(γ * m.planet.R * T_p)
        # Mach = vel_pp_mag / sound_velocity
        # S = sqrt(γ/2) * Mach    # Molecular speed ratio
        # heat_load = in_cond[8] * config.cnf.MU / config.cnf.TU^2 # * 1e4

        # if config.cnf.drag_state == true
        #     ## Check type of fluid and check if this changes for different planets
        #     Kn = 1.26 * sqrt(γ) * Mach / (Re + 1e-5)
        #     if index_phase_aerobraking == 2
        #         if (alt < 80000) && (config.cnf.index_warning_alt == 0)
        #             println("WARNING: Altitude < 80 km!")
        #         end

        #         config.cnf.index_warning_alt = 1
        #     elseif alt > 100000
        #         config.cnf.index_warning_alt = 0
        #     end

        #     if Kn < 0.1 && config.cnf.index_warning_flow == 0
        #         if Bool(args[:print_res])
        #             println("WARNING: Transitional flow passage!")
        #         end
                
        #         config.cnf.index_warning_flow = 1
        #     elseif Kn >= 0.1
        #         config.cnf.index_warning_flow = 0
        #     end
        # end

        # # Heat rate and Control
        # if (index_phase_aerobraking == 2 || index_phase_aerobraking == 1.75 || index_phase_aerobraking == 2.25) && config.cnf.drag_state && length(config.cnf.initial_position_closed_form) != 0
        #     # evaluates the closed form solution the first time at EI km
        #     if abs(pos_ii_mag - m.planet.Rp_e - args[:EI] * 1e3) <= 1e-2 && (args[:control_mode] == 2 || args[:control_mode] == 3) && config.cnf.time_switch_1 == 0
        #         if ip.cm == 3
        #             control_solarpanels_openloop(ip, m, args, [1,0], [T_p, ρ, S], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, 0, true, gram_atmosphere)
        #         elseif ip.cm == 2
        #             control_solarpanels_heatload(ip, m, args, [1,0], [T_p, ρ, S], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, 0, gram_atmosphere)
        #         elseif ip.cm == 1
        #             control_solarpanels_heatrate(ip, m, args, [1,0], [T_p, ρ, S], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form)
        #         elseif ip.cm == 0
        #             no_control(ip, m, args, [1,0], [T_p, ρ, S], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form)
        #         end
        #     end

        #     if index_phase_aerobraking == 2
        #         if Bool(args[:control_in_loop])
        #             config.cnf.state_flesh1 = [[T_p, ρ, S]]
        #             if ip.cm == 3
        #                 config.cnf.α = control_solarpanels_openloop(ip, m, args, [1,1], config.cnf.state_flesh1[1], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE, true, gram_atmosphere)
        #             elseif ip.cm == 2
        #                 config.cnf.α = control_solarpanels_heatload(ip, m, args, [1,1], config.cnf.state_flesh1[1], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE, gram_atmosphere)
        #             elseif ip.cm == 1
        #                 config.cnf.α = control_solarpanels_heatrate(ip, m, args, [1,1], config.cnf.state_flesh1[1], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE)
        #             elseif ip.cm == 0
        #                 config.cnf.α = no_control(ip, m, args, [1,1], config.cnf.state_flesh1[1], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE)
        #             end
        #         elseif args[:control_in_loop] == false && args[:integrator] == "Julia"
        #             if config.controller.count_controller != config.controller.count_prev_controller && config.controller.stored_state == 0 && t0 != config.controller.prev_time
        #                 push!(config.cnf.state_flesh1, [T_p, ρ, S])

        #                 if config.controller.count_controller == 2
        #                     state = config.cnf.state_flesh1[end]
        #                 else
        #                     state = config.cnf.state_flesh1[end-1]
        #                     deleteat!(config.cnf.state_flesh1, 1)
        #                 end

        #                 config.controller.stored_state = 1
        #                 config.controller.prev_time = time_0

        #                 if ip.cm == 3
        #                     config.cnf.α = control_solarpanels_openloop(ip, m, args, [1,1], state, t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE, true, gram_atmosphere)
        #                 elseif ip.cm == 2
        #                     config.cnf.α = control_solarpanels_heatload(ip, m, args, [1,1], state, t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE, gram_atmosphere)
        #                 elseif ip.cm == 1
        #                     config.cnf.α = control_solarpanels_heatrate(ip, m, args, [1,1], state, t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE)
        #                 elseif ip.cm == 0
        #                     config.cnf.α = no_control(ip, m, args, [1,1], state, t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE)
        #                 end
        #             end
        #         end
        #     end

        #     # Heat Rate 
        #     if ip.tm == 1
        #         heat_rate = heatrate_convective_radiative(S, T_p, m, ρ, vel_pp_mag, config.cnf.α)
        #     elseif ip.tm == 2
        #         heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, config.cnf.α)
        #     end
            
        #     cp = m.planet.γ / (m.planet.γ - 1) * m.planet.R

        #     T_r = 0.0
        # else
        #     T_r = 0.0
        #     heat_rate = 0.0
        # end

        # # Convert wind to pp(PCPF) frame
        # wE = wind[1] # positive to the east , m / s
        # wN = wind[2] # positive to the north , m / s
        # wU = wind[3] # positive up , m / s

        # wind_pp = wN * uN + wE * uE - wU * uD         # wind velocity in pp frame, m / s 
        # vel_pp_rw = vel_pp + wind_pp                  # relative wind vector, m / s
        # vel_pp_rw_hat = vel_pp_rw / norm(vel_pp_rw)   # relative wind unit vector 

        # # Dynamic Pressure, CHANGE THE VELOCITY WITH THE WIND VELOCITY
        # q = 0.5 * ρ * norm(vel_pp_rw)^2               # dynamic pressure based on wind, Pa

        # if args[:struct_ctrl] == 1
        #     α_struct = control_struct_load(ip, m, args, S, T_p, q, MonteCarlo)

        #     config.cnf.α = min(config.cnf.α, α_struct) # limit the angle of attack to the structural load control
        # end

        # lambda_switch = (k_cf * 2.0 * m.body.mass * vel_ii_mag) / (area_tot * CD_slope * pi)

        # if lambdav_ii < lambda_switch
        #     config.cnf.α = 0.0
        # else
        #     config.cnf.α = m.aerodynamics.α
        # end

        # α = config.cnf.α

        # # Heat Rate 
        # if (index_phase_aerobraking == 2 || index_phase_aerobraking == 1.75 || index_phase_aerobraking == 2.25) && config.cnf.drag_state && length(config.cnf.initial_position_closed_form) != 0
            
        #     if ip.tm == 1
        #         heat_rate = heatrate_convective_radiative(S, T_p, m, ρ, vel_pp_mag, config.cnf.α)
        #     elseif ip.tm == 2
        #         heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, config.cnf.α)
        #     end
            
        #     cp = m.planet.γ / (m.planet.γ - 1) * m.planet.R

        #     T_r = 0.0
        # else
        #     T_r = 0.0
        #     heat_rate = 0.0
        # end
               
        # # Nominal gravity calculation
        # if ip.gm == 0
        #     gravity_ii = mass * gravity_const(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        # elseif ip.gm == 1
        #     gravity_ii = mass * gravity_invsquared(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        # elseif ip.gm == 2
        #     gravity_ii = mass * (args[:gravity_harmonics] == 1 ? gravity_invsquared(pos_ii_mag, pos_ii, m.planet, mass, vel_ii) : gravity_invsquared_J2(pos_ii_mag, pos_ii, m.planet, mass, vel_ii))
        # elseif ip.gm == 3
        #     gravity_ii = mass * gravity_GRAM(pos_ii, lat, lon, alt, m.planet, mass, vel_ii, el_time, gram_atmosphere, args, gram)
        # end

        # if length(args[:n_bodies]) != 0
        #     for k = 1:length(args[:n_bodies])  
        #         gravity_ii += mass * gravity_n_bodies(config.cnf.et, pos_ii, m.planet, config.cnf.n_bodies_list[k])
        #     end
        # end

        # srp_ii = zeros(3) # solar radiation pressure vector
        # if args[:srp] == true
        #     p_srp_unscaled = 4.56e-6  # N / m ^ 2, solar radiation pressure at 1 AU
        #     srp_ii = mass * srp(m.planet, p_srp_unscaled, m.aerodynamics.reflection_coefficient, m.body.area_tot, m.body.mass, pos_ii, config.cnf.et)
        # end

        # if args[:gravity_harmonics] == 1
        #     gravity_ii += mass * m.planet.L_PI' * acc_gravity_pines!(pos_pp, m.planet.Clm, m.planet.Slm, args[:L], args[:M], m.planet.μ, m.planet.Rp_e, m.planet)
        # end

        # bank_angle = deg2rad(0.0)

        # lift_pp_hat = cross(h_pp_hat, vel_pp_rw_hat)

        # # Vehicle Aerodynamic Forces
        # # CL and CD
        # if ip.am == 0
        #     CL, CD = aerodynamic_coefficient_constant(α, m.body, T_p, S, m.aerodynamics, MonteCarlo)
        # elseif ip.am == 1
        #     CL, CD = aerodynamic_coefficient_fM(α, m.body, T_p, S, m.aerodynamics, MonteCarlo)
        # elseif ip.am == 2
        #     CL, CD = aerodynamic_coefficient_no_ballistic_flight(α, m.body, args, T_p, S, m.aerodynamics, MonteCarlo)
        # end

        # β = mass / (CD*area_tot)    # ballistic coefficient, kg / m ^ 2

        # # Force Calculation
        # drag_pp_hat = -vel_pp_rw_hat    # PLanet relative drag force direction

        # drag_pp = q * CD * area_tot * drag_pp_hat                       # PLanet relative drag force vector
        # lift_pp = q * CL * area_tot * lift_pp_hat * cos(bank_angle)     # PLanet relative lift force vector

        # drag_ii = m.planet.L_PI' * drag_pp   # Inertial drag force vector
        # lift_ii = m.planet.L_PI' * lift_pp   # Inertial lift force vector

        # # Check if propellant mass is greater than 0 kg
        # if config.cnf.index_propellant_mass == 1
        #     if mass - args[:dry_mass] <= 0.5
        #         config.cnf.index_propellant_mass = 0
        #         m.engines.T = 0

        #         if Bool(args[:print_res])
        #             println("WARNING: No fuel left!")
        #         end
        #     end
        # end

        # # Thrust
        # thrust_pp_mag = 0.0
        
        # # Rodrigues rotation formula to rotate thrust vector of angle phi around angular vector from D direction
        # D_L_per_pp_hat = cross(drag_pp_hat, lift_pp_hat)
        # thrust_pp_hat =  drag_pp_hat * cos(args[:phi]) + cross(D_L_per_pp_hat, drag_pp_hat) * sin(args[:phi]) + D_L_per_pp_hat * dot(D_L_per_pp_hat, drag_pp_hat) * (1 - cos(args[:phi]))
        # #these two ways give the same direction
        # thrust_pp = thrust_pp_mag * thrust_pp_hat
        # thrust_ii = m.planet.L_PI' * thrust_pp

        # # Total Force
        # # Total inetrial external force vector on body [N]
        # force_ii = drag_ii + lift_ii + gravity_ii + thrust_ii + srp_ii

        # g_ii = norm(gravity_ii) / mass

        # # EOM
        # lambdav_dot = -3*ρ*vel_ii_mag^2*α / pi + 
        #               lambdav_ii * (ρ* area_tot * CD * vel_ii_mag) / mass - 
        #               lambdagamma_ii * ((ρ * area_tot * CL) / (2 * mass) + g_ii /  vel_ii_mag^2 + 1 / (pos_ii_mag)) - 
        #               lambdah_ii * γ_ii
        
        # lambdag_dot = lambdav_ii * g_ii - lambdah_ii * vel_ii_mag
        
        # lambdah_dot = ρ * vel_ii_mag^3 * α / (pi * m.planet.H) - 
        #               lambdav_ii * ((ρ * area_tot * CD * vel_ii_mag^2) / (2 * mass * m.planet.H) + 2 * g_ii * γ_ii/ (pos_ii_mag)) + 
        #               lambdagamma_ii * (ρ * area_tot * CL * vel_ii_mag / (2 * mass * m.planet.H) - 2 * g_ii / (pos_ii_mag * vel_ii_mag) + vel_ii_mag / (pos_ii_mag)^2)
        
        # y_dot[1:3] = vel_ii * (config.cnf.TU / config.cnf.DU)
        # y_dot[4:6] = force_ii / mass * (config.cnf.TU^2 / config.cnf.DU) 
        # y_dot[7] = lambdav_dot
        # y_dot[8] = lambdag_dot
        # y_dot[9] = lambdah_dot
        # y_dot[10] = heat_rate * config.cnf.TU^3 / config.cnf.MU # * 1e-4
        # energy = (vel_ii_mag^2)/2 - (m.planet.μ / pos_ii_mag)

        return y_dot
    end

    ## EVENTS
    function out_drag_pass_condition(y, t, integrator)
        m = integrator.p[1]
        args = integrator.p[8]

        norm(y[1:3]) - m.planet.Rp_e - args[:AE]*1e3  # upcrossing
    end
    function out_drag_pass_affect!(integrator)
        # println("entered out_drag_passage_affect! in Eoms.jl")
        config.cnf.t_out_drag_passage = integrator.t
        terminate!(integrator)
    end
    out_drag_pass = ContinuousCallback(out_drag_pass_condition, out_drag_pass_affect!, nothing)

    function shooting_residual!(resid, z, p, param)

        r0 = p[1]
        v0 = p[2]

        hf = p[3]
        vf = p[4]
        γf = p[5]

        m = param[1]
        time_0 = param[7]

        lambdav_0 = z[1]
        lambdagamma_0 = z[2]
        lambdah_0 = z[3]

        v_E = p[6]

        # p1 = z[4]
        # p2 = z[5]
        # p3 = z[6]

        # tf = z[7]

        in_cond = [r0[1], r0[2], r0[3], v0[1], v0[2], v0[3], lambdav_0, lambdagamma_0, lambdah_0, 0.0]

        # Time initialization
        initial_time, final_time = time_0, 1500

        method = Tsit5()
        a_tol = 1e-9
        r_tol = 1e-9

        events = out_drag_pass

        # Run simulation
        prob = ODEProblem(f_ctrl!, in_cond, (initial_time, final_time), param)
        sol = solve(prob, method, abstol=a_tol, reltol=r_tol, callback=events)

        lambdav_fin = sol[7,end]
        lambdagamma_fin = sol[8,end]
        lambdah_fin = sol[9,end]

        r_fin = norm(sol[1:3,end])
        v_fin = norm(sol[4:6,end])

        gamma_fin = asin( sol[1:3,end]' * sol[4:6,end] / (r_fin * v_fin) )

        # m = param[1]
        # index_phase_aerobraking = param[2]
        # ip = param[3]
        # aerobraking_phase = param[4]
        # t_prev = param[5]
        # date_initial = param[6]
        # time_0 = param[7]
        # args = param[8]
        # initial_state = param[9]
        # gram_atmosphere = param[10]
        # gram = param[11]
        # k_cf = param[12]

        # t0 = sol.t[end]

        # # Clock
        # current_epoch = date_initial + t0*seconds # Precompute the current epoch
        # time_real = DateTime(current_epoch) # date_initial + Second(t0)
        # timereal = ref_sys.clock(Dates.year(time_real), Dates.month(time_real), Dates.day(time_real), Dates.hour(time_real), Dates.minute(time_real), Dates.second(time_real))

        # # Timing variables
        # el_time = value(seconds(current_epoch - m.initial_condition.DateTimeIC)) # Elapsed time since the beginning of the simulation
        # current_time =  value(seconds(current_epoch - m.initial_condition.DateTimeJ2000)) # current time in seconds since J2000
        # time_real_utc = to_utc(time_real) # Current time in UTC as a DateTime object
        # et = utc2et(time_real_utc) # Current time in Ephemeris Time
        # m.planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_"*uppercase(m.planet.name), et))*m.planet.J2000_to_pci' # Construct a rotation matrix from J2000 (Planet-fixed frame 0.0 seconds past the J2000 epoch) to planet-fixed frame
        
        # pos_ii = SVector{3, Float64}(sol[1:3,end])       # Inertial position 
        # vel_ii = SVector{3, Float64}(sol[4:6,end])       # Inertial velocity
        # mass = m.body.mass                             # Mass kg
        # pos_ii_mag = norm(pos_ii)   # Magnitude of the inertial position
        # vel_ii_mag = norm(vel_ii)   # Magnitude of the inertial velocity
        # lambdav_ii = sol[7,end]

        # # Assign Parameters
        # ω_planet = m.planet.ω
        # γ = m.planet.γ
        # area_tot = m.body.area_tot

        # # TRANSFORM THE STATE
        # # Inertial to planet relative transformation
        # pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, et) # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
        # pos_pp_mag = norm(pos_pp) # Magnitude of the planet relative position
        # vel_pp_mag = norm(vel_pp)

        # # Orbital Elements
        # OE = rvtoorbitalelement(pos_ii, vel_ii, mass, m.planet)

        # # Timing variables
        # el_time = value(seconds((date_initial + t0*seconds) - from_utc(DateTime(args[:year], args[:month], args[:day], args[:hours], args[:minutes], args[:secs])))) # Elapsed time since the beginning of the simulation
        # current_time =  value(seconds(date_initial + t0*seconds - TAIEpoch(2000, 1, 1, 12, 0, 0.0))) # current time in seconds since J2000
        # time_real_utc = to_utc(time_real) # Current time in UTC as a DateTime object
        # et = utc2et(time_real_utc) # Current time in Ephemeris Time

        # # Angular Momentum Calculations
        # h_ii = cross(pos_ii, vel_ii)        # Inertial angular momentum vector[m ^ 2 / s]
        # h_ii_mag = norm(h_ii)               # Magnitude of the inertial angular momentum
        # h_pp = cross(pos_pp, vel_pp)        # Planet relative angular momentum vector
        # h_pp_mag = norm(h_pp)               # Magnitude of the planet relative angular momentum
        # h_pp_hat = h_pp / h_pp_mag          # Unit vector of the planet relative angular momentum

        # # Inertial flight path angle
        # arg = median([-1, 1, h_ii_mag / (pos_ii_mag * vel_ii_mag)])
        # γ_ii = acos(arg)

        # if dot(pos_ii, vel_ii) < 0
        #     γ_ii = -γ_ii
        # end

        # # Relative flight path angle
        # arg = median([-1, 1, h_pp_mag / (pos_pp_mag * vel_pp_mag)])
        # γ_pp = acos(arg)

        # if dot(pos_pp, vel_pp) < 0
        #     γ_pp = -γ_pp
        # end

        # # Derived Quantity Calculations

        # # Compute Latitude and Longitude
        # LatLong = rtolatlong(pos_pp, m.planet)
        # lat = LatLong[2]
        # lon = LatLong[3]
        # alt = LatLong[1]

        # # Compute NED basis unit vectors
        # uDuNuE = latlongtoNED(LatLong)  # nd
        # uD = uDuNuE[1]
        # uE = uDuNuE[3]
        # uN = uDuNuE[2]

        # # Get density, pressure , temperature and winds
        # if ip.dm == 0
        #     ρ, T_p, wind = density_constant(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        # elseif ip.dm == 1
        #     ρ, T_p, wind = density_exp(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        # elseif ip.dm == 2
        #     ρ, T_p, wind = density_no(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        # elseif ip.dm == 3
        #     ρ, T_p, wind = density_gram(alt, m.planet, lat, lon, MonteCarlo, wind_m, args, el_time, gram_atmosphere, gram)
        #     ρ, T_p, wind = pyconvert(Any, ρ), pyconvert(Any, T_p), [pyconvert(Any, wind[1]), pyconvert(Any, wind[2]), pyconvert(Any, wind[3])]
        # end

        # # Mach Number
        # sound_velocity = sqrt(γ * m.planet.R * T_p)
        # Mach = vel_pp_mag / sound_velocity
        # S = sqrt(γ/2) * Mach   # molecular speed ratio

        
        # # lambda_switch = (k_cf * 2.0 * m.body.mass * vel_ii_mag) / (area_tot * CD_slope * pi)

        # # if lambdav_ii < lambda_switch
        # #     aoa = 0.0
        # # else
        # #     aoa = m.aerodynamics.α
        # # end

        # if lambdav_ii > 0
        #     aoa = 0.0
        # else
        #     aoa = m.aerodynamics.α
        # end

        # # Convert wind to pp(PCPF) frame
        # wE = wind[1] # positive to the east , m / s
        # wN = wind[2] # positive to the north , m / s
        # wU = wind[3] # positive up , m / s

        # wind_pp = wN * uN + wE * uE - wU * uD        # wind velocity in pp frame , m / s
        # vel_pp_rw = vel_pp + wind_pp                 # relative wind vector , m / s
        # vel_pp_rw_hat = vel_pp_rw / norm(vel_pp_rw)  # relative wind unit vector , nd

        # # Dynamic pressure, CHANGE THE VELOCITY WITH THE WIND VELOCITY
        # q = 0.5 * ρ * norm(vel_pp_rw)^2            # base on wind - relative velocity

        # # Heat Rate
        # heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, aoa)

        # # Add the control for the heat rate if flash == 3
        # if heat_rate_control == true && heat_rate > args[:max_heat_rate]
        #     state = [T_p, ρ, S]
        #     index_ratio = [1]
        #     aoa_hr = control_solarpanels_heatrate(ip, m, args, index_ratio, state)

        #     if args[:struct_ctrl] == 1
        #         α_struct = control_struct_load(ip, m, args, S, T_p, q, MonteCarlo)

        #         aoa = min(aoa_hr, α_struct) # limit the angle of attack to the structural load control
        #     end

        #     if aoa != aoa_hr
        #         heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, aoa)
        #     else
        #         aoa = aoa_hr
        #         heat_rate = args[:max_heat_rate]  
        #     end

        #     # heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, aoa)
        # end

        # # Rotation Calculation
        # L_PI = pxform("J2000", "IAU_"*uppercase(m.planet.name), current_time)*m.planet.J2000_to_pci'
        # # rot_angle = norm(ω_planet) * t0     # rad
        # # L_PI = [cos(rot_angle)  sin(rot_angle)  0.0;
        # #         -sin(rot_angle) cos(rot_angle)  0.0; 
        # #         0.0             0.0             1.0]    # rotation matrix

        # if ip.gm == 0
        #     gravity_ii = mass * gravity_const(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        # elseif ip.gm == 1
        #     gravity_ii = mass * gravity_invsquared(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        # elseif ip.gm == 2
        #     gravity_ii = mass * gravity_invsquared_J2(pos_ii_mag, pos_ii, m.planet, mass, vel_ii)
        # elseif ip.gm == 3
        #     gravity_ii = mass * gravity_GRAM(pos_ii, lat, lon, alt, m.planet, mass, vel_ii, el_time, gram_atmosphere, args, gram)
        # end

        # if length(args[:n_bodies]) != 0

        #     for k = 1:length(args[:n_bodies])  
        #         gravity_ii += mass * gravity_n_bodies(et, pos_ii, m.planet, config.cnf.n_bodies_list[k])
        #     end
        # end

        # # Vehicle Aerodynamic Forces
        # # CL and CD
        # CL, CD = aerodynamic_coefficient_fM(aoa, m.body, T_p, S, m.aerodynamics, 0)

        # g_ii = norm(gravity_ii) / mass

        # Hf = lambdah_fin * v_fin * gamma_fin - lambdav_fin * (ρ*area_tot*CD*v_fin^2/(2*mass) + g_ii * gamma_fin) + lambdagamma_fin * (ρ*area_tot*CL*v_fin/(2*mass) - g_ii/v_fin + v_fin/r_fin)

        # # Residuals
        # resid[1] = r_fin - m.planet.Rp_e - hf
        # resid[2] = v_fin - vf
        # resid[3] = gamma_fin - γf
        # resid[4] = lambdav_fin - vf
        # resid[5] = lambdagamma_fin - 0.0
        # resid[6] = lambdah_fin - (m.planet.μ / (m.planet.Rp_e + hf)^2)

        # resid[1] = v_fin - vf
        # resid[2] = gamma_fin - γf
        # resid[3] = r_fin - (m.planet.Rp_e + hf)
        # resid[4] = lambdav_fin - (v_fin + p1)
        # resid[5] = (lambdagamma_fin - (0.0 + p2))/1e7
        # resid[6] = lambdah_fin - ((m.planet.μ / (r_fin)^2) + p3)

        resid[1] = lambdav_fin - v_E*v_fin
        resid[2] = lambdagamma_fin - 0.0
        resid[3] = lambdah_fin - v_E*(m.planet.μ / (r_fin)^2)
        # resid[4] = (v_fin^2/2 - (m.planet.μ / r_fin)) - energy_f

        # resid[7] = Hf

        # resid[1] = lambdav_fin - ((v_fin^2/2 - m.planet.μ / r_fin^2 - vf^2/2 + m.planet.μ / (m.planet.Rp_e + hf)^2) *v_fin)
        # resid[2] = lambdagamma_fin - 0.0
        # resid[3] = lambdah_fin - ((v_fin^2/2 - m.planet.μ / r_fin^2 - vf^2/2 + m.planet.μ / (m.planet.Rp_e + hf)^2) * (m.planet.μ / (r_fin)^2))

        return 
    end

    index_phase_aerobraking = nothing
    aerobraking_phase = nothing
    initial_state = nothing

    # z0 = [-550, 80000, 9, -1407.1, -2.30485e7, 7.0444]
    # z0 = [4300, 0.0, 160000]

    z0 = [-500, 80000, 9]

    param = (m, index_phase_aerobraking, ip, aerobraking_phase, t_prev, date_initial, time_0, args, initial_state, gram_atmosphere, gram, k_cf)

    p = (r0, v0, hf, vf, γf, v_E)

    # prob = NonlinearProblem(shooting_residual!, z0, p)
    # sol_NL = solve(prob, NewtonRaphson())

    sol_NL = nlsolve((residuals, z) -> shooting_residual!(residuals, z, p, param), z0, show_trace=true)

    println(sol_NL)

    in_cond = [r0[1], r0[2], r0[3], v0[1], v0[2], v0[3], sol_NL.zero[1], sol_NL.zero[2], sol_NL.zero[3], 0.0]

    # Time initialization
    initial_time, final_time = time_0, time_0 + 1500

    method = Tsit5()
    a_tol = 1e-9
    r_tol = 1e-9

    events = out_drag_pass

    # Run simulation
    prob = ODEProblem(f_ctrl!, in_cond, (initial_time, final_time), param)
    sol = solve(prob, method, abstol=a_tol, reltol=r_tol, callback=events)

    v_ii_mag = [norm(sol[4:6,i]) for i in 1:length(sol.t)]

    # println("v_ii_mag: ", v_ii_mag)
    # println("CD_slope: ", CD_slope)

    lambda_switch_list = (2.0 * m.body.mass * v_ii_mag) ./ (m.body.area_tot * CD_slope * pi)

    # println("lambda_switch_list: ", lambda_switch_list)

    # lambda_switch_list = zeros(size(sol.t))

    # println("lambda_switch_list: ", lambda_switch_list)

    push!(config.cnf.lambda_switch_list, lambda_switch_list...)
    push!(config.cnf.time_switch_list, sol.t...)

    return sol
end