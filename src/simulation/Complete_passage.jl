include("../utils/Reference_system.jl")
include("../integrator/Integrators.jl")
include("../integrator/Events.jl")
# include("../integrator/implicit_midpoint_jacobian.jl")
include("../utils/Save_results.jl")
include("../utils/quaternion_utils.jl")

# include("../physical_models/Gravity_models.jl")
include("../physical_models/Density_models.jl")
include("../physical_models/Aerodynamic_models.jl")
include("../physical_models/Thermal_models.jl")
# include("../physical_models/Perturbations.jl")
# include("../physical_models/DynamicEffectors.jl")

include("../control/Control.jl")
include("../control/utils/Propulsive_maneuvers.jl")
include("../control/targeting_control/targeting.jl")

include("../control/utils/Eom_ctrl.jl")
include("../control/targeting_control/Eom_targeting.jl")
include("../control/targeting_control/sim_targeting.jl")

using LinearAlgebra
using OrdinaryDiffEq
using DiffEqCallbacks
using Dates
using AstroTime
using SPICE
using PythonCall
using StaticArrays
using Quaternions
using Arrow
sys = pyimport("sys")

using .SimulationModel
import .ref_sys
import .quaternion_utils


const R0 = 149597870.7e3 # 1AU, m
const g_e = 9.81 # Gravitational acceleration of Earth at surface, m/s^2
const solarRadFlux = 1368.0 # Solar radiation flux at 1 AU, W/m^2
const speedLight = 299792458.0 # Speed of light, m/s
const AstU = 149597870.7e3 # Astronomical Unit, m



function asim(ip, initial_state, numberofpassage, args, params, gram_atmosphere=nothing, gram=nothing)
    cnf, m, solution = params
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

    T_ijk = SMatrix{3, 3, Float64}([cos(Ω)*cos(ω)-sin(Ω)*sin(ω)*cos(i)   sin(Ω)*cos(ω)+cos(Ω)*sin(ω)*cos(i)    sin(ω)*sin(i);
             -cos(Ω)*sin(ω)-sin(Ω)*cos(ω)*cos(i)  -sin(Ω)*sin(ω)+cos(Ω)*cos(ω)*cos(i)   cos(ω)*sin(i);
             sin(Ω)*sin(i)                        -cos(Ω)*sin(i)                        cos(i)])


    r0, v0 = orbitalelemtorv(OE, m.planet)

    Mass = OE[end]

    # Clock
    date_initial = from_utc(DateTime(m.initial_condition.year, 
                                    m.initial_condition.month,
                                    m.initial_condition.day, 
                                    m.initial_condition.hour, 
                                    m.initial_condition.minute, 
                                    m.initial_condition.second))

    cnf.count_numberofpassage += 1
    t_prev = 0.0
    # if cnf.count_numberofpassage != 1
    #     t_prev = solution.orientation.time[end]
    # else
    #     t_prev = m.initial_condition.time_rot # value(seconds(date_initial - from_utc(DateTime(2000, 1, 1, 12, 0, 0)))) # m.initial_condition.time_rot
    # end

    # sol_lam = asim_ctrl_plot(ip, m, 0, OE, args, 0.03565, true, gram_atmosphere)

    # hf = 160e3 # args[:AE] * 1e3
    # vf = 4195.4809 # 4196.4868
    # γf = 0.10979 # deg2rad(5.874)
    # # energy_f = -3.265e6 
    # energy_f = vf^2 / 2 - (m.planet.μ /(hf + m.planet.Rp_e))

    # println("Targeting energy: ", energy_f)

    # # cnf.targeting = 1

    # # cnf.t_switch_targeting = control_solarpanels_targeting(energy_f, ip, m, 0, OE, args, gram_atmosphere)

    # # println("Targeting switch time: ", cnf.t_switch_targeting)

    # sol_lam = asim_ctrl_targeting_plot(ip, m, 0, OE, args, hf, vf, γf, energy_f, 0.03565, false, gram_atmosphere)
    
    # println("hf: ", norm(sol_lam[1:3,end]) - m.planet.Rp_e)
    # println("vf: ", norm(sol_lam[4:6,end]))
    # println("γf: ", asin(sol_lam[1:3,end]'*sol_lam[4:6,end]/norm(sol_lam[4:6,end]) / norm(sol_lam[1:3,end])))

    # println("Final Energy: ", norm(sol_lam[4:6,end])^2/2 - m.planet.μ/norm(sol_lam[1:3,end]))

    # push!(cnf.time_list, sol_lam.t...)
    # push!(cnf.lamv_list, sol_lam[7,:]...)

    function f!(y_dot, in_cond, param, t0::Float64)
        m = param.m
        index_phase_aerobraking = param.index_phase_aerobraking
        ip = param.ip
        aerobraking_phase = param.aerobraking_phase
        # t_prev = param[5]
        date_initial = param.date_initial
        time_0 = param.time_0
        args = param.args
        initial_state = param.initial_state
        gram_atmosphere = param.gram_atmosphere
        gram = param.gram
        numberofpassage = param.numberofpassage
        orientation_sim = param.orientation_sim
        # args = param[16]

        # Orbital Elements
        pos_ii = SVector{3, Float64}((@view in_cond[1:3]) * cnf.DU)                      # Inertial position 
        vel_ii = SVector{3, Float64}((@view in_cond[4:6]) * cnf.DU / cnf.TU)      # Inertial velocity
        mass = in_cond[7] * cnf.MU                                          # Mass kg
        # println("pos_ii: ", pos_ii)
        OE = rvtoorbitalelement(pos_ii, vel_ii, mass, m.planet)
        vi = OE[6]
        
        ## Counters

        t0 *= cnf.TU

        # Clock
        current_epoch = date_initial + (t0-m.initial_condition.el_time)*seconds # Precompute the current epoch
        time_real = DateTime(current_epoch) # date_initial + Second(t0)
        timereal = ref_sys.clock(Dates.year(time_real), Dates.month(time_real), Dates.day(time_real), Dates.hour(time_real), Dates.minute(time_real), Dates.second(time_real))

        # Timing variables
        el_time = value(seconds(current_epoch - m.initial_condition.DateTimeIC)) # Elapsed time since the beginning of the simulation
        current_time =  value(seconds(current_epoch - m.initial_condition.DateTimeJ2000)) # current time in seconds since J2000
        time_real_utc = to_utc(time_real) # Current time in UTC as a DateTime object
        cnf.et = utc2et(time_real_utc) # Current time in Ephemeris Time
        m.planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_"*uppercase(m.planet.name), cnf.et))*m.planet.J2000_to_pci' # Construct a rotation matrix from J2000 (Planet-fixed frame 0.0 seconds past the J2000 epoch) to planet-fixed frame
        

        # Assign state
        # println(in_cond)
        # quaternion = SVector{4, Float64}(in_cond[1:4]) # Quaternion
        # pos_ii = SVector{3, Float64}(in_cond[1:3] * cnf.DU)                      # Inertial position
        # vel_ii = SVector{3, Float64}(in_cond[4:6] * cnf.DU / cnf.TU)      # Inertial velocity
        # mass = in_cond[7] * cnf.MU                                          # Mass kg
        # ω = SVector{3, Float64}(in_cond[9:11] / cnf.TU)                # Angular velocity vector [rad / s]
        
        quat_idx = 8 + length(m.body.links)
        if orientation_sim
            quaternion = SVector{4, Float64}(@view in_cond[quat_idx:quat_idx+3]) # Quaternion
            ω = SVector{3, Float64}((@view in_cond[quat_idx+4:quat_idx+6]) / cnf.TU)                # Angular velocity vector [rad / s]
            m.body.roots[1].q .= quaternion
            # quaternion = SVector{4, Float64}(m.body.roots[1].q)
            m.body.roots[1].ω .= ω # Body frame angular velocity
        else
            quaternion = orbital_elements_to_lvlh_quaternion(OE[4], OE[3], OE[5], OE[6]) # Quaternion, set to align with LVLH frame if orientation simulation is not enabled
            ω = SVector{3, Float64}(0.0, 0.0, 0.0)                # Angular velocity vector [rad / s]
        end

        rot_body_to_inertial = rot(quaternion)' # Rotation matrix from body frame to inertial frame

        pos_ii_mag = norm(pos_ii)                                  # Magnitude of the inertial position
        vel_ii_mag = norm(vel_ii)                                  # Magnitude of the inertial velocity

        # Define unnormalized state
        state_vector = SVector{14, Float64}(pos_ii..., vel_ii..., mass, quaternion..., ω...) # State vector

        # Assign parameters
        ω_planet = m.planet.ω
        γ = m.planet.γ
        μ_fluid = m.planet.μ_fluid
        bodies, root_index = traverse_bodies(m.body, m.body.roots[1]) # Get all bodies in the simulation

        # TRANSFORM THE STATE
        # Inertial to planet relative transformation
        pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, cnf.et) # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
        pos_pp_mag = norm(pos_pp) # Magnitude of the planet relative position

        vel_pp_mag = norm(vel_pp)

        

        Mars_Gram_recalled_at_periapsis = false

        if vi > 0 && vi < pi/2 && cnf.ascending_phase == false
            cnf.ascending_phase = true
        elseif vi >= pi/2 && vi <= pi && cnf.ascending_phase == true && args[:body_shape] == "Blunted Cone"
            cnf.ascending_phase = false
        end

        if cnf.ascending_phase == true && cnf.MarsGram_recall == false
            cnf.atmospheric_data = Dict()
        end

        # Angular Momentum Calculations 
        h_ii = cross(pos_ii, vel_ii)    # Inertial angular momentum vector [m ^ 2 / s]

        h_ii_mag = norm(h_ii)           # Magnitude of the inertial angular momentum [m ^ 2 / s]
        h_pp = cross(pos_pp, vel_pp)
        
        h_pp_mag = norm(h_pp)
        h_pp_hat = normalize(h_pp) # Unit vector of the planet relative angular momentum
        # param[14][:] .= h_pp_hat # Update the angular momentum unit vector in the parameter array for later use
        # Inertial flight path angle 
        γ_ii = acos(clamp(h_ii_mag / (pos_ii_mag * vel_ii_mag), -1.0, 1.0))     # limit to[-1, 1]
        # γ_ii = acos(arg)    
        if dot(pos_ii, vel_ii) < 0.0
            γ_ii = -γ_ii
        end

        # Relative flight path angle
        γ_pp = acos(clamp(h_pp_mag / (pos_pp_mag * vel_pp_mag), -1.0, 1.0))     # limit to[-1, 1]
        # γ_pp = acos(arg)
        if dot(pos_pp, vel_pp) < 0.0
            γ_pp = -γ_pp
        end

        ## Derived Quantity Calculations

        # Compute latitude and longitude
        alt,lat,lon = rtolatlong(pos_pp, m.planet, args[:topography_model] == "Spherical Harmonics" && norm(pos_ii) - m.planet.Rp_e < args[:EI] * 1e3)
        
        # alt,lat,lon = LatLong
        
        # println(" ")
        # println(" Altitude: ", alt)
        # println(" ")

        if aerobraking_phase == 2 || aerobraking_phase == 0
            if (pos_ii_mag - m.planet.Rp_e - args[:EI] * 1.0e3) <= 0.0 && cnf.drag_state == false && cnf.ascending_phase == false
                cnf.drag_state = true
                cnf.time_IEI = t0
            elseif (pos_ii_mag - m.planet.Rp_e >= args[:EI] * 1.0e3) && cnf.drag_state == true && cnf.ascending_phase
                cnf.drag_state = false
                cnf.time_OEI = t0
            end
        end


        if aerobraking_phase == 2 || aerobraking_phase == 0
            if args[:control_mode] == 1
                x = 120.0
            else
                x = 140.0
            end

            if (any(i -> i > 0.005, cnf.heat_rate_prev) || abs(pos_ii_mag - m.planet.Rp_e <= x*1.0e3)) && cnf.sensible_loads == false && cnf.ascending_phase == false
                cnf.sensible_loads = true
            elseif any(i -> i > 0.005, cnf.heat_rate_prev) && cnf.sensible_loads == true && cnf.ascending_phase
                cnf.sensible_loads = false
            end
        end

        # Compute NED basis unit vectors
        uD, uN, uE = latlongtoNED([alt,lat,lon])

        # copmute azimuth
        vN = dot(vel_pp, uN)
        vE = dot(vel_pp, uE)
        azi_pp = atan(vE, vN)

        # Get density, pressure , temperature and winds
        cnf.Gram_justrecalled = 0
        if ip.dm == 0
            ρ, T_p, wind = density_constant(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 1
            ρ, T_p, wind = density_exp(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 2
            ρ, T_p, wind = density_no(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 3
            ρ, T_p, wind = density_gram(alt, m.planet, lat, lon, MonteCarlo, wind_m, args, el_time, param, gram_atmosphere, gram)
        elseif ip.dm == 4
            ρ, T_p, wind = density_nrlmsise(alt, m.planet, lat, lon, MonteCarlo, wind_m, args, time_real)
        end

        # Define output.txt containing density data
        p = 0.0
        if args[:body_shape] == "Spacecraft"
            length_car = get_spacecraft_length(m.body, m.body.roots[1]) # Length of the spacecraft
        elseif args[:body_shape] == "Blunted Cone"
            length_car = m.body.base_radius * 2.0
        end

        Re = vel_pp_mag * ρ * length_car / μ_fluid  # Reynolds number

        # Mach Number
        sound_velocity = sqrt(γ * m.planet.R * T_p)
        Mach = vel_pp_mag / sound_velocity
        S = sqrt(γ/2.0) * Mach    # Molecular speed ratio
        # param[18] .= [ρ, T_p, S] # Update the density, temperature and speed ratio in the parameter array for later use
        heat_load = in_cond[8:8+length(bodies)-1] * cnf.MU / cnf.TU^2 # * 1e4

        if cnf.drag_state == true
            ## Check type of fluid and check if this changes for different planets
            Kn = 1.26 * sqrt(γ) * Mach / (Re + 1.0e-5)
            if index_phase_aerobraking == 2
                if (alt < 80000.0) && (cnf.index_warning_alt == 0)
                    println("WARNING: Altitude < 80 km!")
                end

                cnf.index_warning_alt = 1
            elseif alt > 100000.0
                cnf.index_warning_alt = 0
            end

            if Kn < 0.1 && cnf.index_warning_flow == 0
                if Bool(args[:print_res])
                    println("WARNING: Transitional flow passage!")
                end
                
                cnf.index_warning_flow = 1
            elseif Kn >= 0.1
                cnf.index_warning_flow = 0
            end
        end

        if isempty(cnf.heat_load_past)
            cnf.heat_load_past = zeros(length(bodies))
        end
        cnf.heat_load_past .= heat_load

        # println([cnf.drag_state, length(cnf.initial_position_closed_form)])

        # Heat rate and Control
        if (index_phase_aerobraking == 2 || index_phase_aerobraking == 1.75 || index_phase_aerobraking == 2.25) && cnf.drag_state && cnf.initial_position_closed_form[1] != 0
            # evaluates the closed form solution the first time at EI km
            if abs(pos_ii_mag - m.planet.Rp_e - args[:EI] * 1.0e3) <= 1.0e-2 && (args[:control_mode] == 2 || args[:control_mode] == 3) && cnf.time_switch_1 == 0
                if ip.cm == 3
                    control_solarpanels_openloop(ip, m, args, [1,0], [T_p, ρ, S], t0 - cnf.time_IEI, cnf.initial_position_closed_form, 0, true, gram_atmosphere)
                elseif ip.cm == 2
                    control_solarpanels_heatload(ip, m, args, [1,0], [T_p, ρ, S], t0 - cnf.time_IEI, cnf.initial_position_closed_form, 0, gram_atmosphere)
                elseif ip.cm == 1
                    control_solarpanels_heatrate(ip, m, args, [1,0], [T_p, ρ, S], t0 - cnf.time_IEI, cnf.initial_position_closed_form)
                elseif ip.cm == 0
                    no_control(ip, m, args, [1,0], [T_p, ρ, S], t0 - cnf.time_IEI, cnf.initial_position_closed_form)
                end
            end

            if index_phase_aerobraking == 2
                if Bool(args[:control_in_loop])
                    cnf.state_flesh1 = [[T_p, ρ, S]]
                    if ip.cm == 3
                        # println("comp pass: ", OE)
                        # println("comp pass init_cf: ", cnf.initial_position_closed_form)
                        cnf.α = control_solarpanels_openloop(ip, m, args, [1,1], cnf.state_flesh1[1], t0 - cnf.time_IEI, cnf.initial_position_closed_form, OE, true, gram_atmosphere)
                    elseif ip.cm == 2
                        cnf.α = control_solarpanels_heatload(ip, m, args, [1,1], cnf.state_flesh1[1], t0 - cnf.time_IEI, cnf.initial_position_closed_form, OE, gram_atmosphere)
                    elseif ip.cm == 1
                        cnf.α = control_solarpanels_heatrate(ip, m, args, [1,1], cnf.state_flesh1[1], t0 - cnf.time_IEI, cnf.initial_position_closed_form, OE)
                    elseif ip.cm == 0
                        cnf.α = no_control(ip, m, args, [1,1], cnf.state_flesh1[1], t0 - cnf.time_IEI, cnf.initial_position_closed_form, OE)
                    end
                elseif args[:control_in_loop] == false && args[:integrator] == "Julia"
                    if controller.count_controller != controller.count_prev_controller && controller.stored_state == 0 && t0 != controller.prev_time
                        push!(cnf.state_flesh1, [T_p, ρ, S])

                        if controller.count_controller == 2
                            state = cnf.state_flesh1[end]
                        else
                            state = cnf.state_flesh1[end-1]
                            deleteat!(cnf.state_flesh1, 1)
                        end

                        controller.stored_state = 1
                        controller.prev_time = time_0

                        if ip.cm == 3
                            cnf.α = control_solarpanels_openloop(ip, m, args, [1,1], state, t0 - cnf.time_IEI, cnf.initial_position_closed_form, OE, true, gram_atmosphere)
                        elseif ip.cm == 2
                            cnf.α = control_solarpanels_heatload(ip, m, args, [1,1], state, t0 - cnf.time_IEI, cnf.initial_position_closed_form, OE, gram_atmosphere)
                        elseif ip.cm == 1
                            cnf.α = control_solarpanels_heatrate(ip, m, args, [1,1], state, t0 - cnf.time_IEI, cnf.initial_position_closed_form, OE)
                        elseif ip.cm == 0
                            cnf.α = no_control(ip, m, args, [1,1], state, t0 - cnf.time_IEI, cnf.initial_position_closed_form, OE)
                        end
                    end
                end
            end

            # Heat Rate 
            heat_rate = MVector{length(bodies), Float64}(zeros(length(bodies))) # Heat rate vector for each body
            for (i, b) in enumerate(bodies)
                if ip.tm == 1
                    heat_rate[i] = heatrate_convective_radiative(S, T_p, m, ρ, vel_pp_mag, b.α)
                elseif ip.tm == 2
                    heat_rate[i] = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, b.α)
                end
            end
            
            cp = m.planet.γ / (m.planet.γ - 1) * m.planet.R

            T_r = 0.0
        else
            T_r = 0.0
            heat_rate = MVector{length(bodies), Float64}(zeros(length(bodies))) # Heat rate vector for each body
        end

        if isempty(cnf.heat_rate_prev)
            cnf.heat_rate_prev = zeros(length(bodies)) # Initialize heat rate vector if it is empty
        end
        cnf.heat_rate_prev .= heat_rate # save current heat rate

        # Convert wind to pp(PCPF) frame
        wE, wN, wU = wind # positive to the east , m / s
        # wN = wind[2] # positive to the north , m / s
        # wU = wind[3] # positive up , m / s

        wind_pp = wN * uN + wE * uE - wU * uD         # wind velocity in pp frame, m / s 
        vel_pp_rw = vel_pp + wind_pp                  # relative wind vector, m / s
        # param[15] .= vel_pp_rw # Update the relative wind vector in the parameter array for later use
        vel_pp_rw_hat = normalize(vel_pp_rw)   # relative wind unit vector 

        # Dynamic Pressure, CHANGE THE VELOCITY WITH THE WIND VELOCITY
        q = 0.5 * ρ * norm(vel_pp_rw)^2               # dynamic pressure based on wind, Pa

        if args[:struct_ctrl] == 1
            α_struct = control_struct_load(ip, m, args, S, T_p, q, MonteCarlo)

            println("Structural limit angle of attack: ", rad2deg(α_struct))
            println("Angle of attack before structural limit: ", rad2deg(cnf.α))

            cnf.α = min(cnf.α, α_struct) # limit the angle of attack to the structural load control

            println("Final angle of attack after structural limit: ", rad2deg(cnf.α))
            println(" ")
        end

        if cnf.targeting == 1
            # if t0 >= cnf.t_switch_targeting
            if t0 >= cnf.ts_targ_1 && t0 <= cnf.ts_targ_2
                cnf.α = 0
            else
                state = [T_p, ρ, S]
                index_ratio = [1,1]
                cnf.α = control_solarpanels_heatrate(ip, m, args, index_ratio, state, t0 - cnf.t_switch_targeting, cnf.initial_position_closed_form, OE)
            end
        end

        α = cnf.α
        
        # Assumes that the spacecraft is the standard 2 panels one bus
        root = m.body.roots[1]
        # bodies, root_index = traverse_bodies(m.body, root)
        for body in bodies
            if !body.root
                axis = SVector{3, Float64}(abs.(body.r))
                # Rotate the solar panel to the angle α
                rotate_link(body, axis, - α + m.body.roots[root_index].α)
            end
        end

        # Heat Rate 
        if (index_phase_aerobraking == 2 || index_phase_aerobraking == 1.75 || index_phase_aerobraking == 2.25) && cnf.drag_state
            
            # Heat Rate 
            heat_rate = MVector{length(bodies), Float64}(zeros(length(bodies))) # Heat rate vector for each body
            for (i, b) in enumerate(bodies)
                if ip.tm == 1
                    heat_rate[i] = heatrate_convective_radiative(S, T_p, m, ρ, vel_pp_mag, b.α)
                elseif ip.tm == 2
                    heat_rate[i] = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, b.α)
                end
            end
            
            cp = m.planet.γ / (m.planet.γ - 1) * m.planet.R

            T_r = 0.0
        else
            T_r = 0.0
            heat_rate = MVector{length(bodies), Float64}(zeros(length(bodies))) # Heat rate vector for each body
        end

        cnf.heat_rate_prev .= heat_rate # save current heat rate
        
        # Update the force on each link on the spacecraft
        gravity_ii = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize gravity vector
        tau_grav = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize gravity torque vector
        for effector in m.body.dynamic_effectors
            force, torque = DynamicEffectors.calcForceTorque(effector, state_vector, param)
            gravity_ii .+= force
            tau_grav .+= torque
            # gravity_ii, tau_grav .+= DynamicEffectors.calcForceTorque(effector, state_vector)
        end
        # Nominal gravity calculation
        # if ip.gm == 0
        #     gravity_ii += mass * gravity_const(pos_ii_mag, pos_ii, m.planet)
        # elseif ip.gm == 1
        #     gravity_ii += mass * gravity_invsquared(pos_ii_mag, pos_ii, m.planet)
        # elseif ip.gm == 2
        #     gravity_ii += mass * (args[:gravity_harmonics] == 1 ? gravity_invsquared(pos_ii_mag, pos_ii, m.planet) : gravity_invsquared_J2(pos_ii_mag, pos_ii, m.planet))
        # elseif ip.gm == 3
        #     gravity_ii += mass * gravity_GRAM(pos_ii, lat, lon, alt, m.planet, mass, vel_ii, el_time, gram_atmosphere, args, gram)
        # end


        # n-body perturbations
        # if length(args[:n_bodies]) != 0
        #     for k = 1:length(args[:n_bodies])
        #         gravity_ii .+= mass * gravity_n_bodies(cnf.et, pos_ii, m.planet, cnf.n_bodies_list[k])
        #     end
        # end

        # Calculate gravitational harmonics using Pines' method
        # if args[:gravity_harmonics] == 1
        #     gravity_ii .+= mass * m.planet.L_PI' * acc_gravity_pines!(pos_pp, m.planet.Clm, m.planet.Slm, args[:L], args[:M], m.planet.μ, m.planet.Rp_e, m.planet)
        # end

        if orientation_sim
            Rot = [MMatrix{3,3,Float64}(zeros(3, 3)) for i in eachindex(bodies)] # Rotation matrix from the root body to the spacecraft link
            @inbounds for (i, b) in enumerate(bodies)
                Rot[i] .= rotate_to_inertial(m.body, b, root_index) # Rotation matrix from the spacecraft link to the inertial frame
            end
        end
        if args[:srp] == true
            r_sun_planet = m.planet.J2000_to_pci * SVector{3, Float64}(spkpos("SUN", cnf.et, "J2000", "NONE", uppercase(m.planet.name))[1])*1e3 # Vector describing the position of the Sun wrt the planet in J2000 frame
            eclipse_ratio = args[:eclipse] ? eclipse_area_calc(pos_ii, r_sun_planet, m.planet.Rp_e) : 1.0
            numAU = AstU / norm(r_sun_planet - pos_ii) # 1 / Number of Astronomical Units from the Sun
            P_srp = (solarRadFlux / speedLight) * numAU * numAU # Solar radiation pressure at the spacecraft location, N/m^2

            @inbounds for (i, b) in enumerate(bodies)
                # Calculate the position of the spacecraft link in inertial frame
                pos_ii_body = pos_ii + rot_body_to_inertial * b.r # Update the position of the spacecraft link in inertial frame
                sun_direction = normalize(r_sun_planet - pos_ii)
                srp!(m.body, root_index, sun_direction, b, P_srp, eclipse_ratio, orientation_sim)
            end
        end

        if args[:magnetic_field]
            B_ii = get_magnetic_field_dipole(pos_pp, m.planet.L_PI) # Magnetic field vector in inertial frame, T
            # B_ii = get_magnetic_field(time_real, lat, lon, alt, m.planet.L_PI) # Magnetic field vector in inertial frame, T
            @inbounds for (i, b) in enumerate(bodies)
                if !isempty(b.magnets)
                    @inbounds for magnet in b.magnets
                        b.net_torque .+= calculate_magnetic_torque(magnet.m, rot_body_to_inertial' * B_ii)
                        # break
                    end
                end
            end
        end
        bank_angle = deg2rad(0.0)

        lift_pp_hat = normalize(cross(h_pp_hat, vel_pp_rw_hat))
        # lift_pp_hat /= norm(lift_pp_hat) # Normalize the lift vector in planet relative frame
        drag_pp_hat = -vel_pp_rw_hat # Planet relative drag force direction
        cross_pp_hat = cross(drag_pp_hat, lift_pp_hat) # Cross product of the drag and lift vectors in planet relative frame
        
        CL, CD = 0.0, 0.0 # Initialize aerodynamic coefficients
        total_area = 0.0 # Initialize total area
        
        lift_ii = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize inertial lift force vector
        drag_ii = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize inertial drag force vector
        drag_pp = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize planet relative drag force vector
        lift_pp = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize planet relative lift force vector
        α = MVector{length(bodies), Float64}(zeros(length(bodies))) # Initialize angle of attack vector
        β = MVector{length(bodies), Float64}(zeros(length(bodies))) # Initialize sideslip angle vector
        R = MMatrix{3, 3, Float64}(zeros(3, 3)) # Rotation matrix from the root body to the spacecraft link
        # Determine angle of attack (α) and sideslip angle (β)
        # Vehicle Aerodynamic Forces
        # CL and CD
        @inbounds for (i, b) in enumerate(bodies)
            if orientation_sim
                R .= Rot[i] # Rotation matrix from the spacecraft link to the inertial frame
                body_frame_velocity = R' * m.planet.L_PI' * vel_pp_rw # Velocity of the spacecraft link in inertial frame
                
                α_body = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack in radians
                β_body = atan(body_frame_velocity[2], norm([body_frame_velocity[1], body_frame_velocity[3]])) # Sideslip angle in radians
                α[i] = α_body # Angle of attack for the spacecraft link
                β[i] = β_body # Sideslip angle for the spacecraft link
                b.α = α_body
                b.β = β_body
                b.θ = acos(clamp(vel_pp_rw[1]/norm(vel_pp_rw), -1.0, 1.0)) # Elevation angle for the spacecraft link
            else
                # TODO: Change this so that it just uses above code even with orientation_sim = false
                if b.root
                    # if the body is the root body, then the angle of attack is 90 degrees
                    α[i] = pi/2
                    b.α = pi/2 # Angle of attack for the root body
                else
                    body_frame_velocity = rot(b.q) * SVector{3, Float64}(1.0, 0.0, 0.0) # Velocity of the spacecraft link in inertial frame
                    α[i] = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack for the spacecraft link
                    # α[i] = pi/2 # Angle of attack for the spacecraft link, temporary hard code for testing
                    b.α = α[i] # Angle of attack for the spacecraft link
                end
            end
            if ip.am == 0
                CL, CD = aerodynamic_coefficient_constant(α, m.body, T_p, S, m.aerodynamics, MonteCarlo)
            elseif ip.am == 1
                if orientation_sim
                    CL_body, CD_body, CS_body, Cl_body, Cm_body, Cn_body = aerodynamic_coefficient_fM(b, T_p, S, m.aerodynamics, MonteCarlo)
                else
                    # CL_body, CD_body = aerodynamic_coefficient_fM(α[i], m.body, T_p, S, m.aerodynamics, MonteCarlo)
                    CL_body, CD_body, CS_body, Cl_body, Cm_body, Cn_body = aerodynamic_coefficient_fM(b, T_p, S, m.aerodynamics, MonteCarlo)
                end
            elseif ip.am == 2
                CL, CD = aerodynamic_coefficient_no_ballistic_flight(α, m.body, args, T_p, S, m.aerodynamics, MonteCarlo)
            end

            drag_pp_body = q * CD_body * b.ref_area * drag_pp_hat                       # Planet relative drag force vector
            lift_pp_body = q * CL_body * b.ref_area * lift_pp_hat * cos(bank_angle)     # Planet relative lift force vector
            if orientation_sim
                cross_pp_body = q * CS_body * b.ref_area * cross_pp_hat # Planet relative cross force vector
                cross_body = m.planet.L_PI' * cross_pp_body # Inertial cross force vector
            else
                cross_pp_body = SVector{3, Float64}(0.0, 0.0, 0.0) # Planet relative cross force vector
                cross_body = SVector{3, Float64}(0.0, 0.0, 0.0) # Inertial cross force vector
            end

            drag_body = m.planet.L_PI' * drag_pp_body   # Inertial drag force vector
            lift_body = m.planet.L_PI' * lift_pp_body   # Inertial lift force vector

            # Update the force on the spacecraft link
            b.net_force .+= drag_body + lift_body + cross_body # Update the force on the spacecraft link, inertial frame
            # aero_torque = q * Cl_body * b.ref_area * b.dims[1] * SVector{3, Float64}(1.0, 0.0, 0.0) + # Aerodynamic roll torque, body frame
            #               q * Cm_body * b.ref_area * b.dims[2] * SVector{3, Float64}(0.0, 1.0, 0.0) + # Aerodynamic pitch torque, body frame
            #               q * Cn_body * b.ref_area * b.dims[3] * SVector{3, Float64}(0.0, 0.0, 1.0)   # Aerodynamic yaw torque, body frame
            # b.net_torque .+= aero_torque # Update the torque on the spacecraft link, body frame
            b.net_torque .+= cross(b.r, rot_body_to_inertial' * (drag_body + lift_body + cross_body)) # Update the torque on the spacecraft link, body frame
            # Update the total CL/CD
            CL += CL_body * b.ref_area
            CD += CD_body * b.ref_area
            total_area += b.ref_area # Update the total area
            drag_ii += drag_body # Update the total drag force
            lift_ii += lift_body # Update the total lift force
            drag_pp += drag_pp_body # Update the total drag force in planet relative frame
            lift_pp += lift_pp_body # Update the total lift force in planet relative frame
        end
        
        # Normalize the aerodynamic coefficients
        CL = CL / total_area
        CD = CD / total_area

        # Check if propellant mass is greater than 0 kg
        if cnf.index_propellant_mass == 1
            if mass - get_spacecraft_mass(m.body, m.body.roots[1], dry=true) <= 0.5
                cnf.index_propellant_mass = 0
                m.engines.T = 0

                if Bool(args[:print_res])
                    println("WARNING: No fuel left!")
                end
            end
        end

        # Thrust
        Δv = g_e * m.engines.Isp * log(initial_state.m / mass)

        if ip.tc == 0
            thrust_pp_mag = no_maneuver(t0, m.engines.T, Δv, args, index_phase_aerobraking)
        elseif ip.tc == 1
            thrust_pp_mag = abms(t0, m.engines.T, Δv, args, index_phase_aerobraking)
        elseif ip.tc == 2
            thrust_pp_mag = deceleration_drag_passage(t0, m.engines.T, Δv, args, index_phase_aerobraking)
        end

        # Rodrigues rotation formula to rotate thrust vector of angle phi around angular vector from D direction
        D_L_per_pp_hat = cross(drag_pp_hat, lift_pp_hat)
        thrust_pp_hat = normalize(drag_pp_hat * cos(args[:phi]) + cross(D_L_per_pp_hat, drag_pp_hat) * sin(args[:phi]) + D_L_per_pp_hat * dot(D_L_per_pp_hat, drag_pp_hat) * (1 - cos(args[:phi])))
        #these two ways give the same direction
        thrust_pp = thrust_pp_mag * thrust_pp_hat
        thrust_ii = m.planet.L_PI' * thrust_pp

        # Attitude control torques
        τ_rw = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize reaction wheel torque vector
        total_rw_h = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize total reaction wheel angular momentum vector
        rw_h = MVector{m.body.n_reaction_wheels, Float64}(zeros(m.body.n_reaction_wheels)) # Initialize vector of reaction wheel angular momentum magnitudes
        rw_τ = MVector{m.body.n_reaction_wheels, Float64}(zeros(m.body.n_reaction_wheels)) # Initialize vector of reaction wheel torque magnitudes
        thruster_forces = MVector{m.body.n_thrusters, Float64}(zeros(m.body.n_thrusters)) # Initialize vector of thruster forces
        thruster_fuel_mass_consumption = 0.0 # Tracks mass consumption rate of all attitude control thrusters
        if orientation_sim
            # ω_wheel_derivatives = MVector{m.body.n_reaction_wheels, Float64}(zeros(m.body.n_reaction_wheels)) # Initialize vector of reaction wheel angular momentum derivatives
            counter = 1 # Counter for reaction wheel angular momentum vector
            counter_thrusters = 1 # Counter for thruster forces vector
            @inbounds for (i, b) in enumerate(bodies)
                n_wheels = b.rw_assembly.n_wheels # Number of reaction wheels on the body
                R .= Rot[i] # Rotation matrix from the spacecraft link to the inertial frame
                if n_wheels != 0 # If the body has reaction wheels
                    # Determine the angular momentum derivatives of the reaction wheels
                    τ = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize reaction wheel torque vector
                    clamp!(b.rw_assembly.h_wheels, -b.rw_assembly.max_wheel_h, b.rw_assembly.max_wheel_h) # Clamp the reaction wheel angular momentum to the maximum angular momentum
                    @inbounds for j in 1:n_wheels
                        if abs(b.rw_assembly.h_wheels[j] - b.rw_assembly.max_wheel_h) == 0.0 && sign(b.ω_wheel_derivatives[j]) == sign(b.rw_assembly.h_wheels[j]) # If the reaction wheel angular momentum is at its maximum and the derivative is in the same direction
                            b.ω_wheel_derivatives[j] = 0.0 # Set the angular momentum derivative to zero if the maximum angular momentum is reached
                        end
                        rw_torque = b.rw_assembly.J_rw[:, j] * b.ω_wheel_derivatives[j] # Update the reaction wheel torque
                        if norm(rw_torque) > b.rw_assembly.max_wheel_torque
                            rw_torque = normalize(rw_torque) * b.rw_assembly.max_wheel_torque # Limit the reaction wheel torque to the maximum torque
                        end
                        τ .+= rw_torque # Sum the reaction wheel torques
                        total_rw_h .+= b.rw_assembly.J_rw[:, j] * b.rw_assembly.h_wheels[j] # Update the total reaction wheel angular momentum
                        rw_h[counter] = b.rw_assembly.h_wheels[j] # Update the reaction wheel angular momentum vector
                        rw_τ[counter] = clamp(b.ω_wheel_derivatives[j], -b.rw_assembly.max_wheel_torque, b.rw_assembly.max_wheel_torque) # Update the reaction wheel torque vector
                        counter += 1 # Increment the counter for the reaction wheel angular momentum vector
                    end
                    b.rw_assembly.tau_body_net .= τ # Save the reaction wheel torque in the body
                    τ_rw .+= τ # Sum the reaction wheel torques in the inertial frame
                    b.net_torque .-= τ # Update the torque on the spacecraft link. Subtract the reaction wheel torque because the reaction torque on the spacecraft is opposite to the reaction wheel torque
                end

                # Attitude control thruster torques and forces
                # Check that the current body has thrusters

                if !isempty(b.thrusters)
                    @inbounds for thrust_idx in eachindex(b.thrusters)
                        thruster = b.thrusters[thrust_idx]
                        # Check if the thruster has a stop firing time and if the current time is greater than or equal to the stop firing time
                        thruster.thrust = clamp(thruster.thrust, 0.0, thruster.max_thrust) # Ensure thruster magnitudes are capped at the max
                        thrust = thruster.thrust # Get the thrust magnitude of the thruster (after thrust factor is applied) thruster.max_thrust * thruster.κ
                        b.net_force .+= R * thruster.direction * thrust # Get the net force in the inertial frame
                        thruster_fuel_mass_consumption -= thrust / (g_e * thruster.Isp)
                        thruster_forces[counter_thrusters] = thrust # Update the thruster forces vector
                        counter_thrusters += 1 # Increment the counter for the thruster forces vector
                        rot_to_body = rotate_to_body(b) # Get the rotation matrix from the link frame to the body frame
                        # println("Body net torque up to thruster $thrust_idx: ", b.net_torque)
                        b.net_torque .+= cross(rot_to_body * thruster.location + b.r, rot_to_body * thruster.direction * thrust)  # Get the net torque in the body frame
                    end
                end
            end
        end

        # Total Force
        # Total inertial external force vector on body [N]
        body_forces = sum([b.net_force for b in bodies]) # Sum of all forces on the spacecraft links
        force_ii = body_forces + gravity_ii + thrust_ii
        
        # Torques
        # τ_body = MVector{3, Float64}(1.0e-3, 2.0e-3, -3.0e-3).*sin(t0/9952) # Initialize torque vector
        τ_body = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize torque vector
        inertia_tensor = get_inertia_tensor(m.body, root_index) # Inertia tensor of the body
        # τ_body = param[end] # Use the saved gravity gradient torque
        if orientation_sim
            # Gravity gradient torque
            # R .= rotate_to_inertial(m.body, m.body.roots[1], root_index) # Rotation matrix from the root body to the spacecraft link
            # println("time: $t0")
            # quat = SVector{4, Float64}([args[:q1_itp](t0), args[:q2_itp](t0), args[:q3_itp](t0), args[:q4_itp](t0)])
            # pos = SVector{3, Float64}([args[:p1_itp](t0), args[:p2_itp](t0), args[:p3_itp](t0)])
            # pos_body = rot(quat)*pos
            rHat = pos_ii / pos_ii_mag # Normalize the position vector in body frame
            rHat_B = rot(quaternion) * rHat # Position of the body in body frame
            τ_body += args[:gravity_gradient] ? cross(rHat_B, 3.0 * m.planet.μ / pos_ii_mag / pos_ii_mag / pos_ii_mag * (inertia_tensor * rHat_B)) : SVector{3, Float64}(0.0, 0.0, 0.0) # Gravity gradient torque
            # println("Time: $t0, GG_Torque: $τ_body")
            # τ_body .+= SVector{3, Float64}([args[:tau_1_itp](t0), args[:tau_2_itp](t0), args[:tau_3_itp](t0)])
            # All other torques

            τ_body += sum([b.net_torque for b in bodies]) # Sum of all torques on the spacecraft links
        end
        
        y_dot[1:3] .= vel_ii * (cnf.TU / cnf.DU) # Position derivative in inertial frame
        y_dot[4:6] .= force_ii / mass * (cnf.TU^2 / cnf.DU) # Velocity derivative in inertial frame
        y_dot[7] = (-norm(thrust_ii) / (g_e * m.engines.Isp) + thruster_fuel_mass_consumption) * cnf.TU / cnf.MU       # mass variation
        y_dot[8:8+length(bodies)-1] .= heat_rate * cnf.TU^3 / cnf.MU # Heat load derivatives
        
        next_index = 8 + length(bodies)

        if orientation_sim
            y_dot[next_index:next_index+3] .= (0.5*Ξ(quaternion)*ω) * cnf.TU  # Quaternion derivative
            y_dot[next_index+4:next_index+6] .= (inertia_tensor\(τ_body - cross(ω, inertia_tensor * ω + total_rw_h))) * cnf.TU^2  # Angular velocity derivative
        end

        energy = (vel_ii_mag^2)/2.0 - (m.planet.μ / pos_ii_mag)

        for b in bodies
            b.net_force .= [0.0, 0.0, 0.0] # Reset the net force on each link
            b.net_torque .= [0.0, 0.0, 0.0] # Reset the net torque on each link
        end
        
        ## SAVE RESULTS
        # TODO: Add intermediate solution to params struct and save in a callback instead of here
        if Bool(cnf.results_save)
            # println("Number of passages in complete_passage: ", numberofpassage)
            param.intermediate_solution = IntermediateSolution(t0, timereal.year, timereal.month, timereal.day, timereal.hour, timereal.minute,
                        timereal.second, numberofpassage, pos_ii, pos_ii_mag, vel_ii, vel_ii_mag, pos_pp, 
                        pos_pp_mag, vel_pp, vel_pp_mag, OE[1:6],
                        lat, lon, alt, γ_ii, γ_pp, h_ii, h_ii_mag, h_pp, h_pp_mag, uD, uE, uN, vN, vE,
                        azi_pp, ρ, T_p, p, wind, CL, CD, S, mass, T_r, 
                        q, gravity_ii, drag_ii, drag_pp, lift_ii, lift_pp, force_ii, τ_body, energy, cnf.index_MonteCarlo, Int64(cnf.drag_state),
                        quaternion, ω, cnf.α, vec(inertia_tensor), τ_rw, α, β, heat_rate, heat_load, rw_h, rw_τ, thruster_forces)
            # if !isempty(cnf.solution_intermediate) && cnf.solution_intermediate[end].time == t0
            #     cnf.solution_intermediate[end] = intermediate_solution
            # else
            #     push!(cnf.solution_intermediate, intermediate_solution)
            # end
        end

        return y_dot
    end

    ## EVENTS: 
    function every_step_condition(y, t, integrator)
        """
        Event function to be run at every step. Used for reaction wheels.
        """
        true
    end

    function save_intermediate_solution_affect!(integrator)
        """
        Event function to save the intermediate solution at every step.
        """
        # println("Saving callback")
        save_results(integrator.p)
    end

    save_intermediate_solution = DiscreteCallback(every_step_condition, save_intermediate_solution_affect!)

    function increment_counters!(integrator)
        # Counter for all along the simulation of all passages
        cnf.count_aerobraking += 1
        # Counter for one entire passage
        cnf.count_dori += 1
        # Counter for one phase
        cnf.count_phase += 1
    end

    increment_counter = DiscreteCallback(true, increment_counters!)
    function time_condition(y, t, integrator)
        """
        Event function to check if the time condition is met.
        """
        cnf = integrator.p.cnf
        # Check if the time is greater than the end time
        if lowercase(args[:type_of_mission]) == "time"
            return t*cnf.TU - args[:mission_time]
        else
            return -1 # Do not terminate if the mission type is not "time"
        end
        # return lowercase(args[:type_of_mission]) == "time" && t >= args[:mission_time] / cnf.TU
    end
    function time_affect!(integrator)
        """
        Event function to end the sim if the time condition is met.
        """
        cnf.time_termination = true
        terminate!(integrator) # Terminate the integrator if the time condition is met
    end

    time_check = ContinuousCallback(time_condition, time_affect!, nothing)

    function save_steps_condition(y, t, integrator)
        """
        Check the length of solution_intermediate to determine if it exceeds the :num_steps_to_save argument
        """
        return length(cnf.solution_intermediate) >= args[:num_steps_to_save]
    end

    function save_steps_affect!(integrator)
        terminate!(integrator) # Stop the integration 
    end

    save_steps_check = DiscreteCallback(save_steps_condition, save_steps_affect!)

    function reaction_wheels_affect!(integrator)
        """
        Event function to update the reaction wheels at every step.
        """
        if args[:orientation_sim]
            m = integrator.p.m
            args = integrator.p.args
            bodies, root_index = traverse_bodies(m.body, m.body.roots[1])
            for b in bodies
                if b.gyro != 0.0
                    reaction_wheel_model!(b, b.rw_τ, integrator.dt*cnf.TU)
                end
            end
        end
    end

    reaction_wheel_update = DiscreteCallback(every_step_condition, reaction_wheels_affect!)

    function run_attitude_controller_condition(y, t, integrator)
        """
        Check if the attitude controller should be run.
        """
        # true
        m = integrator.p.m
        if (m.body.n_reaction_wheels == 0 && m.body.n_thrusters == 0) || !args[:orientation_sim]
            return false # Do not run the attitude controller if there are no reaction wheels or orientation simulation is disabled
        end
        solution_states = integrator.p.intermediate_solution
        wE, wN, wU = solution_states.wind
        uD, uN, uE = solution_states.uD, solution_states.uN, solution_states.uE
        wind_pp = wN * uN + wE * uE - wU * uD
        vel_pp_rw = SVector{3, Float64}(solution_states.vel_pp + wind_pp) # Relative wind vector
        h_pp_hat = SVector{3, Float64}(solution_states.h_pp / solution_states.h_pp_mag) # Relative wind unit vector
        aerobraking_phase = 2 # Aerobraking phase   

        Rot = rotate_to_inertial(m.body, m.body.roots[1], 1)

        # Calculate the wind-relative velocity in the inertial frame
        wind_relative_velocity = m.planet.L_PI' * vel_pp_rw
        # Calculate the orientation quaternion from the inertial x-axis to the wind-relative velocity
        orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), wind_relative_velocity)
        error_quat = error_quaternion(SVector{4, Float64}(m.body.roots[1].q), orientation_quat)
        # Check if the error quaternion is significant
        if acosd(error_quat[4])*2 > 0.01 # Check if the angle is greater than 0.01 degrees
            return true  # Run the attitude controller if the error quaternion is significant
        else
            return false  # Do not run the attitude controller if the error quaternion is small
        end
    end
    function run_attitude_controller!(integrator)
        """
        Event function to run the attitude controller periodically.
        """
        if args[:orientation_sim]
            m = integrator.p.m
            solution_states = integrator.p.intermediate_solution
            wE, wN, wU = solution_states.wind
            uD, uN, uE = solution_states.uD, solution_states.uN, solution_states.uE
            wind_pp = wN * uN + wE * uE - wU * uD
            vel_pp_rw = SVector{3, Float64}(solution_states.vel_pp + wind_pp) # Relative wind vector
            h_pp_hat = SVector{3, Float64}(solution_states.h_pp / solution_states.h_pp_mag) # Relative wind unit vector
            aerobraking_phase = 2 # Aerobraking phase
            bodies, root_index = traverse_bodies(m.body, m.body.roots[1])
            for b in bodies
                if b.rw_assembly.n_wheels != 0 || !isempty(b.thrusters)
                    b.attitude_control_function(m, b, root_index, vel_pp_rw, h_pp_hat, aerobraking_phase, integrator.t * cnf.TU) # Calculate the reaction wheel torque
                end
            end
        end
    end

    attitude_controller = args[:orientation_sim] && (m.body.n_reaction_wheels != 0 || m.body.n_thrusters != 0) ? PeriodicCallback(run_attitude_controller!, m.body.roots[1].attitude_control_rate / cnf.TU) : nothing
    # attitude_controller_orbit = DiscreteCallback(run_attitude_controller_condition, run_attitude_controller!)
    attitude_controller_orbit = args[:orientation_sim] && (m.body.n_reaction_wheels != 0 || m.body.n_thrusters != 0) ? PeriodicCallback(run_attitude_controller!, m.body.roots[1].attitude_control_rate / cnf.TU) : nothing
    function thrust_factor_integrator!(integrator)
        """
        Integrate the thrust factors for the thrusters
        """
        m = integrator.p.m
        bodies, root_index = traverse_bodies(m.body, m.body.roots[1])
        @inbounds for b in bodies
            if !isempty(b.thrusters)
                @inbounds for thruster in b.thrusters
                    if thruster.stop_firing_time == 0.0
                        continue # If there is no stop firing time, skip to the next thruster
                    end
                    if thruster.stop_firing_time < integrator.t * cnf.TU + integrator.dt * cnf.TU
                        thruster.κ *= exp(-thruster.cutoff_frequency * integrator.dt * cnf.TU)
                    elseif thruster.stop_firing_time > integrator.t * cnf.TU + integrator.dt * cnf.TU
                        thruster.κ = 1 + (thruster.κ - 1) * exp(-thruster.cutoff_frequency * integrator.dt * cnf.TU)
                    end
                end
            end
        end
    end
    thrust_factor_integrator = nothing#m.body.n_thrusters != 0 ? DiscreteCallback(every_step_condition, thrust_factor_integrator!) : nothing

    # function run_solar_panel_controller!(integrator)
    #     """
    #     Event function to run the solar panel controller periodically.
    #     """

    #     index_phase_aerobraking = integrator.p[2]  # Index of the aerobraking phase
        
    #     if index_phase_aerobraking == 2 && cnf.drag_state && length(cnf.initial_position_closed_form) != 0
    #         m = integrator.p[1]
    #         t0 = integrator.t * cnf.TU  # Current time in seconds
    #         args = integrator.p[8]  # Arguments passed to the integrator
    #         ρ, T_p, S = integrator.p[18]  # Atmospheric density, temperature, and solar panel area
    #         cnf.state_flesh1 = [[T_p, ρ, S]]
    #         if ip.cm == 3
    #             cnf.α = control_solarpanels_openloop(ip, m, args, [1,1], cnf.state_flesh1[1], t0 - cnf.time_IEI, cnf.initial_position_closed_form, OE, true, gram_atmosphere)
    #         elseif ip.cm == 2
    #             println("Control solar panels with heat load")
    #             cnf.α = control_solarpanels_heatload(ip, m, args, [1,1], cnf.state_flesh1[1], t0 - cnf.time_IEI, cnf.initial_position_closed_form, OE, gram_atmosphere)
    #             # println("control_solarpanels_heatload: ", cnf.α)
    #         elseif ip.cm == 1
    #             cnf.α = control_solarpanels_heatrate(ip, m, args, [1,1], cnf.state_flesh1[1], t0 - cnf.time_IEI, cnf.initial_position_closed_form, OE)
    #         elseif ip.cm == 0
    #             cnf.α = no_control(ip, m, args, [1,1], cnf.state_flesh1[1], t0 - cnf.time_IEI, cnf.initial_position_closed_form, OE)
    #         end
    #     end
    # end

    # solar_panel_controller = ip.cm != 0 ? PeriodicCallback(run_solar_panel_controller!, args[:solar_panel_control_rate] / cnf.TU) : nothing

    # attitude_controller_orbit = PeriodicCallback(run_attitude_controller!, m.body.roots[1].attitude_control_rate / cnf.TU)
    function quaternion_update_affect!(integrator)
        """
        Event function to update the quaternion at every step.
        """
        # quaternion = m.body.roots[1].q
        # ω = integrator.u[13:15] / cnf.TU  # Angular velocity
        # dt = integrator.dt * cnf.TU  # Time step in seconds
        # quaternion_update_function = (du, u, p, t) -> du[:] .= 0.5*Ξ(SVector{4, Float64}(u))*ω
        # # Update the quaternion using the angular velocity
        # prob = ODEProblem(quaternion_update_function, quaternion, (0.0, dt))
        # m.body.roots[1].q .= solve(prob, Tsit5()).u[end]  # Update the quaternion in the body
        # integrator.u[9:12] .= m.body.roots[1].q  # Update the quaternion in the integrator
    end
    quaternion_update = nothing

    function quaternion_normalize_condition(y, t, integrator)
        """
        Event function to be run at every step. Used for quaternion normalization.
        """
        m = integrator.p.m
        quat_idx = length(m.body.links)
        if args[:orientation_sim]
            abs(norm(y[8+quat_idx:8+quat_idx+3]) - 1.0) > args[:a_tol_quaternion]  # Check if the quaternion is not normalized
            # true
        else
            false  # If orientation simulation is not enabled, do not normalize
        end
    end
    function quaternion_normalize_affect!(integrator)
        """
        Event function to update the quaternion at every step.
        """
        m = integrator.p.m
        quat_idx = length(m.body.links)
        normalize!(integrator.u[8+quat_idx:8+quat_idx+3])  # Normalize the quaternion
    end
    quaternion_normalize = args[:orientation_sim] ? DiscreteCallback(quaternion_normalize_condition, quaternion_normalize_affect!) : nothing

    function eventfirststep_condition(y, t, integrator)
        """
        Event function to detect the entry interface downcrossing.
        """
        m = integrator.p.m
        cnf = integrator.p.cnf
        norm(y[1:3]) * cnf.DU - m.planet.Rp_e - args[:EI]*1e3   #  downcrossing
    end
    function eventfirststep_affect!(integrator)
        """
        Event function to terminate the integration at the entry interface downcrossing.
        """
        
        integrator.p.cnf.count_eventfirststep += 1
        terminate!(integrator)
    end
    eventfirststep = ContinuousCallback(eventfirststep_condition, nothing, eventfirststep_affect!)

    function eventfirststep_periapsis_condition(y, t, integrator)
        """
        Event function to detect the periapsis crossing.
        """
        m = integrator.p.m
        cnf = integrator.p.cnf
        pos_ii = SVector{3, Float64}([y[1], y[2], y[3]] * cnf.DU)  # Inertial position
        vel_ii = SVector{3, Float64}([y[4], y[5], y[6]] * cnf.DU / cnf.TU)  # Inertial Velocity

        vi = rvtoorbitalelement(pos_ii, vel_ii, y[7] * cnf.MU, m.planet)[6]

        rad2deg(vi) - 180.0  # downcrossing
    end

    function eventfirststep_periapsis_affect!(integrator)
        """
        Event function to terminate the integration at the periapsis crossing.
        """
        integrator.p.cnf.eventfirststep_periapsis += 1
        terminate!(integrator)
    end
    eventfirststep_periapsis = ContinuousCallback(eventfirststep_periapsis_condition, eventfirststep_periapsis_affect!)

    function eventsecondstep_condition(y, t, integrator)
        """
        Event function to detect the atmospheric exit upcrossing.
        """
        m = integrator.p.m
        norm(y[1:3]) * cnf.DU - m.planet.Rp_e - (args[:AE])*1.0e3   # upcrossing
    end
    function eventsecondstep_affect!(integrator)
        """
        Event function to terminate the integration at the atmospheric exit upcrossing.
        """
        integrator.p.cnf.count_eventsecondstep += 1
        terminate!(integrator)
    end
    eventsecondstep = ContinuousCallback(eventsecondstep_condition, eventsecondstep_affect!, nothing)

    function reached_EI_condition(y, t, integrator)
        """
        Event function to detect the entry interface downcrossing.
        """
        m = integrator.p.m
        args = integrator.p.args
        cnf = integrator.p.cnf
        norm(y[1:3]) * cnf.DU - m.planet.Rp_e - args[:EI]*1.0e3  # downcrossing
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
        m = integrator.p.m
        args = integrator.p.args
        cnf = integrator.p.cnf
        norm(y[1:3]) * cnf.DU - m.planet.Rp_e - args[:AE]*1.0e3  # upcrossing
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
        cnf = integrator.p.cnf
        if abs(norm(y[1:3]) * cnf.DU - m.planet.Rp_e - args[:AE]*1e3) <= 1e-5 
            if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 2
                cnf.α = m.aerodynamics.α
            elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 3
                cnf.α = 0.0
            end
        end

        norm(y[1:3]) * cnf.DU - m.planet.Rp_e - args[:AE]*1e3  # upcrossing
    end
    function out_drag_passage_affect!(integrator)
        integrator.p.cnf.count_out_drag_passage += 1
        terminate!(integrator)
    end
    out_drag_passage = ContinuousCallback(out_drag_passage_condition, out_drag_passage_affect!, nothing)

    function in_drag_passage_condition(y, t, integrator)
        m = integrator.p.m
        args = integrator.p.args
        t_prev = integrator.p.t_prev
        ip = integrator.p.ip
        date_initial = integrator.p.date_initial
        cnf = integrator.p.cnf
        pos_ii = SVector{3, Float64}([y[1], y[2], y[3]] * cnf.DU)  # Inertial position
        vel_ii = SVector{3, Float64}([y[4], y[5], y[6]] * cnf.DU / cnf.TU)  # Inertial velocity
        pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, cnf.et)

        LatLong = rtolatlong(pos_pp, m.planet)

        h0 = LatLong[1]

        if ip.gm == 2
            cond = h0 - args[:EI] * 1e3
            thr = 500
        else
            cond = norm(y[1:3]) * cnf.DU - m.planet.Rp_e - args[:EI]*1e3
            thr = 1e-5
        end

        if abs(cond) <= thr
            controller.guidance_t_eval = collect(t*cnf.TU:1/args[:flash1_rate]:(t*cnf.TU)+1500)
        
            # State definition for control 2, 3 State used by closed-form solution
            pos_ii = SVector{3, Float64}([y[1], y[2], y[3]]) * cnf.DU                     # Inertial position
            vel_ii = SVector{3, Float64}([y[4], y[5], y[6]]) * cnf.DU / cnf.TU     # Inertial velocity
            OE_closedform = rvtoorbitalelement(pos_ii, vel_ii, y[7] * cnf.MU, m.planet)
            cnf.initial_position_closed_form = OE_closedform
        end

        cond  # downcrossing
    end
    function in_drag_passage_affect!(integrator) 
        integrator.p.cnf.count_in_drag_passage += 1
        terminate!(integrator)
    end
    in_drag_passage = ContinuousCallback(in_drag_passage_condition, nothing, in_drag_passage_affect!)

    function in_drag_passage_nt_condition(y, t, integrator)
        m = integrator.p.m
        args = integrator.p.args
        t_prev = integrator.p.t_prev
        ip = integrator.p.ip
        date_initial = integrator.p.date_initial
        cnf = integrator.p.cnf
        pos_ii = SVector{3, Float64}([y[1], y[2], y[3]] * cnf.DU)                    # Inertial position
        vel_ii = SVector{3, Float64}([y[4], y[5], y[6]] * cnf.DU / cnf.TU)    # Inertial velocity
        pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, cnf.et)

        LatLong = rtolatlong(pos_pp, m.planet)#, args[:topography_model] == "Spherical Harmonics" && norm(pos_ii) < m.planet.Rp_e + args[:EI]*1e3)

        h0 = LatLong[1]

        if ip.gm == 2 && !Bool(args[:drag_passage])
            cond = (h0 - args[:EI]*1e3)
            thr = 500
        else
            cond = norm(y[1:3]) * cnf.DU - m.planet.Rp_e - args[:EI]*1e3
            thr = 1e-3
        end

        # println([h0, args[:EI]*1e3])
        # println("In drag passage nt condition check: thr = $thr , cond = $cond and ", length(cnf.initial_position_closed_form))

        if h0 <= args[:EI]*1e3 && cnf.initial_position_closed_form[1] == 0
            # println(h0)
            controller.guidance_t_eval = collect(t*cnf.TU:1/args[:flash1_rate]:(t*cnf.TU)+1500)

            # State definition for control 2, 3 State used by closed-form solution
            OE_closedform = rvtoorbitalelement(pos_ii, vel_ii, y[7] * cnf.MU, m.planet)
            cnf.initial_position_closed_form = OE_closedform
        end

        norm(y[1:3]) * cnf.DU - m.planet.Rp_e - args[:EI]*1e3  # downcrossing
    end
    function in_drag_passage_nt_affect!(integrator)
        integrator.p.cnf.count_in_drag_passage_nt += 1
        nothing
    end
    in_drag_passage_nt = ContinuousCallback(in_drag_passage_nt_condition, nothing, in_drag_passage_nt_affect!)

    function apoapsispoint_condition(y, t, integrator)
        m = integrator.p.m
        pos_ii = SVector{3, Float64}([y[1], y[2], y[3]]) * cnf.DU  # Inertial position
        vel_ii = SVector{3, Float64}([y[4], y[5], y[6]]) * cnf.DU / cnf.TU  # Inertial Velocity

        vi = rvtoorbitalelement(pos_ii, vel_ii, y[7] * cnf.MU, m.planet)[6]
        
        rad2deg(vi) - 180 # upcrossing
    end
    function apoapsispoint_affect!(integrator)
        integrator.p.cnf.count_apoapsispoint += 1
        terminate!(integrator)
    end
    apoapsispoint = ContinuousCallback(apoapsispoint_condition, apoapsispoint_affect!, nothing)

    function periapsispoint_condition(y, t, integrator)
        m = integrator.p.m
        pos_ii = SVector{3, Float64}([y[1], y[2], y[3]]) * cnf.DU  # Inertial position
        vel_ii = SVector{3, Float64}([y[4], y[5], y[6]]) * cnf.DU / cnf.TU  # Inertial Velocity

        vi = rvtoorbitalelement(pos_ii, vel_ii, y[7] * cnf.MU, m.planet)[6]

        rad2deg(vi) - 180  # downcrossing
    end
    function periapsispoint_affect!(integrator)
        cnf = integrator.p.cnf
        cnf.count_periapsispoint += 1
        planet = integrator.p.m.planet
        r_p, _ = r_intor_p!(SVector{3, Float64}(integrator.u[1:3] * cnf.DU), 
                        SVector{3, Float64}(integrator.u[4:6] * cnf.DU / cnf.TU), 
                        planet, 
                        cnf.et)
        r, lat, lon = rtolatlong(r_p, planet, args[:topography_model] == "Spherical Harmonics")
        append!(cnf.altitude_periapsis, r*1e-3)
        append!(cnf.latitude_periapsis, rad2deg(lat))
        append!(cnf.longitude_periapsis, rad2deg(lon))
        nothing
    end
    periapsispoint = ContinuousCallback(periapsispoint_condition, nothing, periapsispoint_affect!)

    function impact_condition(y, t, integrator)
        m = integrator.p.m
        args = integrator.p.args
        cnf = integrator.p.cnf

        if args[:body_shape] == "Blunted Body" || args[:type_of_mission] == "Entry"
            min_alt = 0 * 1e3
        else
            min_alt = 35 * 1e3
        end

        norm(y[1:3]) * cnf.DU - (m.planet.Rp_e + min_alt) # upcrossing and downcrossing
    end
    function impact_affect!(integrator)
        cnf.count_impact += 1
        terminate!(integrator)
    end
    impact = ContinuousCallback(impact_condition, impact_affect!)

    function apoapsisgreaterperiapsis_condition(y, t, integrator)
        m = integrator.p.m
        args = integrator.p.args

        r = y[1:3] * cnf.DU
        v = y[4:6] * cnf.DU / cnf.TU
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
        cnf.count_apoapsisgreaterperiapsis += 1
        terminate!(integrator)
    end
    apoapsisgreaterperiapsis = ContinuousCallback(apoapsisgreaterperiapsis_condition, apoapsisgreaterperiapsis_affect!)

    function stop_firing_condition(y, t, integrator)
        m = integrator.p.m
        args = integrator.p.args
        initial_state = integrator.p.initial_state

        mass = y[7] * cnf.MU
        Δv = (g_e * m.engines.Isp) * log(initial_state.m/mass)
        m.body.prop_mass .= [mass - get_spacecraft_mass(m.body, m.body.roots[1], dry=true)]
        Δv - args[:delta_v]  # upcrossing and downcrossing
    end
    function stop_firing_affect!(integrator)
        cnf.count_stop_firing += 1
        terminate!(integrator)
    end
    stop_firing = ContinuousCallback(stop_firing_condition, stop_firing_affect!)

    # TODO: Update these with some guidance/controller struct
    function guidance_condition(y, t, integrator)
        if controller.stored_state == 1
            controller.count_prev_controller = controller.count_controller
        end

        if t * cnf.TU - controller.t > 1 # t - controller.t > 1
            println("Decrease step size of integration and tolerance")
        end

        if abs(t * cnf.TU - controller.t) <= 1e-8 # abs(t - controller.t) <= 1e-8 
            controller.stored_state = 0
            controller.count_controller += 1
            controller.t = controller.guidance_t_eval[controller.count_controller]
        end

        t * cnf.TU - controller.t  # upcrossing
    end
    function guidance_affect!(integrator)
        nothing
    end
    # guidance = ContinuousCallback(guidance_condition, guidance_affect!, nothing)
    guidance = nothing

    function heat_rate_check_condition(y, t, integrator)
        m = integrator.p.m
        args = integrator.p.args
        cnf = integrator.p.cnf
        if args[:control_mode] == 1
            x = 120
        else
            x = 160
        end

        if abs(norm(y[1:3]) * cnf.DU - m.planet.Rp_e - x*1e3) <= 1e-5
            controller.guidance_t_eval = collect(range(start=t * cnf.TU, stop=(t * cnf.TU)+2500, step=1/args[:flash1_rate]))
        end

        norm(y[1:3]) * cnf.DU - m.planet.Rp_e - x*1e3  # upcrossing and downcrossing
    end
    function heat_rate_check_affect!(integrator)
        integrator.p.cnf.count_heat_rate_check += 1
        terminate!(integrator)
    end
    heat_rate_check = ContinuousCallback(heat_rate_check_condition, heat_rate_check_affect!)

    function heat_load_check_exit_condition(y, t, integrator)
        m = integrator.p.m
        args = integrator.p.args

        if args[:control_mode] == 1
            x = 120
        else
            x = 160
        end

        norm(y[1:3]) * cnf.DU - m.planet.Rp_e - x*1e3  # upcrossing
    end
    function heat_load_check_exit_affect!(integrator)
        integrator.p.cnf.count_heat_load_check_exit += 1
        terminate!(integrator)
    end
    heat_load_check_exit = ContinuousCallback(heat_load_check_exit_condition, heat_load_check_exit_affect!, nothing)

    function final_entry_altitude_reached_condition(y, t, integrator)
        m = integrator.p.m
        args = integrator.p.args
        min_alt = args[:final_altitude]

        r_p, _ = r_intor_p!(SVector{3, Float64}(y[1:3] * cnf.DU), 
                        SVector{3, Float64}(y[4:6] * cnf.DU / cnf.TU), 
                        integrator.p[1].planet, 
                        cnf.et)
        alt, lat, lon = rtolatlong(r_p, integrator.p[1].planet, args[:topography_model] == "Spherical Harmonics")
        alt - min_alt # upcrossing and downcrossing
    end

    function final_entry_altitude_reached_affect!(integrator)
        integrator.p.cnf.count_final_entry_altitude_reached += 1
        terminate!(integrator)
    end
    final_entry_altitude_reached = ContinuousCallback(final_entry_altitude_reached_condition, final_entry_altitude_reached_affect!)

    time_0 = m.initial_condition.el_time
    if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 2
        cnf.α = m.aerodynamics.α
    elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 3
        cnf.α = 0.0
    end

    stop_simulation = false

    # Double check these 2
    save_pre_index = 1
    save_post_index = 1

    cnf.impact = false
    cnf.solution_intermediate = []
    cnf.count_dori = 0
    cnf.atmospheric_data = Dict()
    cnf.previous_atmospheric_data = Dict()
    cnf.ascending_phase = false
    cnf.evaluate_switch_heat_load = false
    cnf.state_inner_boundary_atmosphere = [] # used in Density model for vi def
    cnf.time_IP = cnf.time_OP
    cnf.time_IEI = 0
    cnf.time_OEI = 0
    cnf.time_switch_1 = 0.0
    cnf.time_switch_2 = 0.0

    if args[:heat_load_sol] == 2
        cnf.time_switch_2 = 1000.0
    end

    cnf.timer_revaluation = 0
    cnf.closed_form_solution_off = 1         # used in closed form solution online to run the solution only once
    cnf.initial_position_closed_form = SVector{7, Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)  # used to store the state at IEI for closed-form solution
    cnf.heat_rate_list = []               
    cnf.α_list = []

    # controller.guidance_t_eval = []
    # controller.count_prev_controller = 0
    # controller.count_controller = 1
    # controller.stored_state = 1
    # controller.prev_time = 0
    # controller.t = 0

    cnf.security_mode = false
    cnf.stop_simulation = false
    cnf.results_save = 1
    cnf.drag_state = false
    cnf.α_past = m.aerodynamics.α

    if norm(r0) - m.planet.Rp_e <= args[:EI]*1e3
        cnf.drag_state = true
        cnf.initial_position_closed_form = OE
    end

    cnf.sensible_loads = false
    cnf.counter_integrator = 0

    continue_campaign = false

    # Def initial conditions
    root_body = m.body.roots[1]  # Root body of the mission
    if args[:orientation_sim]
        in_cond = Float64[r0[1], r0[2], r0[3], v0[1], v0[2], v0[3], Mass+1e-10, zeros(length(m.body.links))..., root_body.q[1],
               root_body.q[2], root_body.q[3], root_body.q[4], root_body.ω[1], root_body.ω[2], root_body.ω[3]]
    else
        in_cond = Float64[r0[1], r0[2], r0[3], v0[1], v0[2], v0[3], Mass+1e-10, zeros(length(m.body.links))...]
    end

    # non dimensionalization
    in_cond[1:3] ./= cnf.DU
    in_cond[4:6] .*= cnf.TU / cnf.DU
    in_cond[7] /= cnf.MU
    for i in eachindex(m.body.links)
        in_cond[7 + i] *= cnf.TU^2 / cnf.MU  # Mass of the links
    end
    next_index = 7 + length(m.body.links) + 1
    # in_cond[8] *= cnf.TU^2 / cnf.MU # * 1e4
    if args[:orientation_sim]
        normalize!(in_cond[next_index:next_index+3])  # Quaternion normalization
        in_cond[next_index+4:next_index+6] .*= cnf.TU  # Angular velocity
    end

    # If aerobraking maneuver allowed, add a prephase 0
    range_phase_i = 1
    if args[:thrust_control] == "Aerobraking Maneuver" && args[:keplerian] == false
        range_phase_i = 0
    elseif args[:keplerian] == false && norm(r0) - m.planet.Rp_e <= args[:EI]*1e3
        range_phase_i = 2
    elseif norm(r0) - m.planet.Rp_e >= args[:AE]*1e3 && OE[6] < π - 0.1
        range_phase_i = 3
    end

    index_phase_aerobraking = range_phase_i # for scope reasons
    aerobraking_phase = range_phase_i       # for scope reasons
    method = Tsit5()                        # for scope reasons

    # Solve Equations of Motion 
    for aerobraking_phase in range(range_phase_i, 3)
        index_phase_aerobraking = aerobraking_phase
        # If keplerian, just use step 3 so that it automatically exits the for loop at the end
        if args[:keplerian]
            aerobraking_phase = 3
        end
        if (index_steps_EOM == 1 || Bool(args[:drag_passage])) && (aerobraking_phase == 1 || aerobraking_phase == 3 || aerobraking_phase == 0)
            continue
        elseif (lowercase(args[:type_of_mission]) == "orbits" || lowercase(args[:type_of_mission]) == "time") && args[:keplerian] == true && aerobraking_phase == 2
            continue
        end

        # Definition of eventsecondstep
        if aerobraking_phase == 0
            events = CallbackSet(stop_firing, apoapsisgreaterperiapsis, impact, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller_orbit, time_check, thrust_factor_integrator, save_intermediate_solution)
            t_event_0 = "stop_firing"
            t_event_1 = "apoapsisgreaterperiapsis"
        # elseif aerobraking_phase == 1 && args[:keplerian] == true
        #     events = CallbackSet(eventfirststep_periapsis, apoapsisgreaterperiapsis, impact, periapsispoint, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller_orbit, time_check, thrust_factor_integrator, gravity_gradient)
        #     t_event_0 = "eventfirststep_periapsis"
        #     t_event_1 = "apoapsisgreaterperiapsis"
        elseif aerobraking_phase == 1
            events = CallbackSet(eventfirststep, apoapsisgreaterperiapsis, impact, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller_orbit, time_check, thrust_factor_integrator, save_intermediate_solution)
            t_event_0 = "eventfirststep"
            t_event_1 = "apoapsisgreaterperiapsis"
        elseif aerobraking_phase == 2 && Bool(args[:drag_passage]) && args[:type_of_mission] != "Entry"
            # events = CallbackSet(out_drag_passage, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, impact, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller, solar_panel_controller, time_check, thrust_factor_integrator, gravity_gradient)
            events = CallbackSet(out_drag_passage, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, impact, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller, time_check, thrust_factor_integrator, save_intermediate_solution)
            t_event_0 = "out_drag_passage"
            t_event_1 = "periapsispoint"
        elseif aerobraking_phase == 2 && Bool(args[:drag_passage]) && args[:type_of_mission] == "Entry"
            # events = CallbackSet(out_drag_passage, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, final_entry_altitude_reached, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller, solar_panel_controller, time_check, thrust_factor_integrator, gravity_gradient)
            events = CallbackSet(out_drag_passage, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, final_entry_altitude_reached, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller, time_check, thrust_factor_integrator, save_intermediate_solution)
            t_event_0 = "final_altitude_reached"
            t_event_1 = "out_drag_passage"
        elseif aerobraking_phase == 2 && index_steps_EOM == 1 && args[:body_shape] == "Blunted Cone"
            # events = CallbackSet(out_drag_passage, apoapsispoint, periapsispoint, impact, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller, solar_panel_controller, time_check, thrust_factor_integrator, gravity_gradient)
            events = CallbackSet(out_drag_passage, apoapsispoint, periapsispoint, impact, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller, time_check, thrust_factor_integrator, save_intermediate_solution)
            t_event_0 = "out_drag_passage"
            t_event_1 = "apoapsispoint"
        elseif aerobraking_phase == 2 && index_steps_EOM == 1 && args[:drag_passage] == false
            # events = CallbackSet(apoapsispoint, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, impact, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller, solar_panel_controller, time_check, thrust_factor_integrator, gravity_gradient)
            events = CallbackSet(apoapsispoint, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, impact, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller, time_check, thrust_factor_integrator, save_intermediate_solution)
            t_event_0 = "apoapsispoint"
            t_event_1 = "periapsispoint"
        elseif aerobraking_phase == 2
            # events = CallbackSet(eventsecondstep, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, impact, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller, solar_panel_controller, time_check, thrust_factor_integrator, gravity_gradient)
            events = CallbackSet(eventsecondstep, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis, impact, reaction_wheel_update, quaternion_normalize, quaternion_update, attitude_controller, time_check, thrust_factor_integrator, save_intermediate_solution)
            t_event_0 = "eventsecondstep"
            t_event_1 = "periapsispoint"
        elseif aerobraking_phase == 3 && args[:keplerian] && lowercase(args[:type_of_mission]) == "time" # In cases where the mission type is time and the orbit is keplerian, just integrate for some predetermined number of steps before saving to avoid the CYGNSS true anomaly issues
            events = CallbackSet(impact, attitude_controller_orbit, reaction_wheel_update, quaternion_normalize, time_check, save_steps_check, save_intermediate_solution)
            t_event_0 = ""
            t_event_1 = ""
        elseif aerobraking_phase == 3
            events = CallbackSet(apoapsispoint, periapsispoint, apoapsisgreaterperiapsis, impact, reaction_wheel_update, quaternion_normalize, attitude_controller_orbit, time_check, thrust_factor_integrator, save_intermediate_solution)
            t_event_0 = "apoapsispoint"
            t_event_1 = "periapsispoint"
        end

        if index_phase_aerobraking == 2 || (index_phase_aerobraking == 1 && args[:keplerian] == 1)
            save_pre_index = length(solution.orientation.time) + 1
            simulator = args[:integrator]
        end

        # Definition step size
        if aerobraking_phase == 1 || aerobraking_phase == 3
            step = 5.0
            r_tol = args[:r_tol_orbit] == 0.0 ? args[:r_tol] : args[:r_tol_orbit] #1e-7
            a_tol = args[:a_tol_orbit] == 0.0 ? args[:a_tol] : args[:a_tol_orbit] #1e-9
            simulator = "Julia"
            method = TRBDF2(autodiff=false)
            dt_max = args[:dt_max_orbit] == 0.0 ? args[:dt_max] : args[:dt_max_orbit] #1e-3
            save_ratio = args[:save_rate]
        elseif aerobraking_phase == 0
            step = 1.0
            r_tol = args[:r_tol_orbit] == 0.0 ? args[:r_tol] : args[:r_tol_orbit] #1e-7
            a_tol = args[:a_tol_orbit] == 0.0 ? args[:a_tol] : args[:a_tol_orbit] #1e-9
            simulator = "Julia"
            method = TRBDF2(autodiff=false)
            dt_max = args[:dt_max_orbit] == 0.0 ? args[:dt_max] : args[:dt_max_orbit] #1e-3
            save_ratio = 5
        elseif aerobraking_phase == 2
            if args[:integrator] == "Julia"
                step = 0.1
                r_tol = args[:r_tol_drag] == 0.0 ? args[:r_tol] : args[:r_tol_drag] #1e-7
                a_tol = args[:a_tol_drag] == 0.0 ? args[:a_tol] : args[:a_tol_drag] #1e-9
                dt_max = args[:dt_max_drag] == 0.0 ? args[:dt_max] : args[:dt_max_drag] #1e-3
                if MonteCarlo
                    method = TRBDF2(autodiff=false) # KenCarp4(autodiff=false)
                else
                    method = TRBDF2(autodiff=false) # KenCarp4(autodiff=false)
                end

                save_ratio = args[:save_rate]
            else
                step = 1/args[:trajectory_rate]
                save_ratio = round(args[:trajectory_rate] / args[:save_rate])
            end
        end

        # Definition length simulation
        length_sim = lowercase(args[:type_of_mission]) == "time" ? min(args[:mission_time], 1e8) : 1e8
        i_sim = 0
        time_solution = Float64[]  # Time solution for the simulation

        cnf.continue_simulation = true

        while cnf.continue_simulation
            index_phase_aerobraking = aerobraking_phase
            # if control mode =! 0, redefine sim setting and creates two more phases until reaching EI and out of the AE phase 2: between 120 km alt
            if aerobraking_phase == 2 && (args[:control_mode] != 0 && args[:control_in_loop] == 0 && cnf.drag_state == true && cnf.sensible_loads == true && cnf.ascending_phase == false)
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

                controller.t = controller.guidance_t_eval[controller.count_controller]
            
            #phase 1.75: between EI km alt and 120 kmcontroller.t
            elseif aerobraking_phase == 2 && args[:control_mode] != 0 && args[:control_in_loop] == 0 && cnf.drag_state == true && cnf.ascending_phase == false
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
            elseif aerobraking_phase == 2 && args[:control_mode] != 0 && args[:control_in_loop] == 0 && cnf.drag_state == false && cnf.ascending_phase == false
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
            elseif aerobraking_phase == 2 && args[:control_mode] != 0 && args[:control_in_loop] == 0 && cnf.drag_state == true && cnf.ascending_phase == true
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
            elseif aerobraking_phase == 2 && args[:control_mode] != 0 && args[:control_in_loop] == 0 && cnf.ascending_phase == true
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
                # Counter for events
                cnf.count_eventfirststep = 0
                cnf.eventfirststep_periapsis = 0
                cnf.count_eventsecondstep = 0
                cnf.count_reached_EI = 0
                cnf.count_reached_AE = 0
                cnf.count_out_drag_passage = 0
                cnf.count_in_drag_passage = 0
                cnf.count_in_drag_passage_nt = 0
                cnf.count_apoapsispoint = 0
                cnf.count_periapsispoint = 0
                cnf.count_impact = 0
                cnf.count_apoapsisgreaterperiapsis = 0
                cnf.count_stop_firing = 0
                cnf.count_guidance = 0
                cnf.count_heat_rate_check = 0
                cnf.count_heat_load_check_exit = 0
                cnf.count_final_entry_altitude_reached = 0
                cnf.time_termination = false

                ## Julia Integrator
                # Time initialization
                initial_time, final_time = time_0 / cnf.TU, (time_0 + length_sim) / cnf.TU

                # Parameter Definition
                param = ODEParams(m, cnf, solution, Float64(index_phase_aerobraking), ip, aerobraking_phase, t_prev, date_initial, time_0, initial_state, gram_atmosphere, gram, numberofpassage, Bool(args[:orientation_sim]), args, IntermediateSolution())

                # Energy targeting 
                if args[:targeting_ctrl] == 1 && aerobraking_phase == 2

                    # println("in_cond: ", [(norm(in_cond[1:3])-m.planet.Rp_e)/1e3, norm(in_cond[4:6])])
                    # println("index_phase_aerobraking: ", index_phase_aerobraking)
                    # println("aerobraking_phase: ", aerobraking_phase)
                    # println("config_drag_state: ", cnf.drag_state)

                    # println("m_aero_alpha_a before targeting: ", m.aerodynamics.α)

                    OE_AI = rvtoorbitalelement(SVector{3, Float64}(in_cond[1:3]), SVector{3, Float64}(in_cond[4:6]), in_cond[7], m.planet)

                    cnf.initial_position_closed_form = OE_AI # need this for simulation to enter the control conditional statements in f!

                    energy_f = target_planning(f!, ip, m, args, param, OE_AI, initial_time, final_time, a_tol, r_tol, method, events, in_cond)

                    # println(cnf.targeting)
                    # println("energy_f: ", energy_f)

                    m.aerodynamics.α = deg2rad(args[:α])    
                    cnf.ascending_phase = false
                    cnf.drag_state = true
                    # println("m_aero_alpha_a after targeting: ", m.aerodynamics.α)

                    if cnf.targeting == 1
                        # current_epoch = date_initial + time_0*seconds # Precompute the current epoch
                        # time_real = DateTime(current_epoch) # date_initial + Second(t0)
                        # timereal = ref_sys.clock(Dates.year(time_real), Dates.month(time_real), Dates.day(time_real), Dates.hour(time_real), Dates.minute(time_real), Dates.second(time_real))

                        # # Timing variables
                        # el_time = value(seconds(current_epoch - m.initial_condition.DateTimeIC)) # Elapsed time since the beginning of the simulation
                        # current_time =  value(seconds(current_epoch - m.initial_condition.DateTimeJ2000)) # current time in seconds since J2000
                        # time_real_utc = to_utc(time_real) # Current time in UTC as a DateTime object
                        # cnf.et = utc2et(time_real_utc) # Current time in Ephemeris Time
                        # m.planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_"*uppercase(m.planet.name), cnf.et))*m.planet.J2000_to_pci' # Construct a rotation matrix from J2000 (Planet-fixed frame 0.0 seconds past the J2000 epoch) to planet-fixed frame

                        # Uncomment from here for targeting with shooting
                        v_E = control_solarpanels_targeting_heatload(energy_f, param, OE_AI) # 28.075

                        println("v_E: ", v_E)

                        cnf.lambda_switch_list = []
                        cnf.time_switch_list = []

                        # root finding num int
                        # cnf.t_switch_targeting = control_solarpanels_targeting_num_int(energy_f, param, time_0, in_cond)
                        
                        sol_lam, time_switch = asim_ctrl_rf(ip, m, time_0, OE_AI, args, v_E, 1.0, false, gram_atmosphere)

                        if time_switch[2] == Inf
                            time_switch[2] = 1e50
                        end

                        cnf.ts_targ_1 = time_switch[1]
                        cnf.ts_targ_2 = time_switch[2]

                        # time_switch = control_solarpanels_targeting_closed_form(energy_f, param, OE_AI)

                        # time_switch = control_solarpanels_targeting_closed_form(energy_f, ip, m, OE_AI, args, 0, true, 0, 0)

                        # println("Time Switch: ", time_switch)

                        # cnf.ts_targ_1 = initial_time + time_switch[1]
                        # cnf.ts_targ_2 = initial_time + time_switch[2]

                        println("hf: ", norm(sol_lam[1:3,end]) - m.planet.Rp_e)
                        println("vf: ", norm(sol_lam[4:6,end]))
                        println("γf: ", asin(sol_lam[1:3,end]'*sol_lam[4:6,end]/norm(sol_lam[4:6,end]) / norm(sol_lam[1:3,end])))

                        fin_energy = norm(sol_lam[4:6,end])^2/2 - m.planet.μ/norm(sol_lam[1:3,end])

                        println("Targeting energy: ", energy_f)
                        println("Final Energy: ", fin_energy)

                        push!(cnf.time_list, sol_lam.t...)
                        push!(cnf.lamv_list, sol_lam[7,:]...)

                        # println(cnf.lambda_switch_list)
                        # println(cnf.time_switch_list)
                        # println(cnf.lamv_list)

                        ip.cm = 0

                        # Parameter Definition
                        if solution.simulation.solution_states != 0
                            println("ip: ", ip.cm)
                            param = (m, index_phase_aerobraking, ip, aerobraking_phase, t_prev, date_initial, time_0, args, initial_state, gram_atmosphere, gram, numberofpassage, Bool(args[:orientation_sim]), MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0), args, ip, MVector{3, Float64}(0.0, 0.0, 0.0), zeros(solution.simulation.solution_states))
                        else
                            param = (m, index_phase_aerobraking, ip, aerobraking_phase, t_prev, date_initial, time_0, args, initial_state, gram_atmosphere, gram, numberofpassage, Bool(args[:orientation_sim]), MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0), args, ip, MVector{3, Float64}(0.0, 0.0, 0.0))
                        end
                    end
                end

                println("ip: ", ip.cm)

                # cnf.targeting = 1

                # hf = 160e3 # args[:AE] * 1e3
                # vf = 4195.0 # 4196.4868
                # γf = 0.10979 # deg2rad(5.874)
                # energy_f = -3.2687e6 # -3.265e6 
                # # energy_f = vf^2 / 2 - (m.planet.μ /(hf + m.planet.Rp_e))

                # current_epoch = date_initial # Precompute the current epoch
                # time_real = DateTime(current_epoch) # date_initial + Second(t0)
                # timereal = ref_sys.clock(Dates.year(time_real), Dates.month(time_real), Dates.day(time_real), Dates.hour(time_real), Dates.minute(time_real), Dates.second(time_real))

                # # Timing variables
                # el_time = value(seconds(current_epoch - m.initial_condition.DateTimeIC)) # Elapsed time since the beginning of the simulation
                # current_time =  value(seconds(current_epoch - m.initial_condition.DateTimeJ2000)) # current time in seconds since J2000
                # time_real_utc = to_utc(time_real) # Current time in UTC as a DateTime object
                # cnf.et = utc2et(time_real_utc) # Current time in Ephemeris Time
                # m.planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_"*uppercase(m.planet.name), cnf.et))*m.planet.J2000_to_pci' # Construct a rotation matrix from J2000 (Planet-fixed frame 0.0 seconds past the J2000 epoch) to planet-fixed frame

                # cnf.t_switch_targeting = control_solarpanels_targeting_closed_form(energy_f, param, OE)

                # println("Targeting switch time: ", cnf.t_switch_targeting)

                # cnf.t_switch_targeting = control_solarpanels_targeting_num_int(energy_f, param, time_0, in_cond)

                # println("Targeting switch time: ", cnf.t_switch_targeting)

                # sol_lam = asim_ctrl_targeting_plot(ip, m, 0, OE, args, hf, vf, γf, energy_f, 100, 0.03565, false, gram_atmosphere)

                # v_E = control_solarpanels_targeting_heatload(energy_f, param, OE) # 28.075

                # v_E = 27.892163870200108

                # println("v_E: ", v_E)

                # cnf.lambda_switch_list = []
                # cnf.time_switch_list = []
                
                # sol_lam, time_switch = asim_ctrl_rf(ip, m, 0, OE, args, v_E, 1.0, true, gram_atmosphere)

                # cnf.ts_targ_1 = time_switch[1]
                # cnf.ts_targ_2 = time_switch[2]
                
                # println("hf: ", norm(sol_lam[1:3,end]) - m.planet.Rp_e)
                # println("vf: ", norm(sol_lam[4:6,end]))
                # println("γf: ", asin(sol_lam[1:3,end]'*sol_lam[4:6,end]/norm(sol_lam[4:6,end]) / norm(sol_lam[1:3,end])))

                # println("Targeting energy: ", energy_f)
                # println("Final Energy: ", norm(sol_lam[4:6,end])^2/2 - m.planet.μ/norm(sol_lam[1:3,end]))


                # push!(cnf.time_list, sol_lam.t...)
                # push!(cnf.lamv_list, sol_lam[7,:]...)

                # Run simulation
                # method = TRBDF2(autodiff=false)
                if !cnf.prob_set
                    # cnf.prob = ODEProblem(ODEFunction(f!, jac=f_jac), in_cond, (initial_time, final_time), param)
                    cnf.prob = ODEProblem(f!, in_cond, (initial_time, final_time), param)
                    # prob = complete(modelingtoolkitize(prob))
                    # ModelingToolkit.generate_jacobian(prob)[2]
                    # cnf.prob = prob
                    cnf.prob_set = true
                else
                    cnf.prob = remake(cnf.prob, u0=in_cond, tspan=(initial_time, final_time), p=param)
                end
                # a_tol = 1e-3
                # r_tol = 1e-2
                a_tol_list = a_tol * ones(length(in_cond))
                r_tol_list = r_tol * ones(length(in_cond))
                if args[:orientation_sim]
                    r_tol_quat = args[:r_tol_quaternion] == 0.0 ? args[:r_tol] : args[:r_tol_quaternion] #1e-7
                    a_tol_quat = args[:a_tol_quaternion] == 0.0 ? args[:a_tol] : args[:a_tol_quaternion] #1e-9
                    a_tol_list[next_index:next_index+3] .= a_tol_quat*ones(4)  # Quaternion
                    r_tol_list[next_index:next_index+3] .= r_tol_quat*ones(4)  # Quaternion
                end
                method = Tsit5()
                sol = solve(cnf.prob, method, abstol=a_tol_list, reltol=r_tol_list, callback=events, dtmax=dt_max/cnf.TU)
                cnf.counter_integrator += 1
                in_cond .= sol.u[end]  # Update initial condition for next step
                # if args[:orientation_sim]
                #     in_cond .= [sol[1,end], sol[2,end], sol[3, end], sol[4, end], sol[5, end], sol[6, end], sol[7, end], sol[8, end], 
                #                 sol[9, end], sol[10, end], sol[11, end], sol[12, end], sol[13, end], sol[14, end], sol[15, end]]
                # else
                #     in_cond .= [sol[1,end], sol[2,end], sol[3, end], sol[4, end], sol[5, end], sol[6, end], sol[7, end], sol[8, end]]
                # end
                # Update model parameters
                # m.body.roots[1].q = SVector{4, Float64}(sol[9, end], sol[10, end], sol[11, end], sol[12, end])
                # m.body.roots[1].ω = SVector{3, Float64}(sol[13, end], sol[14, end], sol[15, end])
                # m.body.roots[1].q = normalize(m.body.roots[1].q)  # Quaternion normalization

                # Save results 
                push!(time_solution, (sol.t * cnf.TU)...)
                time_0 = time_solution[end]

                if aerobraking_phase == 2
                    r0 = SVector{3, Float64}(in_cond[1:3])
                    v0 = SVector{3, Float64}(in_cond[4:6])

                    energy_fin = norm(v0)^2/2 - m.planet.μ/norm(r0)

                    a = -m.planet.μ/(2*energy_fin)

                    e = sqrt(1 - (norm(cross(r0, v0))^2)/(m.planet.μ*a))

                    rp = a*(1 - e)

                    println("rp_fin = ", 2*a - rp)
                end

                if aerobraking_phase == 0
                    new_periapsis(m, in_cond[1:3] * cnf.DU, in_cond[4:6] * cnf.DU / cnf.TU, args)
                end
            end

            # Define breaker campaign impact km or apoapsis greater than periapsis
            continue_campaign = event(cnf.count_impact, cnf.count_apoapsisgreaterperiapsis)

            if continue_campaign == false
                cnf.impact = true
                break
            end

            time_ev_0, time_ev_1 = 0, 0

            if t_event_0 == "stop_firing"
                time_ev_0 = cnf.count_stop_firing
            elseif t_event_0 == "eventfirststep"
                time_ev_0 = cnf.count_eventfirststep
            elseif t_event_0 == "eventfirststep_periapsis"
                time_ev_0 = cnf.eventfirststep_periapsis
            elseif t_event_0 == "eventsecondstep"
                time_ev_0 = cnf.count_eventsecondstep
            elseif t_event_0 == "out_drag_passage"
                time_ev_0 = cnf.count_out_drag_passage
            elseif t_event_0 == "in_drag_passage"
                time_ev_0 = cnf.count_in_drag_passage
            elseif t_event_0 == "in_drag_passage_nt"
                time_ev_0 = cnf.count_in_drag_passage_nt
            elseif t_event_0 == "apoapsispoint"
                time_ev_0 = cnf.count_apoapsispoint
            elseif t_event_0 == "periapsispoint"
                time_ev_0 = cnf.count_periapsispoint
            elseif t_event_0 == "apoapsisgreaterperiapsis"
                time_ev_0 = cnf.count_apoapsisgreaterperiapsis
            elseif t_event_0 == "heat_load_check_exit"
                time_ev_0 = cnf.count_heat_load_check_exit
            elseif t_event_0 == "final_altitude_reached"
                time_ev_0 = cnf.count_final_entry_altitude_reached
            end

            if t_event_1 == "stop_firing"
                time_ev_1 = cnf.count_stop_firing
            elseif t_event_1 == "eventfirststep"
                time_ev_1 = cnf.count_eventfirststep
            elseif t_event_1 == "eventfirststep_periapsis"
                time_ev_1 = cnf.eventfirststep_periapsis
            elseif t_event_1 == "eventsecondstep"
                time_ev_1 = cnf.count_eventsecondstep
            elseif t_event_1 == "out_drag_passage"
                time_ev_1 = cnf.count_out_drag_passage
            elseif t_event_1 == "in_drag_passage"
                time_ev_1 = cnf.count_in_drag_passage
            elseif t_event_1 == "in_drag_passage_nt"
                time_ev_1 = cnf.count_in_drag_passage_nt
            elseif t_event_1 == "apoapsispoint"
                time_ev_1 = cnf.count_apoapsispoint
            elseif t_event_1 == "periapsispoint"
                time_ev_1 = cnf.count_periapsispoint
            elseif t_event_1 == "apoapsisgreaterperiapsis"
                time_ev_1 = cnf.count_apoapsisgreaterperiapsis
            elseif t_event_1 == "heat_load_check_exit"
                time_ev_1 = cnf.count_heat_load_check_exit
            elseif t_event_1 == "final_altitude_reached"
                time_ev_1 = cnf.count_final_entry_altitude_reached
            end

            # Breaker conditions
            if simulator == "Julia"
                if time_ev_0 != 0 || (Bool(args[:drag_passage]) && index_phase_aerobraking == 2.25 && time_ev_1 != 0)
                    cnf.continue_simulation = false
                    break
                end

                if lowercase(args[:type_of_mission]) == "time" && (cnf.time_termination == true || time_solution[end] >= args[:mission_time])
                    continue_campaign = false
                    cnf.continue_simulation = false
                    if args[:print_res]
                        println("Setting continue_simulation to false due to mission time condition.")
                    end
                    break
                end

                if args[:keplerian]
                    cnf.continue_simulation = false
                end
            else
                if args[:control_mode] == 0 && stop_simulation == true
                    cnf.continue_simulation = false
                    break
                end
            end
        end

        # Save Results
        time_0 = time_solution[end]
        # time_0 = save_results(time_solution, save_ratio, param)

        if index_phase_aerobraking == 2 || index_phase_aerobraking == 2.5 || (index_phase_aerobraking == 2.25 && Bool(args[:drag_passage])) || (index_phase_aerobraking == 3 && Bool(args[:keplerian]))
            save_post_index = length(solution.orientation.time)
        end

        # Re-Set count index to 0
        cnf.count_phase = 0

        # Define breaker campaign
        if continue_campaign == false
            if save_post_index == 1
                save_post_index = length(solution.orientation.time)
            end

            break
        end
    end

    # Re-set count index to 0
    cnf.count_dori = 0

    # Check tolerance here on true anomaly
    if args[:drag_passage] == false && !args[:keplerian] && (pi - solution.orientation.oe[end][end] > 1e-3) && continue_campaign == true && args[:body_shape] != "Blunted Cone"
        final_conditions_notmet = true
        if lowercase(args[:type_of_mission]) == "time"
            if args[:print_res]
                println("Final time conditions not met, re-running simulation...")
                println("Current time: ", solution.orientation.time[end], " | Mission time: ", args[:mission_time])
            end
            events = CallbackSet(apoapsispoint, time_check, quaternion_normalize, periapsispoint, save_intermediate_solution)
        else
            if args[:print_res]
                println("Final conditions not met, re-running simulation for apoapsis point...")
            end
            events = CallbackSet(apoapsispoint, quaternion_normalize, periapsispoint, save_intermediate_solution)
        end
    else
        final_conditions_notmet = false
    end

    count_temp = 0

    while final_conditions_notmet
        if args[:print_res]
            println("Final conditions not met, re-running simulation...")
        end
        heat_loads = MVector{length(m.body.links), Float64}(zeros(length(m.body.links)))
        for i in 1:length(m.body.links)
            heat_loads[i] = solution.performance.heat_load[i][end]
        end
        if args[:orientation_sim]
            in_cond = MVector{14+length(m.body.links), Float64}([solution.orientation.pos_ii[1][end], solution.orientation.pos_ii[2][end], solution.orientation.pos_ii[3][end], 
                    solution.orientation.vel_ii[1][end], solution.orientation.vel_ii[2][end], solution.orientation.vel_ii[3][end],
                    solution.performance.mass[end], heat_loads..., solution.orientation.quaternion[1][end],
                    solution.orientation.quaternion[2][end], solution.orientation.quaternion[3][end], solution.orientation.quaternion[4][end],
                    solution.orientation.ω[1][end], solution.orientation.ω[2][end], solution.orientation.ω[3][end]])
        else
            in_cond = MVector{7+length(m.body.links), Float64}([solution.orientation.pos_ii[1][end], solution.orientation.pos_ii[2][end], solution.orientation.pos_ii[3][end], 
                    solution.orientation.vel_ii[1][end], solution.orientation.vel_ii[2][end], solution.orientation.vel_ii[3][end],
                    solution.performance.mass[end], heat_loads...])
        end
        # println("In_cond: ", in_cond)
        # non dimensionalization
        in_cond[1:3] ./= cnf.DU
        in_cond[4:6] .*= cnf.TU / cnf.DU
        in_cond[7] /= cnf.MU
        for i in eachindex(m.body.links)
            in_cond[7 + i] *= cnf.TU^2 / cnf.MU  # Heat loads
        end
        next_index = 7 + length(m.body.links) + 1
        # in_cond[8] *= cnf.TU^2 / cnf.MU # * 1e4
        if args[:orientation_sim]
            normalize!(in_cond[next_index:next_index+3])  # Quaternion normalization
            in_cond[next_index+4:next_index+6] .*= cnf.TU  # Angular velocity
        end


        initial_time, final_time = time_0 / cnf.TU, (time_0 + 5000) / cnf.TU
        step = 0.05
        r_tol = 1e-12
        a_tol = 1e-13
        dt_max = args[:dt_max_orbit] == 0.0 ? args[:dt_max]/10.0 : args[:dt_max_orbit]/10.0 #1e-3

        # counter for events
        cnf.count_eventfirststep = 0
        cnf.eventfirststep_periapsis = 0
        cnf.count_eventsecondstep = 0
        cnf.count_reached_EI = 0
        cnf.count_reached_AE = 0
        cnf.count_out_drag_passage = 0
        cnf.count_in_drag_passage = 0
        cnf.count_in_drag_passage_nt = 0
        cnf.count_apoapsispoint = 0
        cnf.count_periapsispoint = 0
        cnf.count_impact = 0
        cnf.count_apoapsisgreaterperiapsis = 0
        cnf.count_stop_firing = 0
        cnf.count_guidance = 0
        cnf.count_heat_rate_check = 0
        cnf.count_heat_load_check_exit = 0

        # Parameter Definition
        if solution.simulation.solution_states != 0
                                param = (m, index_phase_aerobraking, ip, aerobraking_phase, t_prev, date_initial, time_0, args, initial_state, gram_atmosphere, gram, numberofpassage, Bool(args[:orientation_sim]), MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0), args, ip, MVector{3, Float64}(0.0, 0.0, 0.0), zeros(solution.simulation.solution_states))
        else
            param = (m, 
                index_phase_aerobraking, 
                ip, 
                aerobraking_phase, 
                t_prev, 
                date_initial, 
                time_0, 
                args, 
                initial_state, 
                gram_atmosphere, 
                gram, 
                numberofpassage, 
                Bool(args[:orientation_sim]), 
                MVector{3, Float64}(0.0, 0.0, 0.0), 
                MVector{3, Float64}(0.0, 0.0, 0.0), 
                args, 
                ip, 
                MVector{3, Float64}(0.0, 0.0, 0.0))
        end

        # Run simulation
        prob = ODEProblem(f!, in_cond, (initial_time, final_time), param)
        sol = solve(prob, method, abstol=a_tol, reltol=r_tol, callback=events, dtmax=dt_max/cnf.TU)#abstol=a_tol, reltol=r_tol,

        # Update model parameters
        # m.body.roots[1].q = SVector{4, Float64}(sol[9, end], sol[10, end], sol[11, end], sol[12, end])
        # m.body.roots[1].ω = SVector{3, Float64}([sol[13, end], sol[14, end], sol[15, end]] / cnf.TU)
        # m.body.roots[1].q = normalize(m.body.roots[1].q)  # Quaternion normalization

        cnf.counter_integrator += 1
        time_0 = sol.t[end] * cnf.TU
        # time_0 = save_results(sol.t * cnf.TU, args[:save_rate])
        count_temp += 1

        if count_temp > 15
            println("entering")
            break
        end

        if args[:drag_passage] == false && pi - solution.orientation.oe[end][end] > 1e-3 && continue_campaign == true
            final_conditions_notmet = true
            if args[:type_of_mission] == "time"
                events = CallbackSet(apoapsispoint, time_check, quaternion_normalize, periapsispoint, save_intermediate_solution)
            else
                events = CallbackSet(apoapsispoint, quaternion_normalize, periapsispoint, save_intermediate_solution)
            end
        else
            final_conditions_notmet = false
        end

        if lowercase(args[:type_of_mission]) == "time" && (cnf.time_termination == true || sol.t[end] * cnf.TU >= args[:mission_time])
            continue_campaign = false
            cnf.continue_simulation = false
            final_conditions_notmet = false
            if args[:print_res]
                println("Setting continue_simulation to false due to mission time condition.")
            end
            break
        end
    end

    cnf.save_index_heat = length(solution.orientation.time)
    cnf.time_OP = length(solution.orientation.time)

    max_heat_rate = maximum(solution.performance.heat_rate[1][save_pre_index:save_post_index])
    for i in 2:size(solution.performance.heat_rate, 1)
        max_heat_rate = max(max_heat_rate, maximum(solution.performance.heat_rate[i][save_pre_index:save_post_index]))
    end
    if args[:print_res]
        println("Max heat rate is " * string(max_heat_rate) * " W/cm^2")
    end
    append!(cnf.max_heatrate, max_heat_rate)
    cnf.Δv_man = (g_e * m.engines.Isp) * log((get_spacecraft_mass(m.body, m.body.roots[1])) / solution.performance.mass[end])


    if Bool(args[:print_res])
        # Print Actual periapsis altitude and Vacuum periapsis altitude
        if args[:type_of_mission] != "Entry" && !isempty(cnf.altitude_periapsis)
            println("Actual periapsis altitude " * string(cnf.altitude_periapsis[end]) * " km - Vacuum periapsis altitude = " * string((solution.orientation.oe[1][end] * (1 - solution.orientation.oe[2][end]) - m.planet.Rp_e)*1e-3) * " km")

        # Print Ra new (Apoapsis)
            println("Ra new = " * string(solution.orientation.pos_ii_mag[end]/1e3) * " km") # string((solution.orientation.oe[1][end] * (1 + solution.orientation.oe[2][end]))*1e-3) * " km")
        end

        # Print Heat Rate and Heat Load
        max_heat_load = maximum(solution.performance.heat_load[1][save_pre_index:save_post_index])
        for i in 2:size(solution.performance.heat_load, 1)
            max_heat_load = max(max_heat_load, maximum(solution.performance.heat_load[i][save_pre_index:save_post_index]))
        end

        if args[:keplerian] == false
            println("HEAT RATE IS " * string(maximum(cnf.max_heatrate[end])) * " W/cm^2")
            println("HEAT LOAD IS " * string(max_heat_load) * " J/cm^2")
        end

        # Print Fuel Mass
        println("Fuel Mass is " * string(solution.performance.mass[end] - get_spacecraft_mass(m.body, dry=true)) * " kg")

        # Print Total Time
        if args[:keplerian] == false && save_post_index > 1
            println("Total time is " * string(solution.orientation.time[save_post_index-1] - solution.orientation.time[save_pre_index]) * " s")
        else
            println("Total time is " * string(solution.orientation.time[end] - solution.orientation.time[save_pre_index]) * " s")
        end

        # Print Delta-v and Delta-E

        println("Delta-v is " * string(cnf.Δv_man) * " m/s")
        if !isempty(solution.forces.energy)
            ΔE = solution.forces.energy[max(1, save_post_index-1)] - solution.forces.energy[save_pre_index]
        else
            ΔE = 0.0
        end
        println("Delta-E is " * string(ΔE * 1e-3) * " kJ")#Int64(cnf.time_IP)

        # Find periapsis latitude and longitude
        min_index = argmin(solution.orientation.alt[save_pre_index:save_post_index]) + save_pre_index - 1
        println("Latitude of periapsis is " * string(rad2deg(solution.orientation.lat[min_index])) * " deg")
        println("Longitude of periapsis is " * string(rad2deg(solution.orientation.lon[min_index])) * " deg")

        # If the body shape is 'Blunted Cone', print additional information
        if args[:body_shape] == "Blunted Cone"
            # Max Dynamic Pressure
            max_value = maximum(solution.performance.q)
            max_index = argmax(solution.performance.q)
            println("Max Dynamic Pressure " * string(max_value) * " N/m^2 at time " * string(solution.orientation.time[max_index]) * " s")

            # Max Heat Rate
            max_value = maximum(maximum(solution.performance.heat_rate))
            max_index = argmax(maximum(solution.performance.heat_rate))
            println("Max Heat Rate " * string(max_value) * " W/cm^2 at time " * string(solution.orientation.time[max_index]) * " s")

            # Max Heat Load
            max_value = maximum(maximum(solution.performance.heat_load))
            max_index = argmax(maximum(solution.performance.heat_load))
            println("Max Heat Load " * string(max_value) * " J/cm^2 at time " * string(solution.orientation.time[max_index]) * " s")
        end
    end

    append!(cnf.periapsis_list, minimum(solution.orientation.alt[save_pre_index:save_post_index])*1e-3)
    append!(cnf.orbit_number_list, cnf.count_numberofpassage + 1)
    append!(cnf.Δv_list, cnf.Δv_man)
    return continue_campaign
end