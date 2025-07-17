include("../utils/Reference_system.jl")
include("../utils/quaternion_utils.jl")
include("../physical_models/Gravity_models.jl")
include("../physical_models/Density_models.jl")
include("../physical_models/Aerodynamic_models.jl")
include("../physical_models/Thermal_models.jl")
include("../physical_models/Perturbations.jl")

include("../control/Control.jl")
include("../control/Propulsive_maneuvers.jl")

using LinearAlgebra
using OrdinaryDiffEq
using ForwardDiff
using StaticArrays
using SPICE
using Dates
using AstroTime
using PythonCall
using FiniteDiff


import .config
import .ref_sys
import .quaternion_utils

function f_jac(J, u, param, t0)
    q1, q2, q3, q4 = u[9:12]
    ω1, ω2, ω3 = u[13:15]
    m_ref = param[1]
    index_phase_aerobraking = param[2]
    ip = param[3]
    aerobraking_phase = param[4]
    # t_prev = param[5]
    date_initial = param[6]
    time_0 = param[7]
    args = param[8]
    initial_state = param[9]
    gram_atmosphere = param[10]
    gram = param[11]
    numberofpassage = param[12]
    orientation_sim = param[13]
    args = param[16]
    ip = param[17]

    MonteCarlo = Bool(ip.mc)
    wind_m = Bool(ip.wm)
    OE = SVector{7, Float64}([initial_state.a, initial_state.e, initial_state.i, initial_state.Ω, initial_state.ω, initial_state.vi, initial_state.m])

    if (OE[1] > (m_ref.planet.Rp_e*1e-3 + args[:EI])*1e3) && (args[:drag_passage] == false) && (args[:body_shape] == "Spacecraft")
        index_steps_EOM = 3
    else
        index_steps_EOM = 1
    end

    # i = OE[3]
    # Ω = OE[4]
    # ω = OE[5]


    r0, v0 = orbitalelemtorv(OE, m_ref.planet)

    Mass = OE[end]

    # Clock
    date_initial = from_utc(DateTime(m_ref.initial_condition.year, 
                                    m_ref.initial_condition.month,
                                    m_ref.initial_condition.day, 
                                    m_ref.initial_condition.hour, 
                                    m_ref.initial_condition.minute, 
                                    m_ref.initial_condition.second))

    fill!(J[1:3, :], 0.0)
    J[1:3, 4:6] .= I(3)
    fill!(J[7, :], 0.0)
    if orientation_sim
        fill!(J[9:12, 1:8], 0.0)
        J[9:12, 9:end] .= 0.5*Float64[0.0 ω3 -ω2 ω1 q4 -q3 q2;
                            -ω2 0.0 ω1 ω2 q3 q4 -q1;
                            ω2 -ω1 0.0 ω3 -q2 q1 q4;
                            -ω1 -ω2 -ω3 0.0 -q1 -q2 -q3]
    end

    t0 *= config.cnf.TU # Convert time to seconds
    function f_jac_forces!(du::AbstractVector, u::AbstractVector)
        m = deepcopy(m_ref)
        current_epoch = date_initial + (t0-m.initial_condition.el_time)*seconds # Precompute the current epoch
        time_real = DateTime(current_epoch) # date_initial + Second(t0)
        timereal = ref_sys.clock(Dates.year(time_real), Dates.month(time_real), Dates.day(time_real), Dates.hour(time_real), Dates.minute(time_real), Dates.second(time_real))
        el_time = value(seconds(current_epoch - m.initial_condition.DateTimeIC)) # Elapsed time since the beginning of the simulation
 
        pos_ii = SVector{3, Float64}(u[1:3] * config.cnf.DU)   # Inertial position
        vel_ii = SVector{3, Float64}(u[4:6] * config.cnf.DU)   # Inertial velocity
        if orientation_sim
            quaternion = SVector{4, Float64}(@view u[9:12]) # Quaternion
            ω = SVector{3, Float64}((@view u[13:15]) / config.cnf.TU)                # Angular velocity vector [rad / s]
            m.body.roots[1].q .= quaternion
            m.body.roots[1].ω .= ω
        else
            quaternion = SVector{4, Float64}(0.0, 0.0, 0.0, 0.0) # Quaternion, set to all zeros if orientation simulation is not enabled
            ω = SVector{3, Float64}(0.0, 0.0, 0.0)                # Angular velocity vector [rad / s]
        end

        mass = u[7] * config.cnf.MU # Mass of the spacecraft
        pos_ii_mag = norm(pos_ii) # Magnitude of the inertial position vector
        vel_ii_mag = norm(vel_ii) # Magnitude of the inertial velocity vector
        γ = m.planet.γ
        μ_fluid = m.planet.μ_fluid
        bodies, root_index = config.traverse_bodies(m.body, m.body.roots[1])
        m.planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_"*uppercase(m.planet.name), config.cnf.et))*m.planet.J2000_to_pci'
        γ = m.planet.γ
        μ_fluid = m.planet.μ_fluid
        pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, config.cnf.et)
        pos_pp_mag = norm(pos_pp) # Magnitude of the planet-relative position vector
        vel_pp_mag = norm(vel_pp) # Magnitude of the planet-relative velocity vector
        OE = rvtoorbitalelement(pos_ii, vel_ii, mass, m.planet)
        vi = OE[6] # Initial velocity in the inertial frame

        h_ii = cross(pos_ii, vel_ii) # Inertial angular momentum vector
        h_pp = cross(pos_pp, vel_pp) # Planet-relative angular momentum vector
        h_ii_mag = norm(h_ii) # Magnitude of the inertial angular momentum vector
        h_pp_mag = norm(h_pp) # Magnitude of the planet-relative angular momentum vector
        h_pp_hat = h_pp / h_pp_mag # Planet-relative angular momentum unit vector
        h_ii_hat = h_ii / h_ii_mag # Inertial angular momentum unit vector
        γ_ii = acos(clamp(h_ii_mag / (pos_ii_mag * vel_ii_mag), -1.0, 1.0))
        if dot(pos_ii, vel_ii) < 0
            γ_ii = -γ_ii # Adjust the angle if the dot product is negative
        end

        γ_pp = acos(clamp(h_pp_mag / (pos_pp_mag * vel_pp_mag), -1.0, 1.0)) # Planet-relative angle
        if dot(pos_pp, vel_pp) < 0
            γ_pp = -γ_pp # Adjust the angle if the dot product is negative
        end

        alt,lat,lon = rtolatlong(pos_pp, m.planet, args[:topography_model] == "Spherical Harmonics" && norm(pos_ii) - m.planet.Rp_e < args[:EI] * 1e3) # Planet-relative altitude, latitude, and longitude

        uD, uN, uE = latlongtoNED([alt, lat, lon])
        vN = dot(vel_pp, uN)
        vE = dot(vel_pp, uE)
        axi_pp = atan(vE, vN) # Planet-relative azimuth angle

        if ip.dm == 0
            ρ, T_p, wind = density_constant(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 1
            ρ, T_p, wind = density_exp(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 2
            ρ, T_p, wind = density_no(alt, m.planet, lat, lon, timereal, t0, t_prev, MonteCarlo, wind_m, args)
        elseif ip.dm == 3
            ρ, T_p, wind = density_gram(alt, m.planet, lat, lon, MonteCarlo, wind_m, args, el_time, gram_atmosphere, gram)
            ρ, T_p, wind = pyconvert(Float64, ρ), pyconvert(Float64, T_p), SVector{3, Float64}([pyconvert(Float64, wind[1]), pyconvert(Float64, wind[2]), pyconvert(Float64, wind[3])])
        elseif ip.dm == 4
            ρ, T_p, wind = density_nrlmsise(alt, m.planet, lat, lon, MonteCarlo, wind_m, args, time_real)
        end

        p = 0.0
        if args[:body_shape] == "Spacecraft"
            length_car = config.get_spacecraft_length(m.body, m.body.roots[1]) # Length of the spacecraft
        elseif args[:body_shape] == "Blunted Cone"
            length_car = m.body.base_radius * 2.0
        end

        Re = vel_pp_mag * ρ * length_car / μ_fluid  # Reynolds number

        # Mach Number
        sound_velocity = sqrt(γ * m.planet.R * T_p)
        Mach = vel_pp_mag / sound_velocity
        S = sqrt(γ/2.0) * Mach    # Molecular speed ratio
        heat_load = u[8] * config.cnf.MU / config.cnf.TU^2 # * 1e4

        if (index_phase_aerobraking == 2 || index_phase_aerobraking == 1.75 || index_phase_aerobraking == 2.25) && config.cnf.drag_state && length(config.cnf.initial_position_closed_form) != 0
            # evaluates the closed form solution the first time at EI km
            if abs(pos_ii_mag - m.planet.Rp_e - args[:EI] * 1.0e3) <= 1.0e-2 && (args[:control_mode] == 2 || args[:control_mode] == 3) && config.cnf.time_switch_1 == 0
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

            if index_phase_aerobraking == 2
                if Bool(args[:control_in_loop])
                    state_flesh1 = [[T_p, ρ, S]]
                    if ip.cm == 3
                        α = control_solarpanels_openloop(ip, m, args, [1,1], config.cnf.state_flesh1[1], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE, true, gram_atmosphere)
                    elseif ip.cm == 2
                        α = control_solarpanels_heatload(ip, m, args, [1,1], config.cnf.state_flesh1[1], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE, gram_atmosphere)
                        # println("control_solarpanels_heatload: ", config.cnf.α)
                    elseif ip.cm == 1
                        α = control_solarpanels_heatrate(ip, m, args, [1,1], config.cnf.state_flesh1[1], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE)
                    elseif ip.cm == 0
                        α = no_control(ip, m, args, [1,1], config.cnf.state_flesh1[1], t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE)
                    end
                elseif args[:control_in_loop] == false && args[:integrator] == "Julia"
                    if config.controller.count_controller != config.controller.count_prev_controller && config.controller.stored_state == 0 && t0 != config.controller.prev_time
                        # push!(config.cnf.state_flesh1, [T_p, ρ, S]) # might have to change to push!

                        # if config.controller.count_controller == 2
                        #     state = config.cnf.state_flesh1[end]
                        # else
                        #     state = config.cnf.state_flesh1[end-1]
                        #     deleteat!(config.cnf.state_flesh1, 1)
                        # end

                        # config.controller.stored_state = 1
                        # config.controller.prev_time = time_0

                        if ip.cm == 3
                            α = control_solarpanels_openloop(ip, m, args, [1,1], state, t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE, true, gram_atmosphere)
                        elseif ip.cm == 2
                            α = control_solarpanels_heatload(ip, m, args, [1,1], state, t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE, gram_atmosphere)
                        elseif ip.cm == 1
                            α = control_solarpanels_heatrate(ip, m, args, [1,1], state, t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE)
                        elseif ip.cm == 0
                            α = no_control(ip, m, args, [1,1], state, t0 - config.cnf.time_IEI, config.cnf.initial_position_closed_form, OE)
                        end
                    end
                end
            end

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

        wE, wN, wU = wind
        wind_pp = SVector{3, Float64}(wN*uN + wE*uE - wU*uD) # Wind vector in the planet-relative frame
        vel_pp_rw = vel_pp + wind_pp # Wind-relative velocity in the planet-relative frame
        vel_pp_rw_hat = normalize(vel_pp_rw) # Wind-relative velocity unit vector in the planet-relative frame
        q = 0.5 * ρ * norm(vel_pp_rw)^2 # Dynamic pressure in the planet-relative frame

        gravity_ii = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize gravity vector
        # Nominal gravity calculation
        if ip.gm == 0
            gravity_ii += mass * gravity_const(pos_ii_mag, pos_ii, m.planet)
        elseif ip.gm == 1
            gravity_ii += mass * gravity_invsquared(pos_ii_mag, pos_ii, m.planet)
        elseif ip.gm == 2
            gravity_ii += mass * (args[:gravity_harmonics] == 1 ? gravity_invsquared(pos_ii_mag, pos_ii, m.planet) : gravity_invsquared_J2(pos_ii_mag, pos_ii, m.planet))
        elseif ip.gm == 3
            gravity_ii += mass * gravity_GRAM(pos_ii, lat, lon, alt, m.planet, mass, vel_ii, el_time, gram_atmosphere, args, gram)
        end

        if length(args[:n_bodies]) != 0
            for k = 1:length(args[:n_bodies])  
                gravity_ii += mass * gravity_n_bodies(config.cnf.et, pos_ii, m.planet, config.cnf.n_bodies_list[k])
            end
        end
        if args[:gravity_harmonics] == 1
            gravity_ii += mass * m.planet.L_PI' * acc_gravity_pines!(pos_pp, m.planet.Clm, m.planet.Slm, args[:L], args[:M], m.planet.μ, m.planet.Rp_e, m.planet)
        end

        srp_ii = MVector{3, Float64}(zeros(3)) # solar radiation pressure vector
        if orientation_sim
            Rot = [MMatrix{3,3,Float64}(zeros(3, 3)) for i in eachindex(bodies)] # Rotation matrix from the root body to the spacecraft link
            @inbounds for (i, b) in enumerate(bodies)
                Rot[i] .= config.rotate_to_inertial(m.body, b, root_index) # Rotation matrix from the root body to the spacecraft link
            end
        end
        if args[:srp] == true
            # sun_earth_vector = m.planet.J2000_to_pci * SVector{3, Float64}(spkpos("SUN", config.cnf.et, "J2000", "NONE", uppercase(m.planet.name))[1])
            for (i, b) in enumerate(bodies)
                mass_body = b.m # Mass of the spacecraft link
                if b.root
                    mass_body += mass-config.get_spacecraft_mass(m.body, b, dry=true) # Add the mass of the root body
                end
                # Calculate the position of the spacecraft link in inertial frame
                R = Rot[i] # Rotation matrix from the root body to the spacecraft link
                pos_ii_body = pos_ii + R * b.r # Update the position of the spacecraft link in inertial frame
                pos_ii_body_mag = norm(pos_ii_body) # Magnitude of the inertial position of the spacecraft link
                p_srp_unscaled = 4.56e-6  # N / m ^ 2, solar radiation pressure at 1 AU
                srp_ii = mass_body * srp(m.planet, p_srp_unscaled, m.aerodynamics.reflection_coefficient, b.ref_area, b.m, pos_ii_body, config.cnf.et)
                # # Account for angle of incidence of sunlight
                # normal_vector_ii = (rot(root_body.q)'*rot(b.q)'*[1;0;0])
                # normal_vector_ii_hat = normal_vector_ii / norm(normal_vector_ii) # Unit vector of the normal vector
                # # Sun vector
                # sun_vector_ii = sun_earth_vector - pos_ii_body # Vector from the spacecraft link to the Sun
                # sun_vector_ii_hat = sun_vector_ii / norm(sun_vector_ii) # Unit vector of the Sun vector
                # # Calculate the angle of incidence
                # cos_θ = dot(normal_vector_ii_hat, sun_vector_ii_hat) # Cosine of the angle of incidence
                b.net_force += srp_ii # Update the force on the spacecraft link
                if orientation_sim
                    b.net_torque += cross(R*b.r, srp_ii) # Update the torque on the spacecraft link
                end
            end
        end

        lift_pp_hat = normalize(cross(h_pp_hat, vel_pp_rw_hat))
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
                R .= Rot[i] # Rotation matrix from the root body to the spacecraft link
                body_frame_velocity = R' * m.planet.L_PI' * vel_pp_rw # Velocity of the spacecraft link in inertial frame
                
                α_body = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack in radians
                β_body = atan(body_frame_velocity[2], norm([body_frame_velocity[1], body_frame_velocity[3]])) # Sideslip angle in radians
                α[i] = α_body # Angle of attack for the spacecraft link
                β[i] = β_body # Sideslip angle for the spacecraft link
            else
                if b.root
                    # if the body is the root body, then the angle of attack is 90 degrees
                    α[i] = pi/2
                else
                    # body_frame_velocity = rot(b.q) * m.planet.L_PI' * vel_pp_rw # Velocity of the spacecraft link in inertial frame
                    # α[i] = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack for the spacecraft link
                    α[i] = pi/2 # Angle of attack for the spacecraft link, temporary hard code for testing
                end
            end
            if ip.am == 0
                CL, CD = aerodynamic_coefficient_constant(α, m.body, T_p, S, m.aerodynamics, MonteCarlo)
            elseif ip.am == 1
                if orientation_sim
                    CL_body, CD_body, CS_body = aerodynamic_coefficient_fM(α_body, β_body, b, T_p, S, m.aerodynamics, MonteCarlo)
                else
                    CL_body, CD_body = aerodynamic_coefficient_fM(α[i], m.body, T_p, S, m.aerodynamics, MonteCarlo)
                end
            elseif ip.am == 2
                CL, CD = aerodynamic_coefficient_no_ballistic_flight(α, m.body, args, T_p, S, m.aerodynamics, MonteCarlo)
            end

            drag_pp_body = q * CD_body * b.ref_area * drag_pp_hat                       # Planet relative drag force vector
            lift_pp_body = q * CL_body * b.ref_area * lift_pp_hat                       # Planet relative lift force vector
            if orientation_sim
                cross_pp_body = q * CS_body * b.ref_area * cross_pp_hat # Planet relative cross force vector
                cross_ii_body = m.planet.L_PI' * cross_pp_body # Inertial cross force vector
            else
                cross_pp_body = SVector{3, Float64}(0.0, 0.0, 0.0) # Planet relative cross force vector
                cross_ii_body = SVector{3, Float64}(0.0, 0.0, 0.0) # Inertial cross force vector
            end

            drag_ii_body = m.planet.L_PI' * drag_pp_body   # Inertial drag force vector
            lift_ii_body = m.planet.L_PI' * lift_pp_body   # Inertial lift force vector

            # Update the force on the spacecraft link
            b.net_force += drag_ii_body + lift_ii_body + cross_ii_body # Update the force on the spacecraft link
            b.net_torque += cross(R*b.r, drag_ii_body + lift_ii_body + cross_ii_body) # Update the torque on the spacecraft link
            # Update the total CL/CD
            CL += CL_body * b.ref_area
            CD += CD_body * b.ref_area
            total_area += b.ref_area # Update the total area
            drag_ii += drag_ii_body # Update the total drag force
            lift_ii += lift_ii_body # Update the total lift force
            drag_pp += drag_pp_body # Update the total drag force in planet relative frame
            lift_pp += lift_pp_body # Update the total lift force in planet relative frame
        end
        
        # Normalize the aerodynamic coefficients
        CL = CL / total_area
        CD = CD / total_area

        # Check if propellant mass is greater than 0 kg
        if config.cnf.index_propellant_mass == 1
            if mass - config.get_spacecraft_mass(m.body, m.body.roots[1], dry=true) <= 0.5
                config.cnf.index_propellant_mass = 0
                m.engines.T = 0
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
        thrust_pp_hat = normalize(drag_pp_hat * cos(args[:phi]) + cross(D_L_per_pp_hat, drag_pp_hat) * sin(args[:phi]) + D_L_per_pp_hat * dot(D_L_per_pp_hat, drag_pp_hat) * (1 - cos(args[:phi])))
        #these two ways give the same direction
        thrust_pp = thrust_pp_mag * thrust_pp_hat
        thrust_ii = m.planet.L_PI' * thrust_pp

        # Attitude control torques
        τ_rw = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize reaction wheel torque vector
        total_rw_h = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize total reaction wheel angular momentum vector
        rw_h = MVector{m.body.n_reaction_wheels, Float64}(zeros(m.body.n_reaction_wheels)) # Initialize vector of reaction wheel angular momentum magnitudes
        rw_τ = MVector{m.body.n_reaction_wheels, Float64}(zeros(m.body.n_reaction_wheels)) # Initialize vector of reaction wheel torque magnitudes
        if orientation_sim
            counter = 1 # Counter for reaction wheel angular momentum vector
            @inbounds for (i, b) in enumerate(bodies)
                R .= Rot[i] # Rotation matrix from the root body to the spacecraft link
                if b.gyro != 0.0 # If the body has reaction wheels and we are in the aerobraking phase
                    τ = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize reaction wheel torque vector
                    clamp!(b.rw, -b.max_h, b.max_h) # Clamp the reaction wheel angular momentum to the maximum angular momentum
                    @inbounds for j in 1:b.gyro
                        if abs(b.rw[j] - b.max_h) == 0.0 && sign(b.ω_wheel_derivatives[j]) == sign(b.rw[j]) # If the reaction wheel angular momentum is at its maximum and the derivative is in the same direction
                            b.ω_wheel_derivatives[j] = 0.0 # Set the angular momentum derivative to zero if the maximum angular momentum is reached
                        end
                        rw_torque = R*b.J_rw[:, j] * b.ω_wheel_derivatives[j] # Update the reaction wheel torque
                        if norm(rw_torque) > b.max_torque
                            rw_torque = normalize(rw_torque) * b.max_torque # Limit the reaction wheel torque to the maximum torque
                        end
                        τ .+= rw_torque # Sum the reaction wheel torques
                        total_rw_h .+= R * b.J_rw[:, j] * b.rw[j] # Update the total reaction wheel angular momentum
                        rw_h[counter] = b.rw[j] # Update the reaction wheel angular momentum vector
                        rw_τ[counter] = norm(rw_torque) # Update the reaction wheel torque vector
                        counter += 1 # Increment the counter for the reaction wheel angular momentum vector
                    end
                    b.rw_τ .= R'*τ # Save the reaction wheel torque in the body
                    τ_rw .+= b.rw_τ # Sum the reaction wheel torques
                    b.net_torque .+= τ # Update the torque on the spacecraft link
                end
            end
        end

        body_forces = sum([b.net_force for b in bodies]) # Sum of all forces on the spacecraft links
        force_ii = body_forces + gravity_ii + thrust_ii
        # force_ii = aerodynamics_ii + gravity_ii + thrust_ii + srp_ii # Total inertial external force vector on body [N]
        
        # Torques
        if orientation_sim
            τ_ii = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize torque vector
            # Gravity gradient torque
            R .= config.rotate_to_inertial(m.body, m.body.roots[1], root_index) # Rotation matrix from the root body to the spacecraft link
            inertia_tensor = R * config.get_inertia_tensor(m.body, root_index) * R' # Inertia tensor of the body
            τ_ii += 3.0*m.planet.μ * cross(pos_ii, inertia_tensor*pos_ii) / pos_ii_mag^5 # Gravity gradient torque
            τ_ii += sum([b.net_torque for b in bodies]) # Sum of all torques on the spacecraft links
        end
        
        # f = MVector{7, Float64}(zeros(7))
        du[1:3] .= force_ii / mass * (config.cnf.TU^2 / config.cnf.DU)
        du[4] = heat_rate * config.cnf.TU^3 / config.cnf.MU # * 1e-4
        if orientation_sim
            du[5:7] .= (inertia_tensor\(-hat(ω)*(inertia_tensor*ω + total_rw_h)+ τ_ii)) * config.cnf.TU^2  # Angular acceleration
        end
    end
    # output = orientation_sim ? MMatrix{7, length(u), Float64}(zeros(7, length(u))) : MMatrix{4, length(u), Float64}(zeros(4, length(u))) # Initialize output vector
    fd_jac_block = orientation_sim ? MMatrix{7, length(u), Float64}(zeros(7, length(u))) : MMatrix{4, length(u), Float64}(zeros(4, length(u))) # Initialize finite difference Jacobian block
    fd_out_vec = orientation_sim ? MVector{7, Float64}(zeros(7)) : MVector{4, Float64}(zeros(4)) # Initialize finite difference output vector
    FiniteDiff.finite_difference_jacobian!(fd_jac_block, f_jac_forces!, u) # Calculate the Jacobian of the forces
    # Acceleration jacobian
    J[4:6, :] .= fd_jac_block[1:3, :] # Update the Jacobian with the acceleration
    J[8, :] .= fd_jac_block[4, :] # Update the Jacobian with the heat rate
    if orientation_sim
        J[13:15, :] .= fd_jac_block[5:7, :]
    end
end

