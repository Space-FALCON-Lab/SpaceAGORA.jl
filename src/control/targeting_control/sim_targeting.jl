include("../../physical_models/Density_models.jl")
include("../../physical_models/Aerodynamic_models.jl")
# include("../../physical_models/Gravity_models.jl")
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

function asim_ctrl_targeting(t_switch, param, time_0, in_cond)
    ip = param[3]

    wind_m = false
    if ip.wm == 1
        wind_m = true
    end

    MonteCarlo = false
    if ip.mc == 1
        MonteCarlo = true
    end

    function f_ctrl_rf!(y_dot, in_cond, param, t0)

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
        t_switch = param[12]
        
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
        OE = rvtoorbitalelement(pos_ii, vel_ii, mass, m.planet)
        vi = OE[6]

        # if vi > 0 && vi < pi/2 && config.cnf.ascending_phase == false
        #     config.cnf.ascending_phase = true
        # elseif vi >= pi/2 && vi <= pi && config.cnf.ascending_phase == true && args[:body_shape] == "Blunted Cone"
        #     config.cnf.ascending_phase = false
        # end

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

        # if aerobraking_phase == 2 || aerobraking_phase == 0
        #     if (pos_ii_mag - m.planet.Rp_e - args[:EI] * 1e3) <= 0 && config.cnf.drag_state == false && config.cnf.ascending_phase == false
        #         config.cnf.drag_state = true
        #         config.cnf.time_IEI = t0
        #     elseif (pos_ii_mag - m.planet.Rp_e >= args[:EI] * 1e3) && config.cnf.drag_state == true && config.cnf.ascending_phase
        #         config.cnf.drag_state = false
        #         config.cnf.time_OEI = t0
        #     end
        # end

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

        # Heat rate
        if (index_phase_aerobraking == 2 || index_phase_aerobraking == 1.75 || index_phase_aerobraking == 2.25) && config.cnf.drag_state && length(config.cnf.initial_position_closed_form) != 0
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

        # Convert wind to pp(PCPF) frame
        wE = wind[1] # positive to the east , m / s
        wN = wind[2] # positive to the north , m / s
        wU = wind[3] # positive up , m / s

        wind_pp = wN * uN + wE * uE - wU * uD         # wind velocity in pp frame, m / s 
        vel_pp_rw = vel_pp + wind_pp                  # relative wind vector, m / s
        vel_pp_rw_hat = vel_pp_rw / norm(vel_pp_rw)   # relative wind unit vector 

        # Dynamic Pressure, CHANGE THE VELOCITY WITH THE WIND VELOCITY
        q = 0.5 * ρ * norm(vel_pp_rw)^2               # dynamic pressure based on wind, Pa

        if args[:struct_ctrl] == 1
            α_struct = control_struct_load(ip, m, args, S, T_p, q, MonteCarlo)

            config.cnf.α = min(config.cnf.α, α_struct) # limit the angle of attack to the structural load control
        end

        if t0 >= t_switch
            config.cnf.α = 0
        else
            state = [T_p, ρ, S]
            index_ratio = [1,1]
            config.cnf.α = control_solarpanels_heatrate(ip, m, args, index_ratio, state, t0 - config.cnf.t_switch_targeting, config.cnf.initial_position_closed_form, OE)
        end

        α = config.cnf.α

        # Heat Rate 
        if (index_phase_aerobraking == 2 || index_phase_aerobraking == 1.75 || index_phase_aerobraking == 2.25) && config.cnf.drag_state && length(config.cnf.initial_position_closed_form) != 0
            
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

        return y_dot
    end

    ## EVENTS
    function out_drag_passage_condition(y, t, integrator)
        m = integrator.p[1]
        args = integrator.p[8]

        if abs(norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - args[:AE]*1e3) <= 1e-5 
            if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 2
                config.cnf.α = m.aerodynamics.α
            elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 3
                config.cnf.α = 0.0
            end
        end

        norm(y[1:3]) * config.cnf.DU - m.planet.Rp_e - args[:AE]*1e3  # upcrossing
    end
    function out_drag_passage_affect!(integrator)
        terminate!(integrator)
    end
    out_drag_passage = ContinuousCallback(out_drag_passage_condition, out_drag_passage_affect!, nothing)

    # # Time initialization
    initial_time, final_time = time_0 / config.cnf.TU, (time_0 + 1e8) / config.cnf.TU

    param = (param..., t_switch)

    method = Tsit5()
    a_tol = 1e-9
    r_tol = 1e-9

    events = out_drag_passage

    # Run simulation
    prob = ODEProblem(f_ctrl_rf!, in_cond, (initial_time, final_time), param)
    sol = solve(prob, method, abstol=a_tol, reltol=r_tol, callback=events)

    return sol
end

function closed_form_targeting(t0, mission, initialcondition, T, t_cf, t_p, mass, α_profile)
    v0, γ0, h0 = initialcondition[1], initialcondition[2], initialcondition[3]

    # t_cf = collect(range(start=0, stop=Δt, length=step_time))

    cost_3 = v0 * γ0

    h_cf = h0 .+ cost_3*(t_cf - (t_cf.^2/(2*t_p)))

    ρ = density_polyfit(h_cf, mission.planet)[1]

    RT = T * mission.planet.R
    S = v0/sqrt(2*RT)
    CL90, CD90 = aerodynamic_coefficient_fM(pi/2, mission.body, T, S, mission.aerodynamics)
    CL0, CD0 = aerodynamic_coefficient_fM(0, mission.body, T, S, mission.aerodynamics)
    Area_tot = mission.body.area_SC + mission.body.area_SA
    
    Rp = mission.planet.Rp_e

    # if length(α_profile) == 0
    #     α_profile = α * ones(length(t_cf))
    # else
    #     if length(α_profile) > length(t_cf)
    #         α_profile = α_profile[1:length(t_cf)]
    #     elseif length(α_profile) < length(t_cf)
    #         last_α = α_profile[end]
    #         α_profile = α_profile .+ last_α * ones(length(t_cf) - length(α_profile))
    #     end
    # end

    CD_t = CD0 .+ (α_profile * (CD90 - CD0)) / (pi/2)
    CL_t = CL0 .+ (α_profile * (CL90 - CL0)) / (pi/2)

    # CL_t = zeros(length(t_cf))
    # CD_t = zeros(length(t_cf))

    # for i = 1:length(t_cf)
    #     CL_t[i], CD_t[i] = aerodynamic_coefficient_fM(pi/2, mission.body, T, S, mission.aerodynamics)
    # end

    # println(α_profile)
    # println(CD_t)
    # println(" ")

    cost_1 = ρ .* CD_t * Area_tot / (2*mass) .* α_profile
    cost_2 = ρ .* CL_t * Area_tot / (2*mass)

    # a0 = 0.0016
    # c0 = 5e-6
    # mean_a = 3.38
    # mean_c = 2.6
    # mean_b = -8.25
    # mean_d = -0.001

    # f1 = -0.005 * v0 + 27.87
    # f2 = (a0 * (mean_a^(2 * abs(rad2deg(γ0) + 3)) * exp(mean_b * (v0/1000 - 3.7))) + 
    #       c0 * (mean_c^(2 * abs(rad2deg(γ0) + 3)) * exp(mean_d * (v0/1000 - 3.7)))) * (t_cf) / (2 * t_p)

    if mission.planet.name == "mars"
        v0_first = 3900
        γ0_end = -3.0

        mean_b = 2.1858e-4
        mean_d = -0.0036
        
        a2 = -11.7322
        a3 = 1.7060
        b2 = 0.2450
        b3 = -0.4948

        f1 = (-4.894e-11)*v0^4 + (8.678e-7)*v0^3 + (-0.005762)*v0^2 + (16.98)*v0 - 1.871e4

        f2 = exp((a2*exp(b2*(rad2deg(γ0) - γ0_end))*exp(mean_b*(v0 - v0_first)) + 
                  a3*exp(b3*(rad2deg(γ0) - γ0_end))*exp(mean_d*(v0 - v0_first)))) * (t_cf) / (2 * t_p)

        # a0 = 0.0016
        # c0 = 5e-6
        # mean_a = 3.38
        # mean_c = 2.6
        # mean_b = -8.25
        # mean_d = -0.001

        # f1 = -0.005 * v0 + 27.87
        # f2 = (a0 * (mean_a^(2 * abs(rad2deg(γ0) + 3)) * exp(mean_b * (v0/1000 - 3.7))) + 
        #       c0 * (mean_c^(2 * abs(rad2deg(γ0) + 3)) * exp(mean_d * (v0/1000 - 3.7)))) * (t_cf) / (2 * t_p)

    elseif mission.planet.name == "venus"
        # v0_first = 8400
        # γ0_end = -3

        # mean_b = -0.0123
        # mean_d = -0.0039
        
        # a2 = 6.5539e-4
        # a3 = 0.0026
        # b2 = -5.9985
        # b3 = -3.4274

        f1 = (-1.364e-11)*v0^4 + (4.984e-7)*v0^3 + (-0.006825)*v0^2 + (41.51)*v0 - 9.459e4

        # f2 = (a2*exp(b2*(rad2deg(γ0) - γ0_end))*exp(mean_b*(v0 - v0_first)) + 
        #       a3*exp(b3*(rad2deg(γ0) - γ0_end))*exp(mean_d*(v0 - v0_first))) * (t_cf) / (2 * t_p)

        x = rad2deg(γ0)
        y = v0

        p00 =   4274 
        p10 =  -865.8 
        p01 =  -2.216 
        p20 =   30.28
        p11 =   0.2986  
        p02 =   0.0004211 
        p30 =   1.309 
        p21 =  -0.003775  
        p12 =  -3.272e-5  
        p03 =  -3.478e-8  
        p40 =  -0.2126  
        p31 =  -0.0005503 
        p22 =  -2.191e-7
        p13 =   1.063e-9  
        p04 =   1.044e-12

        f2 = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3 + p04*y^4

        f2 = exp(f2) * (t_cf) / (2 * t_p)
    elseif mission.planet.name == "earth"
        # v0_first = 8350
        # γ0_end = -3

        # mean_b = -0.0164
        # mean_d = -7.9792e-4
        
        # a2 = 2.7402e-6
        # a3 = 0.0047
        # b2 = -5.1505
        # b3 = -0.4760

        f1 = (1.271e-12)*v0^4 + (-4.733e-8)*v0^3 + (0.0006621)*v0^2 + (-4.127)*v0 + 9697

        # f2 = (a2*exp(b2*(rad2deg(γ0) - γ0_end))*exp(mean_b*(v0 - v0_first)) + 
        #       a3*exp(b3*(rad2deg(γ0) - γ0_end))*exp(mean_d*(v0 - v0_first))) * (t_cf) / (2 * t_p)

        x = rad2deg(γ0)
        y = v0

        p00 =        4731  
        p10 =      -759.6
        p01 =      -2.402 
        p20 =       47.19  
        p11 =      0.2952  
        p02 =   0.0004592  
        p30 =      -1.377  
        p21 =    -0.01233  
        p12 =  -3.823e-05  
        p03 =  -3.913e-08  
        p40 =     0.00949  
        p31 =   0.0001648  
        p22 =   7.891e-07  
        p13 =   1.644e-09 
        p04 =   1.252e-12  

        f2 = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3 + p04*y^4

        f2 = exp(f2) * (t_cf) / (2 * t_p)
    elseif mission.planet.name == "titan"
        # v0_first = 1900
        # γ0_end = -6

        # mean_b = -0.0188
        # mean_d = -0.0049
        
        # a2 = 2.3512e-6
        # a3 = 1.0144e-4
        # b2 = -1.0153
        # b3 = -0.5682

        f1 = (9.39e-12)*v0^4 + (-8.141e-8)*v0^3 + (0.0002664)*v0^2 + (-0.3915)*v0 + 219.4

        # f2 = (a2*exp(b2*(rad2deg(γ0) - γ0_end))*exp(mean_b*(v0 - v0_first)) + 
        #       a3*exp(b3*(rad2deg(γ0) - γ0_end))*exp(mean_d*(v0 - v0_first))) * (t_cf) / (2 * t_p)

        x = rad2deg(γ0)
        y = v0

        # println("x: ", x)
        # println("y: ", y)

        p00 =       934.1
        p10 =       -33.4  
        p01 =      -1.903  
        p20 =      0.4016 
        p11 =     0.04816  
        p02 =    0.001434  
        p30 =    0.001327  
        p21 =  -0.0003013  
        p12 =  -2.301e-5
        p03 =  -4.798e-7  
        p40 =  -0.0001555 
        p31 =  -4.661e-6
        p22 =   2.175e-8
        p13 =   3.506e-9
        p04 =   5.988e-11

        f2 = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3 + p04*y^4

        # println("f2: ", exp(f2))

        f2 = exp(f2) * (t_cf) / (2 * t_p)
    end

    f2_solar_panels = f2 .* α_profile * mission.body.area_SA / Area_tot
    f2_spacecraft = f2 * pi / 2 * mission.body.area_SC / Area_tot

    ϵ = f1 .+ f2_solar_panels .+ f2_spacecraft

    k1 = (cost_2 .+ (1 ./ (Rp .+ h_cf)))
    k2 = (cost_1 * cost_3) .* (1 .- t_cf/t_p)
    k3 = -mission.planet.g_ref .- ϵ

    cost = v0 - (k2[1]/k1[1] - sqrt((k2[1]/k1[1])^2 - 4 * (k3[1]/k1[1]))) / 2

    v_cf = ((k2 ./ k1) .- sqrt.((k2 ./ k1).^2 - 4*(k3 ./ k1))) / 2 .+ cost
    γ_cf = cost_3 * (1 .- t_cf./t_p) ./ v_cf
    t_cf = [item + t0 for item in t_cf]

    return t_cf, h_cf, γ_cf, v_cf
end