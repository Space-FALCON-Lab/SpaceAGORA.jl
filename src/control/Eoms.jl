include("../physical_models/Density_models.jl")
include("../physical_models/Aerodynamic_models.jl")
include("../physical_models/Gravity_models.jl")
include("../physical_models/Thermal_models.jl")

include("../utils/Reference_system.jl")
include("../utils/Closed_form_solution.jl")

# include("Control.jl")
# include("heatload_control/Utils_timeswitch.jl")

using LinearAlgebra
using OrdinaryDiffEq
using Dates
using AstroTime
using SPICE
using PythonCall
sys = pyimport("sys")


function asim_ctrl(ip, m, time_0, OE, args, k_cf, heat_rate_control, time_switch_eval=false, gram_atmosphere=nothing, time_switch_2=0, reevaluation_mode=1)
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

    CL_90, CD_90 = aerodynamic_coefficient_fM(pi/2, m.body, T, S, m.aerodynamics, 0)
    CL_0, CD_0 = aerodynamic_coefficient_fM(0, m.body, T, S, m.aerodynamics, 0)
    CD_slope = (CD_90 - CD_0) / (pi/2)

    function f_ctrl!(y_dot, in_cond, param, t0)
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

        if time_switch_eval == true

            lambda_switch = (k_cf * 2.0 * m.body.mass * vel_ii_mag) ./ (area_tot * CD_slope * pi)

            if args[:heat_load_sol] == 0
                if lambdav_ii < lambda_switch
                    aoa = 0.0001
                else
                    aoa = m.aerodynamics.α
                end
            elseif args[:heat_load_sol] == 1
                if lambdav_ii < lambda_switch
                    aoa = m.aerodynamics.α
                else
                    aoa = 0.0001
                end
            end
        else
            if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 3
                if t0 >= config.cnf.time_switch_1 && t0 <= time_switch_2
                    aoa = 0.0001
                else
                    aoa = m.aerodynamics.α
                end
            elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 2
                if t0 >= config.cnf.time_switch_1 && t0 <= time_switch_2
                    aoa = m.aerodynamics.α
                else
                    aoa = 0.0001
                end
            end
        end

        # Heat Rate
        heat_rate = heatrate_convective_maxwellian(S, T_p, m, ρ, vel_pp_mag, aoa)

        # Add the control for the heat rate if flash == 3
        if heat_rate_control == true && heat_rate > args[:max_heat_rate]
            state = [T_p, ρ, S]
            index_ratio = [1]
            aoa = control_solarpanels_heatrate(ip, m, args, index_ratio, state)
            heat_rate = args[:max_heat_rate]
        end

        # Convert wind to pp(PCPF) frame
        wE = wind[1] # positive to the east , m / s
        wN = wind[2] # positive to the north , m / s
        wU = wind[3] # positive up , m / s

        wind_pp = wN * uN + wE * uE - wU * uD        # wind velocity in pp frame , m / s
        vel_pp_rw = vel_pp + wind_pp                 # relative wind vector , m / s
        vel_pp_rw_hat = vel_pp_rw / norm(vel_pp_rw)  # relative wind unit vector , nd

        # Dynamic pressure, CHANGE THE VELOCITY WITH THE WIND VELOCITY
        q = 0.5 * ρ * norm(vel_pp_rw)^2            # base on wind - relative velocity

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
        lambdav_dot = -3 * k_cf * ρ * vel_ii_mag^2 * aoa / pi + 
                      lambdav_ii * (ρ* area_tot * CD * vel_ii_mag) / mass - 
                      lambdagamma_ii * ((ρ * area_tot * CL) / (2 * mass) + g_ii /  vel_ii_mag^2 + 1 / (pos_ii_mag)) - 
                      lambdah_ii * γ_ii
        
        lambdag_dot = lambdav_ii * g_ii - lambdah_ii * vel_ii_mag
        
        lambdah_dot = k_cf * ρ * vel_ii_mag^3 * aoa/ (pi * m.planet.H) - 
                      lambdav_ii * ((ρ * area_tot *CD * vel_ii_mag^2) / (2 * mass * m.planet.H) + 2 * g_ii * γ_ii/ (pos_ii_mag)) + 
                      lambdagamma_ii * (ρ * area_tot * CL * vel_ii_mag / (2 * mass * m.planet.H) - 2 * g_ii / ((pos_ii_mag) * vel_ii_mag) + vel_ii_mag / (pos_ii_mag)^2)

        y_dot[1:3] = vel_ii
        y_dot[4:6] = force_ii / mass
        y_dot[7] = lambdav_dot
        y_dot[8] = lambdag_dot
        y_dot[9] = lambdah_dot
        y_dot[10] = heat_rate

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

    function time_switch_func_condition(y, t, integrator)
        m = integrator.p[1]
        k_cf = integrator.p[12]

        # CL_90, CD_90 = aerodynamic_coefficient_fM(pi/2, m.body, T, S, m.aerodynamics, 0)
        # CL_0, CD_0 = aerodynamic_coefficient_fM(0, m.body, T, S, m.aerodynamics, 0)
        # CD_slope = (CD_90 - CD_0) / (pi/2)

        vel_ii = y[4:6]
        vel_ii_mag = norm(vel_ii)

        lambda_switch = (k_cf * 2 * m.body.mass * vel_ii_mag) ./ (m.body.area_tot * CD_slope * pi)
        lambda_switch - y[7]
    end
    function time_switch_func_affect!(integrator)
        # println("entered time_switch_func_affect! in Eoms.jl")
        append!(config.cnf.t_time_switch_func, integrator.t)
        nothing
    end
    time_switch_func = ContinuousCallback(time_switch_func_condition, time_switch_func_affect!)

    if time_switch_eval == true
        # SOLVE EQUATIONS OF MOTIONS - 1 steps
        # USE CLOSED FORM SOLUTION TO DEFINE lambda_zero:
        T = m.planet.T  # fixed temperature
        t_cf, h_cf, γ_cf, v_cf =closed_form(args, m, OE, T, true, m.aerodynamics.α)  # define closed-form solution

        lambdav = v_cf[end]
        lambdag = 0.0
        lambdah = m.planet.μ / (m.planet.Rp_e + h_cf[end])^2

        lambda_v_fin = 10000
        lambda_γ_fin = 10000
        lambda_h_fin = 10000

        lambda_v_fin_actual = 0
        lambda_γ_fin_actual = 0
        lambda_h_fin_actual = 0

        count = 0

        index_phase_aerobraking = nothing
        aerobraking_phase = nothing
        initial_state = nothing

        step = 5

        while abs(lambda_v_fin_actual - lambda_v_fin) > 0.1 || abs(lambda_γ_fin_actual - lambda_γ_fin) > 0.1 || abs(lambda_h_fin_actual - lambda_h_fin) > 0.01
            count += 1

            in_cond = [r0[1], r0[2], r0[3], v0[1], v0[2], v0[3], lambdav, lambdag, lambdah, 0.0]

            # println("in_cond 1: ", in_cond)

            # Time initialization
            initial_time, final_time = time_0, time_0 + 1500

            # Parameter Definition
            param = (m, index_phase_aerobraking, ip, aerobraking_phase, t_prev, date_initial, time_0, args, initial_state, gram_atmosphere, gram, k_cf)

            method = Tsit5()
            a_tol = 1e-7
            r_tol = 1e-7

            events = out_drag_pass

            # Run simulation
            prob = ODEProblem(f_ctrl!, in_cond, (initial_time, final_time), param)
            sol = solve(prob, method, abstol=a_tol, reltol=r_tol, dtmax=step, callback=events)

            r_fin = norm(sol[1:3,end])
            v_fin = norm(sol[4:6,end])
            Q_fin = sol[end,end]
            lambda_v_fin = v_fin
            lambda_γ_fin = 0.0
            lambda_h_fin = m.planet.μ / r_fin^2
            lambda_v_fin_actual = sol[7,end]
            lambda_γ_fin_actual = sol[8,end]
            lambda_h_fin_actual = sol[9,end]

            # println("lambda 2: ", [lambda_v_fin_actual, lambda_γ_fin_actual, lambda_h_fin_actual])

            in_cond = [sol[1,end], sol[2,end], sol[3,end], sol[4,end], sol[5,end], sol[6,end], lambda_v_fin, lambda_γ_fin, lambda_h_fin, Q_fin]

            # println("in_cond 2: ", in_cond)
            
            if (abs(lambda_v_fin_actual - lambda_v_fin) < 0.1 && abs(lambda_γ_fin_actual - lambda_γ_fin) < 0.1 && abs(lambda_h_fin_actual - lambda_h_fin) < 0.01) || count > 4
                break
            end

            # println("time: ", config.cnf.t_out_drag_passage)

            prob = ODEProblem(f_ctrl!, in_cond, (sol.t[end], -10), param)
            sol = solve(prob, method, abstol=a_tol, reltol=r_tol, dtmax=step, callback=events)

            # println("time rev: ", sol.t[end])

            lambdav = sol[7,end]
            lambdag = sol[8,end]
            lambdah = sol[9,end]

            # println("lambda 3: ", [lambdav, lambdag, lambdah])
        end

        # Rerun the simulation with smaller step-size and the right lambda zero
        # Initial condition initialization
        in_cond = [r0[1], r0[2], r0[3], v0[1], v0[2], v0[3], lambdav, lambdag, lambdah, 0.0]

        # println("in_cond: ", in_cond)

        # Time initialization
        initial_time, final_time = time_0, time_0 + 1500

        # Parameter Definition
        param = (m, index_phase_aerobraking, ip, aerobraking_phase, t_prev, date_initial, time_0, args, initial_state, gram_atmosphere, gram, k_cf)

        method = Tsit5()
        a_tol = 1e-9
        r_tol = 1e-9

        events = CallbackSet(out_drag_pass, time_switch_func)

        # Run simulation
        prob = ODEProblem(f_ctrl!, in_cond, (initial_time, final_time), param)
        sol = solve(prob, method, abstol=a_tol, reltol=r_tol, dtmax=step, callback=events)

        temp = config.cnf.t_time_switch_func

        ## Time switch definition
        time_switch = [0.0, 0.0]

        # println("length of temp: ", temp)

        if length(temp) == 2
            # time_switch = temp
            time_switch = temp
            # time_switch[2] = temp[end]
        elseif length(temp) == 1
            time_switch[1] = temp[1]
            time_switch[2] = sol.t[end]
        end

        config.cnf.t_time_switch_func = []

    else  # second time evaluation
        temp_0 = 0
        tp = 1000

        index_phase_aerobraking = nothing
        aerobraking_phase = nothing
        initial_state = nothing

        # Initial Condition Initialization
        in_cond = [r0[1], r0[2], r0[3], v0[1], v0[2], v0[3], 0.0, 0.0, 0.0, config.cnf.heat_load_past]

        step = 1

        if reevaluation_mode == 1  # bigger step for the first times revaluation is performed
            # Time initialization
            initial_time, final_time = time_0, time_0 + 1000

            # Parameter Definition
            param = (m, index_phase_aerobraking, ip, aerobraking_phase, t_prev, date_initial, time_0, args, initial_state, gram_atmosphere, gram, k_cf)

            method = Tsit5()
            a_tol = 1e-9
            r_tol = 1e-9

            events = out_drag_pass

            # Run simulation
            prob = ODEProblem(f_ctrl!, in_cond, (initial_time, final_time), param)
            sol = solve(prob, method, dtmax=step, callback=events)
        elseif reevaluation_mode == 2  # stricter conditions for the last times revaluation is performed
            # Time initialization
            initial_time, final_time = time_0, time_0 + 1000

            # Parameter Definition
            param = (m, index_phase_aerobraking, ip, aerobraking_phase, t_prev, date_initial, time_0, args, initial_state, gram_atmosphere, gram, k_cf)

            method = Tsit5()
            a_tol = 1e-9
            r_tol = 1e-9

            events = out_drag_pass

            # Run simulation
            prob = ODEProblem(f_ctrl!, in_cond, (initial_time, final_time), param)
            sol = solve(prob, method, dtmax=step, callback=events)
        end

        return sol
    end

    # println("time_switch: ", time_switch)

    return sol, time_switch
end