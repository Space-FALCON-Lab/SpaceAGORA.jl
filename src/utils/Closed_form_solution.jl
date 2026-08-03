include("../physical_models/Density_models.jl")
include("../physical_models/Aerodynamic_models.jl")
include("../physical_models/MonteCarlo_pertrubations.jl")
include("../utils/Reference_system.jl")
include("Misc.jl")

using LinearAlgebra
using Statistics
using AstroTime

function closed_form(args, mission, initialcondition = 0, T = 0, online = false, α=0, α_profile = [])
    date_initial = from_utc(DateTime(mission.initial_condition.year, mission.initial_condition.month, mission.initial_condition.day, mission.initial_condition.hour, mission.initial_condition.minute, mission.initial_condition.second))

    if args[:body_shape] == "Blunted Cone"
        len_sol = length(config.solution.orientation.time)
        results(zeros(len_sol), zeros(len_sol), zeros(len_sol), zeros(len_sol))

        return zeros(len_sol), zeros(len_sol), zeros(len_sol), zeros(len_sol)
    end

    if online == false
        if args[:type_of_mission] == "Drag Passage"
            step_time = length(config.solution.orientation.time)
            initialcondition = [config.solution.orientation.oe[1][1], config.solution.orientation.oe[2][1], config.solution.orientation.oe[3][1], config.solution.orientation.oe[4][1], config.solution.orientation.oe[5][1], config.solution.orientation.oe[6][1], config.solution.performance.mass[1]]
            T = config.solution.physical_properties.T[1]
            α = config.solution.physical_properties.α[1]
            t0 = config.solution.orientation.time[1]

            t_cf, h_cf, γ_cf, v_cf = closed_form_calculation(args, t0, mission, initialcondition, α, T, date_initial, step_time)
            results(t_cf, h_cf, γ_cf, v_cf)
        elseif args[:type_of_mission] != "Drag Passage"
            # Calculate the number of orbit
            number_orbits = config.solution.orientation.number_of_passage[end]
            length_solution = length(config.solution.orientation.number_of_passage)

            t, h, γ, v = zeros(length_solution), zeros(length_solution), zeros(length_solution), zeros(length_solution)
            cnt = 0
            for i in range(1,ceil(number_orbits))

                # idx_orbit = findall(val -> val == i, config.solution.orientation.number_of_passage)
                # # idx_orbit = [idx for idx, val in enumerate(solution.orientation.numberofpassage) if val == i]
                
                # alt = [(config.solution.orientation.pos_ii_mag[item] - mission.planet.Rp_e) for item in idx_orbit]
                # # alt_index = [idx for idx, val in enumerate(alt) if val <= 160*1e3]

                # # TODO: CHANGE TO ARGS[:EI]
                # alt_index = findall(val -> val <= 160*1e3, alt)

                idx_orbit = findall(x -> x == i, config.solution.orientation.number_of_passage) # [idx for (idx, val) in enumerate(config.solution.orientation.number_of_passage) if val == i]

                alt = [config.solution.orientation.pos_ii_mag[item] - mission.planet.Rp_e for item in idx_orbit]
                alt_index = findall(x -> x < args[:EI]*1e3, alt) # [idx for (idx, val) in enumerate(alt) if val <= 160e3]

                if length(alt_index) == 0
                    len_sol = length(config.solution.orientation.time)
                    results(zeros(len_sol), zeros(len_sol), zeros(len_sol), zeros(len_sol))
                    
                    return zeros(len_sol), zeros(len_sol), zeros(len_sol), zeros(len_sol)
                end

                index = alt_index[1] + idx_orbit[1]
                step_time = length(alt_index)

                initialcondition = [config.solution.orientation.oe[1][index], config.solution.orientation.oe[2][index], config.solution.orientation.oe[3][index], config.solution.orientation.oe[4][index], config.solution.orientation.oe[5][index], config.solution.orientation.oe[6][index], config.solution.performance.mass[index]]

                T = config.solution.physical_properties.T[index]
                α = config.solution.physical_properties.α[index]
                t0 = config.solution.orientation.time[index]

                t_cf, h_cf, γ_cf, v_cf = closed_form_calculation(args, t0, mission, initialcondition, α, T, date_initial, step_time)
                # println("Length of interval: ", length(t[(alt_index[1]+idx_orbit[1]):(alt_index[1]+step_time+idx_orbit[1])]))
                step_time = step_time - 1 # -1 because we start from 0
                t[(alt_index[1]+idx_orbit[1]):(alt_index[1]+step_time+idx_orbit[1])] = t_cf
                h[(alt_index[1]+idx_orbit[1]):(alt_index[1]+step_time+idx_orbit[1])] = h_cf
                γ[(alt_index[1]+idx_orbit[1]):(alt_index[1]+step_time+idx_orbit[1])] = γ_cf
                v[(alt_index[1]+idx_orbit[1]):(alt_index[1]+step_time+idx_orbit[1])] = v_cf
            end
            # For loop for the number of orbits

            results(t, h, γ, v)
        end
    else # online for control
        if args[:montecarlo] == true && Bool(config.cnf.closed_form_solution_off)
            config.cnf.closed_form_solution_off = 0
            state = Dict(:ra => 0, :rp => 0, :i => 0, :Ω => 0, :ω => 0, :vi => 0)
            state[:ra], state[:rp], state[:i], state[:Ω], state[:ω], state[:vi] = initialcondition[1]*(1+initialcondition[2]), initialcondition[1]*(1-initialcondition[2]), initialcondition[3], initialcondition[4], initialcondition[5], initialcondition[6]
            state = monte_carlo_guidance_closedform(state, args)
            initialcondition[1], initialcondition[2], initialcondition[3], initialcondition[4], initialcondition[5], initialcondition[6] = (state[:ra]+state[:rp])/2, (state[:ra]-state[:rp])/(state[:ra]+state[:rp]), state[:i], state[:Ω], state[:ω], state[:vi]
        end

        t_cf, h_cf, γ_cf, v_cf = closed_form_calculation(args, 0, mission, initialcondition, α, T, date_initial, 0, α_profile)
    end

    return t_cf, h_cf, γ_cf, v_cf

    # drag_passage -> save results and initial conditions given
    # all_passage -> save results and initial conditions not given
    # online -> not save results and initial conditions given
end

function closed_form_calculation(args, t0, mission, initialcondition, α, T, date_initial, step_time = 0, α_profile = [], online = 0)
    
    if config.cnf.count_numberofpassage != 1
        t_prev = config.solution.orientation.time[end]
    else
        t_prev = mission.initial_condition.time_rot # value(seconds(date_initial - from_utc(DateTime(2000, 1, 1, 12, 0, 0)))) # mission.initial_condition.time_rot
    end

    pos_ii_org, vel_ii_org = orbitalelemtorv(initialcondition, mission.planet)
    pos_ii = SVector{3, Float64}(pos_ii_org)
    vel_ii = SVector{3, Float64}(vel_ii_org)

    r0 = norm(pos_ii) # Inertial position magnitude
    v0 = norm(vel_ii) # Inertial velocity magnitude
    h0 = r0 - mission.planet.Rp_e

    pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, mission.planet, config.cnf.et)

    LatLong = rtolatlong(pos_pp, mission.planet)
    lat = LatLong[2]
    lon = LatLong[3]
    h0 = LatLong[1]

    h_ii = cross(pos_ii, vel_ii)
    arg = median([-1, 1, norm(h_ii)/(r0*v0)])   # limit to[-1, 1]
    γ0 = acos(arg)

    if dot(pos_ii, vel_ii) < 0
        γ0 = -γ0
    end

    initial_state_angle = initialcondition[6]
    e = initialcondition[2]
    a = initialcondition[1]
    final_state_angle = -initial_state_angle
    E_initialstate = 2 * atan(sqrt((1-e)/(1+e)) * tan(initial_state_angle/2))
    E_finalstate = 2 * atan(sqrt((1-e)/(1+e)) * tan(final_state_angle/2))

    # Evaluate time to reach next state
    Δt = sqrt(a^3 / mission.planet.μ) * ((E_finalstate - e*sin(E_finalstate)) - (E_initialstate - e*sin(E_initialstate)))
    t_p = Δt/2

    mass = initialcondition[end]

    if h0 < args[:EI]*1e3 #if initial condition are lower than drag passage initial condition # this happens only running MC cases
        # let's calculate pos_ii,v_ii for the point of trajectory corresponding to h = 160 km
        h0 = args[:EI]*1e3
        r = mission.planet.Rp_e + h0
        OE = rvtoorbitalelement(pos_ii, vel_ii, mass, mission.planet)
        a, e, i, Ω, ω, vi = OE[1], OE[2], OE[3], OE[4], OE[5], OE[6]
        vi = 2*pi - acos(((a*(1 - e^2)/r)-1)/e)
        E_real_finalstate = 2 * atan(sqrt((1-e)/(1+e)) * tan(-vi/2)) # eccentric anomaly
        Δt = sqrt(a^3 / mission.planet.μ) * ((E_real_finalstate - e*sin(E_real_finalstate)) - (E_initialstate - e*sin(E_initialstate)))
        t_p = Δt/2
    end

    if step_time == 0
        temp = Δt * args[:trajectory_rate]/10

        if temp > length(config.cnf.heat_rate_list)
            step_time = ceil(Int, temp)
        else
            step_time = length(config.cnf.heat_rate_list)
        end
    end
    
    t_cf = collect(range(start=0, stop=Δt, length=step_time))

    # println("t_cf length: ", t_cf)
    # println("step_time: ", step_time)
    # println("Δt: ", Δt)

    cost_3 = v0 * γ0

    h_cf = h0 .+ cost_3*(t_cf - (t_cf.^2/(2*t_p)))

    ρ = density_cheby_norm(h_cf, mission.planet)[1]

    RT = T * mission.planet.R
    S = v0/sqrt(2*RT)
    CL90, CD90 = aerodynamic_coefficient_fM(pi/2, mission.body, T, S, mission.aerodynamics)
    CL0, CD0 = aerodynamic_coefficient_fM(0, mission.body, T, S, mission.aerodynamics)
    Area_tot = mission.body.area_SC + mission.body.area_SA
    
    Rp = mission.planet.Rp_e

    if length(α_profile) == 0
        α_profile = α * ones(length(t_cf))
    else
        if length(α_profile) > length(t_cf)
            α_profile = α_profile[1:length(t_cf)]
        elseif length(α_profile) < length(t_cf)
            last_α = α_profile[end]
            α_profile = α_profile .+ last_α * ones(length(t_cf) - length(α_profile))
        end
    end

    CD_t = CD0 .+ (α_profile * (CD90 - CD0)) / (pi/2)
    CL_t = CL0 .+ (α_profile * (CL90 - CL0)) / (pi/2)

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

    """
    Return the first five Chebyshev polynomials evaluated at z:

        T₀(z), T₁(z), ..., T₄(z)
    """
    function chebyshev_T0_to_T4(z)
        T0 = one(z)
        T1 = z
        T2 = 2z^2 - 1
        T3 = 4z^3 - 3z
        T4 = 8z^4 - 8z^2 + 1

        return (T0, T1, T2, T3, T4)
    end


    """
    Evaluate a fourth-degree univariate Chebyshev expansion:

        a₀T₀(x) + a₁T₁(x) + ... + a₄T₄(x)
    """
    function evaluate_chebyshev4(x, a)
        T = chebyshev_T0_to_T4(x)

        return a[1] * T[1] +
            a[2] * T[2] +
            a[3] * T[3] +
            a[4] * T[4] +
            a[5] * T[5]
    end


    """
    Evaluate the total-degree-four bivariate Chebyshev expansion

        Σ d_{ℓn} T_ℓ(γ̂) T_n(V̂),    ℓ+n ≤ 4.

    The coefficient tuple must use the ordering

        d00, d10, d01,
        d20, d11, d02,
        d30, d21, d12, d03,
        d40, d31, d22, d13, d04.
    """
    function evaluate_chebyshev44(γhat, Vhat, d)
        Tγ = chebyshev_T0_to_T4(γhat)
        TV = chebyshev_T0_to_T4(Vhat)

        d00, d10, d01,
        d20, d11, d02,
        d30, d21, d12, d03,
        d40, d31, d22, d13, d04 = d

        return d00 * Tγ[1] * TV[1] +
            d10 * Tγ[2] * TV[1] +
            d01 * Tγ[1] * TV[2] +

            d20 * Tγ[3] * TV[1] +
            d11 * Tγ[2] * TV[2] +
            d02 * Tγ[1] * TV[3] +

            d30 * Tγ[4] * TV[1] +
            d21 * Tγ[3] * TV[2] +
            d12 * Tγ[2] * TV[3] +
            d03 * Tγ[1] * TV[4] +

            d40 * Tγ[5] * TV[1] +
            d31 * Tγ[4] * TV[2] +
            d22 * Tγ[3] * TV[3] +
            d13 * Tγ[2] * TV[4] +
            d04 * Tγ[1] * TV[5]
    end

    planet_name = lowercase(mission.planet.name)

    if planet_name == "mars"

        # Calibration-domain normalization
        V_center = 4350.0
        ΔV       = 450.0

        γ_center = -5.25       # deg
        Δγ       = 2.25        # deg

        # f₁ Chebyshev coefficients:
        # a₀, a₁, a₂, a₃, a₄
        f1_coeffs = (
            4.560850,
            -1.715925,
            -0.375150,
            0.376025,
            -0.254500
        )

        # f₂ Chebyshev coefficients:
        # d00, d10, d01,
        # d20, d11, d02,
        # d30, d21, d12, d03,
        # d40, d31, d22, d13, d04
        f2_coeffs = (
            -4.74628450,
            -3.69824500,
            -6.06348525,

            1.07688200,
            3.54740275,
            1.27053200,

            -0.26357900,
            -1.12508600,
            -0.91570600,
            -0.17288575,

            0.045299125,
            0.280380250,
            0.224278500,
            -0.004521000,
            -0.025270125
        )

    elseif planet_name == "venus"

        # Calibration-domain normalization
        V_center = 9000.0
        ΔV       = 600.0

        γ_center = -4.375      # deg
        Δγ       = 1.375       # deg

        f1_coeffs = (
            16.3939125,
            -0.6316500,
            -0.3870500,
            0.3974500,
            -0.2209625
        )

        f2_coeffs = (
            -11.578333625,
            -6.131747500,
            -3.504510250,

            1.440936250,
            2.801654000,
            0.536784250,

            0.049912500,
            -0.281389000,
            -0.520214000,
            -0.099717750,

            -0.094991375,
            -0.214584750,
            -0.037281250,
            0.078927750,
            0.016912750
        )

    elseif planet_name == "earth"

        # Calibration-domain normalization
        V_center = 9075.0
        ΔV       = 725.0

        γ_center = -5.125      # deg
        Δγ       = 2.125       # deg

        f1_coeffs = (
            19.807700,
            -3.164075,
            0.592100,
            -0.113625,
            0.043900
        )

        f2_coeffs = (
            -1.989471500,
            0.523520500,
            -1.992902250,

            0.462889750,
            1.965380250,
            0.591760750,

            -0.182282500,
            -0.886610000,
            -0.871433000,
            -0.200818750,

            0.024188625,
            0.286623750,
            0.468237250,
            0.332824000,
            0.043238125
        )

    elseif planet_name == "titan"

        # Calibration-domain normalization
        V_center = 2050.0
        ΔV       = 150.0

        γ_center = -12.5       # deg
        Δγ       = 6.5         # deg

        f1_coeffs = (
            148.480400,
            10.480925,
            0.030500,
            -0.003725,
            0.000600
        )

        f2_coeffs = (
            -8.397098625,
            -4.739801250,
            -1.473748250,

            0.671454000,
            1.024378750,
            0.184619500,

            -0.031104750,
            -0.118313000,
            -0.145654000,
            -0.027513750,

            -0.034697125,
            -0.048001000,
            0.005169000,
            0.019228250,
            0.003789250
        )

    else
        throw(ArgumentError(
            "No f₁/f₂ correction fit is defined for planet: " *
            string(mission.planet.name)
        ))
    end

    # Normalize the dimensional inputs.
    Vhat = (v0 - V_center) / ΔV
    γ0_deg = rad2deg(γ0)
    γhat = (γ0_deg - γ_center) / Δγ

    # Evaluate f₁ directly in m/s².
    f1 = evaluate_chebyshev4(Vhat, f1_coeffs)

    # Evaluate the dimensionless logarithmic f₂ polynomial.
    log_f2 = evaluate_chebyshev44(γhat, Vhat, f2_coeffs)

    # Recover the positive correction and apply the original time scaling.
    f2 = exp(log_f2) * t_cf / (2 * t_p)

    # if mission.planet.name == "mars"
    #     v0_first = 3900
    #     γ0_end = -3.0

    #     mean_b = 2.1858e-4
    #     mean_d = -0.0036
        
    #     a2 = -11.7322
    #     a3 = 1.7060
    #     b2 = 0.2450
    #     b3 = -0.4948

    #     f1 = (-4.894e-11)*v0^4 + (8.678e-7)*v0^3 + (-0.005762)*v0^2 + (16.98)*v0 - 1.871e4

    #     f2 = exp((a2*exp(b2*(rad2deg(γ0) - γ0_end))*exp(mean_b*(v0 - v0_first)) + 
    #               a3*exp(b3*(rad2deg(γ0) - γ0_end))*exp(mean_d*(v0 - v0_first)))) * (t_cf) / (2 * t_p)

    #     # a0 = 0.0016
    #     # c0 = 5e-6
    #     # mean_a = 3.38
    #     # mean_c = 2.6
    #     # mean_b = -8.25
    #     # mean_d = -0.001

    #     # f1 = -0.005 * v0 + 27.87
    #     # f2 = (a0 * (mean_a^(2 * abs(rad2deg(γ0) + 3)) * exp(mean_b * (v0/1000 - 3.7))) + 
    #     #       c0 * (mean_c^(2 * abs(rad2deg(γ0) + 3)) * exp(mean_d * (v0/1000 - 3.7)))) * (t_cf) / (2 * t_p)

    # elseif mission.planet.name == "venus"
    #     # v0_first = 8400
    #     # γ0_end = -3

    #     # mean_b = -0.0123
    #     # mean_d = -0.0039
        
    #     # a2 = 6.5539e-4
    #     # a3 = 0.0026
    #     # b2 = -5.9985
    #     # b3 = -3.4274

    #     f1 = (-1.364e-11)*v0^4 + (4.984e-7)*v0^3 + (-0.006825)*v0^2 + (41.51)*v0 - 9.459e4

    #     # f2 = (a2*exp(b2*(rad2deg(γ0) - γ0_end))*exp(mean_b*(v0 - v0_first)) + 
    #     #       a3*exp(b3*(rad2deg(γ0) - γ0_end))*exp(mean_d*(v0 - v0_first))) * (t_cf) / (2 * t_p)

    #     x = rad2deg(γ0)
    #     y = v0

    #     p00 =   4274 
    #     p10 =  -865.8 
    #     p01 =  -2.216 
    #     p20 =   30.28
    #     p11 =   0.2986  
    #     p02 =   0.0004211 
    #     p30 =   1.309 
    #     p21 =  -0.003775  
    #     p12 =  -3.272e-5  
    #     p03 =  -3.478e-8  
    #     p40 =  -0.2126  
    #     p31 =  -0.0005503 
    #     p22 =  -2.191e-7
    #     p13 =   1.063e-9  
    #     p04 =   1.044e-12

    #     f2 = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3 + p04*y^4

    #     f2 = exp(f2) * (t_cf) / (2 * t_p)
    # elseif mission.planet.name == "earth"
    #     # v0_first = 8350
    #     # γ0_end = -3

    #     # mean_b = -0.0164
    #     # mean_d = -7.9792e-4
        
    #     # a2 = 2.7402e-6
    #     # a3 = 0.0047
    #     # b2 = -5.1505
    #     # b3 = -0.4760

    #     f1 = (1.271e-12)*v0^4 + (-4.733e-8)*v0^3 + (0.0006621)*v0^2 + (-4.127)*v0 + 9697

    #     # f2 = (a2*exp(b2*(rad2deg(γ0) - γ0_end))*exp(mean_b*(v0 - v0_first)) + 
    #     #       a3*exp(b3*(rad2deg(γ0) - γ0_end))*exp(mean_d*(v0 - v0_first))) * (t_cf) / (2 * t_p)

    #     x = rad2deg(γ0)
    #     y = v0

    #     p00 =        4731  
    #     p10 =      -759.6
    #     p01 =      -2.402 
    #     p20 =       47.19  
    #     p11 =      0.2952  
    #     p02 =   0.0004592  
    #     p30 =      -1.377  
    #     p21 =    -0.01233  
    #     p12 =  -3.823e-05  
    #     p03 =  -3.913e-08  
    #     p40 =     0.00949  
    #     p31 =   0.0001648  
    #     p22 =   7.891e-07  
    #     p13 =   1.644e-09 
    #     p04 =   1.252e-12  

    #     f2 = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3 + p04*y^4

    #     f2 = exp(f2) * (t_cf) / (2 * t_p)
    # elseif mission.planet.name == "titan"
    #     # v0_first = 1900
    #     # γ0_end = -6

    #     # mean_b = -0.0188
    #     # mean_d = -0.0049
        
    #     # a2 = 2.3512e-6
    #     # a3 = 1.0144e-4
    #     # b2 = -1.0153
    #     # b3 = -0.5682

    #     f1 = (9.39e-12)*v0^4 + (-8.141e-8)*v0^3 + (0.0002664)*v0^2 + (-0.3915)*v0 + 219.4

    #     # f2 = (a2*exp(b2*(rad2deg(γ0) - γ0_end))*exp(mean_b*(v0 - v0_first)) + 
    #     #       a3*exp(b3*(rad2deg(γ0) - γ0_end))*exp(mean_d*(v0 - v0_first))) * (t_cf) / (2 * t_p)

    #     x = rad2deg(γ0)
    #     y = v0

    #     # println("x: ", x)
    #     # println("y: ", y)

    #     p00 =       934.1
    #     p10 =       -33.4  
    #     p01 =      -1.903  
    #     p20 =      0.4016 
    #     p11 =     0.04816  
    #     p02 =    0.001434  
    #     p30 =    0.001327  
    #     p21 =  -0.0003013  
    #     p12 =  -2.301e-5
    #     p03 =  -4.798e-7  
    #     p40 =  -0.0001555 
    #     p31 =  -4.661e-6
    #     p22 =   2.175e-8
    #     p13 =   3.506e-9
    #     p04 =   5.988e-11

    #     f2 = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3 + p04*y^4

    #     # println("f2: ", exp(f2))

    #     f2 = exp(f2) * (t_cf) / (2 * t_p)
    # end

    f2_solar_panels = f2 .* α * mission.body.area_SA / Area_tot
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

function results(t_cf, h_cf, γ_cf, v_cf)
    # Save results
    append!(config.solution.closed_form.t_cf, t_cf)
    append!(config.solution.closed_form.h_cf, h_cf)
    append!(config.solution.closed_form.γ_cf, γ_cf)
    append!(config.solution.closed_form.v_cf, v_cf)
end