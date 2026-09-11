"""
This module provides high-level functions to run common simulation scenarios.
"""

using StaticArrays

# ---------- 1. Δv stop: callback that integrates Δv from non-grav accel ----------
# 1.1. Discrete event
"""
    Create a discrete callback to stop integration when a Δv target is reached.

    Inputs:
        p: parameter dictionary (must include :masses and any keys required by laser_forces())
        sat: index of satellite to monitor (1 to N)
        target: target Δv value (m/s)
        mode: :magnitude | :component | :direction
        axis: for :component mode, choose :R | :T | :N
        dirECI: for :direction mode, specify direction in ECI frame (3-vector)

    Returns:
        DiscreteCallback object for use with DifferentialEquations.jl solvers
            that stops integration when the Δv target is reached.
"""
function make_dv_stop_callback(p; sat::Int=1, target::Float64=15.0,
                            mode::Symbol=:magnitude, axis::Symbol=:T,
                            dirECI::SVector{3,Float64}=SVector(1.0,0.0,0.0)) # this function creates a discrete callback that stops the integration when a certain Δv target is reached
    masses = p[:masses]
    dv_accum = zeros(3)
    last_t   = Ref{Float64}(NaN)

    idx6(i,off) = 6*(i-1) + off
    function rtn_basis_from_state(u, j)
        r = @SVector [u[idx6(j,1)], u[idx6(j,2)], u[idx6(j,3)]] # r is the position vector of satellite j centered at the center of the Earth
        v = @SVector [u[idx6(j,4)], u[idx6(j,5)], u[idx6(j,6)]]
        er = r/(norm(r)+1e-12); en = cross(r,v); en = en/(norm(en)+1e-12); et = cross(en,er)
        return er,et,en # return RTN unit vectors
    end 

    condition(u,t,integrator) = begin # this function checks if the condition for stopping the integration is metric
        # condition() is called at every step; return true to trigger.
        # Writing condition() inside make_dv_stop_callback makes it a closure: It can access and modify variables (dv_accum, last_t, masses, p, etc.) that are local to make_dv_stop_callback.
        if isnan(last_t[]) # isnan checks if last_t is NaN # last_t[] is initialized on line 679
            last_t[] = t
            return false
        end
        dt = t - last_t[]
        last_t[] = t

        # Sum non-grav force on 'sat'
        Fdict = laser_forces(u, p)
        a = zeros(3)
        for ((_, j), Fj) in Fdict #many Fj from different sat i onto one sat j
            if j == sat
                a .+= Fj ./ masses[sat] # adding acceleration caused by forces from different sast i onto one sat j
            end
        end
        dv_accum .+= a .* dt      # integrate Δv

        # Choose metric
        metric = 0.0 # initialize metric
        if mode == :magnitude
            metric = norm(dv_accum) # total Δv magnitude
        elseif mode == :component
            er, et, en = rtn_basis_from_state(u, sat)
            comp = axis == :R ? dot(dv_accum, er) : axis == :T ? dot(dv_accum, et) : dot(dv_accum, en) # component along chosen RTN axis
            # "?" is a ternary operator, equivalent to if-elseif-else
            metric = comp
        elseif mode == :direction
            û = dirECI/(norm(dirECI)+1e-12) # unit vector in specified ECI direction, ECI means Earth-Centered Inertial
            metric = dot(dv_accum, û)
        else
            error("Unknown dv target mode: $mode")
        end
        return metric >= target # return true if the metric exceeds the target
    end

    affect!(integrator) = terminate!(integrator)
    # affect!() is called when condition() returns true; it stops the integration.
    return DiscreteCallback(condition, affect!; save_positions=(false,false)) # DiscreteCallback function from DifferentialEquations.jl library 
    # Create and return a discrete event callback for the ODE solver that will monitor the simulation at every step.
    # If the condition function returns true, the callback will execute the affect! function (which stops the simulation).
    # The option save_positions=(false,false) tells the solver not to save the state before or after the event.
    # https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#SciMLBase.DiscreteCallback
end

# 1.2. Continuous event: g(u,t) = (chosen Δv metric) - target  crosses 0
"""
    This function creates a continuous callback that stops the integration when a certain Δv target is reached
"""
function make_dv_stop_callback_cont(p; sat::Int=1, target::Float64=15.0,
                                    mode::Symbol=:magnitude, axis::Symbol=:T,
                                    dirECI::SVector{3,Float64}=SVector(1.0,0.0,0.0)) 
    N = p[:N]
    n6 = 6N

    g(u, t, integrator) = begin
        # Δv state is appended as last 3 entries
        dv = @SVector [u[n6+1], u[n6+2], u[n6+3]]
        if mode == :magnitude
            return norm(dv) - target
        elseif mode == :direction
            û = dirECI/(norm(dirECI)+1e-12)
            return dot(dv, û) - target
        elseif mode == :component
            # # RTN axis taken from instantaneous (r,v) of the tracked sat
            # r = @SVector [u[idx(sat,1)], u[idx(sat,2)], u[idx(sat,3)]]
            # v = @SVector [u[idx(sat,4)], u[idx(sat,5)], u[idx(sat,6)]]
            # er = r/(norm(r)+1e-12); en = cross(r,v); en = en/(norm(en)+1e-12); et = cross(en,er)
            # eaxis = axis == :R ? er : axis == :T ? et : en
            # println("Current Δv ", mode, " for sat ", sat, ": ", dot(dv, eaxis), " m/s")
            # return dot(dv, eaxis) - target # this is code is wrong for the same reason showing in ORACLE 2025-10-06 meeting
            dv_component = axis == :R ? dv[1] : axis == :T ? dv[2] : dv[3]
            #println("Current Δv ", mode, " for sat ", sat, ": ", dv_component, " m/s")
            return dv_component - target
        else
            error("Unknown dv target mode: $mode")
        end
    end

    affect!(integrator) = terminate!(integrator) # stop the integration
    return ContinuousCallback(g, affect!; save_positions=(false,false))
    # https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#SciMLBase.ContinuousCallback
end

# 1.3. Continuous event: stop when RTN position metric reaches target
"""
    Create a continuous callback that stops integration when the RTN position
    metric for a satellite reaches a target value.

    Inputs:
        p: parameter dictionary (must include :N)
        sat: satellite index to monitor (1..N)
        target: target value for the chosen metric (meters)
        mode: :magnitude | :component
        axis: for :component mode, one of :R | :T | :N

    Returns:
        ContinuousCallback that terminates the integrator when the metric crosses the target.
"""
function make_r_RTN_stop_callback_cont(p; sat::Int=1, target::Float64=1.0,
                                       mode::Symbol=:component, axis::Symbol=:R)
    N = p[:N]
    n6 = 6N

    # Helper to compute RTN basis from current state for satellite `sat`
    idx6(i, off) = 6*(i-1) + off
    rtn_basis_from_state(u, j) = begin
        r = @SVector [u[idx6(j,1)], u[idx6(j,2)], u[idx6(j,3)]]
        v = @SVector [u[idx6(j,4)], u[idx6(j,5)], u[idx6(j,6)]]
        er = r/(norm(r)+1e-12)
        en = cross(r, v); en = en/(norm(en)+1e-12)
        et = cross(en, er)
        return er, et, en
    end

    g(u, t, integrator) = begin
        # Current ECI position of satellite `sat`
        r = @SVector [u[idx6(sat,1)], u[idx6(sat,2)], u[idx6(sat,3)]]
        er, et, en = rtn_basis_from_state(u, sat)
        # Project position into RTN
        r_R = dot(r, er)
        r_T = dot(r, et)
        r_N = dot(r, en)
        if mode == :magnitude
            return sqrt(r_R^2 + r_T^2 + r_N^2) - target
        elseif mode == :component
            comp = axis == :R ? r_R : axis == :T ? r_T : r_N
            return comp - target
        else
            error("Unknown r RTN target mode: $(mode)")
        end
    end

    affect!(integrator) = terminate!(integrator)
    return ContinuousCallback(g, affect!; save_positions=(false,false))
end

# ----------- 2. OE runners: common scenarios with OE inputs -----------
# 2.1. Open-cavity runner
"""
    This function runs a two-satellite OPEN-CAVITY case with OE inputs.
        Open-cavity with laser cavity between two satellites with power buildup and potential leakage.
            Fint = (2*B*Pin / c) * û    # Internal cavity forces (enhanced)
            Fleak = (f*B*Pin / c) * û   # Additional leakage thrusts
        The first satellite (sat 1) has the laser, and the second satellite (sat 2) has the retroreflector.
        The laser power is Pin (W), and the cavity gain is B (unitless).
    
    Inputs:
    - Orbital elements for two satellites (oe1, oe2)
    - Laser/cavity parameters (Pin, B, mass_kg, etc.)
    - Simulation parameters (use_los, R_atm, atm_clearance, min_range, max_range)
    - Δv target parameters (stop_on_dv, dv_target_sat, dv_target_value, dv_target_mode, etc.)
    - Plotting and verbosity options (verbose, plots)
    - Additional options (e.g., time limits, error tolerances)

    Returns:
    - sol: the solution object from the ODE solver containing the time history of the simulation
    - p: the parameters dictionary used in the simulation
    - (ΔPdict, ΔEdict): dictionaries containing total mechanical impulse and work done by laser/cavity forces
        keyed by (i,j) tuples, where ΔPdict[(i,j)] is the impulse on satellite j due to satellite i

    Defaults per request:
    - mass_kg = 227.0
    - B = 2000.0
    - max_range = 300 km
    - atm_clearance = 5 km above R_atm (default R_E + 100 km)
    - target_dv_mps = 15 m/s (for runtime estimate)
    - Optional stop_on_dv to terminate when |Δv| (or component/direction) hits target.
    - oe1: (a_m=R_EARTH+600e3, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=0.0)
    - oe2: (a_m=R_EARTH+600e3, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=2.0)
    - Pin = 100 kW
"""
function run_open_cavity_oe(; mass_kg=227.0,
                            Pin=1.0e5, B=2000.0, # # laser parameter (100 kW input, 2000 baseline cavity)
                            use_los=true, # use_los – If true, apply line-of-sight and clearance checks before applying forces.
                            R_atm=R_ATMDEF, 
                            atm_clearance=5_000.0, 
                            min_range=0.0, max_range=Inf,
                            target_dv_mps=50.0,       # for runtime estimate only
                            stop_on_dv=false,         # [contain Δv stop] terminate when Δv reaches target
                            dv_target_sat=2,          # [contain Δv stop] which satellite to monitor
                            dv_target_value=50.0,     # [contain Δv stop] target Δv [m/s]
                            dv_target_mode=:magnitude,# [contain Δv stop] :magnitude | :component | :direction
                            dv_target_axis=:T,        # [contain Δv stop] for :component → :R | :T | :N
                            dv_target_dirECI=SVector(1.0,0.0,0.0), # [contain Δv stop] for :direction
                            verbose=true, plots=true,
                            oe1=(a_m=R_EARTH+600e3, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=0.0), # orbital elements for satellite 1 initial location
                            oe2=(a_m=R_EARTH+600e3, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=2.0)) # orbital elements for satellite 2 initial location
                            #oe means orbital elements

    # 1. Obtain initial state u0 from initial OE of satellites
    s1 = state_from_OE(oe1.a_m; e=oe1.e, i_deg=oe1.i_deg, Ω_deg=oe1.Ω_deg, ω_deg=oe1.ω_deg, ν_deg=oe1.ν_deg)
    s2 = state_from_OE(oe2.a_m; e=oe2.e, i_deg=oe2.i_deg, Ω_deg=oe2.Ω_deg, ω_deg=oe2.ω_deg, ν_deg=oe2.ν_deg)
    u0 = vcat(collect(s1), collect(s2))

    # 2. Input parameters
    # 2.1 Common parameters
    N = 2 # number of satellites = 2
    masses = fill(mass_kg, N) # masses is a vector of length N with each entry equal to mass_kg

    # 2.2 Single pass parameters (Deactivate)
    Pm = zeros(N,N)  # power matrix # no single-pass for the cavity pair

    # 2.3. Open-cavity parameters (Active)
    cavity = Dict{Tuple{Int,Int},Dict{Symbol,Any}}()
    cavity[(1,2)] = Dict(:B=>B, :Pin=>Pin,
                        :leak_i_frac=>0.01, :leak_i_along=>true,
                        :leak_j_frac=>0.01, :leak_j_along=>false) # leak fractions and directions
                        # leak_i_frac=>0.01 means 1% of the power is leaked from satellite i
                        # leak_i_along=>true means the leakage from satellite i is along the beam direction

    # 2.4. Build parameters dictionary
    p = Dict(:mu=>MU, :c=>C, :N=>N, :masses=>masses, :Pmatrix=>Pm,
            :use_los=>use_los, :cavity=>cavity,
            :R_atm=>R_atm, :atm_clearance=>atm_clearance,
            :min_range=>min_range, :max_range=>max_range)


    # 3. Runtime estimate and simulation time T
    # Runtime estimate (upper bound if link isn't always on) # runtime estimate is based on nominal acceleration and target Δv
    a_nom = (2*B*Pin/C) / mass_kg # nominal acceleration if always on
    T_est = dv_target_value / max(a_nom, 1e-12) # estimated time to reach target Δv at nominal acceleration

    a1 = oe1.a_m # semi-major axis of the first satellite
    T_orbit = 2π*sqrt(a1^3/MU) # orbital period of the first satellite
    
    T = max(0.25*T_orbit, 1.5*T_est) # run at least 1/4 orbit or 1.5x estimated time (1.5 is the safety factor)
    println("\nEstimated runtime for target Δv ≈ $(dv_target_value) m/s with B=$(B), Pin=$(Pin/1e3) kW, m=$(mass_kg) kg:")
    @printf("  a_nominal ≈ %.6e m/s²,  T_est ≈ %.1f s, orbit ≈ %.1f s  → using T = %.1f s\n",
                a_nom, T_est, T_orbit, T)
    # T_est is the estimated time to reach the target Δv at nominal acceleration
    # T is is the simulation time that it plan to run (1.5 is the safety factor)

    # 4. Solve ODE
    # Build problem (augmented if we need the Δv stop)
    if stop_on_dv # if stop_on_dv is true, we need to augment the state vector to include Δv
        u0_aug = vcat(u0, 0.0, 0.0, 0.0)            # Δv starts at zero
        p_aug  = merge(p, Dict(:track_dv_sat=>dv_target_sat))
        prob   = ODEProblem(nbody_photon_aug!, u0_aug, (0.0, T), p_aug) # use augmented dynamics
        cb     = make_dv_stop_callback_cont(p_aug; sat=dv_target_sat,
                                            target=dv_target_value,
                                            mode=dv_target_mode,
                                            axis=dv_target_axis,
                                            dirECI=dv_target_dirECI)
        sol = solve(prob, Vern9(); reltol=1e-12, abstol=1e-12, callback=cb)
    else
        prob = ODEProblem(nbody_photon!, u0, (0.0, T), p)
        sol  = solve(prob, Vern9(); reltol=1e-12, abstol=1e-12)
    end

    # 5. Plots & reports
    if plots
        plot_orbits(sol; fn="orbits_OC.png")
        plot_angmom(sol, masses; fn="H_OC.png")
        plot_momentum(sol, masses; fn="P_OC.png")
        plot_link_range(sol, p; i=1, j=2, show_clearance=true, units=:km, fn="link_range_1-2.png")
    end

    Eorb = [sum(orbital_energy(u, masses, MU)) for u in sol.u]
    ΔPdict, ΔEdict = evaluate_laser_exchanges(sol, p)
    ΔE_mech = sum(values(ΔEdict))
    println("\nEnergy audit:")
    println("  ΔE_orb_total (MJ)        = ", (Eorb[end]-Eorb[1])/1e6) # total change in orbital energy for all satellites
    println("  ∑ laser/cavity work (MJ) = ",  ΔE_mech/1e6) # total mechanical energy added/removed by lasers/cavity for all satellites
    println("  Balance residual (MJ)    = ", ((Eorb[end]-Eorb[1]) - ΔE_mech)/1e6) 

    print_initial_final_elements(sol, MU; degrees=true)
    for s in 1:N
        report_and_plot_dv_RTN(sol, p; sat=s, fn_prefix="dv_RTN_sat")
    end
    print_los_summary(sol; R_atm=R_atm)
    return sol, p, (ΔPdict, ΔEdict)
end

# 2.2. Single-pass runner
"""
    Optional: SINGLE-PASS runner with OE inputs (for A/B comparisons).
    Note: This does not include any cavity effects, just single-pass laser forces.
    Single-pass beams are direct laser beam from satellite i to satellite j
        F = (P/c) * met.direction  # Simple radiation pressure

    Inputs:
        mass_kg: mass of each satellite (kg)
        P12: laser power from sat 1 to sat 2 (W)
        use_los: whether to use line-of-sight checks
        R_atm: radius of atmosphere (m)
        atm_clearance: minimum clearance above atmosphere (m)
        min_range: minimum range for link (m)
        max_range: maximum range for link (m)
        verbose: whether to print detailed output
        plots: whether to generate plots
        oe1: orbital elements for satellite 1
        oe2: orbital elements for satellite 2
        T_seconds: total simulation time (s)

    Returns:
        sol: the solution object from the ODE solver containing the time history of the simulation
        p: parameters dictionary used in the simulation
        (ΔPdict, ΔEdict): dictionaries of momentum and energy exchanges due to laser forces
"""
function run_single_pass_oe(; mass_kg=227.0,
                             P12=1.0e5, # laser parameter (laser power from sat 1 to sat 2 (W))
                             use_los=true, # whether to use line-of-sight checks
                             R_atm=R_ATMDEF, # radius of atmosphere (m)
                             atm_clearance=5_000.0,
                             min_range=0.0, max_range=300_000.0,
                             verbose=true, 
                             plots=true,
                             oe1=(a_m=R_EARTH+600e3, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=0.0), # orbital elements for satellite 1 initial location
                             oe2=(a_m=R_EARTH+600e3, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=2.0), # orbital elements for satellite 2 initial location
                             T_seconds=2.0*2π*sqrt((R_EARTH+600e3)^3/MU)) # total simulation time (s); defaults to two orbits of the first satellite if not provided

    # 1. Obtain initial state u0 from initial OE of satellites
    s1 = state_from_OE(oe1.a_m; e=oe1.e, i_deg=oe1.i_deg, Ω_deg=oe1.Ω_deg, ω_deg=oe1.ω_deg, ν_deg=oe1.ν_deg)
    s2 = state_from_OE(oe2.a_m; e=oe2.e, i_deg=oe2.i_deg, Ω_deg=oe2.Ω_deg, ω_deg=oe2.ω_deg, ν_deg=oe2.ν_deg)
    u0 = vcat(collect(s1), collect(s2))

    # 2. Input parameters
    # 2.1 Common parameters
    N = 2 
    masses = fill(mass_kg, N)

    # 2.2 Single pass parameters (Active)
    Pm = zeros(N,N) # power matrix for single-pass links
    Pm[1,2] = P12 # laser power from sat 1 to sat 2 (W)

    # 2.3. Open-cavity parameters (Deactivate)
    cavity = Dict{Tuple{Int,Int},Dict{Symbol,Any}}()

    # 2.4. Build parameters dictionary
    p = Dict(:mu=>MU, :c=>C, :N=>N, :masses=>masses, :Pmatrix=>Pm,
            :use_los=>use_los, :cavity=>cavity, :R_atm=>R_atm,
            :atm_clearance=>atm_clearance, :min_range=>min_range, :max_range=>max_range)

    # 3. Solve ODE
    #sol = solve((ODEProblem(nbody_photon!, u0, (0.0, T_seconds), p)), Vern9(); reltol=1e-10, abstol=1e-12)
    prob = (ODEProblem(nbody_photon!, u0, (0.0, T_seconds), p))
    sol = solve(prob, Vern9(); reltol=1e-12, abstol=1e-12)

    # 5. Plots & reports
    if plots
        plot_orbits(sol; fn="orbits_SP.png")
        plot_angmom(sol, masses; fn="H_SP.png")
        plot_momentum(sol, masses; fn="P_SP.png")
    end
    print_initial_final_elements(sol, MU)
    for s in 1:N
        report_and_plot_dv_RTN(sol, p; sat=s, fn_prefix="dv_RTN_SP_sat")
    end
    print_los_summary(sol; R_atm=R_atm)
    return sol, p
end

# 2.3. Multi-satellite OPEN-CAVITY runner
"""
    This function runs a multi-satellite OPEN-CAVITY case with OE inputs.
        This extends run_open_cavity_oe to N>2 satellites with arbitrary cavity and power matrix.

    Inputs:
        - oe_list: list of orbital element named-tuples for each satellite
        - mass_kg: mass of each satellite (kg)
        - masses: optional vector of masses for each satellite (overrides mass_kg if provided)
        - Pm: NxN power matrix (W) for single-pass links (zeros if none)
        - cavity: dictionary defining cavity pairs and their parameters
        - use_los: whether to use line-of-sight checks
        - R_atm: radius of atmosphere (m)
        - atm_clearance: minimum clearance above atmosphere (m)
        - min_range: minimum range for links (m)
        - max_range: maximum range for links (m)
        - stop_on_dv: whether to stop simulation when a Δv target is reached
        - dv_target_sat: which satellite to monitor for Δv
        - dv_target_value: target Δv value (m/s)
        - dv_target_mode: mode for Δv target (:magnitude, :component, :direction)
        - dv_target_axis: axis for component mode (:R, :T, :N)
        - dv_target_dirECI: direction vector for direction mode
        - T_seconds: total simulation time (s); defaults to two orbits of the first satellite if not provided
        - verbose: whether to print detailed output
        - plots: whether to generate plots
        - IMG_DIR: directory to save images

    Returns:
        - sol: the solution object from the ODE solver containing the time history of the simulation
        - p: parameters dictionary used in the simulation
        - (ΔPdict, ΔEdict): dictionaries of momentum and energy exchanges due to laser forces
        - masses: vector of masses for each satellite
""" 
function run_open_cavity_multi(oe_list::AbstractVector;
                               mass_kg::Real=227.0, masses::Union{Nothing,AbstractVector}=nothing, # optional vector of masses for each satellite (overrides mass_kg if provided)
                               Pm::Union{Nothing,AbstractMatrix}=nothing,
                               cavity::Dict{Tuple{Int,Int},Dict{Symbol,Any}}=Dict{Tuple{Int,Int},Dict{Symbol,Any}}(),
                               use_los::Bool=true, R_atm::Real=R_ATMDEF, atm_clearance::Real=5_000.0,
                               min_range::Real=0.0, max_range::Real=Inf,
                               stop_on_dv::Bool=false, dv_target_sat::Int=0, dv_target_value::Real=50.0,
                               dv_target_mode::Symbol=:magnitude, dv_target_axis::Symbol=:T, dv_target_dirECI::SVector{3,Float64}=SVector(1.0,0.0,0.0),
                               stop_on_r::Bool=false, r_target_sat::Int=0, r_target_value::Real=100.0,
                               r_target_mode::Symbol=:component, r_target_axis::Symbol=:R,
                               T_seconds::Union{Nothing,Real}=nothing,
                               verbose::Bool=true, result_plots::Bool=true, target_only::Bool=false,
                               IMG_DIR::String="Kuang's Prototype Code/output/images/",
                               helper_num::Int=10,
                               use_J2::Bool=true, useDrag::Bool=false,
                               gve_schedule::Symbol=:none,
                               gve_target_idx::Union{Nothing,Int}=nothing)

    # 1. Obtain initial state u0 from initial OE of satellites
    u0 = build_u_from_oe(oe_list; μ=MU)

    # 2. Input parameters
    # 2.1 Common parameters
    N = length(oe_list)
    N < 1 && error("oe_list must contain at least one satellite")

    # Helper IDs: helpers are the first `helper_num` IDs in 1:N
    if helper_num > N
        error("helper_num ($(helper_num)) cannot exceed N ($(N))")
    end
    helper_ids = collect(1:helper_num)

    # Target IDs: the remaining IDs after helper_num
    target_ids = collect((helper_num + 1):N)

    mvec = masses === nothing ? fill(mass_kg, N) : collect(masses) # Masses # if masses is not provided, use mass_kg for all satellites
    # length(mvec) == N || error("length(masses) must equal N=$(N)") # check that the length of masses matches the number of satellites
    
    # 2.2 Single pass parameters
    Pm_use = Pm === nothing ? zeros(N, N) : Array(Pm) # if Pm is not provided, use a zero matrix # zero means no single-pass
    # power matrix is a matrix that defines the power of the laser beams between satellites
    # size(Pm_use,1) == N && size(Pm_use,2) == N || error("Pmatrix must be N×N with N=$(N)") # check that the size of Pm matches the number of satellites

    # 2.3. Open-cavity parameters
    # cavity in directly provided in the function input, so don't need to be created here

    # 2.4 Dynamic helper map: key = target sat index, value holds partner and link_type (or nothing)
    #current_helpers = Dict{Int,Union{Nothing,Dict{Symbol,Any}}}(i => nothing for i in 1:N) # Dict{Symbol,Any} is a dictionary with keys of type Symbol and values of any type
    
    # 2.5. Build parameters dictionary
    # gve_target_idx defaults to the first target satellite if not provided explicitly
    _gve_tgt = gve_target_idx === nothing ? (isempty(target_ids) ? N : target_ids[1]) : gve_target_idx
    p = Dict(:mu=>MU, :c=>C, :N=>N, :masses=>mvec, :Pmatrix=>Pm_use,
             :use_los=>use_los, :cavity=>cavity,
             :R_atm=>R_atm, :atm_clearance=>atm_clearance,
             :min_range=>min_range, :max_range=>max_range,
             :helper_ids=>helper_ids, :target_ids=>target_ids,
             :use_J2=>use_J2, :useDrag=>useDrag,
             :gve_schedule=>String(gve_schedule), :gve_target_idx=>_gve_tgt)
            #  :current_helpers=>current_helpers,
            #  :helper_num=>helper_num)

    # 3. Runtime estimate and simulation time T
    # Time horizon: default to two orbits of sat 1 if not provided
    T_orbit = 2π*sqrt(oe_list[1].a_m^3/MU)
    T = T_seconds === nothing ? 2.0*T_orbit : float(T_seconds)
    if verbose
        @printf("\nRuntime horizon: T ≈ %.1f s (%.2f orbits of sat 1)\n", T, T/T_orbit)
    end

    # 4. Solve ODE
    if stop_on_dv && 1 <= dv_target_sat <= N # augmented if dv stop requested on a valid index
        u0_aug = vcat(u0, 0.0, 0.0, 0.0) # Δv starts at zero

        #p_aug  = merge(p, Dict(:track_dv_sat=>dv_target_sat)) # add the target satellite to track into the parameters dictionary
                
        # pass tracking settings into p
        p_aug = merge(p, Dict(
            :track_dv_sat     => dv_target_sat,
            :dv_target_mode    => String(dv_target_mode),  # "component" | "magnitude" | "direction"
            :dv_target_axis    => String(dv_target_axis)  # "R" | "T" | "N"
            # :dv_target_mode   => dv_target_mode,   # :component in your Main5.jl
            # :dv_target_axis   => dv_target_axis,   # :R / :T / :N
            # :dv_target_dirECI => dv_target_dirECI # if you use :direction
        ))
        
        prob   = ODEProblem(nbody_photon_aug!, u0_aug, (0.0, T), p_aug)
        cb     = make_dv_stop_callback_cont(p_aug; sat=dv_target_sat,
                                            target=dv_target_value,
                                            mode=dv_target_mode,
                                            axis=dv_target_axis,
                                            dirECI=dv_target_dirECI)
        sol = solve(prob, Vern9(); reltol=1e-12, abstol=1e-12, callback=cb)
        #p = p_aug  # keep for logs/plots
    
    elseif stop_on_r && 1 <= r_target_sat <= N # augmented if dv stop requested on a valid index
        #u0_aug = vcat(u0, 0.0, 0.0, 0.0) # dummy augmentation to match function signature

        # p_aug  = merge(p, Dict(:track_r_sat=>r_target_sat)) # add the target satellite to track into the parameters dictionary
                
        # prob   = ODEProblem(nbody_photon_aug!, u0_aug, (0.0, T), p_aug)
        prob = ODEProblem(nbody_photon!, u0, (0.0, T), p)
        cb     = make_r_RTN_stop_callback_cont(p; sat=r_target_sat,
                                               target=r_target_value,
                                               mode=r_target_mode,
                                               axis=r_target_axis)
        sol = solve(prob, Vern9(); reltol=1e-12, abstol=1e-12, callback=cb)
        #p = p_aug  # keep for logs/plots
    else
        prob = ODEProblem(nbody_photon!, u0, (0.0, T), p)
        sol  = solve(prob, Vern9(); reltol=1e-12, abstol=1e-12)
    end

    # 5. Reports & plots
    # delete all existing PNG files in IMG_DIR
    # for f in readdir(IMG_DIR)
    #     endswith(f, ".png") && rm(joinpath(IMG_DIR, f))
    # end

    mkpath(IMG_DIR)
    for entry in readdir(IMG_DIR; join=true)
        rm(entry; recursive=true, force=true)
    end

    if result_plots
        plot_orbits(sol; IMG_DIR = IMG_DIR, fn="satellite_orbits.png")
        plot_angmom_two_axes(sol, mvec; IMG_DIR = IMG_DIR, fn="angular_momentum.png")
        plot_momentum_two_axes(sol, mvec; IMG_DIR = IMG_DIR, fn="linear_momentum.png")
        plot_orbit_energy(sol, p; IMG_DIR = IMG_DIR, fn="orbital_energy.png", subplot = true)
        plot_orbit_energy_individual_satellites(sol, p; IMG_DIR = IMG_DIR, fn="orbit_energy_individual.png")
        plot_orbit_energy_total(sol, p; IMG_DIR = IMG_DIR, fn="orbital_energy_total.png")
        # Optional: user can call plot_link_range for specific pairs
    end

    Eorb = [sum(orbital_energy(u, mvec, MU)) for u in sol.u]
    ΔPdict, ΔEdict = evaluate_laser_exchanges(sol, p)
    ΔE_mech = sum(values(ΔEdict))
    if result_plots
        println("\n=================== Energy audit ==================")
        #println("\nEnergy audit:")
        println("  ΔE_orb_total (MJ)        = ", (Eorb[end]-Eorb[1])/1e6)
        println("  ∑ laser/cavity work (MJ) = ",  ΔE_mech/1e6)
        println("  Balance residual (MJ)    = ", ((Eorb[end]-Eorb[1]) - ΔE_mech)/1e6) # should be close to zero
        println("\n")

        println("\n=============== Initial & Final Elements ===============")
        print_initial_final_elements(sol, MU; degrees=true)

        println("\n=============== Historic Data Plots ===============")
        mkpath(IMG_DIR)
        for sub in ("r_RTN_sat","v_RTN_sat","a_RTN_sat","dv_from_laser_in_RTN_for_sat","F_from_laser_in_RTN_for_sat",
            "delta_P_from_laser_in_RTN_for_sat")
            mkpath(joinpath(IMG_DIR, sub))
        end


        if target_only
            for s in target_ids
                mkpath(joinpath(IMG_DIR, "orbital_elements_sat/orbital_elements_sat_$s"))
            end
            for s in target_ids
                report_and_plot_r_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="r_RTN_sat/r_RTN_sat_$s", R_only=true)
                report_and_plot_v_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="v_RTN_sat/v_RTN_sat_$s", R_only=true)
                report_and_plot_a_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="a_RTN_sat/a_RTN_sat_$s", show_a=false, show_a_gravity=false, show_a_laser=true, R_only=true)
                report_and_plot_dv_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="dv_from_laser_in_RTN_for_sat/dv_from_laser_in_RTN_for_sat_$s", RTN_separate=true)
                report_and_plot_F_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="F_from_laser_in_RTN_for_sat/F_from_laser_in_RTN_for_sat_$s", RTN_separate=true)
                report_and_plot_delta_P_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="delta_P_from_laser_in_RTN_for_sat/delta_P_from_laser_in_RTN_for_sat_$s", R_only=true)
                report_and_plot_OE(sol, MU; sat=s, IMG_DIR=IMG_DIR, fn_prefix="orbital_elements_sat/orbital_elements_sat_$s")
                report_and_plot_rp_ra(sol, MU, R_EARTH; sat=s, IMG_DIR=IMG_DIR, fn_prefix="apogee_perigee")
            end
            if N >= 2
                report_and_plot_OE_diff(sol, MU; sat1=1, sat2=2, IMG_DIR=IMG_DIR, fn_prefix="orbital_elements_diff")
            end
        else
            for s in 1:N
                mkpath(joinpath(IMG_DIR, "orbital_elements_sat/orbital_elements_sat_$s"))
            end
            for s in 1:N
                report_and_plot_r_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="r_RTN_sat/r_RTN_sat_$s", R_only=true)
                report_and_plot_v_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="v_RTN_sat/v_RTN_sat_$s", R_only=true)
                report_and_plot_a_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="a_RTN_sat/a_RTN_sat_$s", show_a=false, show_a_gravity=false, show_a_laser=true, R_only=true)
                report_and_plot_dv_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="dv_from_laser_in_RTN_for_sat/dv_from_laser_in_RTN_for_sat_$s", R_only=true)
                report_and_plot_F_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="F_from_laser_in_RTN_for_sat/F_from_laser_in_RTN_for_sat_$s", R_only=true)
                report_and_plot_delta_P_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="delta_P_from_laser_in_RTN_for_sat/delta_P_from_laser_in_RTN_for_sat_$s", R_only=true)
                report_and_plot_OE(sol, MU; sat=s, IMG_DIR=IMG_DIR, fn_prefix="orbital_elements_sat/orbital_elements_sat_$s")
                report_and_plot_rp_ra(sol, MU, R_EARTH; sat=s, IMG_DIR=IMG_DIR, fn_prefix="apogee_perigee")
            end
            if N >= 2
                report_and_plot_OE_diff(sol, MU; sat1=1, sat2=2, IMG_DIR=IMG_DIR, fn_prefix="orbital_elements_diff")
            end
        end

        # println("=============== Link Summary ===============")
        # print_los_summary(sol; R_atm=R_atm)
    end

    # if result_plots == false
    #     println("\n=================== Energy audit ==================")
    #     #println("\nEnergy audit:")
    #     println("  ΔE_orb_total (MJ)        = ", (Eorb[end]-Eorb[1])/1e6)
    #     println("  ∑ laser/cavity work (MJ) = ",  ΔE_mech/1e6)
    #     println("  Balance residual (MJ)    = ", ((Eorb[end]-Eorb[1]) - ΔE_mech)/1e6) # should be close to zero
    #     println("\n")

    #     println("\n=============== Initial & Final Elements ===============")
    #     print_initial_final_elements(sol, MU; degrees=true)

    #     println("\n=============== Historic Data Plots ===============")
    #     mkpath(IMG_DIR)
    #     for sub in ("r_RTN_sat","v_RTN_sat","a_RTN_sat","dv_from_laser_in_RTN_for_sat","F_from_laser_in_RTN_for_sat",
    #         "delta_P_from_laser_in_RTN_for_sat")
    #         mkpath(joinpath(IMG_DIR, sub))
    #     end
    #     for s in 1:N
    #         mkpath(joinpath(IMG_DIR, "orbital_elements_sat/orbital_elements_sat_$s"))
    #     end

    #     for s in 1:N
    #         report_and_plot_dv_RTN(sol, p; sat=s, IMG_DIR=IMG_DIR, fn_prefix="dv_from_laser_in_RTN_for_sat/dv_from_laser_in_RTN_for_sat_$s", R_only=true)
    #     end
    # end


    return sol, p, (ΔPdict, ΔEdict), masses
end
