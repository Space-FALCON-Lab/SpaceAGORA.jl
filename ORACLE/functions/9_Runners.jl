# Run the simulation for the given options, returning a summary, time-series DataFrame, and the full solution object.
function run_open_cavity_case(opts::OracleCase2Options)
    # --- Block 1: build config ---
    # Build the SimulationConfiguration, laser model, and initial orbital period from opts
    args, laser_model, target_period_s = build_case_config(opts)

    # --- Block 2: initial orbital elements ---
    u0 = SimulationEngine.build_initial_conditions(args)  # convert args → initial state vector
    r0 = SVector{3, Float64}(u0.sc[1].pos)                # target initial position [m]
    v0 = SVector{3, Float64}(u0.sc[1].vel)                # target initial velocity [m/s]
    oe0 = _rv_to_elements(r0, v0, args.environment_model.planet.μ)  # initial orbital elements (a,e,i,raan)

    # --- Block 3: callbacks ---
    # impulse_cb must come FIRST so it reads the active_helper_idx that was set
    # for the step just completed, before the scheduler updates it for the next.
    impulse_tracker = _LaserImpulseTracker()               # accumulates dV in RTN at every ODE step
    impulse_cb  = _make_laser_impulse_callback(laser_model, impulse_tracker, opts.mass_kg)
    scheduler_cb = laser_link_scheduler_callback(laser_model)  # updates which helper fires next

    # --- Block 4: run the ODE solver ---
    result = run_simulation(
        args;
        isolate_state=false,       # no deep copy — laser_model IS the live model, so the callback can use it directly
        return_solution=true,      # keep the full solution object (dense output + saved snapshots)
        return_solver_metadata=true,  # keep solver name and step stats for the summary
        extra_callbacks=(impulse_cb, scheduler_cb),
    )

    # --- Block 5: post-simulation quantities ---
    sol = result.solution
    orbit_counts   = _orbit_count_from_sol(sol, args.environment_model.planet.μ)  # varying-period orbit counter
    orbits_elapsed = orbit_counts[end]   # total orbits completed (accounts for changing period)
    final_state = sol.u[end].sc[1]       # target state at the last saved time step
    rf = SVector{3, Float64}(final_state.pos)  # final position [m]
    vf = SVector{3, Float64}(final_state.vel)  # final velocity [m/s]
    oef = _rv_to_elements(rf, vf, args.environment_model.planet.μ)  # final orbital elements
    # laser_model is the live model (isolate_state=false), so its final scheduler state is already up to date

    # --- Block 6: pack results into a summary named tuple ---
    # Each field becomes one column when written to the summary CSV
    summary = (
        case_id=_case_id(opts),                          # unique string label for this run
        helpers=opts.helpers,                            # number of helper satellites
        helper_altitude_km=opts.helper_altitude_km,      # helper orbital altitude
        schedule=opts.schedule,                          # which helper-selection algorithm was used
        target_period_s=target_period_s,                 # initial orbital period of the target [s]
        active_helper=laser_model.active_helper_idx,     # which helper was firing at end of sim (0 = none)
        activations=laser_model.link_activation_count,   # total number of helper switches
        active_steps=laser_model.active_link_step_count, # total ODE steps where the laser was on
        target_altitude_km=opts.target_altitude_km,
        target_inclination_deg=opts.target_inclination_deg,
        helper_inclination_deg=opts.helper_inclination_deg,
        orbits=opts.orbits,                              # requested number of orbits (input)
        orbits_elapsed=orbits_elapsed,                   # actual orbits completed (varying period)
        laser_range_km=opts.laser_range_km,
        laser_power_w=opts.laser_power_w,
        magnification=opts.magnification,
        beta=opts.beta,
        eta=opts.eta,
        mass_kg=opts.mass_kg,
        dt_max_s=opts.dt_max_s,
        dv_r_mps=impulse_tracker.dv_R,   # total radial ΔV from laser [m/s]
        dv_t_mps=impulse_tracker.dv_T,   # total along-track ΔV from laser [m/s]
        dv_n_mps=impulse_tracker.dv_N,   # total cross-track ΔV from laser [m/s]
        da_m=oef.a - oe0.a,              # change in semi-major axis [m]
        de=oef.e - oe0.e,                # change in eccentricity
        di_deg=rad2deg(oef.i - oe0.i),   # change in inclination [deg]
        draan_deg=rad2deg(oef.raan - oe0.raan),  # change in RAAN [deg]
        retcode=sol.retcode,             # ODE solver exit status (Success, MaxIters, etc.)
        solver=result.solver_trace[end].solver,  # name of the solver used
    )

    # --- Block 7: build the time-series DataFrame and return ---
    timeseries = _build_timeseries_dataframe(
        opts,
        sol,
        oe0,
        orbit_counts,
        args.environment_model.planet.μ,
        impulse_tracker,
    )
    return (summary=summary, timeseries=timeseries, sol=sol,
            helper_num=opts.helpers, impulse_tracker=impulse_tracker,
            mu=args.environment_model.planet.μ)
end

# =============================================================================
# SaveField helpers for laser-specific quantities (used by native output pipeline)
# =============================================================================

# Returns extra SaveField entries that capture cumulative laser ΔV and the active
# helper index at every ODE save step.  These are merged with SpaceAGORA's default
# save fields so the standard output bundle (feather + csv + manifest) contains both
# the usual trajectory columns AND these laser-specific columns.
function _build_laser_save_fields(impulse_tracker::_LaserImpulseTracker, laser_model)
    return SaveField[
        SaveField(
            :dv_r_accumulated,
            (u, t, integrator) -> impulse_tracker.dv_R;
            per_satellite=false, column_prefix="dv_r_accumulated",
        ),
        SaveField(
            :dv_t_accumulated,
            (u, t, integrator) -> impulse_tracker.dv_T;
            per_satellite=false, column_prefix="dv_t_accumulated",
        ),
        SaveField(
            :dv_n_accumulated,
            (u, t, integrator) -> impulse_tracker.dv_N;
            per_satellite=false, column_prefix="dv_n_accumulated",
        ),
        SaveField(
            :laser_active_helper,
            (u, t, integrator) -> laser_model.active_helper_idx;
            per_satellite=false, column_prefix="laser_active_helper",
        ),
    ]
end

# =============================================================================
# Native-output runner — uses SpaceAGORA's output pipeline (Option B)
# =============================================================================

# Like run_open_cavity_case, but delegates timeseries output to SpaceAGORA's native
# pipeline (simulation_results.feather + .csv + .manifest.toml in a per-scenario
# subfolder).  The custom post-run summary (one row with totals) is still written as
# summary.csv alongside the SpaceAGORA files.
#
# Activate with --native-output on the command line.
# The return value has no `timeseries` field; use `results_dir` to find the output.
function run_open_cavity_case_native(opts::OracleCase2Options)
    # --- Compute per-scenario output directory BEFORE building args ---
    # (We need target_period_s to build the stem, and that requires knowing the planet radius.)
    _validate_options(opts)
    _planet_pre        = make_no_gram_planet(:earth)
    _target_radius_m   = _planet_pre.Rp_e + opts.target_altitude_km * 1e3
    _target_period_pre = 2pi * sqrt(_target_radius_m^3 / _planet_pre.μ)
    _T_s = round(Int, _target_period_pre * opts.orbits)
    results_dir = joinpath(opts.output_dir,
        @sprintf("h%.0fkm_t%.0fkm",        opts.helper_altitude_km, opts.target_altitude_km),
        @sprintf("ih%.1fdeg_it%.1fdeg",     opts.helper_inclination_deg, opts.target_inclination_deg),
        @sprintf("N%d",                     opts.helpers),
        @sprintf("T%ds",                    _T_s),
    )

    # --- Block 1: build config with native output enabled ---
    args, laser_model, target_period_s = build_case_config_native(opts, results_dir)

    # --- Block 2: initial orbital elements ---
    u0  = SimulationEngine.build_initial_conditions(args)
    r0  = SVector{3, Float64}(u0.sc[1].pos)
    v0  = SVector{3, Float64}(u0.sc[1].vel)
    oe0 = _rv_to_elements(r0, v0, args.environment_model.planet.μ)

    # --- Block 3: callbacks ---
    impulse_tracker = _LaserImpulseTracker()
    impulse_cb   = _make_laser_impulse_callback(laser_model, impulse_tracker, opts.mass_kg)
    scheduler_cb = laser_link_scheduler_callback(laser_model)

    # --- Block 4: build save fields (SpaceAGORA defaults + laser extras) ---
    laser_fields = _build_laser_save_fields(impulse_tracker, laser_model)
    all_save_fields = vcat(SimulationModel.default_save_fields(args), laser_fields)

    # --- Block 5: run — SpaceAGORA writes feather + csv + manifest automatically ---
    result = run_simulation(
        args;
        isolate_state=false,
        return_solution=true,
        return_solver_metadata=true,
        extra_callbacks=(impulse_cb, scheduler_cb),
        save_fields=all_save_fields,
    )

    # --- Block 6: post-simulation quantities (same as run_open_cavity_case) ---
    sol            = result.solution
    orbit_counts   = _orbit_count_from_sol(sol, args.environment_model.planet.μ)
    orbits_elapsed = orbit_counts[end]
    final_state    = sol.u[end].sc[1]
    rf  = SVector{3, Float64}(final_state.pos)
    vf  = SVector{3, Float64}(final_state.vel)
    oef = _rv_to_elements(rf, vf, args.environment_model.planet.μ)

    # --- Block 7: summary row (SpaceAGORA has no concept of a one-row summary) ---
    summary = (
        case_id=_case_id(opts),
        helpers=opts.helpers,
        helper_altitude_km=opts.helper_altitude_km,
        schedule=opts.schedule,
        target_period_s=target_period_s,
        active_helper=laser_model.active_helper_idx,
        activations=laser_model.link_activation_count,
        active_steps=laser_model.active_link_step_count,
        target_altitude_km=opts.target_altitude_km,
        target_inclination_deg=opts.target_inclination_deg,
        helper_inclination_deg=opts.helper_inclination_deg,
        orbits=opts.orbits,
        orbits_elapsed=orbits_elapsed,
        laser_range_km=opts.laser_range_km,
        laser_power_w=opts.laser_power_w,
        magnification=opts.magnification,
        beta=opts.beta,
        eta=opts.eta,
        mass_kg=opts.mass_kg,
        dt_max_s=opts.dt_max_s,
        dv_r_mps=impulse_tracker.dv_R,
        dv_t_mps=impulse_tracker.dv_T,
        dv_n_mps=impulse_tracker.dv_N,
        da_m=oef.a - oe0.a,
        de=oef.e - oe0.e,
        di_deg=rad2deg(oef.i - oe0.i),
        draan_deg=rad2deg(oef.raan - oe0.raan),
        retcode=sol.retcode,
        solver=result.solver_trace[end].solver,
    )

    # Write the one-row summary CSV (skipped when feather_only=true)
    opts.feather_only || _write_csv!(joinpath(results_dir, "summary.csv"), [summary])

    return (summary=summary, sol=sol,
            helper_num=opts.helpers, impulse_tracker=impulse_tracker,
            mu=args.environment_model.planet.μ, results_dir=results_dir)
end

# =============================================================================
# Plotting helper: mirrors the prototype's run_open_cavity_multi plot block
# but adapts to the SpaceAGORA solution format.
#
# Usage (in main / paper-grid loop):
#   result = run_open_cavity_case(opts)
#   plot_open_cavity_results(result, opts)
# =============================================================================
"""
    plot_open_cavity_results(result, opts;
                             IMG_DIR = "output/oracle_case2_laser_links/images/",
                             target_only = false)

Generate the full diagnostic figure suite for a completed Case-2 simulation.

The function:
1. Wraps SpaceAGORA's `result.sol` in a flat-vector `_FlatSol` adapter so all
   ported 7_Plots.jl functions work without modification.
2. Constructs a compatible parameter dictionary `p` that includes `:tracker`
   (the impulse tracker) so laser-force plots are driven by the accumulated dV
   history instead of re-evaluating the ODE right-hand side.
3. Replicates the prototype `run_open_cavity_multi` result_plots block:
   orbits, energy, momentum, orbital elements, RTN kinematics, dV, force, ΔP.

`target_only = true` limits the per-satellite series to the target (sc[1]).
"""
function plot_open_cavity_results(result, opts::OracleCase2Options;
                                  IMG_DIR::String = joinpath(
                                      normpath(joinpath(@__DIR__, "..", "..", "output",
                                               "oracle_case2_laser_links")), "images"),
                                  target_only::Bool = false)

    sa_sol  = result.sol           # SpaceAGORA structured solution
    tracker = result.impulse_tracker
    mu      = result.mu            # gravitational parameter [m³/s²]
    N       = opts.helpers + 1     # total spacecraft: 1 target + N helpers

    # -----------------------------------------------------------------
    # 1. Build flat-vector sol adapter (N total spacecraft)
    # -----------------------------------------------------------------
    flat_sol = _make_flat_sol(sa_sol, N)

    # -----------------------------------------------------------------
    # 2. Build a plotting-compatible parameter dictionary
    # -----------------------------------------------------------------
    R_earth  = 6_378_137.0       # [m]  WGS-84 mean equatorial radius
    R_atm_def = R_earth + 100e3  # [m]  100 km atmosphere top
    mvec = fill(opts.mass_kg, N) # all spacecraft share the same mass
    p = Dict{Symbol,Any}(
        :N            => N,
        :masses       => mvec,
        :mu           => mu,
        :c            => 299_792_458.0,
        :eta          => opts.eta,
        :helper_ids   => collect(2:N),
        :target_ids   => [1],
        :R_atm        => R_atm_def,
        :atm_clearance => 5_000.0,
        :min_range    => 0.0,
        :max_range    => opts.laser_range_km * 1e3,
        :tracker      => tracker,    # used by delta_v / force / laser-exchange functions
        :sa_sol       => sa_sol,     # SpaceAGORA structured sol for velocity lookup in energy calc
    )

    target_ids = [1]
    helper_ids = collect(2:N)
    sats_to_plot = target_only ? target_ids : collect(1:N)

    # -----------------------------------------------------------------
    # 3. Clear output directory and recreate sub-folders
    # -----------------------------------------------------------------
    mkpath(IMG_DIR)
    for entry in readdir(IMG_DIR; join=true)
        rm(entry; recursive=true, force=true)
    end
    for sub in ("r_RTN_sat", "v_RTN_sat", "a_RTN_sat",
                "dv_from_laser_in_RTN_for_sat",
                "F_from_laser_in_RTN_for_sat",
                "delta_P_from_laser_in_RTN_for_sat")
        mkpath(joinpath(IMG_DIR, sub))
    end
    for s in sats_to_plot
        mkpath(joinpath(IMG_DIR, "orbital_elements_sat", "orbital_elements_sat_$s"))
    end

    # -----------------------------------------------------------------
    # 4. Overview plots (same calls as prototype's run_open_cavity_multi)
    # -----------------------------------------------------------------
    plot_orbits(flat_sol; IMG_DIR=IMG_DIR, fn="satellite_orbits.png")
    plot_angmom_two_axes(flat_sol, mvec; IMG_DIR=IMG_DIR, fn="angular_momentum.png")
    plot_momentum_two_axes(flat_sol, mvec; IMG_DIR=IMG_DIR, fn="linear_momentum.png")
    plot_orbit_energy(flat_sol, p; IMG_DIR=IMG_DIR, fn="orbital_energy.png", subplot=true)
    plot_orbit_energy_individual_satellites(flat_sol, p; IMG_DIR=IMG_DIR, fn="orbit_energy_individual.png")
    plot_orbit_energy_total(flat_sol, p; IMG_DIR=IMG_DIR, fn="orbital_energy_total.png")
    plot_orbital_energy_change(flat_sol, p; IMG_DIR=IMG_DIR, fn="orbital_energy_change.png", total_only=true)
    plot_laser_dE_time_series(flat_sol, p; IMG_DIR=IMG_DIR, fn="laser_dE_timeseries.png", total_only=true)
    plot_laser_dP_time_series(flat_sol, p; IMG_DIR=IMG_DIR, fn="laser_dP_timeseries.png", total_only=true)

    # -----------------------------------------------------------------
    # 5. Energy audit (print to console)
    # -----------------------------------------------------------------
    Eorb    = [sum(orbital_energy(u, mvec, mu)) for u in flat_sol.u]
    ΔPdict, ΔEdict = evaluate_laser_exchanges(sa_sol, p)
    ΔE_mech = sum(values(ΔEdict))
    println("\n=================== Energy audit ==================")
    println("  ΔE_orb_total (MJ)        = ", (Eorb[end]-Eorb[1])/1e6)
    println("  ∑ laser/cavity work (MJ) = ",  ΔE_mech/1e6)
    println("  Balance residual (MJ)    = ", ((Eorb[end]-Eorb[1]) - ΔE_mech)/1e6)

    # -----------------------------------------------------------------
    # 6. ΔV RTN start→end summary
    # -----------------------------------------------------------------
    println("\n=============== Δv RTN Start → End ===============")
    let (_t_dv, _Δv_hist) = delta_v_RTN_time_series(flat_sol, p)
        for s in sats_to_plot
            dv = _Δv_hist[s]
            @printf("  Sat %d  Δv_R (m/s): %+.4f → %+.4f\n", s, dv[1,1],   dv[1,end])
            @printf("  Sat %d  Δv_T (m/s): %+.4f → %+.4f\n", s, dv[2,1],   dv[2,end])
            @printf("  Sat %d  Δv_N (m/s): %+.4f → %+.4f\n", s, dv[3,1],   dv[3,end])
            @printf("  Sat %d  |Δv| (m/s): %+.4f → %+.4f\n", s,
                    norm(dv[:,1]), norm(dv[:,end]))
        end
    end

    # -----------------------------------------------------------------
    # 7. Initial & final orbital elements
    # -----------------------------------------------------------------
    println("\n=============== Initial & Final Elements ===============")
    print_initial_final_elements(flat_sol, mu; degrees=true)

    # -----------------------------------------------------------------
    # 8. Per-satellite time-series plots
    # -----------------------------------------------------------------
    println("\n=============== Historic Data Plots ===============")
    for s in sats_to_plot
        report_and_plot_r_RTN(flat_sol, p; sat=s, IMG_DIR=IMG_DIR,
                              fn_prefix="r_RTN_sat/r_RTN_sat_$s", R_only=true)
        report_and_plot_v_RTN(flat_sol, p; sat=s, IMG_DIR=IMG_DIR,
                              fn_prefix="v_RTN_sat/v_RTN_sat_$s", R_only=true)
        report_and_plot_a_RTN(flat_sol, p; sat=s, IMG_DIR=IMG_DIR,
                              fn_prefix="a_RTN_sat/a_RTN_sat_$s",
                              show_a=false, show_a_gravity=false,
                              show_a_laser=(s==1), R_only=true)
        if s == 1  # tracker data only covers the target
            report_and_plot_dv_RTN(flat_sol, p; sat=s, IMG_DIR=IMG_DIR,
                                   fn_prefix="dv_from_laser_in_RTN_for_sat/dv_from_laser_in_RTN_for_sat_$s",
                                   RTN_separate=true)
            report_and_plot_F_RTN(flat_sol, p; sat=s, IMG_DIR=IMG_DIR,
                                  fn_prefix="F_from_laser_in_RTN_for_sat/F_from_laser_in_RTN_for_sat_$s",
                                  RTN_separate=true)
            report_and_plot_delta_P_RTN(flat_sol, p; sat=s, IMG_DIR=IMG_DIR,
                                        fn_prefix="delta_P_from_laser_in_RTN_for_sat/delta_P_from_laser_in_RTN_for_sat_$s",
                                        R_only=true)
        end
        report_and_plot_OE(flat_sol, mu; sat=s, IMG_DIR=IMG_DIR,
                           fn_prefix="orbital_elements_sat/orbital_elements_sat_$s")
        report_and_plot_rp_ra(flat_sol, mu, R_earth; sat=s, IMG_DIR=IMG_DIR,
                              fn_prefix="apogee_perigee")
    end
    if N >= 2
        report_and_plot_OE_diff(flat_sol, mu; sat1=1, sat2=2, IMG_DIR=IMG_DIR,
                                fn_prefix="orbital_elements_diff")
    end

    println("\nPlots saved to: $IMG_DIR")
    return nothing
end

