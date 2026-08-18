# ORACLE Case 2 scenario builder and runner.
# Included by run_oracle_laser_links.jl; uses SpaceAGORA + OracleAnalysis symbols
# already brought into scope by that script.

using LinearAlgebra
using StaticArrays
using DataFrames
using Arrow
using Printf

# Returns extra SaveField entries for cumulative laser ΔV and active helper index.
function _build_laser_save_fields(impulse_tracker::LaserImpulseTracker, laser_model)
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

function _oracle_spacecraft(id::Int, mass_kg::Float64, ic::InitialCondition)
    bus = Link(root=true, m=mass_kg, ref_area=1.0)
    return SpacecraftModel(
        Joint[],
        [bus],
        bus,
        true,
        mass_kg,
        0.0,
        bus.inertia,
        0,
        0,
        ic,
        id,
    )
end

"""
    build_oracle_case_config(opts::OracleOptions, results_directory=nothing)

Build the `SimulationConfiguration` for the ORACLE open-cavity Case 2.
Returns `(args, laser_model, target_period_s)`.
"""
function build_oracle_case_config(opts::OracleOptions, results_directory::Union{String,Nothing}=nothing)
    _validate_options(opts)

    planet          = make_no_gram_planet(opts.planet)
    target_radius_m = planet.Rp_e + opts.target_altitude_km * 1e3
    helper_radius_m = planet.Rp_e + opts.helper_altitude_km * 1e3

    spacecraft = SpacecraftModel[]
    push!(spacecraft, _oracle_spacecraft(1, opts.mass_kg, InitialCondition(
        target_radius_m, opts.target_ecc, opts.target_inclination_deg, 0.0, 0.0, opts.target_nu_deg
    )))
    for helper in 1:opts.helpers
        nu_deg = 360.0 * (helper - 1) / opts.helpers
        push!(spacecraft, _oracle_spacecraft(helper + 1, opts.mass_kg, InitialCondition(
            helper_radius_m, 0.0, opts.helper_inclination_deg, 0.0, 0.0, nu_deg
        )))
    end

    laser_model = OpenCavityLaserLinkModel(
        1,
        collect(2:(opts.helpers + 1));
        range_m=opts.laser_range_km * 1e3,
        power_w=opts.laser_power_w,
        magnification=opts.magnification,
        beta=opts.beta,
        eta=opts.eta,
        schedule=opts.schedule,
    )

    target_period_s = 2pi * sqrt(target_radius_m^3 / planet.μ)
    mission_time_s  = opts.orbits * target_period_s

    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=results_directory !== nothing,
            verbose=false,
            results_directory=results_directory !== nothing ? results_directory : "output",
            generate_plots=false,
            save_csv=results_directory !== nothing && !opts.feather_only,
            normalize=false,
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=mission_time_s,
            orientation_sim=false,
            num_steps_to_save=1000,
            data_rate=max(10.0, mission_time_s / 1000.0),
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false,
        ),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(), laser_model)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-12,
            abstol_orbit=1e-12,
            dt_max_orbit=opts.dt_max_s,
        ),
        solver_config=SolverConfig(solver_mode=:tsit5),
    )
    return args, laser_model, target_period_s
end

"""
    run_oracle_open_cavity_case(opts::OracleOptions)

Run the ORACLE open-cavity Case 2 using SpaceAGORA's native output pipeline.
Returns a named tuple with fields `summary`, `flat_sol`, `sol`, `helper_num`,
`impulse_tracker`, `mu`, `results_dir`.
"""
function run_oracle_open_cavity_case(opts::OracleOptions)
    _validate_options(opts)
    _planet_pre       = make_no_gram_planet(opts.planet)
    _target_radius_m  = _planet_pre.Rp_e + opts.target_altitude_km * 1e3
    _target_period_pre = 2pi * sqrt(_target_radius_m^3 / _planet_pre.μ)
    _T_s = round(Int, _target_period_pre * opts.orbits)
    results_dir = joinpath(opts.output_dir,
        @sprintf("h%.0fkm_t%.0fkm",        opts.helper_altitude_km, opts.target_altitude_km),
        @sprintf("ih%.1fdeg_it%.1fdeg",     opts.helper_inclination_deg, opts.target_inclination_deg),
        @sprintf("N%d",                     opts.helpers),
        @sprintf("T%ds",                    _T_s),
        @sprintf("%s_e%.4f_nu%.4f",         string(opts.schedule), opts.target_ecc, opts.target_nu_deg),
    )

    args, laser_model, target_period_s = build_oracle_case_config(opts, results_dir)

    # Initial OE read directly from InitialCondition (avoids build_initial_conditions).
    ic     = args.dynamics_model.spacecraft[1].initial_condition
    oe0_a  = ic.a; oe0_e = ic.e; oe0_i = ic.i; oe0_raan = ic.Ω

    impulse_tracker = LaserImpulseTracker()
    impulse_cb   = laser_impulse_callback(laser_model, impulse_tracker, opts.mass_kg)
    scheduler_cb = laser_link_scheduler_callback(laser_model)

    all_save_fields = vcat(default_save_fields(args), _build_laser_save_fields(impulse_tracker, laser_model))

    result = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false") do
        run_simulation(
            args;
            isolate_state=false,
            return_solution=true,
            return_solver_metadata=true,
            extra_callbacks=(impulse_cb, scheduler_cb),
            save_fields=all_save_fields,
        )
    end

    sol = result.solution

    feather_df = DataFrame(Arrow.Table(joinpath(results_dir, "simulation_results.feather")))
    flat_sol   = _make_flat_sol_from_feather(feather_df, opts.helpers + 1)

    # Pass planet object so rvtoorbitalelement is used for SMA (consistent with rest of Kuang-ver2).
    orbit_counts   = _orbit_count_from_flat_sol(flat_sol, args.environment_model.planet)
    orbits_elapsed = orbit_counts[end]

    # Use rvtoorbitalelement for the final OE delta (SVector{6}: a,e,i,Ω,ω,ν).
    final_sc = sol.u[end].sc[1]
    oef = rvtoorbitalelement(
        SVector{3,Float64}(final_sc.pos), SVector{3,Float64}(final_sc.vel),
        args.environment_model.planet,
    )

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
        da_m=oef[1] - oe0_a,
        de=oef[2] - oe0_e,
        di_deg=rad2deg(oef[3] - oe0_i),
        draan_deg=rad2deg(oef[4] - oe0_raan),
        retcode=sol.retcode,
        solver=result.solver_trace[end].solver,
    )

    opts.feather_only || _write_csv!(joinpath(results_dir, "summary.csv"), [summary])

    return (summary=summary, flat_sol=flat_sol, sol=sol,
            helper_num=opts.helpers, impulse_tracker=impulse_tracker,
            mu=args.environment_model.planet.μ, results_dir=results_dir)
end

"""
    plot_open_cavity_results(result, opts; IMG_DIR, target_only)

Generate the full diagnostic figure suite for a completed Case-2 simulation.
"""
function plot_open_cavity_results(result, opts::OracleOptions;
                                  IMG_DIR::String=joinpath(result.results_dir, "images"),
                                  target_only::Bool=false)
    flat_sol = result.flat_sol
    tracker  = result.impulse_tracker
    mu       = result.mu
    N        = opts.helpers + 1

    R_earth   = 6_378_137.0
    R_atm_def = R_earth + 100e3
    mvec = fill(opts.mass_kg, N)
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
        :tracker      => tracker,
        :sa_sol       => nothing,
    )

    sats_to_plot = target_only ? [1] : collect(1:N)

    mkpath(IMG_DIR)
    for entry in readdir(IMG_DIR; join=true); rm(entry; recursive=true, force=true); end
    for sub in ("r_RTN_sat", "v_RTN_sat", "a_RTN_sat",
                "dv_from_laser_in_RTN_for_sat",
                "F_from_laser_in_RTN_for_sat",
                "delta_P_from_laser_in_RTN_for_sat")
        mkpath(joinpath(IMG_DIR, sub))
    end
    for s in sats_to_plot
        mkpath(joinpath(IMG_DIR, "orbital_elements_sat", "orbital_elements_sat_$s"))
    end

    plot_orbits(flat_sol; IMG_DIR=IMG_DIR, fn="satellite_orbits.png")
    plot_angmom_two_axes(flat_sol, mvec; IMG_DIR=IMG_DIR, fn="angular_momentum.png")
    plot_momentum_two_axes(flat_sol, mvec; IMG_DIR=IMG_DIR, fn="linear_momentum.png")
    plot_orbit_energy(flat_sol, p; IMG_DIR=IMG_DIR, fn="orbital_energy.png", subplot=true)
    plot_orbit_energy_individual_satellites(flat_sol, p; IMG_DIR=IMG_DIR, fn="orbit_energy_individual.png")
    plot_orbit_energy_total(flat_sol, p; IMG_DIR=IMG_DIR, fn="orbital_energy_total.png")
    plot_orbital_energy_change(flat_sol, p; IMG_DIR=IMG_DIR, fn="orbital_energy_change.png", total_only=true)
    plot_laser_dE_time_series(flat_sol, p; IMG_DIR=IMG_DIR, fn="laser_dE_timeseries.png", total_only=true)
    plot_laser_dP_time_series(flat_sol, p; IMG_DIR=IMG_DIR, fn="laser_dP_timeseries.png", total_only=true)

    Eorb  = [sum(orbital_energy(u, mvec, mu)) for u in flat_sol.u]
    ΔPdict, ΔEdict = evaluate_laser_exchanges(flat_sol, p)
    ΔE_mech = sum(values(ΔEdict))
    println("\n=================== Energy audit ==================")
    println("  ΔE_orb_total (MJ)        = ", (Eorb[end]-Eorb[1])/1e6)
    println("  ∑ laser/cavity work (MJ) = ",  ΔE_mech/1e6)
    println("  Balance residual (MJ)    = ", ((Eorb[end]-Eorb[1]) - ΔE_mech)/1e6)

    println("\n=============== Δv RTN Start → End ===============")
    let (_t_dv, _Δv_hist) = delta_v_RTN_time_series(flat_sol, p)
        for s in sats_to_plot
            dv = _Δv_hist[s]
            @printf("  Sat %d  Δv_R (m/s): %+.4f → %+.4f\n", s, dv[1,1], dv[1,end])
            @printf("  Sat %d  Δv_T (m/s): %+.4f → %+.4f\n", s, dv[2,1], dv[2,end])
            @printf("  Sat %d  Δv_N (m/s): %+.4f → %+.4f\n", s, dv[3,1], dv[3,end])
            @printf("  Sat %d  |Δv| (m/s): %+.4f → %+.4f\n", s, norm(dv[:,1]), norm(dv[:,end]))
        end
    end

    println("\n=============== Initial & Final Elements ===============")
    print_initial_final_elements(flat_sol, mu; degrees=true)

    println("\n=============== Historic Data Plots ===============")
    for s in sats_to_plot
        report_and_plot_r_RTN(flat_sol, p; sat=s, IMG_DIR=IMG_DIR,
                              fn_prefix="r_RTN_sat/r_RTN_sat_$s", R_only=true)
        report_and_plot_v_RTN(flat_sol, p; sat=s, IMG_DIR=IMG_DIR,
                              fn_prefix="v_RTN_sat/v_RTN_sat_$s", R_only=true)
        report_and_plot_a_RTN(flat_sol, p; sat=s, IMG_DIR=IMG_DIR,
                              fn_prefix="a_RTN_sat/a_RTN_sat_$s",
                              show_a=false, show_a_gravity=false, show_a_laser=(s==1), R_only=true)
        if s == 1
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
