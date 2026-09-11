#!/usr/bin/env julia
# ORACLE Case 2 — laser-link deorbit simulation
# Uses the build_constellation API added to 4_SpaceAGORA.jl-main.
# Equivalent to 2_SpaceAGORA.jl/ORACLE/run_case2_laser_links.jl (single-run mode).
# Run:  julia --project=. examples/oracle_laser_links.jl [options]

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "oracle_laser_plots.jl"))
using StaticArrays
using LinearAlgebra
using Printf
using SpaceAGORA.TelemetryVerification: rvtoorbitalelement

# ── Options ────────────────────────────────────────────────────────────────────

Base.@kwdef struct OracleLaserOptions
    helpers::Int               = 10
    helper_altitude_km::Float64 = 1050.0
    target_altitude_km::Float64 = 1000.0
    target_inclination_deg::Float64 = 0.0
    helper_inclination_deg::Float64 = 0.0
    target_nu_deg::Float64     = 0.0
    target_ecc::Float64        = 0.0
    orbits::Float64            = 10.0
    schedule::Symbol           = :naive_next_entering
    laser_range_km::Float64    = 200.0
    laser_power_w::Float64     = 10_000.0
    magnification::Float64     = 100.0
    beta::Float64              = 1.0
    eta::Float64               = 1.0
    mass_kg::Float64           = 227.0
    dt_max_s::Float64          = 10.0
    output_dir::String         = joinpath(REPO_ROOT, "output", "oracle_laser_links")
    timeseries_points::Int     = 1001
end

function _usage()
    return """
    Usage: julia --project=. examples/oracle_laser_links.jl [options]

    Options:
      --helpers N                        (default: 10)
      --helper-altitude-km KM            (default: 1050.0)
      --target-altitude-km KM            (default: 1000.0)
      --target-inclination-deg DEG       (default: 0.0)
      --helper-inclination-deg DEG       (default: 0.0)
      --target-nu-deg DEG                (default: 0.0)
      --target-ecc VALUE                 (default: 0.0)
      --orbits N                         (default: 10.0)
      --schedule SYMBOL                  naive_next_entering|positive_along_track|gve_sma|... (default: naive_next_entering)
      --laser-range-km KM                (default: 200.0)
      --laser-power-w W                  (default: 10000.0)
      --magnification VALUE              (default: 100.0)
      --beta VALUE                       (default: 1.0)
      --eta VALUE                        (default: 1.0)
      --mass-kg KG                       (default: 227.0)
      --dt-max-s SEC                     (default: 10.0)
      --output-dir PATH
      --timeseries-points N              (default: 1001)
    """
end

const _INT_OPTS_ORACLE    = (:helpers, :timeseries_points)                      # parsed as Int
const _SYMBOL_OPTS_ORACLE = (:schedule,)                                        # parsed as Symbol
const _PATH_OPTS_ORACLE   = (:output_dir,)                                      # parsed as abspath
const _FLOAT_OPTS_ORACLE  = (                                                   # parsed as Float64
    :helper_altitude_km, :target_altitude_km, :target_inclination_deg,
    :helper_inclination_deg, :target_nu_deg, :target_ecc, :orbits,
    :laser_range_km, :laser_power_w, :magnification, :beta, :eta,
    :mass_kg, :dt_max_s,
)

function _parse_options(argv)::OracleLaserOptions
    opts = Dict{Symbol, Any}()
    i = 1
    while i <= length(argv)
        arg = argv[i]
        arg in ("--help", "-h") && (println(_usage()); exit(0))
        startswith(arg, "--") || throw(ArgumentError("Unexpected argument '$arg'."))
        key = Symbol(replace(arg[3:end], '-' => '_'))
        i < length(argv) || throw(ArgumentError("Missing value for $arg."))
        val = argv[i + 1]
        if     key in _INT_OPTS_ORACLE    opts[key] = parse(Int, val)
        elseif key in _SYMBOL_OPTS_ORACLE opts[key] = Symbol(val)
        elseif key in _PATH_OPTS_ORACLE   opts[key] = abspath(val)
        elseif key in _FLOAT_OPTS_ORACLE  opts[key] = parse(Float64, val)
        else   throw(ArgumentError("Unknown option $arg."))
        end
        i += 2
    end
    return OracleLaserOptions(; opts...)
end

# ── Spacecraft factory ─────────────────────────────────────────────────────────

# Minimal single-body spacecraft matching the ORACLE bus: one rigid box, no panels.
# id is a placeholder — build_constellation reassigns it by position.
function _make_oracle_sc(initial_conditions, mass_kg::Float64)
    bus = Link(root=true, m=mass_kg)
    return SpacecraftModel(
        joints            = Joint[],
        links             = [bus],
        root              = bus,
        instant_actuation = true,
        # dry_mass here is computed as the sum of link masses (= mass_kg)
        prop_mass         = 0.0,
        inertia_tensor    = bus.inertia,
        n_reaction_wheels = 0,
        n_thrusters       = 0,
        initial_condition = initial_conditions,
        id                = 0,
    )
end
# ── Main runner ────────────────────────────────────────────────────────────────

function run_oracle_laser(opts::OracleLaserOptions)
    planet = make_no_gram_planet(:earth)

    # 1. Build constellation: one target + N evenly spaced helpers
    target_radius_m = planet.Rp_e + opts.target_altitude_km * 1e3
    helper_radius_m = planet.Rp_e + opts.helper_altitude_km * 1e3

    target_ic = InitialCondition(
        ra = target_radius_m * (1.0 + opts.target_ecc), rp = target_radius_m * (1.0 - opts.target_ecc),
        i  = opts.target_inclination_deg, ω = 0.0, Ω = 0.0, ν = opts.target_nu_deg
    )

    helper_ics = [
        InitialCondition(
            ra = helper_radius_m, rp = helper_radius_m,
            i  = opts.helper_inclination_deg, ω = 0.0, Ω = 0.0,
            ν  = 360.0 * (k - 1) / opts.helpers
        )
        for k in 1:opts.helpers
    ]

    constellation = build_constellation(_make_oracle_sc.(vcat([target_ic], helper_ics), opts.mass_kg))

    # 2. Laser links: one LaserLinkModel per (target, helper) pair (uses positional indices = satellite array position)
    laser_links = [
        build_LaserLinkModel(
            sat_i          = 1,
            sat_j          = h,
            range_m        = opts.laser_range_km * 1e3,
            schedule       = opts.schedule,
            laser_type     = :thruster,
            power_w        = opts.laser_power_w,
            magnification  = opts.magnification,
            beta           = opts.beta,
            eta            = opts.eta,
        )
        for h in 2:(opts.helpers + 1)
    ]

    # 3. Mission time
    target_period_s = 2π * sqrt(target_radius_m^3 / planet.μ)
    mission_time_s  = opts.orbits * target_period_s

    # 4. SimulationConfiguration
    results_dir = joinpath(opts.output_dir,
        @sprintf("h%.0fkm_t%.0fkm_N%d_%s",
            opts.helper_altitude_km, opts.target_altitude_km,
            opts.helpers, string(opts.schedule))
    )
    mkpath(results_dir)

    args = SimulationConfiguration(
        simulation_settings = SimulationSettings(
            results           = true,
            verbose           = false,
            results_directory = results_dir,
            generate_plots    = false,
            save_csv          = true,
            normalize         = false,
        ),
        mission_configuration = MissionConfiguration(
            mission_type       = MissionTime,
            keplerian          = true,
            number_of_orbits   = 1,
            mission_time       = mission_time_s,
            orientation_sim    = false,
            num_steps_to_save  = opts.timeseries_points,
            data_rate          = max(10.0, mission_time_s / opts.timeseries_points),
        ),
        environment_model = EnvironmentModel(
            planet           = planet,
            EI               = 120.0,
            density_model    = NoAtmosphereModel(),
            ephemerides_model = SimpleEphemeridesModel(),
            thermal_model    = MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography       = false,
            wind             = false,
        ),
        dynamics_model = DynamicsModel(
            constellation.spacecraft,
            (InverseSquaredJ2GravityModel(), laser_links...)
        ),
        guidance_model   = GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model = NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model    = ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time     = InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances = IntegrationTolerances(
            reltol_orbit = 1e-12,
            abstol_orbit = 1e-12,
            dt_max_orbit = opts.dt_max_s,
        ),
        solver_config = SolverConfig(solver_mode=:tsit5),
    )

    # 5. Callbacks
    impulse_tracker = LaserImpulseTracker()
    impulse_cb      = laser_impulse_callback(constellation, impulse_tracker, opts.mass_kg)
    scheduler_cb    = laser_link_scheduler_callback(constellation)

    # 6. Initial state for post-analysis
    u0  = SimulationEngine.build_initial_conditions(args)
    r0  = SVector{3,Float64}(u0.sc[1].pos)
    v0  = SVector{3,Float64}(u0.sc[1].vel)
    oe0 = rvtoorbitalelement(r0, v0, planet)

    # 7. Run
    result = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false") do
        run_simulation(
            args;
            isolate_state       = false,
            return_solution     = true,
            return_solver_metadata = true,
            extra_callbacks     = (impulse_cb, scheduler_cb),
        )
    end

    # 8. Post-analysis
    sol          = result.solution
    final_state  = sol.u[end].sc[1]
    rf  = SVector{3,Float64}(final_state.pos)
    vf  = SVector{3,Float64}(final_state.vel)
    oef = rvtoorbitalelement(rf, vf, planet)

    # Locate the feather file written by run_simulation
    feather_path = joinpath(results_dir, "simulation_results.feather")

    return (
        dv_r_mps        = sum(values(impulse_tracker.dv_R); init=0.0),
        dv_t_mps        = sum(values(impulse_tracker.dv_T); init=0.0),
        dv_n_mps        = sum(values(impulse_tracker.dv_N); init=0.0),
        active_steps    = impulse_tracker.active_link_steps,
        da_m            = oef[1] - oe0[1],
        de              = oef[2] - oe0[2],
        di_deg          = rad2deg(oef[3] - oe0[3]),
        target_period_s = target_period_s,
        retcode         = sol.retcode,
        results_dir     = results_dir,
        feather_path    = feather_path,
        impulse_tracker = impulse_tracker,
        n_helpers       = opts.helpers,
    )
end

# ── Entry point ────────────────────────────────────────────────────────────────

function main(argv=ARGS)
    opts = _parse_options(argv)
    println("Running ORACLE laser-link simulation...")
    println("  helpers=$(opts.helpers)  schedule=$(opts.schedule)  orbits=$(opts.orbits)")

    elapsed = @elapsed s = run_oracle_laser(opts)

    println("\n── Results ──────────────────────────────────────────")
    @printf("  dv_R  = %+.6e m/s\n", s.dv_r_mps)
    @printf("  dv_T  = %+.6e m/s\n", s.dv_t_mps)
    @printf("  dv_N  = %+.6e m/s\n", s.dv_n_mps)
    @printf("  Δa    = %+.3f m\n",   s.da_m)
    @printf("  Δe    = %+.3e\n",     s.de)
    @printf("  Δi    = %+.6f deg\n", s.di_deg)
    @printf("  active steps     : %d\n", s.active_steps)
    @printf("  retcode          : %s\n", s.retcode)
    @printf("  run time         : %.1f s\n", elapsed)
    println("  output → $(s.results_dir)")

    # Generate diagnostic plots
    println("\nGenerating plots...")
    plot_oracle_laser_results(s.feather_path, s.impulse_tracker, s.results_dir;
                              n_helpers=s.n_helpers)
    return s
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
