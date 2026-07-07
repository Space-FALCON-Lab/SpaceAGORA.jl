#!/usr/bin/env julia

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "examples", "common.jl"))

using CSV
using DataFrames
using LinearAlgebra
using Printf
using StaticArrays
using .SimulationModel

const DEFAULT_SUMMARY_CSV = joinpath(REPO_ROOT, "output", "oracle_case2_laser_links", "case2_laser_summary.csv")
const DEFAULT_TIMESERIES_CSV = joinpath(REPO_ROOT, "output", "oracle_case2_laser_links", "case2_laser_timeseries.csv")
const PAPER_HELPER_ALTITUDES_KM = (1150.0, 1050.0, 1000.0, 950.0, 850.0)
const PAPER_HELPER_INCLINATION_DELTAS_DEG = (0.0, 0.5, 1.0)

Base.@kwdef struct OracleCase2Options
    helpers::Int = 200
    helper_altitude_km::Float64 = 1050.0
    target_altitude_km::Float64 = 1000.0
    inclination_deg::Float64 = 0.0
    helper_inclination_delta_deg::Float64 = 0.0
    orbits::Float64 = 80.0
    schedule::Symbol = :naive_next_entering
    laser_range_km::Float64 = 200.0
    laser_power_w::Float64 = 10_000.0
    magnification::Float64 = 100.0
    beta::Float64 = 1.0
    eta::Float64 = 1.0
    mass_kg::Float64 = 227.0
    dt_max_s::Float64 = 10.0
    paper_grid::Bool = false
    output_csv::String = DEFAULT_SUMMARY_CSV
    timeseries_csv::String = DEFAULT_TIMESERIES_CSV
    append_output::Bool = false
    timeseries_points::Int = 1001
end

function _usage()
    return """
    Usage:
      julia --project=. ORACLE/run_case2_laser_links.jl [options]

    Options:
      --helpers N
      --helper-altitude-km KM
      --target-altitude-km KM
      --inclination-deg DEG
      --helper-inclination-delta-deg DEG
      --orbits N
      --schedule naive_next_entering|positive_along_track
      --laser-range-km KM
      --laser-power-w W
      --magnification B
      --beta VALUE
      --eta VALUE
      --mass-kg KG
      --dt-max-s SEC
      --paper-grid
      --output-csv PATH
      --timeseries-csv PATH
      --append-output
      --timeseries-points N
    """
end

function _parse_options(argv)::OracleCase2Options
    opts = Dict{Symbol, Any}()
    i = 1
    while i <= length(argv)
        arg = argv[i]
        if arg == "--help" || arg == "-h"
            println(_usage())
            exit(0)
        elseif arg == "--paper-grid"
            opts[:paper_grid] = true
            i += 1
            continue
        elseif arg == "--append-output"
            opts[:append_output] = true
            i += 1
            continue
        end
        startswith(arg, "--") || throw(ArgumentError("Unexpected argument '$arg'.\n$(_usage())"))
        i < length(argv) || throw(ArgumentError("Missing value for $arg."))
        key = Symbol(replace(arg[3:end], '-' => '_'))
        val = argv[i + 1]
        if key in (:helpers, :timeseries_points)
            opts[key] = parse(Int, val)
        elseif key == :schedule
            opts[key] = Symbol(val)
        elseif key in (:output_csv, :timeseries_csv)
            opts[key] = abspath(val)
        elseif key in (
            :helper_altitude_km, :target_altitude_km, :inclination_deg,
            :helper_inclination_delta_deg, :orbits,
            :laser_range_km, :laser_power_w, :magnification, :beta, :eta,
            :mass_kg, :dt_max_s
        )
            opts[key] = parse(Float64, val)
        else
            throw(ArgumentError("Unknown option $arg.\n$(_usage())"))
        end
        i += 2
    end
    return OracleCase2Options(; opts...)
end

function _validate_options(opts::OracleCase2Options)
    opts.helpers >= 1 || throw(ArgumentError("--helpers must be >= 1."))
    opts.helper_altitude_km > 0.0 || throw(ArgumentError("--helper-altitude-km must be positive."))
    opts.target_altitude_km > 0.0 || throw(ArgumentError("--target-altitude-km must be positive."))
    opts.orbits > 0.0 || throw(ArgumentError("--orbits must be positive."))
    opts.schedule in (:naive_next_entering, :positive_along_track) ||
        throw(ArgumentError("--schedule must be naive_next_entering or positive_along_track."))
    opts.laser_range_km >= 0.0 || throw(ArgumentError("--laser-range-km must be nonnegative."))
    opts.laser_power_w >= 0.0 || throw(ArgumentError("--laser-power-w must be nonnegative."))
    opts.magnification >= 0.0 || throw(ArgumentError("--magnification must be nonnegative."))
    opts.beta >= 0.0 || throw(ArgumentError("--beta must be nonnegative."))
    opts.eta >= 0.0 || throw(ArgumentError("--eta must be nonnegative."))
    opts.mass_kg > 0.0 || throw(ArgumentError("--mass-kg must be positive."))
    opts.dt_max_s > 0.0 || throw(ArgumentError("--dt-max-s must be positive."))
    opts.timeseries_points >= 2 || throw(ArgumentError("--timeseries-points must be >= 2."))
    return nothing
end

################################################################ DEFINE SPACECRAFT HERE ################################################################
function _spacecraft(id::Int, mass_kg::Float64, ic::InitialCondition)
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
########################################################################### END SPACECRAFT DEFINITION ###############################################################

@inline function _rtn_basis(r::SVector{3, Float64}, v::SVector{3, Float64})
    rhat = r / norm(r)
    nhat = cross(r, v)
    nhat = nhat / norm(nhat)
    that = cross(nhat, rhat)
    return rhat, that, nhat
end

function _rv_to_elements(r::SVector{3, Float64}, v::SVector{3, Float64}, mu::Float64)
    rmag = norm(r)
    v2 = dot(v, v)
    h = cross(r, v)
    hmag = norm(h)
    n = cross(SVector{3, Float64}(0.0, 0.0, 1.0), h)
    nmag = norm(n)
    evec = cross(v, h) / mu - r / rmag
    e = norm(evec)
    a = -mu / (2.0 * (0.5 * v2 - mu / rmag))
    inc = acos(clamp(h[3] / hmag, -1.0, 1.0))
    raan = nmag <= 1e-12 ? 0.0 : mod(atan(n[2], n[1]), 2pi)
    argp = if e <= 1e-12 || nmag <= 1e-12
        0.0
    else
        angle = acos(clamp(dot(n, evec) / (nmag * e), -1.0, 1.0))
        evec[3] < 0.0 ? 2pi - angle : angle
    end
    ta = if e <= 1e-12
        mod(atan(r[2], r[1]), 2pi)
    else
        angle = acos(clamp(dot(evec, r) / (e * rmag), -1.0, 1.0))
        dot(r, v) < 0.0 ? 2pi - angle : angle
    end
    return (a=a, e=e, i=inc, raan=raan, argp=argp, ta=ta)
end

function build_case2_config(opts::OracleCase2Options)
    _validate_options(opts)
    planet = make_no_gram_planet(:earth)
    target_radius_m = planet.Rp_e + opts.target_altitude_km * 1e3
    helper_radius_m = planet.Rp_e + opts.helper_altitude_km * 1e3

    spacecraft = SpacecraftModel[]
    push!(spacecraft, _spacecraft(1, opts.mass_kg, InitialCondition(
        target_radius_m, 0.0, opts.inclination_deg, 0.0, 0.0, 0.0
    )))
    helper_inclination_deg = opts.inclination_deg + opts.helper_inclination_delta_deg
    for helper in 1:opts.helpers
        nu_deg = 360.0 * (helper - 1) / opts.helpers
        push!(spacecraft, _spacecraft(helper + 1, opts.mass_kg, InitialCondition(
            helper_radius_m, 0.0, helper_inclination_deg, 0.0, 0.0, nu_deg
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
    mission_time_s = opts.orbits * target_period_s
    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
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

@inline function _target_rv_at(sol, t::Float64)
    u = if isapprox(t, Float64(sol.t[1]); atol=0.0, rtol=0.0)
        sol.u[1]
    elseif isapprox(t, Float64(sol.t[end]); atol=0.0, rtol=0.0)
        sol.u[end]
    else
        sol(t)
    end
    sc = u.sc[1]
    return SVector{3, Float64}(sc.pos), SVector{3, Float64}(sc.vel)
end

function _case_id(opts::OracleCase2Options)::String
    return @sprintf(
        "Nh%d_h%.0fkm_dinc%.1fdeg",
        opts.helpers,
        opts.helper_altitude_km,
        opts.helper_inclination_delta_deg,
    )
end

function _build_timeseries_dataframe(
    opts::OracleCase2Options,
    sol,
    baseline_sol,
    oe0,
    target_period_s::Float64,
    mu::Float64,
)::DataFrame
    tf = Float64(sol.t[end])
    times = collect(range(0.0, tf; length=opts.timeseries_points))
    rows = NamedTuple[]
    case_id = _case_id(opts)
    helper_inclination_deg = opts.inclination_deg + opts.helper_inclination_delta_deg
    for t in times
        r, v = _target_rv_at(sol, Float64(t))
        rb, vb = _target_rv_at(baseline_sol, Float64(t))
        rhat_b, that_b, nhat_b = _rtn_basis(rb, vb)
        dv = v - vb
        oe = _rv_to_elements(r, v, mu)
        oeb = _rv_to_elements(rb, vb, mu)
        push!(rows, (
            case_id=case_id,
            time_s=Float64(t),
            orbit=Float64(t) / target_period_s,
            helpers=opts.helpers,
            helper_altitude_km=opts.helper_altitude_km,
            target_altitude_km=opts.target_altitude_km,
            target_inclination_deg=opts.inclination_deg,
            helper_inclination_delta_deg=opts.helper_inclination_delta_deg,
            helper_inclination_deg=helper_inclination_deg,
            schedule=opts.schedule,
            laser_range_km=opts.laser_range_km,
            laser_power_w=opts.laser_power_w,
            magnification=opts.magnification,
            beta=opts.beta,
            eta=opts.eta,
            mass_kg=opts.mass_kg,
            dv_r_mps=dot(dv, rhat_b),
            dv_t_mps=dot(dv, that_b),
            dv_n_mps=dot(dv, nhat_b),
            da_m=oe.a - oe0.a,
            de=oe.e - oe0.e,
            di_deg=rad2deg(oe.i - oe0.i),
            draan_deg=rad2deg(oe.raan - oe0.raan),
            laser_da_m=oe.a - oeb.a,
            laser_de=oe.e - oeb.e,
            laser_di_deg=rad2deg(oe.i - oeb.i),
            laser_draan_deg=rad2deg(oe.raan - oeb.raan),
        ))
    end
    return DataFrame(rows)
end

function run_case2(opts::OracleCase2Options)
    args, laser_model, target_period_s = build_case2_config(opts)
    u0 = SimulationEngine.build_initial_conditions(args)
    r0 = SVector{3, Float64}(u0.sc[1].pos)
    v0 = SVector{3, Float64}(u0.sc[1].vel)
    oe0 = _rv_to_elements(r0, v0, args.environment_model.planet.μ)

    baseline_args = SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=args.simulation_settings,
        mission_configuration=args.mission_configuration,
        environment_model=args.environment_model,
        dynamics_model=DynamicsModel(SpacecraftModel[args.dynamics_model.spacecraft[1]], (InverseSquaredJ2GravityModel(),)),
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances,
        solver_config=args.solver_config,
    )
    baseline_sol = run_simulation(baseline_args; return_solution=true)
    baseline_final = baseline_sol.u[end].sc[1]
    rb = SVector{3, Float64}(baseline_final.pos)
    vb = SVector{3, Float64}(baseline_final.vel)
    rhat_b, that_b, nhat_b = _rtn_basis(rb, vb)

    callback = laser_link_scheduler_callback(laser_model)
    result = run_simulation(
        args;
        return_solution=true,
        return_solver_metadata=true,
        extra_callbacks=(callback,),
    )
    sol = result.solution
    final_state = sol.u[end].sc[1]
    rf = SVector{3, Float64}(final_state.pos)
    vf = SVector{3, Float64}(final_state.vel)
    dv = vf - vb
    oef = _rv_to_elements(rf, vf, args.environment_model.planet.μ)
    oeb = _rv_to_elements(rb, vb, args.environment_model.planet.μ)
    model_after = result.solution.prob.p.args.dynamics_model.dynamic_effectors[2]

    summary = (
        case_id=_case_id(opts),
        helpers=opts.helpers,
        helper_altitude_km=opts.helper_altitude_km,
        schedule=opts.schedule,
        target_period_s=target_period_s,
        active_helper=model_after.active_helper_idx,
        activations=model_after.link_activation_count,
        active_steps=model_after.active_link_step_count,
        target_altitude_km=opts.target_altitude_km,
        inclination_deg=opts.inclination_deg,
        target_inclination_deg=opts.inclination_deg,
        helper_inclination_delta_deg=opts.helper_inclination_delta_deg,
        helper_inclination_deg=opts.inclination_deg + opts.helper_inclination_delta_deg,
        orbits=opts.orbits,
        laser_range_km=opts.laser_range_km,
        laser_power_w=opts.laser_power_w,
        magnification=opts.magnification,
        beta=opts.beta,
        eta=opts.eta,
        mass_kg=opts.mass_kg,
        dt_max_s=opts.dt_max_s,
        dv_r_mps=dot(dv, rhat_b),
        dv_t_mps=dot(dv, that_b),
        dv_n_mps=dot(dv, nhat_b),
        da_m=oef.a - oe0.a,
        de=oef.e - oe0.e,
        di_deg=rad2deg(oef.i - oe0.i),
        draan_deg=rad2deg(oef.raan - oe0.raan),
        laser_da_m=oef.a - oeb.a,
        laser_de=oef.e - oeb.e,
        retcode=sol.retcode,
        solver=result.solver_trace[end].solver,
    )
    timeseries = _build_timeseries_dataframe(
        opts,
        sol,
        baseline_sol,
        oe0,
        target_period_s,
        args.environment_model.planet.μ,
    )
    return (summary=summary, timeseries=timeseries)
end

function _summary_dataframe(rows)::DataFrame
    return DataFrame(rows)
end

function _write_summary_rows!(path::String, rows; append::Bool=false)
    isempty(path) && return ""
    mkpath(dirname(path))
    df = _summary_dataframe(rows)
    do_append = append && isfile(path)
    CSV.write(path, df; append=do_append, header=!do_append)
    return path
end

function _write_timeseries!(path::String, df::DataFrame; append::Bool=false)
    isempty(path) && return ""
    mkpath(dirname(path))
    do_append = append && isfile(path)
    CSV.write(path, df; append=do_append, header=!do_append)
    return path
end

function _print_summary(s)
    println("ORACLE Case 2 laser-link run")
    println("  helpers              : $(s.helpers)")
    println("  helper altitude       : $(s.helper_altitude_km) km")
    println("  schedule              : $(s.schedule)")
    println("  solver                : $(s.solver)")
    println("  retcode               : $(s.retcode)")
    @printf("  mission time          : %.2f s (%.4f orbits)\n", s.orbits * s.target_period_s, s.orbits)
    println("  link activations      : $(s.activations)")
    println("  accepted active steps : $(s.active_steps)")
    @printf("  dV_RTN [m/s]          : R=% .6e  T=% .6e  N=% .6e\n", s.dv_r_mps, s.dv_t_mps, s.dv_n_mps)
    @printf("  da [m]                : % .6e\n", s.da_m)
    @printf("  de                    : % .6e\n", s.de)
    @printf("  di [deg]              : % .6e\n", s.di_deg)
    @printf("  dRAAN [deg]           : % .6e\n", s.draan_deg)
    @printf("  laser-only da [m]     : % .6e\n", s.laser_da_m)
    @printf("  laser-only de         : % .6e\n", s.laser_de)
end

function main(argv=ARGS)
    opts = _parse_options(argv)
    if opts.paper_grid
        opts.append_output || (isfile(opts.output_csv) && rm(opts.output_csv))
        opts.append_output || (isfile(opts.timeseries_csv) && rm(opts.timeseries_csv))
        for helper_alt in PAPER_HELPER_ALTITUDES_KM
            for helper_inclination_delta in PAPER_HELPER_INCLINATION_DELTAS_DEG
                result = run_case2(OracleCase2Options(
                    helpers=opts.helpers,
                    helper_altitude_km=helper_alt,
                    target_altitude_km=opts.target_altitude_km,
                    inclination_deg=opts.inclination_deg,
                    helper_inclination_delta_deg=helper_inclination_delta,
                    orbits=opts.orbits,
                    schedule=opts.schedule,
                    laser_range_km=opts.laser_range_km,
                    laser_power_w=opts.laser_power_w,
                    magnification=opts.magnification,
                    beta=opts.beta,
                    eta=opts.eta,
                    mass_kg=opts.mass_kg,
                    dt_max_s=opts.dt_max_s,
                    output_csv=opts.output_csv,
                    timeseries_csv=opts.timeseries_csv,
                    append_output=true,
                    timeseries_points=opts.timeseries_points,
                ))
                s = result.summary
                _write_summary_rows!(opts.output_csv, [s]; append=true)
                _write_timeseries!(opts.timeseries_csv, result.timeseries; append=true)
                @printf(
                    "helpers=%d helper_alt_km=%.1f helper_dinc_deg=%.1f dv_R=%.6e dv_T=%.6e dv_N=%.6e da_m=%.6e activations=%d\n",
                    s.helpers, s.helper_altitude_km, s.helper_inclination_delta_deg,
                    s.dv_r_mps, s.dv_t_mps, s.dv_n_mps, s.laser_da_m, s.activations
                )
            end
        end
        println("Summary CSV: $(opts.output_csv)")
        println("Time-series CSV: $(opts.timeseries_csv)")
    else
        result = run_case2(opts)
        s = result.summary
        _print_summary(s)
        out = _write_summary_rows!(opts.output_csv, [s]; append=opts.append_output)
        isempty(out) || println("Summary CSV: $out")
        ts_out = _write_timeseries!(opts.timeseries_csv, result.timeseries; append=opts.append_output)
        isempty(ts_out) || println("Time-series CSV: $ts_out")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
