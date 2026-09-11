# Backfill density_kg_m3, dynamic_pressure_pa, altitude_m, and in_atmosphere
# columns in existing aerobraking perturbation MC outputs without rerunning dynamics.
#
# Usage:
#   julia --project=. benchmarks/studies/aerobraking_perturbation_mc/backfill_density_history.jl RUN_DIR [--force] [--dry-run]

using Arrow
using DataFrames
using LinearAlgebra
using Printf
using StaticArrays

include(joinpath(@__DIR__, "study.jl"))
using .AerobrakingPerturbationMC

const MC = AerobrakingPerturbationMC
const DENSITY_COLUMNS = (:altitude_m, :density_kg_m3, :dynamic_pressure_pa, :in_atmosphere)

function _backfill_args(info, planet_cache, density_cache)
    planet = get!(planet_cache, info.planet) do
        MC._planet(info.planet)
    end
    density_key = (info.planet, info.dynamics_case, info.density_scale)
    density_model = get!(density_cache, density_key) do
        MC._density_model(info.planet, info.dynamics_case, info.density_scale)
    end
    rp_alt_m = MC._periapsis_altitude_m(info.planet, info.periapsis_regime)
    environment_model = (
        planet=planet,
        ephemerides_model=MC.SM.SimpleEphemeridesModel(),
        density_model=density_model,
        EI=max(220.0, rp_alt_m / 1e3 + 80.0),
    )
    mission_configuration = (keplerian=false,)
    return (
        environment_model=environment_model,
        mission_configuration=mission_configuration,
        initial_time=MC.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
    )
end

function _parse_case_dir(dirname::String)
    m = match(
        r"^case_\d+_(mars|earth|venus|titan)_(shallow|nominal|deep)_apo_(\d+)km_(.+)_density_(none|low|nominal|high)_ms(\d+p\d+)_inc(\d+)_aop(\d+)$",
        dirname,
    )
    m === nothing && return nothing
    density_case = Symbol(m[5])
    density_scale = density_case == :high ? 1.1 :
                    density_case == :low ? 0.9 :
                    density_case == :nominal ? 1.0 : NaN
    return (
        planet=Symbol(m[1]),
        periapsis_regime=Symbol(m[2]),
        apoapsis_alt_km=parse(Float64, m[3]),
        dynamics_case=Symbol(m[4]),
        density_case=density_case,
        density_scale=density_scale,
        spacecraft_mass_scale=parse(Float64, replace(m[6], "p" => ".")),
        inclination_deg=parse(Float64, m[7]),
        argp_deg=parse(Float64, m[8]),
        norbits=1,
        aero_solver=MC.DEFAULT_AERO_SOLVER,
        aero_auto_maxiters=MC.DEFAULT_AERO_AUTO_MAXITERS,
        aero_stiff_maxiters=MC.DEFAULT_AERO_STIFF_MAXITERS,
        aero_auto_stiff_switch_max=MC.DEFAULT_AERO_AUTO_STIFF_SWITCH_MAX,
        dt_max_orbit_s=MC.DEFAULT_DT_MAX_ORBIT_S,
        dt_max_atmosphere_s=MC.DEFAULT_DT_MAX_ATMOSPHERE_S,
    )
end

function _trajectory_paths(case_dir::String)
    paths = String[]
    for name in ("trajectory_with_active_force.feather", "simulation_results.feather")
        path = joinpath(case_dir, name)
        isfile(path) && push!(paths, path)
    end
    orbit_paths = filter(sort(filter(isfile, readdir(case_dir; join=true)))) do path
        occursin(r"trajectory_with_active_force_orbit_\d+\.feather$", basename(path))
    end
    append!(paths, orbit_paths)
    return unique(paths)
end

function _density_needs_backfill(df::DataFrame)::Bool
    all(name -> name in propertynames(df), DENSITY_COLUMNS) || return true
    for name in (:density_kg_m3, :dynamic_pressure_pa)
        vals = Float64.(df[!, name])
        any(isfinite, vals) && return false
    end
    return true
end

function _state_from_row(row)
    pos = SVector{3, Float64}(
        Float64(row.sc1_pos_1),
        Float64(row.sc1_pos_2),
        Float64(row.sc1_pos_3),
    )
    vel = SVector{3, Float64}(
        Float64(row.sc1_vel_1),
        Float64(row.sc1_vel_2),
        Float64(row.sc1_vel_3),
    )
    mass = Float64(row.sc1_mass)
    return MC.SM.StateSample(pos, vel, mass), pos, vel
end

function _sample_density_history(args, df::DataFrame)
    n = nrow(df)
    altitude = Vector{Float64}(undef, n)
    density = Vector{Float64}(undef, n)
    qdyn = Vector{Float64}(undef, n)
    in_atmosphere = Vector{Bool}(undef, n)
    density_model = args.environment_model.density_model
    p_stub = (; args=args)
    EI_m = args.environment_model.EI * 1e3

    for (i, row) in enumerate(eachrow(df))
        t = Float64(row.time)
        state, _, _ = _state_from_row(row)
        frame = MC._planet_frame_sample(args, state, t)
        altitude[i] = frame.alt_m
        in_atmosphere[i] = frame.alt_m <= EI_m

        if density_model isa MC.NoAtmosphereModel
            density[i] = 0.0
            qdyn[i] = 0.0
        else
            rho, _, wind = MC.getDensity(
                density_model,
                Float64(frame.alt_m),
                Float64(frame.lat_rad),
                Float64(frame.lon_rad),
                t,
                true,
                p_stub,
            )
            density[i] = rho
            qdyn[i] = 0.5 * rho * norm(frame.vel_pp - wind)^2
        end
    end

    return (; altitude, density, qdyn, in_atmosphere)
end

function _write_density_columns!(path::String, args; dry_run::Bool=false, force::Bool=false)
    df = DataFrame(Arrow.Table(path); copycols=true)
    required = (:time, :sc1_pos_1, :sc1_pos_2, :sc1_pos_3, :sc1_vel_1, :sc1_vel_2, :sc1_vel_3, :sc1_mass)
    all(name -> name in propertynames(df), required) || return nothing
    if !force && !_density_needs_backfill(df)
        return nothing
    end

    hist = _sample_density_history(args, df)
    df[!, :altitude_m] = hist.altitude
    df[!, :density_kg_m3] = hist.density
    df[!, :dynamic_pressure_pa] = hist.qdyn
    df[!, :in_atmosphere] = hist.in_atmosphere

    dry_run || MC._write_feather(path, df)
    finite_density = filter(isfinite, hist.density)
    finite_qdyn = filter(isfinite, hist.qdyn)
    return (
        path=path,
        rows=nrow(df),
        peak_density=isempty(finite_density) ? NaN : maximum(finite_density),
        max_dynamic_pressure=isempty(finite_qdyn) ? NaN : maximum(finite_qdyn),
        dataframe=df,
    )
end

function _copy_density_columns!(
    source_path::String,
    target_path::String;
    source_df=nothing,
    dry_run::Bool=false,
    force::Bool=false,
)
    source_path == target_path && return nothing
    source = source_df === nothing ? DataFrame(Arrow.Table(source_path)) : source_df
    target = DataFrame(Arrow.Table(target_path))
    all(name -> name in propertynames(source), (:time, DENSITY_COLUMNS...)) || return nothing
    :time in propertynames(target) || return nothing
    if !force && !_density_needs_backfill(target)
        return nothing
    end

    source_index = Dict(Float64(source.time[i]) => i for i in 1:nrow(source))
    idxs = Vector{Int}(undef, nrow(target))
    for i in 1:nrow(target)
        t = Float64(target.time[i])
        idx = get(source_index, t, 0)
        idx == 0 && return nothing
        idxs[i] = idx
    end
    for name in DENSITY_COLUMNS
        target[!, name] = source[idxs, name]
    end
    dry_run || MC._write_feather(target_path, target)

    finite_density = filter(isfinite, Float64.(target.density_kg_m3))
    finite_qdyn = filter(isfinite, Float64.(target.dynamic_pressure_pa))
    return (
        path=target_path,
        rows=nrow(target),
        peak_density=isempty(finite_density) ? NaN : maximum(finite_density),
        max_dynamic_pressure=isempty(finite_qdyn) ? NaN : maximum(finite_qdyn),
    )
end

function _density_metrics_from_trajectory(path::String)
    df = DataFrame(Arrow.Table(path); copycols=true)
    all(name -> name in propertynames(df), (:time, :density_kg_m3, :dynamic_pressure_pa, :in_atmosphere)) || return nothing
    nrow(df) == 0 && return nothing
    sort!(df, :time)
    rho = Float64.(df.density_kg_m3)
    q = Float64.(df.dynamic_pressure_pa)
    inside = Bool.(df.in_atmosphere)
    t = Float64.(df.time)
    finite_rho = filter(isfinite, rho)
    finite_q = filter(isfinite, q)
    (isempty(finite_rho) || isempty(finite_q)) && return nothing

    integrated_density = 0.0
    integrated_dynamic_pressure = 0.0
    time_below_interface_s = 0.0
    for i in 2:length(t)
        dt = max(0.0, t[i] - t[i - 1])
        if isfinite(rho[i]) && isfinite(rho[i - 1])
            integrated_density += 0.5 * (rho[i] + rho[i - 1]) * dt
        end
        if isfinite(q[i]) && isfinite(q[i - 1])
            integrated_dynamic_pressure += 0.5 * (q[i] + q[i - 1]) * dt
        end
        if inside[i] || inside[i - 1]
            time_below_interface_s += dt
        end
    end

    return (
        peak_density=maximum(finite_rho),
        integrated_density=integrated_density,
        max_dynamic_pressure=maximum(finite_q),
        integrated_dynamic_pressure=integrated_dynamic_pressure,
        time_below_interface_s=time_below_interface_s,
    )
end

function _update_result_row!(case_dir::String, metrics; dry_run::Bool=false)
    path = joinpath(case_dir, "result_row.feather")
    isfile(path) || return false
    df = DataFrame(Arrow.Table(path); copycols=true)
    nrow(df) == 0 && return false
    for key in keys(metrics)
        df[!, key] = fill(getfield(metrics, key), nrow(df))
    end
    dry_run || MC._write_feather(path, df)
    return true
end

function _update_run_results!(run_dir::String, case_metrics::Dict{String, NamedTuple}; dry_run::Bool=false)
    path = joinpath(run_dir, "results.feather")
    isfile(path) || return 0
    df = DataFrame(Arrow.Table(path); copycols=true)
    :case_dir in propertynames(df) || return 0
    updated = 0
    for i in 1:nrow(df)
        case_dir = abspath(String(df.case_dir[i]))
        metrics = get(case_metrics, case_dir, nothing)
        metrics === nothing && continue
        for key in keys(metrics)
            key in propertynames(df) || continue
            df[i, key] = getfield(metrics, key)
        end
        updated += 1
    end
    if updated > 0 && !dry_run
        MC._write_feather(path, df)
        delta_df = MC.paired_deltas(df)
        aggregate_df = MC.aggregate_deltas(delta_df)
        MC._write_feather(joinpath(run_dir, "paired_deltas.feather"), delta_df)
        MC._write_feather(joinpath(run_dir, "aggregates.feather"), aggregate_df)
    end
    return updated
end

function backfill_run!(run_dir::String; dry_run::Bool=false, force::Bool=false)
    case_dirs = sort(filter(isdir, readdir(run_dir; join=true)))
    updated_files = 0
    updated_cases = 0
    case_metrics = Dict{String, NamedTuple}()
    planet_cache = Dict{Symbol, Any}()
    density_cache = Dict{Tuple{Symbol, Symbol, Float64}, Any}()
    processed_cases = 0
    progress_interval = 5000
    for case_dir in case_dirs
        info = _parse_case_dir(basename(case_dir))
        info === nothing && continue
        paths = _trajectory_paths(case_dir)
        isempty(paths) && continue
        processed_cases += 1
        args = _backfill_args(info, planet_cache, density_cache)

        case_updated = false
        primary_path = joinpath(case_dir, "trajectory_with_active_force.feather")
        sample_paths = isfile(primary_path) ? [primary_path] : paths
        copy_paths = isfile(primary_path) ? [path for path in paths if path != primary_path] : String[]
        primary_df = nothing

        for path in sample_paths
            result = _write_density_columns!(path, args; dry_run=dry_run, force=force)
            result === nothing && continue
            updated_files += 1
            case_updated = true
            path == primary_path && (primary_df = result.dataframe)
        end

        for path in copy_paths
            result = _copy_density_columns!(primary_path, path; source_df=primary_df, dry_run=dry_run, force=force)
            result === nothing && continue
            updated_files += 1
            case_updated = true
        end

        combined_path = joinpath(case_dir, "trajectory_with_active_force.feather")
        if isfile(combined_path) && (case_updated || force)
            metrics = _density_metrics_from_trajectory(combined_path)
            if metrics !== nothing
                case_metrics[abspath(case_dir)] = metrics
                _update_result_row!(case_dir, metrics; dry_run=dry_run)
            end
        end
        case_updated && (updated_cases += 1)
        if processed_cases % progress_interval == 0
            @printf("  processed=%d updated_cases=%d updated_files=%d\n",
                processed_cases, updated_cases, updated_files)
            flush(stdout)
        end
    end

    updated_result_rows = _update_run_results!(run_dir, case_metrics; dry_run=dry_run)
    @printf("Backfill complete: cases=%d files=%d run_result_rows=%d%s\n",
        updated_cases, updated_files, updated_result_rows, dry_run ? " (dry run)" : "")
    return (; updated_cases, updated_files, updated_result_rows)
end

function _parse_args(args)
    isempty(args) && error("Usage: julia --project=. $(PROGRAM_FILE) RUN_DIR [--force] [--dry-run]")
    run_dir = ""
    force = false
    dry_run = false
    for arg in args
        if arg == "--force"
            force = true
        elseif arg == "--dry-run"
            dry_run = true
        elseif startswith(arg, "-")
            error("Unknown option: $arg")
        else
            isempty(run_dir) || error("Only one RUN_DIR is supported.")
            run_dir = arg
        end
    end
    isempty(run_dir) && error("RUN_DIR is required.")
    isdir(run_dir) || error("RUN_DIR does not exist: $run_dir")
    return (; run_dir=abspath(run_dir), force, dry_run)
end

function main(args=ARGS)
    opts = _parse_args(args)
    println("Run directory : $(opts.run_dir)")
    println("Mode          : $(opts.dry_run ? "dry-run" : "rewrite")")
    println("Force         : $(opts.force)")
    backfill_run!(opts.run_dir; dry_run=opts.dry_run, force=opts.force)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
