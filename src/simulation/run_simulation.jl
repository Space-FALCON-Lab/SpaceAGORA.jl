include("../utils/Reference_system.jl")
using LinearAlgebra
using StaticArrays
using LoopVectorization
using ComponentArrays
using DifferentialEquations
using CSV
using DataFrames
using Polyester
using Arrow
using Dates
using Serialization
using SHA
using TOML

const _normalize_warning_emitted = Ref(false)
const RESULTS_BUNDLE_SCHEMA_VERSION = "1"
const CHECKPOINT_SCHEMA_VERSION = "1"

@inline _typed_normalize_warning_enabled() = get(ENV, "SPACEAGORA_WARN_NORMALIZE", "1") == "1"
@inline _typed_allow_legacy_normalize() = get(ENV, "SPACEAGORA_ALLOW_TYPED_NORMALIZE", "0") == "1"
@inline _typed_save_bundle_enabled() = get(ENV, "SPACEAGORA_SAVE_BUNDLE", "1") == "1"
@inline _typed_checkpoint_enabled(args) = args.simulation_settings.checkpoint_enabled || args.simulation_settings.resume_from_checkpoint

function _validate_orientation_inertia!(args)
    if !args.mission_configuration.orientation_sim
        return nothing
    end
    for (i, sc) in enumerate(args.dynamics_model.spacecraft)
        inertia_tensor = sc.inertia_tensor
        inertia_matrix = Matrix(inertia_tensor)
        if !all(isfinite, inertia_matrix)
            throw(ArgumentError("Spacecraft index $i has non-finite inertia tensor entries, required for orientation_sim=true."))
        end
        if !issymmetric(inertia_matrix)
            throw(ArgumentError("Spacecraft index $i has non-symmetric inertia tensor, required for orientation_sim=true."))
        end
        if !isposdef(Symmetric(inertia_matrix))
            throw(ArgumentError("Spacecraft index $i has non-positive-definite inertia tensor, required for orientation_sim=true."))
        end
    end
    return nothing
end

function _validate_thermal_model_support!(args)
    thermal_model = args.environment_model.thermal_model
    if !hasmethod(SimulationModel.getHeatRate, Tuple{typeof(thermal_model), Float64, Float64, Float64, Float64, Float64})
        throw(ArgumentError(
            "Thermal model $(nameof(typeof(thermal_model))) must implement " *
            "getHeatRate(model, S, T, ρ, v, α)::Float64."
        ))
    end
    return nothing
end

function _initialize_heat_rate_buffers!(p)
    n_sats = length(p.args.dynamics_model.spacecraft)
    if length(p.shared_buffers.heat_rates) != n_sats
        resize!(p.shared_buffers.heat_rates, n_sats)
    end
    @inbounds for i in 1:n_sats
        n_links = length(p.args.dynamics_model.spacecraft[i].links)
        rates = p.shared_buffers.heat_rates[i]
        if !(rates isa Vector{Float64})
            rates = Float64[]
            p.shared_buffers.heat_rates[i] = rates
        end
        if length(rates) != n_links
            resize!(rates, n_links)
        end
        fill!(rates, 0.0)
    end
    return nothing
end

@inline function _gram_per_sat_instances_enabled()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_GRAM_PER_SAT_INSTANCES", "0")))
    if raw in ("1", "true", "yes", "on")
        return true
    elseif raw in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError(
        "Invalid SPACEAGORA_GRAM_PER_SAT_INSTANCES='$raw'. Use one of: 1/0, true/false, yes/no, on/off."
    ))
end

@inline function _parse_positive_float_env(name::String, default::Float64)::Float64
    raw = strip(get(ENV, name, string(default)))
    parsed = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a floating-point value, got '$raw'"))
    end
    parsed > 0.0 || throw(ArgumentError("$name must be > 0.0, got $parsed"))
    return parsed
end

@inline function _effector_parallel_mode()::Symbol
    return SimulationModel.ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_EFFECTOR_PARALLEL"; default="auto")
end

@inline function _effector_thread_threshold()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_THREAD_THRESHOLD", 2)
end

@inline function _effector_max_threads()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_MAX_THREADS", 4)
end

@inline function _effector_outer_parallel_hint()::Bool
    return SimulationModel.ParallelPolicy.outer_parallel_active()
end

@inline function _effector_allow_with_outer()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_EFFECTOR_PARALLEL_ALLOW_WITH_OUTER", false)
end

@inline function _effector_heavy_only()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY", true)
end

@inline function _effector_long_mission_threshold_s()::Float64
    return _parse_positive_float_env("SPACEAGORA_EFFECTOR_LONG_MISSION_THRESHOLD_S", 5400.0)
end

@inline function _effector_long_orbit_threshold()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_LONG_ORBIT_THRESHOLD", 8)
end

@inline function _dynamic_effector_threadsafe(::Any)::Bool
    return false
end
@inline _dynamic_effector_threadsafe(::SimulationModel.InverseSquaredGravityModel)::Bool = true
@inline _dynamic_effector_threadsafe(::SimulationModel.InverseSquaredJ2GravityModel)::Bool = true
@inline _dynamic_effector_threadsafe(::SimulationModel.NBodyGravityModel)::Bool = true
@inline _dynamic_effector_threadsafe(::SimulationModel.GravitationalHarmonicsModel)::Bool = true
@inline _dynamic_effector_threadsafe(::SimulationModel.SolarRadiationPressureModel)::Bool = true
@inline _dynamic_effector_threadsafe(::SimulationModel.AerodynamicCoefficientfM)::Bool = true

@inline function _dynamic_effectors_parallel_supported(dynamic_effectors::Tuple)::Bool
    aero_fm_count = 0
    @inbounds for effector in dynamic_effectors
        if effector isa SimulationModel.AerodynamicCoefficientfM
            aero_fm_count += 1
        end
        _dynamic_effector_threadsafe(effector) || return false
    end
    return aero_fm_count <= 1
end

@inline function _mission_is_long_for_effector_threads(args)::Bool
    mission_cfg = args.mission_configuration
    if mission_cfg.mission_type == SimulationModel.MissionOrbits
        return mission_cfg.number_of_orbits >= _effector_long_orbit_threshold()
    end
    return mission_cfg.mission_time >= _effector_long_mission_threshold_s()
end

@inline function _dynamic_effector_thread_decision(args::SimulationConfiguration, dynamic_effectors::Tuple, num_sats::Int)
    mode = _effector_parallel_mode()
    n_effectors = length(dynamic_effectors)
    if n_effectors <= 1
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end
    if !_dynamic_effectors_parallel_supported(dynamic_effectors)
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end

    single_sat = num_sats == 1
    single_body = single_sat && length(args.dynamics_model.spacecraft[1].links) <= 1
    heavy_work = _mission_is_long_for_effector_threads(args)
    if mode == :auto && (!single_sat || !single_body || !heavy_work)
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end

    policy = SimulationModel.ParallelPolicy.thread_policy_decision(
        n_effectors;
        mode=mode,
        threshold=_effector_thread_threshold(),
        heavy_work=heavy_work,
        heavy_only=_effector_heavy_only(),
        outer_active=_effector_outer_parallel_hint(),
        allow_with_outer=_effector_allow_with_outer(),
        source=:dynamic_effectors
    )
    allotment = min(policy.allotment, _effector_max_threads())
    use_threads = policy.use_threads && allotment > 1
    return (use_threads=use_threads, allotment=use_threads ? allotment : 1, mode=mode, policy_applied=true)
end

function _initialize_density_model_instances!(p)
    instances = p.shared_buffers.density_models
    empty!(instances)

    density_model = p.args.environment_model.density_model
    if !_gram_per_sat_instances_enabled() || !(
        density_model isa SimulationModel.EnvironmentModels.GRAMAtmosphereModel ||
        density_model isa SimulationModel.EnvironmentModels.GRAMAtmosphereModelSurrogate
    )
        return nothing
    end

    n_sats = length(p.args.dynamics_model.spacecraft)
    sizehint!(instances, n_sats)
    @inbounds for _ in 1:n_sats
        # One GRAM handle per satellite avoids sharing mutable native model state.
        push!(instances, deepcopy(density_model))
    end
    return nothing
end

function _initialize_density_cache_buffers!(p)
    n_sats = length(p.args.dynamics_model.spacecraft)
    caches = p.shared_buffers.gram_density_cache
    if length(caches) != n_sats
        resize!(caches, n_sats)
    end
    fill!(caches, nothing)
    return nothing
end

@inline function _resolve_component_tolerance(component_tol::Float64, fallback_tol::Float64, name::String)::Float64
    if component_tol < 0.0
        throw(ArgumentError("$name must be >= 0.0, got $component_tol"))
    end
    return component_tol == 0.0 ? fallback_tol : component_tol
end

@inline function _requires_componentwise_tolerances(args)::Bool
    tol = args.integration_tolerances
    return args.mission_configuration.orientation_sim ||
           tol.reltol_mass != 0.0 || tol.abstol_mass != 0.0 ||
           tol.reltol_heat_load != 0.0 || tol.abstol_heat_load != 0.0
end

function _build_solver_tolerances(u_state::ComponentVector, args)
    tol = args.integration_tolerances
    if !_requires_componentwise_tolerances(args)
        return tol.reltol_orbit, tol.abstol_orbit
    end

    reltol_mass = _resolve_component_tolerance(tol.reltol_mass, tol.reltol_orbit, "reltol_mass")
    abstol_mass = _resolve_component_tolerance(tol.abstol_mass, tol.abstol_orbit, "abstol_mass")
    reltol_heat = _resolve_component_tolerance(tol.reltol_heat_load, tol.reltol_orbit, "reltol_heat_load")
    abstol_heat = _resolve_component_tolerance(tol.abstol_heat_load, tol.abstol_orbit, "abstol_heat_load")
    reltol_ω = _resolve_component_tolerance(tol.reltol_angular_rate, tol.reltol_orbit, "reltol_angular_rate")
    abstol_ω = _resolve_component_tolerance(tol.abstol_angular_rate, tol.abstol_orbit, "abstol_angular_rate")

    reltol_state = copy(u_state)
    abstol_state = copy(u_state)
    reltol_state .= tol.reltol_orbit
    abstol_state .= tol.abstol_orbit
    @inbounds for i in eachindex(reltol_state.sc)
        reltol_state.sc[i].mass = reltol_mass
        abstol_state.sc[i].mass = abstol_mass
        reltol_state.sc[i].heat_loads .= reltol_heat
        abstol_state.sc[i].heat_loads .= abstol_heat
    end

    if args.mission_configuration.orientation_sim
        @inbounds for i in eachindex(reltol_state.sc)
            reltol_state.sc[i].ω .= reltol_ω
            abstol_state.sc[i].ω .= abstol_ω
        end
    end

    if args.mission_configuration.orientation_sim
        @inbounds for i in eachindex(reltol_state.sc)
        reltol_state.sc[i].q .= tol.reltol_quaternion
        abstol_state.sc[i].q .= tol.abstol_quaternion
        end
    end
    return reltol_state, abstol_state
end

@inline function _solver_policy_mode()::Symbol
    mode = lowercase(strip(get(ENV, "SPACEAGORA_SOLVER_MODE", "tsit5")))
    if mode in ("tsit5", "default")
        return :tsit5
    elseif mode in ("auto_stiff", "auto-stiff", "autostiff", "auto")
        return :auto_stiff
    elseif mode in ("rodas5p", "rodas", "stiff")
        return :rodas5p
    end
    throw(ArgumentError(
        "Unsupported SPACEAGORA_SOLVER_MODE='$mode'. Use one of: tsit5, auto_stiff, rodas5p."
    ))
end

@inline function _retcode_is_stiff_symptom(retcode)::Bool
    rc = string(retcode)
    return rc in ("Unstable", "DtLessThanMin", "MaxIters", "InitialFailure")
end

@inline function _solver_maxiters()::Union{Nothing, Int}
    raw = strip(get(ENV, "SPACEAGORA_SOLVER_MAXITERS", ""))
    isempty(raw) && return nothing
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_SOLVER_MAXITERS must be an integer, got '$raw'."))
    end
    parsed > 0 || throw(ArgumentError("SPACEAGORA_SOLVER_MAXITERS must be > 0, got $parsed."))
    return parsed
end

@inline function _solve_with_explicit_solver(prob, args, alg, reltol_tol, abstol_tol)
    maxiters = _solver_maxiters()
    if maxiters === nothing
        return solve(
            prob,
            alg;
            reltol=reltol_tol,
            abstol=abstol_tol,
            dtmax=args.integration_tolerances.dt_max_orbit
        )
    end
    return solve(
        prob,
        alg;
        reltol=reltol_tol,
        abstol=abstol_tol,
        dtmax=args.integration_tolerances.dt_max_orbit,
        maxiters=maxiters
    )
end

function _solve_with_solver_policy(prob, args, reltol_tol, abstol_tol)
    mode = _solver_policy_mode()
    if mode == :rodas5p
        sol = _solve_with_explicit_solver(prob, args, Rodas5P(autodiff=false), reltol_tol, abstol_tol)
        return sol, (
            solver="Rodas5P",
            initial_solver="Rodas5P",
            fallback_used=false,
            trigger_retcode=missing
        )
    end

    tsit_sol = _solve_with_explicit_solver(prob, args, Tsit5(), reltol_tol, abstol_tol)
    if mode == :auto_stiff &&
       !SciMLBase.successful_retcode(tsit_sol.retcode) &&
       _retcode_is_stiff_symptom(tsit_sol.retcode)
        stiff_sol = _solve_with_explicit_solver(prob, args, Rodas5P(autodiff=false), reltol_tol, abstol_tol)
        return stiff_sol, (
            solver="Rodas5P",
            initial_solver="Tsit5",
            fallback_used=true,
            trigger_retcode=string(tsit_sol.retcode)
        )
    end
    return tsit_sol, (
        solver="Tsit5",
        initial_solver="Tsit5",
        fallback_used=false,
        trigger_retcode=missing
    )
end

function _warn_legacy_normalize_flag!(args)
    if !args.simulation_settings.normalize || !_typed_normalize_warning_enabled()
        return nothing
    end
    if _normalize_warning_emitted[]
        return nothing
    end
    _normalize_warning_emitted[] = true
    @warn "SimulationSettings.normalize=true is legacy-only in typed run_simulation; propagation is always SI-native (m, s, kg). Set normalize=false to silence this warning."
    return nothing
end

function _enforce_typed_normalize_policy!(args)
    if !args.simulation_settings.normalize
        return nothing
    end
    if _typed_allow_legacy_normalize()
        _warn_legacy_normalize_flag!(args)
        return nothing
    end
    throw(ArgumentError(
        "SimulationSettings.normalize=true is unsupported in typed run_simulation. " *
        "Typed propagation is SI-native (m, s, kg). Set normalize=false, or set " *
        "SPACEAGORA_ALLOW_TYPED_NORMALIZE=1 only for legacy transition checks."
    ))
end

function _debug_print_nan_parameter_paths!(x, path::AbstractString="p")
    if x isa Number
        if isnan(x)
            println("NaN found in parameter: $path")
        end
        return nothing
    end

    if x isa Base.RefValue{<:Number}
        xv = x[]
        if isnan(xv)
            println("NaN found in parameter: $path[]")
        end
        return nothing
    end

    if x isa AbstractArray{<:Number}
        for (idx, xv) in pairs(x)
            if isnan(xv)
                println("NaN found in parameter: $path[$idx]")
            end
        end
        return nothing
    end

    # Skip generic arrays of non-numeric types to keep debug scans bounded.
    if x isa AbstractArray
        return nothing
    end

    T = typeof(x)
    if isstructtype(T)
        for field in fieldnames(T)
            val = getfield(x, field)
            _debug_print_nan_parameter_paths!(val, string(path, ".", field))
        end
    end
    return nothing
end

@inline function _results_bundle_prefix(args)::String
    return joinpath(args.simulation_settings.results_directory, "simulation_results")
end

function _atomic_write_file(path::String, writer::Function)
    dir = dirname(path)
    mkpath(dir)
    tmp = joinpath(
        dir,
        string(".", basename(path), ".tmp.", getpid(), ".", Threads.threadid(), ".", rand(UInt))
    )
    try
        writer(tmp)
        mv(tmp, path; force=true)
    finally
        if isfile(tmp)
            rm(tmp; force=true)
        end
    end
    return path
end

function _sha256_hex(path::String)::String
    open(path, "r") do io
        return bytes2hex(SHA.sha256(read(io)))
    end
end

@inline function _checkpoint_directory(args)::String
    if isempty(args.simulation_settings.checkpoint_directory)
        return joinpath(args.simulation_settings.results_directory, "checkpoints")
    end
    return args.simulation_settings.checkpoint_directory
end

@inline function _checkpoint_paths(args)
    ckpt_dir = _checkpoint_directory(args)
    return (
        data=joinpath(ckpt_dir, "simulation_checkpoint.bin"),
        manifest=joinpath(ckpt_dir, "simulation_checkpoint.manifest.toml")
    )
end

function _write_checkpoint!(args, t::Float64, u_state)
    paths = _checkpoint_paths(args)
    payload = (
        schema_version=CHECKPOINT_SCHEMA_VERSION,
        created_utc=string(now(UTC)),
        t=t,
        u=deepcopy(u_state)
    )
    _atomic_write_file(paths.data, tmp -> open(tmp, "w") do io
        serialize(io, payload)
    end)

    manifest = Dict{String, Any}(
        "schema_version" => CHECKPOINT_SCHEMA_VERSION,
        "created_utc" => string(now(UTC)),
        "time_s" => t,
        "data_path" => paths.data,
        "data_size_bytes" => filesize(paths.data),
        "data_sha256" => _sha256_hex(paths.data)
    )
    _atomic_write_file(paths.manifest, tmp -> open(tmp, "w") do io
        TOML.print(io, manifest)
    end)
    return nothing
end

function _load_checkpoint(args)
    paths = _checkpoint_paths(args)
    if !isfile(paths.data)
        return nothing
    end
    payload = open(paths.data, "r") do io
        deserialize(io)
    end
    if !haskey(payload, :t) || !haskey(payload, :u)
        throw(ArgumentError("Checkpoint payload missing required keys (:t, :u)."))
    end
    return (t=Float64(payload[:t]), u=payload[:u], data_path=paths.data, manifest_path=paths.manifest)
end

function _clear_checkpoint!(args)
    paths = _checkpoint_paths(args)
    isfile(paths.data) && rm(paths.data; force=true)
    isfile(paths.manifest) && rm(paths.manifest; force=true)
    return nothing
end

function _append_saved_segment!(times_acc::Vector{Float64}, data_acc::Vector, saved_values)
    seg_len = length(saved_values.t)
    seg_len == 0 && return nothing
    start_idx = 1
    if !isempty(times_acc) && isapprox(times_acc[end], saved_values.t[1]; atol=0.0, rtol=0.0)
        start_idx = 2
    end
    if start_idx <= seg_len
        append!(times_acc, @view saved_values.t[start_idx:seg_len])
        append!(data_acc, @view saved_values.saveval[start_idx:seg_len])
    end
    return nothing
end

@inline function _is_flat_scalar(value)::Bool
    return value === missing || value === nothing || value isa Number || value isa AbstractString || value isa Symbol || value isa Bool
end

function _find_sample_value(series)
    for value in series
        if value !== nothing
            return value
        end
    end
    return nothing
end

function _append_series_columns!(results_df::DataFrame, prefix::String, series)
    sample = _find_sample_value(series)
    if sample === nothing || _is_flat_scalar(sample)
        results_df[!, prefix] = collect(series)
        return nothing
    end

    if sample isa NamedTuple
        for key in keys(sample)
            child_series = [value === nothing ? nothing : getproperty(value, key) for value in series]
            _append_series_columns!(results_df, string(prefix, "_", key), child_series)
        end
        return nothing
    end

    if sample isa AbstractDict
        for key in sort!(collect(keys(sample)); by=string)
            child_series = [value === nothing ? nothing : value[key] for value in series]
            _append_series_columns!(results_df, string(prefix, "_", key), child_series)
        end
        return nothing
    end

    if sample isa Tuple || sample isa AbstractArray
        for idx in eachindex(sample)
            child_series = [value === nothing ? nothing : value[idx] for value in series]
            _append_series_columns!(results_df, string(prefix, "_", idx), child_series)
        end
        return nothing
    end

    results_df[!, prefix] = collect(series)
    return nothing
end

function _append_save_field_columns!(results_df::DataFrame, field, saved_data::Vector, num_sats::Int)
    field_series = [snapshot[field.name] for snapshot in saved_data]
    if field.per_satellite
        for sat_idx in 1:num_sats
            sat_series = [value[sat_idx] for value in field_series]
            _append_series_columns!(results_df, "sc$(sat_idx)_$(field.column_prefix)", sat_series)
        end
        return nothing
    end
    _append_series_columns!(results_df, field.column_prefix, field_series)
    return nothing
end

function _build_results_dataframe(times::Vector{Float64}, saved_data::Vector, save_fields, args)::DataFrame
    results_df = DataFrame(time=times)
    num_sats = length(args.dynamics_model.spacecraft)
    for field in save_fields
        _append_save_field_columns!(results_df, field, saved_data, num_sats)
    end
    return results_df
end

function _write_results_bundle!(results_df::DataFrame, times::Vector{Float64}, args)
    prefix = _results_bundle_prefix(args)
    feather_path = prefix * ".feather"
    manifest_path = prefix * ".manifest.toml"

    _atomic_write_file(feather_path, tmp -> Arrow.write(tmp, results_df))

    files = Dict{String, Any}()
    files["feather"] = Dict(
        "path" => feather_path,
        "size_bytes" => filesize(feather_path),
        "sha256" => _sha256_hex(feather_path)
    )

    if args.simulation_settings.save_csv
        csv_path = prefix * ".csv"
        _atomic_write_file(csv_path, tmp -> CSV.write(tmp, results_df))
        files["csv"] = Dict(
            "path" => csv_path,
            "size_bytes" => filesize(csv_path),
            "sha256" => _sha256_hex(csv_path)
        )
    end

    manifest = Dict{String, Any}(
        "schema_version" => RESULTS_BUNDLE_SCHEMA_VERSION,
        "created_utc" => string(now(UTC)),
        "mission_time_s" => args.mission_configuration.mission_time,
        "steps" => length(times),
        "spacecraft_count" => length(args.dynamics_model.spacecraft),
        "orientation_sim" => args.mission_configuration.orientation_sim,
        "files" => files
    )

    _atomic_write_file(manifest_path, tmp -> begin
        open(tmp, "w") do io
            TOML.print(io, manifest)
        end
    end)

    return nothing
end

function run_simulation(
    args;
    isolate_state::Bool=true,
    return_solution::Bool=false,
    return_solver_metadata::Bool=false,
    save_fields=nothing
)
    return SimulationModel.ParallelPolicy.with_policy_context() do
    # Isolate mutable campaign/model state by default so repeated/concurrent runs
    # do not alias shared in-memory objects.
    args = isolate_state ? deepcopy(args) : args

    # Typed pipeline is SI-native (meters, seconds, kilograms). The
    # `simulation_settings.normalize` field is legacy-only and rejected by default.
    _enforce_typed_normalize_policy!(args)
    _validate_orientation_inertia!(args)
    _validate_thermal_model_support!(args)
    try
        SimulationModel.ParallelPolicy.reset_policy_telemetry!()
    catch
    end

    # Set up the model and initial conditions
    initial_conditions = build_initial_conditions(args)
    if args.simulation_settings.verbose
        println("Initial conditions:")
        println(initial_conditions)
    end

    # Define the ODE parameters and callbacks
    p = SimulationModel.ODEParams{length(args.dynamics_model.spacecraft)}(args=args) # Define the parameters for the ODE problem, including the shared buffers for the callbacks
    _initialize_heat_rate_buffers!(p)
    _initialize_density_model_instances!(p)
    _initialize_density_cache_buffers!(p)
    p.shared_buffers.debug_control[] = get(ENV, "SPACEAGORA_DEBUG_CONTROL", "0") == "1"
    p.shared_buffers.debug_initial_derivative[] = get(ENV, "SPACEAGORA_DEBUG_INITIAL_DERIVATIVE", "0") == "1"
    save_fields_resolved = isnothing(save_fields) ? SimulationModel.default_save_fields(args) : collect(save_fields)
    save_field_names = Symbol[field.name for field in save_fields_resolved]
    length(unique(save_field_names)) == length(save_field_names) || throw(ArgumentError("save_fields names must be unique. Got $(save_field_names)."))
    saved_values = SavedValues(Float64, SimulationModel.SaveData)
    callbacks = SimulationModel.get_callbacks(
        length(args.dynamics_model.spacecraft),
        args.dynamics_model.dynamic_effectors,
        args;
        saved_values=saved_values,
        save_fields=save_fields_resolved
    ) # Get the callbacks based on the number of satellites and the dynamic effectors being used in the simulation
    initial_time = args.initial_time
    start_epoch = from_utc(DateTime(
            initial_time.year,
            initial_time.month,
            initial_time.day,
            initial_time.hour,
            initial_time.minute,
            initial_time.second
        ))
    et_start = lock(SimulationModel.SPICE_LOCK) do
        utc2et(to_utc(start_epoch))
    end
    p.shared_buffers.et_start[] = et_start
    lock(SimulationModel.SPICE_LOCK) do
        args.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_$(args.environment_model.planet.name)", et_start)) * args.environment_model.planet.J2000_to_pci' # Initialize the planet frame at the start of the simulation (will be updated in the callback)
    end
    mission_end = args.mission_configuration.mission_time
    checkpoint_active = _typed_checkpoint_enabled(args)
    if checkpoint_active && args.simulation_settings.checkpoint_interval_s <= 0.0
        throw(ArgumentError("SimulationSettings.checkpoint_interval_s must be > 0 when checkpointing is enabled."))
    end

    t_start = 0.0
    u_start = initial_conditions
    if args.simulation_settings.resume_from_checkpoint
        ckpt = _load_checkpoint(args)
        if ckpt === nothing
            @warn "resume_from_checkpoint=true but no checkpoint file was found; starting from initial conditions."
        else
            t_start = ckpt.t
            u_start = ckpt.u
            if args.simulation_settings.verbose
                println("Resuming simulation from checkpoint at t=$(round(t_start, digits=6)) s")
            end
        end
    end

    # println("Initial conditions:")
    # println(initial_conditions)
    # println("ODE parameters:")
    # println(p)
    # println("args.mission_configuration.mission_time: $(args.mission_configuration.mission_time)")
    prob_debug = ODEProblem(spacecraft_dynamics!, u_start, (t_start, mission_end), p, callback=callbacks)
    if p.shared_buffers.debug_initial_derivative[]
        # 1. Manually evaluate the derivative at the start
        du_test = copy(prob_debug.u0)
        try
            prob_debug.f(du_test, prob_debug.u0, prob_debug.p, prob_debug.tspan[1])
        catch e
            @error "The derivative function itself crashed!" exception=e
            throw(ErrorException("Initial derivative evaluation failed; aborting solve in debug mode."))
        end

        # 2. Check for NaNs and print exactly where they are
        if any(isnan, du_test)
            println("--- INITIAL NaN DETECTED ---")

            # Check global parameters in p
            _debug_print_nan_parameter_paths!(prob_debug.p)

            # Check the state vector (u)
            # Assuming your u has a .sc field for satellites
            for (i, sat) in enumerate(du_test.sc)
                if any(isnan, sat.pos) || any(isnan, sat.vel)
                    println("NaN found in Satellite $i derivative!")
                    println("  Pos: $(sat.pos)")
                    println("  Vel: $(sat.vel)")
                end
            end

            throw(ErrorException("Initial derivative contains NaN; aborting solve in debug mode."))
        end
    end

    reltol_tol, abstol_tol = _build_solver_tolerances(u_start, args)
    last_sol = nothing
    solver_trace = NamedTuple[]
    checkpoint_saved_times = Float64[]
    checkpoint_saved_data = SimulationModel.SaveData[]

    if t_start < mission_end && checkpoint_active
        interval = args.simulation_settings.checkpoint_interval_s
        t_cursor = t_start
        u_cursor = deepcopy(u_start)

        while t_cursor < mission_end
            t_next = min(t_cursor + interval, mission_end)
            empty!(saved_values.t)
            empty!(saved_values.saveval)
            prob = ODEProblem(spacecraft_dynamics!, u_cursor, (t_cursor, t_next), p, callback=callbacks)
            seg_sol, solve_meta = _solve_with_solver_policy(prob, args, reltol_tol, abstol_tol)
            push!(solver_trace, solve_meta)
            if !SciMLBase.successful_retcode(seg_sol.retcode)
                throw(ErrorException("Checkpointed solve failed with retcode=$(seg_sol.retcode)."))
            end
            _append_saved_segment!(checkpoint_saved_times, checkpoint_saved_data, saved_values)
            last_sol = seg_sol
            t_cursor = Float64(seg_sol.t[end])
            u_cursor = deepcopy(seg_sol.u[end])
            _write_checkpoint!(args, t_cursor, u_cursor)
        end

    elseif t_start < mission_end
        prob = ODEProblem(spacecraft_dynamics!, u_start, (t_start, mission_end), p, callback=callbacks)
        sol, solve_meta = _solve_with_solver_policy(prob, args, reltol_tol, abstol_tol)
        push!(solver_trace, solve_meta)
        if !SciMLBase.successful_retcode(sol.retcode)
            throw(ErrorException("Solve failed with retcode=$(sol.retcode)."))
        end
        last_sol = sol
    end

    # Process and save results
    if args.simulation_settings.results
        results_times = checkpoint_active ? checkpoint_saved_times : saved_values.t
        results_data = checkpoint_active ? checkpoint_saved_data : saved_values.saveval
        results_df = _build_results_dataframe(results_times, results_data, save_fields_resolved, args)
        # Keep backwards-compatible CSV contract used by existing scripts/tests.
        CSV.write("simulation_results.csv", results_df)
        if _typed_save_bundle_enabled()
            _write_results_bundle!(results_df, results_times, args)
        end
    end

    if return_solution
        if checkpoint_active && args.simulation_settings.checkpoint_interval_s < mission_end
            @warn "return_solution=true with checkpointed integration returns the final segment ODESolution, not a stitched full-history ODESolution."
        end
        if return_solver_metadata
            return (
                solution=last_sol,
                solver_mode=string(_solver_policy_mode()),
                solver_trace=solver_trace
            )
        end
        return last_sol
    end
    return nothing
    end
end

function spacecraft_dynamics!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    # Unpack the state vector
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    debug_control = p.shared_buffers.debug_control[]
    p.shared_buffers.current_time[] = t
    effector_decision = _dynamic_effector_thread_decision(p.args, dynamic_effectors, length(spacecraft))
    minbatch = Int(ceil(length(spacecraft) / Polyester.num_cores())) # Determine the batch size for LoopVectorization based on the number of spacecraft and available CPU cores
    # Loop over each spacecraft and compute its dynamics
    @batch minbatch=minbatch for i in eachindex(sc_state)
        if !p.is_active[i]
            sc_du[i] .= 0.0 # Set the derivatives to zero for inactive spacecraft
            continue
        end
        @views begin
            sc_view = sc_state[i]
            du_view = sc_du[i]
            # println("Computing dynamics for spacecraft $i at time $t seconds...")
            # Compute forces and torques using the dynamic effectors
            forces = MVector{3, Float64}(0.0, 0.0, 0.0)
            torques = MVector{3, Float64}(0.0, 0.0, 0.0)
            mass_rate = 0.0
            effector_started_ns = time_ns()
            if effector_decision.use_threads
                fx = Threads.Atomic{Float64}(0.0)
                fy = Threads.Atomic{Float64}(0.0)
                fz = Threads.Atomic{Float64}(0.0)
                tx = Threads.Atomic{Float64}(0.0)
                ty = Threads.Atomic{Float64}(0.0)
                tz = Threads.Atomic{Float64}(0.0)
                SimulationModel.ParallelPolicy.threaded_foreach(length(dynamic_effectors), effector_decision.allotment) do eff_idx
                    effector = dynamic_effectors[eff_idx]
                    force, torque = SimulationModel.calcForceTorque(effector, sc_view, p, i)
                    Threads.atomic_add!(fx, force[1])
                    Threads.atomic_add!(fy, force[2])
                    Threads.atomic_add!(fz, force[3])
                    Threads.atomic_add!(tx, torque[1])
                    Threads.atomic_add!(ty, torque[2])
                    Threads.atomic_add!(tz, torque[3])
                end
                forces .= SVector{3, Float64}(fx[], fy[], fz[])
                torques .= SVector{3, Float64}(tx[], ty[], tz[])
            else
                @inbounds for effector in dynamic_effectors
                    force, torque = SimulationModel.calcForceTorque(effector, sc_view, p, i)
                    forces .+= force
                    torques .+= torque
                end
            end
            if effector_decision.policy_applied
                SimulationModel.ParallelPolicy.record_policy_observation!(
                    :dynamic_effectors;
                    mode=effector_decision.mode,
                    num_items=length(dynamic_effectors),
                    use_threads=effector_decision.use_threads,
                    elapsed_ns=(time_ns() - effector_started_ns)
                )
            end

            # Compute control forces and torques using the control effectors (if any)
                @inbounds for control_effector in p.args.control_model.control_effectors
                    control_force, control_torque = SimulationModel.calcControlForceTorque(control_effector, sc_view, p, i, t)
                    control_mass_rate = SimulationModel.calcControlMassFlowRate(control_effector, sc_view, p, i, t)
                    if debug_control && (norm(control_force) > 0.0 || norm(control_torque) > 0.0)
                        println("Applying control effect for spacecraft $i at time $t seconds:")
                        println("  Control force: $control_force")
                        # println("  Control torque: $control_torque")
                    end
                # println("Control force for spacecraft $i at time $t seconds: $control_force")
                forces .+= control_force
                torques .+= control_torque
                mass_rate += isfinite(control_mass_rate) ? control_mass_rate : 0.0
            end

            # Update the derivatives of position and velocity
            du_view.pos .= sc_view.vel
            du_view.vel .= forces / sc_view.mass
            du_view.mass = mass_rate

            if p.args.mission_configuration.orientation_sim
                # Update the derivatives of orientation (quaternion) and angular velocity
                ω_body = SVector{3, Float64}(sc_view.ω)
                inertia_tensor = spacecraft[i].inertia_tensor
                τ_body = SVector{3, Float64}(torques)

                du_view.q .= 0.5 * quat_mult(SVector{4, Float64}(ω_body..., 0.0), sc_view.q)
                du_view.ω .= inertia_tensor \ (τ_body - cross(ω_body, inertia_tensor * ω_body))
            end

            # Thermal callback precomputes per-link heat-rate derivatives once per step.
            heat_rates = p.shared_buffers.heat_rates[i]
            if length(heat_rates) == length(du_view.heat_loads)
                du_view.heat_loads .= heat_rates
            else
                du_view.heat_loads .= 0.0
                n_copy = min(length(heat_rates), length(du_view.heat_loads))
                @inbounds for j in 1:n_copy
                    du_view.heat_loads[j] = heat_rates[j]
                end
            end
        end
    end
end # function spacecraft_dynamics!

function build_initial_conditions(args)::ComponentVector
    # 1. Build the structure (Axis) based on each spacecraft's unique body count
    # This identifies exactly how many heat_load slots each SC needs
    sc_shapes = map(args.dynamics_model.spacecraft) do sc
        # Get the number of bodies for this specific spacecraft
        n_bodies = length(sc.links)
        mass = sc.dry_mass + sc.prop_mass
        # Create the initial state for this spacecraft with the correct size for heat_loads
        if args.mission_configuration.orientation_sim
            return (
                pos = zeros(3), 
                vel = zeros(3), 
                mass = mass, 
                heat_loads = zeros(n_bodies), # Variable size!
                q = Float64[0.0, 0.0, 0.0, 1.0], 
                ω = zeros(3)
            )
        else
            return (
                pos = zeros(3), 
                vel = zeros(3), 
                mass = mass, 
                heat_loads = zeros(n_bodies)
            )
        end
    end

    # 2. Pack everything into one ComponentVector
    # Julia allocates ONE flat array and calculates all offsets automatically
    state = ComponentVector(sc = sc_shapes) # Add more components here as needed in the future (e.g., separate orientation state if not using quaternions, etc.)

    # 3. Fill the values (The logic remains the same)
    for i in eachindex(args.dynamics_model.spacecraft)
        spacecraft = args.dynamics_model.spacecraft[i]
        sc_view = state.sc[i]
        r0, v0 = orbitalelemtorv(spacecraft.initial_condition, args.environment_model.planet)
        sc_view.pos .= r0
        sc_view.vel .= v0
        # sc_view.mass .= spacecraft.dry_mass + spacecraft.prop_mass
        # Note: heat_loads is already the correct size for this specific i!
        sc_view.heat_loads .= 0.0  
        
        if args.mission_configuration.orientation_sim
            sc_view.q .= spacecraft.initial_condition.q
            sc_view.ω .= spacecraft.initial_condition.ang_vel
        end
    end

    return state
end
