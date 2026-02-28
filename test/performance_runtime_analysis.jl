const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_OUTPUT_DIR = joinpath(REPO_ROOT, "output", "performance")
const SPICE_PATH = joinpath(REPO_ROOT, "GRAM Suite 2.0", "SPICE")
const EARTH_HARMONICS_FILE = joinpath(REPO_ROOT, "Gravity_harmonics_data", "EarthGGM05C.csv")
const _PERF_POLICY_ENV_NAMES = (
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE",
    "SPACEAGORA_DENSITY_CALLBACK_PARALLEL",
    "SPACEAGORA_CONTROL_CALLBACK_PARALLEL",
    "SPACEAGORA_MULTIBODY_PARALLEL",
)
const _PERF_POLICY_ENV_BASELINE = Dict{String, Union{Nothing, String}}(
    name => (haskey(ENV, name) ? String(ENV[name]) : nothing)
    for name in _PERF_POLICY_ENV_NAMES
)

using CSV
using DataFrames
using Dates
using Distributed
using LinearAlgebra
using Random
using SPICE
using StaticArrays
using Statistics

if myid() == 1
    ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
    using Plots
end

include(joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
using .SimulationModel

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
include(joinpath(REPO_ROOT, "src", "simulation", "run_simulation.jl"))

Base.@kwdef struct ProfileSpec
    name::String
    repeats::Int
    warmup::Int
    max_attempts::Int
    mission_short_s::Float64
    mission_long_s::Float64
    montecarlo_samples::Int
    montecarlo_mission_s::Float64
end

Base.@kwdef struct BenchmarkCase
    name::String
    category::String
    description::String
    args_template::SimulationConfiguration
    run_in_quick::Bool = true
end

@inline _safe_div(num::Float64, den::Float64) = den > 0.0 ? num / den : NaN

@inline function _alloc_calls(gcstats::Base.GC_Diff)::Int
    return gcstats.malloc + gcstats.realloc + gcstats.poolalloc + gcstats.bigalloc
end

@inline function _destat_int(sol, field::Symbol)
    if hasproperty(sol.destats, field)
        value = getproperty(sol.destats, field)
        return value isa Integer ? Int(value) : try
            Int(value)
        catch
            missing
        end
    end
    return missing
end

@inline function _solve_success(sol)::Bool
    retcode = string(sol.retcode)
    # For paper/runtime benchmarking, only full mission completion is valid.
    # "Terminated" usually means an early-stop callback path (e.g., impact),
    # which makes timing samples non-comparable.
    return retcode == "Success"
end

@inline function _effector_signature(effectors::Tuple)
    isempty(effectors) && return "none"
    return join([string(nameof(typeof(e))) for e in effectors], "+")
end

@inline function _safe_unique_join(values; delimiter::String=",")
    vec = collect(skipmissing(values))
    isempty(vec) && return missing
    return join(sort(unique(string.(vec))), delimiter)
end

@inline _perf_error_text(err) = sprint(showerror, err)

@inline function _perf_strict_errors()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_PERF_STRICT_ERRORS", "0")))
    return raw in ("1", "true", "yes", "on")
end

@inline function _perf_stack_head(bt)::String
    bt === nothing && return "stack=unavailable"
    st = stacktrace(bt)
    if isempty(st)
        return "stack=empty"
    end
    for frame in st
        file = String(frame.file)
        if !(occursin("/julia/", file) || occursin("boot.jl", file) || occursin("task.jl", file))
            return string(file, ":", frame.line, " in ", frame.func)
        end
    end
    frame = st[1]
    return string(String(frame.file), ":", frame.line, " in ", frame.func)
end

@inline function _solve_retcode_from_error(err)::Union{Missing, String}
    msg = _perf_error_text(err)
    m = match(r"retcode=([A-Za-z0-9_]+)", msg)
    return m === nothing ? missing : String(m.captures[1])
end

@inline function _profile_from_name(name::String)::ProfileSpec
    if name == "full"
        return ProfileSpec(
            name="full",
            repeats=5,
            warmup=2,
            max_attempts=4,
            mission_short_s=3600.0,
            mission_long_s=14400.0,
            montecarlo_samples=20,
            montecarlo_mission_s=3600.0
        )
    elseif name == "quick"
        return ProfileSpec(
            name="quick",
            repeats=3,
            warmup=1,
            max_attempts=4,
            mission_short_s=1800.0,
            mission_long_s=7200.0,
            montecarlo_samples=8,
            montecarlo_mission_s=1800.0
        )
    else
        throw(ArgumentError("Unsupported profile '$name'. Valid values: quick, full."))
    end
end

function parse_cli()
    profile_name = lowercase(get(ENV, "SPACEAGORA_PERF_PROFILE", "quick"))
    outdir = get(ENV, "SPACEAGORA_PERF_OUTDIR", DEFAULT_OUTPUT_DIR)

    for arg in ARGS
        if arg in ("quick", "full")
            profile_name = arg
        elseif startswith(arg, "--profile=")
            profile_name = lowercase(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]
        else
            throw(ArgumentError("Unknown argument '$arg'. Supported: [quick|full], --profile=..., --outdir=..."))
        end
    end

    return _profile_from_name(profile_name), abspath(outdir)
end

@inline function perf_parallel_enabled()::Bool
    mode = lowercase(strip(get(ENV, "SPACEAGORA_PERF_PARALLEL", "auto")))
    if mode in ("0", "false", "off", "no")
        return false
    end
    return Threads.nthreads() > 1
end


@inline function perf_parallel_backend()::Symbol
    mode = lowercase(strip(get(ENV, "SPACEAGORA_PERF_PARALLEL_BACKEND", "threads")))
    if mode in ("none", "serial", "off", "0", "false", "no")
        return :none
    elseif mode in ("threads", "thread")
        return perf_parallel_enabled() ? :threads : :none
    elseif mode in ("process", "processes", "distributed", "pmap")
        return :process
    elseif mode == "auto"
        return :auto
    else
        throw(ArgumentError("Unsupported SPACEAGORA_PERF_PARALLEL_BACKEND='$mode'. Use one of: threads, process, none, auto."))
    end
end

@inline _threads_or_none_backend()::Symbol = perf_parallel_enabled() ? :threads : :none

Base.@kwdef struct ParallelPriorityPlan
    outer_route::Symbol = :none
    density_mode::String = "off"
    control_mode::String = "off"
    multibody_mode::String = "off"
end

@inline function _parse_positive_int_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, string(default)))
    value = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be an integer, got '$raw'"))
    end
    return max(1, value)
end

@inline function _parse_nonnegative_float_env(name::String, default::Float64)::Float64
    raw = strip(get(ENV, name, string(default)))
    value = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a number, got '$raw'"))
    end
    return max(0.0, value)
end

@inline function _priority_inner_sat_threshold()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_PRIORITY_INNER_SAT_THRESHOLD", 8)
end

@inline function _priority_inner_link_threshold()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_PRIORITY_INNER_LINK_THRESHOLD", 12)
end

@inline function _priority_outer_light_sat_threshold()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_PRIORITY_OUTER_LIGHT_SAT_THRESHOLD", 2)
end

@inline function _priority_outer_light_link_threshold()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_PRIORITY_OUTER_LIGHT_LINK_THRESHOLD", 4)
end

@inline function _priority_outer_light_mission_threshold_s()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_PRIORITY_OUTER_LIGHT_MISSION_THRESHOLD_S", 14_400.0)
end

@inline function _machine_parallel_class()::Symbol
    cpu_threads = Sys.CPU_THREADS
    if cpu_threads >= 24
        return :large
    elseif cpu_threads >= 12
        return :medium
    end
    return :small
end

@inline function _total_links(case::BenchmarkCase)::Int
    total = 0
    @inbounds for sc in case.args_template.dynamics_model.spacecraft
        total += length(sc.links)
    end
    return total
end

@inline function _has_atmo_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa AerodynamicCoefficientConstant || effector isa AerodynamicCoefficientfM || effector isa AerodynamicCoefficientNoBallisticFlight
            return true
        end
    end
    return false
end

@inline function _has_nbody_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa NBodyGravityModel
            return true
        end
    end
    return false
end

@inline function _max_harmonics_degree(effectors::Tuple)::Int
    degree = 0
    @inbounds for effector in effectors
        if effector isa GravitationalHarmonicsModel
            degree = max(degree, effector.L)
        end
    end
    return degree
end

@inline function _priority_outer_route(case::BenchmarkCase)::Symbol
    if !perf_parallel_enabled()
        return :none
    end

    n_sats = length(case.args_template.dynamics_model.spacecraft)
    n_links = _total_links(case)
    mission_time_s = case.args_template.mission_configuration.mission_time
    dynamic_effectors = case.args_template.dynamics_model.dynamic_effectors
    has_nbody = _has_nbody_dynamic_effector(dynamic_effectors)
    harmonics_degree = _max_harmonics_degree(dynamic_effectors)

    machine = _machine_parallel_class()
    if case.category == "montecarlo"
        return machine in (:large, :medium) ? :process : _threads_or_none_backend()
    end

    # Lightweight scenarios are usually overhead-bound, so keep them serial.
    if n_sats <= _priority_outer_light_sat_threshold() &&
       n_links <= _priority_outer_light_link_threshold() &&
       mission_time_s <= _priority_outer_light_mission_threshold_s() &&
       !has_nbody &&
       harmonics_degree == 0
        return :none
    end

    sat_threshold = _priority_inner_sat_threshold()
    link_threshold = _priority_inner_link_threshold()

    # Keep one dominant layer: if the case is large enough for intra-sim parallelism,
    # avoid also parallelizing the outer loop.
    if n_sats >= sat_threshold || n_links >= link_threshold
        return :none
    end

    if case.category == "satellite_scaling" && n_sats >= 4
        return _threads_or_none_backend()
    end
    if has_nbody || harmonics_degree >= 20
        return machine in (:large, :medium) ? :process : _threads_or_none_backend()
    end
    if machine in (:large, :medium)
        return :process
    end
    return _threads_or_none_backend()
end

@inline function parallel_priority_plan(case::BenchmarkCase, outer_route::Symbol)::ParallelPriorityPlan
    resolved_outer = outer_route == :auto ? _priority_outer_route(case) : outer_route
    n_sats = length(case.args_template.dynamics_model.spacecraft)
    n_links = _total_links(case)
    dynamic_effectors = case.args_template.dynamics_model.dynamic_effectors
    has_atmo = _has_atmo_dynamic_effector(dynamic_effectors)
    has_nbody = _has_nbody_dynamic_effector(dynamic_effectors)
    has_density = has_atmo || !(case.args_template.environment_model.density_model isa NoAtmosphereModel)
    has_control = !isempty(case.args_template.control_model.control_effectors)
    sat_threshold = _priority_inner_sat_threshold()
    link_threshold = _priority_inner_link_threshold()

    if resolved_outer != :none
        return ParallelPriorityPlan(
            outer_route=resolved_outer,
            density_mode="off",
            control_mode="off",
            multibody_mode="off"
        )
    end

    density_mode = has_density && n_sats >= sat_threshold ? "auto" : "off"
    control_mode = has_control && n_sats >= sat_threshold ? "auto" : "off"
    multibody_mode = (has_atmo || has_nbody) && n_links >= link_threshold ? "auto" : "off"

    return ParallelPriorityPlan(
        outer_route=resolved_outer,
        density_mode=density_mode,
        control_mode=control_mode,
        multibody_mode=multibody_mode
    )
end

@inline function _existing_or_policy_env(name::String, value::String)::String
    # Only respect user-provided baseline overrides captured at startup.
    # This avoids inheriting transient values left by other benchmark phases.
    baseline = get(_PERF_POLICY_ENV_BASELINE, name, nothing)
    return baseline === nothing ? value : baseline
end

@inline function parallel_priority_env_pairs(plan::ParallelPriorityPlan)::Vector{Pair{String, String}}
    outer_flag = plan.outer_route == :none ? "0" : "1"
    return Pair{String, String}[
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => _existing_or_policy_env("SPACEAGORA_OUTER_PARALLEL_ACTIVE", outer_flag),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => _existing_or_policy_env("SPACEAGORA_DENSITY_CALLBACK_PARALLEL", plan.density_mode),
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => _existing_or_policy_env("SPACEAGORA_CONTROL_CALLBACK_PARALLEL", plan.control_mode),
        "SPACEAGORA_MULTIBODY_PARALLEL" => _existing_or_policy_env("SPACEAGORA_MULTIBODY_PARALLEL", plan.multibody_mode),
    ]
end

@inline function _env_pairs_key(env_pairs::Vector{Pair{String, String}})::String
    isempty(env_pairs) && return ""
    return join([string(first(p), "=", last(p)) for p in env_pairs], ";")
end

function _thread_plan_groups(
    cases::Vector{BenchmarkCase},
    indices::Vector{Int},
    outer_route::Symbol
)::Vector{Tuple{Vector{Pair{String, String}}, Vector{Tuple{Int, ParallelPriorityPlan}}}}
    grouped_pairs = Dict{String, Vector{Pair{String, String}}}()
    grouped_payload = Dict{String, Vector{Tuple{Int, ParallelPriorityPlan}}}()
    ordered_keys = String[]
    for idx in indices
        plan = parallel_priority_plan(cases[idx], outer_route)
        env_pairs = parallel_priority_env_pairs(plan)
        key = _env_pairs_key(env_pairs)
        if !haskey(grouped_payload, key)
            grouped_pairs[key] = env_pairs
            grouped_payload[key] = Tuple{Int, ParallelPriorityPlan}[]
            push!(ordered_keys, key)
        end
        push!(grouped_payload[key], (idx, plan))
    end
    return [(grouped_pairs[key], grouped_payload[key]) for key in ordered_keys]
end

@inline function auto_backend_for_case(case::BenchmarkCase)::Symbol
    return _priority_outer_route(case)
end

@inline function perf_process_workers_target()::Int
    raw = strip(get(ENV, "SPACEAGORA_PERF_PROCS", "0"))
    if isempty(raw) || raw == "0"
        return max(1, Sys.CPU_THREADS - 1)
    end
    n = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_PERF_PROCS must be an integer, got '$raw'"))
    end
    return max(1, n)
end

const _perf_workers_initialized = Ref(false)
const _perf_worker_planet_cache = Ref{Any}(nothing)

function ensure_perf_workers!()
    _perf_workers_initialized[] && return nothing
    target_workers = perf_process_workers_target()
    missing_workers = target_workers - nworkers()
    if missing_workers > 0
        addprocs(
            missing_workers;
            exeflags=`--startup-file=no --project=$(joinpath(REPO_ROOT, ".AGORA"))`
        )
    end
    script_path = abspath(@__FILE__)
    @everywhere workers() include($script_path)
    for w in workers()
        remotecall_wait(perf_worker_planet, w)
    end
    _perf_workers_initialized[] = true
    return nothing
end

function perf_worker_planet()::Earth
    cached = _perf_worker_planet_cache[]
    if cached === nothing
        cached = Earth("", SPICE_PATH)
        _perf_worker_planet_cache[] = cached
    end
    return cached::Earth
end

@inline function orbital_period_seconds(spacecraft::SpacecraftModel, planet::AbstractPlanet)::Float64
    a = spacecraft.initial_condition.a
    if !isfinite(a) || a <= 0.0
        throw(ArgumentError("Invalid semimajor axis for orbital period calculation: $a"))
    end
    return 2π * sqrt(a^3 / planet.μ)
end

function make_spacecraft(
    planet::AbstractPlanet;
    id::Int=1,
    ra_alt_m::Float64=550e3,
    rp_alt_m::Float64=500e3,
    i_deg::Float64=35.0,
    ω_deg::Float64=40.0,
    Ω_deg::Float64=10.0,
    ν_deg::Float64=170.0,
    with_panel::Bool=true,
    orientation_state::Union{Nothing, Tuple{SVector{4, Float64}, SVector{3, Float64}}}=nothing,
    root_mass::Float64=500.0,
    root_area::Float64=12.0,
    panel_mass::Float64=30.0,
    panel_area::Float64=6.0,
    panel_offset_y::Float64=1.2,
    prop_mass::Float64=0.0
)::SpacecraftModel
    root = Link{0}(root=true, m=root_mass, ref_area=root_area)
    links = Link[root]
    if with_panel
        panel = Link{0}(root=false, m=panel_mass, ref_area=panel_area, r=MVector{3, Float64}(0.0, panel_offset_y, 0.0))
        push!(links, panel)
    end

    ra = planet.Rp_e + ra_alt_m
    rp = planet.Rp_e + rp_alt_m

    ic = if isnothing(orientation_state)
        InitialCondition(ra=ra, rp=rp, i=i_deg, ω=ω_deg, Ω=Ω_deg, ν=ν_deg)
    else
        q0, w0 = orientation_state
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        InitialCondition(a, e, i_deg, ω_deg, Ω_deg, ν_deg, q0, w0)
    end

    dry_mass = sum(link.m for link in links)
    return SpacecraftModel(
        Joint[],
        links,
        root,
        true,
        dry_mass,
        prop_mass,
        root.inertia,
        0,
        0,
        ic,
        id
    )
end

function make_constellation(planet::AbstractPlanet, n::Int; with_panel::Bool=false)::Vector{SpacecraftModel}
    sats = SpacecraftModel[]
    for i in 1:n
        ra_alt = 540e3 + 20e3 * (i - 1)
        rp_alt = 470e3 + 15e3 * (i - 1)
        if rp_alt >= ra_alt
            ra_alt = rp_alt + 8e3
        end
        ν = 140.0 + 180.0 * (i - 1) / n
        push!(sats, make_spacecraft(planet; id=i, ra_alt_m=ra_alt, rp_alt_m=rp_alt, ν_deg=ν, with_panel=with_panel))
    end
    return sats
end

function build_config(;
    planet::AbstractPlanet,
    spacecraft::Vector{SpacecraftModel},
    mission_time_s::Float64,
    orientation_sim::Bool,
    dynamic_effectors::Tuple,
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    dt_max_orbit::Float64=1.0,
    reltol_orbit::Float64=1e-9,
    abstol_orbit::Float64=1e-9
)::SimulationConfiguration
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false,
            save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=mission_time_s,
            orientation_sim=orientation_sim,
            num_steps_to_save=400
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel(spacecraft, dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=reltol_orbit,
            abstol_orbit=abstol_orbit,
            dt_max_orbit=dt_max_orbit
        )
    )
end

function make_montecarlo_config(seed::Int, planet::Earth, mission_time_s::Float64)::SimulationConfiguration
    rng = MersenneTwister(seed)
    ra_alt = 530e3 + randn(rng) * 20e3
    rp_alt = 490e3 + randn(rng) * 20e3
    rp_alt = max(rp_alt, 120e3)
    if rp_alt >= ra_alt
        ra_alt = rp_alt + 8e3
    end

    sc = make_spacecraft(
        planet;
        id=1,
        with_panel=false,
        root_mass=160.0,
        root_area=1.5,
        prop_mass=20.0,
        ra_alt_m=ra_alt,
        rp_alt_m=rp_alt,
        i_deg=28.0 + randn(rng) * 0.8,
        ω_deg=15.0 + randn(rng) * 2.0,
        Ω_deg=20.0 + randn(rng) * 2.0,
        ν_deg=160.0 + randn(rng) * 8.0
    )

    thruster = BaseThrusterModel(
        thrust=[0.6 + rand(rng) * 0.5],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[0.0],
        stop_burn_time=[180.0 + rand(rng) * 40.0],
        Isp=[280.0 + rand(rng) * 40.0]
    )

    return build_config(
        planet=planet,
        spacecraft=[sc],
        mission_time_s=mission_time_s,
        orientation_sim=false,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thruster,),
        control_rates=[1.0],
        dt_max_orbit=1.0,
        reltol_orbit=1e-9,
        abstol_orbit=1e-9
    )
end

function build_cases(spec::ProfileSpec, planet::Earth)::Vector{BenchmarkCase}
    harmonics20 = GravitationalHarmonicsModel(20, 20, EARTH_HARMONICS_FILE, planet)
    harmonics50 = GravitationalHarmonicsModel(50, 50, EARTH_HARMONICS_FILE, planet)
    nbody_sun_moon = NBodyGravityModel(body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet)

    q0 = normalize(SVector{4, Float64}(0.15, -0.05, 0.2, 0.96))
    w0 = SVector{3, Float64}(0.01, -0.02, 0.015)

    sc_baseline = [make_spacecraft(planet; id=1, with_panel=false)]
    sc_orientation = [make_spacecraft(planet; id=1, with_panel=true, orientation_state=(q0, w0))]
    sc_thruster = [make_spacecraft(planet; id=1, with_panel=false, prop_mass=30.0)]
    sc_super_constellation = make_constellation(planet, 8; with_panel=false)
    thruster = BaseThrusterModel(
        thrust=[0.8],
        direction=[0.0],
        Δv=[35.0],
        start_burn_time=[-1.0],
        stop_burn_time=[-1.0],
        Isp=[300.0]
    )
    super_constellation_thruster = BaseThrusterModel(
        thrust=fill(0.18, 8),
        direction=fill(0.0, 8),
        Δv=fill(4.0, 8),
        start_burn_time=fill(-1.0, 8),
        stop_burn_time=fill(-1.0, 8),
        Isp=fill(300.0, 8)
    )

    cases = BenchmarkCase[
        BenchmarkCase(
            name="single_baseline_gravity",
            category="core",
            description="1 spacecraft, position-only, inverse-square gravity",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(),)
            )
        ),
        BenchmarkCase(
            name="single_orientation_aero",
            category="orientation",
            description="1 spacecraft, orientation dynamics on, aerodynamic model active",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_orientation,
                mission_time_s=spec.mission_short_s,
                orientation_sim=true,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM())
            )
        ),
        BenchmarkCase(
            name="multi_2_gravity",
            category="satellite_scaling",
            description="2 spacecraft, position-only, inverse-square gravity",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 2; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(),)
            )
        ),
        BenchmarkCase(
            name="multi_4_gravity",
            category="satellite_scaling",
            description="4 spacecraft, position-only, inverse-square gravity",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 4; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(),)
            )
        ),
        BenchmarkCase(
            name="multi_8_gravity",
            category="satellite_scaling",
            description="8 spacecraft, position-only, inverse-square gravity",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 8; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(),)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="single_j2",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square + J2 gravity",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredJ2GravityModel(),)
            )
        ),
        BenchmarkCase(
            name="single_nbody_sun_moon",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square gravity + N-body Sun/Moon perturbations",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), nbody_sun_moon)
            )
        ),
        BenchmarkCase(
            name="single_harmonics_l20",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square gravity + spherical harmonics L=M=20",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics20)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="single_harmonics_l50",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square gravity + spherical harmonics L=M=50",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics50)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="single_thruster_control",
            category="control",
            description="1 spacecraft, inverse-square gravity + BaseThrusterModel control callback",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_thruster,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(),),
                control_effectors=(thruster,),
                control_rates=[1.0]
            )
        ),
        BenchmarkCase(
            name="super_constellation_8sat_l20_control",
            category="control_stress",
            description="8 spacecraft, L20 harmonics + BaseThrusterModel control callback stress case",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_super_constellation,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics20),
                control_effectors=(super_constellation_thruster,),
                control_rates=[1.0]
            )
        ),
        BenchmarkCase(
            name="single_baseline_long_mission",
            category="mission_length",
            description="1 spacecraft baseline gravity with longer mission duration",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_long_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(),)
            )
        )
    ]

    return cases
end

function measure_case(case::BenchmarkCase, profile_name::String, repeat_idx::Int; seed::Union{Missing, Int}=missing, attempt::Int=1)
    timestamp_utc = string(now(UTC))
    args_meta = case.args_template
    n_sats = length(args_meta.dynamics_model.spacecraft)
    mission_time_s = args_meta.mission_configuration.mission_time

    GC.gc()
    copy_timed = @timed deepcopy(case.args_template)
    args_run = copy_timed.value

    GC.gc()
    solve_timed = @timed begin
        try
            result = run_simulation(
                args_run;
                isolate_state=false,
                return_solution=true,
                return_solver_metadata=true
            )
            (ok=true, result=result, err=nothing, bt=nothing)
        catch err
            if err isa InterruptException
                rethrow()
            end
            if _perf_strict_errors()
                rethrow(err)
            end
            (ok=false, result=nothing, err=err, bt=catch_backtrace())
        end
    end

    copy_time_s = copy_timed.time
    solve_time_s = solve_timed.time
    total_time_s = copy_time_s + solve_time_s

    copy_bytes_mb = copy_timed.bytes / 1024^2
    solve_bytes_mb = solve_timed.bytes / 1024^2
    total_bytes_mb = copy_bytes_mb + solve_bytes_mb

    copy_alloc_calls = _alloc_calls(copy_timed.gcstats)
    solve_alloc_calls = _alloc_calls(solve_timed.gcstats)

    solver_mode = missing
    solver_sequence = missing
    solver_fallback_used = missing
    solver_fallback_count = missing
    solver_fallback_trigger = missing
    solve_retcode = missing
    solve_success = false
    saved_points = missing
    accepted_steps = missing
    rejected_steps = missing
    sim_seconds_per_wall_second = missing
    satellite_sim_seconds_per_wall_second = missing

    solve_payload = solve_timed.value
    if solve_payload.ok
        solve_result = solve_payload.result
        sol = solve_result.solution
        solver_mode = solve_result.solver_mode
        solver_trace = solve_result.solver_trace
        solver_sequence = isempty(solver_trace) ? missing : join([meta.solver for meta in solver_trace], "->")
        solver_fallback_count = count(meta -> meta.fallback_used, solver_trace)
        solver_fallback_used = solver_fallback_count > 0
        fallback_triggers = [meta.trigger_retcode for meta in solver_trace if !(meta.trigger_retcode isa Missing)]
        solver_fallback_trigger = isempty(fallback_triggers) ? missing : _safe_unique_join(fallback_triggers; delimiter="|")
        solve_retcode = string(sol.retcode)
        solve_success = _solve_success(sol)
        saved_points = length(sol.t)
        accepted_steps = _destat_int(sol, :naccept)
        rejected_steps = _destat_int(sol, :nreject)
        sim_seconds_per_wall_second = _safe_div(mission_time_s, total_time_s)
        satellite_sim_seconds_per_wall_second = _safe_div(mission_time_s * n_sats, total_time_s)
    else
        solve_err = solve_payload.err
        solve_bt = solve_payload.bt
        solve_retcode = _solve_retcode_from_error(solve_err)
        if ismissing(solve_retcode)
            solve_retcode = "Exception"
            @warn "[perf] case=$(case.name) repeat=$(repeat_idx) attempt=$(attempt) threw $(typeof(solve_err)): $(_perf_error_text(solve_err)) @ $(_perf_stack_head(solve_bt))"
        end
    end

    return (
        profile=profile_name,
        category=case.category,
        scenario=case.name,
        description=case.description,
        seed=seed,
        repeat=repeat_idx,
        attempt=attempt,
        satellites=n_sats,
        orientation=args_meta.mission_configuration.orientation_sim,
        mission_time_s=mission_time_s,
        dynamic_effectors=_effector_signature(args_meta.dynamics_model.dynamic_effectors),
        control_effectors=_effector_signature(args_meta.control_model.control_effectors),
        copy_time_s=copy_time_s,
        solve_time_s=solve_time_s,
        total_time_s=total_time_s,
        copy_compile_time_s=copy_timed.compile_time,
        solve_compile_time_s=solve_timed.compile_time,
        copy_gctime_s=copy_timed.gctime,
        solve_gctime_s=solve_timed.gctime,
        solver_mode=solver_mode,
        solver_sequence=solver_sequence,
        solver_fallback_used=solver_fallback_used,
        solver_fallback_count=solver_fallback_count,
        solver_fallback_trigger=solver_fallback_trigger,
        solve_retcode=solve_retcode,
        solve_success=solve_success,
        copy_bytes_mb=copy_bytes_mb,
        solve_bytes_mb=solve_bytes_mb,
        total_bytes_mb=total_bytes_mb,
        copy_alloc_calls=copy_alloc_calls,
        solve_alloc_calls=solve_alloc_calls,
        saved_points=saved_points,
        accepted_steps=accepted_steps,
        rejected_steps=rejected_steps,
        sim_seconds_per_wall_second=sim_seconds_per_wall_second,
        satellite_sim_seconds_per_wall_second=satellite_sim_seconds_per_wall_second,
        timestamp_utc=timestamp_utc
    )
end

function run_warmup(case::BenchmarkCase, warmup::Int)
    for i in 1:warmup
        args_run = deepcopy(case.args_template)
        try
            run_simulation(args_run; isolate_state=false, return_solution=false)
        catch err
            if err isa InterruptException
                rethrow()
            end
            if _perf_strict_errors()
                rethrow(err)
            end
            err_bt = catch_backtrace()
            solve_retcode = _solve_retcode_from_error(err)
            if ismissing(solve_retcode)
                @warn "[warmup] $(case.name) $(i)/$(warmup) threw $(typeof(err)): $(_perf_error_text(err)) @ $(_perf_stack_head(err_bt)); continuing"
            else
                println("  warmup $(i)/$(warmup): failed retcode=$(solve_retcode), continuing")
            end
        end
    end
    return nothing
end

function _run_case_batch_core!(
    rows::Vector{NamedTuple},
    case::BenchmarkCase,
    spec::ProfileSpec,
    idx::Int,
    total::Int,
    plan::ParallelPriorityPlan
)
    println(
        "[$idx/$total] $(case.name) :: warmup x$(spec.warmup), repeats x$(spec.repeats), " *
        "outer=$(plan.outer_route), density=$(plan.density_mode), control=$(plan.control_mode), multibody=$(plan.multibody_mode)"
    )
    run_warmup(case, spec.warmup)
    for rep in 1:spec.repeats
        last_row = nothing
        for attempt in 1:spec.max_attempts
            row = measure_case(case, spec.name, rep; attempt=attempt)
            last_row = row
            if row.solve_success
                push!(rows, row)
                println("  repeat $(rep)/$(spec.repeats): total=$(round(row.total_time_s; digits=3)) s, solve=$(round(row.solve_time_s; digits=3)) s")
                break
            end
            println("  repeat $(rep)/$(spec.repeats) attempt $(attempt)/$(spec.max_attempts): failed retcode=$(row.solve_retcode), retrying")
        end
        if !(last_row === nothing) && !last_row.solve_success
            push!(rows, last_row)
            println("  repeat $(rep)/$(spec.repeats): failed after $(spec.max_attempts) attempts, retcode=$(last_row.solve_retcode)")
        end
    end
    return nothing
end

function run_case_batch!(
    rows::Vector{NamedTuple},
    case::BenchmarkCase,
    spec::ProfileSpec,
    idx::Int,
    total::Int;
    outer_route::Symbol=:none,
    plan::Union{Nothing, ParallelPriorityPlan}=nothing,
    apply_env::Bool=true
)
    resolved_plan = isnothing(plan) ? parallel_priority_plan(case, outer_route) : plan
    if apply_env
        env_pairs = parallel_priority_env_pairs(resolved_plan)
        return withenv(env_pairs...) do
            _run_case_batch_core!(rows, case, spec, idx, total, resolved_plan)
        end
    end
    return _run_case_batch_core!(rows, case, spec, idx, total, resolved_plan)
end

function measure_montecarlo_seed(
    spec::ProfileSpec,
    planet::Earth,
    mission_time_s::Float64,
    seed::Int;
    outer_route::Symbol=:none,
    plan::Union{Nothing, ParallelPriorityPlan}=nothing,
    apply_env::Bool=true
)
    case = BenchmarkCase(
        name="montecarlo_randomized",
        category="montecarlo",
        description="Randomized initial conditions + thruster, one run per seed",
        args_template=make_montecarlo_config(seed, planet, mission_time_s),
        run_in_quick=true
    )
    resolved_plan = isnothing(plan) ? parallel_priority_plan(case, outer_route) : plan
    run_seed = () -> begin
        last_row = nothing
        for attempt in 1:spec.max_attempts
            row = measure_case(case, spec.name, 1; seed=seed, attempt=attempt)
            last_row = row
            if row.solve_success
                return row, nothing
            end
        end
        return last_row, last_row === nothing ? "failed without attempt data" : "failed after $(spec.max_attempts) attempts, retcode=$(last_row.solve_retcode)"
    end
    if apply_env
        env_pairs = parallel_priority_env_pairs(resolved_plan)
        return withenv(env_pairs...) do
            run_seed()
        end
    end
    return run_seed()
end

function perf_worker_montecarlo_warmup(
    spec::ProfileSpec,
    mission_time_s::Float64,
    seed::Int,
    outer_route::Symbol=:process
)
    planet = perf_worker_planet()
    warmup_case = BenchmarkCase(
        name="montecarlo_randomized",
        category="montecarlo",
        description="Randomized initial conditions + thruster, one run per seed",
        args_template=make_montecarlo_config(seed, planet, mission_time_s),
        run_in_quick=true
    )
    plan = parallel_priority_plan(warmup_case, outer_route)
    env_pairs = parallel_priority_env_pairs(plan)
    withenv(env_pairs...) do
        run_warmup(warmup_case, spec.warmup)
    end
    return nothing
end

function perf_worker_measure_montecarlo_seed(
    spec::ProfileSpec,
    mission_time_s::Float64,
    seed::Int,
    outer_route::Symbol=:process
)
    planet = perf_worker_planet()
    return measure_montecarlo_seed(spec, planet, mission_time_s, seed; outer_route=outer_route)
end

function perf_worker_run_case_batch(
    case::BenchmarkCase,
    spec::ProfileSpec,
    idx::Int,
    total::Int,
    outer_route::Symbol=:process
)
    local_rows = NamedTuple[]
    run_case_batch!(local_rows, case, spec, idx, total; outer_route=outer_route)
    return local_rows
end

function perf_worker_measure_per_orbit_scenario(
    base_case::BenchmarkCase,
    spec::ProfileSpec,
    period_s::Float64,
    orbit_counts::Vector{Int},
    outer_route::Symbol=:process
)
    return measure_per_orbit_scenario(base_case, spec, period_s, orbit_counts; outer_route=outer_route)
end

function run_montecarlo_batch!(rows::Vector{NamedTuple}, spec::ProfileSpec, planet::Earth)
    seeds = collect(1001:(1000 + spec.montecarlo_samples))
    warmup_case = BenchmarkCase(
        name="montecarlo_randomized",
        category="montecarlo",
        description="Randomized initial conditions + thruster, one run per seed",
        args_template=make_montecarlo_config(first(seeds), planet, spec.montecarlo_mission_s),
        run_in_quick=true
    )
    println("[montecarlo] warmup x$(spec.warmup), seeds=$(length(seeds))")

    backend = perf_parallel_backend()
    mc_backend = backend == :auto ? :process : backend
    if mc_backend == :process
        ensure_perf_workers!()
        warmup_seed = first(seeds)
        for w in workers()
            remotecall_wait(perf_worker_montecarlo_warmup, w, spec, spec.montecarlo_mission_s, warmup_seed, :process)
        end
    else
        plan = parallel_priority_plan(warmup_case, mc_backend)
        env_pairs = parallel_priority_env_pairs(plan)
        withenv(env_pairs...) do
            run_warmup(warmup_case, spec.warmup)
        end
    end

    seed_rows = Vector{NamedTuple}(undef, length(seeds))
    seed_msgs = Vector{String}(undef, length(seeds))

    if mc_backend == :process
        seed_results = pmap(seed -> perf_worker_measure_montecarlo_seed(spec, spec.montecarlo_mission_s, seed, :process), seeds)
        for i in eachindex(seeds)
            seed = seeds[i]
            row, err = seed_results[i]
            seed_rows[i] = row
            if row.solve_success
                seed_msgs[i] = "  seed $(i)/$(length(seeds))=$(seed): total=$(round(row.total_time_s; digits=3)) s"
            else
                seed_msgs[i] = "  seed $(i)/$(length(seeds))=$(seed): $(err)"
            end
        end
    elseif mc_backend == :threads
        threaded_plan = parallel_priority_plan(warmup_case, :threads)
        threaded_env = parallel_priority_env_pairs(threaded_plan)
        withenv(threaded_env...) do
            Threads.@threads for i in eachindex(seeds)
                seed = seeds[i]
                row, err = measure_montecarlo_seed(
                    spec,
                    planet,
                    spec.montecarlo_mission_s,
                    seed;
                    outer_route=:threads,
                    plan=threaded_plan,
                    apply_env=false
                )
                seed_rows[i] = row
                if row.solve_success
                    seed_msgs[i] = "  seed $(i)/$(length(seeds))=$(seed): total=$(round(row.total_time_s; digits=3)) s"
                else
                    seed_msgs[i] = "  seed $(i)/$(length(seeds))=$(seed): $(err)"
                end
            end
        end
    else
        for i in eachindex(seeds)
            seed = seeds[i]
            row, err = measure_montecarlo_seed(spec, planet, spec.montecarlo_mission_s, seed; outer_route=:none)
            seed_rows[i] = row
            if row.solve_success
                seed_msgs[i] = "  seed $(i)/$(length(seeds))=$(seed): total=$(round(row.total_time_s; digits=3)) s"
            else
                seed_msgs[i] = "  seed $(i)/$(length(seeds))=$(seed): $(err)"
            end
        end
    end

    for i in eachindex(seeds)
        push!(rows, seed_rows[i])
        println(seed_msgs[i])
    end
    return nothing
end

function run_benchmarks(spec::ProfileSpec, cases::Vector{BenchmarkCase}, planet::Earth)::DataFrame
    selected = spec.name == "full" ? cases : [c for c in cases if c.run_in_quick]
    rows = NamedTuple[]
    total = length(selected)
    backend = perf_parallel_backend()
    if backend == :auto
        case_rows = Vector{Vector{NamedTuple}}(undef, total)
        process_tasks = Tuple{Int, BenchmarkCase}[]
        thread_indices = Int[]
        serial_indices = Int[]

        for (idx, case) in enumerate(selected)
            route = auto_backend_for_case(case)
            if route == :process
                push!(process_tasks, (idx, case))
            elseif route == :threads
                push!(thread_indices, idx)
            else
                push!(serial_indices, idx)
            end
        end

        if !isempty(process_tasks)
            ensure_perf_workers!()
            process_rows = pmap(process_tasks) do task
                idx, case = task
                perf_worker_run_case_batch(case, spec, idx, total, :process)
            end
            for (k, task) in enumerate(process_tasks)
                idx = task[1]
                case_rows[idx] = process_rows[k]
            end
        end

        if !isempty(thread_indices)
            for (env_pairs, payload) in _thread_plan_groups(selected, thread_indices, :threads)
                withenv(env_pairs...) do
                    Threads.@threads for j in eachindex(payload)
                        idx, plan = payload[j]
                        local_rows = NamedTuple[]
                        run_case_batch!(
                            local_rows,
                            selected[idx],
                            spec,
                            idx,
                            total;
                            outer_route=:threads,
                            plan=plan,
                            apply_env=false
                        )
                        case_rows[idx] = local_rows
                    end
                end
            end
        end

        for idx in serial_indices
            local_rows = NamedTuple[]
            run_case_batch!(local_rows, selected[idx], spec, idx, total; outer_route=:none)
            case_rows[idx] = local_rows
        end

        for idx in eachindex(case_rows)
            append!(rows, case_rows[idx])
        end
    elseif backend == :process
        ensure_perf_workers!()
        tasks = collect(enumerate(selected))
        case_rows = pmap(tasks) do task
            idx, case = task
            perf_worker_run_case_batch(case, spec, idx, total, :process)
        end
        for idx in eachindex(case_rows)
            append!(rows, case_rows[idx])
        end
    elseif backend == :threads
        case_rows = Vector{Vector{NamedTuple}}(undef, total)
        thread_indices = collect(eachindex(selected))
        for (env_pairs, payload) in _thread_plan_groups(selected, thread_indices, :threads)
            withenv(env_pairs...) do
                Threads.@threads for j in eachindex(payload)
                    idx, plan = payload[j]
                    local_rows = NamedTuple[]
                    run_case_batch!(
                        local_rows,
                        selected[idx],
                        spec,
                        idx,
                        total;
                        outer_route=:threads,
                        plan=plan,
                        apply_env=false
                    )
                    case_rows[idx] = local_rows
                end
            end
        end
        for idx in eachindex(case_rows)
            append!(rows, case_rows[idx])
        end
    else
        for (idx, case) in enumerate(selected)
            run_case_batch!(rows, case, spec, idx, total; outer_route=:none)
        end
    end

    run_montecarlo_batch!(rows, spec, planet)
    return DataFrame(rows)
end
@inline function selected_cases(spec::ProfileSpec, cases::Vector{BenchmarkCase})::Vector{BenchmarkCase}
    return spec.name == "full" ? cases : [c for c in cases if c.run_in_quick]
end

function measure_per_orbit_scenario(
    base_case::BenchmarkCase,
    spec::ProfileSpec,
    period_s::Float64,
    orbit_counts::Vector{Int};
    outer_route::Symbol=:none,
    apply_env::Bool=true
)
    rows = NamedTuple[]
    logs = String[]
    for orbit_count in orbit_counts
        mission_time = orbit_count * period_s
        args_template = deepcopy(base_case.args_template)
        args_template = SimulationConfiguration(
            file_paths=args_template.file_paths,
            simulation_settings=args_template.simulation_settings,
            mission_configuration=MissionConfiguration(
                mission_type=args_template.mission_configuration.mission_type,
                keplerian=args_template.mission_configuration.keplerian,
                number_of_orbits=args_template.mission_configuration.number_of_orbits,
                mission_time=mission_time,
                orientation_sim=args_template.mission_configuration.orientation_sim,
                num_steps_to_save=args_template.mission_configuration.num_steps_to_save
            ),
            environment_model=args_template.environment_model,
            dynamics_model=args_template.dynamics_model,
            guidance_model=args_template.guidance_model,
            navigation_model=args_template.navigation_model,
            control_model=args_template.control_model,
            initial_time=args_template.initial_time,
            integration_tolerances=args_template.integration_tolerances
        )

        case = BenchmarkCase(
            name=base_case.name,
            category=base_case.category,
            description=base_case.description,
            args_template=args_template,
            run_in_quick=base_case.run_in_quick
        )

        plan = parallel_priority_plan(case, outer_route)
        run_case = () -> begin
            run_warmup(case, spec.warmup)
            for rep in 1:spec.repeats
                last_row = nothing
                for attempt in 1:spec.max_attempts
                    row = measure_case(case, spec.name, rep; attempt=attempt)
                    row_orbit = merge(
                        row,
                        (
                            orbit_count=orbit_count,
                            orbital_period_s=period_s
                        )
                    )
                    last_row = row_orbit
                    if row_orbit.solve_success
                        push!(rows, row_orbit)
                        push!(logs, "    orbit=$(orbit_count) repeat $(rep)/$(spec.repeats): total=$(round(row_orbit.total_time_s; digits=3)) s")
                        break
                    end
                end
                if !(last_row === nothing) && !last_row.solve_success
                    push!(rows, last_row)
                    push!(logs, "    orbit=$(orbit_count) repeat $(rep)/$(spec.repeats): failed after $(spec.max_attempts) attempts, retcode=$(last_row.solve_retcode)")
                end
            end
        end
        if apply_env
            env_pairs = parallel_priority_env_pairs(plan)
            withenv(env_pairs...) do
                run_case()
            end
        else
            run_case()
        end
    end
    return rows, logs
end

function run_montecarlo_per_orbit!(
    rows::Vector{NamedTuple},
    spec::ProfileSpec,
    planet::Earth,
    period_s::Float64,
    orbit_counts::Vector{Int}
)
    seeds = collect(1:spec.montecarlo_samples)
    println("  scenario montecarlo_randomized (per-orbit, seeds=$(length(seeds)))")
    backend = perf_parallel_backend()
    mc_backend = backend == :auto ? :process : backend

    if mc_backend == :process
        ensure_perf_workers!()
    end

    for orbit_count in orbit_counts
        mission_time = orbit_count * period_s
        println("    orbit=$(orbit_count)")
        orbit_rows = Vector{NamedTuple}(undef, length(seeds))
        orbit_msgs = Vector{String}(undef, length(seeds))

        if mc_backend == :process
            seed_results = pmap(seed -> perf_worker_measure_montecarlo_seed(spec, mission_time, seed, :process), seeds)
            for i in eachindex(seeds)
                seed = seeds[i]
                row, err = seed_results[i]
                row_orbit = merge(row, (orbit_count=orbit_count, orbital_period_s=period_s))
                orbit_rows[i] = row_orbit
                if row_orbit.solve_success
                    orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): total=$(round(row_orbit.total_time_s; digits=3)) s"
                else
                    orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): $(err)"
                end
            end
        elseif mc_backend == :threads
            warmup_case = BenchmarkCase(
                name="montecarlo_randomized",
                category="montecarlo",
                description="Randomized initial conditions + thruster, one run per seed",
                args_template=make_montecarlo_config(first(seeds), planet, mission_time),
                run_in_quick=true
            )
            threaded_plan = parallel_priority_plan(warmup_case, :threads)
            threaded_env = parallel_priority_env_pairs(threaded_plan)
            withenv(threaded_env...) do
                Threads.@threads for i in eachindex(seeds)
                    seed = seeds[i]
                    row, err = measure_montecarlo_seed(
                        spec,
                        planet,
                        mission_time,
                        seed;
                        outer_route=:threads,
                        plan=threaded_plan,
                        apply_env=false
                    )
                    row_orbit = merge(row, (orbit_count=orbit_count, orbital_period_s=period_s))
                    orbit_rows[i] = row_orbit
                    if row_orbit.solve_success
                        orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): total=$(round(row_orbit.total_time_s; digits=3)) s"
                    else
                        orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): $(err)"
                    end
                end
            end
        else
            for i in eachindex(seeds)
                seed = seeds[i]
                row, err = measure_montecarlo_seed(spec, planet, mission_time, seed; outer_route=:none)
                row_orbit = merge(row, (orbit_count=orbit_count, orbital_period_s=period_s))
                orbit_rows[i] = row_orbit
                if row_orbit.solve_success
                    orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): total=$(round(row_orbit.total_time_s; digits=3)) s"
                else
                    orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): $(err)"
                end
            end
        end

        for i in eachindex(seeds)
            push!(rows, orbit_rows[i])
            println(orbit_msgs[i])
        end
    end
    return nothing
end
function run_per_orbit_for_scenarios(spec::ProfileSpec, cases::Vector{BenchmarkCase}, planet::Earth)::DataFrame
    baseline_sc = make_spacecraft(planet; id=1, with_panel=false)
    period_s = orbital_period_seconds(baseline_sc, planet)
    orbit_counts = spec.name == "full" ? collect(1:5) : collect(1:3)
    selected = selected_cases(spec, cases)

    println("[per-orbit] scenarios=$(length(selected)), baseline period=$(round(period_s; digits=3)) s, orbit counts=$(first(orbit_counts)):$(last(orbit_counts))")
    rows = NamedTuple[]
    backend = perf_parallel_backend()

    if backend == :auto
        scenario_rows = Vector{Vector{NamedTuple}}(undef, length(selected))
        scenario_logs = Vector{Vector{String}}(undef, length(selected))
        process_tasks = Tuple{Int, BenchmarkCase}[]
        thread_indices = Int[]
        serial_indices = Int[]

        for (idx, base_case) in enumerate(selected)
            route = auto_backend_for_case(base_case)
            if route == :process
                push!(process_tasks, (idx, base_case))
            elseif route == :threads
                push!(thread_indices, idx)
            else
                push!(serial_indices, idx)
            end
        end

        if !isempty(process_tasks)
            ensure_perf_workers!()
            process_results = pmap(process_tasks) do task
                idx, base_case = task
                local_rows, local_logs = perf_worker_measure_per_orbit_scenario(base_case, spec, period_s, orbit_counts, :process)
                return (rows=local_rows, logs=local_logs)
            end
            for (k, task) in enumerate(process_tasks)
                idx = task[1]
                scenario_rows[idx] = process_results[k].rows
                scenario_logs[idx] = process_results[k].logs
            end
        end

        if !isempty(thread_indices)
            for (env_pairs, payload) in _thread_plan_groups(selected, thread_indices, :threads)
                withenv(env_pairs...) do
                    Threads.@threads for j in eachindex(payload)
                        idx, _ = payload[j]
                        local_rows, local_logs = measure_per_orbit_scenario(
                            selected[idx],
                            spec,
                            period_s,
                            orbit_counts;
                            outer_route=:threads,
                            apply_env=false
                        )
                        scenario_rows[idx] = local_rows
                        scenario_logs[idx] = local_logs
                    end
                end
            end
        end

        for idx in serial_indices
            local_rows, local_logs = measure_per_orbit_scenario(selected[idx], spec, period_s, orbit_counts; outer_route=:none)
            scenario_rows[idx] = local_rows
            scenario_logs[idx] = local_logs
        end

        for (idx, base_case) in enumerate(selected)
            println("  scenario $(idx)/$(length(selected)) = $(base_case.name)")
            append!(rows, scenario_rows[idx])
            for log_line in scenario_logs[idx]
                println(log_line)
            end
        end
    elseif backend == :process
        ensure_perf_workers!()
        tasks = collect(enumerate(selected))
        scenario_results = pmap(tasks) do task
            idx, base_case = task
            local_rows, local_logs = perf_worker_measure_per_orbit_scenario(base_case, spec, period_s, orbit_counts, :process)
            return (rows=local_rows, logs=local_logs)
        end
        for (idx, base_case) in enumerate(selected)
            println("  scenario $(idx)/$(length(selected)) = $(base_case.name)")
            append!(rows, scenario_results[idx].rows)
            for log_line in scenario_results[idx].logs
                println(log_line)
            end
        end
    elseif backend == :threads
        scenario_rows = Vector{Vector{NamedTuple}}(undef, length(selected))
        scenario_logs = Vector{Vector{String}}(undef, length(selected))
        thread_indices = collect(eachindex(selected))
        for (env_pairs, payload) in _thread_plan_groups(selected, thread_indices, :threads)
            withenv(env_pairs...) do
                Threads.@threads for j in eachindex(payload)
                    idx, _ = payload[j]
                    local_rows, local_logs = measure_per_orbit_scenario(
                        selected[idx],
                        spec,
                        period_s,
                        orbit_counts;
                        outer_route=:threads,
                        apply_env=false
                    )
                    scenario_rows[idx] = local_rows
                    scenario_logs[idx] = local_logs
                end
            end
        end
        for (idx, base_case) in enumerate(selected)
            println("  scenario $(idx)/$(length(selected)) = $(base_case.name)")
            append!(rows, scenario_rows[idx])
            for log_line in scenario_logs[idx]
                println(log_line)
            end
        end
    else
        for (idx, base_case) in enumerate(selected)
            println("  scenario $(idx)/$(length(selected)) = $(base_case.name)")
            local_rows, local_logs = measure_per_orbit_scenario(base_case, spec, period_s, orbit_counts; outer_route=:none)
            append!(rows, local_rows)
            for log_line in local_logs
                println(log_line)
            end
        end
    end

    run_montecarlo_per_orbit!(rows, spec, planet, period_s, orbit_counts)
    return DataFrame(rows)
end
function _safe_stat(values, op::Function)
    vec = collect(skipmissing(values))
    isempty(vec) && return missing
    return op(vec)
end

function summarize_per_orbit_results(orbit_raw_df::DataFrame)::DataFrame
    keys = [:category, :scenario, :description, :orbit_count, :orbital_period_s]
    counts = combine(
        groupby(orbit_raw_df, keys),
        nrow => :samples_total,
        :solve_success => (v -> count(identity, v)) => :samples_success
    )
    counts[!, :samples_failed] = counts.samples_total .- counts.samples_success
    counts[!, :success_rate] = Float64.(counts.samples_success) ./ Float64.(counts.samples_total)

    success_df = orbit_raw_df[orbit_raw_df.solve_success .== true, :]
    summary = counts
    if nrow(success_df) > 0
        success_summary = combine(
            groupby(success_df, keys),
            nrow => :samples,
            :mission_time_s => (v -> _safe_stat(v, mean)) => :mission_time_mean_s,
            :total_time_s => (v -> _safe_stat(v, mean)) => :total_time_mean_s,
            :total_time_s => (v -> _safe_stat(v, x -> quantile(x, 0.9))) => :total_time_p90_s,
            :solve_time_s => (v -> _safe_stat(v, mean)) => :solve_time_mean_s,
            :total_bytes_mb => (v -> _safe_stat(v, mean)) => :total_bytes_mean_mb,
            :sim_seconds_per_wall_second => (v -> _safe_stat(v, mean)) => :sim_seconds_per_wall_second_mean
        )
        summary = leftjoin(counts, success_summary, on=keys)
    else
        summary[!, :samples] = fill(missing, nrow(summary))
        summary[!, :mission_time_mean_s] = fill(missing, nrow(summary))
        summary[!, :total_time_mean_s] = fill(missing, nrow(summary))
        summary[!, :total_time_p90_s] = fill(missing, nrow(summary))
        summary[!, :solve_time_mean_s] = fill(missing, nrow(summary))
        summary[!, :total_bytes_mean_mb] = fill(missing, nrow(summary))
        summary[!, :sim_seconds_per_wall_second_mean] = fill(missing, nrow(summary))
    end

    summary[!, :time_per_orbit_mean_s] = [
        (ismissing(tt) || orbit <= 0) ? missing : tt / orbit
        for (tt, orbit) in zip(summary.total_time_mean_s, summary.orbit_count)
    ]
    summary[!, :orbits_per_wall_second_mean] = [
        (ismissing(tt) || tt <= 0.0) ? missing : orbit / tt
        for (tt, orbit) in zip(summary.total_time_mean_s, summary.orbit_count)
    ]

    sort!(summary, [:scenario, :orbit_count])
    return summary
end

function summarize_results(raw_df::DataFrame)::DataFrame
    keys = [:category, :scenario, :description, :satellites, :orientation, :mission_time_s, :dynamic_effectors, :control_effectors]
    counts = combine(
        groupby(raw_df, keys),
        nrow => :samples_total,
        :solve_success => (v -> count(identity, v)) => :samples_success
    )
    counts[!, :samples_failed] = counts.samples_total .- counts.samples_success
    counts[!, :success_rate] = Float64.(counts.samples_success) ./ Float64.(counts.samples_total)

    success_df = raw_df[raw_df.solve_success .== true, :]
    metric_cols = [
        :samples,
        :copy_time_mean_s,
        :solve_time_mean_s,
        :total_time_mean_s,
        :total_time_median_s,
        :total_time_std_s,
        :total_time_min_s,
        :total_time_max_s,
        :total_time_p90_s,
        :total_bytes_mean_mb,
        :copy_alloc_mean,
        :solve_alloc_mean,
        :saved_points_mean,
        :accepted_steps_mean,
        :rejected_steps_mean,
        :solver_modes,
        :solver_sequences,
        :solver_fallback_any,
        :solver_fallback_count_mean,
        :solver_fallback_triggers,
        :sim_seconds_per_wall_second_mean,
        :satellite_sim_seconds_per_wall_second_mean
    ]

    summary = counts
    if nrow(success_df) > 0
        success_summary = combine(
            groupby(success_df, keys),
            nrow => :samples,
            :copy_time_s => (v -> _safe_stat(v, mean)) => :copy_time_mean_s,
            :solve_time_s => (v -> _safe_stat(v, mean)) => :solve_time_mean_s,
            :total_time_s => (v -> _safe_stat(v, mean)) => :total_time_mean_s,
            :total_time_s => (v -> _safe_stat(v, median)) => :total_time_median_s,
            :total_time_s => (v -> _safe_stat(v, x -> std(x; corrected=false))) => :total_time_std_s,
            :total_time_s => (v -> _safe_stat(v, minimum)) => :total_time_min_s,
            :total_time_s => (v -> _safe_stat(v, maximum)) => :total_time_max_s,
            :total_time_s => (v -> _safe_stat(v, x -> quantile(x, 0.9))) => :total_time_p90_s,
            :total_bytes_mb => (v -> _safe_stat(v, mean)) => :total_bytes_mean_mb,
            :copy_alloc_calls => (v -> _safe_stat(v, mean)) => :copy_alloc_mean,
            :solve_alloc_calls => (v -> _safe_stat(v, mean)) => :solve_alloc_mean,
            :saved_points => (v -> _safe_stat(v, mean)) => :saved_points_mean,
            :accepted_steps => (v -> _safe_stat(v, mean)) => :accepted_steps_mean,
            :rejected_steps => (v -> _safe_stat(v, mean)) => :rejected_steps_mean,
            :solver_mode => (v -> _safe_unique_join(v)) => :solver_modes,
            :solver_sequence => (v -> _safe_unique_join(v)) => :solver_sequences,
            :solver_fallback_used => (v -> any(skipmissing(v))) => :solver_fallback_any,
            :solver_fallback_count => (v -> _safe_stat(v, mean)) => :solver_fallback_count_mean,
            :solver_fallback_trigger => (v -> _safe_unique_join(v; delimiter="|")) => :solver_fallback_triggers,
            :sim_seconds_per_wall_second => (v -> _safe_stat(v, mean)) => :sim_seconds_per_wall_second_mean,
            :satellite_sim_seconds_per_wall_second => (v -> _safe_stat(v, mean)) => :satellite_sim_seconds_per_wall_second_mean
        )
        summary = leftjoin(counts, success_summary, on=keys)
    else
        for col in metric_cols
            summary[!, col] = fill(missing, nrow(summary))
        end
    end

    baseline_idx = findfirst(==("single_baseline_gravity"), summary.scenario)
    if baseline_idx === nothing || ismissing(summary.total_time_mean_s[baseline_idx]) || summary.total_time_mean_s[baseline_idx] <= 0.0
        summary[!, :relative_to_baseline] = fill(missing, nrow(summary))
        summary[!, :speedup_vs_baseline] = fill(missing, nrow(summary))
    else
        baseline = summary.total_time_mean_s[baseline_idx]
        summary[!, :relative_to_baseline] = [
            ismissing(row_total) ? missing : row_total / baseline
            for row_total in summary.total_time_mean_s
        ]
        summary[!, :speedup_vs_baseline] = [
            ismissing(row_total) || row_total <= 0.0 ? missing : baseline / row_total
            for row_total in summary.total_time_mean_s
        ]
    end

    summary[!, :_sort_key] = [ismissing(x) ? -Inf : x for x in summary.total_time_mean_s]
    sort!(summary, :_sort_key, rev=true)
    select!(summary, Not(:_sort_key))
    return summary
end

@inline _plot_label(name::AbstractString) = replace(name, "_" => " ")

@inline function _plot_ready()::Bool
    return myid() == 1 && isdefined(Main, :Plots)
end

function _plot_metric_pairs(df::DataFrame, label_col::Symbol, metric_col::Symbol)
    labels = String[]
    values = Float64[]
    for row in eachrow(df)
        value = row[metric_col]
        if !(value isa Missing)
            f = Float64(value)
            if isfinite(f)
                push!(labels, _plot_label(string(row[label_col])))
                push!(values, f)
            end
        end
    end
    return labels, values
end

@inline function _has_row_fields(row, fields::Vector{Symbol})::Bool
    for field in fields
        if row[field] isa Missing
            return false
        end
    end
    return true
end

function _save_runtime_plot!(
    artifacts::Vector{String},
    plt,
    outdir::String,
    basename::String,
    spec::ProfileSpec,
    stamp::String
)
    path = joinpath(outdir, "$(basename)_$(spec.name)_$(stamp).png")
    try
        Plots.savefig(plt, path)
        push!(artifacts, path)
    catch err
        @warn "[perf] failed to save plot $basename: $(_perf_error_text(err))"
    end
    return nothing
end

function generate_runtime_plots(
    outdir::String,
    spec::ProfileSpec,
    stamp::String,
    raw_df::DataFrame,
    summary_df::DataFrame,
    orbit_summary_df::DataFrame
)::Vector{String}
    !_plot_ready() && return String[]

    plot_artifacts = String[]
    default_size = (1280, 720)
    success_summary = summary_df[summary_df.samples_success .> 0, :]

    # 1) Mean total runtime per scenario.
    labels_totals, values_totals = _plot_metric_pairs(success_summary, :scenario, :total_time_mean_s)
    if !isempty(values_totals)
        plt = Plots.bar(
            labels_totals,
            values_totals;
            legend=false,
            title="Mean Total Runtime by Scenario",
            xlabel="Scenario",
            ylabel="Seconds",
            xrotation=30,
            size=default_size
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_totals", spec, stamp)
    end

    # 2) Speedup against single_baseline_gravity.
    labels_speedup, values_speedup = _plot_metric_pairs(success_summary, :scenario, :speedup_vs_baseline)
    if !isempty(values_speedup)
        plt = Plots.bar(
            labels_speedup,
            values_speedup;
            legend=false,
            title="Speedup vs single_baseline_gravity",
            xlabel="Scenario",
            ylabel="Speedup (x)",
            xrotation=30,
            size=default_size
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_speedup", spec, stamp)
    end

    # 3) Runtime mean/std variability map.
    variability_df = success_summary[
        .!ismissing.(success_summary.total_time_mean_s) .&
        .!ismissing.(success_summary.total_time_std_s), :
    ]
    if nrow(variability_df) > 0
        x = Float64.(variability_df.total_time_mean_s)
        y = Float64.(variability_df.total_time_std_s)
        plt = Plots.scatter(
            x,
            y;
            legend=false,
            title="Runtime Variability Map",
            xlabel="Mean total time (s)",
            ylabel="Std total time (s)",
            size=default_size
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_variability", spec, stamp)
    end

    # 4) Copy vs solve breakdown.
    labels_breakdown = String[]
    copy_vals = Float64[]
    solve_vals = Float64[]
    for row in eachrow(success_summary)
        if _has_row_fields(row, [:copy_time_mean_s, :solve_time_mean_s])
            push!(labels_breakdown, _plot_label(row.scenario))
            push!(copy_vals, Float64(row.copy_time_mean_s))
            push!(solve_vals, Float64(row.solve_time_mean_s))
        end
    end
    if !isempty(labels_breakdown)
        plt = Plots.bar(
            labels_breakdown,
            hcat(copy_vals, solve_vals);
            label=["copy" "solve"],
            title="Copy vs Solve Runtime Breakdown",
            xlabel="Scenario",
            ylabel="Seconds",
            xrotation=30,
            size=default_size
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_breakdown_copy_solve", spec, stamp)
    end

    # 5) Memory and allocations proxies.
    mem_df = success_summary[
        .!ismissing.(success_summary.total_bytes_mean_mb) .&
        .!ismissing.(success_summary.solve_alloc_mean), :
    ]
    if nrow(mem_df) > 0
        x = Float64.(mem_df.total_bytes_mean_mb)
        y = Float64.(mem_df.solve_alloc_mean)
        plt = Plots.scatter(
            x,
            y;
            legend=false,
            title="Memory vs Allocation Calls",
            xlabel="Mean allocated bytes (MB)",
            ylabel="Mean solve allocation calls",
            size=default_size
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_memory_alloc", spec, stamp)
    end

    # 6) Solver accepted/rejected steps.
    labels_solver = String[]
    accepted_vals = Float64[]
    rejected_vals = Float64[]
    for row in eachrow(success_summary)
        if _has_row_fields(row, [:accepted_steps_mean, :rejected_steps_mean])
            push!(labels_solver, _plot_label(row.scenario))
            push!(accepted_vals, Float64(row.accepted_steps_mean))
            push!(rejected_vals, Float64(row.rejected_steps_mean))
        end
    end
    if !isempty(labels_solver)
        plt = Plots.bar(
            labels_solver,
            hcat(accepted_vals, rejected_vals);
            label=["accepted" "rejected"],
            title="Solver Workload (Accepted/Rejected Steps)",
            xlabel="Scenario",
            ylabel="Step count",
            xrotation=30,
            size=default_size
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_solver_workload", spec, stamp)
    end

    # 7) Throughput by scenario.
    labels_throughput, values_throughput = _plot_metric_pairs(success_summary, :scenario, :sim_seconds_per_wall_second_mean)
    if !isempty(values_throughput)
        plt = Plots.bar(
            labels_throughput,
            values_throughput;
            legend=false,
            title="Throughput by Scenario",
            xlabel="Scenario",
            ylabel="Sim seconds / wall second",
            xrotation=30,
            size=default_size
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_throughput", spec, stamp)
    end

    # 8) Satellite scaling curve.
    sat_df = success_summary[
        (success_summary.category .== "satellite_scaling") .&
        .!ismissing.(success_summary.satellites) .&
        .!ismissing.(success_summary.total_time_mean_s), :
    ]
    if nrow(sat_df) > 0
        sort!(sat_df, :satellites)
        x = Int.(sat_df.satellites)
        y = Float64.(sat_df.total_time_mean_s)
        plt = Plots.plot(
            x,
            y;
            marker=:circle,
            legend=false,
            title="Satellite Scaling",
            xlabel="Satellites",
            ylabel="Mean total time (s)",
            size=default_size
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_satellite_scaling", spec, stamp)
    end

    # 9) Dynamics fidelity ladder.
    fidelity_order = [
        "single_baseline_gravity",
        "single_j2",
        "single_nbody_sun_moon",
        "single_harmonics_l20",
        "single_harmonics_l50",
    ]
    fidelity_labels = String[]
    fidelity_values = Float64[]
    for scenario in fidelity_order
        idx = findfirst(==(scenario), success_summary.scenario)
        if idx !== nothing
            value = success_summary.total_time_mean_s[idx]
            if !(value isa Missing)
                push!(fidelity_labels, _plot_label(scenario))
                push!(fidelity_values, Float64(value))
            end
        end
    end
    if length(fidelity_values) >= 2
        plt = Plots.plot(
            fidelity_labels,
            fidelity_values;
            marker=:circle,
            legend=false,
            title="Dynamics Fidelity Ladder",
            xlabel="Scenario",
            ylabel="Mean total time (s)",
            xrotation=20,
            size=default_size
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_fidelity_ladder", spec, stamp)
    end

    # 10) Monte Carlo runtime histogram.
    mc_df = raw_df[(raw_df.category .== "montecarlo") .& (raw_df.solve_success .== true), :]
    mc_times = [Float64(v) for v in mc_df.total_time_s if !(v isa Missing)]
    if !isempty(mc_times)
        plt = Plots.histogram(
            mc_times;
            bins=min(20, max(5, round(Int, sqrt(length(mc_times))))),
            legend=false,
            title="Monte Carlo Runtime Distribution",
            xlabel="Total time (s)",
            ylabel="Count",
            size=default_size
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_montecarlo_hist", spec, stamp)
    end

    # 11) Monte Carlo seed trace.
    mc_seed_df = mc_df[.!ismissing.(mc_df.seed), :]
    if nrow(mc_seed_df) > 0
        sort!(mc_seed_df, :seed)
        seeds = Int.(mc_seed_df.seed)
        totals = Float64.(mc_seed_df.total_time_s)
        plt = Plots.plot(
            seeds,
            totals;
            marker=:circle,
            legend=false,
            title="Monte Carlo Runtime by Seed",
            xlabel="Seed",
            ylabel="Total time (s)",
            size=default_size
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_montecarlo_seed_trace", spec, stamp)
    end

    orbit_valid = orbit_summary_df[
        (orbit_summary_df.samples_success .> 0) .&
        .!ismissing.(orbit_summary_df.orbit_count), :
    ]

    # 12) Per-orbit total scaling curves by scenario.
    orbit_scaling_df = orbit_valid[.!ismissing.(orbit_valid.total_time_mean_s), :]
    if nrow(orbit_scaling_df) > 0
        plt = Plots.plot(
            title="Per-Orbit Total Runtime Scaling",
            xlabel="Orbit count",
            ylabel="Mean total time (s)",
            size=default_size
        )
        for grp in groupby(orbit_scaling_df, :scenario)
            local_df = DataFrame(grp)
            sort!(local_df, :orbit_count)
            x = Int.(local_df.orbit_count)
            y = Float64.(local_df.total_time_mean_s)
            Plots.plot!(plt, x, y; marker=:circle, label=_plot_label(local_df.scenario[1]))
        end
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_per_orbit_scaling", spec, stamp)
    end

    # 13) Per-orbit efficiency curves by scenario.
    orbit_eff_df = orbit_valid[.!ismissing.(orbit_valid.orbits_per_wall_second_mean), :]
    if nrow(orbit_eff_df) > 0
        plt = Plots.plot(
            title="Per-Orbit Efficiency",
            xlabel="Orbit count",
            ylabel="Orbits per wall second",
            size=default_size
        )
        for grp in groupby(orbit_eff_df, :scenario)
            local_df = DataFrame(grp)
            sort!(local_df, :orbit_count)
            x = Int.(local_df.orbit_count)
            y = Float64.(local_df.orbits_per_wall_second_mean)
            Plots.plot!(plt, x, y; marker=:circle, label=_plot_label(local_df.scenario[1]))
        end
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_per_orbit_efficiency", spec, stamp)
    end

    # 14) Per-orbit time-per-orbit heatmap.
    heat_df = orbit_valid[.!ismissing.(orbit_valid.time_per_orbit_mean_s), :]
    if nrow(heat_df) > 0
        scenarios = sort(unique(String.(heat_df.scenario)))
        orbits = sort(unique(Int.(heat_df.orbit_count)))
        z = fill(NaN, length(scenarios), length(orbits))
        for row in eachrow(heat_df)
            si = findfirst(==(String(row.scenario)), scenarios)
            oi = findfirst(==(Int(row.orbit_count)), orbits)
            if si !== nothing && oi !== nothing
                z[si, oi] = Float64(row.time_per_orbit_mean_s)
            end
        end
        plt = Plots.heatmap(
            orbits,
            1:length(scenarios),
            z;
            colorbar_title="s/orbit",
            yticks=(1:length(scenarios), _plot_label.(scenarios)),
            xlabel="Orbit count",
            ylabel="Scenario",
            title="Per-Orbit Time Heatmap",
            size=default_size
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_per_orbit_heatmap", spec, stamp)
    end

    return plot_artifacts
end

@inline function _fmt(v; digits::Int=3)
    if v isa Missing
        return "n/a"
    elseif v isa AbstractFloat
        return isfinite(v) ? string(round(v; digits=digits)) : "n/a"
    else
        return string(v)
    end
end

function _scenario_metric(summary_df::DataFrame, scenario::String, metric::Symbol)
    idx = findfirst(==(scenario), summary_df.scenario)
    idx === nothing && return nothing
    return summary_df[idx, metric]
end

function write_report(
    path::String,
    spec::ProfileSpec,
    raw_df::DataFrame,
    summary_df::DataFrame,
    orbit_summary_df::DataFrame;
    plot_paths::Vector{String}=String[]
)
    generated = string(now(UTC))
    julia_ver = string(VERSION)
    nthreads = Threads.nthreads()

    valid_rows = findall(x -> !ismissing(x), summary_df.total_time_mean_s)
    fastest = nothing
    slowest = nothing
    if !isempty(valid_rows)
        vals = summary_df.total_time_mean_s[valid_rows]
        fastest = summary_df[valid_rows[argmin(vals)], :]
        slowest = summary_df[valid_rows[argmax(vals)], :]
    end

    baseline = _scenario_metric(summary_df, "single_baseline_gravity", :total_time_mean_s)
    orientation = _scenario_metric(summary_df, "single_orientation_aero", :total_time_mean_s)
    multi2 = _scenario_metric(summary_df, "multi_2_gravity", :total_time_mean_s)
    multi4 = _scenario_metric(summary_df, "multi_4_gravity", :total_time_mean_s)
    harmonics20 = _scenario_metric(summary_df, "single_harmonics_l20", :total_time_mean_s)
    harmonics50 = _scenario_metric(summary_df, "single_harmonics_l50", :total_time_mean_s)

    mc_rows = raw_df[(raw_df.category .== "montecarlo") .& (raw_df.solve_success .== true), :]
    mc_mean = nrow(mc_rows) > 0 ? mean(mc_rows.total_time_s) : missing
    mc_p90 = nrow(mc_rows) > 0 ? quantile(mc_rows.total_time_s, 0.9) : missing
    mc_std = nrow(mc_rows) > 0 ? std(mc_rows.total_time_s; corrected=false) : missing
    fallback_rows = raw_df[(raw_df.solve_success .== true) .& (raw_df.solver_fallback_used .== true), :]
    solver_modes = _safe_unique_join(raw_df.solver_mode)

    open(path, "w") do io
        println(io, "# SpaceAGORA Computational Time Analysis (`$(spec.name)` profile)")
        println(io)
        println(io, "- Generated (UTC): $generated")
        println(io, "- Julia: `$julia_ver`")
        println(io, "- Threads: `$nthreads`")
        println(io, "- Repeats per deterministic scenario: `$(spec.repeats)`")
        println(io, "- Warmup runs per scenario: `$(spec.warmup)`")
        println(io, "- Monte Carlo seeds: `$(spec.montecarlo_samples)`")
        println(io, "- Solver mode(s) observed: `$(_fmt(solver_modes))`")
        println(io)
        println(io, "## Key Findings")
        println(io)
        if fastest === nothing || slowest === nothing
            println(io, "- No successful runs were recorded.")
        else
            println(io, "- Slowest successful scenario: `$(slowest.scenario)` with mean total time `$(round(slowest.total_time_mean_s; digits=3)) s`.")
            println(io, "- Fastest successful scenario: `$(fastest.scenario)` with mean total time `$(round(fastest.total_time_mean_s; digits=3)) s`.")
        end
        if baseline !== nothing && orientation !== nothing && !ismissing(baseline) && !ismissing(orientation) && baseline > 0.0
            println(io, "- Orientation + aerodynamic run vs baseline: `$(round(orientation / baseline; digits=2))x` runtime.")
        end
        if baseline !== nothing && multi2 !== nothing && multi4 !== nothing && !ismissing(baseline) && !ismissing(multi2) && !ismissing(multi4) && baseline > 0.0
            println(io, "- Multi-satellite scaling: `2-sat=$(round(multi2 / baseline; digits=2))x`, `4-sat=$(round(multi4 / baseline; digits=2))x` relative to single-sat baseline.")
        end
        if harmonics20 !== nothing && harmonics50 !== nothing && !ismissing(harmonics20) && !ismissing(harmonics50) && harmonics20 > 0.0
            println(io, "- Harmonics scaling: `L=50` is `$(round(harmonics50 / harmonics20; digits=2))x` relative to `L=20`.")
        end
        if !(mc_mean isa Missing)
            println(io, "- Monte Carlo runtime spread: mean `$(round(mc_mean; digits=3)) s`, p90 `$(round(mc_p90; digits=3)) s`, std `$(round(mc_std; digits=3)) s`.")
        end
        println(io, "- Auto-stiff fallback activations (successful runs): `$(nrow(fallback_rows))`.")
        println(io, "- Plot artifacts generated: `$(length(plot_paths))`.")
        failed_groups = summary_df[summary_df.samples_failed .> 0, :]
        if nrow(failed_groups) > 0
            println(io, "- Solver failures detected in `$(nrow(failed_groups))` scenario groups; timings only use successful runs.")
        end
        if !isempty(plot_paths)
            println(io)
            println(io, "## Plot Artifacts")
            println(io)
            for plot_path in plot_paths
                println(io, "- `$(plot_path)`")
            end
            println(io)
        end
        println(io)
        println(io, "## Scenario Summary")
        println(io)
        println(io, "| Scenario | Category | Success/Total | Solver(s) | Fallback Any | Mean Total (s) | P90 (s) | Mean Solve (s) | Mean Copy (s) | Sim sec / wall sec | Rel. Baseline |")
        println(io, "|---|---|---:|---|---|---:|---:|---:|---:|---:|---:|")
        for row in eachrow(summary_df)
            println(
                io,
                "| $(row.scenario) | $(row.category) | $(row.samples_success)/$(row.samples_total) | $(_fmt(row.solver_sequences)) | $(_fmt(row.solver_fallback_any)) | $(_fmt(row.total_time_mean_s)) | $(_fmt(row.total_time_p90_s)) | $(_fmt(row.solve_time_mean_s)) | $(_fmt(row.copy_time_mean_s)) | $(_fmt(row.sim_seconds_per_wall_second_mean)) | $(_fmt(row.relative_to_baseline)) |"
            )
        end
        println(io)
        println(io, "## Per-Orbit Results (All Scenarios)")
        println(io)
        println(io, "| Scenario | Category | Orbit Count | Success/Total | Mission Time (s) | Mean Total (s) | P90 (s) | Time / Orbit (s) | Orbits / Wall-sec |")
        println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|")
        for row in eachrow(orbit_summary_df)
            println(
                io,
                "| $(row.scenario) | $(row.category) | $(row.orbit_count) | $(row.samples_success)/$(row.samples_total) | $(_fmt(row.mission_time_mean_s)) | $(_fmt(row.total_time_mean_s)) | $(_fmt(row.total_time_p90_s)) | $(_fmt(row.time_per_orbit_mean_s)) | $(_fmt(row.orbits_per_wall_second_mean)) |"
            )
        end
    end
end

function main()
    spec, outdir = parse_cli()
    mkpath(outdir)

    println("Performance runtime analysis profile: $(spec.name)")
    println("Output directory: $outdir")

    planet = Earth("", SPICE_PATH)
    cases = build_cases(spec, planet)
    raw_df = run_benchmarks(spec, cases, planet)
    summary_df = summarize_results(raw_df)
    orbit_raw_df = run_per_orbit_for_scenarios(spec, cases, planet)
    orbit_summary_df = summarize_per_orbit_results(orbit_raw_df)

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    raw_path = joinpath(outdir, "runtime_raw_$(spec.name)_$(stamp).csv")
    summary_path = joinpath(outdir, "runtime_summary_$(spec.name)_$(stamp).csv")
    orbit_raw_path = joinpath(outdir, "runtime_per_orbit_raw_$(spec.name)_$(stamp).csv")
    orbit_summary_path = joinpath(outdir, "runtime_per_orbit_summary_$(spec.name)_$(stamp).csv")
    report_path = joinpath(outdir, "runtime_report_$(spec.name)_$(stamp).md")
    plot_paths = generate_runtime_plots(outdir, spec, stamp, raw_df, summary_df, orbit_summary_df)

    CSV.write(raw_path, raw_df)
    CSV.write(summary_path, summary_df)
    CSV.write(orbit_raw_path, orbit_raw_df)
    CSV.write(orbit_summary_path, orbit_summary_df)
    write_report(report_path, spec, raw_df, summary_df, orbit_summary_df; plot_paths=plot_paths)

    println("Analysis complete.")
    println("Raw results: $raw_path")
    println("Summary: $summary_path")
    println("Per-orbit raw: $orbit_raw_path")
    println("Per-orbit summary: $orbit_summary_path")
    println("Plots generated: $(length(plot_paths))")
    println("Report: $report_path")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
