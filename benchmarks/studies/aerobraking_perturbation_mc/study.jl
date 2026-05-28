module AerobrakingPerturbationMC

using Arrow
using CSV
using DataFrames
using Dates
using DiffEqCallbacks: PresetTimeCallback
using Distributed
using DiffEqBase: DiscreteCallback, terminate!
using LinearAlgebra
using Printf
using Statistics
using StaticArrays
using TOML

include(joinpath(@__DIR__, "..", "..", "..", "examples", "common.jl"))

import SpaceAGORA: getDensity, getDensityBatch!

const DEFAULT_PLANETS = [:mars, :earth, :venus, :titan]
const DEFAULT_PERIAPSIS_REGIMES = [:shallow, :nominal, :deep]
const MIN_APOAPSIS_ALTITUDE_KM = 1_000.0
const DEFAULT_APOAPSIS_ALTITUDES_KM = Dict{Symbol, Vector{Float64}}(
    :mars  => [1_000.0, 1_250.0, 1_500.0, 1_750.0, 2_000.0, 2_500.0, 3_000.0, 4_000.0, 5_000.0, 7_500.0, 10_000.0, 12_000.0, 15_000.0, 20_000.0, 30_000.0],
    :earth => [1_000.0, 1_500.0, 2_000.0, 3_000.0, 5_000.0, 7_500.0, 10_000.0, 15_000.0, 20_000.0, 30_000.0, 36_000.0, 45_000.0, 60_000.0],
    :venus => [1_000.0, 1_500.0, 2_000.0, 3_000.0, 5_000.0, 7_500.0, 10_000.0, 15_000.0, 20_000.0, 30_000.0, 40_000.0],
    :titan => [1_000.0, 1_500.0, 2_000.0, 3_000.0, 5_000.0, 7_500.0, 10_000.0, 15_000.0, 20_000.0, 30_000.0, 40_000.0],
)
const DEFAULT_DYNAMICS_CASES = [
    :two_body,
    :j2,
    :harmonics_low,
    :srp,
    :third_body_sun,
    :gram_aero,
    :full_environment,
]
const DEFAULT_DENSITY_SCALES = Dict(:low => 0.9, :nominal => 1.0, :high => 1.1)
const OUTSIDE_ATMOSPHERE_SAVE_DT_S = 100.0
const INSIDE_ATMOSPHERE_SAVE_DT_S = 5.0
const DEFAULT_WORKER_MAX_CASES = 16
const DEFAULT_WORKER_MAX_RSS_GB = 6.0
const DEFAULT_AERO_SOLVER = :auto_stiff_then_rodas5p
const DEFAULT_DT_MAX_ORBIT_S = OUTSIDE_ATMOSPHERE_SAVE_DT_S
const DEFAULT_DT_MAX_ATMOSPHERE_S = INSIDE_ATMOSPHERE_SAVE_DT_S
const _WORKER_POOL_LOCK = ReentrantLock()

@inline _aero_remote_workers() = [w for w in workers() if w != myid()]

const NOMINAL_SPACECRAFT_MASS_SCALE = 1.0
const NOMINAL_INCLINATION_DEG = 93.0
const NOMINAL_ARGP_DEG = 80.0

const DEFAULT_SPACECRAFT_MASS_SCALES = [0.25, 0.5, 1.0, 2.0, 4.0]
const DEFAULT_INCLINATIONS_DEG = [15.0, 30.0, 45.0, 60.0, 75.0, 93.0, 105.0, 120.0, 135.0, 150.0, 165.0]
const DEFAULT_ARGP_DEGS = [0.0, 30.0, 45.0, 60.0, 80.0, 90.0, 120.0, 135.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0]
const DEFAULT_AERO_AUTO_MAXITERS = 500_000
const DEFAULT_AERO_STIFF_MAXITERS = 20_000_000
const DEFAULT_AERO_AUTO_STIFF_SWITCH_MAX = 5
const DEFAULT_AERO_STIFF_TOL_SCALE = 10.0
const DEFAULT_CASE_TIMEOUT_MIN = 30.0
const DEFAULT_DEORBIT_BAILOUT = true
const DEFAULT_DEORBIT_BAILOUT_MARGIN_KM = 25.0
const DEFAULT_DEORBIT_BAILOUT_CHECK_DT_S = 30.0

@inline _timestamp() = Dates.format(now(), "HH:MM:SS")

@inline function _sample_field(sample::NamedTuple, key::Symbol, default)
    return haskey(sample, key) ? getfield(sample, key) : default
end

function _sample_label(sample::NamedTuple)::String
    density = _sample_field(sample, :density_case, :none)
    density_part = density === :none ? "" : "/density=$(density)"
    apo_str = @sprintf("%.0fkm", sample.apoapsis_alt_km)
    return "$(sample.planet)/$(sample.periapsis_regime)/apo=$(apo_str)/ms=$(sample.spacecraft_mass_scale)/inc=$(sample.inclination_deg)/aop=$(sample.argp_deg)/$(sample.dynamics_case)$(density_part)"
end

function _progress_prefix(sample::NamedTuple)::String
    idx = _sample_field(sample, :sample_index, 0)
    total = _sample_field(sample, :sample_count, 0)
    counter = idx > 0 && total > 0 ? "[$idx/$total] " : ""
    return "[aero-perturb $(_timestamp())] $(counter)worker=$(myid())"
end

function _sample_slug(sample::NamedTuple)::String
    idx = _sample_field(sample, :sample_index, 0)
    prefix = idx > 0 ? @sprintf("case_%03d", idx) : "case"
    density = _sample_field(sample, :density_case, :none)
    ms_str = replace(@sprintf("%.2f", sample.spacecraft_mass_scale), "." => "p")
    return join((
        prefix,
        String(sample.planet),
        String(sample.periapsis_regime),
        @sprintf("apo_%05.0fkm", sample.apoapsis_alt_km),
        String(sample.dynamics_case),
        "density_" * String(density),
        "ms" * ms_str,
        @sprintf("inc%03.0f", sample.inclination_deg),
        @sprintf("aop%03.0f", sample.argp_deg),
    ), "_")
end

function _log_sample(sample::NamedTuple, message::AbstractString)
    println("$(_progress_prefix(sample)) $message")
    flush(stdout)
    return nothing
end

function _arrow_safe_dataframe(df::DataFrame)::DataFrame
    out = copy(df)
    for name in names(out)
        col = out[!, name]
        if any(value -> value isa Symbol, col)
            out[!, name] = [
                value === missing ? missing : value isa Symbol ? String(value) : value
                for value in col
            ]
        end
    end
    return out
end

function _write_feather(path::String, df::DataFrame)::String
    Arrow.write(path, _arrow_safe_dataframe(df))
    return path
end

function _write_feather_table(base_path_without_ext::String, df::DataFrame)
    feather_path = base_path_without_ext * ".feather"
    _write_feather(feather_path, df)
    return (; feather_path)
end

Base.@kwdef struct StudySpec
    planets::Vector{Symbol} = copy(DEFAULT_PLANETS)
    periapsis_regimes::Vector{Symbol} = copy(DEFAULT_PERIAPSIS_REGIMES)
    apoapsis_altitudes_km::Dict{Symbol, Vector{Float64}} = deepcopy(DEFAULT_APOAPSIS_ALTITUDES_KM)
    dynamics_cases::Vector{Symbol} = copy(DEFAULT_DYNAMICS_CASES)
    density_scales::Dict{Symbol, Float64} = copy(DEFAULT_DENSITY_SCALES)
    spacecraft_mass_scales::Vector{Float64} = copy(DEFAULT_SPACECRAFT_MASS_SCALES)
    inclinations_deg::Vector{Float64} = copy(DEFAULT_INCLINATIONS_DEG)
    argp_degs::Vector{Float64} = copy(DEFAULT_ARGP_DEGS)
    nominal_mass_scale::Float64 = NOMINAL_SPACECRAFT_MASS_SCALE
    nominal_inclination_deg::Float64 = NOMINAL_INCLINATION_DEG
    nominal_argp_deg::Float64 = NOMINAL_ARGP_DEG
    norbits::Int = 1
    procs::Int = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_PROCS", "0"))
    worker_max_cases::Int = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_WORKER_MAX_CASES", string(DEFAULT_WORKER_MAX_CASES)))
    worker_max_rss_gb::Float64 = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_WORKER_MAX_RSS_GB", string(DEFAULT_WORKER_MAX_RSS_GB)))
    aero_solver::Symbol = Symbol(get(ENV, "SPACEAGORA_AERO_PERTURB_AERO_SOLVER", string(DEFAULT_AERO_SOLVER)))
    aero_auto_maxiters::Int = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_AERO_AUTO_MAXITERS", string(DEFAULT_AERO_AUTO_MAXITERS)))
    aero_stiff_maxiters::Int = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_AERO_STIFF_MAXITERS", string(DEFAULT_AERO_STIFF_MAXITERS)))
    aero_auto_stiff_switch_max::Int = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_AERO_AUTO_STIFF_SWITCH_MAX", string(DEFAULT_AERO_AUTO_STIFF_SWITCH_MAX)))
    aero_stiff_tol_scale::Float64 = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_AERO_STIFF_TOL_SCALE", string(DEFAULT_AERO_STIFF_TOL_SCALE)))
    case_timeout_min::Float64 = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_CASE_TIMEOUT_MIN", string(DEFAULT_CASE_TIMEOUT_MIN)))
    deorbit_bailout::Bool = get(ENV, "SPACEAGORA_AERO_PERTURB_DEORBIT_BAILOUT", DEFAULT_DEORBIT_BAILOUT ? "1" : "0") in ("1", "true", "TRUE", "yes", "YES", "on", "ON")
    deorbit_bailout_margin_km::Float64 = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_DEORBIT_BAILOUT_MARGIN_KM", string(DEFAULT_DEORBIT_BAILOUT_MARGIN_KM)))
    deorbit_bailout_check_dt_s::Float64 = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_DEORBIT_BAILOUT_CHECK_DT_S", string(DEFAULT_DEORBIT_BAILOUT_CHECK_DT_S)))
    dt_max_orbit_s::Float64 = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_DT_MAX_ORBIT_S", string(DEFAULT_DT_MAX_ORBIT_S)))
    dt_max_atmosphere_s::Float64 = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_DT_MAX_ATMOSPHERE_S", string(DEFAULT_DT_MAX_ATMOSPHERE_S)))
    output_dir::String = joinpath(REPO_ROOT, "output", "aerobraking_perturbation_mc")
    results::Bool = false
    smoke::Bool = get(ENV, "SPACEAGORA_AERO_PERTURB_SMOKE", "0") == "1"
    resume::Bool = false
    resume_dir::String = ""
end

const GRAM_TRAJECTORY_DENSITY_ENV = Dict(
    "SPACEAGORA_GRAM_TRACK_CACHE" => "off",
    "SPACEAGORA_VACUUM_GRAM_CACHE" => "1",
    "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "24",
    "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => "900",
    "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "5000",
)

function configure_gram_trajectory_density!()
    applied = Dict{String, String}()
    for (name, default_value) in GRAM_TRAJECTORY_DENSITY_ENV
        override_name = "SPACEAGORA_AERO_PERTURB_" * replace(name, "SPACEAGORA_" => "")
        value = get(ENV, override_name, default_value)
        ENV[name] = value
        applied[name] = value
    end
    return applied
end

function _aero_worker_init_expr()
    study_dir = @__DIR__
    project = REPO_ROOT
    return quote
        import Pkg
        Pkg.activate($project; io=devnull)
        include(joinpath($study_dir, "study.jl"))
        nothing
    end
end

function _start_aero_perturb_worker!()::Int
    project = REPO_ROOT
    worker = only(addprocs(1; exeflags="--project=$(project) --compiled-modules=existing"))
    try
        remotecall_wait(Core.eval, worker, Main, _aero_worker_init_expr())
        println("[aero-perturb] worker=$worker ready")
        flush(stdout)
        return worker
    catch
        try
            rmprocs(worker; waitfor=10.0)
        catch
        end
        rethrow()
    end
end

function ensure_aero_perturb_workers!(target_workers::Int)
    target_workers <= 0 && return nothing
    lock(_WORKER_POOL_LOCK) do
        missing = target_workers - length(_aero_remote_workers())
        missing <= 0 && return nothing
        println("[aero-perturb] starting $(missing) worker process(es); this can take a bit on first run")
        flush(stdout)
        project = REPO_ROOT
        new_workers = addprocs(missing; exeflags="--project=$(project) --compiled-modules=existing")
        println("[aero-perturb] initializing $(length(new_workers)) new worker(s) in parallel...")
        flush(stdout)
        init_ex = _aero_worker_init_expr()
        init_tasks = map(new_workers) do worker
            @async begin
                try
                    remotecall_wait(Core.eval, worker, Main, init_ex)
                    println("[aero-perturb] worker=$worker ready")
                    flush(stdout)
                catch
                    try
                        rmprocs(worker; waitfor=10.0)
                    catch
                    end
                    rethrow()
                end
            end
        end
        foreach(wait, init_tasks)
    end
    return nothing
end

function _recycle_aero_perturb_worker!(worker::Int; reason::AbstractString="case limit")::Int
    lock(_WORKER_POOL_LOCK) do
        println("[aero-perturb] recycling worker=$worker after $(reason)")
        flush(stdout)
        try
            rmprocs(worker; waitfor=10.0)
        catch err
            @warn "[aero-perturb] failed to remove worker=$worker before recycling" exception=(err, catch_backtrace())
        end
        return _start_aero_perturb_worker!()
    end
end

function _malloc_trim!()
    Sys.islinux() || return nothing
    try
        ccall(:malloc_trim, Cint, (Cint,), 0)
    catch
    end
    return nothing
end

function _post_sample_gc!()
    GC.gc(true)
    _malloc_trim!()
    return nothing
end

function _current_rss_bytes()::Int64
    if Sys.islinux() && isfile("/proc/self/status")
        rss_kb = open("/proc/self/status", "r") do io
            for line in eachline(io)
                if startswith(line, "VmRSS:")
                    parts = split(line)
                    length(parts) >= 2 && return parse(Int64, parts[2])
                end
            end
            return nothing
        end
        rss_kb === nothing || return rss_kb * 1024
    end
    return Int64(Sys.maxrss())
end

_current_rss_mib() = _current_rss_bytes() / 2.0^20

struct DensityScaleModel{M} <: SM.AbstractDensityModel
    base::M
    scale::Float64
end

function Base.getproperty(model::DensityScaleModel, name::Symbol)
    if name === :base || name === :scale
        return getfield(model, name)
    end
    base = getfield(model, :base)
    if hasproperty(base, name)
        return getproperty(base, name)
    end
    return getfield(model, name)
end

SM.SimulationCallbacks._is_gram_density_model(model::DensityScaleModel) =
    SM.SimulationCallbacks._is_gram_density_model(model.base)

SM.SimulationCallbacks._gram_track_trajectory_supported(model::DensityScaleModel) =
    SM.SimulationCallbacks._gram_track_trajectory_supported(model.base)

SM.SimulationCallbacks.density_model_threadsafe(model::DensityScaleModel)::Bool =
    SM.SimulationCallbacks.density_model_threadsafe(model.base)

function getDensity(model::DensityScaleModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p)
    rho, temp, wind_vec = Base.invokelatest(getDensity, model.base, h, lat, lon, el_time, wind, p)
    return model.scale * rho, temp, wind_vec
end

function getDensity(model::DensityScaleModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)
    rho, temp, wind_vec = Base.invokelatest(getDensity, model.base, h, lat, lon, el_time, wind)
    return model.scale * rho, temp, wind_vec
end

function getDensityBatch!(
    rhos,
    temps,
    winds,
    model::DensityScaleModel,
    hs,
    lats,
    lons,
    el_time,
    wind::Bool,
    p,
)
    Base.invokelatest(getDensityBatch!, rhos, temps, winds, model.base, hs, lats, lons, el_time, wind, p)
    @inbounds for i in eachindex(rhos)
        rhos[i] *= model.scale
    end
    return nothing
end

@inline _is_aero_case(case::Symbol)::Bool = case in (:gram_aero, :full_environment)

function _validate_aero_solver(solver::Symbol)::Symbol
    solver in (:auto_stiff, :rodas5p, :auto_stiff_then_rodas5p) && return solver
    throw(ArgumentError("Unsupported aero solver '$(solver)'. Use one of: auto_stiff, rodas5p, auto_stiff_then_rodas5p."))
end

function _aero_solver_for_config(sample::NamedTuple)::Symbol
    !_is_aero_case(sample.dynamics_case) && return :default
    requested = _validate_aero_solver(Symbol(_sample_field(sample, :aero_solver, DEFAULT_AERO_SOLVER)))
    requested === :auto_stiff_then_rodas5p && return :auto_stiff
    return requested
end

function _aero_solver_maxiters(sample::NamedTuple)::Int
    solver = _aero_solver_for_config(sample)
    solver === :rodas5p && return Int(_sample_field(sample, :aero_stiff_maxiters, DEFAULT_AERO_STIFF_MAXITERS))
    return Int(_sample_field(sample, :aero_auto_maxiters, DEFAULT_AERO_AUTO_MAXITERS))
end

function _aero_auto_stiff_switch_max(sample::NamedTuple)::Int
    return Int(_sample_field(sample, :aero_auto_stiff_switch_max, DEFAULT_AERO_AUTO_STIFF_SWITCH_MAX))
end

function _aero_tolerance_scale(sample::NamedTuple)::Float64
    solver = _aero_solver_for_config(sample)
    solver === :rodas5p || return 1.0
    scale = Float64(_sample_field(sample, :aero_stiff_tol_scale, DEFAULT_AERO_STIFF_TOL_SCALE))
    scale > 0.0 || throw(ArgumentError("aero_stiff_tol_scale must be > 0; got $(scale)."))
    return scale
end

function _planet(planet_id::Symbol)
    if planet_id == :mars
        return Mars("", SPICE_PATH)
    elseif planet_id == :earth
        return Earth("", SPICE_PATH)
    elseif planet_id == :venus
        return Venus("", SPICE_PATH)
    elseif planet_id == :titan
        return Titan("", SPICE_PATH)
    end
    throw(ArgumentError("Unsupported planet '$planet_id'. Supported planets: mars, earth, venus, titan."))
end

function _planet_name(planet_id::Symbol)::String
    planet_id == :mars && return "mars"
    planet_id == :earth && return "earth"
    planet_id == :venus && return "venus"
    planet_id == :titan && return "titan"
    return String(planet_id)
end

function _periapsis_altitude_m(planet_id::Symbol, regime::Symbol)::Float64
    table = if planet_id == :mars
        Dict(:shallow => 150e3, :nominal => 125e3, :deep => 110e3)
    elseif planet_id == :venus
        Dict(:shallow => 180e3, :nominal => 150e3, :deep => 135e3)
    elseif planet_id == :titan
        Dict(:shallow => 900e3, :nominal => 720e3, :deep => 650e3)
    else
        Dict(:shallow => 180e3, :nominal => 145e3, :deep => 120e3)
    end
    return table[regime]
end


function _harmonics_file(planet_id::Symbol)::String
    if planet_id == :mars
        return joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "Mars50c.csv")
    elseif planet_id == :earth
        return joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    elseif planet_id == :venus
        return joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "MGNP180U.csv")
    elseif planet_id == :titan
        return joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "titan5.csv")
    end
    throw(ArgumentError("No harmonics file configured for '$planet_id'."))
end

function _harmonics_degree_order(planet_id::Symbol)::Tuple{Int, Int}
    planet_id == :titan && return (5, 5)
    return (20, 20)
end

function _third_body_names(planet_id::Symbol)
    planet_id == :earth && return ("Sun", "Moon")
    planet_id == :titan && return ("Saturn", "Sun")
    return ("Sun",)
end

function _density_model(planet_id::Symbol, dynamics_case::Symbol, density_scale::Float64)
    if !_is_aero_case(dynamics_case)
        return NoAtmosphereModel()
    end
    setup_gram_example!(@__MODULE__)
    base = Base.invokelatest(GRAMAtmosphereModel; planet_name=_planet_name(planet_id))
    return DensityScaleModel(base, density_scale)
end

function _spacecraft(planet, rp_alt_m::Float64, ra_alt_m::Float64;
        mass_scale::Float64=NOMINAL_SPACECRAFT_MASS_SCALE,
        inclination_deg::Float64=NOMINAL_INCLINATION_DEG,
        argp_deg::Float64=NOMINAL_ARGP_DEG)
    return make_three_body_spacecraft(
        bus_dims=(2.2, 2.6, 1.7),
        panel_dims=(0.01, 3.0, 1.5),
        bus_mass=350.0 * mass_scale,
        panel_mass_each=10.0 * mass_scale,
        panel_offset_y=2.6,
        ic=InitialCondition(
            ra=planet.Rp_e + ra_alt_m,
            rp=planet.Rp_e + rp_alt_m,
            i=inclination_deg,
            ω=argp_deg,
            Ω=30.0,
            ν=180.0,
        ),
        reflection_coefficient=0.9,
        prop_mass=50.0 * mass_scale,
        id=1,
    )
end

function _dynamic_effectors(planet_id::Symbol, planet, spacecraft, dynamics_case::Symbol)
    if dynamics_case == :two_body
        return (InverseSquaredGravityModel(),)
    elseif dynamics_case == :j2
        return (InverseSquaredJ2GravityModel(),)
    elseif dynamics_case == :harmonics_low
        degree, order = _harmonics_degree_order(planet_id)
        return (GravitationalHarmonicsModel(degree, order, _harmonics_file(planet_id), planet),)
    elseif dynamics_case == :srp
        return (
            InverseSquaredGravityModel(),
            SolarRadiationPressureModel(spacecraft.root.reflection_coefficient, spacecraft.root.ref_area),
        )
    elseif dynamics_case == :third_body_sun
        body_names = _third_body_names(planet_id)
        return (
            InverseSquaredGravityModel(),
            NBodyGravityModel(body_names=body_names, primary_body_name=planet.name, planet=planet),
        )
    elseif dynamics_case == :gram_aero
        return (InverseSquaredGravityModel(), AerodynamicCoefficientfM())
    elseif dynamics_case == :full_environment
        body_names = _third_body_names(planet_id)
        degree, order = _harmonics_degree_order(planet_id)
        return (
            NBodyGravityModel(body_names=body_names, primary_body_name=planet.name, planet=planet),
            GravitationalHarmonicsModel(degree, order, _harmonics_file(planet_id), planet),
            SolarRadiationPressureModel(spacecraft.root.reflection_coefficient, spacecraft.root.ref_area),
            AerodynamicCoefficientfM(),
        )
    end
    throw(ArgumentError("Unsupported dynamics case '$dynamics_case'."))
end

@inline function _period_s(planet, rp_alt_m::Float64, ra_alt_m::Float64)::Float64
    a = planet.Rp_e + 0.5 * (rp_alt_m + ra_alt_m)
    return 2π * sqrt(a^3 / planet.μ)
end

@inline function _orbit_index_for_time(t::Float64, period_s::Float64, norbits::Int)::Int
    t <= 0.0 && return 1
    return clamp(ceil(Int, t / period_s), 1, max(1, norbits))
end

function _eccentric_anomaly_from_true_anomaly(ν::Float64, e::Float64)::Float64
    return 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * sin(ν / 2.0), cos(ν / 2.0))
end

function _atmosphere_half_window_s(planet, rp_alt_m::Float64, ra_alt_m::Float64, ei_alt_m::Float64)::Float64
    rp_r = planet.Rp_e + rp_alt_m
    ra_r = planet.Rp_e + ra_alt_m
    interface_r = planet.Rp_e + ei_alt_m
    period_s = _period_s(planet, rp_alt_m, ra_alt_m)
    interface_r <= rp_r && return 0.0
    interface_r >= ra_r && return 0.5 * period_s

    a = 0.5 * (rp_r + ra_r)
    e = (ra_r - rp_r) / (ra_r + rp_r)
    e <= eps(Float64) && return 0.0
    p = a * (1.0 - e^2)
    cosν = clamp((p / interface_r - 1.0) / e, -1.0, 1.0)
    ν = acos(cosν)
    E = _eccentric_anomaly_from_true_anomaly(ν, e)
    M = E - e * sin(E)
    n = 2π / period_s
    return abs(M / n)
end

function _aero_save_times(
    planet,
    rp_alt_m::Float64,
    ra_alt_m::Float64,
    ei_alt_m::Float64,
    norbits::Int,
)::Vector{Float64}
    period_s = _period_s(planet, rp_alt_m, ra_alt_m)
    mission_time = norbits * period_s
    half_window = _atmosphere_half_window_s(planet, rp_alt_m, ra_alt_m, ei_alt_m)
    times = Float64[0.0, mission_time]
    for orbit in 1:(norbits - 1)
        push!(times, orbit * period_s)
    end

    t = OUTSIDE_ATMOSPHERE_SAVE_DT_S
    while t < mission_time
        push!(times, t)
        t += OUTSIDE_ATMOSPHERE_SAVE_DT_S
    end

    if half_window > 0.0
        for orbit in 0:(norbits - 1)
            peri_t = (orbit + 0.5) * period_s
            t_inside = max(orbit * period_s, peri_t - half_window)
            t_exit = min((orbit + 1) * period_s, peri_t + half_window)
            while t_inside <= t_exit
                push!(times, t_inside)
                t_inside += INSIDE_ATMOSPHERE_SAVE_DT_S
            end
            push!(times, t_exit)
        end
    end

    sort!(times)
    unique_times = Float64[]
    last_t = -Inf
    for t in times
        t_clamped = clamp(t, 0.0, mission_time)
        if isempty(unique_times) || abs(t_clamped - last_t) > 1e-7
            push!(unique_times, t_clamped)
            last_t = t_clamped
        end
    end
    return unique_times
end

function _study_tolerances(
    base::SM.IntegrationTolerances,
    dt_max_orbit_s::Float64,
    dt_max_atmosphere_s::Float64,
    tolerance_scale::Float64=1.0,
)::SM.IntegrationTolerances
    tolerance_scale > 0.0 || throw(ArgumentError("tolerance_scale must be > 0; got $(tolerance_scale)."))
    scaled(value) = value == 0.0 ? 0.0 : value * tolerance_scale
    return SM.IntegrationTolerances(
        reltol=scaled(base.reltol),
        abstol=scaled(base.abstol),
        reltol_orbit=scaled(base.reltol_orbit),
        abstol_orbit=scaled(base.abstol_orbit),
        reltol_atmosphere=scaled(base.reltol_atmosphere),
        abstol_atmosphere=scaled(base.abstol_atmosphere),
        reltol_quaternion=scaled(base.reltol_quaternion),
        abstol_quaternion=scaled(base.abstol_quaternion),
        reltol_mass=scaled(base.reltol_mass),
        abstol_mass=scaled(base.abstol_mass),
        reltol_heat_load=scaled(base.reltol_heat_load),
        abstol_heat_load=scaled(base.abstol_heat_load),
        reltol_angular_rate=scaled(base.reltol_angular_rate),
        abstol_angular_rate=scaled(base.abstol_angular_rate),
        dt_max=base.dt_max,
        dt_max_orbit=dt_max_orbit_s,
        dt_max_atmosphere=dt_max_atmosphere_s,
    )
end

function _with_study_dtmax(
    args::SM.SimulationConfiguration,
    dt_max_orbit_s::Float64,
    dt_max_atmosphere_s::Float64,
    tolerance_scale::Float64=1.0,
)::SM.SimulationConfiguration
    dt_max_orbit_s > 0.0 || throw(ArgumentError("dt_max_orbit_s must be > 0; got $(dt_max_orbit_s)."))
    dt_max_atmosphere_s > 0.0 || throw(ArgumentError("dt_max_atmosphere_s must be > 0; got $(dt_max_atmosphere_s)."))
    return SM.SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=args.simulation_settings,
        mission_configuration=args.mission_configuration,
        environment_model=args.environment_model,
        dynamics_model=args.dynamics_model,
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=_study_tolerances(args.integration_tolerances, dt_max_orbit_s, dt_max_atmosphere_s, tolerance_scale),
        solver_config=args.solver_config,
    )
end

function _make_config(sample::NamedTuple, results_directory::String, results::Bool)
    planet = _planet(sample.planet)
    rp_alt_m = _periapsis_altitude_m(sample.planet, sample.periapsis_regime)
    ra_alt_m = sample.apoapsis_alt_km * 1e3
    spacecraft = _spacecraft(planet, rp_alt_m, ra_alt_m;
        mass_scale=sample.spacecraft_mass_scale,
        inclination_deg=sample.inclination_deg,
        argp_deg=sample.argp_deg)
    mission_time = sample.norbits * _period_s(planet, rp_alt_m, ra_alt_m)
    dynamic_effectors = _dynamic_effectors(sample.planet, planet, spacecraft, sample.dynamics_case)
    density_model = _density_model(sample.planet, sample.dynamics_case, sample.density_scale)
    aero_solver = _aero_solver_for_config(sample)
    solver_config = _is_aero_case(sample.dynamics_case) ? SM.SolverConfig(
        solver_mode=aero_solver,
        maxiters=_aero_solver_maxiters(sample),
        auto_stiff_switch_max=_aero_auto_stiff_switch_max(sample),
    ) : nothing

    args = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=mission_time,
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=dynamic_effectors,
        density_model=density_model,
        orientation_sim=false,
        keplerian=false,
        EI_km=max(220.0, rp_alt_m / 1e3 + 80.0),
        verbose=false,
        results=results,
        results_directory=results_directory,
        solver_config=solver_config,
    )
    dt_max_orbit_s = Float64(_sample_field(sample, :dt_max_orbit_s, DEFAULT_DT_MAX_ORBIT_S))
    dt_max_atmosphere_s = Float64(_sample_field(sample, :dt_max_atmosphere_s, DEFAULT_DT_MAX_ATMOSPHERE_S))
    return _with_study_dtmax(args, dt_max_orbit_s, dt_max_atmosphere_s, _aero_tolerance_scale(sample))
end

function _ofat_spacecraft_combos(spec::StudySpec)::Vector{NTuple{3, Float64}}
    seen = Set{NTuple{3, Float64}}()
    combos = NTuple{3, Float64}[]
    function _add!(m, i, a)
        key = (Float64(m), Float64(i), Float64(a))
        key in seen && return
        push!(seen, key)
        push!(combos, key)
    end
    nm, ni, na = spec.nominal_mass_scale, spec.nominal_inclination_deg, spec.nominal_argp_deg
    _add!(nm, ni, na)
    for m in spec.spacecraft_mass_scales; _add!(m, ni, na); end
    for i in spec.inclinations_deg;      _add!(nm, i, na); end
    for a in spec.argp_degs;             _add!(nm, ni, a); end
    return combos
end

function make_samples(spec::StudySpec)
    samples = NamedTuple[]
    periapsis_regimes = spec.smoke ? spec.periapsis_regimes[1:min(end, 1)] : spec.periapsis_regimes
    dynamics_cases = spec.smoke ? [:two_body, :gram_aero] : spec.dynamics_cases
    density_items = sort(collect(spec.density_scales); by=x -> String(x.first))
    nominal_combo = (spec.nominal_mass_scale, spec.nominal_inclination_deg, spec.nominal_argp_deg)
    spacecraft_combos = spec.smoke ? [nominal_combo] : _ofat_spacecraft_combos(spec)
    for planet in spec.planets
        all_apo = get(spec.apoapsis_altitudes_km, planet, DEFAULT_APOAPSIS_ALTITUDES_KM[planet])
        apo_alts = spec.smoke ? all_apo[1:min(end, 1)] : all_apo
        for periapsis_regime in periapsis_regimes
            for apoapsis_alt_km in apo_alts
                for (mass_scale, inclination_deg, argp_deg) in spacecraft_combos
                    for dynamics_case in dynamics_cases
                        density_iter = _is_aero_case(dynamics_case) ? density_items : [(:none => NaN)]
                        for density_pair in density_iter
                            push!(samples, (
                                planet=planet,
                                periapsis_regime=periapsis_regime,
                                apoapsis_alt_km=apoapsis_alt_km,
                                spacecraft_mass_scale=mass_scale,
                                inclination_deg=inclination_deg,
                                argp_deg=argp_deg,
                                dynamics_case=dynamics_case,
                                density_case=Symbol(density_pair.first),
                                density_scale=Float64(density_pair.second),
                                norbits=spec.smoke ? 1 : spec.norbits,
                                aero_solver=spec.aero_solver,
                                aero_auto_maxiters=spec.aero_auto_maxiters,
                                aero_stiff_maxiters=spec.aero_stiff_maxiters,
                                aero_auto_stiff_switch_max=spec.aero_auto_stiff_switch_max,
                                aero_stiff_tol_scale=spec.aero_stiff_tol_scale,
                                case_timeout_min=spec.case_timeout_min,
                                deorbit_bailout=spec.deorbit_bailout,
                                deorbit_bailout_margin_km=spec.deorbit_bailout_margin_km,
                                deorbit_bailout_check_dt_s=spec.deorbit_bailout_check_dt_s,
                                dt_max_orbit_s=spec.dt_max_orbit_s,
                                dt_max_atmosphere_s=spec.dt_max_atmosphere_s,
                            ))
                        end
                    end
                end
            end
        end
    end
    return samples
end

function _state_pos_vel(u)
    sc = getproperty(u, :sc)[1]
    pos = SVector{3, Float64}(getproperty(sc, :pos))
    vel = SVector{3, Float64}(getproperty(sc, :vel))
    return pos, vel
end

function _elements_from_rv(pos::SVector{3, Float64}, vel::SVector{3, Float64}, μ::Float64)
    r = norm(pos)
    v2 = dot(vel, vel)
    hvec = cross(pos, vel)
    h = norm(hvec)
    nvec = cross(SVector{3, Float64}(0.0, 0.0, 1.0), hvec)
    n = norm(nvec)
    evec = cross(vel, hvec) / μ - pos / r
    e = norm(evec)
    energy = 0.5 * v2 - μ / r
    a = -μ / (2.0 * energy)
    inc = acos(clamp(hvec[3] / h, -1.0, 1.0))
    raan = n > 0.0 ? atan(nvec[2], nvec[1]) : 0.0
    argp = (n > 0.0 && e > 0.0) ? acos(clamp(dot(nvec, evec) / (n * e), -1.0, 1.0)) : 0.0
    if e > 0.0 && evec[3] < 0.0
        argp = 2π - argp
    end
    return (; a, e, inc, raan, argp, energy, rp=a * (1.0 - e), ra=a * (1.0 + e))
end

function _density_metrics(args, sol)
    density_model = args.environment_model.density_model
    if density_model isa NoAtmosphereModel
        return (peak_density=0.0, integrated_density=0.0, max_dynamic_pressure=0.0, integrated_dynamic_pressure=0.0, time_below_interface_s=0.0)
    end

    planet = args.environment_model.planet
    EI_m = args.environment_model.EI * 1e3
    peak_density = 0.0
    integrated_density = 0.0
    max_dynamic_pressure = 0.0
    integrated_dynamic_pressure = 0.0
    time_below_interface_s = 0.0

    for i in eachindex(sol.t)
        pos, vel = _state_pos_vel(sol.u[i])
        r = norm(pos)
        alt = r - planet.Rp_e
        lat = asind(clamp(pos[3] / r, -1.0, 1.0))
        lon = atan(pos[2], pos[1]) * 180.0 / π
        density_state = try
            getDensity(density_model, Float64(alt), Float64(lat), Float64(lon), Float64(sol.t[i]), true)
        catch
            return (peak_density=NaN, integrated_density=NaN, max_dynamic_pressure=NaN, integrated_dynamic_pressure=NaN, time_below_interface_s=time_below_interface_s)
        end
        rho, _, wind = density_state
        q = 0.5 * rho * norm(vel - wind)^2
        peak_density = max(peak_density, rho)
        max_dynamic_pressure = max(max_dynamic_pressure, q)
        if i > 1
            dt = sol.t[i] - sol.t[i - 1]
            integrated_density += rho * dt
            integrated_dynamic_pressure += q * dt
            if alt <= EI_m
                time_below_interface_s += dt
            end
        end
    end

    return (; peak_density, integrated_density, max_dynamic_pressure, integrated_dynamic_pressure, time_below_interface_s)
end

@inline function _state_mass_kg(u)::Float64
    sc = getproperty(u, :sc)[1]
    return Float64(getproperty(sc, :mass))
end

function _trajectory_dataframe_from_solution(sol)
    times = Float64[]
    pos1 = Float64[]
    pos2 = Float64[]
    pos3 = Float64[]
    vel1 = Float64[]
    vel2 = Float64[]
    vel3 = Float64[]
    mass = Float64[]
    for i in eachindex(sol.t)
        pos, vel = _state_pos_vel(sol.u[i])
        push!(times, Float64(sol.t[i]))
        push!(pos1, pos[1]); push!(pos2, pos[2]); push!(pos3, pos[3])
        push!(vel1, vel[1]); push!(vel2, vel[2]); push!(vel3, vel[3])
        push!(mass, _state_mass_kg(sol.u[i]))
    end
    return DataFrame(
        time=times,
        sc1_pos_1=pos1,
        sc1_pos_2=pos2,
        sc1_pos_3=pos3,
        sc1_vel_1=vel1,
        sc1_vel_2=vel2,
        sc1_vel_3=vel3,
        sc1_mass=mass,
    )
end

@inline function _row_float(row, name::Symbol, default::Float64=NaN)::Float64
    hasproperty(row, name) || return default
    value = getproperty(row, name)
    value === missing && return default
    return Float64(value)
end

function _trajectory_dataframe(args, sol, trajectory_csv_path::String)
    if isfile(trajectory_csv_path)
        try
            df = CSV.read(trajectory_csv_path, DataFrame)
            required = (:time, :sc1_pos_1, :sc1_pos_2, :sc1_pos_3, :sc1_vel_1, :sc1_vel_2, :sc1_vel_3, :sc1_mass)
            if all(name -> name in Symbol.(names(df)), required)
                return df
            end
            @warn "Ignoring stale simulation_results.csv without trajectory columns; reconstructing from solution" path=trajectory_csv_path
        catch err
            @warn "Ignoring unreadable simulation_results.csv; reconstructing from solution" path=trajectory_csv_path exception=(err, catch_backtrace())
        end
    end
    return _trajectory_dataframe_from_solution(sol)
end

function _state_sample_from_row(row, spacecraft)
    pos = SVector{3, Float64}(
        _row_float(row, :sc1_pos_1),
        _row_float(row, :sc1_pos_2),
        _row_float(row, :sc1_pos_3),
    )
    vel = SVector{3, Float64}(
        _row_float(row, :sc1_vel_1),
        _row_float(row, :sc1_vel_2),
        _row_float(row, :sc1_vel_3),
    )
    fallback_mass = Float64(spacecraft.dry_mass + spacecraft.prop_mass)
    mass = _row_float(row, :sc1_mass, fallback_mass)
    return SM.StateSample(pos, vel, mass; spacecraft=spacecraft)
end

function _planet_frame_sample(args, state::SM.StateSample, t::Float64)
    planet = args.environment_model.planet
    ephemerides_model = args.environment_model.ephemerides_model
    et0 = SM.ephemerides_time_seconds(args.initial_time, ephemerides_model)
    et = et0 + t
    l_pi = SM.planet_frame_lpi(planet, et, ephemerides_model)
    pos_pp = SVector{3, Float64}(l_pi * state.pos_ii)
    vel_pp = SVector{3, Float64}(l_pi * (state.vel_ii - cross(planet.ω, state.pos_ii)))
    alt, lat, lon = SM.SimulationCallbacks.rtolatlong(pos_pp, planet)
    return SM.PlanetFrameSample(l_pi, pos_pp, vel_pp, alt, lat, lon)
end

function _solar_sample(args, t::Float64)
    planet = args.environment_model.planet
    ephemerides_model = args.environment_model.ephemerides_model
    et0 = SM.ephemerides_time_seconds(args.initial_time, ephemerides_model)
    et = et0 + t
    primary = SM.DynamicEffectors._spice_query_name(planet.name)
    sun_pos = SM.EphemeridesModels.spice_position_j2000_m("sun", et, primary)
    return SM.SolarEphemerisSample(SVector{3, Float64}(sun_pos))
end

function _third_body_sample(args, model::SM.NBodyGravityModel, t::Float64)
    ephemerides_model = args.environment_model.ephemerides_model
    et0 = SM.ephemerides_time_seconds(args.initial_time, ephemerides_model)
    et = et0 + t
    primary = SM.DynamicEffectors._spice_query_name(model.primary_body_name)
    positions = ntuple(length(model.body_names)) do k
        body = SM.DynamicEffectors._spice_query_name(model.body_names[k])
        SVector{3, Float64}(SM.EphemeridesModels.spice_position_j2000_m(body, et, primary))
    end
    return SM.ThirdBodyEphemerisSample(model.body_names, positions)
end

function _aero_force_from_row(row)::SVector{3, Float64}
    force = SVector{3, Float64}(0.0, 0.0, 0.0)
    for prefix in (:sc1_drag, :sc1_lift, :sc1_cross)
        force += SVector{3, Float64}(
            _row_float(row, Symbol(prefix, "_1"), 0.0),
            _row_float(row, Symbol(prefix, "_2"), 0.0),
            _row_float(row, Symbol(prefix, "_3"), 0.0),
        )
    end
    return force
end

function _central_force(args, state::SM.StateSample, env)::SVector{3, Float64}
    force, _ = SM.wrench(InverseSquaredGravityModel(), state, env, 0.0)
    return force
end

function _c20_harmonics_force(model::GravitationalHarmonicsModel, state::SM.StateSample, planet_frame)::SVector{3, Float64}
    model.L >= 2 || return SVector{3, Float64}(0.0, 0.0, 0.0)
    C20 = Float64(model.C[3, 1])
    isfinite(C20) && C20 != 0.0 || return SVector{3, Float64}(0.0, 0.0, 0.0)
    pos_pp = planet_frame.pos_pp
    r = norm(pos_pp)
    r > 0.0 || return SVector{3, Float64}(0.0, 0.0, 0.0)
    x, y, z = pos_pp
    r2 = r * r
    z2 = z * z
    J2 = -sqrt(5.0) * C20
    common = 1.5 * model.planet.μ * J2 * model.reference_radius_m^2 / (r2 * r2)
    accel_pp = SVector{3, Float64}(
        common * (x / r) * (5.0 * z2 / r2 - 1.0),
        common * (y / r) * (5.0 * z2 / r2 - 1.0),
        common * (z / r) * (5.0 * z2 / r2 - 3.0),
    )
    return state.mass_kg * (planet_frame.l_pi' * accel_pp)
end

function _force_for_effector(args, state::SM.StateSample, row, effector, t::Float64)::SVector{3, Float64}
    planet_frame = nothing
    solar = nothing
    third_bodies = nothing
    if effector isa InverseSquaredJ2GravityModel || effector isa GravitationalHarmonicsModel
        planet_frame = _planet_frame_sample(args, state, t)
    elseif effector isa SolarRadiationPressureModel
        solar = _solar_sample(args, t)
    elseif effector isa NBodyGravityModel
        third_bodies = _third_body_sample(args, effector, t)
    elseif effector isa AerodynamicCoefficientfM
        return _aero_force_from_row(row)
    end
    env = SM.EnvironmentSample(args.environment_model.planet; planet_frame=planet_frame, solar=solar, third_bodies=third_bodies)
    force, _ = SM.wrench(effector, state, env, t)
    if effector isa InverseSquaredJ2GravityModel
        return force - _central_force(args, state, env)
    elseif effector isa GravitationalHarmonicsModel
        return force - _central_force(args, state, env) - _c20_harmonics_force(effector, state, planet_frame)
    end
    return force
end

function _active_perturbation_force(args, sample::NamedTuple, state::SM.StateSample, row, t::Float64)::SVector{3, Float64}
    sample.dynamics_case == :two_body && return SVector{3, Float64}(0.0, 0.0, 0.0)

    total = SVector{3, Float64}(0.0, 0.0, 0.0)
    for effector in args.dynamics_model.dynamic_effectors
        if effector isa InverseSquaredGravityModel
            continue
        end
        total += _force_for_effector(args, state, row, effector, t)
    end
    return total
end

mutable struct OrbitChunkWriter
    args::Any
    sample::NamedTuple
    case_dir::String
    period_s::Float64
    norbits::Int
    current_orbit::Int
    buffer::Vector{NamedTuple}
    chunks::Vector{NamedTuple}
    samples_saved::Int
    peak_density::Float64
    integrated_density::Float64
    max_dynamic_pressure::Float64
    integrated_dynamic_pressure::Float64
    time_below_interface_s::Float64
    previous_t::Float64
    previous_rho::Float64
    previous_q::Float64
    previous_inside::Bool
    first_pos::Any
    first_vel::Any
    last_pos::Any
    last_vel::Any
end

function OrbitChunkWriter(args, sample::NamedTuple, case_dir::String, period_s::Float64)
    return OrbitChunkWriter(
        args,
        sample,
        case_dir,
        period_s,
        Int(sample.norbits),
        1,
        NamedTuple[],
        NamedTuple[],
        0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        NaN,
        0.0,
        0.0,
        false,
        nothing,
        nothing,
        nothing,
        nothing,
    )
end

@inline function _cache_force(cache::Vector{SVector{3, Float64}}, idx::Int=1)::SVector{3, Float64}
    return idx <= length(cache) ? cache[idx] : SVector{3, Float64}(0.0, 0.0, 0.0)
end

function _stream_density_update!(
    writer::OrbitChunkWriter,
    state::SM.StateSample,
    p,
    t::Float64,
)
    args = writer.args
    EI_m = args.environment_model.EI * 1e3
    frame = _planet_frame_sample(args, state, t)
    alt = frame.alt_m
    inside = alt <= EI_m
    rho = 0.0
    q = 0.0
    density_model = SM.SimulationCallbacks._density_model_for_sat(p, 1)
    if !(density_model isa NoAtmosphereModel)
        density_state = try
            getDensity(
                density_model,
                Float64(frame.alt_m),
                Float64(frame.lat_rad),
                Float64(frame.lon_rad),
                Float64(t),
                true,
                p,
            )
        catch
            writer.peak_density = NaN
            writer.integrated_density = NaN
            writer.max_dynamic_pressure = NaN
            writer.integrated_dynamic_pressure = NaN
            writer.previous_t = t
            writer.previous_inside = inside
            return (; alt, rho=NaN, q=NaN, inside)
        end
        rho, _, wind = density_state
        q = 0.5 * rho * norm(frame.vel_pp - wind)^2
        writer.peak_density = max(writer.peak_density, rho)
        writer.max_dynamic_pressure = max(writer.max_dynamic_pressure, q)
    end
    if isfinite(writer.previous_t)
        dt = max(0.0, t - writer.previous_t)
        if isfinite(writer.integrated_density)
            writer.integrated_density += 0.5 * (writer.previous_rho + rho) * dt
            writer.integrated_dynamic_pressure += 0.5 * (writer.previous_q + q) * dt
        end
        if inside || writer.previous_inside
            writer.time_below_interface_s += dt
        end
    end
    writer.previous_t = t
    writer.previous_rho = rho
    writer.previous_q = q
    writer.previous_inside = inside
    return (; alt, rho, q, inside)
end

function _flush_orbit_chunk!(writer::OrbitChunkWriter)
    isempty(writer.buffer) && return nothing
    orbit = writer.current_orbit
    chunk_path = joinpath(writer.case_dir, @sprintf("trajectory_with_active_force_orbit_%03d.feather", orbit))
    _write_feather(chunk_path, DataFrame(writer.buffer))
    push!(writer.chunks, (
        orbit=orbit,
        path=chunk_path,
        rows=length(writer.buffer),
        t_start=Float64(first(writer.buffer).time),
        t_end=Float64(last(writer.buffer).time),
    ))
    empty!(writer.buffer)
    return nothing
end

function _stream_sample!(writer::OrbitChunkWriter, u, t_raw, integrator)
    t = Float64(t_raw)
    orbit = _orbit_index_for_time(t, writer.period_s, writer.norbits)
    if orbit != writer.current_orbit
        _flush_orbit_chunk!(writer)
        writer.current_orbit = orbit
    end

    args = writer.args
    p = integrator.p
    spacecraft = args.dynamics_model.spacecraft[1]
    pos = SimulationEngine._state_position_ii(u, 1)
    vel = SimulationEngine._state_velocity_ii(u, 1)
    mass = SimulationEngine._state_mass_kg(u, args, 1)
    if writer.samples_saved == 0
        writer.first_pos = pos
        writer.first_vel = vel
    end
    writer.last_pos = pos
    writer.last_vel = vel
    state = SM.StateSample(pos, vel, mass; spacecraft=spacecraft)
    drag = _cache_force(p.save_cache.drag_cache)
    lift = _cache_force(p.save_cache.lift_cache)
    cross_force = _cache_force(p.save_cache.cross_cache)
    row_for_force = (
        time=t,
        sc1_pos_1=pos[1], sc1_pos_2=pos[2], sc1_pos_3=pos[3],
        sc1_vel_1=vel[1], sc1_vel_2=vel[2], sc1_vel_3=vel[3],
        sc1_mass=mass,
        sc1_drag_1=drag[1], sc1_drag_2=drag[2], sc1_drag_3=drag[3],
        sc1_lift_1=lift[1], sc1_lift_2=lift[2], sc1_lift_3=lift[3],
        sc1_cross_1=cross_force[1], sc1_cross_2=cross_force[2], sc1_cross_3=cross_force[3],
    )
    force = _active_perturbation_force(args, writer.sample, state, row_for_force, t)
    density = _stream_density_update!(writer, state, p, t)
    push!(writer.buffer, merge(row_for_force, (
        orbit=orbit,
        altitude_m=density.alt,
        density_kg_m3=density.rho,
        dynamic_pressure_pa=density.q,
        in_atmosphere=density.inside,
        active_perturbation_force_ii_1=force[1],
        active_perturbation_force_ii_2=force[2],
        active_perturbation_force_ii_3=force[3],
        active_perturbation_force_mag=norm(force),
    )))
    writer.samples_saved += 1
    return 0
end

function _finish_streaming_writer!(writer::OrbitChunkWriter)
    _flush_orbit_chunk!(writer)
    manifest_path = joinpath(writer.case_dir, "orbit_chunk_manifest.feather")
    _write_feather(manifest_path, DataFrame(writer.chunks))
    combined_paths = _combine_orbit_chunks!(writer, manifest_path)
    return (
        simulation_results_csv="",
        simulation_results_feather=combined_paths.simulation_results_feather,
        active_force_history_csv="",
        active_force_history_feather=combined_paths.active_force_history_feather,
        trajectory_with_active_force_csv="",
        trajectory_with_active_force_feather=combined_paths.trajectory_with_active_force_feather,
        trajectory_chunk_manifest=manifest_path,
        orbit_chunk_count=length(writer.chunks),
        samples_saved=writer.samples_saved,
        peak_density=writer.peak_density,
        integrated_density=writer.integrated_density,
        max_dynamic_pressure=writer.max_dynamic_pressure,
        integrated_dynamic_pressure=writer.integrated_dynamic_pressure,
        time_below_interface_s=writer.time_below_interface_s,
    )
end

function _combine_orbit_chunks!(writer::OrbitChunkWriter, manifest_path::String)
    chunks = DataFrame(Arrow.Table(manifest_path))
    if nrow(chunks) == 0
        empty_df = DataFrame()
        simulation_results_path = joinpath(writer.case_dir, "simulation_results.feather")
        active_force_path = joinpath(writer.case_dir, "active_force_history.feather")
        combined_path = joinpath(writer.case_dir, "trajectory_with_active_force.feather")
        _write_feather(simulation_results_path, empty_df)
        _write_feather(active_force_path, empty_df)
        _write_feather(combined_path, empty_df)
        return (
            simulation_results_feather=simulation_results_path,
            active_force_history_feather=active_force_path,
            trajectory_with_active_force_feather=combined_path,
        )
    end

    sort!(chunks, :orbit)
    chunk_dfs = DataFrame[]
    for path in String.(chunks.path)
        push!(chunk_dfs, DataFrame(Arrow.Table(path)))
    end
    combined_df = vcat(chunk_dfs...; cols=:union)
    sort!(combined_df, :time)

    force_cols = [
        :time,
        :active_perturbation_force_ii_1,
        :active_perturbation_force_ii_2,
        :active_perturbation_force_ii_3,
        :active_perturbation_force_mag,
    ]
    active_force_df = combined_df[:, [c for c in force_cols if c in propertynames(combined_df)]]

    active_force_names = Set(force_cols[2:end])
    trajectory_names = [name for name in propertynames(combined_df) if !(name in active_force_names)]
    simulation_df = combined_df[:, trajectory_names]

    simulation_results_path = joinpath(writer.case_dir, "simulation_results.feather")
    active_force_path = joinpath(writer.case_dir, "active_force_history.feather")
    combined_path = joinpath(writer.case_dir, "trajectory_with_active_force.feather")
    _write_feather(simulation_results_path, simulation_df)
    _write_feather(active_force_path, active_force_df)
    _write_feather(combined_path, combined_df)
    empty!(chunk_dfs)
    GC.gc(false)
    return (
        simulation_results_feather=simulation_results_path,
        active_force_history_feather=active_force_path,
        trajectory_with_active_force_feather=combined_path,
    )
end

function _streaming_save_callback(writer::OrbitChunkWriter, save_times::Vector{Float64})
    return PresetTimeCallback(save_times, integrator -> _stream_sample!(writer, integrator.u, integrator.t, integrator))
end

mutable struct StudyGuardState
    start_wall_s::Float64
    timeout_s::Float64
    deorbit_enabled::Bool
    deorbit_ra_alt_limit_m::Float64
    deorbit_check_dt_s::Float64
    next_deorbit_check_t::Float64
    termination_reason::String
    termination_time_s::Float64
    termination_ra_alt_m::Float64
    termination_alt_m::Float64
end

function StudyGuardState(sample::NamedTuple, args)::StudyGuardState
    timeout_min = Float64(_sample_field(sample, :case_timeout_min, DEFAULT_CASE_TIMEOUT_MIN))
    timeout_s = timeout_min > 0.0 ? timeout_min * 60.0 : Inf
    deorbit_enabled = Bool(_sample_field(sample, :deorbit_bailout, DEFAULT_DEORBIT_BAILOUT)) && _is_aero_case(sample.dynamics_case)
    margin_m = Float64(_sample_field(sample, :deorbit_bailout_margin_km, DEFAULT_DEORBIT_BAILOUT_MARGIN_KM)) * 1e3
    check_dt_s = max(1.0, Float64(_sample_field(sample, :deorbit_bailout_check_dt_s, DEFAULT_DEORBIT_BAILOUT_CHECK_DT_S)))
    return StudyGuardState(
        time(),
        timeout_s,
        deorbit_enabled,
        args.environment_model.EI * 1e3 + margin_m,
        check_dt_s,
        check_dt_s,
        "",
        NaN,
        NaN,
        NaN,
    )
end

function _study_guard_callback(guard::StudyGuardState, writer::OrbitChunkWriter)
    condition(u, t, integrator) = true

    function affect!(integrator)
        t = Float64(integrator.t)
        if isempty(guard.termination_reason) && isfinite(guard.timeout_s) && (time() - guard.start_wall_s) >= guard.timeout_s
            guard.termination_reason = "wall_clock_timeout"
            guard.termination_time_s = t
            terminate!(integrator)
            return nothing
        end

        if guard.deorbit_enabled && isempty(guard.termination_reason) && t >= guard.next_deorbit_check_t
            guard.next_deorbit_check_t = t + guard.deorbit_check_dt_s
            args = writer.args
            planet = args.environment_model.planet
            pos = SimulationEngine._state_position_ii(integrator.u, 1)
            vel = SimulationEngine._state_velocity_ii(integrator.u, 1)
            el = _elements_from_rv(pos, vel, planet.μ)
            ra_alt_m = el.ra - planet.Rp_e
            alt_m = norm(pos) - planet.Rp_e
            if isfinite(ra_alt_m) && ra_alt_m <= guard.deorbit_ra_alt_limit_m
                guard.termination_reason = "deorbit_bailout"
                guard.termination_time_s = t
                guard.termination_ra_alt_m = ra_alt_m
                guard.termination_alt_m = alt_m
                _stream_sample!(writer, integrator.u, integrator.t, integrator)
                terminate!(integrator)
            end
        end
        return nothing
    end

    return DiscreteCallback(condition, affect!)
end

function _write_active_force_history!(args, sample::NamedTuple, sol, case_dir::String)
    trajectory_csv_path = joinpath(case_dir, "simulation_results.csv")
    trajectory_feather_path = joinpath(case_dir, "simulation_results.feather")
    trajectory_df = _trajectory_dataframe(args, sol, trajectory_csv_path)
    _write_feather(trajectory_feather_path, trajectory_df)
    spacecraft = args.dynamics_model.spacecraft[1]

    fx = Vector{Float64}(undef, nrow(trajectory_df))
    fy = similar(fx)
    fz = similar(fx)
    fmag = similar(fx)
    for (i, row) in enumerate(eachrow(trajectory_df))
        t = _row_float(row, :time)
        state = _state_sample_from_row(row, spacecraft)
        force = _active_perturbation_force(args, sample, state, row, t)
        fx[i] = force[1]
        fy[i] = force[2]
        fz[i] = force[3]
        fmag[i] = norm(force)
    end

    force_df = DataFrame(
        time=trajectory_df.time,
        active_perturbation_force_ii_1=fx,
        active_perturbation_force_ii_2=fy,
        active_perturbation_force_ii_3=fz,
        active_perturbation_force_mag=fmag,
    )
    force_paths = _write_feather_table(joinpath(case_dir, "active_force_history"), force_df)

    combined_df = copy(trajectory_df)
    combined_df[!, :active_perturbation_force_ii_1] = fx
    combined_df[!, :active_perturbation_force_ii_2] = fy
    combined_df[!, :active_perturbation_force_ii_3] = fz
    combined_df[!, :active_perturbation_force_mag] = fmag
    combined_paths = _write_feather_table(joinpath(case_dir, "trajectory_with_active_force"), combined_df)

    return (
        simulation_results_csv=isfile(trajectory_csv_path) ? trajectory_csv_path : "",
        simulation_results_feather=trajectory_feather_path,
        active_force_history_csv="",
        active_force_history_feather=force_paths.feather_path,
        trajectory_with_active_force_csv="",
        trajectory_with_active_force_feather=combined_paths.feather_path,
    )
end

const _RESULT_ROW_SYMBOL_FIELDS = Set([:planet, :periapsis_regime, :dynamics_case, :density_case])

function _is_case_complete(case_dir::String)::Bool
    return isfile(joinpath(case_dir, "orbit_chunk_manifest.feather")) ||
           isfile(joinpath(case_dir, "active_force_history.feather"))
end

function _result_row_dataframe(row::NamedTuple)::DataFrame
    pairs = Pair{Symbol, Any}[]
    for key in keys(row)
        val = getfield(row, key)
        if val isa Symbol
            push!(pairs, key => String(val))
        elseif val isa Number || val isa Bool || val isa AbstractString
            push!(pairs, key => val)
        end
    end
    return DataFrame(pairs)
end

function _write_result_row(row::NamedTuple, case_dir::String)
    _write_feather(joinpath(case_dir, "result_row.feather"), _result_row_dataframe(row))
    return nothing
end

function _read_result_row(case_dir::String)::Union{NamedTuple, Nothing}
    feather_path = joinpath(case_dir, "result_row.feather")
    csv_path = joinpath(case_dir, "result_row.csv")
    df = if isfile(feather_path)
        DataFrame(Arrow.Table(feather_path))
    elseif isfile(csv_path)
        CSV.read(csv_path, DataFrame)
    else
        return nothing
    end
    nrow(df) == 0 && return nothing
    row = first(eachrow(df))
    pairs = []
    for name in names(df)
        sym = Symbol(name)
        val = getproperty(row, sym)
        val_converted = (sym in _RESULT_ROW_SYMBOL_FIELDS && val isa AbstractString) ? Symbol(val) : val
        push!(pairs, sym => val_converted)
    end
    return NamedTuple(pairs)
end

# Reconstruct a result row for a case that completed without a result_row file (e.g. runs
# predating this feature). Reads position/velocity from the trajectory feather to recompute
# orbital elements. Density metrics are unavailable and are set to NaN.
# Pass a pre-created `planet` object to avoid repeated SPICE kernel loads across many cases.
function _reconstruct_result_row(sample::NamedTuple, case_dir::String; planet=nothing)::NamedTuple
    manifest_path = joinpath(case_dir, "orbit_chunk_manifest.feather")
    legacy_traj_path = joinpath(case_dir, "simulation_results.feather")
    planet = planet !== nothing ? planet : _planet(sample.planet)
    rp_alt_m = _periapsis_altitude_m(sample.planet, sample.periapsis_regime)
    ra_alt_m = sample.apoapsis_alt_km * 1e3
    mission_time = sample.norbits * _period_s(planet, rp_alt_m, ra_alt_m)
    if !isfile(manifest_path) && !isfile(legacy_traj_path)
        return _sample_failure_row(sample; output_dir=dirname(case_dir), stage="resume_missing_trajectory")
    end
    n = 0
    first_tbl = nothing
    last_tbl = nothing
    first_idx = 1
    last_idx = 1
    if isfile(manifest_path)
        manifest = DataFrame(Arrow.Table(manifest_path))
        sort!(manifest, :orbit)
        n = sum(Int.(manifest.rows))
        first_path = String(manifest.path[1])
        last_path = String(manifest.path[end])
        first_tbl = Arrow.Table(first_path)
        last_tbl = Arrow.Table(last_path)
        first_idx = 1
        last_idx = length(last_tbl.time)
    else
        first_tbl = Arrow.Table(legacy_traj_path)
        last_tbl = first_tbl
        n = length(first_tbl.time)
        first_idx = 1
        last_idx = n
    end
    first_pos = SVector{3, Float64}(first_tbl.sc1_pos_1[first_idx], first_tbl.sc1_pos_2[first_idx], first_tbl.sc1_pos_3[first_idx])
    first_vel = SVector{3, Float64}(first_tbl.sc1_vel_1[first_idx], first_tbl.sc1_vel_2[first_idx], first_tbl.sc1_vel_3[first_idx])
    last_pos  = SVector{3, Float64}(last_tbl.sc1_pos_1[last_idx], last_tbl.sc1_pos_2[last_idx], last_tbl.sc1_pos_3[last_idx])
    last_vel  = SVector{3, Float64}(last_tbl.sc1_vel_1[last_idx], last_tbl.sc1_vel_2[last_idx], last_tbl.sc1_vel_3[last_idx])
    initial_el = _elements_from_rv(first_pos, first_vel, planet.μ)
    final_el   = _elements_from_rv(last_pos,  last_vel,  planet.μ)
    return merge(sample, (
        requested_rp_alt_m=rp_alt_m,
        requested_ra_alt_m=ra_alt_m,
        mission_time_s=mission_time,
        elapsed_s=NaN,
        retcode="unknown",
        success=true,
        aero_solver_requested=String(Symbol(_sample_field(sample, :aero_solver, :unknown))),
        aero_solver_used="unknown",
        aero_solver_maxiters=0,
        aero_auto_stiff_switch_max=0,
        aero_tolerance_scale=NaN,
        aero_solver_fallback_used=false,
        study_termination_reason="",
        study_termination_time_s=NaN,
        study_termination_ra_alt_m=NaN,
        study_termination_alt_m=NaN,
        initial_a_m=initial_el.a,
        initial_e=initial_el.e,
        initial_rp_alt_m=initial_el.rp - planet.Rp_e,
        initial_ra_alt_m=initial_el.ra - planet.Rp_e,
        final_a_m=final_el.a,
        final_e=final_el.e,
        final_rp_alt_m=final_el.rp - planet.Rp_e,
        final_ra_alt_m=final_el.ra - planet.Rp_e,
        final_i_deg=rad2deg(final_el.inc),
        final_raan_deg=rad2deg(final_el.raan),
        final_argp_deg=rad2deg(final_el.argp),
        final_energy_j_kg=final_el.energy,
        delta_a_m=final_el.a - initial_el.a,
        delta_e=final_el.e - initial_el.e,
        delta_rp_alt_m=(final_el.rp - initial_el.rp),
        delta_ra_alt_m=(final_el.ra - initial_el.ra),
        samples_saved=n,
        case_dir=case_dir,
        peak_density=NaN,
        integrated_density=NaN,
        max_dynamic_pressure=NaN,
        integrated_dynamic_pressure=NaN,
        time_below_interface_s=NaN,
        simulation_results_csv="",
        simulation_results_feather=isfile(joinpath(case_dir, "simulation_results.feather")) ? joinpath(case_dir, "simulation_results.feather") : (isfile(manifest_path) ? manifest_path : legacy_traj_path),
        active_force_history_csv="",
        active_force_history_feather=isfile(joinpath(case_dir, "active_force_history.feather")) ? joinpath(case_dir, "active_force_history.feather") : (isfile(manifest_path) ? manifest_path : joinpath(case_dir, "active_force_history.feather")),
        trajectory_with_active_force_csv="",
        trajectory_with_active_force_feather=isfile(joinpath(case_dir, "trajectory_with_active_force.feather")) ? joinpath(case_dir, "trajectory_with_active_force.feather") : (isfile(manifest_path) ? manifest_path : joinpath(case_dir, "trajectory_with_active_force.feather")),
        trajectory_chunk_manifest=isfile(manifest_path) ? manifest_path : "",
        orbit_chunk_count=isfile(manifest_path) ? nrow(DataFrame(Arrow.Table(manifest_path))) : 0,
        worker_id=myid(),
        worker_rss_after_gc_mb=NaN,
        failure_stage="",
        failure_error="",
        failure_log="",
    ))
end

function _run_sample_once(sample::NamedTuple; output_dir::String=joinpath(REPO_ROOT, "output"), results::Bool=false)
    total_elapsed = time()
    _log_sample(sample, "START $(_sample_label(sample))")
    case_dir = joinpath(output_dir, _sample_slug(sample))
    mkpath(case_dir)
    args = _make_config(sample, case_dir, results)
    planet = args.environment_model.planet
    rp_alt_m = _periapsis_altitude_m(sample.planet, sample.periapsis_regime)
    ra_alt_m = sample.apoapsis_alt_km * 1e3
    period_s = _period_s(planet, rp_alt_m, ra_alt_m)
    save_times = _aero_save_times(planet, rp_alt_m, ra_alt_m, args.environment_model.EI * 1e3, sample.norbits)
    stream_writer = OrbitChunkWriter(args, sample, case_dir, period_s)
    stream_callback = _streaming_save_callback(stream_writer, save_times)
    guard_state = StudyGuardState(sample, args)
    guard_callback = _study_guard_callback(guard_state, stream_writer)
    _log_sample(
        sample,
        @sprintf(
            "CONFIG ready mission_time=%.1f s rp_alt=%.1f km ra_alt=%.1f km save_points=%d timeout=%.1f min deorbit_limit=%.1f km",
            args.mission_configuration.mission_time,
            rp_alt_m / 1e3,
            ra_alt_m / 1e3,
            length(save_times),
            isfinite(guard_state.timeout_s) ? guard_state.timeout_s / 60.0 : 0.0,
            guard_state.deorbit_ra_alt_limit_m / 1e3,
        ),
    )
    _log_sample(sample, "SOLVE start")
    elapsed_s = withenv(
        "SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false",
        "SPACEAGORA_SOLVER_SAVE_ON" => "false",
        "SPACEAGORA_SOLVER_SAVE_START" => "false",
        "SPACEAGORA_SOLVER_SAVE_END" => "false",
    ) do
        @elapsed begin
            run_simulation(
                args;
                return_solution=false,
                isolate_state=false,
                extra_callbacks=(stream_callback, guard_callback),
            )
        end
    end
    if guard_state.termination_reason == "wall_clock_timeout"
        error("Case exceeded wall-clock timeout of $(guard_state.timeout_s / 60.0) minutes.")
    end
    retcode = "Success"
    if !isempty(guard_state.termination_reason)
        retcode = "Terminated:$(guard_state.termination_reason)"
    end
    force_paths = _finish_streaming_writer!(stream_writer)
    _log_sample(sample, @sprintf("SOLVE done elapsed=%.2f s retcode=%s streamed=%d chunks=%d", elapsed_s, retcode, force_paths.samples_saved, force_paths.orbit_chunk_count))
    stream_writer.first_pos === nothing && error("Streaming writer did not capture an initial state.")
    stream_writer.last_pos === nothing && error("Streaming writer did not capture a final state.")
    first_pos = stream_writer.first_pos
    first_vel = stream_writer.first_vel
    final_pos = stream_writer.last_pos
    final_vel = stream_writer.last_vel
    initial_el = _elements_from_rv(first_pos, first_vel, planet.μ)
    final_el = _elements_from_rv(final_pos, final_vel, planet.μ)
    _log_sample(sample, "ORBIT chunks written $(force_paths.trajectory_chunk_manifest)")
    _log_sample(
        sample,
        @sprintf(
            "DONE total=%.2f s final_ra_alt=%.3f km final_rp_alt=%.3f km",
            time() - total_elapsed,
            (final_el.ra - planet.Rp_e) / 1e3,
            (final_el.rp - planet.Rp_e) / 1e3,
        ),
    )

    return merge(sample, (
        requested_rp_alt_m=rp_alt_m,
        requested_ra_alt_m=ra_alt_m,
        mission_time_s=args.mission_configuration.mission_time,
        elapsed_s=elapsed_s,
        retcode=retcode,
        success=true,
        aero_solver_requested=String(Symbol(_sample_field(sample, :aero_solver_requested, _sample_field(sample, :aero_solver, :default)))),
        aero_solver_used=String(_aero_solver_for_config(sample)),
        aero_solver_maxiters=_aero_solver_maxiters(sample),
        aero_auto_stiff_switch_max=_aero_auto_stiff_switch_max(sample),
        aero_tolerance_scale=_aero_tolerance_scale(sample),
        aero_solver_fallback_used=Bool(_sample_field(sample, :aero_solver_fallback_used, false)),
        study_termination_reason=guard_state.termination_reason,
        study_termination_time_s=guard_state.termination_time_s,
        study_termination_ra_alt_m=guard_state.termination_ra_alt_m,
        study_termination_alt_m=guard_state.termination_alt_m,
        initial_a_m=initial_el.a,
        initial_e=initial_el.e,
        initial_rp_alt_m=initial_el.rp - planet.Rp_e,
        initial_ra_alt_m=initial_el.ra - planet.Rp_e,
        final_a_m=final_el.a,
        final_e=final_el.e,
        final_rp_alt_m=final_el.rp - planet.Rp_e,
        final_ra_alt_m=final_el.ra - planet.Rp_e,
        final_i_deg=rad2deg(final_el.inc),
        final_raan_deg=rad2deg(final_el.raan),
        final_argp_deg=rad2deg(final_el.argp),
        final_energy_j_kg=final_el.energy,
        delta_a_m=final_el.a - initial_el.a,
        delta_e=final_el.e - initial_el.e,
        delta_rp_alt_m=(final_el.rp - initial_el.rp),
        delta_ra_alt_m=(final_el.ra - initial_el.ra),
        samples_saved=force_paths.samples_saved,
        case_dir=case_dir,
    ), force_paths)
end

function _stiff_retry_error(err, bt)::Bool
    err_text = sprint(showerror, err, bt)
    return occursin("retcode=MaxIters", err_text) ||
           occursin("retcode=DtLessThanMin", err_text) ||
           occursin("retcode=Unstable", err_text) ||
           occursin("retcode=InitialFailure", err_text)
end

function _write_solver_retry_log(case_dir::String, auto_err, auto_bt; stiff_err=nothing, stiff_bt=nothing)
    mkpath(case_dir)
    open(joinpath(case_dir, "solver_retry.log"), "w") do io
        println(io, "auto_stiff_failed=true")
        println(io, "stiff_fallback_attempted=true")
        println(io, "auto_stiff_error_begin")
        println(io, sprint(showerror, auto_err, auto_bt))
        println(io, "auto_stiff_error_end")
        if stiff_err !== nothing
            println(io, "stiff_fallback_failed=true")
            println(io, "stiff_fallback_error_begin")
            println(io, sprint(showerror, stiff_err, stiff_bt))
            println(io, "stiff_fallback_error_end")
        end
    end
    return nothing
end

function run_sample(sample::NamedTuple; output_dir::String=joinpath(REPO_ROOT, "output"), results::Bool=false)
    requested_solver = Symbol(_sample_field(sample, :aero_solver, DEFAULT_AERO_SOLVER))
    if _is_aero_case(sample.dynamics_case) && requested_solver === :auto_stiff_then_rodas5p
        case_dir = joinpath(output_dir, _sample_slug(sample))
        auto_sample = merge(sample, (
            aero_solver=:auto_stiff,
            aero_solver_requested=:auto_stiff_then_rodas5p,
            aero_solver_fallback_used=false,
        ))
        try
            return _run_sample_once(auto_sample; output_dir=output_dir, results=results)
        catch err
            bt = catch_backtrace()
            _stiff_retry_error(err, bt) || rethrow()
            _log_sample(sample, "SOLVE retry with Rodas5P after auto-stiff failure: $(sprint(showerror, err))")
            _write_solver_retry_log(case_dir, err, bt)
            _post_sample_gc!()
            stiff_sample = merge(sample, (
                aero_solver=:rodas5p,
                aero_solver_requested=:auto_stiff_then_rodas5p,
                aero_solver_fallback_used=true,
            ))
            try
                return _run_sample_once(stiff_sample; output_dir=output_dir, results=results)
            catch stiff_err
                stiff_bt = catch_backtrace()
                _write_solver_retry_log(case_dir, err, bt; stiff_err=stiff_err, stiff_bt=stiff_bt)
                throw(ErrorException(
                    "Rodas5P fallback failed after auto-stiff failure. " *
                    "auto_stiff_error=$(sprint(showerror, err)); " *
                    "rodas5p_error=$(sprint(showerror, stiff_err))"
                ))
            end
        end
    end

    _validate_aero_solver(requested_solver)
    return _run_sample_once(merge(sample, (
        aero_solver=requested_solver,
        aero_solver_requested=requested_solver,
        aero_solver_fallback_used=false,
    )); output_dir=output_dir, results=results)
end

function _sample_failure_row(
    sample::NamedTuple;
    output_dir::String,
    stage::AbstractString,
    err=nothing,
    bt=nothing,
)
    case_dir = joinpath(output_dir, _sample_slug(sample))
    mkpath(case_dir)
    err_text = err === nothing ? "" : sprint(showerror, err, bt)
    failure_path = joinpath(case_dir, "failure.log")
    open(failure_path, "w") do io
        println(io, "stage=$(stage)")
        println(io, "worker=$(myid())")
        println(io, "sample=$sample")
        isempty(err_text) || println(io, err_text)
    end
    _log_sample(sample, "FAILED stage=$(stage) $(failure_path)")
    return merge(sample, (
        requested_rp_alt_m=NaN,
        requested_ra_alt_m=NaN,
        mission_time_s=NaN,
        elapsed_s=NaN,
        retcode="failed",
        success=false,
        aero_solver_requested=String(Symbol(_sample_field(sample, :aero_solver, :unknown))),
        aero_solver_used="unknown",
        aero_solver_maxiters=0,
        aero_auto_stiff_switch_max=0,
        aero_tolerance_scale=NaN,
        aero_solver_fallback_used=false,
        study_termination_reason="",
        study_termination_time_s=NaN,
        study_termination_ra_alt_m=NaN,
        study_termination_alt_m=NaN,
        failure_stage=String(stage),
        failure_error=err_text,
        failure_log=failure_path,
        initial_a_m=NaN,
        initial_e=NaN,
        initial_rp_alt_m=NaN,
        initial_ra_alt_m=NaN,
        final_a_m=NaN,
        final_e=NaN,
        final_rp_alt_m=NaN,
        final_ra_alt_m=NaN,
        final_i_deg=NaN,
        final_raan_deg=NaN,
        final_argp_deg=NaN,
        final_energy_j_kg=NaN,
        delta_a_m=NaN,
        delta_e=NaN,
        delta_rp_alt_m=NaN,
        delta_ra_alt_m=NaN,
        samples_saved=0,
        case_dir=case_dir,
        peak_density=NaN,
        integrated_density=NaN,
        max_dynamic_pressure=NaN,
        integrated_dynamic_pressure=NaN,
        time_below_interface_s=NaN,
        simulation_results_csv="",
        simulation_results_feather="",
        active_force_history_csv="",
        active_force_history_feather="",
        trajectory_with_active_force_csv="",
        trajectory_with_active_force_feather="",
        trajectory_chunk_manifest="",
        orbit_chunk_count=0,
        worker_id=myid(),
        worker_rss_after_gc_mb=NaN,
    ))
end

function run_sample_guarded(sample::NamedTuple; output_dir::String=joinpath(REPO_ROOT, "output"), results::Bool=false)
    local out
    try
        row = run_sample(sample; output_dir=output_dir, results=results)
        out = merge(row, (failure_stage="", failure_error="", failure_log=""))
    catch err
        out = _sample_failure_row(sample; output_dir=output_dir, stage="sample_exception", err=err, bt=catch_backtrace())
    finally
        _post_sample_gc!()
    end
    out = merge(out, (worker_id=myid(), worker_rss_after_gc_mb=_current_rss_mib()))
    _write_result_row(out, out.case_dir)
    return out
end

function _run_indexed_sample_remote(sample::NamedTuple, outdir::String, results::Bool)
    return run_sample_guarded(sample; output_dir=outdir, results=results)
end

function _sample_timeout_s(sample::NamedTuple)::Float64
    timeout_min = Float64(_sample_field(sample, :case_timeout_min, DEFAULT_CASE_TIMEOUT_MIN))
    return timeout_min > 0.0 ? timeout_min * 60.0 : Inf
end

function _remote_sample_fetch_with_timeout(worker::Int, sample::NamedTuple, outdir::String, results::Bool)
    fut = remotecall(_run_indexed_sample_remote, worker, sample, outdir, results)
    timeout_s = _sample_timeout_s(sample)
    start_s = time()
    while !isready(fut)
        if isfinite(timeout_s) && (time() - start_s) >= timeout_s
            return (:timeout, nothing)
        end
        sleep(isfinite(timeout_s) ? min(5.0, max(0.1, timeout_s - (time() - start_s))) : 5.0)
    end
    return (:ok, fetch(fut))
end

function _run_samples_distributed(indexed_samples, outdir::String, results::Bool, worker_max_cases::Int, worker_max_rss_gb::Float64)
    available = _aero_remote_workers()
    isempty(available) && return map(sample -> run_sample_guarded(sample; output_dir=outdir, results=results), indexed_samples)

    rows = Vector{Any}(undef, length(indexed_samples))
    next_sample = Threads.Atomic{Int}(1)
    worker_tasks = Task[]
    case_recycle_enabled = worker_max_cases > 0
    rss_recycle_enabled = worker_max_rss_gb > 0
    worker_max_rss_mb = worker_max_rss_gb * 1024.0

    for worker in available
        task = @async begin
            current_worker = worker
            cases_on_worker = 0
            while true
                idx = Threads.atomic_add!(next_sample, 1)
                idx > length(indexed_samples) && break
                sample = indexed_samples[idx]
                try
                    status, row = _remote_sample_fetch_with_timeout(current_worker, sample, outdir, results)
                    if status === :timeout
                        timeout_min = _sample_timeout_s(sample) / 60.0
                        rows[idx] = _sample_failure_row(
                            sample;
                            output_dir=outdir,
                            stage=@sprintf("worker_%d_timeout_after_%.2f_min", current_worker, timeout_min),
                        )
                        current_worker = _recycle_aero_perturb_worker!(current_worker; reason=@sprintf("case timeout %.2f min", timeout_min))
                        cases_on_worker = 0
                        continue
                    end
                    rows[idx] = row
                    cases_on_worker += 1
                    row = rows[idx]
                    rss_mb = haskey(row, :worker_rss_after_gc_mb) ? row.worker_rss_after_gc_mb : NaN
                    should_recycle_rss = rss_recycle_enabled && isfinite(rss_mb) && rss_mb >= worker_max_rss_mb
                    should_recycle_cases = case_recycle_enabled && cases_on_worker >= worker_max_cases
                    if (should_recycle_rss || should_recycle_cases) && next_sample[] <= length(indexed_samples)
                        reason = should_recycle_rss ?
                            @sprintf("RSS %.2f GiB >= %.2f GiB", rss_mb / 1024.0, worker_max_rss_gb) :
                            "case limit $(worker_max_cases)"
                        current_worker = _recycle_aero_perturb_worker!(current_worker; reason=reason)
                        cases_on_worker = 0
                    end
                catch err
                    rows[idx] = _sample_failure_row(
                        sample;
                        output_dir=outdir,
                        stage="worker_$(current_worker)_terminated",
                        err=err,
                        bt=catch_backtrace(),
                    )
                    try
                        current_worker = _recycle_aero_perturb_worker!(current_worker; reason="worker failure")
                        cases_on_worker = 0
                    catch recycle_err
                        @warn "[aero-perturb] unable to recycle failed worker=$(current_worker)" exception=(recycle_err, catch_backtrace())
                    end
                    break
                end
            end
        end
        push!(worker_tasks, task)
    end

    foreach(wait, worker_tasks)
    for idx in eachindex(rows)
        if !isassigned(rows, idx)
            rows[idx] = _sample_failure_row(
                indexed_samples[idx];
                output_dir=outdir,
                stage="not_started_after_worker_exit",
            )
        end
    end
    return rows
end

function paired_deltas(results_df::DataFrame)
    baselines = results_df[
        (results_df.dynamics_case .== :two_body) .&
        (results_df.density_case .== :none),
        :,
    ]
    base_cols = [
        :planet,
        :periapsis_regime,
        :apoapsis_alt_km,
        :spacecraft_mass_scale,
        :inclination_deg,
        :argp_deg,
        :final_a_m,
        :final_e,
        :final_rp_alt_m,
        :final_ra_alt_m,
        :final_i_deg,
        :final_raan_deg,
        :final_argp_deg,
        :final_energy_j_kg,
        :max_dynamic_pressure,
        :integrated_dynamic_pressure,
        :time_below_interface_s,
    ]
    base = baselines[:, base_cols]
    join_keys = [:planet, :periapsis_regime, :apoapsis_alt_km, :spacecraft_mass_scale, :inclination_deg, :argp_deg]
    rename!(base, Dict(c => Symbol("baseline_", c) for c in names(base) if !(Symbol(c) in join_keys)))
    joined = leftjoin(results_df, base; on=join_keys)
    for col in (:final_a_m, :final_e, :final_rp_alt_m, :final_ra_alt_m, :final_i_deg, :final_raan_deg, :final_argp_deg, :final_energy_j_kg, :max_dynamic_pressure, :integrated_dynamic_pressure, :time_below_interface_s)
        joined[!, Symbol("delta_vs_baseline_", col)] = joined[!, col] .- joined[!, Symbol("baseline_", col)]
    end
    return joined
end

function aggregate_deltas(delta_df::DataFrame)
    grouped = groupby(delta_df, [:planet, :periapsis_regime, :apoapsis_alt_km, :spacecraft_mass_scale, :inclination_deg, :argp_deg, :dynamics_case, :density_case])
    return combine(
        grouped,
        :delta_vs_baseline_final_ra_alt_m => mean => :mean_delta_final_ra_alt_m,
        :delta_vs_baseline_final_rp_alt_m => mean => :mean_delta_final_rp_alt_m,
        :delta_vs_baseline_final_a_m => mean => :mean_delta_final_a_m,
        :delta_vs_baseline_final_e => mean => :mean_delta_final_e,
        :delta_vs_baseline_max_dynamic_pressure => mean => :mean_delta_max_dynamic_pressure,
        :delta_vs_baseline_integrated_dynamic_pressure => mean => :mean_delta_integrated_dynamic_pressure,
        :elapsed_s => mean => :mean_elapsed_s,
        :success => (x -> count(identity, x) / length(x)) => :success_rate,
    )
end

function _manifest(spec::StudySpec, samples, outdir::String)
    gram_density_env = Dict(
        name => get(ENV, name, "")
        for name in sort(collect(keys(GRAM_TRAJECTORY_DENSITY_ENV)))
    )
    return Dict(
        "created_at" => string(now()),
        "repo_root" => REPO_ROOT,
        "julia_version" => string(VERSION),
        "workers" => length(_aero_remote_workers()),
        "threads" => Threads.nthreads(),
        "output_dir" => outdir,
        "spec" => Dict(
            "planets" => string.(spec.planets),
            "periapsis_regimes" => string.(spec.periapsis_regimes),
            "apoapsis_altitudes_km" => Dict(String(k) => v for (k, v) in spec.apoapsis_altitudes_km),
            "spacecraft_mass_scales" => spec.spacecraft_mass_scales,
            "inclinations_deg" => spec.inclinations_deg,
            "argp_degs" => spec.argp_degs,
            "nominal_mass_scale" => spec.nominal_mass_scale,
            "nominal_inclination_deg" => spec.nominal_inclination_deg,
            "nominal_argp_deg" => spec.nominal_argp_deg,
            "dynamics_cases" => string.(spec.dynamics_cases),
            "density_scales" => Dict(String(k) => v for (k, v) in spec.density_scales),
            "min_apoapsis_altitude_km" => MIN_APOAPSIS_ALTITUDE_KM,
            "norbits" => spec.norbits,
            "procs" => spec.procs,
            "worker_max_cases" => spec.worker_max_cases,
            "worker_max_rss_gb" => spec.worker_max_rss_gb,
            "aero_solver" => String(spec.aero_solver),
            "aero_auto_maxiters" => spec.aero_auto_maxiters,
            "aero_stiff_maxiters" => spec.aero_stiff_maxiters,
            "aero_auto_stiff_switch_max" => spec.aero_auto_stiff_switch_max,
            "aero_stiff_tol_scale" => spec.aero_stiff_tol_scale,
            "case_timeout_min" => spec.case_timeout_min,
            "deorbit_bailout" => spec.deorbit_bailout,
            "deorbit_bailout_margin_km" => spec.deorbit_bailout_margin_km,
            "deorbit_bailout_check_dt_s" => spec.deorbit_bailout_check_dt_s,
            "dt_max_orbit_s" => spec.dt_max_orbit_s,
            "dt_max_atmosphere_s" => spec.dt_max_atmosphere_s,
            "save_simulation_csv" => spec.results,
            "save_feather" => true,
            "active_force_history" => "orbit_chunked",
            "outside_atmosphere_save_dt_s" => OUTSIDE_ATMOSPHERE_SAVE_DT_S,
            "inside_atmosphere_save_dt_s" => INSIDE_ATMOSPHERE_SAVE_DT_S,
            "orbit_chunked_outputs" => true,
            "gram_density_method" => "vacuum_predicted_gram_cache",
            "gram_density_env" => gram_density_env,
            "smoke" => spec.smoke,
            "resume" => spec.resume,
            "resume_dir" => spec.resume_dir,
            "sample_count" => length(samples),
        ),
    )
end

function _latest_run_dir(output_dir::String)::String
    isdir(output_dir) || error("[aero-perturb] output directory does not exist for --resume: $(output_dir)")
    run_dirs = filter(readdir(output_dir; join=true)) do d
        isdir(d) && occursin(r"^\d{8}_\d{6}$", basename(d))
    end
    isempty(run_dirs) && error("[aero-perturb] no timestamped runs found under $(output_dir)")
    return last(sort(run_dirs))
end

function run_study(spec::StudySpec=StudySpec())
    gram_density_env = configure_gram_trajectory_density!()
    samples = make_samples(spec)
    indexed_samples = [
        merge(sample, (sample_index=i, sample_count=length(samples)))
        for (i, sample) in enumerate(samples)
    ]

    resume_dir = if !isempty(spec.resume_dir)
        spec.resume_dir
    elseif spec.resume
        _latest_run_dir(spec.output_dir)
    else
        ""
    end

    outdir = if !isempty(resume_dir)
        isdir(resume_dir) || error("[aero-perturb] resume directory does not exist: $(resume_dir)")
        resume_dir
    else
        stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
        d = joinpath(spec.output_dir, stamp)
        mkpath(d)
        d
    end

    # When resuming, partition samples into already-complete and still-pending.
    completed_rows = Any[]
    pending_samples = eltype(indexed_samples)[]
    if !isempty(resume_dir)
        planet_cache = Dict{Symbol, Any}()
        for sample in indexed_samples
            case_dir = joinpath(outdir, _sample_slug(sample))
            if _is_case_complete(case_dir)
                row = _read_result_row(case_dir)
                if row !== nothing
                    push!(completed_rows, row)
                else
                    planet = get!(planet_cache, sample.planet) do; _planet(sample.planet); end
                    push!(completed_rows, _reconstruct_result_row(sample, case_dir; planet=planet))
                end
            else
                push!(pending_samples, sample)
            end
        end
        println("[aero-perturb] RESUME $(resume_dir)")
        println("[aero-perturb] completed=$(length(completed_rows)) pending=$(length(pending_samples)) total=$(length(samples)) workers=$(length(_aero_remote_workers())) worker_max_cases=$(spec.worker_max_cases) worker_max_rss_gb=$(spec.worker_max_rss_gb) aero_solver=$(spec.aero_solver) output=$(outdir)")
    else
        pending_samples = indexed_samples
        println("[aero-perturb] samples=$(length(samples)), norbits=$(spec.smoke ? 1 : spec.norbits), workers=$(length(_aero_remote_workers())), worker_max_cases=$(spec.worker_max_cases), worker_max_rss_gb=$(spec.worker_max_rss_gb), aero_solver=$(spec.aero_solver), output=$(outdir)")
    end

    println("[aero-perturb] GRAM density method=vacuum_predicted_gram_cache env=$(gram_density_env)")
    println("[aero-perturb] progress logging is per simulation: START -> CONFIG -> SOLVE start -> SOLVE done -> DONE")
    flush(stdout)

    new_rows = if isempty(pending_samples)
        Any[]
    elseif !isempty(_aero_remote_workers())
        _run_samples_distributed(pending_samples, outdir, spec.results, spec.worker_max_cases, spec.worker_max_rss_gb)
    else
        map(sample -> run_sample_guarded(sample; output_dir=outdir, results=spec.results), pending_samples)
    end

    all_rows = vcat(completed_rows, new_rows)
    sort!(all_rows, by = r -> r.sample_index)

    results_df = DataFrame(all_rows)
    delta_df = paired_deltas(results_df)
    aggregate_df = aggregate_deltas(delta_df)

    _write_feather_table(joinpath(outdir, "results"), results_df)
    _write_feather_table(joinpath(outdir, "paired_deltas"), delta_df)
    _write_feather_table(joinpath(outdir, "aggregates"), aggregate_df)
    open(joinpath(outdir, "manifest.toml"), "w") do io
        TOML.print(io, _manifest(spec, samples, outdir))
    end

    println("[aero-perturb] wrote results/paired_deltas/aggregates as .feather, plus manifest.toml")
    return (; outdir, results=results_df, deltas=delta_df, aggregates=aggregate_df)
end

function _parse_symbols(raw::String)::Vector{Symbol}
    isempty(strip(raw)) && return Symbol[]
    return Symbol.(strip.(split(raw, ",")))
end

function _filter_apoapsis_altitudes_km(alts::Vector{Float64})::Vector{Float64}
    filtered = [alt for alt in alts if alt >= MIN_APOAPSIS_ALTITUDE_KM]
    isempty(filtered) && throw(ArgumentError("At least one apoapsis altitude must be >= $(MIN_APOAPSIS_ALTITUDE_KM) km."))
    dropped = length(alts) - length(filtered)
    if dropped > 0
        @warn "[aero-perturb] dropped $(dropped) apoapsis altitude(s) below $(MIN_APOAPSIS_ALTITUDE_KM) km"
    end
    return filtered
end

function spec_from_args(args::Vector{String})
    planets = copy(DEFAULT_PLANETS)
    dynamics_cases = copy(DEFAULT_DYNAMICS_CASES)
    apoapsis_altitudes_km = deepcopy(DEFAULT_APOAPSIS_ALTITUDES_KM)
    spacecraft_mass_scales = copy(DEFAULT_SPACECRAFT_MASS_SCALES)
    inclinations_deg = copy(DEFAULT_INCLINATIONS_DEG)
    argp_degs = copy(DEFAULT_ARGP_DEGS)
    procs = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_PROCS", "0"))
    worker_max_cases = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_WORKER_MAX_CASES", string(DEFAULT_WORKER_MAX_CASES)))
    worker_max_rss_gb = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_WORKER_MAX_RSS_GB", string(DEFAULT_WORKER_MAX_RSS_GB)))
    aero_solver = _validate_aero_solver(Symbol(get(ENV, "SPACEAGORA_AERO_PERTURB_AERO_SOLVER", string(DEFAULT_AERO_SOLVER))))
    aero_auto_maxiters = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_AERO_AUTO_MAXITERS", string(DEFAULT_AERO_AUTO_MAXITERS)))
    aero_stiff_maxiters = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_AERO_STIFF_MAXITERS", string(DEFAULT_AERO_STIFF_MAXITERS)))
    aero_auto_stiff_switch_max = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_AERO_AUTO_STIFF_SWITCH_MAX", string(DEFAULT_AERO_AUTO_STIFF_SWITCH_MAX)))
    aero_stiff_tol_scale = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_AERO_STIFF_TOL_SCALE", string(DEFAULT_AERO_STIFF_TOL_SCALE)))
    case_timeout_min = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_CASE_TIMEOUT_MIN", string(DEFAULT_CASE_TIMEOUT_MIN)))
    deorbit_bailout = get(ENV, "SPACEAGORA_AERO_PERTURB_DEORBIT_BAILOUT", DEFAULT_DEORBIT_BAILOUT ? "1" : "0") in ("1", "true", "TRUE", "yes", "YES", "on", "ON")
    deorbit_bailout_margin_km = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_DEORBIT_BAILOUT_MARGIN_KM", string(DEFAULT_DEORBIT_BAILOUT_MARGIN_KM)))
    deorbit_bailout_check_dt_s = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_DEORBIT_BAILOUT_CHECK_DT_S", string(DEFAULT_DEORBIT_BAILOUT_CHECK_DT_S)))
    dt_max_orbit_s = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_DT_MAX_ORBIT_S", string(DEFAULT_DT_MAX_ORBIT_S)))
    dt_max_atmosphere_s = parse(Float64, get(ENV, "SPACEAGORA_AERO_PERTURB_DT_MAX_ATMOSPHERE_S", string(DEFAULT_DT_MAX_ATMOSPHERE_S)))
    norbits = 1
    output_dir = joinpath(REPO_ROOT, "output", "aerobraking_perturbation_mc")
    smoke = get(ENV, "SPACEAGORA_AERO_PERTURB_SMOKE", "0") == "1"
    results = false
    resume = false
    resume_dir = ""

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--planets"
            i += 1; planets = _parse_symbols(args[i])
        elseif arg == "--dynamics"
            i += 1; dynamics_cases = _parse_symbols(args[i])
        elseif arg == "--apoapsis-alts"
            i += 1
            alts = _filter_apoapsis_altitudes_km(parse.(Float64, strip.(split(args[i], ","))))
            for p in keys(apoapsis_altitudes_km)
                apoapsis_altitudes_km[p] = alts
            end
        elseif arg == "--mass-scales"
            i += 1; spacecraft_mass_scales = parse.(Float64, strip.(split(args[i], ",")))
        elseif arg == "--inclinations"
            i += 1; inclinations_deg = parse.(Float64, strip.(split(args[i], ",")))
        elseif arg == "--argps"
            i += 1; argp_degs = parse.(Float64, strip.(split(args[i], ",")))
        elseif arg == "--norbits"
            i += 1; norbits = parse(Int, args[i])
        elseif arg == "--procs"
            i += 1; procs = parse(Int, args[i])
        elseif arg == "--worker-max-cases"
            i += 1; worker_max_cases = parse(Int, args[i])
        elseif arg == "--worker-max-rss-gb"
            i += 1; worker_max_rss_gb = parse(Float64, args[i])
        elseif arg == "--aero-solver"
            i += 1; aero_solver = _validate_aero_solver(Symbol(args[i]))
        elseif arg == "--aero-auto-maxiters"
            i += 1; aero_auto_maxiters = parse(Int, args[i])
        elseif arg == "--aero-stiff-maxiters"
            i += 1; aero_stiff_maxiters = parse(Int, args[i])
        elseif arg == "--aero-auto-stiff-switch-max"
            i += 1; aero_auto_stiff_switch_max = parse(Int, args[i])
        elseif arg == "--aero-stiff-tol-scale"
            i += 1; aero_stiff_tol_scale = parse(Float64, args[i])
        elseif arg == "--case-timeout-min"
            i += 1; case_timeout_min = parse(Float64, args[i])
        elseif arg == "--deorbit-bailout"
            deorbit_bailout = true
        elseif arg == "--no-deorbit-bailout"
            deorbit_bailout = false
        elseif arg == "--deorbit-bailout-margin-km"
            i += 1; deorbit_bailout_margin_km = parse(Float64, args[i])
        elseif arg == "--deorbit-bailout-check-dt"
            i += 1; deorbit_bailout_check_dt_s = parse(Float64, args[i])
        elseif arg == "--dt-max-orbit"
            i += 1; dt_max_orbit_s = parse(Float64, args[i])
        elseif arg == "--dt-max-atmosphere"
            i += 1; dt_max_atmosphere_s = parse(Float64, args[i])
        elseif arg == "--output-dir"
            i += 1; output_dir = abspath(args[i])
        elseif arg == "--resume"
            resume = true
        elseif arg == "--resume-dir"
            i += 1; resume_dir = abspath(args[i])
        elseif arg == "--smoke"
            smoke = true
        elseif arg == "--save-simulation-csv"
            @warn "--save-simulation-csv is ignored; aerobraking perturbation MC now writes Feather outputs only"
            results = false
        elseif arg == "--no-save-simulation-csv"
            results = false
        elseif arg in ("-h", "--help")
            println("""
            Usage:
              julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl [options]

            Options:
              --planets mars,earth,venus,titan
              --dynamics two_body,j2,harmonics_low,srp,third_body_sun,gram_aero,full_environment
              --apoapsis-alts km1,km2,...    override apoapsis altitude grid (applied to all planets)
              --mass-scales s1,s2,...        spacecraft mass scale factors (default: 0.25,0.5,1.0,2.0,4.0)
              --inclinations deg1,deg2,...   orbital inclinations in degrees (default: 15,30,45,60,75,93,105,120,135,150,165)
              --argps deg1,deg2,...          arguments of perigee in degrees (default: 0,30,45,60,80,90,120,135,150,180,210,240,270,300,330)
              --norbits N
              --procs N                    number of worker processes
              --worker-max-cases N         recycle each worker after N cases (default: $(DEFAULT_WORKER_MAX_CASES), 0 disables)
              --worker-max-rss-gb GB       recycle workers above post-GC RSS (default: $(DEFAULT_WORKER_MAX_RSS_GB), 0 disables)
              --aero-solver MODE           auto_stiff, rodas5p, or auto_stiff_then_rodas5p (default: $(DEFAULT_AERO_SOLVER))
              --aero-auto-maxiters N       max iterations for auto-stiff aero attempt (default: $(DEFAULT_AERO_AUTO_MAXITERS))
              --aero-stiff-maxiters N      max iterations for Rodas5P aero attempt (default: $(DEFAULT_AERO_STIFF_MAXITERS))
              --aero-auto-stiff-switch-max N  AutoTsit5 switch_max for aero cases (default: $(DEFAULT_AERO_AUTO_STIFF_SWITCH_MAX))
              --aero-stiff-tol-scale X     multiply tolerances for Rodas5P aero attempts (default: $(DEFAULT_AERO_STIFF_TOL_SCALE))
              --case-timeout-min MINUTES   terminate a case after this wall-clock time (default: $(DEFAULT_CASE_TIMEOUT_MIN), 0 disables)
              --deorbit-bailout / --no-deorbit-bailout
              --deorbit-bailout-margin-km KM  terminate aero cases when osculating apoapsis <= EI + margin (default: $(DEFAULT_DEORBIT_BAILOUT_MARGIN_KM))
              --deorbit-bailout-check-dt SECONDS  simulation-time check cadence (default: $(DEFAULT_DEORBIT_BAILOUT_CHECK_DT_S))
              --dt-max-orbit SECONDS       max coast-arc solver step (default: $(DEFAULT_DT_MAX_ORBIT_S))
              --dt-max-atmosphere SECONDS  max atmosphere solver step (default: $(DEFAULT_DT_MAX_ATMOSPHERE_S))
              --output-dir PATH
              --resume                       resume from the latest run under --output-dir
              --resume-dir PATH              resume from a specific run directory
              --smoke
              --save-simulation-csv       ignored; this study writes Feather outputs only
              --no-save-simulation-csv    default; retained for compatibility
            """)
            exit(0)
        else
            throw(ArgumentError("Unknown argument '$arg'. Use --help for options."))
        end
        i += 1
    end

    return StudySpec(
        planets=planets,
        dynamics_cases=dynamics_cases,
        apoapsis_altitudes_km=apoapsis_altitudes_km,
        spacecraft_mass_scales=spacecraft_mass_scales,
        inclinations_deg=inclinations_deg,
        argp_degs=argp_degs,
        norbits=norbits,
        procs=procs,
        worker_max_cases=worker_max_cases,
        worker_max_rss_gb=worker_max_rss_gb,
        aero_solver=aero_solver,
        aero_auto_maxiters=aero_auto_maxiters,
        aero_stiff_maxiters=aero_stiff_maxiters,
        aero_auto_stiff_switch_max=aero_auto_stiff_switch_max,
        aero_stiff_tol_scale=aero_stiff_tol_scale,
        case_timeout_min=case_timeout_min,
        deorbit_bailout=deorbit_bailout,
        deorbit_bailout_margin_km=deorbit_bailout_margin_km,
        deorbit_bailout_check_dt_s=deorbit_bailout_check_dt_s,
        dt_max_orbit_s=dt_max_orbit_s,
        dt_max_atmosphere_s=dt_max_atmosphere_s,
        output_dir=output_dir,
        results=results,
        smoke=smoke,
        resume=resume,
        resume_dir=resume_dir,
    )
end

end
