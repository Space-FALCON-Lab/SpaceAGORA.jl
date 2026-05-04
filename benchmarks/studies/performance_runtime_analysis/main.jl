const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_OUTPUT_DIR = joinpath(REPO_ROOT, "output", "performance")
const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH_HARMONICS_FILE = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
const EARTH_GRAM_SURROGATE_FILE = joinpath(REPO_ROOT, "data", "GRAM_surrogate", "earth_surrogate.jls")
const PERF_BASELINE_SCENARIO = "single_j2"
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
const _PERF_THREADS_BACKEND_WARNING_EMITTED = Ref(false)

include(joinpath(REPO_ROOT, "src", "parallel", "routing", "parallel_profiles.jl"))
using .ParallelProfiles

using CSV
using DataFrames
using Dates
using Distributed
using LinearAlgebra
using Random
using Sockets
using SPICE
using StaticArrays
using Statistics

if myid() == 1
    ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
    using Plots
end

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel

function _ensure_perf_gramsuite_loaded!()
    vendored_gramsuite = joinpath(REPO_ROOT, "data", "GRAMSuite.jl")
    if Base.find_package("GRAMSuite") === nothing && isdir(vendored_gramsuite)
        pushfirst!(LOAD_PATH, vendored_gramsuite)
    end
    @eval import GRAMSuite
    return nothing
end

_ensure_perf_gramsuite_loaded!()

const EM = SimulationModel.EnvironmentModels
const PERF_GRAM_LOCK = RuntimeServices.GRAM_LOCK

function EM.GRAMAtmosphereModel(; kwargs...)
    return EM.GRAMAtmosphereModel(GRAMSuite.GRAMAtmosphereModel(; kwargs...))
end

function EM.GRAMAtmosphereModelSurrogate(;
    surrogate_file::String="",
    point_fallback_below_m::Union{Nothing, Real}=nothing,
    kwargs...
)
    base_model = EM.GRAMAtmosphereModel(; kwargs...)
    file = isempty(strip(surrogate_file)) ?
        GRAMSuite.gram_default_surrogate_file(base_model.planet_name; gram_root=base_model.gram_root) :
        GRAMSuite.resolve_path(surrogate_file)
    point_fallback = if point_fallback_below_m === nothing
        GRAMSuite.gram_default_point_fallback_below_m(base_model.planet_name)
    else
        value = Float64(point_fallback_below_m)
        value >= 0.0 || throw(ArgumentError("point_fallback_below_m must be >= 0.0 m, got $value"))
        value
    end
    return EM.GRAMAtmosphereModelSurrogate(base_model, file, point_fallback)
end

function Base.deepcopy_internal(model::EM.GRAMAtmosphereModel, stackdict::IdDict)
    haskey(stackdict, model) && return stackdict[model]
    copied = lock(PERF_GRAM_LOCK) do
        EM.GRAMAtmosphereModel(deepcopy(model.core))
    end
    stackdict[model] = copied
    return copied
end

function Base.deepcopy_internal(model::EM.GRAMAtmosphereModelSurrogate, stackdict::IdDict)
    haskey(stackdict, model) && return stackdict[model]
    copied = lock(PERF_GRAM_LOCK) do
        EM.GRAMAtmosphereModelSurrogate(
            deepcopy(model.base_model),
            model.surrogate_file,
            model.point_fallback_below_m
        )
    end
    stackdict[model] = copied
    return copied
end

function EM._gram_core_density_state(
    core::GRAMSuite.GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    lock_obj,
    vacuum_temperature::Float64
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return GRAMSuite.density_state(
        core, h, lat, lon, el_time, wind;
        lock_obj=lock_obj,
        vacuum_temperature=vacuum_temperature
    )
end

@inline function EM._gram_point_density(
    model::EM.GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    h_gram = max(h, -30.0)
    return GRAMSuite.point_density_state(model.core, h_gram, lat, lon, el_time, wind; lock_obj=PERF_GRAM_LOCK)
end

function EM.getDensity(
    model::EM.GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p::params
)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0

    if h > 2000.0e3
        return 0.0, p.args.environment_model.planet.T_ref, SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif !drag_state && !p.args.mission_configuration.keplerian
        return EM.density_polyfit(h, p)
    end

    h_gram = max(h, -30.0)
    return GRAMSuite.density_state(
        model.core,
        h_gram,
        lat,
        lon,
        el_time,
        wind;
        lock_obj=PERF_GRAM_LOCK,
        vacuum_temperature=p.args.environment_model.planet.T_ref
    )
end

function EM.getDensity(
    model::EM.GRAMAtmosphereModelSurrogate,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p::params
)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0

    if h > 2000.0e3
        return 0.0, p.args.environment_model.planet.T_ref, SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif !drag_state && !p.args.mission_configuration.keplerian
        return EM.density_polyfit(h, p)
    end

    base_model = model.base_model isa EM.GRAMAtmosphereModel ? model.base_model.core : model.base_model
    point_fallback = model.base_model isa EM.GRAMAtmosphereModel ? nothing :
        (m, h_i, lat_i, lon_i, t_i, w_i) -> EM._gram_point_density(m, h_i, lat_i, lon_i, t_i, w_i)
    h_gram = max(h, -30.0)

    return GRAMSuite.surrogate_density_state(
        base_model,
        model.surrogate_file,
        model.point_fallback_below_m,
        h_gram,
        lat,
        lon,
        el_time,
        wind;
        lock_obj=PERF_GRAM_LOCK,
        point_density_fallback=point_fallback,
        vacuum_temperature=p.args.environment_model.planet.T_ref
    )
end

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :SimulationEngine)
    include(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SimulationEngine.run_simulation
end


include(joinpath(@__DIR__, "case_catalog.jl"))
include(joinpath(@__DIR__, "measurement.jl"))
include(joinpath(@__DIR__, "reporting.jl"))
include(joinpath(@__DIR__, "cli.jl"))
