##
# GRAM/offline-surrogate wiring for the include-world drivers in this folder.
# The SpaceAGORAGRAMSuiteExt package extension only loads when SpaceAGORA is
# loaded as a package; scripts that `include` the source tree must wire the
# EnvironmentModels hooks themselves (same pattern as test/integration/runtests.jl),
# including the two getDensity methods the extension normally provides.
# Requires harness_48hr_common.jl to be included first (REPO_ROOT, modules).
##
if !isdefined(@__MODULE__, :_CYGNSS_GRAM_WIRING_LOADED)
const _CYGNSS_GRAM_WIRING_LOADED = true

"""
Prepare a minimal offline GRAM root so `GRAMAtmosphereModel` construction can
resolve it and fall back to the offline surrogate when no native build exists:
an empty `Build/` directory plus the Earth surrogate symlinked from
`data/GRAM_surrogate/`. No-ops when everything is already in place.
"""
function ensure_offline_gram_root!()
    root = joinpath(REPO_ROOT, "data", "GRAMSuite.jl", "GRAM Suite 2.0")
    isdir(root) || error("GRAMSuite submodule not populated at $(root)")
    mkpath(joinpath(root, "Build"))
    mkpath(joinpath(root, "Earth"))
    link = joinpath(root, "Earth", "earth_surrogate.jls")
    surrogate = joinpath(REPO_ROOT, "data", "GRAM_surrogate", "earth_surrogate.jls")
    if !isfile(link) && isfile(surrogate)
        symlink(surrogate, link)
    end
    return nothing
end
ensure_offline_gram_root!()

const HAS_GRAMSUITE = let
    vendored_gramsuite = joinpath(REPO_ROOT, "data", "GRAMSuite.jl")
    try
        if Base.find_package("GRAMSuite") === nothing && isdir(vendored_gramsuite)
            pushfirst!(LOAD_PATH, vendored_gramsuite)
        end
        @eval import GRAMSuite
        true
    catch err
        @info "GRAMSuite unavailable" exception=(err,)
        false
    end
end
HAS_GRAMSUITE || error("drag rung requires GRAMSuite")
const _GRAM_WIRE_LOCK = RuntimeServices.GRAM_LOCK
const EM = SimulationModel.EnvironmentModels
EM._GRAM_USE_GLOBAL_LOCK_FN[] = () -> GRAMSuite.gram_use_global_lock()
EM._GRAM_DEFAULT_SURROGATE_FILE_FN[] = planet -> GRAMSuite.gram_default_surrogate_file(planet)
EM._CLEAR_GRAM_STATIC_GRID_CACHE_FN[] = () -> GRAMSuite.clear_gram_static_grid_cache!()
EM._CLEAR_GRAM_OFFLINE_SURROGATE_CACHE_FN[] = () -> GRAMSuite.clear_gram_offline_surrogate_cache!()
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
        Float64(point_fallback_below_m)
    end
    return EM.GRAMAtmosphereModelSurrogate(base_model, file, point_fallback)
end
function Base.deepcopy_internal(model::EM.GRAMAtmosphereModel, stackdict::IdDict)
    haskey(stackdict, model) && return stackdict[model]
    copied = lock(_GRAM_WIRE_LOCK) do
        EM.GRAMAtmosphereModel(deepcopy(model.core))
    end
    stackdict[model] = copied
    return copied
end

@inline function _gram_call_lock_wire(model::EM.GRAMAtmosphereModel)::ReentrantLock
    return EM._gram_lock_scope() === :model ? model.instance_lock : _GRAM_WIRE_LOCK
end
# EXPERIMENT: density from NOAA WAM-IPE assimilated data sampled along the
# flight track (60 s cadence, WamIPEDensity.jl), time-interpolated.
using DelimitedFiles
const _WAM_TABLE = let
    raw, _ = readdlm("/private/tmp/claude-501/-Users-josephine--julia-dev-SpaceAGORA-jl/dd935b09-42e7-437b-8706-6710375e6226/scratchpad/wam_density_table.csv", ',', header=true)
    (t = Float64.(raw[:, 1]), rho = Float64.(raw[:, 2]))
end
@inline function _wam_rho(el_time::Float64)::Float64
    t, rho = _WAM_TABLE.t, _WAM_TABLE.rho
    el_time <= t[1] && return rho[1]
    el_time >= t[end] && return rho[end]
    i = searchsortedfirst(t, el_time)
    w = (el_time - t[i-1]) / (t[i] - t[i-1])
    return rho[i-1] * (1 - w) + rho[i] * w
end
function EM.getDensity(
    model::EM.GRAMAtmosphereModelSurrogate,
    h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params
)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0
    if h > 2000.0e3
        return 0.0, p.args.environment_model.planet.T_ref, SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif !drag_state && !p.args.mission_configuration.keplerian
        return EM.density_polyfit(h, p)
    end
    return _wam_rho(el_time), 900.0, SVector{3, Float64}(0.0, 0.0, 0.0)
end
function EM.getDensity(
    model::EM.GRAMAtmosphereModel,
    h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params
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
        model.core, h_gram, lat, lon, el_time, wind;
        lock_obj=_gram_call_lock_wire(model),
        vacuum_temperature=p.args.environment_model.planet.T_ref)
end
# --- end GRAMSuite wiring ---

end # _CYGNSS_GRAM_WIRING_LOADED guard
