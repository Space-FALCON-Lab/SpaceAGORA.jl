# Native-free near-surface payloads: loading and evaluation belong to GRAMSuite. Resolve the optional API at
# construction time so an older GRAMSuite still loads the extension for its other model families.
function EM.GRAMNearSurfaceAtmosphereModel(; kwargs...)
    isdefined(GRAMSuite, :GRAMNearSurfaceAtmosphereModel) || throw(ArgumentError(
        "The loaded GRAMSuite version does not provide the native-free near-surface API. " *
        "Use a GRAMSuite version with GRAMNearSurfaceAtmosphereModel, restart Julia, " *
        "and load GRAMSuite before constructing SpaceAGORA.GRAMNearSurfaceAtmosphereModel."
    ))
    core = GRAMSuite.GRAMNearSurfaceAtmosphereModel(; kwargs...)
    return EM.GRAMNearSurfaceAtmosphereModel(core)
end

@inline function EM.getDensity(
    model::EM.GRAMNearSurfaceAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return GRAMSuite.density_state(model.core, h, lat, lon, el_time, wind)
end

@inline function EM.getDensity(
    model::EM.GRAMNearSurfaceAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p,
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return EM.getDensity(model, h, lat, lon, el_time, wind)
end
