# Kept separate from the native adapter: this path only loads and evaluates an
# owned offline snapshot. Resolve the optional API at construction time so an
# older GRAMSuite can still load the extension for its existing model families.
if isdefined(GRAMSuite, :GRAMGridAtmosphereModel)
    function SpaceAGORA.SimulationModel.SceneVisualization._atmosphere_grid_bounds(
        model::EM.GRAMGridAtmosphereModel{GRAMSuite.GRAMGridAtmosphereModel}
    )
        grid = model.core.surrogate
        return (alt_min_m=first(grid.alt_nodes_m), alt_max_m=last(grid.alt_nodes_m),
            lat_min_rad=first(grid.lat_nodes_rad), lat_max_rad=last(grid.lat_nodes_rad))
    end
end

function EM.GRAMGridAtmosphereModel(; kwargs...)
    isdefined(GRAMSuite, :GRAMGridAtmosphereModel) || throw(ArgumentError(
        "The loaded GRAMSuite version does not provide the native-free grid API. " *
        "Use a GRAMSuite version with GRAMGridAtmosphereModel, restart Julia, " *
        "and load GRAMSuite before constructing SpaceAGORA.GRAMGridAtmosphereModel."
    ))
    core = GRAMSuite.GRAMGridAtmosphereModel(; kwargs...)
    return EM.GRAMGridAtmosphereModel(core)
end

@inline function EM.getDensity(
    model::EM.GRAMGridAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return GRAMSuite.density_state(model.core, h, lat, lon, el_time, wind)
end

@inline function EM.getDensity(
    model::EM.GRAMGridAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p,
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return EM.getDensity(model, h, lat, lon, el_time, wind)
end
