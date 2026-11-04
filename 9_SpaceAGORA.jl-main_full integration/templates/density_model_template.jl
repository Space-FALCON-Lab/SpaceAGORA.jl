module ExampleDensityModelTemplate

using StaticArrays
using SpaceAGORA

struct ExampleDensityModel <: SpaceAGORA.AbstractDensityModel
    rho0_kg_m3::Float64
    scale_height_m::Float64
    temperature_k::Float64
end

function SpaceAGORA.getDensity(
    model::ExampleDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
)
    rho = model.rho0_kg_m3 * exp(-h / model.scale_height_m)
    wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    return rho, model.temperature_k, wind_vec
end

function SpaceAGORA.getDensity(
    model::ExampleDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p,
)
    return SpaceAGORA.getDensity(model, h, lat, lon, el_time, wind)
end

end # module ExampleDensityModelTemplate
