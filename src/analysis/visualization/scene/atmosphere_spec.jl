# Atmosphere description for the viewer: an entry-interface shell and
# a density profile for the verified pure analytic atmosphere models.
# Scene export never queries native or arbitrary user-defined models.
# Their trajectory density comes from the existing saved density field.
# Advanced callers can explicitly sample another model through
# atmosphere_spec(...; sample_model=true), outside scene/run export.

using ..EnvironmentModels: getDensity, NoAtmosphereModel, ExponentialAtmosphereModel, PiecewiseExponentialAtmosphereModel
using ...RuntimeServices: GRAM_LOCK

const ATMOSPHERE_PROFILE_POINTS = 64
const ATMOSPHERE_MAP_STEP_DEG = 5.0

"""
    AtmosphereSpec

Atmosphere layer of a scene: `model` (type name), `ei_altitude_m`, a density
profile (`profile_altitude_m`, `profile_density_kg_m3`, empty when the model
could not be sampled), and an optional latitude/longitude density map at
`map_altitude_m` (`map_lat_deg` x `map_lon_deg`, row-major by latitude).
"""
struct AtmosphereSpec
    model::String
    ei_altitude_m::Float64
    profile_altitude_m::Vector{Float64}
    profile_density_kg_m3::Vector{Float64}
    map_altitude_m::Float64
    map_lat_deg::Vector{Float64}
    map_lon_deg::Vector{Float64}
    map_density_kg_m3::Vector{Float64}
end

Base.:(==)(a::AtmosphereSpec, b::AtmosphereSpec) = all(getfield(a, f) == getfield(b, f) for f in fieldnames(AtmosphereSpec))

@inline _altitude_only_model(model)::Bool = model isa ExponentialAtmosphereModel || model isa PiecewiseExponentialAtmosphereModel

# Density at one point through whichever getDensity form is available; `nothing` when none is.
function _density_sample(model, h_m::Float64, lat_rad::Float64, lon_rad::Float64, t_s::Float64, density_params)
    try
        if density_params !== nothing
            return lock(GRAM_LOCK) do
                Float64(getDensity(model, h_m, lat_rad, lon_rad, t_s, false, density_params)[1])
            end
        elseif hasmethod(getDensity, Tuple{typeof(model), Float64, Float64, Float64, Float64, Bool})
            return Float64(getDensity(model, h_m, lat_rad, lon_rad, t_s, false)[1])
        end
    catch
    end
    return nothing
end

"""
    atmosphere_spec(args; density_params=nothing, sample_model=false, profile_points=64, map_altitude_m=nothing, map_step_deg=5.0) -> Union{Nothing, AtmosphereSpec}

Describe the configuration's atmosphere: `nothing` for `NoAtmosphereModel`.
By default only `ExponentialAtmosphereModel` and
`PiecewiseExponentialAtmosphereModel` receive a sampled profile. Other
models retain the entry-interface shell with an empty profile and map;
no `getDensity` call is made, even when `density_params` is provided.

Advanced standalone callers may set `sample_model=true` to allow live
queries of another model. Such calls can mutate native or user model
state and are never enabled by scene export or `run_simulation`.
The profile spans the surface to 1.5 x the entry interface at latitude,
longitude and elapsed time zero. Horizontally varying models also get
a map at `map_altitude_m` (default 0.6 x EI) on a `map_step_deg` grid.
"""
function atmosphere_spec(
    args::SimulationConfiguration;
    density_params=nothing,
    sample_model::Bool=false,
    profile_points::Integer=ATMOSPHERE_PROFILE_POINTS,
    map_altitude_m::Union{Nothing, Real}=nothing,
    map_step_deg::Real=ATMOSPHERE_MAP_STEP_DEG
)::Union{Nothing, AtmosphereSpec}
    model = args.environment_model.density_model
    model isa NoAtmosphereModel && return nothing
    ei_m = Float64(args.environment_model.EI) * 1000.0
    ei_m > 0.0 || return nothing
    name = String(nameof(typeof(model)))
    map_alt = map_altitude_m === nothing ? 0.6 * ei_m : Float64(map_altitude_m)
    if !_altitude_only_model(model) && !sample_model
        return AtmosphereSpec(name, ei_m, Float64[], Float64[], map_alt, Float64[], Float64[], Float64[])
    end

    altitudes = collect(range(0.0, 1.5 * ei_m; length=max(Int(profile_points), 2)))
    profile = Float64[]
    ok = true
    for h in altitudes
        ρ = _density_sample(model, h, 0.0, 0.0, 0.0, density_params)
        if ρ === nothing || !isfinite(ρ)
            ok = false
            break
        end
        push!(profile, max(ρ, 0.0))
    end
    ok || (altitudes = Float64[]; profile = Float64[])

    lats = Float64[]
    lons = Float64[]
    grid = Float64[]
    if !_altitude_only_model(model) && ok
        step = Float64(map_step_deg)
        lats = collect(range(-90.0 + step / 2, 90.0 - step / 2; step=step))
        lons = collect(range(-180.0 + step / 2, 180.0 - step / 2; step=step))
        grid = Vector{Float64}(undef, length(lats) * length(lons))
        map_ok = true
        for (i, lat) in enumerate(lats), (j, lon) in enumerate(lons)
            ρ = _density_sample(model, map_alt, deg2rad(lat), deg2rad(lon), 0.0, density_params)
            if ρ === nothing || !isfinite(ρ)
                map_ok = false
                break
            end
            grid[(i - 1) * length(lons) + j] = max(ρ, 0.0)
        end
        if !map_ok
            lats = Float64[]; lons = Float64[]; grid = Float64[]
        end
    end
    return AtmosphereSpec(name, ei_m, altitudes, profile, map_alt, lats, lons, grid)
end

function _scene_dict(a::AtmosphereSpec)
    return Dict{String, Any}(
        "model" => a.model,
        "ei_altitude_m" => a.ei_altitude_m,
        "profile" => Dict{String, Any}("altitude_m" => copy(a.profile_altitude_m), "density_kg_m3" => copy(a.profile_density_kg_m3)),
        "map_altitude_m" => a.map_altitude_m,
        "map" => isempty(a.map_density_kg_m3) ? nothing : Dict{String, Any}(
            "altitude_m" => a.map_altitude_m,
            "lat_deg" => copy(a.map_lat_deg),
            "lon_deg" => copy(a.map_lon_deg),
            "density_kg_m3" => copy(a.map_density_kg_m3),
        ),
    )
end

function _atmosphere_from(d)::Union{Nothing, AtmosphereSpec}
    d === nothing && return nothing
    profile = get(d, "profile", Dict{String, Any}())
    m = get(d, "map", nothing)
    return AtmosphereSpec(
        String(d["model"]),
        Float64(d["ei_altitude_m"]),
        Float64[Float64(x) for x in get(profile, "altitude_m", Any[])],
        Float64[Float64(x) for x in get(profile, "density_kg_m3", Any[])],
        Float64(get(d, "map_altitude_m", m === nothing ? 0.0 : m["altitude_m"])),
        m === nothing ? Float64[] : Float64[Float64(x) for x in m["lat_deg"]],
        m === nothing ? Float64[] : Float64[Float64(x) for x in m["lon_deg"]],
        m === nothing ? Float64[] : Float64[Float64(x) for x in m["density_kg_m3"]],
    )
end
