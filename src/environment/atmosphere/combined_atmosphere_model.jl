"""
    CombinedAtmosphereModel(lower, upper; handover_height_m)

Height-based composition of two frozen, native-free atmosphere snapshots. Queries below `handover_height_m` (metres
above the reference ellipsoid) go to `lower`; queries at or above it go to `upper`. Both components must be
`GRAMNearSurfaceAtmosphereModel` or `GRAMGridAtmosphereModel`.

The composition covers only what its components cover. Each keeps its own domain, refusals and winds: there is no
blending, extrapolation or fallback between them, so a query the selected component refuses throws that component's
`DomainError`. With a near-surface lower component, coverage below the handover stops 5 m above supported terrain and
excludes that component's refused regions; for the published preset these are latitudes beyond 85 degrees, volcano
flanks and positions where a component is unavailable. Density, temperature and wind can change abruptly at the
handover. A format-1 near-surface component stores no winds, so its wind is zero below the handover. An explicitly
supplied format-2 near-surface component returns its own stored winds; the composition does not blend either side.

The constructor checks compatibility where both components record it: the planet, the frozen instant and the
reference-ellipsoid radii must agree, and a grid component's height range must contain the handover. The caller must
ensure that the lower component serves heights up to the handover (a near-surface component's top is an areoid
height and is not checked here) and that each component's documented validation suits the intended use.

The composition has no validated accuracy claim of its own and adds no time evolution. `atmosphere_provenance`
records the rule and both components' provenance. The engine treats the composition like its components: every
current coordinate reaches the evaluator, and concurrent queries share the read-only snapshots. A propagation that
descends below 50 km needs a touchdown specification; otherwise the engine's default 50 km impact stop applies.

For the published Mars presets:

```julia
lower = surrogate_preset_model("mars_global_near_surface_p20_frozen_v1"; version="1.1.0")
upper = surrogate_preset_model("mars_global_upper_p20_frozen_v1"; version="1.0.0")
density_model = CombinedAtmosphereModel(lower, upper; handover_height_m=80e3)
```

Near-surface version 1.1.0 reaches 81 km areoid height, above the upper preset's 80 km ellipsoidal floor everywhere.
"""
struct CombinedAtmosphereModel{L<:_NativeFreeSnapshotComponent, U<:_NativeFreeSnapshotComponent} <: AbstractDensityModel
    lower::L
    upper::U
    handover_height_m::Float64
    function CombinedAtmosphereModel(lower::L, upper::U; handover_height_m::Real) where {
            L<:_NativeFreeSnapshotComponent, U<:_NativeFreeSnapshotComponent}
        handover_height_m isa Bool && throw(ArgumentError("handover_height_m must be a finite height in metres."))
        h = Float64(handover_height_m)
        isfinite(h) || throw(ArgumentError("handover_height_m must be a finite height in metres; got $handover_height_m."))
        _check_combined_components(_component_facts(lower), _component_facts(upper), h)
        return new{L, U}(lower, upper, h)
    end
end

function CombinedAtmosphereModel(lower, upper; handover_height_m=nothing)
    throw(ArgumentError(
        "CombinedAtmosphereModel joins native-free snapshot models (GRAMNearSurfaceAtmosphereModel or " *
        "GRAMGridAtmosphereModel); got $(typeof(lower)) below the handover and $(typeof(upper)) above it."))
end

# Snapshot models whose spatial domain must see every current coordinate: the components and their combinations.
const _NativeFreeSnapshotModel = Union{_NativeFreeSnapshotComponent, CombinedAtmosphereModel}

# What a component records about itself in its payload metadata; `nothing` for a fact it does not record.
function _component_facts(model::_NativeFreeSnapshotComponent)
    core = model.core
    metadata = hasproperty(core, :metadata) ? getproperty(core, :metadata) : nothing
    md = metadata isa AbstractDict ? metadata : Dict{String, Any}()
    planet = get(md, "planet", nothing)
    heights = nothing
    if model isa GRAMGridAtmosphereModel && hasproperty(core, :surrogate) && hasproperty(core.surrogate, :alt_nodes_m)
        nodes = core.surrogate.alt_nodes_m
        isempty(nodes) || (heights = (Float64(first(nodes)), Float64(last(nodes))))
    end
    return (planet = planet isa AbstractString ? lowercase(planet) : nothing, instant = _frozen_instant(md),
            radii_km = _ellipsoid_radii_km(md), heights_m = heights)
end

# The frozen instant as (year, month, day, hour, minute, second), from `epoch_utc` or `initial_time` metadata.
function _frozen_instant(md::AbstractDict)
    epoch = get(md, "epoch_utc", nothing)
    if epoch isa AbstractString
        m = match(r"^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)Z?$", epoch)
        m === nothing && return nothing
        return (parse.(Int, m.captures[1:5])..., parse(Float64, m.captures[6]))
    end
    t = get(md, "initial_time", nothing)
    if t isa AbstractDict && all(k -> haskey(t, k) && t[k] isa Real, ("year", "month", "day", "hour", "minute", "second"))
        return (Int(t["year"]), Int(t["month"]), Int(t["day"]), Int(t["hour"]), Int(t["minute"]), Float64(t["second"]))
    end
    return nothing
end

function _ellipsoid_radii_km(md::AbstractDict)
    config = get(md, "generation_config", nothing)
    config isa AbstractDict || return nothing
    a, b = get(config, "equatorial_radius_km", nothing), get(config, "polar_radius_km", nothing)
    return a isa Real && b isa Real ? (Float64(a), Float64(b)) : nothing
end

function _check_combined_components(lower, upper, handover_m::Float64)
    if lower.planet !== nothing && upper.planet !== nothing && lower.planet != upper.planet
        throw(ArgumentError("CombinedAtmosphereModel components are for different planets: $(lower.planet) below " *
            "the handover and $(upper.planet) above it."))
    end
    if lower.instant !== nothing && upper.instant !== nothing &&
            (lower.instant[1:5] != upper.instant[1:5] || abs(lower.instant[6] - upper.instant[6]) > 1e-6)
        throw(ArgumentError("CombinedAtmosphereModel components are frozen at different instants, $(lower.instant) " *
            "below the handover and $(upper.instant) above it; combining them would mix atmospheric states."))
    end
    if lower.radii_km !== nothing && upper.radii_km !== nothing && !all(isapprox.(lower.radii_km, upper.radii_km; rtol=1e-12))
        throw(ArgumentError("CombinedAtmosphereModel components use different reference ellipsoids: equatorial and " *
            "polar radii $(lower.radii_km) km below the handover and $(upper.radii_km) km above it."))
    end
    tol = 1e-12 * max(1.0, abs(handover_m))
    if upper.heights_m !== nothing && !(upper.heights_m[1] - tol <= handover_m < upper.heights_m[2])
        throw(ArgumentError("The upper component's grid covers $(upper.heights_m[1]) to $(upper.heights_m[2]) m above " *
            "the ellipsoid; the handover height $handover_m m must lie at or above its floor and below its ceiling."))
    end
    if lower.heights_m !== nothing && !(lower.heights_m[1] < handover_m <= lower.heights_m[2] + tol)
        throw(ArgumentError("The lower component's grid covers $(lower.heights_m[1]) to $(lower.heights_m[2]) m above " *
            "the ellipsoid; the handover height $handover_m m must lie above its floor and at or below its ceiling."))
    end
    return nothing
end

@inline function getDensity(
    model::CombinedAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
)
    return h < model.handover_height_m ? getDensity(model.lower, h, lat, lon, el_time, wind) :
                                          getDensity(model.upper, h, lat, lon, el_time, wind)
end

@inline function getDensity(
    model::CombinedAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p,
)
    return getDensity(model, h, lat, lon, el_time, wind)
end
