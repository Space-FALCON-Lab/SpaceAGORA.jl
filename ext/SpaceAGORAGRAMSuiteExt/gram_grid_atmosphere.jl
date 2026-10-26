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

# SpaceAGORA keeps its surrogate grids under data/planets/gram_surrogate. The
# pinned GRAMSuite wrapper still probes the former adjacent data/GRAM_surrogate
# location, between its own package-local grid and the bundled GRAM Suite grids.
# Passing the relocated directory as a caller search root would move it ahead
# of the wrapper's package-local grid, so it is supplied only when that grid is
# absent. The selection order therefore stays: package-local, SpaceAGORA,
# bundled. Explicit files and caller-supplied roots are never altered.
const _SPACEAGORA_GRID_DIR = normpath(joinpath(pkgdir(SpaceAGORA), "data", "planets", "gram_surrogate"))

function _grid_planet_filename(planet)::String
    key = isdefined(GRAMSuite, :_gram_normalize_planet_name) ?
        GRAMSuite._gram_normalize_planet_name(String(planet)) : lowercase(strip(String(planet)))
    return "$(key)_surrogate.jls"
end

function _grid_search_roots(planet, search_roots;
                            package_root::AbstractString=GRAMSuite._package_root(),
                            spaceagora_dir::AbstractString=_SPACEAGORA_GRID_DIR)
    (search_roots isa AbstractVector && isempty(search_roots)) || return search_roots
    filename = _grid_planet_filename(planet)
    package_local = joinpath(package_root, "data", "GRAM_surrogate", filename)
    (ispath(package_local) || islink(package_local)) && return search_roots
    relocated = joinpath(spaceagora_dir, filename)
    (ispath(relocated) || islink(relocated)) || return search_roots
    return [String(spaceagora_dir)]
end

function EM.GRAMGridAtmosphereModel(; kwargs...)
    isdefined(GRAMSuite, :GRAMGridAtmosphereModel) || throw(ArgumentError(
        "The loaded GRAMSuite version does not provide the native-free grid API. " *
        "Use a GRAMSuite version with GRAMGridAtmosphereModel, restart Julia, " *
        "and load GRAMSuite before constructing SpaceAGORA.GRAMGridAtmosphereModel."
    ))
    kw = Dict{Symbol, Any}(kwargs)
    explicit = strip(String(get(kw, :surrogate_file, "")))
    if isempty(explicit) && haskey(kw, :planet) && kw[:planet] isa AbstractString
        kw[:search_roots] = _grid_search_roots(kw[:planet], get(kw, :search_roots, String[]))
    end
    core = GRAMSuite.GRAMGridAtmosphereModel(; kw...)
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
