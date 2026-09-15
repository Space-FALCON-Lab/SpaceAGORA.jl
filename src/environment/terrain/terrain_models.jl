module TerrainModels
    # Surface terrain: height of the ground above a body's reference sphere at a
    # latitude and longitude, from digital elevation maps (DEMs). The landing
    # guidance uses it as the radar altimeter, the touchdown event stops the
    # integration on it, and the viewer displaces its globe with the same grids.
    using ..AbstractTypes: AbstractTerrainModel
    using StaticArrays
    using JSON

    export AbstractTerrainModel, NoTerrainModel, DEMGrid, DEMTerrainModel
    export terrain_height, terrain_radius, load_dem_grid, load_site_terrain, dem_grid_covers

    """
        NoTerrainModel()

    Flat reference sphere: the terrain height is zero everywhere.
    """
    struct NoTerrainModel <: AbstractTerrainModel end

    """
        DEMGrid(heights, lat_min, lat_max, lon_min, lon_max; name="", source="")

    A regular latitude/longitude grid of heights (meters above the reference
    sphere). `heights` is rows x cols, row 1 at `lat_max` (north) and column 1
    at `lon_min` (west), each cell spanning an equal angle; the bounds are the
    outer edges of the grid. Longitudes are east-positive and compared modulo
    360, so a grid written with 0..360 longitudes answers -180..180 queries.
    """
    struct DEMGrid
        heights::Matrix{Float32}
        lat_min::Float64
        lat_max::Float64
        lon_min::Float64
        lon_max::Float64
        name::String
        source::String
        function DEMGrid(heights::AbstractMatrix{<:Real}, lat_min::Real, lat_max::Real, lon_min::Real, lon_max::Real; name::AbstractString="", source::AbstractString="")
            size(heights, 1) >= 2 && size(heights, 2) >= 2 || throw(ArgumentError("a DEM grid needs at least 2 x 2 samples"))
            lat_max > lat_min || throw(ArgumentError("DEM grid: lat_max must exceed lat_min"))
            lon_max > lon_min || throw(ArgumentError("DEM grid: lon_max must exceed lon_min"))
            return new(Matrix{Float32}(heights), Float64(lat_min), Float64(lat_max), Float64(lon_min), Float64(lon_max), String(name), String(source))
        end
    end

    Base.size(g::DEMGrid) = size(g.heights)

    @inline function _lon_into(g::DEMGrid, lon_deg::Float64)::Float64
        # bring the query longitude into the grid's own longitude window
        lon = lon_deg
        while lon < g.lon_min - 1e-12 && lon + 360.0 <= g.lon_max + 1e-9
            lon += 360.0
        end
        while lon > g.lon_max + 1e-12 && lon - 360.0 >= g.lon_min - 1e-9
            lon -= 360.0
        end
        return lon
    end

    "Whether the grid covers the point (edges included)."
    function dem_grid_covers(g::DEMGrid, lat_deg::Real, lon_deg::Real)::Bool
        lon = _lon_into(g, Float64(lon_deg))
        return g.lat_min <= lat_deg <= g.lat_max && g.lon_min <= lon <= g.lon_max
    end

    "Bilinear height from one grid (meters above the sphere); NaN outside it."
    function grid_height(g::DEMGrid, lat_deg::Real, lon_deg::Real)::Float64
        dem_grid_covers(g, lat_deg, lon_deg) || return NaN
        rows, cols = size(g.heights)
        lon = _lon_into(g, Float64(lon_deg))
        # cell centers: row r spans [lat_max - r*dlat, lat_max - (r-1)*dlat]
        dlat = (g.lat_max - g.lat_min) / rows
        dlon = (g.lon_max - g.lon_min) / cols
        fr = (g.lat_max - Float64(lat_deg)) / dlat - 0.5   # fractional row index, 0 at the first center
        fc = (lon - g.lon_min) / dlon - 0.5
        fr = clamp(fr, 0.0, rows - 1.0)
        fc = clamp(fc, 0.0, cols - 1.0)
        r0 = clamp(floor(Int, fr), 0, rows - 2); c0 = clamp(floor(Int, fc), 0, cols - 2)
        tr = fr - r0; tc = fc - c0
        h00 = Float64(@inbounds g.heights[r0 + 1, c0 + 1]); h01 = Float64(@inbounds g.heights[r0 + 1, c0 + 2])
        h10 = Float64(@inbounds g.heights[r0 + 2, c0 + 1]); h11 = Float64(@inbounds g.heights[r0 + 2, c0 + 2])
        return (1 - tr) * ((1 - tc) * h00 + tc * h01) + tr * ((1 - tc) * h10 + tc * h11)
    end

    """
        DEMTerrainModel(grids; reference_radius_m, fallback_height_m=0.0)

    Terrain from one or more DEM grids tried in order (put the finest first);
    the first grid covering the point answers, and `fallback_height_m` answers
    where no grid does. `reference_radius_m` is the sphere the heights refer to.
    """
    struct DEMTerrainModel <: AbstractTerrainModel
        grids::Vector{DEMGrid}
        reference_radius_m::Float64
        fallback_height_m::Float64
        function DEMTerrainModel(grids::AbstractVector{DEMGrid}; reference_radius_m::Real, fallback_height_m::Real=0.0)
            isempty(grids) && throw(ArgumentError("DEMTerrainModel needs at least one grid"))
            reference_radius_m > 0 || throw(ArgumentError("reference_radius_m must be positive"))
            return new(collect(grids), Float64(reference_radius_m), Float64(fallback_height_m))
        end
    end

    "Height of the ground above the reference sphere (m) at planetocentric latitude and east longitude (degrees)."
    terrain_height(::NoTerrainModel, lat_deg::Real, lon_deg::Real)::Float64 = 0.0
    function terrain_height(model::DEMTerrainModel, lat_deg::Real, lon_deg::Real)::Float64
        for g in model.grids
            h = grid_height(g, lat_deg, lon_deg)
            isnan(h) || return h
        end
        return model.fallback_height_m
    end

    "Distance from the body center to the ground (m) at the point, for a body of reference radius `reference_radius_m`."
    terrain_radius(model::AbstractTerrainModel, lat_deg::Real, lon_deg::Real, reference_radius_m::Real)::Float64 =
        Float64(reference_radius_m) + terrain_height(model, lat_deg, lon_deg)
    terrain_radius(model::DEMTerrainModel, lat_deg::Real, lon_deg::Real)::Float64 = model.reference_radius_m + terrain_height(model, lat_deg, lon_deg)

    """
        load_dem_grid(json_path) -> DEMGrid

    Read a grid written by `scripts/dev/terrain/fetch_moon_site.py`: a JSON
    header (rows, cols, lat_min, lat_max, lon_min, lon_max, source) beside a
    `.f32` file of little-endian Float32 heights, row-major from north to south.
    """
    function load_dem_grid(json_path::AbstractString)::DEMGrid
        meta = JSON.parsefile(String(json_path))
        rows, cols = Int(meta["rows"]), Int(meta["cols"])
        raw = String(json_path)[1:end - 5] * ".f32"
        isfile(raw) || throw(ArgumentError("DEM heights not found beside $(json_path): $(raw)"))
        data = Vector{Float32}(undef, rows * cols)
        read!(raw, data)
        heights = permutedims(reshape(data, cols, rows))   # file is row-major
        return DEMGrid(heights, meta["lat_min"], meta["lat_max"], meta["lon_min"], meta["lon_max"];
            name=get(meta, "name", basename(json_path)), source=get(meta, "source", ""))
    end

    """
        load_site_terrain(site_json) -> (model::DEMTerrainModel, site::NamedTuple)

    Read a site directory's `site.json` (from `fetch_moon_site.py`): every DEM it
    lists, finest first, plus the site latitude, longitude and its terrain height.
    """
    function load_site_terrain(site_json::AbstractString)
        meta = JSON.parsefile(String(site_json))
        dir = dirname(String(site_json))
        grids = DEMGrid[]
        for d in meta["dem"]
            g = load_dem_grid(joinpath(dir, String(d["name"]) * ".json"))
            push!(grids, DEMGrid(g.heights, g.lat_min, g.lat_max, g.lon_min, g.lon_max; name=String(d["name"]), source=g.source))
        end
        isempty(grids) && throw(ArgumentError("$(site_json) lists no DEM"))
        radius = Float64(get(meta["dem"][1], "reference_radius_m", 1737400.0))
        model = DEMTerrainModel(grids; reference_radius_m=radius)
        lat = Float64(meta["site"]["lat_deg"]); lon = Float64(meta["site"]["lon_deg"])
        return model, (lat_deg=lat, lon_deg=lon, name=String(get(meta["site"], "name", "site")), height_m=terrain_height(model, lat, lon), reference_radius_m=radius)
    end
end # module TerrainModels
