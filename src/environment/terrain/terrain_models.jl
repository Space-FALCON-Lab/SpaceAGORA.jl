module TerrainModels
    # Extracted from Evan's terrain query implementation in PR #121 (cc891c76).
    # This module only loads and queries regional, reference-sphere DEMs. It does
    # not connect terrain to propagation, touchdown, guidance, or visualization.
    using ..AbstractTypes: AbstractTerrainModel
    using JSON

    export AbstractTerrainModel, NoTerrainModel, DEMGrid, DEMTerrainModel
    export terrain_height, terrain_radius, load_dem_grid, load_site_terrain, dem_grid_covers

    function _finite_number(value, label)::Float64
        value isa Real || throw(ArgumentError("$(label) must be a real number"))
        result = try
            Float64(value)
        catch
            throw(ArgumentError("$(label) must be representable as Float64"))
        end
        isfinite(result) || throw(ArgumentError("$(label) must be finite"))
        return result
    end

    function _reference_radius(value)::Float64
        radius = _finite_number(value, "reference_radius_m")
        radius > 0 || throw(ArgumentError("reference_radius_m must be positive"))
        return radius
    end

    function _latitude(value)::Float64
        lat = _finite_number(value, "latitude")
        -90 <= lat <= 90 || throw(ArgumentError("latitude must be within [-90, 90] degrees"))
        return lat
    end

    _coordinates(lat, lon) = (_latitude(lat), _finite_number(lon, "longitude"))

    """
        NoTerrainModel()

    A reference sphere with zero terrain height. Queries use finite planetocentric
    latitude in [-90, 90] degrees and finite east-positive longitude in degrees.
    """
    struct NoTerrainModel <: AbstractTerrainModel end

    """
        DEMGrid(heights, lat_min, lat_max, lon_min, lon_max;
                name="", source="", reference_radius_m=nothing)

    A regional grid of finite heights in metres above a reference sphere. Bounds
    are outer cell edges, and samples are at cell centres: rows run north to south
    and columns west to east. Coordinates are planetocentric latitude and
    east-positive longitude in degrees. Latitude bounds must lie in [-90, 90].

    Longitude bounds must describe an increasing, unwrapped span smaller than
    360 degrees. The western bound is stored modulo 360 in [0, 360), preserving
    the span; the stored eastern bound can exceed 360. Queries wrap modulo 360.
    Full-globe grids and missing/nodata samples are not supported.

    Heights are copied into an owned Matrix{Float32} and validated after conversion.
    A loaded grid records its positive reference radius; an in-memory grid may
    omit it and receive its datum from DEMTerrainModel. Stored arrays must not be
    mutated after construction, since mutation bypasses these checks.
    """
    struct DEMGrid
        heights::Matrix{Float32}
        lat_min::Float64
        lat_max::Float64
        lon_min::Float64
        lon_max::Float64
        name::String
        source::String
        reference_radius_m::Union{Nothing,Float64}
        function DEMGrid(heights::AbstractMatrix{<:Real}, lat_min::Real, lat_max::Real,
                         lon_min::Real, lon_max::Real; name::AbstractString="",
                         source::AbstractString="", reference_radius_m::Union{Nothing,Real}=nothing)
            rows, cols = size(heights)
            rows >= 2 && cols >= 2 || throw(ArgumentError("a DEM grid needs at least 2 x 2 samples"))
            south, north = _latitude(lat_min), _latitude(lat_max)
            north > south || throw(ArgumentError("DEM grid: lat_max must exceed lat_min"))
            west_input = _finite_number(lon_min, "lon_min")
            east_input = _finite_number(lon_max, "lon_max")
            span = east_input - west_input
            0 < span < 360 || throw(ArgumentError("DEM longitude span must be positive and smaller than 360 degrees"))
            west = mod(west_input, 360.0)
            # Floating-point mod can round a tiny negative input to 360.
            west = west == 360.0 || iszero(west) ? 0.0 : west
            east = west + span
            0 < east - west < 360 || throw(ArgumentError("normalized DEM longitude span is not representable"))
            dlat, dlon = (north - south) / rows, (east - west) / cols
            dlat > 0 && dlon > 0 || throw(ArgumentError("DEM cell spacing is not representable"))
            # Leave room for half-cell rounding: one-ULP spacing can put two
            # neighbouring centres on the same Float64 value.
            dlat >= 2 * eps(max(abs(south), abs(north))) && dlon >= 2 * eps(east) ||
                throw(ArgumentError("adjacent DEM cell centres are not distinct at Float64 precision"))
            north - dlat / 2 > south + dlat / 2 &&
                west + dlon / 2 < east - dlon / 2 ||
                throw(ArgumentError("DEM cell centres are not distinct at Float64 precision"))
            radius = isnothing(reference_radius_m) ? nothing : _reference_radius(reference_radius_m)
            owned = Matrix{Float32}(heights)
            all(isfinite, owned) || throw(ArgumentError("DEM heights must be finite after Float32 conversion; nodata is unsupported"))
            return new(owned, south, north, west, east, String(name), String(source), radius)
        end
    end

    Base.size(g::DEMGrid) = size(g.heights)

    @inline function _lon_into(g::DEMGrid, lon_deg::Float64)::Float64
        # Test canonical intervals before shifting. Subtracting/adding the west
        # can round a point just outside a crossing-zero edge onto that edge.
        lon = mod(lon_deg, 360.0)
        if g.lon_min <= lon <= min(g.lon_max, 360.0)
            return lon
        elseif g.lon_max >= 360.0 && lon <= g.lon_max - 360.0
            return lon + 360.0
        end
        # Retain a mod result rounded to 360 as outside a grid starting at zero.
        return NaN
    end

    """
        dem_grid_covers(grid, lat_deg, lon_deg) -> Bool

    Whether a regional grid covers valid planetocentric coordinates, including
    its outer edges. Latitude is in [-90, 90] degrees; longitude wraps modulo 360.
    Invalid or nonfinite coordinates raise ArgumentError.
    """
    function dem_grid_covers(g::DEMGrid, lat_deg::Real, lon_deg::Real)::Bool
        lat, lon = _coordinates(lat_deg, lon_deg)
        wrapped = _lon_into(g, lon)
        return g.lat_min <= lat <= g.lat_max && g.lon_min <= wrapped <= g.lon_max
    end

    # Bilinear cell-centre interpolation; clamp edge half-cells to edge samples.
    # Only an outside-grid query returns NaN. Public queries validate coordinates.
    function _grid_height(g::DEMGrid, lat::Float64, lon::Float64)::Float64
        wrapped = _lon_into(g, lon)
        g.lat_min <= lat <= g.lat_max && g.lon_min <= wrapped <= g.lon_max || return NaN
        rows, cols = size(g.heights)
        dlat = (g.lat_max - g.lat_min) / rows
        dlon = (g.lon_max - g.lon_min) / cols
        fr = clamp((g.lat_max - lat) / dlat - 0.5, 0.0, rows - 1.0)
        fc = clamp((wrapped - g.lon_min) / dlon - 0.5, 0.0, cols - 1.0)
        r0 = clamp(floor(Int, fr), 0, rows - 2)
        c0 = clamp(floor(Int, fc), 0, cols - 2)
        tr, tc = fr - r0, fc - c0
        h00, h01 = Float64(g.heights[r0 + 1, c0 + 1]), Float64(g.heights[r0 + 1, c0 + 2])
        h10, h11 = Float64(g.heights[r0 + 2, c0 + 1]), Float64(g.heights[r0 + 2, c0 + 2])
        return (1 - tr) * ((1 - tc) * h00 + tc * h01) + tr * ((1 - tc) * h10 + tc * h11)
    end

    """
        DEMTerrainModel(grids; reference_radius_m, fallback_height_m=0.0)

    Copy regional DEM grids in priority order, usually finest first. The first
    covering grid supplies the height in metres; elsewhere return the finite
    fallback height. All known grid reference radii must exactly match the finite,
    positive model radius after Float64 conversion. No datum conversion is done.
    Grid arrays are copied again so the model owns its data independently.
    """
    struct DEMTerrainModel <: AbstractTerrainModel
        grids::Vector{DEMGrid}
        reference_radius_m::Float64
        fallback_height_m::Float64
        function DEMTerrainModel(grids::AbstractVector{DEMGrid}; reference_radius_m::Real,
                                 fallback_height_m::Real=0.0)
            isempty(grids) && throw(ArgumentError("DEMTerrainModel needs at least one grid"))
            radius = _reference_radius(reference_radius_m)
            fallback = _finite_number(fallback_height_m, "fallback_height_m")
            owned = DEMGrid[]
            for g in grids
                isnothing(g.reference_radius_m) || g.reference_radius_m == radius ||
                    throw(ArgumentError("DEM grid reference radius does not match the model"))
                push!(owned, DEMGrid(g.heights, g.lat_min, g.lat_max, g.lon_min, g.lon_max;
                    name=g.name, source=g.source, reference_radius_m=radius))
            end
            return new(owned, radius, fallback)
        end
    end

    """
        terrain_height(model, lat_deg, lon_deg) -> Float64

    Height in metres above the model's reference sphere at finite planetocentric
    latitude in [-90, 90] and east-positive longitude, both in degrees. DEM models
    use the first covering grid, then their fallback; NoTerrainModel returns zero.
    Custom AbstractTerrainModel subtypes may extend this function.
    """
    function terrain_height(::NoTerrainModel, lat_deg::Real, lon_deg::Real)::Float64
        _coordinates(lat_deg, lon_deg)
        return 0.0
    end

    function terrain_height(model::DEMTerrainModel, lat_deg::Real, lon_deg::Real)::Float64
        lat, lon = _coordinates(lat_deg, lon_deg)
        for g in model.grids
            h = _grid_height(g, lat, lon)
            isnan(h) || return h
        end
        return model.fallback_height_m
    end

    """
        terrain_radius(model, lat_deg, lon_deg, reference_radius_m) -> Float64
        terrain_radius(model::DEMTerrainModel, lat_deg, lon_deg) -> Float64

    Distance from the body centre in metres: reference-sphere radius plus terrain
    height. Coordinates follow terrain_height. The reference radius and resulting
    distance must be finite and positive. A DEM model requires an explicitly
    supplied radius to match its own; the three-argument form uses its stored
    radius. The four-argument form also supports custom AbstractTerrainModel types.
    """
    function terrain_radius(model::AbstractTerrainModel, lat_deg::Real, lon_deg::Real,
                            reference_radius_m::Real)::Float64
        lat, lon = _coordinates(lat_deg, lon_deg)
        radius = _reference_radius(reference_radius_m)
        model isa DEMTerrainModel && model.reference_radius_m != radius &&
            throw(ArgumentError("reference radius does not match the DEM terrain model"))
        height = _finite_number(terrain_height(model, lat, lon), "terrain height")
        result = radius + height
        isfinite(result) && result > 0 || throw(ArgumentError("terrain radius must be finite and positive"))
        return result
    end

    terrain_radius(model::DEMTerrainModel, lat_deg::Real, lon_deg::Real)::Float64 =
        terrain_radius(model, lat_deg, lon_deg, model.reference_radius_m)

    function _required(meta::AbstractDict, key::String)
        haskey(meta, key) || throw(ArgumentError("terrain metadata is missing $(key)"))
        return meta[key]
    end

    function _metadata(path::AbstractString)
        meta = JSON.parsefile(String(path))
        meta isa AbstractDict || throw(ArgumentError("terrain metadata must be a JSON object"))
        return meta
    end

    function _dimension(meta::AbstractDict, key::String)::Int
        value = _required(meta, key)
        value isa Integer && !(value isa Bool) || throw(ArgumentError("$(key) must be an integer sample count"))
        count = try
            Int(value)
        catch
            throw(ArgumentError("$(key) is not representable as Int"))
        end
        count >= 2 || throw(ArgumentError("$(key) must be at least 2"))
        return count
    end

    function _metadata_string(value, label)::String
        value isa AbstractString || throw(ArgumentError("$(label) must be a string"))
        return String(value)
    end

    function _validate_format(meta::AbstractDict, radius::Float64)
        if haskey(meta, "layout")
            layout = lowercase(strip(_metadata_string(meta["layout"], "layout")))
            parts = strip.(split(layout, ','))
            parts == ["row-major", "north to south", "west to east", "little-endian float32"] ||
                throw(ArgumentError("unsupported DEM layout; expected north-to-south row-major little-endian Float32"))
        end
        if haskey(meta, "format")
            format = lowercase(strip(_metadata_string(meta["format"], "format")))
            format in ("<f4", "float32", "little-endian float32") ||
                throw(ArgumentError("unsupported DEM format; expected little-endian Float32"))
        end
        if haskey(meta, "units")
            units = lowercase(strip(_metadata_string(meta["units"], "units")))
            if !(units in ("m", "meter", "meters", "metre", "metres"))
                match_result = match(r"^m above the ([0-9]+(?:\.[0-9]*)?(?:e[+-]?[0-9]+)?) (m|km) sphere$", units)
                isnothing(match_result) && throw(ArgumentError("unsupported DEM units; heights must be metres above the reference sphere"))
                units_radius = parse(Float64, match_result.captures[1]) * (match_result.captures[2] == "km" ? 1000.0 : 1.0)
                units_radius == radius || throw(ArgumentError("DEM units name a different reference sphere"))
            end
        end
        return nothing
    end

    """
        load_dem_grid(json_path) -> DEMGrid

    Read a local JSON header and sibling .f32 file. Required fields are integer
    rows and cols (both at least two), lat_min, lat_max, lon_min, lon_max, and an
    explicit reference_radius_m. The binary format is always little-endian Float32,
    row-major, north to south and west to east, with exactly rows*cols samples at
    cell centres. Heights must be finite metres above that reference sphere.

    Optional layout, format, and units declarations must agree with this format.
    The producer's units string, "m above the 1737.4 km sphere", is supported and
    its radius is checked. No nodata handling, reprojection, or datum conversion
    is performed. This function only reads local files.
    """
    function load_dem_grid(json_path::AbstractString)::DEMGrid
        stem, extension = splitext(String(json_path))
        lowercase(extension) == ".json" || throw(ArgumentError("DEM header must have a .json extension"))
        meta = _metadata(json_path)
        rows, cols = _dimension(meta, "rows"), _dimension(meta, "cols")
        byte_count = try
            Base.Checked.checked_mul(Base.Checked.checked_mul(rows, cols), 4)
        catch err
            err isa OverflowError || rethrow()
            throw(ArgumentError("DEM dimensions overflow the byte count"))
        end
        radius = _reference_radius(_required(meta, "reference_radius_m"))
        _validate_format(meta, radius)
        south, north = _latitude(_required(meta, "lat_min")), _latitude(_required(meta, "lat_max"))
        west = _finite_number(_required(meta, "lon_min"), "lon_min")
        east = _finite_number(_required(meta, "lon_max"), "lon_max")
        north > south || throw(ArgumentError("DEM grid: lat_max must exceed lat_min"))
        0 < east - west < 360 || throw(ArgumentError("DEM longitude span must be positive and smaller than 360 degrees"))
        raw = stem * ".f32"
        isfile(raw) || throw(ArgumentError("DEM heights not found beside $(json_path): $(raw)"))
        filesize(raw) == byte_count || throw(ArgumentError("DEM binary size does not match rows*cols*4 bytes"))
        bytes = open(raw, "r") do io
            data = read(io, byte_count)
            length(data) == byte_count && eof(io) || throw(ArgumentError("DEM binary size changed while reading"))
            data
        end
        heights = Matrix{Float32}(undef, rows, cols)
        k = 1
        for row in 1:rows, col in 1:cols
            word = UInt32(bytes[k]) | (UInt32(bytes[k + 1]) << 8) |
                (UInt32(bytes[k + 2]) << 16) | (UInt32(bytes[k + 3]) << 24)
            heights[row, col] = reinterpret(Float32, word)
            k += 4
        end
        return DEMGrid(heights, south, north, west, east;
            name=_metadata_string(get(meta, "name", basename(json_path)), "name"),
            source=_metadata_string(get(meta, "source", ""), "source"), reference_radius_m=radius)
    end

    """
        load_site_terrain(site_json) -> (model::DEMTerrainModel, site::NamedTuple)

    Read a local site.json containing site coordinates and a nonempty dem array
    in priority order. Each DEM entry names a sibling grid header (name + .json)
    and explicitly supplies reference_radius_m. Every entry and header must name
    the same sphere. No implicit lunar radius or top-level radius is inferred.
    Return the terrain model and site coordinates, name, height_m, and radius.
    """
    function load_site_terrain(site_json::AbstractString)
        meta = _metadata(site_json)
        entries = _required(meta, "dem")
        entries isa AbstractVector && !isempty(entries) || throw(ArgumentError("site must list at least one DEM"))
        site = _required(meta, "site")
        site isa AbstractDict || throw(ArgumentError("site metadata must be an object"))
        lat, lon = _coordinates(_required(site, "lat_deg"), _required(site, "lon_deg"))
        grids = DEMGrid[]
        for entry in entries
            entry isa AbstractDict || throw(ArgumentError("each DEM entry must be an object"))
            name = _metadata_string(_required(entry, "name"), "DEM name")
            isempty(name) && throw(ArgumentError("DEM name must not be empty"))
            radius = _reference_radius(_required(entry, "reference_radius_m"))
            g = load_dem_grid(joinpath(dirname(String(site_json)), name * ".json"))
            g.reference_radius_m == radius || throw(ArgumentError("site DEM reference radius does not match its header"))
            push!(grids, DEMGrid(g.heights, g.lat_min, g.lat_max, g.lon_min, g.lon_max;
                name=name, source=g.source, reference_radius_m=radius))
        end
        radius = something(first(grids).reference_radius_m)
        model = DEMTerrainModel(grids; reference_radius_m=radius)
        return model, (lat_deg=lat, lon_deg=lon,
            name=_metadata_string(get(site, "name", "site"), "site name"),
            height_m=terrain_height(model, lat, lon), reference_radius_m=radius)
    end
end # module TerrainModels
