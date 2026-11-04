# Local site export uses the canonical terrain loader and interpolation contract.
const TERRAIN_JSON_MAX_BYTES = 1 << 20
const TERRAIN_INPUT_MAX_BYTES = 128 << 20
const TERRAIN_TILE_MAX_BYTES = 8 << 20
const TERRAIN_TILES_MAX_BYTES = 128 << 20
const TERRAIN_MAX_GRIDS = 64
const TERRAIN_MAX_TILES = 8192
const TERRAIN_MAX_TILE_LEVEL = 20

function _terrain_number(value, label; positive=false)
    value isa Real && !(value isa Bool) || throw(ArgumentError("$(label) must be a number"))
    result = Float64(value)
    isfinite(result) && (!positive || result > 0) || throw(ArgumentError("$(label) must be finite$(positive ? " and positive" : "")"))
    return result
end

function _terrain_integer(value, label, low, high)
    value isa Integer && !(value isa Bool) && low <= value <= high ||
        throw(ArgumentError("$(label) must be an integer in $(low):$(high)"))
    return Int(value)
end

function _terrain_string(value, label; limit=16384)
    value isa AbstractString && ncodeunits(value) <= limit || throw(ArgumentError("$(label) must be a string of at most $(limit) bytes"))
    return String(value)
end

function _terrain_json(path)
    isfile(path) || throw(ArgumentError("terrain metadata not found: $(path)"))
    filesize(path) <= TERRAIN_JSON_MAX_BYTES || throw(ArgumentError("terrain JSON exceeds 1 MiB: $(path)"))
    bytes = open(io -> read(io, TERRAIN_JSON_MAX_BYTES + 1), path)
    length(bytes) <= TERRAIN_JSON_MAX_BYTES || throw(ArgumentError("terrain JSON exceeds 1 MiB"))
    value = JSON.parse(String(bytes))
    value isa AbstractDict || throw(ArgumentError("terrain metadata must be a JSON object"))
    return value
end

function _terrain_required(meta, key)
    haskey(meta, key) || throw(ArgumentError("terrain metadata is missing $(key)"))
    return meta[key]
end

# Site inputs are a portable directory bundle. Resolve symlinks before enforcing
# that every declared file stays inside the bundle rather than embedding outside files.
function _terrain_file(base, relative)
    name = _terrain_string(relative, "terrain relative path"; limit=4096)
    !isempty(name) && !isabspath(name) && !occursin('\0', name) && !(".." in splitpath(name)) ||
        throw(ArgumentError("terrain paths must stay inside their directory"))
    path = joinpath(base, name)
    isfile(path) || throw(ArgumentError("declared terrain file is missing: $(path)"))
    resolved = realpath(path)
    back = relpath(resolved, realpath(base))
    !isabspath(back) && first(splitpath(back)) != ".." || throw(ArgumentError("terrain file resolves outside its directory"))
    return path
end

function _terrain_bounds(meta)
    meta isa AbstractDict || throw(ArgumentError("terrain bounds must be an object"))
    south, north = (_terrain_number(_terrain_required(meta, k), k) for k in ("lat_min", "lat_max"))
    west, east = (_terrain_number(_terrain_required(meta, k), k) for k in ("lon_min", "lon_max"))
    -90 <= south < north <= 90 || throw(ArgumentError("terrain latitude bounds must increase within [-90, 90]"))
    span = east - west
    0 < span < 360 || throw(ArgumentError("terrain longitude span must be positive and less than 360 degrees"))
    west = mod(west, 360.0)
    (west == 360.0 || iszero(west)) && (west = 0.0)
    east = west + span
    0 < east - west < 360 || throw(ArgumentError("terrain longitude bounds are not representable"))
    return Dict{String,Any}("lat_min"=>south, "lat_max"=>north, "lon_min"=>west, "lon_max"=>east)
end

function _terrain_jpeg_size(bytes)
    length(bytes) >= 4 && bytes[1:2] == UInt8[0xff,0xd8] && bytes[end-1:end] == UInt8[0xff,0xd9] ||
        throw(ArgumentError("terrain imagery must be a complete JPEG"))
    i = 3
    while i + 3 <= length(bytes)
        bytes[i] == 0xff || throw(ArgumentError("invalid JPEG marker"))
        while i <= length(bytes) && bytes[i] == 0xff; i += 1; end
        i + 2 <= length(bytes) || break
        marker = bytes[i]; i += 1
        marker in (0xd9, 0xda) && break
        marker in (0x01, 0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7) && continue
        size = (Int(bytes[i]) << 8) | Int(bytes[i+1])
        size >= 2 && i + size - 1 <= length(bytes) || throw(ArgumentError("invalid JPEG segment size"))
        if marker in (0xc0, 0xc1, 0xc2)
            size >= 8 || throw(ArgumentError("invalid JPEG frame"))
            return ((Int(bytes[i+5]) << 8) | Int(bytes[i+6]), (Int(bytes[i+3]) << 8) | Int(bytes[i+4]))
        end
        i += size
    end
    throw(ArgumentError("terrain JPEG has no supported frame dimensions"))
end

"""
    terrain_tiles_payload(dir, tiles_rel) -> Union{Nothing,Dict{String,Any}}

Embed the declared local JPEG quadtree, or return nothing only when tiles_rel
is nothing. A missing declared index or tile is an error. Roots use regional
planetocentric bounds, and nodes use north-to-south y and eastward x indices.
The root tile is required; sparse deeper levels inherit the nearest ancestor.
Files must stay within their index directory. Limits are 8192 tiles, level 20,
4096 pixels per side, 8 MiB per JPEG and 128 MiB total encoded-image bytes.
"""
function terrain_tiles_payload(dir::AbstractString, tiles_rel)::Union{Nothing,Dict{String,Any}}
    tiles_rel === nothing && return nothing
    path = _terrain_file(dir, tiles_rel)
    meta = _terrain_json(path)
    get(meta, "scheme", "quadtree") == "quadtree" || throw(ArgumentError("unsupported terrain tile scheme"))
    root = _terrain_bounds(_terrain_required(meta, "root"))
    tile_px = _terrain_integer(_terrain_required(meta, "tile_px"), "tile_px", 1, 4096)
    listed = _terrain_required(meta, "nodes")
    listed isa AbstractVector && 1 <= length(listed) <= TERRAIN_MAX_TILES || throw(ArgumentError("terrain index must contain 1:8192 nodes"))
    nodes = Dict{String,Any}[]; seen = Set{Tuple{Int,Int,Int}}(); total = 0
    for node in listed
        node isa AbstractDict || throw(ArgumentError("terrain tile node must be an object"))
        level = _terrain_integer(_terrain_required(node, "level"), "tile level", 0, TERRAIN_MAX_TILE_LEVEL)
        x = _terrain_integer(_terrain_required(node, "x"), "tile x", 0, (1 << level)-1)
        y = _terrain_integer(_terrain_required(node, "y"), "tile y", 0, (1 << level)-1)
        key = (level,x,y)
        key in seen && throw(ArgumentError("duplicate terrain tile $(key)")); push!(seen,key)
        mpp = _terrain_number(_terrain_required(node, "m_per_px"), "tile m_per_px"; positive=true)
        file = _terrain_file(dirname(path), _terrain_required(node,"file"))
        bytes_count = filesize(file)
        0 < bytes_count <= TERRAIN_TILE_MAX_BYTES || throw(ArgumentError("terrain tile exceeds 8 MiB or is empty"))
        total += bytes_count
        total <= TERRAIN_TILES_MAX_BYTES || throw(ArgumentError("terrain imagery exceeds 128 MiB"))
        bytes = read(file)
        length(bytes) == bytes_count || throw(ArgumentError("terrain tile size changed while reading"))
        _terrain_jpeg_size(bytes) == (tile_px,tile_px) || throw(ArgumentError("terrain JPEG dimensions do not match tile_px"))
        push!(nodes, Dict{String,Any}("level"=>level,"x"=>x,"y"=>y,"m_per_px"=>mpp,"url"=>_data_url(bytes,"image/jpeg")))
    end
    (0,0,0) in seen || throw(ArgumentError("terrain imagery needs a level-zero root tile"))
    max_level = maximum(n->n["level"],nodes)
    if haskey(meta,"max_level")
        _terrain_integer(meta["max_level"],"max_level",0,TERRAIN_MAX_TILE_LEVEL) == max_level || throw(ArgumentError("max_level does not match the declared tiles"))
    end
    # A root must still have distinct sample centres at the deepest relief level.
    for (low,high) in ((root["lat_min"],root["lat_max"]),(root["lon_min"],root["lon_max"]))
        (high-low)/(2.0^(max_level+3)*16) >= 2*eps(max(abs(low),abs(high))) || throw(ArgumentError("terrain root is too narrow for the declared tile depth"))
    end
    out = Dict{String,Any}("scheme"=>"quadtree","root"=>root,"tile_px"=>tile_px,"max_level"=>max_level,"nodes"=>nodes)
    if haskey(meta,"source"); out["source"] = _terrain_string(meta["source"],"tile source"); end
    if haskey(meta,"attribution")
        credits=meta["attribution"]
        credits isa AbstractVector && length(credits)<=64 || throw(ArgumentError("tile attribution must be a list of at most 64 strings"))
        out["attribution"] = [_terrain_string(v,"tile attribution") for v in credits]
    end
    if haskey(meta,"resolution")
        input=meta["resolution"]; input isa AbstractDict || throw(ArgumentError("tile resolution must be an object"))
        resolution=Dict{String,Any}()
        for key in ("finest_m_per_px","source_grid_m_per_px","feature_scale_m")
            haskey(input,key) && (resolution[key] = input[key] === nothing ? nothing : _terrain_number(input[key],key;positive=true))
        end
        if haskey(input,"detector_sample_m")
            samples=input["detector_sample_m"]
            samples === nothing || (samples isa AbstractVector && length(samples)==2) || throw(ArgumentError("detector_sample_m must contain two samples"))
            resolution["detector_sample_m"] = samples === nothing ? nothing : [_terrain_number(v,"detector sample";positive=true) for v in samples]
        end
        haskey(input,"note") && (resolution["note"]=_terrain_string(input["note"],"resolution note"))
        out["resolution"]=resolution
    end
    return out
end

"""
    terrain_payload(site_json; max_grid=512) -> Dict{String,Any}

Load a regional site with TerrainModels.load_site_terrain and embed its DEMs
in their declared priority order. Heights remain metres above the explicitly
named reference sphere; no implicit lunar radius or datum conversion is used.
Bounds remain the original outer cell edges. When reducing a grid, sample its
canonical bilinear surface at the centres of a uniformly spaced target grid.
Both dimensions use the same reduction factor, round down, and remain at least
one. Thus exported grids may have singleton axes. Float32 values are encoded
little-endian, row-major, north to south and west to east.

max_grid must be an integer in 1:2048. Inputs are limited to 64 grids, 128 MiB
of DEM samples and 1 MiB per JSON document. Imagery is optional; declared files
must exist. Without imagery, the viewer retains height/radar queries but draws
its normal globe. Reduced grids approximate the original bilinear surface;
site.height_m records the original source value. The viewer samples the
exported grid for its marker and radar queries, exposing that approximation.
"""
function terrain_payload(site_json::AbstractString; max_grid=512)::Dict{String,Any}
    bound = _terrain_integer(max_grid,"max_grid",1,2048)
    path=abspath(String(site_json)); meta=_terrain_json(path); dir=dirname(path)
    location=_terrain_required(meta,"site")
    location isa AbstractDict || throw(ArgumentError("terrain site must be an object"))
    lat=_terrain_number(_terrain_required(location,"lat_deg"),"site latitude")
    -90 <= lat <= 90 || throw(ArgumentError("site latitude must lie within [-90,90]"))
    _terrain_number(_terrain_required(location,"lon_deg"),"site longitude")
    entries=_terrain_required(meta,"dem")
    entries isa AbstractVector && 1<=length(entries)<=TERRAIN_MAX_GRIDS || throw(ArgumentError("terrain site must contain 1:64 DEMs"))
    seen=Set{String}(); total=0
    for entry in entries
        entry isa AbstractDict || throw(ArgumentError("terrain DEM entry must be an object"))
        name=_terrain_string(_terrain_required(entry,"name"),"DEM name";limit=4096)
        name in seen && throw(ArgumentError("duplicate terrain DEM name")); push!(seen,name)
        header=_terrain_file(dir,name*".json"); raw=_terrain_file(dir,name*".f32")
        header_meta = _terrain_json(header)
        _terrain_bounds(header_meta)
        _terrain_number(_terrain_required(header_meta,"reference_radius_m"),"reference_radius_m";positive=true)
        count=filesize(raw)
        0 < count <= TERRAIN_INPUT_MAX_BYTES-total || throw(ArgumentError("terrain DEM samples exceed 128 MiB or are empty"))
        total += count
    end
    model, site = TerrainModels.load_site_terrain(path)
    grids=Dict{String,Any}[]
    radius=model.reference_radius_m
    for g in model.grids
        radius + minimum(g.heights) > 0 || throw(ArgumentError("terrain height gives a nonpositive radial distance"))
        rows,cols=size(g); factor=max(1.0,max(rows,cols)/bound)
        r2=max(1,floor(Int,rows/factor)); c2=max(1,floor(Int,cols/factor))
        values=Float32[]; sizehint!(values,r2*c2)
        for r in 1:r2, c in 1:c2
            lat=g.lat_max-(r-0.5)*(g.lat_max-g.lat_min)/r2
            lon=g.lon_min+(c-0.5)*(g.lon_max-g.lon_min)/c2
            value=Float32(TerrainModels._grid_height(g,lat,lon))
            isfinite(value) || throw(ArgumentError("resampled terrain height is not finite"))
            push!(values,value)
        end
        push!(grids,Dict{String,Any}("name"=>g.name,"rows"=>r2,"cols"=>c2,
            "source_rows"=>rows,"source_cols"=>cols,"lat_min"=>g.lat_min,"lat_max"=>g.lat_max,
            "lon_min"=>g.lon_min,"lon_max"=>g.lon_max,"heights"=>_float32_base64(values),"source"=>g.source))
    end
    # Keep site height at full source fidelity. JS computes display/radar queries
    # from the exported approximation; callers can compare both explicitly.
    return Dict{String,Any}(
        "site"=>Dict{String,Any}("lat_deg"=>site.lat_deg,"lon_deg"=>site.lon_deg,"name"=>site.name,"height_m"=>site.height_m),
        "reference_radius_m"=>radius,"fallback_height_m"=>model.fallback_height_m,
        "grids"=>grids,"tiles"=>terrain_tiles_payload(dir,get(meta,"tiles",nothing)),
        "sampling"=>"cell-centres; bilinear resampling; original outer edges",
    )
end
