# Reader for the MERRA-2 binary grids vendored with Earth-GRAM.
#
# Layout is taken from `Earth/source/MERRA2Data.cpp::readM2File`: nineteen 3-D
# blocks of `(pressure, latitude, longitude)` little-endian Float64, longitude
# fastest, followed by twenty 2-D surface blocks. The 3-D block order is
#
#   dens, hgt, uwnd, temp, vwnd, uvcor, dewp, vprs, rhum, sden,
#   spres, suwd, svwd, stmp, sdewp, svprs, srhum, spdav, spdsd
#
# Grid dimensions come from `MERRA2info.txt` (42 x 90 x 180 for the vendored
# 2.0-degree set) and the horizontal indexing follows
# `MERRA2::updateIndices`. Verified against the file: 19*42*90*180 + 20*90*180
# Float64 values is exactly the 106,012,800-byte file size, and `p / (rho * T)`
# recovers 287.0 J/(kg K) at every probed level.
#
# These are monthly climatologies by time-of-day slot, not specific-day
# reanalysis: `MERRA2_3hr_18Z_03.bin` is the March 18Z climatology and
# `MERRA2All_03.bin` is the March all-hours mean. Vertical coverage is the
# 1000 mb to 0.1 mb pressure range, about 0.1 km to 64 km.

const MERRA2_PRESSURES_PA = [
    100000.0, 97500.0, 95000.0, 92500.0, 90000.0, 87500.0, 85000.0, 82500.0, 80000.0,
    77500.0, 75000.0, 72500.0, 70000.0, 65000.0, 60000.0, 55000.0, 50000.0, 45000.0, 40000.0,
    35000.0, 30000.0, 25000.0, 20000.0, 15000.0, 10000.0, 7000.0, 5000.0, 4000.0,
    3000.0, 2000.0, 1000.0, 700.0, 500.0, 400.0, 300.0, 200.0, 100.0,
    70.0, 50.0, 40.0, 30.0, 10.0,
]

# Zero-based index of each needed field in the 3-D block sequence.
const MERRA2_BLOCK_DENS = 0
const MERRA2_BLOCK_HGT = 1
const MERRA2_BLOCK_SDEN = 9

# Heights below terrain are stored as this sentinel, and the corresponding
# density entries as NaN.
const MERRA2_MISSING_HEIGHT = -999.99

struct Merra2Grid
    month::Int
    hour_code::Int
    n_pres::Int
    n_lat::Int
    n_lon::Int
    dens::Array{Float64, 3}   # kg/m^3, (pressure, lat, lon)
    hgt::Array{Float64, 3}    # m
    sden::Array{Float64, 3}   # kg/m^3, interannual density standard deviation
end

"""
    merra2_data_dir() -> String

Directory holding the vendored MERRA-2 grids. Override with
`SPACEAGORA_MERRA2_PATH`.
"""
function merra2_data_dir()::String
    haskey(ENV, "SPACEAGORA_MERRA2_PATH") && return ENV["SPACEAGORA_MERRA2_PATH"]
    return normpath(joinpath(
        STUDY_ROOT, "..", "..", "..",
        "data", "GRAMSuite.jl", "GRAM Suite 2.0", "Earth", "data", "MERRA2data",
    ))
end

"""
    merra2_hour_code(hour) -> Int

GRAM's time-of-day code for a UTC hour: `1..8` select `00Z..21Z`, matching
`MERRA2::buildDataFileName`. Code `9` selects the all-hours monthly mean.
"""
function merra2_hour_code(hour::Integer)::Int
    code = 1 + round(Int, hour / 3.0)
    return code > 8 ? 1 : code
end

function merra2_file_path(month::Integer, hour_code::Integer)::String
    1 <= month <= 12 || throw(ArgumentError("MERRA-2 month must be in 1:12, got $month."))
    root = merra2_data_dir()
    mm = lpad(month, 2, '0')
    if hour_code == 9
        return joinpath(root, "All Mean", "MERRA2All_$mm.bin")
    end
    1 <= hour_code <= 8 ||
        throw(ArgumentError("MERRA-2 hour code must be in 1:9, got $hour_code."))
    zz = lpad(3 * (hour_code - 1), 2, '0') * "Z"
    return joinpath(root, zz, "MERRA2_3hr_$(zz)_$mm.bin")
end

function _read_merra2_info()::Tuple{Int, Int, Int}
    path = joinpath(merra2_data_dir(), "MERRA2info.txt")
    isfile(path) || throw(ErrorException("Missing MERRA-2 info file at $path."))
    n_pres = n_lat = n_lon = 0
    for line in eachline(path)
        m = match(r"^\s*(\w+)\s*=\s*(\d+)", line)
        m === nothing && continue
        name, value = m.captures[1], parse(Int, m.captures[2])
        name == "Pressures" && (n_pres = value)
        name == "Latitudes" && (n_lat = value)
        name == "Longitudes" && (n_lon = value)
    end
    (n_pres > 0 && n_lat > 0 && n_lon > 0) ||
        throw(ErrorException("Could not parse grid sizes from $path."))
    n_pres == length(MERRA2_PRESSURES_PA) ||
        throw(ErrorException("MERRA-2 info reports $n_pres pressure levels, expected $(length(MERRA2_PRESSURES_PA))."))
    return n_pres, n_lat, n_lon
end

function _read_block!(io::IO, dest::Array{Float64, 3}, block_index::Int, block_elems::Int)
    seek(io, block_index * block_elems * sizeof(Float64))
    read!(io, dest)
    return dest
end

const _MERRA2_CACHE = Dict{Tuple{Int, Int}, Merra2Grid}()

"""
    load_merra2_grid(month, hour_code) -> Merra2Grid

Load the density, level-height, and density-standard-deviation blocks for one
month and time-of-day slot. Cached per process; only three of the nineteen 3-D
blocks are read (about 16 MB) rather than the whole 101 MB file.
"""
function load_merra2_grid(month::Integer, hour_code::Integer)::Merra2Grid
    key = (Int(month), Int(hour_code))
    haskey(_MERRA2_CACHE, key) && return _MERRA2_CACHE[key]

    path = merra2_file_path(month, hour_code)
    isfile(path) || throw(ErrorException(
        "Missing MERRA-2 data file at $path. The vendored grids live under " *
        "data/GRAMSuite.jl/GRAM Suite 2.0/Earth/data/MERRA2data; set " *
        "SPACEAGORA_MERRA2_PATH to override."
    ))
    n_pres, n_lat, n_lon = _read_merra2_info()
    block_elems = n_pres * n_lat * n_lon
    expected = 19 * block_elems + 20 * n_lat * n_lon
    filesize(path) == expected * sizeof(Float64) || throw(ErrorException(
        "MERRA-2 file $path is $(filesize(path)) bytes, expected $(expected * sizeof(Float64)) " *
        "for a $(n_pres)x$(n_lat)x$(n_lon) grid."
    ))

    dens = Array{Float64, 3}(undef, n_lon, n_lat, n_pres)
    hgt = similar(dens)
    sden = similar(dens)
    open(path, "r") do io
        _read_block!(io, dens, MERRA2_BLOCK_DENS, block_elems)
        _read_block!(io, hgt, MERRA2_BLOCK_HGT, block_elems)
        _read_block!(io, sden, MERRA2_BLOCK_SDEN, block_elems)
    end

    # The file is C-ordered with longitude fastest, so the Julia read above puts
    # longitude first; permute to (pressure, lat, lon) to match the C++ indexing.
    grid = Merra2Grid(
        Int(month), Int(hour_code), n_pres, n_lat, n_lon,
        permutedims(dens, (3, 2, 1)),
        permutedims(hgt, (3, 2, 1)),
        permutedims(sden, (3, 2, 1)),
    )
    _MERRA2_CACHE[key] = grid
    return grid
end

@inline function _merra2_horizontal_indices(grid::Merra2Grid, lat_deg::Float64, lon_deg::Float64)
    # Follows MERRA2::updateIndices.
    grid_lat = (grid.n_lat - 1) / 180.0
    shift_lat = clamp(lat_deg, -90.0, 90.0) + 90.0
    ilat = min(floor(Int, shift_lat * grid_lat), grid.n_lat - 2)
    flat = shift_lat * grid_lat - ilat

    grid_lon = grid.n_lon / 360.0
    lon360 = mod(lon_deg, 360.0)
    ilon = min(floor(Int, lon360 * grid_lon), grid.n_lon - 1)
    flon = lon360 * grid_lon - ilon
    return ilat, flat, ilon, flon
end

"""
    merra2_column(grid, ilat, ilon, alt_m) -> (log_density, rel_sigma)

Log-density and relative density standard deviation in one grid column at a
geometric altitude, interpolated linearly in height between the bracketing
pressure levels. Returns `(NaN, NaN)` when the altitude falls outside the
column's valid range, which happens below terrain and above the 0.1 mb top.
"""
function merra2_column(grid::Merra2Grid, ilat::Int, ilon::Int, alt_m::Float64)::Tuple{Float64, Float64}
    lo = 0
    hi = 0
    @inbounds for k in 1:grid.n_pres
        h = grid.hgt[k, ilat + 1, ilon + 1]
        (isfinite(h) && h > MERRA2_MISSING_HEIGHT + 1.0) || continue
        isfinite(grid.dens[k, ilat + 1, ilon + 1]) || continue
        if h <= alt_m
            lo = k
        elseif hi == 0
            hi = k
            break
        end
    end
    (lo == 0 || hi == 0) && return (NaN, NaN)

    @inbounds begin
        h0 = grid.hgt[lo, ilat + 1, ilon + 1]
        h1 = grid.hgt[hi, ilat + 1, ilon + 1]
        d0 = grid.dens[lo, ilat + 1, ilon + 1]
        d1 = grid.dens[hi, ilat + 1, ilon + 1]
        s0 = grid.sden[lo, ilat + 1, ilon + 1]
        s1 = grid.sden[hi, ilat + 1, ilon + 1]
    end
    t = h1 > h0 ? (alt_m - h0) / (h1 - h0) : 0.0
    log_rho = muladd(t, log(d1) - log(d0), log(d0))
    rel0 = d0 > 0.0 ? s0 / d0 : 0.0
    rel1 = d1 > 0.0 ? s1 / d1 : 0.0
    rel = muladd(t, rel1 - rel0, rel0)
    return log_rho, max(rel, 0.0)
end

"""
    merra2_density(grid, lat_deg, lon_deg, alt_m) -> (density, rel_sigma)

Bilinear-in-latitude/longitude, log-linear-in-height density in kg/m^3 and the
relative interannual standard deviation. Returns `(NaN, NaN)` when any of the
four surrounding columns has no data at this altitude, which is the signal that
the point is outside MERRA-2's vertical domain.
"""
function merra2_density(grid::Merra2Grid, lat_deg::Float64, lon_deg::Float64, alt_m::Float64)::Tuple{Float64, Float64}
    ilat, flat, ilon, flon = _merra2_horizontal_indices(grid, lat_deg, lon_deg)
    ilon_next = (ilon + 1) % grid.n_lon

    log_rho = 0.0
    rel = 0.0
    for (jlat, wlat) in ((ilat, 1.0 - flat), (ilat + 1, flat))
        wlat == 0.0 && continue
        for (jlon, wlon) in ((ilon, 1.0 - flon), (ilon_next, flon))
            wlon == 0.0 && continue
            lr, rs = merra2_column(grid, jlat, jlon, alt_m)
            isfinite(lr) || return (NaN, NaN)
            w = wlat * wlon
            log_rho += w * lr
            rel += w * rs
        end
    end
    return exp(log_rho), rel
end

"""
    merra2_ceiling_m(grid, lat_deg, lon_deg) -> Float64

Highest altitude with MERRA-2 data in the column containing this point, about
62-65 km depending on location.
"""
function merra2_ceiling_m(grid::Merra2Grid, lat_deg::Float64, lon_deg::Float64)::Float64
    ilat, _, ilon, _ = _merra2_horizontal_indices(grid, lat_deg, lon_deg)
    ilon_next = (ilon + 1) % grid.n_lon
    top = Inf
    for jlat in (ilat, ilat + 1), jlon in (ilon, ilon_next)
        col_top = -Inf
        @inbounds for k in 1:grid.n_pres
            h = grid.hgt[k, jlat + 1, jlon + 1]
            if isfinite(h) && h > MERRA2_MISSING_HEIGHT + 1.0 && isfinite(grid.dens[k, jlat + 1, jlon + 1])
                col_top = max(col_top, h)
            end
        end
        top = min(top, col_top)
    end
    return top
end
