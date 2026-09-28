# Reader for MERRA-2 native model-level granules (collection M2I3NVASM,
# `inst3_3d_asm_Nv`), the day-specific reanalysis rather than the climatology
# vendored with Earth-GRAM.
#
# Why this collection. GRAM ships the MERRA-2 *pressure-level* grid — the 42
# levels of Table 4.1 of the MERRA-2 File Specification, 1000 hPa down to
# 0.1 hPa — averaged into monthly, time-of-day climatologies. That grid is what
# `merra2.jl` reads, and its 0.1 mb ceiling near 64 km is the source of the
# GRAM/MERRA-2 blending seam that dominates this study's weighted metrics.
#
# The native grid is 72 hybrid sigma-p model layers with PTOP fixed at
# 0.01 hPa, roughly 80 km, moving the seam about 16 km above the scored band.
# The data are also per-day and 3-hourly rather than climatological, so the
# truth field carries real weather and real time evolution: neither the
# synthetic dispersion term nor the AR(1) surrogate in `truth_sources.jl` is
# needed with this source.
#
# Schema, from the MERRA-2 File Specification (GMAO Office Note No. 9):
#
#   Dimensions: longitude=576, latitude=361, level=72, time=8
#   PL   mid level pressure     Pa     (t, z, y, x)
#   T    air temperature        K
#   H    mid layer heights      m
#   QV   specific humidity      kg/kg
#   _FillValue = 1e15
#
# Two conventions matter and are easy to get backwards:
#
#  1. "the indexing for the GEOS-5 vertical coordinate system in the vertical is
#     top to bottom, i.e., layer 1 is the top layer of the atmosphere, while
#     layer LM is adjacent to the earth's surface" — the opposite of the
#     pressure-level product `merra2.jl` reads. Levels are flipped to ascending
#     height at load so both readers present the same convention.
#  2. netCDF declares these `tzyx`, and NCDatasets presents dimensions in
#     reverse, so a variable indexes as `[lon, lat, lev, time]`. The loader
#     asserts the sizes so a transposed file fails loudly instead of silently
#     producing plausible nonsense.

using NCDatasets

const MERRA2_NATIVE_FILL = 1.0e14      # spec says 1e15; anything this large is fill
const MERRA2_NATIVE_LEVELS = 72
const MERRA2_NATIVE_R_DRY = 287.0      # J/(kg K); matches p/(rho T) measured on the GRAM grid
const MERRA2_NATIVE_PTOP_PA = 1.0      # 0.01 hPa
const EARTH_RADIUS_M = 6.371e6

"""
    Merra2NativeWindow

A lat/lon/time subset of one or more native granules. Levels run bottom-up
(ascending height) regardless of the on-disk ordering.
"""
struct Merra2NativeWindow
    times::Vector{DateTime}
    lats::Vector{Float64}          # ascending
    lons::Vector{Float64}          # ascending, -180..180
    density::Array{Float64, 4}     # (lev, lat, lon, time), kg/m^3
    height::Array{Float64, 4}      # (lev, lat, lon, time), geometric m
end

"""
    merra2_native_stream(year) -> String

MERRA-2 splits its record into production streams, and the stream number is part
of every filename.
"""
function merra2_native_stream(year::Integer)::String
    year < 1980 && throw(ArgumentError("MERRA-2 begins in 1980, got $year."))
    year <= 1991 && return "100"
    year <= 2000 && return "200"
    year <= 2010 && return "300"
    return "400"
end

"""
    merra2_native_filename(date) -> String

`MERRA2_<stream>.inst3_3d_asm_Nv.<yyyymmdd>.nc4`.
"""
function merra2_native_filename(date::Date)::String
    stamp = Dates.format(date, dateformat"yyyymmdd")
    return "MERRA2_$(merra2_native_stream(Dates.year(date))).inst3_3d_asm_Nv.$stamp.nc4"
end

"""
    merra2_native_dir() -> String

Where granules are expected. Override with `SPACEAGORA_MERRA2_NATIVE_PATH`.
"""
function merra2_native_dir()::String
    haskey(ENV, "SPACEAGORA_MERRA2_NATIVE_PATH") && return ENV["SPACEAGORA_MERRA2_NATIVE_PATH"]
    return normpath(joinpath(STUDY_ROOT, "..", "..", "..", "data", "merra2_native"))
end

function merra2_native_granule_path(date::Date; dir::String=merra2_native_dir())::String
    return joinpath(dir, merra2_native_filename(date))
end

@inline _is_fill(v::Real) = !isfinite(v) || abs(v) >= MERRA2_NATIVE_FILL

"""
    _geometric_height_m(geopotential_m) -> Float64

MERRA-2's `H` is geopotential height. The geometric altitude GRAM is queried at
differs by about `1 km` at `80 km`, which is not negligible when the vertical
correlation scale is `18 km`.

If a check against real granules shows `H` is already geometric, drop this — `verify_merra2_native.jl` prints both interpretations against the GRAM profile
so the choice can be settled on data.
"""
@inline function _geometric_height_m(h::Float64)::Float64
    return EARTH_RADIUS_M * h / (EARTH_RADIUS_M - h)
end

@inline function _index_window(values::Vector{Float64}, lo::Float64, hi::Float64)
    i0 = findlast(<=(lo), values)
    i1 = findfirst(>=(hi), values)
    return (i0 === nothing ? 1 : i0):(i1 === nothing ? length(values) : i1)
end

"""
    load_merra2_native(paths, t0, t1; lat_min, lat_max, lon_min, lon_max) -> Merra2NativeWindow

Read the lat/lon/time hyperslab covering `[t0, t1]` from the given granules. A
full granule is about `2.1 GB`; a corridor-sized window is a few tens of MB, so
the subset is taken during the read rather than after.

`paths` must be in chronological order and contiguous in time.
"""
function load_merra2_native(
    paths::Vector{String},
    t0::DateTime,
    t1::DateTime;
    lat_min::Float64=-90.0,
    lat_max::Float64=90.0,
    lon_min::Float64=-180.0,
    lon_max::Float64=180.0,
)::Merra2NativeWindow
    isempty(paths) && throw(ArgumentError("No MERRA-2 native granules given."))
    for p in paths
        isfile(p) || throw(ErrorException(
            "Missing MERRA-2 native granule at $p.\n" *
            "       These are not vendored: fetch them with benchmarks/studies/edl_aerocapture_gp_uncertainty/fetch_merra2_native.sh\n" *
            "       (needs a NASA Earthdata login), or set SPACEAGORA_MERRA2_NATIVE_PATH."
        ))
    end

    times = DateTime[]
    lats = Float64[]
    lons = Float64[]
    dens_blocks = Array{Float64, 4}[]
    hgt_blocks = Array{Float64, 4}[]

    for path in paths
        NCDataset(path, "r") do ds
            all_lats = map(v -> Float64(v), Array(ds["lat"]))
            all_lons = map(v -> Float64(v), Array(ds["lon"]))
            all_times = collect(DateTime.(ds["time"][:]))

            lat_idx = _index_window(all_lats, lat_min, lat_max)
            lon_idx = _index_window(all_lons, lon_min, lon_max)
            time_idx = findall(t -> t0 - Hour(3) <= t <= t1 + Hour(3), all_times)
            isempty(time_idx) && return

            var = ds["PL"]
            nd = ndims(var)
            nd == 4 || throw(ErrorException("Expected a 4-D PL in $path, found $nd dimensions."))
            sz = size(var)
            sz[3] == MERRA2_NATIVE_LEVELS || throw(ErrorException(
                "Expected $(MERRA2_NATIVE_LEVELS) model levels on axis 3 of $path, found $(sz[3]). " *
                "NCDatasets presents `tzyx` reversed as (lon, lat, lev, time); a different layout " *
                "means this granule is not the native M2I3NVASM product."
            ))
            sz[1] == length(all_lons) && sz[2] == length(all_lats) ||
                throw(ErrorException("PL horizontal dimensions in $path do not match lon/lat."))

            # NCDatasets applies `_FillValue` masking and hands back `missing`;
            # fold that into NaN so one code path covers both masked files and
            # files that carry the raw 1e15 sentinel.
            slab(name) = map(
                v -> v === missing ? NaN : Float64(v),
                Array(ds[name][lon_idx, lat_idx, :, time_idx]),
            )
            pl = slab("PL")
            t_air = slab("T")
            h_geo = slab("H")
            qv = slab("QV")

            nlon, nlat, nlev, ntime = size(pl)
            dens = Array{Float64, 4}(undef, nlev, nlat, nlon, ntime)
            hgt = similar(dens)
            @inbounds for it in 1:ntime, k in 1:nlev, j in 1:nlat, i in 1:nlon
                # Levels arrive top-down; flip to ascending height.
                kk = nlev - k + 1
                p = pl[i, j, k, it]
                tt = t_air[i, j, k, it]
                hh = h_geo[i, j, k, it]
                q = qv[i, j, k, it]
                if _is_fill(p) || _is_fill(tt) || _is_fill(hh) || tt <= 0.0 || p <= 0.0
                    dens[kk, j, i, it] = NaN
                    hgt[kk, j, i, it] = NaN
                else
                    tv = tt * (1.0 + 0.608 * (_is_fill(q) ? 0.0 : q))
                    dens[kk, j, i, it] = p / (MERRA2_NATIVE_R_DRY * tv)
                    hgt[kk, j, i, it] = _geometric_height_m(hh)
                end
            end

            if isempty(lats)
                append!(lats, all_lats[lat_idx])
                append!(lons, all_lons[lon_idx])
            end
            append!(times, all_times[time_idx])
            push!(dens_blocks, dens)
            push!(hgt_blocks, hgt)
        end
    end

    isempty(times) && throw(ErrorException(
        "No MERRA-2 native timesteps cover $t0 to $t1 in the given granules."
    ))
    order = sortperm(times)
    density = cat(dens_blocks...; dims=4)[:, :, :, order]
    height = cat(hgt_blocks...; dims=4)[:, :, :, order]
    return Merra2NativeWindow(times[order], lats, lons, density, height)
end

"""
    load_merra2_native_span(t0, t1; dir, lat/lon bounds) -> Merra2NativeWindow

Resolve and load every granule needed to cover `[t0, t1]`. A `6 hr`
aerocapture-to-EDL gap from an `18:00 UTC` anchor crosses midnight, so this
routinely spans two days.
"""
function load_merra2_native_span(
    t0::DateTime,
    t1::DateTime;
    dir::String=merra2_native_dir(),
    lat_min::Float64=-90.0,
    lat_max::Float64=90.0,
    lon_min::Float64=-180.0,
    lon_max::Float64=180.0,
)::Merra2NativeWindow
    days = Date(t0 - Hour(3)):Day(1):Date(t1 + Hour(3))
    paths = [merra2_native_granule_path(d; dir) for d in days]
    return load_merra2_native(paths, t0, t1; lat_min, lat_max, lon_min, lon_max)
end

@inline function _bracket(values::Vector{Float64}, x::Float64)
    n = length(values)
    n == 1 && return (1, 1, 0.0)
    j = searchsortedfirst(values, x)
    j <= 1 && return (1, 1, 0.0)
    j > n && return (n, n, 0.0)
    lo = j - 1
    return (lo, j, (x - values[lo]) / (values[j] - values[lo]))
end

"""
    _column_log_density(w, ilat, ilon, itime, alt_m) -> Float64

Log-density in one grid column at one timestep, linear in geometric height
between the bracketing model layers. `NaN` outside the column's valid range.
"""
function _column_log_density(w::Merra2NativeWindow, ilat::Int, ilon::Int, itime::Int, alt_m::Float64)::Float64
    nlev = size(w.density, 1)
    lo = 0
    hi = 0
    @inbounds for k in 1:nlev
        h = w.height[k, ilat, ilon, itime]
        (isfinite(h) && isfinite(w.density[k, ilat, ilon, itime])) || continue
        if h <= alt_m
            lo = k
        elseif hi == 0
            hi = k
            break
        end
    end
    (lo == 0 || hi == 0) && return NaN
    @inbounds begin
        h0 = w.height[lo, ilat, ilon, itime]
        h1 = w.height[hi, ilat, ilon, itime]
        d0 = w.density[lo, ilat, ilon, itime]
        d1 = w.density[hi, ilat, ilon, itime]
    end
    t = h1 > h0 ? (alt_m - h0) / (h1 - h0) : 0.0
    return muladd(t, log(d1) - log(d0), log(d0))
end

"""
    merra2_native_density(w, lat_deg, lon_deg, alt_m, dt) -> Float64

Density in kg/m^3, bilinear in latitude/longitude, log-linear in geometric
height, and linear in time between the `3-hourly` analyses. `NaN` when the point
is outside the vertical domain or the loaded window.

The time interpolation is the point of this reader: it gives the study real
atmospheric evolution across the aerocapture-to-EDL gap, replacing the AR(1)
surrogate that the climatology forced.
"""
function merra2_native_density(
    w::Merra2NativeWindow, lat_deg::Float64, lon_deg::Float64, alt_m::Float64, dt::DateTime
)::Float64
    lon = mod(lon_deg + 180.0, 360.0) - 180.0
    (lat_deg < first(w.lats) || lat_deg > last(w.lats)) && return NaN
    (lon < first(w.lons) || lon > last(w.lons)) && return NaN

    ilat0, ilat1, flat = _bracket(w.lats, lat_deg)
    ilon0, ilon1, flon = _bracket(w.lons, lon)
    secs = Float64[Dates.value(t - first(w.times)) * 1.0e-3 for t in w.times]
    target = Dates.value(dt - first(w.times)) * 1.0e-3
    (target < first(secs) || target > last(secs)) && return NaN
    it0, it1, ftime = _bracket(secs, target)

    acc = 0.0
    for (it, wt) in ((it0, 1.0 - ftime), (it1, ftime))
        wt == 0.0 && continue
        for (jlat, wlat) in ((ilat0, 1.0 - flat), (ilat1, flat))
            wlat == 0.0 && continue
            for (jlon, wlon) in ((ilon0, 1.0 - flon), (ilon1, flon))
                wlon == 0.0 && continue
                lr = _column_log_density(w, jlat, jlon, it, alt_m)
                isfinite(lr) || return NaN
                acc += wt * wlat * wlon * lr
            end
        end
    end
    return exp(acc)
end

"""
    merra2_native_ceiling_m(w, lat_deg, lon_deg) -> Float64

Highest altitude with data in the columns surrounding this point, about `80 km`
for the native grid against roughly `64 km` for the pressure-level climatology.
"""
function merra2_native_ceiling_m(w::Merra2NativeWindow, lat_deg::Float64, lon_deg::Float64)::Float64
    lon = mod(lon_deg + 180.0, 360.0) - 180.0
    ilat0, ilat1, _ = _bracket(w.lats, lat_deg)
    ilon0, ilon1, _ = _bracket(w.lons, lon)
    top = Inf
    for jlat in (ilat0, ilat1), jlon in (ilon0, ilon1), it in eachindex(w.times)
        col = -Inf
        @inbounds for k in 1:size(w.height, 1)
            h = w.height[k, jlat, jlon, it]
            if isfinite(h) && isfinite(w.density[k, jlat, jlon, it])
                col = max(col, h)
            end
        end
        top = min(top, col)
    end
    return top
end
