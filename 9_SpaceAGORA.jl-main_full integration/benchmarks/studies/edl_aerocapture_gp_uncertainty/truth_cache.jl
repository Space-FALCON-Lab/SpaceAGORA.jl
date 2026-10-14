struct TruthCache
    times::Vector{DateTime}
    lats_deg::Vector{Float64}
    lons_deg::Vector{Float64}
    alts_m::Vector{Float64}
    densities::Array{Float64, 4}
end

@inline function _axis_bounds(points_a::Vector{TrajectoryPoint}, points_b::Vector{TrajectoryPoint}, getter)
    vals = Float64[getter(p) for p in points_a]
    append!(vals, getter.(points_b))
    return minimum(vals), maximum(vals)
end

function _uniform_axis(lo::Float64, hi::Float64, n::Int)
    n >= 2 || throw(ArgumentError("Axis sample count must be at least 2."))
    return collect(range(lo, hi; length=n))
end

function _time_axis(lo::DateTime, hi::DateTime, n::Int)
    n >= 2 || throw(ArgumentError("Time axis sample count must be at least 2."))
    total_ms = Dates.value(hi - lo)
    return [lo + Millisecond(round(Int, total_ms * (k - 1) / max(n - 1, 1))) for k in 1:n]
end

function build_truth_cache(
    wam_interpolator,
    pair::TrajectoryPair;
    lat_pad_deg::Float64=1.25,
    lon_pad_deg::Float64=1.25,
    alt_pad_m::Float64=5.0e3,
    time_pad_s::Float64=600.0,
    n_time::Int=7,
    n_lat::Int=9,
    n_lon::Int=9,
    n_alt::Int=40
)::TruthCache
    t0 = min(first(pair.aerocapture).dt, first(pair.edl).dt) - Millisecond(round(Int, 1000 * time_pad_s))
    t1 = max(last(pair.aerocapture).dt, last(pair.edl).dt) + Millisecond(round(Int, 1000 * time_pad_s))
    lat_lo, lat_hi = _axis_bounds(pair.aerocapture, pair.edl, p -> p.lat_deg)
    lon_lo, lon_hi = _axis_bounds(pair.aerocapture, pair.edl, p -> p.lon_deg)
    alt_lo, alt_hi = _axis_bounds(pair.aerocapture, pair.edl, p -> p.alt_m)

    times = _time_axis(t0, t1, n_time)
    lats = _uniform_axis(lat_lo - lat_pad_deg, lat_hi + lat_pad_deg, n_lat)
    lons = _uniform_axis(lon_lo - lon_pad_deg, lon_hi + lon_pad_deg, n_lon)
    alts = _uniform_axis(max(0.0, alt_lo - alt_pad_m), alt_hi + alt_pad_m, n_alt)

    n_samples = length(times) * length(lats) * length(lons) * length(alts)
    query_times = Vector{DateTime}(undef, n_samples)
    query_lats = Vector{Float64}(undef, n_samples)
    query_lons = Vector{Float64}(undef, n_samples)
    query_alts = Vector{Float64}(undef, n_samples)
    idx = 1
    @inbounds for it in eachindex(times), ila in eachindex(lats), ilo in eachindex(lons), ialt in eachindex(alts)
        query_times[idx] = times[it]
        query_lats[idx] = lats[ila]
        query_lons[idx] = lons[ilo]
        query_alts[idx] = alts[ialt]
        idx += 1
    end
    density = reshape(
        sample_wam_ipe(wam_interpolator, query_times, query_lats, query_lons, query_alts),
        length(times), length(lats), length(lons), length(alts),
    )
    return TruthCache(times, lats, lons, alts, density)
end

@inline function _search_axis(axis, value)
    value <= axis[1] && return 1, 1, 0.0
    value >= axis[end] && return length(axis), length(axis), 0.0
    hi = searchsortedfirst(axis, value)
    lo = hi - 1
    frac = (value - axis[lo]) / (axis[hi] - axis[lo])
    return lo, hi, frac
end

@inline function _lerp_val(a::Float64, b::Float64, frac::Float64)::Float64
    return muladd(frac, b - a, a)
end

function query_truth(cache::TruthCache, dt::DateTime, lat_deg::Float64, lon_deg::Float64, alt_m::Float64)::Float64
    it0, it1, ft = _search_axis(cache.times, dt)
    ila0, ila1, fla = _search_axis(cache.lats_deg, lat_deg)
    ilo0, ilo1, flo = _search_axis(cache.lons_deg, lon_deg)
    ialt0, ialt1, falt = _search_axis(cache.alts_m, alt_m)

    function _sample(it, ila, ilo, ialt)
        return cache.densities[it, ila, ilo, ialt]
    end

    c0000 = _sample(it0, ila0, ilo0, ialt0)
    c0001 = _sample(it0, ila0, ilo0, ialt1)
    c0010 = _sample(it0, ila0, ilo1, ialt0)
    c0011 = _sample(it0, ila0, ilo1, ialt1)
    c0100 = _sample(it0, ila1, ilo0, ialt0)
    c0101 = _sample(it0, ila1, ilo0, ialt1)
    c0110 = _sample(it0, ila1, ilo1, ialt0)
    c0111 = _sample(it0, ila1, ilo1, ialt1)

    c1000 = _sample(it1, ila0, ilo0, ialt0)
    c1001 = _sample(it1, ila0, ilo0, ialt1)
    c1010 = _sample(it1, ila0, ilo1, ialt0)
    c1011 = _sample(it1, ila0, ilo1, ialt1)
    c1100 = _sample(it1, ila1, ilo0, ialt0)
    c1101 = _sample(it1, ila1, ilo0, ialt1)
    c1110 = _sample(it1, ila1, ilo1, ialt0)
    c1111 = _sample(it1, ila1, ilo1, ialt1)

    v000 = _lerp_val(c0000, c0001, falt)
    v001 = _lerp_val(c0010, c0011, falt)
    v010 = _lerp_val(c0100, c0101, falt)
    v011 = _lerp_val(c0110, c0111, falt)
    v100 = _lerp_val(c1000, c1001, falt)
    v101 = _lerp_val(c1010, c1011, falt)
    v110 = _lerp_val(c1100, c1101, falt)
    v111 = _lerp_val(c1110, c1111, falt)

    w00 = _lerp_val(v000, v001, flo)
    w01 = _lerp_val(v010, v011, flo)
    w10 = _lerp_val(v100, v101, flo)
    w11 = _lerp_val(v110, v111, flo)

    u0 = _lerp_val(w00, w01, fla)
    u1 = _lerp_val(w10, w11, fla)
    return _lerp_val(u0, u1, ft)
end

function truth_along(cache::TruthCache, points::Vector{TrajectoryPoint})::Vector{Float64}
    out = Vector{Float64}(undef, length(points))
    @inbounds for i in eachindex(points)
        p = points[i]
        out[i] = query_truth(cache, p.dt, p.lat_deg, p.lon_deg, p.alt_m)
    end
    return out
end
