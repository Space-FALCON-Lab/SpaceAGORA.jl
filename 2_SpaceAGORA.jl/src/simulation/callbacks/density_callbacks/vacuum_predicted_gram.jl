# Vacuum-predicted GRAM density cache.
#
# On atmosphere entry (or when the cache is stale), propagates a drag-free
# (two-body + J2) trajectory forward for a configurable horizon, queries the
# density model at N evenly-spaced knots, and fits natural cubic splines to
# log(ρ) and T as functions of time.  Subsequent density queries within a
# configurable position deviation of the vacuum prediction are served from the
# splines instead of calling the density model directly.  When the actual
# trajectory deviates beyond the threshold the cache is rebuilt from the
# current state.
#
# Env vars:
#   SPACEAGORA_VACUUM_GRAM_CACHE             – on/off (default: off)
#   SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS     – spline knot count (default: 20)
#   SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S   – look-ahead in seconds (default: 600)
#   SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M – position deviation threshold m (default: 5000)

@inline function _vacuum_gram_cache_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_VACUUM_GRAM_CACHE", false)
end

@inline function _vacuum_gram_cache_npoints()::Int
    raw = strip(get(ENV, "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS", "20"))
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS must be an integer, got '$raw'"))
    end
    return max(4, parsed)
end

@inline _vacuum_gram_cache_horizon_s()::Float64 =
    max(10.0, _parse_float_env("SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S", 600.0))

@inline _vacuum_gram_cache_deviation_m()::Float64 =
    max(100.0, _parse_float_env("SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M", 5000.0))

@inline function VacuumPredictedGRAMCache()
    return VacuumPredictedGRAMCache(
        false, 0.0, 0.0, 1.0,
        Float64[], Float64[], Float64[], Float64[],
        SVector{3, Float64}[], Float64[], SVector{3, Float64}[]
    )
end

@inline function _vacuum_gram_cache_for_sat!(
    caches::Vector{Union{Nothing, VacuumPredictedGRAMCache}},
    sat_idx::Int
)::VacuumPredictedGRAMCache
    if sat_idx <= length(caches)
        cache = @inbounds caches[sat_idx]
        if cache === nothing
            cache = VacuumPredictedGRAMCache()
            @inbounds caches[sat_idx] = cache
        end
        return cache
    end
    return VacuumPredictedGRAMCache()
end

# Two-body + J2 gravitational acceleration in the inertial frame.
# Inlined here to avoid importing a private symbol from GravityEffectors.
@inline function _vacuum_j2_accel(pos::SVector{3, Float64}, planet)::SVector{3, Float64}
    r = norm(pos)
    μ  = Float64(planet.μ)
    J2 = Float64(planet.J2)
    # J2 reference radius is equatorial, consistent with the gravity effector.
    Re = Float64(planet.Rp_e)
    x, y, z = pos
    r2 = r * r
    a_sph = (-μ / r2) * normalize(pos)
    j2_scale = 1.5 * J2 * μ * Re^2 / (r2 * r2)
    j2_term = SVector{3, Float64}(
        x / r * (5.0 * z^2 / r2 - 1.0),
        y / r * (5.0 * z^2 / r2 - 1.0),
        z / r * (5.0 * z^2 / r2 - 3.0),
    )
    return a_sph + j2_scale * j2_term
end

@inline function _vacuum_rk4_step(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet,
    dt::Float64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    k1r = vel
    k1v = _vacuum_j2_accel(pos, planet)
    k2r = vel + 0.5 * dt * k1v
    k2v = _vacuum_j2_accel(pos + 0.5 * dt * k1r, planet)
    k3r = vel + 0.5 * dt * k2v
    k3v = _vacuum_j2_accel(pos + 0.5 * dt * k2r, planet)
    k4r = vel + dt * k3v
    k4v = _vacuum_j2_accel(pos + dt * k3r, planet)
    new_pos = pos + (dt / 6.0) * (k1r + 2.0 * k2r + 2.0 * k3r + k4r)
    new_vel = vel + (dt / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v)
    return new_pos, new_vel
end

# Build natural cubic spline second-derivative coefficients (Ms) for function
# values ys at uniformly-spaced knots with spacing h, using the Thomas algorithm
# for the tridiagonal natural-boundary system [1 4 1] * M = 6/h² * Δ²y.
function _natural_cubic_spline_build!(
    Ms::Vector{Float64},
    ys::Vector{Float64},
    h::Float64
)
    n = length(ys)
    resize!(Ms, n)
    fill!(Ms, 0.0)
    n <= 2 && return
    m = n - 2  # number of interior unknowns M[2]..M[n-1]
    m < 1 && return

    scale = 6.0 / (h * h)

    cp = Vector{Float64}(undef, m)
    dp = Vector{Float64}(undef, m)

    # Forward sweep
    cp[1] = 1.0 / 4.0
    dp[1] = scale * (ys[3] - 2.0 * ys[2] + ys[1]) / 4.0
    for i in 2:m
        rhs_i = scale * (ys[i + 2] - 2.0 * ys[i + 1] + ys[i])
        denom = 4.0 - cp[i - 1]
        cp[i] = 1.0 / denom
        dp[i] = (rhs_i - dp[i - 1]) / denom
    end

    # Back substitution: M[2]..M[n-1] stored at Ms[2]..Ms[n-1]
    Ms[m + 1] = dp[m]
    for i in (m - 1):-1:1
        Ms[i + 1] = dp[i] - cp[i] * Ms[i + 2]
    end
    return nothing
end

@inline function _eval_natural_cubic_spline(
    t::Float64,
    t0::Float64,
    h::Float64,
    ys::Vector{Float64},
    Ms::Vector{Float64}
)::Float64
    n = length(ys)
    idx = clamp(floor(Int, (t - t0) / h) + 1, 1, n - 1)
    dx = t - (t0 + (idx - 1) * h)
    @inbounds yi  = ys[idx]
    @inbounds yi1 = ys[idx + 1]
    @inbounds Mi  = Ms[idx]
    @inbounds Mi1 = Ms[idx + 1]
    bi = (yi1 - yi) / h - h * (Mi1 + 2.0 * Mi) / 6.0
    ci = Mi / 2.0
    di = (Mi1 - Mi) / (6.0 * h)
    return yi + dx * (bi + dx * (ci + dx * di))
end

@inline function _interp_vacuum_alt(cache::VacuumPredictedGRAMCache, t::Float64)::Float64
    n = length(cache.vac_alts)
    idx = clamp(floor(Int, (t - cache.t0) / cache.h) + 1, 1, n - 1)
    x = (t - (cache.t0 + (idx - 1) * cache.h)) / cache.h
    @inbounds return (1.0 - x) * cache.vac_alts[idx] + x * cache.vac_alts[idx + 1]
end

@inline function _interp_vacuum_position(cache::VacuumPredictedGRAMCache, t::Float64)::SVector{3, Float64}
    n = length(cache.vac_positions)
    idx = clamp(floor(Int, (t - cache.t0) / cache.h) + 1, 1, n - 1)
    x = (t - (cache.t0 + (idx - 1) * cache.h)) / cache.h
    @inbounds return (1.0 - x) * cache.vac_positions[idx] + x * cache.vac_positions[idx + 1]
end

@inline function _interp_vacuum_wind(cache::VacuumPredictedGRAMCache, t::Float64)::SVector{3, Float64}
    n = length(cache.winds)
    idx = clamp(floor(Int, (t - cache.t0) / cache.h) + 1, 1, n - 1)
    x = (t - (cache.t0 + (idx - 1) * cache.h)) / cache.h
    @inbounds w0 = cache.winds[idx]
    @inbounds w1 = cache.winds[idx + 1]
    return (1.0 - x) * w0 + x * w1
end

function _build_vacuum_gram_cache!(
    cache::VacuumPredictedGRAMCache,
    density_model,
    p,
    pos_ii::SVector{3, Float64},
    vel_ii::SVector{3, Float64},
    t::Float64,
    n_pts::Int,
    horizon_s::Float64
)
    cache.valid = false
    n_pts < 2 && return

    planet = p.args.environment_model.planet
    h = horizon_s / (n_pts - 1)

    resize!(cache.log_rhos, n_pts)
    resize!(cache.Ms_rho, n_pts)
    resize!(cache.Ts, n_pts)
    resize!(cache.Ms_T, n_pts)
    resize!(cache.winds, n_pts)
    resize!(cache.vac_alts, n_pts)
    resize!(cache.vac_positions, n_pts)

    pos = pos_ii
    vel = vel_ii

    for i in 1:n_pts
        t_i = t + (i - 1) * h

        l_pi = _planet_lpi_at(p, t_i)
        pos_pp = SVector{3, Float64}(l_pi * pos)
        alt_i, lat_i, lon_i = rtolatlong(pos_pp, planet)

        @inbounds cache.vac_alts[i] = alt_i
        @inbounds cache.vac_positions[i] = pos

        rho, T, wind = getDensity(density_model, alt_i, lat_i, lon_i, t_i, true, p)
        @inbounds cache.log_rhos[i] = log(max(rho, 1e-40))
        @inbounds cache.Ts[i] = T
        @inbounds cache.winds[i] = wind

        if i < n_pts
            pos, vel = _vacuum_rk4_step(pos, vel, planet, h)
        end
    end

    _natural_cubic_spline_build!(cache.Ms_rho, cache.log_rhos, h)
    _natural_cubic_spline_build!(cache.Ms_T, cache.Ts, h)

    cache.t0 = t
    cache.t1 = t + (n_pts - 1) * h
    cache.h = h
    cache.valid = true
    return nothing
end

function _query_vacuum_gram_cache!(
    cache::VacuumPredictedGRAMCache,
    density_model,
    p,
    pos_ii::SVector{3, Float64},
    vel_ii::SVector{3, Float64},
    alt::Float64,
    t::Float64,
    n_pts::Int,
    horizon_s::Float64,
    deviation_m::Float64
)::Tuple{Float64, Float64, SVector{3, Float64}}
    if cache.valid && t >= cache.t0 && t <= cache.t1
        vac_pos = _interp_vacuum_position(cache, t)
        if norm(pos_ii - vac_pos) <= deviation_m
            log_rho = _eval_natural_cubic_spline(t, cache.t0, cache.h, cache.log_rhos, cache.Ms_rho)
            T_val   = _eval_natural_cubic_spline(t, cache.t0, cache.h, cache.Ts, cache.Ms_T)
            wind    = _interp_vacuum_wind(cache, t)
            return exp(log_rho), T_val, wind
        end
    end

    # Cache miss or trajectory deviation — rebuild from current state.
    _build_vacuum_gram_cache!(cache, density_model, p, pos_ii, vel_ii, t, n_pts, horizon_s)

    if cache.valid && !isempty(cache.log_rhos)
        # At rebuild time t == cache.t0, so the first knot is the current point.
        @inbounds return exp(cache.log_rhos[1]), cache.Ts[1], cache.winds[1]
    end

    # Fallback: direct density query (should not normally reach here).
    planet = p.args.environment_model.planet
    l_pi = _planet_lpi_at(p, t)
    pos_pp = SVector{3, Float64}(l_pi * pos_ii)
    alt_fb, lat_fb, lon_fb = rtolatlong(pos_pp, planet)
    return getDensity(density_model, alt_fb, lat_fb, lon_fb, t, true, p)
end
