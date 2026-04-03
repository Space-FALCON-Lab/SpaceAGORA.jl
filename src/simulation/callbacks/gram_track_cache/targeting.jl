@inline _gram_track_cache_periapsis_split_enabled() = _parse_bool_env("SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT", true)

@inline function _gram_track_cache_max_npos()::Int
    raw = strip(get(ENV, "SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS", "512"))
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS must be an integer value, got '$raw'"))
    end
    return max(2, parsed)
end

@inline _angle_delta_rad(a::Float64, b::Float64) = atan(sin(b - a), cos(b - a))

@inline function _gram_expected_track_length_m(
    planet,
    alt0::Float64,
    lat0::Float64,
    lon0::Float64,
    alt1::Float64,
    lat1::Float64,
    lon1::Float64
)::Tuple{Float64, Float64}
    radius_m = max(1.0, planet.Rp_m + 0.5 * (alt0 + alt1))
    cos_angle = clamp(
        sin(lat0) * sin(lat1) + cos(lat0) * cos(lat1) * cos(_angle_delta_rad(lon0, lon1)),
        -1.0,
        1.0
    )
    central_angle = acos(cos_angle)
    horizontal_m = radius_m * central_angle
    vertical_m = alt1 - alt0
    length_m = hypot(horizontal_m, vertical_m)
    return max(length_m, 1.0), radius_m
end

@inline function _gram_track_cache_target_spacing_m(
    alt_tol_m::Float64,
    ang_tol_rad::Float64,
    radius_m::Float64
)::Float64
    # Keep interpolation spacing tighter than cache acceptance tolerances.
    ang_tol_m = max(1.0, radius_m * max(ang_tol_rad, 1e-9))
    tol_scale_m = min(max(1.0, alt_tol_m), ang_tol_m)
    return max(1.0, 0.5 * tol_scale_m)
end

@inline function _solve_kepler_elliptic(M::Float64, e::Float64)::Float64
    M2π = mod(M, 2pi)
    E = e < 0.8 ? M2π : pi
    @inbounds for _ in 1:20
        f = E - e * sin(E) - M2π
        fp = 1.0 - e * cos(E)
        abs(fp) < 1e-14 && break
        dE = f / fp
        E -= dE
        abs(dE) <= 1e-12 && break
    end
    return E
end

@inline function _true_to_eccentric_anomaly(ν::Float64, e::Float64)::Float64
    E = atan(sqrt(max(0.0, 1.0 - e^2)) * sin(ν), e + cos(ν))
    return mod(E, 2pi)
end

@inline function _eccentric_to_true_anomaly(E::Float64, e::Float64)::Float64
    ν = atan(sqrt(max(0.0, 1.0 - e^2)) * sin(E), cos(E) - e)
    return mod(ν, 2pi)
end

@inline function _gram_j2_rates(a::Float64, e::Float64, i::Float64, n::Float64, planet)::Tuple{Float64, Float64}
    if !isfinite(planet.J2) || planet.J2 == 0.0
        return 0.0, 0.0
    end
    p = a * (1.0 - e^2)
    if !isfinite(p) || p <= 0.0
        return 0.0, 0.0
    end
    scale = planet.J2 * (planet.Rp_e / p)^2
    Ωdot = -1.5 * n * scale * cos(i)
    ωdot = 0.75 * n * scale * (5.0 * cos(i)^2 - 1.0)
    return Ωdot, ωdot
end

function _gram_kepler_target(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet,
    dt::Float64;
    include_j2::Bool = true
)::Union{Nothing, Tuple{Float64, Float64, Float64}}
    try
        # State is in J2000; convert to PCI (body-equatorial) for rvtoorbitalelement
        pos_pci = planet.J2000_to_pci * pos
        vel_pci = planet.J2000_to_pci * vel
        oe = rvtoorbitalelement(pos_pci, vel_pci, planet)
        a, e, i, Ω, ω, ν = Float64(oe[1]), Float64(oe[2]), Float64(oe[3]), Float64(oe[4]), Float64(oe[5]), Float64(oe[6])
        if !isfinite(dt) || dt <= 1e-6 || !isfinite(a) || !isfinite(e) || !isfinite(ν) || a <= 0.0 || e < 0.0 || e >= 1.0
            return nothing
        end
        n = sqrt(planet.μ / a^3)
        if !isfinite(n) || n <= 0.0
            return nothing
        end
        E0 = _true_to_eccentric_anomaly(ν, e)
        M0 = E0 - e * sin(E0)
        E1 = _solve_kepler_elliptic(M0 + n * dt, e)
        ν1 = _eccentric_to_true_anomaly(E1, e)

        Ω1 = Ω
        ω1 = ω
        if include_j2
            Ωdot, ωdot = _gram_j2_rates(a, e, i, n, planet)
            Ω1 += Ωdot * dt
            ω1 += ωdot * dt
        end

        oe_target = SVector{7, Float64}(a, e, i, Ω1, ω1, ν1, 0.0)
        # orbitalelemtorv gives PCI; convert to J2000 before r_intor_p! (which expects J2000)
        r_pci_vec, v_pci_vec = orbitalelemtorv(oe_target, planet)
        r_target_i = planet.J2000_to_pci' * SVector{3, Float64}(Float64.(r_pci_vec))
        v_target_i = planet.J2000_to_pci' * SVector{3, Float64}(Float64.(v_pci_vec))
        r_target_p, _ = r_intor_p!(r_target_i, v_target_i, planet)
        alt_target, lat_target, lon_target = rtolatlong(r_target_p, planet)
        return alt_target, lat_target, lon_target
    catch
        return nothing
    end
end

function _gram_periapsis_target(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet;
    include_j2::Bool = true
)::Union{Nothing, Tuple{Float64, Float64, Float64, Float64}}
    try
        pos_pci = planet.J2000_to_pci * pos
        vel_pci = planet.J2000_to_pci * vel
        oe = rvtoorbitalelement(pos_pci, vel_pci, planet)
        a, e, ν = Float64(oe[1]), Float64(oe[2]), Float64(oe[6])
        if !isfinite(a) || !isfinite(e) || !isfinite(ν) || a <= 0.0 || e < 0.0 || e >= 1.0
            return nothing
        end

        n = sqrt(planet.μ / a^3)
        if !isfinite(n) || n <= 0.0
            return nothing
        end
        E = _true_to_eccentric_anomaly(ν, e)
        M = mod(E - e * sin(E), 2pi)
        dt_peri = (2pi - M) / n
        if !isfinite(dt_peri) || dt_peri <= 1e-6
            return nothing
        end
        target = _gram_kepler_target(pos, vel, planet, dt_peri; include_j2=include_j2)
        target === nothing && return nothing
        alt_peri, lat_peri, lon_peri = target
        return dt_peri, alt_peri, lat_peri, lon_peri
    catch
        return nothing
    end
end

@inline _gram_descending_periapsis_target(pos::SVector{3, Float64}, vel::SVector{3, Float64}, planet; include_j2::Bool = true) =
    _gram_periapsis_target(pos, vel, planet; include_j2=include_j2)

function _gram_orbit_period_target(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet;
    include_j2::Bool = true
)::Union{Nothing, Tuple{Float64, Float64, Float64, Float64}}
    try
        pos_pci = planet.J2000_to_pci * pos
        vel_pci = planet.J2000_to_pci * vel
        oe = rvtoorbitalelement(pos_pci, vel_pci, planet)
        a, e = Float64(oe[1]), Float64(oe[2])
        if !isfinite(a) || !isfinite(e) || a <= 0.0 || e < 0.0 || e >= 1.0
            return nothing
        end
        n = sqrt(planet.μ / a^3)
        if !isfinite(n) || n <= 0.0
            return nothing
        end
        dt_orbit = 2pi / n
        if !isfinite(dt_orbit) || dt_orbit <= 1e-6
            return nothing
        end
        target = _gram_kepler_target(pos, vel, planet, dt_orbit; include_j2=include_j2)
        target === nothing && return nothing
        alt_end, lat_end, lon_end = target
        return dt_orbit, alt_end, lat_end, lon_end
    catch
        return nothing
    end
end

@inline function _gram_linear_target(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet,
    dt::Float64
)::Tuple{Float64, Float64, Float64}
    pos_target = pos + vel * dt
    rp_target, _ = r_intor_p!(pos_target, vel, planet)
    target_alt, target_lat, target_lon = rtolatlong(rp_target, planet)
    return target_alt, target_lat, target_lon
end

@inline function _gram_entry_reference_area_m2(p, sat_index::Int)::Float64
    try
        spacecraft = p.args.dynamics_model.spacecraft[sat_index]
        area = 0.0
        @inbounds for link in spacecraft.links
            a = Float64(link.ref_area)
            if isfinite(a) && a > 0.0
                area += a
            end
        end
        if area > 0.0
            return area
        end
        a0 = Float64(spacecraft.root.ref_area)
        return (isfinite(a0) && a0 > 0.0) ? a0 : 1.0
    catch
        return 1.0
    end
end

@inline function _gram_entry_mass_kg(p, sat_index::Int, current_mass_kg::Float64)::Float64
    if isfinite(current_mass_kg) && current_mass_kg > 0.0
        return current_mass_kg
    end
    try
        spacecraft = p.args.dynamics_model.spacecraft[sat_index]
        m = Float64(spacecraft.dry_mass + spacecraft.prop_mass)
        return (isfinite(m) && m > 0.0) ? m : 100.0
    catch
        return 100.0
    end
end

@inline function _gram_entry_reference_density(planet, h::Float64)::Float64
    H = max(1.0, planet.H)
    ρ = planet.ρ_ref * exp((planet.h_ref - h) / H)
    return (isfinite(ρ) && ρ > 0.0) ? ρ : 0.0
end

function _gram_entry_target_allen_eggers(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet,
    dt::Float64,
    mass_kg::Float64,
    ref_area_m2::Float64;
    cd::Float64 = _gram_entry_target_cd()
)::Union{Nothing, Tuple{Float64, Float64, Float64}}
    if !isfinite(dt) || dt <= 1e-6 || !isfinite(mass_kg) || mass_kg <= 0.0 || !isfinite(ref_area_m2) || ref_area_m2 <= 0.0
        return nothing
    end
    if !isfinite(planet.ρ_ref) || planet.ρ_ref <= 0.0
        return nothing
    end

    rp, vp = r_intor_p!(pos, vel, planet)
    alt0, lat0, lon0 = rtolatlong(rp, planet)
    uD, uN, uE = latlongtoNED(SVector{3, Float64}(alt0, lat0, lon0))

    vN0 = dot(vp, uN)
    vE0 = dot(vp, uE)
    vU0 = -dot(vp, uD)
    vh0 = hypot(vN0, vE0)
    v = norm(vp)
    if !isfinite(v) || v <= 1e-3
        return nothing
    end

    γ = atan(vU0, max(vh0, 1e-9))
    χ = atan(vE0, vN0)
    h = alt0
    lat = lat0
    lon = lon0
    β = mass_kg / max(1e-6, cd * ref_area_m2)

    n_steps = Int(ceil(dt / _gram_entry_target_dt()))
    n_steps = clamp(n_steps, 2, _gram_entry_target_max_steps())
    Δτ = dt / n_steps

    @inbounds for _ in 1:n_steps
        r = max(1.0, planet.Rp_e + h)
        cosγ = cos(γ)
        sinγ = sin(γ)
        cosχ = cos(χ)
        sinχ = sin(χ)
        coslat = max(1e-6, abs(cos(lat)))
        tanlat = clamp(tan(lat), -100.0, 100.0)

        ρ = _gram_entry_reference_density(planet, h)
        g = planet.μ / (r * r)
        drag_acc = 0.5 * ρ * v * v / β

        v_dot = -drag_acc - g * sinγ
        γ_dot = (v / r - g / max(v, 1e-3)) * cosγ
        χ_dot = v * cosγ * sinχ * tanlat / r
        h_dot = v * sinγ
        lat_dot = v * cosγ * cosχ / r
        lon_dot = v * cosγ * sinχ / (r * coslat) - planet.ω[3]

        v = max(1e-3, v + v_dot * Δτ)
        γ = clamp(γ + γ_dot * Δτ, -deg2rad(89.0), deg2rad(89.0))
        χ = atan(sin(χ + χ_dot * Δτ), cos(χ + χ_dot * Δτ))
        h += h_dot * Δτ
        lat = clamp(lat + lat_dot * Δτ, -deg2rad(89.9), deg2rad(89.9))
        lon = atan(sin(lon + lon_dot * Δτ), cos(lon + lon_dot * Δτ))
    end

    return h, lat, lon
end
