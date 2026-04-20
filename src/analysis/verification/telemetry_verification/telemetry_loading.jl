@inline function _to_float_vector(values, context::String)::Vector{Float64}
    out = Vector{Float64}(undef, length(values))
    @inbounds for i in eachindex(values)
        v = values[i]
        v === missing && throw(ArgumentError("Missing value in $context at index $i"))
        out[i] = Float64(v)
    end
    return out
end

@inline function _planet_fixed_frame_name(planet_name::String)::String
    return lowercase(strip(planet_name)) == "earth" ? "ITRF93" : "IAU_" * uppercase(strip(planet_name))
end

function _initial_time_et(initial_time::InitialTime)::Float64
    utc = @sprintf(
        "%04d-%02d-%02dT%02d:%02d:%09.6f",
        Int(initial_time.year),
        Int(initial_time.month),
        Int(initial_time.day),
        Int(initial_time.hour),
        Int(initial_time.minute),
        Float64(initial_time.second)
    )
    return utc2et(utc)
end

function _transform_state(from_frame::String, to_frame::String, et::Float64, r_m::SVector{3, Float64}, v_mps::SVector{3, Float64})
    xform = SMatrix{6, 6, Float64}(sxform(from_frame, to_frame, et))
    state = xform * SVector{6, Float64}(r_m[1], r_m[2], r_m[3], v_mps[1], v_mps[2], v_mps[3])
    return SVector{3, Float64}(state[1], state[2], state[3]), SVector{3, Float64}(state[4], state[5], state[6])
end

@inline function _planet_fixed_to_j2000_state(planet_name::String, et::Float64, r_m::SVector{3, Float64}, v_mps::SVector{3, Float64})
    return _transform_state(_planet_fixed_frame_name(planet_name), "J2000", et, r_m, v_mps)
end

@inline function _j2000_to_planet_fixed_state(planet_name::String, et::Float64, r_m::SVector{3, Float64}, v_mps::SVector{3, Float64})
    return _transform_state("J2000", _planet_fixed_frame_name(planet_name), et, r_m, v_mps)
end

@inline function _require_column(df::DataFrame, candidates::Vector{String}, context::String)::Vector{Float64}
    for col in candidates
        if col in names(df)
            return _to_float_vector(df[!, col], "$context:$col")
        end
    end
    throw(ArgumentError("Missing required column for $context. Tried: $(join(candidates, ", "))"))
end

function _extract_extrema_series(df::DataFrame, planet, altitude_mode::Symbol)
    x = _require_column(df, ["sc1_pos_1", "sc1_position_1"], "position-x")
    y = _require_column(df, ["sc1_pos_2", "sc1_position_2"], "position-y")
    z = _require_column(df, ["sc1_pos_3", "sc1_position_3"], "position-z")
    vx = _require_column(df, ["sc1_vel_1", "sc1_velocity_1"], "velocity-x")
    vy = _require_column(df, ["sc1_vel_2", "sc1_velocity_2"], "velocity-y")
    vz = _require_column(df, ["sc1_vel_3", "sc1_velocity_3"], "velocity-z")
    t = _to_float_vector(df.time, "simulation-time")

    n = length(t)
    n >= 3 || throw(ArgumentError("Not enough rows to extract extrema (need >=3, got $n)."))

    r = sqrt.(x .^ 2 .+ y .^ 2 .+ z .^ 2)
    speed_kmps = sqrt.(vx .^ 2 .+ vy .^ 2 .+ vz .^ 2) .* 1e-3
    alt_km = Vector{Float64}(undef, n)
    if altitude_mode == :vacuum
        @inbounds for i in eachindex(r)
            alt_km[i] = (r[i] - planet.Rp_e) * 1e-3
        end
    elseif altitude_mode == :oblate
        @inbounds for i in eachindex(r)
            pos_i = SVector{3, Float64}(x[i], y[i], z[i])
            vel_i = SVector{3, Float64}(vx[i], vy[i], vz[i])
            pos_p, _ = r_intor_p!(pos_i, vel_i, planet)
            alt_km[i] = rtolatlong(pos_p, planet)[1] * 1e-3
        end
    else
        throw(ArgumentError("Unsupported altitude_mode='$altitude_mode' in _extract_extrema_series"))
    end
    rdot = (x .* vx .+ y .* vy .+ z .* vz) ./ r

    peri_t = Float64[]
    peri_alt = Float64[]
    peri_speed = Float64[]
    apo_t = Float64[]
    apo_alt = Float64[]
    apo_speed = Float64[]
    last_peri_t = -Inf
    last_apo_t = -Inf
    min_sep_s = 500.0

    for i in 1:(n - 1)
        s0 = rdot[i]
        s1 = rdot[i + 1]
        if !(isfinite(s0) && isfinite(s1))
            continue
        end

        denom = abs(s0) + abs(s1)
        w = denom == 0.0 ? 0.5 : abs(s0) / denom
        t_evt = (1.0 - w) * t[i] + w * t[i + 1]
        alt_evt = (1.0 - w) * alt_km[i] + w * alt_km[i + 1]
        speed_evt = (1.0 - w) * speed_kmps[i] + w * speed_kmps[i + 1]

        if s0 < 0.0 && s1 >= 0.0
            if isempty(peri_t) || (t_evt - last_peri_t) >= min_sep_s
                push!(peri_t, t_evt)
                push!(peri_alt, alt_evt)
                push!(peri_speed, speed_evt)
                last_peri_t = t_evt
            elseif alt_evt < peri_alt[end]
                peri_t[end] = t_evt
                peri_alt[end] = alt_evt
                peri_speed[end] = speed_evt
                last_peri_t = t_evt
            end
        elseif s0 > 0.0 && s1 <= 0.0
            if isempty(apo_t) || (t_evt - last_apo_t) >= min_sep_s
                push!(apo_t, t_evt)
                push!(apo_alt, alt_evt)
                push!(apo_speed, speed_evt)
                last_apo_t = t_evt
            elseif alt_evt > apo_alt[end]
                apo_t[end] = t_evt
                apo_alt[end] = alt_evt
                apo_speed[end] = speed_evt
                last_apo_t = t_evt
            end
        end
    end

    if isempty(peri_alt) || isempty(apo_alt)
        # Fallback for coarse/edge trajectories where radial-velocity sign changes
        # are sparse in sampled output; detect local extrema directly in altitude.
        last_peri_t = -Inf
        last_apo_t = -Inf
        for i in 2:(n - 1)
            a0 = alt_km[i - 1]
            a1 = alt_km[i]
            a2 = alt_km[i + 1]
            ti = t[i]
            vi = speed_kmps[i]
            if a1 <= a0 && a1 < a2
                if isempty(peri_t) || (ti - last_peri_t) >= min_sep_s
                    push!(peri_t, ti)
                    push!(peri_alt, a1)
                    push!(peri_speed, vi)
                    last_peri_t = ti
                elseif a1 < peri_alt[end]
                    peri_t[end] = ti
                    peri_alt[end] = a1
                    peri_speed[end] = vi
                    last_peri_t = ti
                end
            elseif a1 >= a0 && a1 > a2
                if isempty(apo_t) || (ti - last_apo_t) >= min_sep_s
                    push!(apo_t, ti)
                    push!(apo_alt, a1)
                    push!(apo_speed, vi)
                    last_apo_t = ti
                elseif a1 > apo_alt[end]
                    apo_t[end] = ti
                    apo_alt[end] = a1
                    apo_speed[end] = vi
                    last_apo_t = ti
                end
            end
        end
    end

    peri_orbit = collect(1.0:length(peri_alt))
    apo_orbit = collect(1.0:length(apo_alt))
    return (
        peri=(orbit=peri_orbit, altitude=peri_alt, time_s=peri_t, speed_kmps=peri_speed),
        apo=(orbit=apo_orbit, altitude=apo_alt, time_s=apo_t, speed_kmps=apo_speed)
    )
end

function _load_telemetry_curve(path::String, max_points::Int)
    df = DataFrame(Arrow.Table(path))
    @assert all(name in names(df) for name in ["orbit", "altitude"])
    sort!(df, :orbit)
    if max_points > 0 && nrow(df) > max_points
        df = first(df, max_points)
    end
    return (
        orbit=_to_float_vector(df.orbit, "telemetry-orbit"),
        altitude=_to_float_vector(df.altitude, "telemetry-altitude")
    )
end

function _load_time_aligned_telemetry(cfg::TimeAlignedScenarioConfig, max_points::Int)
    df = DataFrame(Arrow.Table(cfg.telemetry_path))

    time_s = _require_column(df, [cfg.telemetry_time_col], "telemetry-time")
    altitude_km = _require_column(df, [cfg.telemetry_altitude_col], "telemetry-altitude")
    x_km = _require_column(df, [cfg.telemetry_x_col], "telemetry-x")
    y_km = _require_column(df, [cfg.telemetry_y_col], "telemetry-y")
    z_km = _require_column(df, [cfg.telemetry_z_col], "telemetry-z")
    n_rows = nrow(df)

    has_keplerian_ic = !any(isnothing, (
        cfg.telemetry_sma_col, cfg.telemetry_ecc_col, cfg.telemetry_inc_col,
        cfg.telemetry_aop_col, cfg.telemetry_raan_col, cfg.telemetry_ta_col
    ))
    has_partial_keplerian_ic = any(!isnothing, (
        cfg.telemetry_sma_col, cfg.telemetry_ecc_col, cfg.telemetry_inc_col,
        cfg.telemetry_aop_col, cfg.telemetry_raan_col, cfg.telemetry_ta_col
    ))
    if has_partial_keplerian_ic && !has_keplerian_ic
        throw(ArgumentError(
            "Time-aligned scenario $(cfg.name) must provide either all six Keplerian telemetry columns or none of them."
        ))
    end

    sma_km = has_keplerian_ic ? _require_column(df, [cfg.telemetry_sma_col], "telemetry-sma") : fill(NaN, n_rows)
    ecc = has_keplerian_ic ? _require_column(df, [cfg.telemetry_ecc_col], "telemetry-ecc") : fill(NaN, n_rows)
    inc_deg = has_keplerian_ic ? _require_column(df, [cfg.telemetry_inc_col], "telemetry-inc") : fill(NaN, n_rows)
    aop_deg = has_keplerian_ic ? _require_column(df, [cfg.telemetry_aop_col], "telemetry-aop") : fill(NaN, n_rows)
    raan_deg = has_keplerian_ic ? _require_column(df, [cfg.telemetry_raan_col], "telemetry-raan") : fill(NaN, n_rows)
    ta_deg = has_keplerian_ic ? _require_column(df, [cfg.telemetry_ta_col], "telemetry-ta") : fill(NaN, n_rows)

    # Optional Cartesian IC columns (read before sorting, IC always from row 1 of sorted data)
    has_cartesian_ic = !any(isnothing, (
        cfg.telemetry_x_ic_col, cfg.telemetry_y_ic_col, cfg.telemetry_z_ic_col,
        cfg.telemetry_vx_ic_col, cfg.telemetry_vy_ic_col, cfg.telemetry_vz_ic_col
    ))
    has_partial_cartesian_ic = any(!isnothing, (
        cfg.telemetry_x_ic_col, cfg.telemetry_y_ic_col, cfg.telemetry_z_ic_col,
        cfg.telemetry_vx_ic_col, cfg.telemetry_vy_ic_col, cfg.telemetry_vz_ic_col
    ))
    if has_partial_cartesian_ic && !has_cartesian_ic
        throw(ArgumentError(
            "Time-aligned scenario $(cfg.name) must provide either all six Cartesian IC telemetry columns or none of them."
        ))
    end
    if !(has_cartesian_ic || has_keplerian_ic)
        throw(ArgumentError(
            "Time-aligned scenario $(cfg.name) must provide either Cartesian IC telemetry columns or Keplerian telemetry columns."
        ))
    end

    x_ic_km_col  = has_cartesian_ic ? _require_column(df, [cfg.telemetry_x_ic_col],  "telemetry-x_ic")  : nothing
    y_ic_km_col  = has_cartesian_ic ? _require_column(df, [cfg.telemetry_y_ic_col],  "telemetry-y_ic")  : nothing
    z_ic_km_col  = has_cartesian_ic ? _require_column(df, [cfg.telemetry_z_ic_col],  "telemetry-z_ic")  : nothing
    vx_ic_col    = has_cartesian_ic ? _require_column(df, [cfg.telemetry_vx_ic_col], "telemetry-vx_ic") : nothing
    vy_ic_col    = has_cartesian_ic ? _require_column(df, [cfg.telemetry_vy_ic_col], "telemetry-vy_ic") : nothing
    vz_ic_col    = has_cartesian_ic ? _require_column(df, [cfg.telemetry_vz_ic_col], "telemetry-vz_ic") : nothing

    perm = sortperm(time_s)
    time_s = time_s[perm]
    altitude_km = altitude_km[perm]
    x_km = x_km[perm]
    y_km = y_km[perm]
    z_km = z_km[perm]
    sma_km = sma_km[perm]
    ecc = ecc[perm]
    inc_deg = inc_deg[perm]
    aop_deg = aop_deg[perm]
    raan_deg = raan_deg[perm]
    ta_deg = ta_deg[perm]

    if max_points > 0 && length(time_s) > max_points
        keep = 1:max_points
        time_s = time_s[keep]
        altitude_km = altitude_km[keep]
        x_km = x_km[keep]
        y_km = y_km[keep]
        z_km = z_km[keep]
        sma_km = sma_km[keep]
        ecc = ecc[keep]
        inc_deg = inc_deg[keep]
        aop_deg = aop_deg[keep]
        raan_deg = raan_deg[keep]
        ta_deg = ta_deg[keep]
    end

    length(time_s) >= 2 || throw(ArgumentError("Need at least 2 telemetry samples for $(cfg.name), got $(length(time_s))."))
    t0 = time_s[1]
    time_s = time_s .- t0

    return (
        time_s=time_s,
        altitude_km=altitude_km,
        x_km=x_km,
        y_km=y_km,
        z_km=z_km,
        vx_kmps=_differentiate_series(x_km, time_s),
        vy_kmps=_differentiate_series(y_km, time_s),
        vz_kmps=_differentiate_series(z_km, time_s),
        sma_km=has_keplerian_ic ? sma_km[1] : NaN,
        ecc=has_keplerian_ic ? ecc[1] : NaN,
        inc_deg=has_keplerian_ic ? inc_deg[1] : NaN,
        aop_deg=has_keplerian_ic ? aop_deg[1] : NaN,
        raan_deg=has_keplerian_ic ? raan_deg[1] : NaN,
        ta_deg=has_keplerian_ic ? ta_deg[1] : NaN,
        x_ic_km=has_cartesian_ic ? Float64(x_ic_km_col[perm[1]]) : nothing,
        y_ic_km=has_cartesian_ic ? Float64(y_ic_km_col[perm[1]]) : nothing,
        z_ic_km=has_cartesian_ic ? Float64(z_ic_km_col[perm[1]]) : nothing,
        vx_ic_kmps=has_cartesian_ic ? Float64(vx_ic_col[perm[1]]) : nothing,
        vy_ic_kmps=has_cartesian_ic ? Float64(vy_ic_col[perm[1]]) : nothing,
        vz_ic_kmps=has_cartesian_ic ? Float64(vz_ic_col[perm[1]]) : nothing
    )
end

function _differentiate_series(values::Vector{Float64}, time_s::Vector{Float64})::Vector{Float64}
    n = length(values)
    n == length(time_s) || throw(ArgumentError("values/time length mismatch in _differentiate_series: $n vs $(length(time_s))"))
    n >= 2 || throw(ArgumentError("Need at least 2 samples in _differentiate_series (got $n)."))
    dv = Vector{Float64}(undef, n)
    @inbounds begin
        dt0 = max(time_s[2] - time_s[1], eps(Float64))
        dv[1] = (values[2] - values[1]) / dt0
        for i in 2:(n - 1)
            dt = max(time_s[i + 1] - time_s[i - 1], eps(Float64))
            dv[i] = (values[i + 1] - values[i - 1]) / dt
        end
        dtn = max(time_s[n] - time_s[n - 1], eps(Float64))
        dv[n] = (values[n] - values[n - 1]) / dtn
    end
    return dv
end

function _extract_extrema_from_time_aligned_telemetry(
    time_s::Vector{Float64},
    altitude_km::Vector{Float64},
    min_sep_s::Float64;
    speed_kmps::Union{Nothing, Vector{Float64}}=nothing
)
    n = length(time_s)
    n == length(altitude_km) || throw(ArgumentError("time/altitude length mismatch: $n vs $(length(altitude_km))"))
    n >= 3 || throw(ArgumentError("Need at least 3 telemetry samples to derive peri/apo extrema (got $n)."))
    speed_series = if speed_kmps === nothing
        fill(NaN, n)
    else
        length(speed_kmps) == n || throw(ArgumentError(
            "time/speed length mismatch in _extract_extrema_from_time_aligned_telemetry: $n vs $(length(speed_kmps))"
        ))
        speed_kmps
    end

    peri_t = Float64[]
    peri_alt = Float64[]
    peri_speed = Float64[]
    apo_t = Float64[]
    apo_alt = Float64[]
    apo_speed = Float64[]
    last_peri_t = -Inf
    last_apo_t = -Inf

    @inbounds for i in 2:(n - 1)
        a0 = altitude_km[i - 1]
        a1 = altitude_km[i]
        a2 = altitude_km[i + 1]
        ti = time_s[i]
        vi = speed_series[i]
        if a1 <= a0 && a1 < a2
            if isempty(peri_t) || (ti - last_peri_t) >= min_sep_s
                push!(peri_t, ti)
                push!(peri_alt, a1)
                push!(peri_speed, vi)
                last_peri_t = ti
            elseif a1 < peri_alt[end]
                peri_t[end] = ti
                peri_alt[end] = a1
                peri_speed[end] = vi
                last_peri_t = ti
            end
        elseif a1 >= a0 && a1 > a2
            if isempty(apo_t) || (ti - last_apo_t) >= min_sep_s
                push!(apo_t, ti)
                push!(apo_alt, a1)
                push!(apo_speed, vi)
                last_apo_t = ti
            elseif a1 > apo_alt[end]
                apo_t[end] = ti
                apo_alt[end] = a1
                apo_speed[end] = vi
                last_apo_t = ti
            end
        end
    end

    isempty(peri_alt) && throw(ArgumentError("No periapsis extrema derived from telemetry altitude history."))
    isempty(apo_alt) && throw(ArgumentError("No apoapsis extrema derived from telemetry altitude history."))

    return (
        peri=(orbit=collect(1.0:length(peri_alt)), altitude=peri_alt, time_s=peri_t, speed_kmps=peri_speed),
        apo=(orbit=collect(1.0:length(apo_alt)), altitude=apo_alt, time_s=apo_t, speed_kmps=apo_speed)
    )
end
