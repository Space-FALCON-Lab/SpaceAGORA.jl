"""Compute path-length and distance references used to normalize RPO cost terms."""
function rpo_path_cost_normalization_refs(points, cfg::RPOPSOConfig)
    pts = Matrix{Float64}(points)
    dx = pts[1, end] - pts[1, 1]
    dy = pts[2, end] - pts[2, 1]
    dz = pts[3, end] - pts[3, 1]
    straight_len = sqrt(dx * dx + dy * dy + dz * dz)
    len_ref = cfg.cost_ref_distance_m > 0.0 ? cfg.cost_ref_distance_m : straight_len
    len_ref = max(len_ref, cfg.sample_ds_m, 1.0e-6)
    v_ref = len_ref / max(cfg.tf_s, 1.0e-6)
    fuel_ref = cfg.mass_kg * v_ref / max(cfg.isp_s * cfg.g0_mps2, 1.0e-9)
    return (straight_len=straight_len, len_ref=len_ref, fuel_ref=max(fuel_ref, 1.0e-12))
end

"""Approximate fuel demand from sampled path increments."""
function rpo_fuel_proxy_from_samples(samples, cfg::RPOPSOConfig)
    pts = Matrix{Float64}(samples)
    size(pts, 2) < 3 && return 0.0
    dt = max(cfg.tf_s / max(size(pts, 2) - 1, 1), 1.0e-6)
    fuel = 0.0
    @inbounds for j in 1:(size(pts, 2) - 2)
        ax = (pts[1, j + 2] - 2.0 * pts[1, j + 1] + pts[1, j]) / (dt * dt)
        ay = (pts[2, j + 2] - 2.0 * pts[2, j + 1] + pts[2, j]) / (dt * dt)
        az = (pts[3, j + 2] - 2.0 * pts[3, j + 1] + pts[3, j]) / (dt * dt)
        fuel += cfg.mass_kg * sqrt(ax * ax + ay * ay + az * az) / max(cfg.isp_s * cfg.g0_mps2, 1.0e-9) * dt
    end
    return fuel
end

"""
Clearance at which the obstacle sigmoid is centred. Eq. 6 of the HyPR
manuscript (`threshold_mode = :manuscript`) uses d_safe + τ_tol; the legacy
objective uses d_safe - τ_tol.
"""
@inline function rpo_obstacle_sigmoid_threshold(safe_distance_m::Real, tol_m::Real, threshold_mode::Symbol=:legacy)
    threshold_mode === :manuscript && return Float64(safe_distance_m) + Float64(tol_m)
    return Float64(safe_distance_m) - Float64(tol_m)
end

"""
HCW feedforward acceleration of Sec. III.B for three consecutive reference
positions at fixed step `dt`, with forward-difference velocities
v_k = (r_{k+1} - r_k)/Δt and v_{k+1} = (r_{k+2} - r_{k+1})/Δt:
u = [(v_{x,k+1} - v_{x,k})/Δt - 3n² x_k - 2n v_{y,k};
     (v_{y,k+1} - v_{y,k})/Δt + 2n v_{x,k};
     (v_{z,k+1} - v_{z,k})/Δt + n² z_k].
"""
@inline function rpo_hcw_feedforward_accel(r0::SVector{3, Float64}, r1::SVector{3, Float64}, r2::SVector{3, Float64}, dt::Float64, n::Float64)
    v0 = (r1 - r0) / dt
    v1 = (r2 - r1) / dt
    return SVector{3, Float64}(
        (v1[1] - v0[1]) / dt - 3.0 * n * n * r0[1] - 2.0 * n * v0[2],
        (v1[2] - v0[2]) / dt + 2.0 * n * v0[1],
        (v1[3] - v0[3]) / dt + n * n * r0[3],
    )
end

"""
    rpo_hcw_fuel_proxy(positions, dt, mean_motion, mass_kg, isp_s, g0_mps2)

Fuel proxy of Sec. III.B on a reference sampled at the fixed step `dt`
(3 x M positions in RTN): Δv_eq = Σ_k ||u_k|| Δt over every k with r_{k+2}
available, u_k from `rpo_hcw_feedforward_accel`, and
J_fuel = m Δv_eq / (Isp g0). Returns `(J_fuel, delta_v_eq_mps)`.
"""
function rpo_hcw_fuel_proxy(positions, dt::Real, mean_motion::Real, mass_kg::Real, isp_s::Real, g0_mps2::Real)
    pts = positions
    m = size(pts, 2)
    step = Float64(dt)
    n = Float64(mean_motion)
    dv = 0.0
    @inbounds for k in 1:(m - 2)
        r0 = SVector{3, Float64}(pts[1, k], pts[2, k], pts[3, k])
        r1 = SVector{3, Float64}(pts[1, k + 1], pts[2, k + 1], pts[3, k + 1])
        r2 = SVector{3, Float64}(pts[1, k + 2], pts[2, k + 2], pts[3, k + 2])
        dv += norm(rpo_hcw_feedforward_accel(r0, r1, r2, step, n)) * step
    end
    return (J_fuel=Float64(mass_kg) * dv / (Float64(isp_s) * Float64(g0_mps2)), delta_v_eq_mps=dv)
end

"""Fixed step at which `:manuscript` mode evaluates each candidate's retimed reference."""
rpo_fuel_proxy_dt_s(cfg::RPOPSOConfig) = cfg.fuel_proxy_dt_s > 0.0 ? cfg.fuel_proxy_dt_s : cfg.retime_dt_s

"""
HCW fuel proxy of an acceleration-limited profile, streamed at step `dt`
without storing the reference: positions r_0..r_K from
`rpo_retimed_reference_from_profile`, then r_{K+1} = r_K because the
reference holds the goal at rest, so the last term includes arrival.
"""
function rpo_profile_hcw_fuel_proxy(profile, dt::Real, mean_motion::Real, mass_kg::Real, isp_s::Real, g0_mps2::Real)
    step = Float64(dt)
    n = Float64(mean_motion)
    goal = SVector{3, Float64}(profile.samples[1, end], profile.samples[2, end], profile.samples[3, end])
    if length(profile.s) < 2
        return (J_fuel=0.0, delta_v_eq_mps=0.0, steps=0, duration_s=0.0)
    end
    K = _rpo_profile_step_count(profile.duration_s, step)
    j = 1
    j, sq, _, _ = _rpo_profile_state_at_time(profile, 0.0, j)
    r0 = K == 0 ? goal : first(_rpo_profile_point_tangent(profile, j, sq))
    r1 = goal
    if K > 1
        j, sq, _, _ = _rpo_profile_state_at_time(profile, step, j)
        r1 = first(_rpo_profile_point_tangent(profile, j, sq))
    end
    dv = 0.0
    for k in 0:(K - 1)
        r2 = goal
        if k + 2 < K
            j, sq, _, _ = _rpo_profile_state_at_time(profile, (k + 2) * step, j)
            r2 = first(_rpo_profile_point_tangent(profile, j, sq))
        end
        dv += norm(rpo_hcw_feedforward_accel(r0, r1, r2, step, n)) * step
        r0 = r1
        r1 = r2
    end
    return (
        J_fuel=Float64(mass_kg) * dv / (Float64(isp_s) * Float64(g0_mps2)),
        delta_v_eq_mps=dv,
        steps=K,
        duration_s=profile.duration_s,
    )
end

"""Summarize obstacle clearance, violations, and penalties along sampled RPO path points."""
function rpo_clearance_stats_from_samples(
    samples,
    geometry,
    safe_distance_m::Real;
    cost_cutoff::Real=Inf,
    w_obs::Real=0.0,
    obstacle_sigmoid_k::Real=1.0e6,
    obstacle_sigmoid_tol_m::Real=0.0,
    threshold_mode::Symbol=:legacy,
)
    min_clearance = Inf
    violation_count = 0
    obstacle_score = 0.0
    safe = Float64(safe_distance_m)
    threshold = rpo_obstacle_sigmoid_threshold(safe, obstacle_sigmoid_tol_m, threshold_mode)
    k = Float64(obstacle_sigmoid_k)
    @inbounds for j in 1:size(samples, 2)
        p = SVector{3, Float64}(samples[1, j], samples[2, j], samples[3, j])
        clearance = rpo_clearance_distance_to_station(p, geometry)
        min_clearance = min(min_clearance, clearance)
        beta = rpo_obstacle_sigmoid_penalty(clearance, threshold, k)
        obstacle_score += beta
        if clearance < 0.0
            violation_count += 1
        end
        if w_obs > 0.0 && isfinite(cost_cutoff) && Float64(w_obs) * obstacle_score > Float64(cost_cutoff)
            return (
                min_clearance=min_clearance,
                violation_count=violation_count,
                violation_fraction=violation_count / max(size(samples, 2), 1),
                obstacle_score=obstacle_score,
                cutoff_exceeded=true,
            )
        end
    end
    return (
        min_clearance=min_clearance,
        violation_count=violation_count,
        violation_fraction=violation_count / max(size(samples, 2), 1),
        obstacle_score=obstacle_score,
        cutoff_exceeded=false,
    )
end

"""Evaluate the smooth obstacle penalty for a clearance value."""
@inline function rpo_obstacle_sigmoid_penalty(clearance::Real, threshold::Real, k::Real)
    x = Float64(k) * (Float64(clearance) - Float64(threshold))
    if x >= 0.0
        y = exp(-x)
        return y / (1.0 + y)
    else
        return 1.0 / (1.0 + exp(x))
    end
end

"""
Objective components of the `:manuscript` mode for one candidate: Eq. 5,
J = w_obs J_obs + w_fuel J_fuel, with J_obs from Eq. 6 (sigmoid centred at
d_safe + τ_tol) over the adaptive samples and J_fuel the HCW fuel proxy of
Sec. III.B on the candidate's retimed reference at `rpo_fuel_proxy_dt_s`.
The candidate is retimed by the same retimer as the final reference; with
`retime_accel_limit_enable` the samples used for J_obs are reused for it.
There is no length term: `J_len` and the normalized values are reported only.
A violation is a sample below the Eq. 6 threshold; `keepout_violation_count`
counts samples inside the keep-out surface itself.
"""
function rpo_manuscript_path_cost_components(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    cost_cutoff::Real=Inf,
)
    samples, params, clearances = rpo_sample_path_with_params(
        points,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=rpo_hypr_sampling_density_m(cfg, safe_distance_m),
        curve_type=cfg.curve_type,
    )
    threshold = rpo_obstacle_sigmoid_threshold(safe_distance_m, cfg.obstacle_sigmoid_tol_m, :manuscript)
    k = Float64(cfg.obstacle_sigmoid_k)
    w_obs = Float64(cfg.w_obs)
    cutoff = Float64(cost_cutoff)
    min_clearance = Inf
    violation_count = 0
    keepout_violation_count = 0
    J_obs = 0.0
    cut = false
    @inbounds for j in 1:size(samples, 2)
        c = clearances[j]
        if isnan(c)
            c = rpo_clearance_distance_to_station(SVector{3, Float64}(samples[1, j], samples[2, j], samples[3, j]), geometry)
            clearances[j] = c
        end
        min_clearance = min(min_clearance, c)
        J_obs += rpo_obstacle_sigmoid_penalty(c, threshold, k)
        c < threshold && (violation_count += 1)
        c < 0.0 && (keepout_violation_count += 1)
        if w_obs > 0.0 && isfinite(cutoff) && w_obs * J_obs > cutoff
            cut = true
            break
        end
    end
    refs = rpo_path_cost_normalization_refs(points, cfg)
    J_len = rpo_path_length(samples)
    if cut
        return (
            total=Inf,
            J_len=J_len,
            J_len_norm=J_len / refs.len_ref,
            J_obs=J_obs,
            J_fuel=0.0,
            J_fuel_norm=0.0,
            min_clearance=min_clearance,
            violation_count=violation_count,
            len_ref=refs.len_ref,
            fuel_ref=refs.fuel_ref,
            keepout_violation_count=keepout_violation_count,
            delta_v_eq_mps=0.0,
            reference_duration_s=0.0,
        )
    end
    dt = rpo_fuel_proxy_dt_s(cfg)
    fuel = if cfg.retime_accel_limit_enable
        profile = rpo_retime_profile(
            RPORetimeCurve(points, cfg.curve_type),
            samples,
            params,
            clearances,
            geometry,
            cfg;
            safe_distance_m=safe_distance_m,
            warn=false,
        )
        rpo_profile_hcw_fuel_proxy(profile, dt, cfg.mean_motion_radps, cfg.mass_kg, cfg.isp_s, cfg.g0_mps2)
    else
        step_cfg = dt == cfg.retime_dt_s ? cfg : rpo_pso_config(cfg; retime_dt_s=dt)
        r_ref, _, _ = Logging.with_logger(Logging.NullLogger()) do
            rpo_retime_path(points, geometry, step_cfg; safe_distance_m=safe_distance_m)
        end
        held = hcat(r_ref, r_ref[:, end])
        proxy = rpo_hcw_fuel_proxy(held, dt, cfg.mean_motion_radps, cfg.mass_kg, cfg.isp_s, cfg.g0_mps2)
        (J_fuel=proxy.J_fuel, delta_v_eq_mps=proxy.delta_v_eq_mps, steps=size(r_ref, 2) - 1, duration_s=(size(r_ref, 2) - 1) * dt)
    end
    return (
        total=w_obs * J_obs + Float64(cfg.w_fuel) * fuel.J_fuel,
        J_len=J_len,
        J_len_norm=J_len / refs.len_ref,
        J_obs=J_obs,
        J_fuel=fuel.J_fuel,
        J_fuel_norm=fuel.J_fuel / refs.fuel_ref,
        min_clearance=min_clearance,
        violation_count=violation_count,
        len_ref=refs.len_ref,
        fuel_ref=refs.fuel_ref,
        keepout_violation_count=keepout_violation_count,
        delta_v_eq_mps=fuel.delta_v_eq_mps,
        reference_duration_s=fuel.duration_s,
    )
end

"""Compute normalized RPO objective components for a candidate path."""
function rpo_normalized_path_cost_components(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    cost_cutoff::Real=Inf,
)
    if cfg.hypr_mode === :manuscript
        return rpo_manuscript_path_cost_components(
            points,
            geometry,
            cfg;
            safe_distance_m=safe_distance_m,
            cost_cutoff=cost_cutoff,
        )
    end
    samples = rpo_sample_path(
        points,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=rpo_hypr_sampling_density_m(cfg, safe_distance_m),
        curve_type=cfg.curve_type,
    )
    stats = rpo_clearance_stats_from_samples(
        samples,
        geometry,
        safe_distance_m;
        cost_cutoff=cost_cutoff,
        w_obs=cfg.w_obs,
        obstacle_sigmoid_k=cfg.obstacle_sigmoid_k,
        obstacle_sigmoid_tol_m=cfg.obstacle_sigmoid_tol_m,
    )
    J_obs = stats.obstacle_score
    if stats.cutoff_exceeded
        return (
            total=Inf,
            J_len=0.0,
            J_len_norm=0.0,
            J_obs=J_obs,
            J_fuel=0.0,
            J_fuel_norm=0.0,
            min_clearance=stats.min_clearance,
            violation_count=stats.violation_count,
            len_ref=0.0,
            fuel_ref=0.0,
        )
    end
    refs = rpo_path_cost_normalization_refs(points, cfg)
    J_len = rpo_path_length(samples)
    J_len_norm = J_len / refs.len_ref
    partial_cost = cfg.w_obs * J_obs + cfg.w_len * J_len_norm^2
    if isfinite(cost_cutoff) && partial_cost > Float64(cost_cutoff)
        return (
            total=Inf,
            J_len=J_len,
            J_len_norm=J_len_norm,
            J_obs=J_obs,
            J_fuel=0.0,
            J_fuel_norm=0.0,
            min_clearance=stats.min_clearance,
            violation_count=stats.violation_count,
            len_ref=refs.len_ref,
            fuel_ref=refs.fuel_ref,
        )
    end
    J_fuel = rpo_fuel_proxy_from_samples(samples, cfg)
    J_fuel_norm = J_fuel / refs.fuel_ref
    return (
        total=cfg.w_len * J_len_norm^2 + cfg.w_obs * J_obs + cfg.w_fuel * J_fuel_norm^2,
        J_len=J_len,
        J_len_norm=J_len_norm,
        J_obs=J_obs,
        J_fuel=J_fuel,
        J_fuel_norm=J_fuel_norm,
        min_clearance=stats.min_clearance,
        violation_count=stats.violation_count,
        len_ref=refs.len_ref,
        fuel_ref=refs.fuel_ref,
    )
end

"""Return the scalar RPO path objective, optionally cutting off expensive candidates early."""
function rpo_path_cost(points, geometry, cfg::RPOPSOConfig; safe_distance_m::Real=0.0, cost_cutoff::Real=Inf)
    return rpo_normalized_path_cost_components(points, geometry, cfg; safe_distance_m=safe_distance_m, cost_cutoff=cost_cutoff).total
end
