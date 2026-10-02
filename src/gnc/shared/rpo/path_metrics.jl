# RPO path metrics shared by HYPR, comparison objectives and reference checks.
# Mathematical inputs are explicit; HYPR owns configuration-based forwarding.
"""Compute path-length and distance references used to normalize RPO cost terms."""
function rpo_path_cost_normalization_refs(points;
    cost_ref_distance_m, sample_ds_m, tf_s, mass_kg, isp_s, g0_mps2)
    pts = Matrix{Float64}(points)
    dx = pts[1, end] - pts[1, 1]
    dy = pts[2, end] - pts[2, 1]
    dz = pts[3, end] - pts[3, 1]
    straight_len = sqrt(dx * dx + dy * dy + dz * dz)
    len_ref = cost_ref_distance_m > 0.0 ? cost_ref_distance_m : straight_len
    len_ref = max(len_ref, sample_ds_m, 1.0e-6)
    v_ref = len_ref / max(tf_s, 1.0e-6)
    fuel_ref = mass_kg * v_ref / max(isp_s * g0_mps2, 1.0e-9)
    return (straight_len=straight_len, len_ref=len_ref, fuel_ref=max(fuel_ref, 1.0e-12))
end

"""Approximate fuel demand from sampled path increments."""
function rpo_fuel_proxy_from_samples(samples; tf_s, mass_kg, isp_s, g0_mps2)
    pts = Matrix{Float64}(samples)
    size(pts, 2) < 3 && return 0.0
    dt = max(tf_s / max(size(pts, 2) - 1, 1), 1.0e-6)
    fuel = 0.0
    @inbounds for j in 1:(size(pts, 2) - 2)
        ax = (pts[1, j + 2] - 2.0 * pts[1, j + 1] + pts[1, j]) / (dt * dt)
        ay = (pts[2, j + 2] - 2.0 * pts[2, j + 1] + pts[2, j]) / (dt * dt)
        az = (pts[3, j + 2] - 2.0 * pts[3, j + 1] + pts[3, j]) / (dt * dt)
        fuel += mass_kg * sqrt(ax * ax + ay * ay + az * az) / max(isp_s * g0_mps2, 1.0e-9) * dt
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
(3 x M positions r_1..r_M in RTN) with zero relative velocity at departure
and arrival: Δv_eq = Σ_k ||u_k|| Δt with u_k from `rpo_hcw_feedforward_accel`
over the positions held at rest before r_1 and after r_M, so the first term
is the departure from rest, v_1/Δt, and the last the braking, -v_{M-1}/Δt,
each with its HCW terms. J_fuel = m Δv_eq / (Isp g0). Returns
`(J_fuel, delta_v_eq_mps)`.
"""
function rpo_hcw_fuel_proxy(positions, dt::Real, mean_motion::Real, mass_kg::Real, isp_s::Real, g0_mps2::Real)
    pts = positions
    m = size(pts, 2)
    step = Float64(dt)
    n = Float64(mean_motion)
    point(k) = SVector{3, Float64}(pts[1, clamp(k, 1, m)], pts[2, clamp(k, 1, m)], pts[3, clamp(k, 1, m)])
    dv = 0.0
    if m >= 2
        @inbounds for k in 0:(m - 1)
            dv += norm(rpo_hcw_feedforward_accel(point(k), point(k + 1), point(k + 2), step, n)) * step
        end
    end
    return (J_fuel=Float64(mass_kg) * dv / (Float64(isp_s) * Float64(g0_mps2)), delta_v_eq_mps=dv)
end

"""
HCW fuel proxy of an acceleration-limited profile, streamed at step `dt`
without storing the reference: positions r_0..r_K from
`rpo_retimed_reference_from_profile`, held at rest before r_0 and after r_K,
so the first term is the departure from rest and the last the arrival, as in
`rpo_hcw_fuel_proxy` on the stored reference.
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
    # Departure from rest: the reference holds r_0 before it starts.
    dv = K > 0 ? norm(rpo_hcw_feedforward_accel(r0, r0, r1, step, n)) * step : 0.0
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
    rpo_reference_accel_demand(ref, mean_motion)

Acceleration an acceleration-limited reference (`rpo_retimed_reference`) asks
of the actuators at each sample: the tangential acceleration along the curve,
the centripetal term v²κ toward the centre of curvature, and the HCW terms,
u = r̈ - [3n² x + 2n ẏ, -2n ẋ, -n² z]. Returns per-sample magnitudes
`tangential_mps2`, `centripetal_mps2` and `hcw_mps2`, and `u_rtn` (3 x K).
A polyline has no curvature between its vertices, so its centripetal term is
zero here.
"""
function rpo_reference_accel_demand(ref, mean_motion::Real)
    profile = ref.profile
    K = length(ref.t_s)
    n = Float64(mean_motion)
    u_rtn = zeros(3, K)
    tangential = zeros(K)
    centripetal = zeros(K)
    hcw = zeros(K)
    ns = length(profile.s)
    j = 1
    @inbounds for k in 1:K
        speed = ref.speed_mps[k]
        a_t = ref.accel_tangential_mps2[k]
        x, y, z = ref.r_rtn[1, k], ref.r_rtn[2, k], ref.r_rtn[3, k]
        vx, vy, vz = ref.v_rtn[1, k], ref.v_rtn[2, k], ref.v_rtn[3, k]
        tangent = SVector{3, Float64}(0.0, 0.0, 0.0)
        κvec = SVector{3, Float64}(0.0, 0.0, 0.0)
        if ns >= 2
            sq = ref.s_m[k]
            while j < ns - 1 && profile.s[j + 1] < sq
                j += 1
            end
            if profile.bezier
                uq = _rpo_profile_param(profile, j, sq)
                d1 = _rpo_curve_d1(profile.curve, uq)
                d1n = norm(d1)
                if d1n > 1.0e-12
                    tangent = d1 / d1n
                    if size(profile.curve.d2, 2) > 0
                        d2 = _rpo_curve_d2(profile.curve, uq)
                        κvec = (d2 - dot(d2, tangent) * tangent) / d1n^2
                    end
                end
            else
                _, tangent = _rpo_profile_point_tangent(profile, j, sq)
            end
        end
        acc = a_t * tangent + speed^2 * κvec
        g = SVector{3, Float64}(3.0 * n * n * x + 2.0 * n * vy, -2.0 * n * vx, -n * n * z)
        u = acc - g
        u_rtn[1, k] = u[1]
        u_rtn[2, k] = u[2]
        u_rtn[3, k] = u[3]
        tangential[k] = abs(a_t)
        centripetal[k] = speed^2 * norm(κvec)
        hcw[k] = norm(g)
    end
    return (tangential_mps2=tangential, centripetal_mps2=centripetal, hcw_mps2=hcw, u_rtn=u_rtn)
end
