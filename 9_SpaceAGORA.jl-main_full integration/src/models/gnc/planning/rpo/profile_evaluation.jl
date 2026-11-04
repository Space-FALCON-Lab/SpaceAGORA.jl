# Evaluate a supplied speed profile without choosing a planner objective.
# Definitions retain their GuidanceHooks identity.

"""
    rpo_accel_limited_speeds(s, v_point, a_max; v_start=0.0, v_end=0.0)

Limit pointwise speeds `v_point` at arc lengths `s` so the tangential
acceleration never exceeds `a_max`: the profile starts at `v_start`, ends at
`v_end`, and a forward and a backward pass bound the change of v² between
neighbouring samples by 2 a_max Δs.
"""
function rpo_accel_limited_speeds(s, v_point, a_max::Real; v_start::Real=0.0, v_end::Real=0.0)
    n = length(s)
    v = Vector{Float64}(v_point)
    n == 0 && return v
    a = Float64(a_max)
    v[1] = Float64(v_start)
    n == 1 && return v
    v[n] = Float64(v_end)
    @inbounds for j in 1:(n - 1)
        ds = max(s[j + 1] - s[j], 0.0)
        v[j + 1] = min(v[j + 1], sqrt(v[j]^2 + 2.0 * a * ds))
    end
    @inbounds for j in (n - 1):-1:1
        ds = max(s[j + 1] - s[j], 0.0)
        v[j] = min(v[j], sqrt(v[j + 1]^2 + 2.0 * a * ds))
    end
    return v
end

"""Curve parameter at arc length `sq` inside profile interval `j` (cubic Hermite with du/ds = 1/|r'(u)| at the samples)."""
@inline function _rpo_profile_param(profile, j::Int, sq::Float64)
    s = profile.s
    u = profile.params
    h = s[j + 1] - s[j]
    h <= 0.0 && return u[j]
    τ = clamp((sq - s[j]) / h, 0.0, 1.0)
    du = u[j + 1] - u[j]
    rp0 = profile.rprime[j]
    rp1 = profile.rprime[j + 1]
    m0 = rp0 > 1.0e-12 ? h / rp0 : du
    m1 = rp1 > 1.0e-12 ? h / rp1 : du
    τ2 = τ * τ
    τ3 = τ2 * τ
    uq = (2.0 * τ3 - 3.0 * τ2 + 1.0) * u[j] + (τ3 - 2.0 * τ2 + τ) * m0 + (-2.0 * τ3 + 3.0 * τ2) * u[j + 1] + (τ3 - τ2) * m1
    return clamp(uq, u[j], u[j + 1])
end

"""Position and unit tangent at arc length `sq` inside profile interval `j`."""
@inline function _rpo_profile_point_tangent(profile, j::Int, sq::Float64)
    if profile.bezier
        uq = _rpo_profile_param(profile, j, sq)
        p = _rpo_curve_point(profile.curve, uq)
        d = _rpo_curve_d1(profile.curve, uq)
        dn = norm(d)
        return p, dn > 1.0e-12 ? d / dn : SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    pts = profile.samples
    s = profile.s
    a = @inbounds SVector{3, Float64}(pts[1, j], pts[2, j], pts[3, j])
    b = @inbounds SVector{3, Float64}(pts[1, j + 1], pts[2, j + 1], pts[3, j + 1])
    h = s[j + 1] - s[j]
    h <= 0.0 && return a, SVector{3, Float64}(0.0, 0.0, 0.0)
    α = clamp((sq - s[j]) / h, 0.0, 1.0)
    return (1.0 - α) * a + α * b, (b - a) / h
end

"""Profile interval, arc length, speed and tangential acceleration at time `tq`, scanning forward from `j`."""
@inline function _rpo_profile_state_at_time(profile, tq::Float64, j::Int)
    t = profile.t
    n = length(t)
    if tq >= t[n]
        return n - 1, profile.s[n], profile.v[n], 0.0
    end
    @inbounds while j < n - 1 && t[j + 1] <= tq
        j += 1
    end
    τ = tq - t[j]
    a = profile.a_seg[j]
    sq = min(profile.s[j] + profile.v[j] * τ + 0.5 * a * τ * τ, profile.s[j + 1])
    return j, sq, max(0.0, profile.v[j] + a * τ), a
end

"""Number of fixed steps that brings the reference to the goal: t_K = K Δt ≥ duration."""
@inline function _rpo_profile_step_count(duration_s::Float64, dt::Float64)
    duration_s <= 0.0 && return 0
    return max(1, ceil(Int, duration_s / dt - 1.0e-9))
end

"""
    rpo_retimed_reference_from_profile(profile, dt)

Sample an acceleration-limited profile at the fixed step `dt`: t_k = k dt for
k = 0..K with t_K the first step at or after arrival. Positions lie on the
curve, velocities are speed times the unit tangent, and the last sample is
the goal at rest.
"""
function rpo_retimed_reference_from_profile(profile, dt::Real)
    step = Float64(dt)
    n = length(profile.s)
    K = _rpo_profile_step_count(profile.duration_s, step)
    t_s = collect(0:K) .* step
    r = zeros(3, K + 1)
    v = zeros(3, K + 1)
    s_ref = zeros(K + 1)
    speed = zeros(K + 1)
    accel = zeros(K + 1)
    j = 1
    @inbounds for k in 0:K
        if n == 1
            r[:, k + 1] .= profile.samples[:, 1]
            continue
        end
        j, sq, vq, a = _rpo_profile_state_at_time(profile, t_s[k + 1], j)
        p, tangent = _rpo_profile_point_tangent(profile, j, sq)
        r[1, k + 1] = p[1]
        r[2, k + 1] = p[2]
        r[3, k + 1] = p[3]
        v[1, k + 1] = vq * tangent[1]
        v[2, k + 1] = vq * tangent[2]
        v[3, k + 1] = vq * tangent[3]
        s_ref[k + 1] = sq
        speed[k + 1] = vq
        accel[k + 1] = a
    end
    r[:, end] .= profile.samples[:, end]
    v[:, end] .= 0.0
    s_ref[end] = profile.s[end]
    speed[end] = 0.0
    accel[end] = 0.0
    return (t_s=t_s, r_rtn=r, v_rtn=v, s_m=s_ref, speed_mps=speed, accel_tangential_mps2=accel, profile=profile)
end
