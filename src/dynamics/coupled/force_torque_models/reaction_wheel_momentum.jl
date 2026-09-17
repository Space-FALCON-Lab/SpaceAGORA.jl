# The reaction torque a set of momentum wheels puts on the body when the wheel
# speeds are PRESCRIBED rather than commanded: the wheel-speed history is given
# as a table (a flight telemetry channel, say), and the effector returns the
# torque that history implies for the body it is mounted in.
#
# For a rigid body carrying wheels whose total momentum in the body frame is
# `H_w(t)`, Euler's equation in the body frame reads
#
#     I ω̇ + ω × (I ω + H_w) + Ḣ_w = τ_external,
#
# so the wheels act on the body as the torque
#
#     τ_wheel = -Ḣ_w - ω × H_w,
#
# with the remaining `ω × I ω` term supplied by the integrator's own rigid-body
# right-hand side. With no external torque this is exactly equivalent to the
# algebraic conservation law `I ω(t) + H_w(t) = constant in inertial space`, and
# it is the form that lets the same statement be integrated alongside every
# other torque the run carries.
#
# `Ḣ_w` is taken analytically from a natural cubic spline through the wheel-speed
# samples, not from differencing them: a flight wheel tachometer is commonly
# quantized to whole revolutions per minute, and differencing a quantized signal
# produces a derivative dominated by the quantization step.
module ReactionWheelMomentum

using StaticArrays
using LinearAlgebra
using ...AbstractTypes: AbstractForceTorqueModel
using ...EffectorSampling: StateSample, EnvironmentSample, EffectorEnvironmentRequirements
import ..DynamicEffectors: wrench, wrench_caching!, environment_requirements

export WheelSpeedSpline, ReactionWheelMomentumModel, ReactionWheelMomentumState
export wheel_speed_spline, wheel_spline_value, wheel_spline_derivative
export wheel_momentum_body, wheel_momentum_rate_body, wheel_reaction_torque
export wheel_speeds_rad_s

"""
    WheelSpeedSpline(t_s, y)

Natural cubic spline through `(t_s, y)`, with `t_s` strictly increasing and not
necessarily uniformly spaced. "Natural" means the second derivative is zero at
both ends; the interior is the standard tridiagonal spline, so the curve is C²
and its derivative is available in closed form.

Queries outside `[t_s[1], t_s[end]]` are clamped to the end knot rather than
extrapolated: a wheel-speed table has no meaning outside the arc it was
measured over, and a cubic extrapolation of one diverges quickly.

Evaluation allocates nothing.
"""
struct WheelSpeedSpline
    t_s::Vector{Float64}
    y::Vector{Float64}
    d2::Vector{Float64}
end

"""
    wheel_speed_spline(t_s, y) -> WheelSpeedSpline

Fit the natural cubic spline. Solves the tridiagonal system for the knot second
derivatives once, at construction.
"""
function wheel_speed_spline(t_s::AbstractVector{<:Real}, y::AbstractVector{<:Real})::WheelSpeedSpline
    x = Vector{Float64}(t_s)
    v = Vector{Float64}(y)
    n = length(x)
    n == length(v) || throw(ArgumentError("wheel_speed_spline: t_s and y must be the same length, got $(n) and $(length(v))."))
    n >= 3 || throw(ArgumentError("wheel_speed_spline: a cubic spline needs at least three samples, got $(n)."))
    all(isfinite, x) && all(isfinite, v) || throw(ArgumentError("wheel_speed_spline: t_s and y must be finite."))
    for k in 2:n
        x[k] > x[k - 1] || throw(ArgumentError("wheel_speed_spline: t_s must be strictly increasing (sample $(k))."))
    end

    h = diff(x)
    lower = zeros(Float64, n)
    diagonal = ones(Float64, n)
    upper = zeros(Float64, n)
    rhs = zeros(Float64, n)
    @inbounds for k in 2:(n - 1)
        lower[k] = h[k - 1] / 6
        diagonal[k] = (h[k - 1] + h[k]) / 3
        upper[k] = h[k] / 6
        rhs[k] = (v[k + 1] - v[k]) / h[k] - (v[k] - v[k - 1]) / h[k - 1]
    end

    # Thomas algorithm; the first and last rows are the natural boundary
    # condition d2 = 0, already in place as a unit diagonal with a zero right
    # hand side.
    c = zeros(Float64, n)
    d = zeros(Float64, n)
    c[1] = upper[1] / diagonal[1]
    d[1] = rhs[1] / diagonal[1]
    @inbounds for k in 2:n
        denom = diagonal[k] - lower[k] * c[k - 1]
        c[k] = upper[k] / denom
        d[k] = (rhs[k] - lower[k] * d[k - 1]) / denom
    end
    d2 = zeros(Float64, n)
    d2[n] = d[n]
    @inbounds for k in (n - 1):-1:1
        d2[k] = d[k] - c[k] * d2[k + 1]
    end
    return WheelSpeedSpline(x, v, d2)
end

@inline function _spline_segment(s::WheelSpeedSpline, t::Float64)::Tuple{Int, Float64}
    n = length(s.t_s)
    tc = clamp(t, s.t_s[1], s.t_s[n])
    i = clamp(searchsortedlast(s.t_s, tc), 1, n - 1)
    return i, tc
end

"""
    wheel_spline_value(spline, t) -> Float64

The spline's value at `t`, clamped to the table's own time span.
"""
@inline function wheel_spline_value(s::WheelSpeedSpline, t::Float64)::Float64
    i, tc = _spline_segment(s, t)
    @inbounds begin
        h = s.t_s[i + 1] - s.t_s[i]
        a = (s.t_s[i + 1] - tc) / h
        b = (tc - s.t_s[i]) / h
        return a * s.y[i] + b * s.y[i + 1] +
            ((a^3 - a) * s.d2[i] + (b^3 - b) * s.d2[i + 1]) * h * h / 6
    end
end

"""
    wheel_spline_derivative(spline, t) -> Float64

The spline's analytic first derivative at `t`. Zero outside the table's time
span, matching the clamped value there.
"""
@inline function wheel_spline_derivative(s::WheelSpeedSpline, t::Float64)::Float64
    n = length(s.t_s)
    (t < s.t_s[1] || t > s.t_s[n]) && return 0.0
    i, tc = _spline_segment(s, t)
    @inbounds begin
        h = s.t_s[i + 1] - s.t_s[i]
        a = (s.t_s[i + 1] - tc) / h
        b = (tc - s.t_s[i]) / h
        return (s.y[i + 1] - s.y[i]) / h +
            ((1 - 3 * a^2) * s.d2[i] + (3 * b^2 - 1) * s.d2[i + 1]) * h / 6
    end
end

"""
    ReactionWheelMomentumState(num_spacecraft)

Per-spacecraft record of what the effector last applied: the wheel momentum in
the body frame (N m s) and the reaction torque it produced (N m). A run with
`isolate_state=false` can read it afterwards; the values are diagnostics, and a
save field that wants the exact value at a save time should evaluate
[`wheel_momentum_body`](@ref) at that time instead of reading the last stage.
"""
mutable struct ReactionWheelMomentumState
    momentum_body::Vector{SVector{3, Float64}}
    torque_body::Vector{SVector{3, Float64}}
end

function ReactionWheelMomentumState(num_spacecraft::Int)
    zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    return ReactionWheelMomentumState(fill(zero3, num_spacecraft), fill(zero3, num_spacecraft))
end

"""
    ReactionWheelMomentumModel(t_s, speeds_rad_s, spin_axes, wheel_inertia_kg_m2; kwargs...)

Force-torque effector for a wheel assembly whose speeds are a prescribed
function of time.

- `t_s` is the sample time of each row of `speeds_rad_s`, in the same clock the
  table was recorded on, strictly increasing.
- `speeds_rad_s` is `N x K`: one column per wheel, signed angular speed about
  that wheel's spin axis, in radians per second.
- `spin_axes` is `3 x K`: wheel `k`'s spin axis in the body frame, in the same
  column order. The columns need not be unit vectors; whatever scale they carry
  multiplies through into the momentum, which is what a matrix regressed from
  flight data will give.
- `wheel_inertia_kg_m2` is the spin-axis moment of inertia of one wheel, common
  to all of them.

so that `H_w(t) = spin_axes * (wheel_inertia * speeds(t))`.

Keyword arguments:

- `spacecraft_index`: which spacecraft of the run the assembly belongs to; the
  effector returns a zero wrench for every other one. Default 1.
- `time_offset_s`: the table's clock minus the run's clock, so the effector
  reads the table at `t + time_offset_s` when the run starts partway into the
  table. Default 0.
- `num_spacecraft`: how many entries the diagnostic state carries. Defaults to
  `spacecraft_index`.
- `state`: an existing [`ReactionWheelMomentumState`](@ref) to write into.
"""
struct ReactionWheelMomentumModel <: AbstractForceTorqueModel
    spacecraft_index::Int
    "wheel_inertia * spin axis, one entry per wheel: N m s per rad/s"
    momentum_axes::Vector{SVector{3, Float64}}
    speeds::Vector{WheelSpeedSpline}
    time_offset_s::Float64
    state::ReactionWheelMomentumState
end

function ReactionWheelMomentumModel(
    t_s::AbstractVector{<:Real},
    speeds_rad_s::AbstractMatrix{<:Real},
    spin_axes::AbstractMatrix{<:Real},
    wheel_inertia_kg_m2::Real;
    spacecraft_index::Int=1,
    time_offset_s::Real=0.0,
    num_spacecraft::Int=spacecraft_index,
    state::ReactionWheelMomentumState=ReactionWheelMomentumState(num_spacecraft),
)
    n_wheels = size(speeds_rad_s, 2)
    size(spin_axes) == (3, n_wheels) || throw(ArgumentError(
        "ReactionWheelMomentumModel: spin_axes must be 3 x $(n_wheels) to match speeds_rad_s, got $(size(spin_axes))."))
    size(speeds_rad_s, 1) == length(t_s) || throw(ArgumentError(
        "ReactionWheelMomentumModel: speeds_rad_s must have one row per entry of t_s, got $(size(speeds_rad_s, 1)) and $(length(t_s))."))
    n_wheels >= 1 || throw(ArgumentError("ReactionWheelMomentumModel: at least one wheel is required."))
    spacecraft_index >= 1 || throw(ArgumentError("ReactionWheelMomentumModel: spacecraft_index must be positive."))
    isfinite(wheel_inertia_kg_m2) && wheel_inertia_kg_m2 > 0.0 || throw(ArgumentError(
        "ReactionWheelMomentumModel: wheel_inertia_kg_m2 must be finite and positive, got $(wheel_inertia_kg_m2)."))
    length(state.momentum_body) >= spacecraft_index || throw(ArgumentError(
        "ReactionWheelMomentumModel: the state carries $(length(state.momentum_body)) spacecraft, fewer than spacecraft_index=$(spacecraft_index)."))

    axes = [SVector{3, Float64}(
        Float64(wheel_inertia_kg_m2) * Float64(spin_axes[1, k]),
        Float64(wheel_inertia_kg_m2) * Float64(spin_axes[2, k]),
        Float64(wheel_inertia_kg_m2) * Float64(spin_axes[3, k])) for k in 1:n_wheels]
    splines = [wheel_speed_spline(t_s, @view speeds_rad_s[:, k]) for k in 1:n_wheels]
    return ReactionWheelMomentumModel(spacecraft_index, axes, splines, Float64(time_offset_s), state)
end

"""
    wheel_speeds_rad_s(model, t) -> SVector

Interpolated wheel speeds (rad/s) at run time `t`, for at most three wheels;
the general case is `wheel_spline_value(model.speeds[k], t + model.time_offset_s)`.
"""
@inline function wheel_speeds_rad_s(model::ReactionWheelMomentumModel, t::Float64)::SVector{3, Float64}
    tt = t + model.time_offset_s
    n = length(model.speeds)
    return SVector{3, Float64}(
        n >= 1 ? wheel_spline_value(model.speeds[1], tt) : 0.0,
        n >= 2 ? wheel_spline_value(model.speeds[2], tt) : 0.0,
        n >= 3 ? wheel_spline_value(model.speeds[3], tt) : 0.0,
    )
end

"""
    wheel_momentum_body(model, t) -> SVector{3, Float64}

`H_w(t)`: the total wheel angular momentum in the body frame, N m s.
"""
@inline function wheel_momentum_body(model::ReactionWheelMomentumModel, t::Float64)::SVector{3, Float64}
    tt = t + model.time_offset_s
    h = SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for k in eachindex(model.speeds)
        h = h + wheel_spline_value(model.speeds[k], tt) * model.momentum_axes[k]
    end
    return h
end

"""
    wheel_momentum_rate_body(model, t) -> SVector{3, Float64}

`Ḣ_w(t)`: the body-frame rate of change of the wheel momentum, N m, taken from
the analytic derivative of the wheel-speed splines.
"""
@inline function wheel_momentum_rate_body(model::ReactionWheelMomentumModel, t::Float64)::SVector{3, Float64}
    tt = t + model.time_offset_s
    hdot = SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for k in eachindex(model.speeds)
        hdot = hdot + wheel_spline_derivative(model.speeds[k], tt) * model.momentum_axes[k]
    end
    return hdot
end

"""
    wheel_reaction_torque(h_wheel_body, h_wheel_rate_body, omega_body) -> SVector{3, Float64}

`τ = -Ḣ_w - ω × H_w`, the torque a wheel assembly of body-frame momentum `H_w`
applies to the body it is mounted in. Constant wheel momentum in a non-rotating
body gives exactly zero.
"""
@inline function wheel_reaction_torque(
    h_wheel_body::SVector{3, Float64},
    h_wheel_rate_body::SVector{3, Float64},
    omega_body::SVector{3, Float64},
)::SVector{3, Float64}
    return -h_wheel_rate_body - cross(omega_body, h_wheel_body)
end

@inline environment_requirements(::ReactionWheelMomentumModel) = EffectorEnvironmentRequirements()

function wrench(
    model::ReactionWheelMomentumModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return _reaction_wheel_wrench(model, x, t, model.spacecraft_index)
end

function wrench_caching!(
    model::ReactionWheelMomentumModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
    p,
    sat_idx::Int,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return _reaction_wheel_wrench(model, x, t, sat_idx)
end

@inline function _reaction_wheel_wrench(
    model::ReactionWheelMomentumModel,
    x::StateSample,
    t::Float64,
    sat_idx::Int,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    sat_idx == model.spacecraft_index || return (zero3, zero3)
    omega = x.ω_body
    omega === nothing && return (zero3, zero3)
    h = wheel_momentum_body(model, t)
    hdot = wheel_momentum_rate_body(model, t)
    torque = wheel_reaction_torque(h, hdot, omega)
    @inbounds if sat_idx <= length(model.state.momentum_body)
        model.state.momentum_body[sat_idx] = h
        model.state.torque_body[sat_idx] = torque
    end
    return (zero3, torque)
end

end # module ReactionWheelMomentum
