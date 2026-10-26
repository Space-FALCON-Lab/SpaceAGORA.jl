# The built-in sensor kinds: rate gyro, star tracker, magnetometer. Included by
# `Components` (src/models/vehicle/spacecraft/components.jl), which owns the
# records; the contract they implement is the `Sensors` module (sensors.jl).

const _ZERO_V3 = SVector{3, Float64}(0.0, 0.0, 0.0)
const _NAN_V3 = SVector{3, Float64}(NaN, NaN, NaN)

@inline _gaussian3(rng) = SVector{3, Float64}(randn(rng), randn(rng), randn(rng))

"""
    RateGyro(; period_s=0.1, angle_random_walk=0.0, rate_random_walk=0.0, initial_bias=zero, seed=1)

A three-axis rate gyro: it measures the body rate in its link's frame, plus a
bias and white noise. `angle_random_walk` (rad/√s) sets the white noise,
`angle_random_walk / √period_s` per axis and sample; the bias starts at
`initial_bias` (rad/s) and walks with `rate_random_walk` (rad/s/√s). `seed`
selects its random stream. Zero noise measures the rate exactly.
"""
mutable struct RateGyro <: AbstractSensor
    period_s::Float64
    angle_random_walk::Float64
    rate_random_walk::Float64
    initial_bias::SVector{3, Float64}
    seed::Int
    # Run state: the current bias, the latest measurement and its time, the random stream.
    bias::SVector{3, Float64}
    latest::SVector{3, Float64}
    latest_time_s::Float64
    rng::Any
end

function RateGyro(; period_s::Real=0.1, angle_random_walk::Real=0.0, rate_random_walk::Real=0.0,
                  initial_bias=_ZERO_V3, seed::Integer=1)
    (isfinite(period_s) && period_s > 0.0) || throw(ArgumentError("RateGyro period_s must be finite and positive, got $(period_s)."))
    (isfinite(angle_random_walk) && angle_random_walk >= 0.0) || throw(ArgumentError("RateGyro angle_random_walk must be finite and nonnegative."))
    (isfinite(rate_random_walk) && rate_random_walk >= 0.0) || throw(ArgumentError("RateGyro rate_random_walk must be finite and nonnegative."))
    b = SVector{3, Float64}(initial_bias)
    return RateGyro(Float64(period_s), Float64(angle_random_walk), Float64(rate_random_walk), b, Int(seed), b, _NAN_V3, NaN, nothing)
end

sensor_dimension(::RateGyro) = 3

function sensor_measure(g::RateGyro, truth::SensorTruth, ctx::ActuatorContext, rng)
    rate = ctx.link_to_body' * truth.ω_body
    noise = g.angle_random_walk > 0.0 ? (g.angle_random_walk / sqrt(g.period_s)) * _gaussian3(rng) : _ZERO_V3
    measurement = rate + g.bias + noise
    if g.rate_random_walk > 0.0
        g.bias = g.bias + (g.rate_random_walk * sqrt(g.period_s)) * _gaussian3(rng)
    end
    return measurement
end

"""
    StarTracker(; period_s=1.0, cross_boresight_noise_rad=0.0, about_boresight_noise_rad=0.0, seed=1)

A star tracker whose boresight is its link's +z axis: it measures the
attitude quaternion (inertial to body, scalar last), perturbed by a small
rotation with `cross_boresight_noise_rad` about the two axes across the
boresight and `about_boresight_noise_rad` about it (one standard deviation,
rad). Zero noise measures the attitude exactly.
"""
mutable struct StarTracker <: AbstractSensor
    period_s::Float64
    cross_boresight_noise_rad::Float64
    about_boresight_noise_rad::Float64
    seed::Int
    latest::SVector{4, Float64}
    latest_time_s::Float64
    rng::Any
end

function StarTracker(; period_s::Real=1.0, cross_boresight_noise_rad::Real=0.0, about_boresight_noise_rad::Real=0.0, seed::Integer=1)
    (isfinite(period_s) && period_s > 0.0) || throw(ArgumentError("StarTracker period_s must be finite and positive, got $(period_s)."))
    (isfinite(cross_boresight_noise_rad) && cross_boresight_noise_rad >= 0.0) || throw(ArgumentError("StarTracker cross_boresight_noise_rad must be finite and nonnegative."))
    (isfinite(about_boresight_noise_rad) && about_boresight_noise_rad >= 0.0) || throw(ArgumentError("StarTracker about_boresight_noise_rad must be finite and nonnegative."))
    return StarTracker(Float64(period_s), Float64(cross_boresight_noise_rad), Float64(about_boresight_noise_rad), Int(seed),
                       SVector{4, Float64}(NaN, NaN, NaN, NaN), NaN, nothing)
end

sensor_dimension(::StarTracker) = 4

# The rotation matrix of the rotation vector θ (Rodrigues).
@inline function _rotation_matrix(θ::SVector{3, Float64})::SMatrix{3, 3, Float64, 9}
    a = norm(θ)
    a < 1e-300 && return SMatrix{3, 3, Float64, 9}(1, 0, 0, 0, 1, 0, 0, 0, 1)
    k = θ / a
    K = SMatrix{3, 3, Float64, 9}(0.0, k[3], -k[2], -k[3], 0.0, k[1], k[2], -k[1], 0.0)
    return SMatrix{3, 3, Float64, 9}(1, 0, 0, 0, 1, 0, 0, 0, 1) + sin(a) * K + (1 - cos(a)) * (K * K)
end

function sensor_measure(st::StarTracker, truth::SensorTruth, ctx::ActuatorContext, rng)
    q = truth.q_ib / norm(truth.q_ib)
    (st.cross_boresight_noise_rad == 0.0 && st.about_boresight_noise_rad == 0.0) && return q
    to_body = ctx.link_to_body
    n = _gaussian3(rng)
    # A small rotation of the body frame: across the boresight (link x, y) and about it (link z).
    θ = to_body * SVector{3, Float64}(st.cross_boresight_noise_rad * n[1], st.cross_boresight_noise_rad * n[2],
                                      st.about_boresight_noise_rad * n[3])
    measured = _rotation_matrix(θ)' * rot(q)        # rot(q) maps inertial to body
    out = dcm_to_quaternion(SMatrix{3, 3, Float64, 9}(measured))
    return dot(out, q) < 0.0 ? -out : out           # the same hemisphere as the truth
end

"""
    Magnetometer(; period_s=0.1, noise_t=0.0, bias_t=zero, seed=1)

A three-axis magnetometer: it measures the local magnetic field in its link's
frame (T), from the run's torque-rod field model (dipole or IGRF, Earth
only), plus a constant `bias_t` and white noise `noise_t` per axis and sample.
"""
mutable struct Magnetometer <: AbstractSensor
    period_s::Float64
    noise_t::Float64
    bias_t::SVector{3, Float64}
    seed::Int
    latest::SVector{3, Float64}
    latest_time_s::Float64
    rng::Any
end

function Magnetometer(; period_s::Real=0.1, noise_t::Real=0.0, bias_t=_ZERO_V3, seed::Integer=1)
    (isfinite(period_s) && period_s > 0.0) || throw(ArgumentError("Magnetometer period_s must be finite and positive, got $(period_s)."))
    (isfinite(noise_t) && noise_t >= 0.0) || throw(ArgumentError("Magnetometer noise_t must be finite and nonnegative."))
    return Magnetometer(Float64(period_s), Float64(noise_t), SVector{3, Float64}(bias_t), Int(seed), _NAN_V3, NaN, nothing)
end

sensor_dimension(::Magnetometer) = 3
sensor_requirements(::Magnetometer) = ActuatorRequirements(magnetic_field=true)

function sensor_measure(m::Magnetometer, truth::SensorTruth, ctx::ActuatorContext, rng)
    field = ctx.link_to_body' * truth.magnetic_field_body
    noise = m.noise_t > 0.0 ? m.noise_t * _gaussian3(rng) : _ZERO_V3
    return field + m.bias_t + noise
end

"""
    GNSSReceiver(; period_s=1.0, position_noise_m=0.0, velocity_noise_mps=0.0, frame=:inertial, max_altitude_m=Inf, seed=1)

A GNSS receiver: it measures the position and velocity (m, m/s; six
components, position first) in the run's inertial frame, or with
`frame=:planet_fixed` in the planet-fixed frame, with white noise
`position_noise_m` and `velocity_noise_mps` per axis and sample. Above
`max_altitude_m` over the planet it has no fix (NaN), as a receiver above the
navigation constellation's coverage. The antenna is taken at the centre of
mass.
"""
mutable struct GNSSReceiver <: AbstractSensor
    period_s::Float64
    position_noise_m::Float64
    velocity_noise_mps::Float64
    frame::Symbol
    max_altitude_m::Float64
    seed::Int
    latest::SVector{6, Float64}
    latest_time_s::Float64
    rng::Any
end

function GNSSReceiver(; period_s::Real=1.0, position_noise_m::Real=0.0, velocity_noise_mps::Real=0.0, frame::Symbol=:inertial,
                      max_altitude_m::Real=Inf, seed::Integer=1)
    (isfinite(period_s) && period_s > 0.0) || throw(ArgumentError("GNSSReceiver period_s must be finite and positive, got $(period_s)."))
    (isfinite(position_noise_m) && position_noise_m >= 0.0) || throw(ArgumentError("GNSSReceiver position_noise_m must be finite and nonnegative."))
    (isfinite(velocity_noise_mps) && velocity_noise_mps >= 0.0) || throw(ArgumentError("GNSSReceiver velocity_noise_mps must be finite and nonnegative."))
    frame in (:inertial, :planet_fixed) || throw(ArgumentError("GNSSReceiver frame must be :inertial or :planet_fixed, got $(repr(frame))."))
    max_altitude_m > 0.0 || throw(ArgumentError("GNSSReceiver max_altitude_m must be positive (Inf for no cutoff), got $(max_altitude_m)."))
    return GNSSReceiver(Float64(period_s), Float64(position_noise_m), Float64(velocity_noise_mps), frame, Float64(max_altitude_m), Int(seed),
                        SVector{6, Float64}(ntuple(_ -> NaN, 6)), NaN, nothing)
end

sensor_dimension(::GNSSReceiver) = 6
sensor_requirements(g::GNSSReceiver) = ActuatorRequirements(planet_frame=(g.frame === :planet_fixed || isfinite(g.max_altitude_m)))

function sensor_measure(g::GNSSReceiver, truth::SensorTruth, ctx::ActuatorContext, rng)
    isfinite(g.max_altitude_m) && !(truth.altitude_m <= g.max_altitude_m) && return SVector{6, Float64}(ntuple(_ -> NaN, 6))
    pos0, vel0 = g.frame === :planet_fixed ? (truth.pos_pp, truth.vel_pp) : (truth.pos_ii, truth.vel_ii)
    pos = g.position_noise_m > 0.0 ? pos0 + g.position_noise_m * _gaussian3(rng) : pos0
    vel = g.velocity_noise_mps > 0.0 ? vel0 + g.velocity_noise_mps * _gaussian3(rng) : vel0
    return SVector{6, Float64}(pos[1], pos[2], pos[3], vel[1], vel[2], vel[3])
end

"""
    SunSensor(; period_s=0.1, noise_rad=0.0, field_of_view_rad=π/2, min_lit_fraction=0.5, albedo=0.0, seed=1)

A sun sensor whose boresight is its link's +z axis: it measures the direction
to the Sun as a unit vector in its link's frame, turned by a small rotation
with `noise_rad` (one standard deviation) about each axis. It reports no
direction (NaN) when the Sun is more than `field_of_view_rad` from the
boresight or when less than `min_lit_fraction` of the Sun's disk is lit
(eclipse). The Sun's direction comes from the run's ephemerides, the lit
fraction from the eclipse model of solar radiation pressure.

`albedo` (0 to 1, about 0.3 for the Earth; default 0, none) adds the sunlight
the planet reflects, as a coarse sun sensor sees it: while the planet's disk
is within the field of view, the reading is the normalized sum of the Sun's
direction and the nadir direction weighted by `η = albedo (2/3) (R/r)² Φ(α)`,
the irradiance of a Lambertian sphere of radius `R` at distance `r` relative
to the Sun's, with the Lambert phase law `Φ(α) = (sin α + (π - α) cos α)/π`
and `α` the angle at the planet between the Sun and the sensor (Φ is 1 over
the subsolar point, 0 over the night side). The reflected light is taken from
the planet's centre; reflectivity maps and the planet's infrared are not
modelled.
"""
mutable struct SunSensor <: AbstractSensor
    period_s::Float64
    noise_rad::Float64
    field_of_view_rad::Float64
    min_lit_fraction::Float64
    albedo::Float64
    seed::Int
    latest::SVector{3, Float64}
    latest_time_s::Float64
    rng::Any
end

function SunSensor(; period_s::Real=0.1, noise_rad::Real=0.0, field_of_view_rad::Real=π / 2, min_lit_fraction::Real=0.5,
                   albedo::Real=0.0, seed::Integer=1)
    (isfinite(period_s) && period_s > 0.0) || throw(ArgumentError("SunSensor period_s must be finite and positive, got $(period_s)."))
    (isfinite(noise_rad) && noise_rad >= 0.0) || throw(ArgumentError("SunSensor noise_rad must be finite and nonnegative."))
    (0.0 < field_of_view_rad <= π) || throw(ArgumentError("SunSensor field_of_view_rad must be in (0, π], got $(field_of_view_rad)."))
    (0.0 <= min_lit_fraction <= 1.0) || throw(ArgumentError("SunSensor min_lit_fraction must be in [0, 1], got $(min_lit_fraction)."))
    (0.0 <= albedo <= 1.0) || throw(ArgumentError("SunSensor albedo must be in [0, 1], got $(albedo)."))
    return SunSensor(Float64(period_s), Float64(noise_rad), Float64(field_of_view_rad), Float64(min_lit_fraction), Float64(albedo),
                     Int(seed), _NAN_V3, NaN, nothing)
end

sensor_dimension(::SunSensor) = 3
sensor_requirements(::SunSensor) = ActuatorRequirements(solar=true)

function sensor_measure(s::SunSensor, truth::SensorTruth, ctx::ActuatorContext, rng)
    d = ctx.link_to_body' * truth.sun_direction_body
    (all(isfinite, d) && truth.sun_lit_fraction >= s.min_lit_fraction) || return _NAN_V3
    acos(clamp(d[3], -1.0, 1.0)) <= s.field_of_view_rad || return _NAN_V3
    s.albedo > 0.0 && (d = _with_albedo(s, d, truth, ctx))
    s.noise_rad > 0.0 || return d
    turned = _rotation_matrix(s.noise_rad * _gaussian3(rng)) * d
    return turned / norm(turned)
end

# The Sun's direction `d` (sensor frame) with the planet's reflected light
# added along the nadir while the planet's disk is in view (see `SunSensor`).
function _with_albedo(s::SunSensor, d::SVector{3, Float64}, truth::SensorTruth, ctx::ActuatorContext)::SVector{3, Float64}
    R = truth.planet_radius_m
    r = norm(truth.pos_ii)
    (isfinite(R) && r > R) || return d
    up = truth.pos_ii / r
    sun_ii = rot(truth.q_ib)' * truth.sun_direction_body
    α = acos(clamp(dot(up, sun_ii), -1.0, 1.0))
    η = s.albedo * (2.0 / 3.0) * (R / r)^2 * (sin(α) + (π - α) * cos(α)) / π
    η > 0.0 || return d
    nadir = ctx.link_to_body' * (rot(truth.q_ib) * (-up))
    acos(clamp(nadir[3], -1.0, 1.0)) <= s.field_of_view_rad + asin(R / r) || return d
    v = d + η * nadir
    return v / norm(v)
end

"""
    Accelerometer(; period_s=0.1, noise_mps2=0.0, bias_mps2=zero, position_m=zero, seed=1)

A three-axis accelerometer at `position_m` in its link's frame: it measures
the specific force there, the acceleration less gravity (m/s²), in its link's
frame, plus a constant `bias_mps2` and white noise `noise_mps2` per axis and
sample. Away from the centre of mass (the body origin), at `r`, the rotation
adds `ω̇ × r + ω × (ω × r)` to the centre of mass's specific force.
Gravity is the force of the run's effectors for which `is_gravitational` is
true, so in free fall it reads zero and under thrust, drag or radiation
pressure it reads their acceleration.
"""
mutable struct Accelerometer <: AbstractSensor
    period_s::Float64
    noise_mps2::Float64
    bias_mps2::SVector{3, Float64}
    position_m::SVector{3, Float64}
    seed::Int
    latest::SVector{3, Float64}
    latest_time_s::Float64
    rng::Any
end

function Accelerometer(; period_s::Real=0.1, noise_mps2::Real=0.0, bias_mps2=_ZERO_V3, position_m=_ZERO_V3, seed::Integer=1)
    (isfinite(period_s) && period_s > 0.0) || throw(ArgumentError("Accelerometer period_s must be finite and positive, got $(period_s)."))
    (isfinite(noise_mps2) && noise_mps2 >= 0.0) || throw(ArgumentError("Accelerometer noise_mps2 must be finite and nonnegative."))
    all(isfinite, position_m) || throw(ArgumentError("Accelerometer position_m must be finite."))
    return Accelerometer(Float64(period_s), Float64(noise_mps2), SVector{3, Float64}(bias_mps2), SVector{3, Float64}(position_m),
                         Int(seed), _NAN_V3, NaN, nothing)
end

sensor_dimension(::Accelerometer) = 3
sensor_needs_specific_force(::Accelerometer) = true

function sensor_measure(a::Accelerometer, truth::SensorTruth, ctx::ActuatorContext, rng)
    r = ctx.link_r_body + ctx.link_to_body * a.position_m
    force = truth.specific_force_body
    iszero(r) || (force = force + cross(truth.angular_acceleration_body, r) + cross(truth.ω_body, cross(truth.ω_body, r)))
    f = ctx.link_to_body' * force + a.bias_mps2
    return a.noise_mps2 > 0.0 ? f + a.noise_mps2 * _gaussian3(rng) : f
end
