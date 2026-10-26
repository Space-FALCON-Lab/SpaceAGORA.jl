# The reaction-wheel actuator kind. Included by `Components`
# (src/models/vehicle/spacecraft/components.jl), which owns and exports the
# records, so serialized spacecraft keep their type identity; the actuator
# contract it implements is in `actuators.jl`.

# ---- Record ---------------------------------------------------------------

"""
    ReactionWheel(; spin_axis, max_torque_nm=Inf, max_momentum_nms=Inf,
                  spin_inertia_kgm2=NaN, base_power_w=0.0, motor_efficiency=1.0, regenerative_efficiency=0.0)

One reaction wheel, mounted on the spacecraft's root link with
`add_actuator!`. Its level is the torque it puts on the vehicle about
its spin axis, in newton-metres; the wheel's own angular momentum along that
axis changes at minus the same rate. `spin_axis` is in the root-link (body)
frame. `max_torque_nm` bounds the level and `max_momentum_nms` stops a torque
that would store more momentum; both default to unlimited.

Power (see `wheel_power`): the wheel draws `base_power_w` (electronics, bearing
losses) plus the mechanical power its motor puts into the wheel, divided by
`motor_efficiency`. The mechanical power is the motor torque on the wheel
(minus the level) times the wheel speed, its stored momentum over
`spin_inertia_kgm2`; without a spin inertia (the default) only the base power
counts. While the motor slows the wheel, `regenerative_efficiency` of the
mechanical power returns to the bus (zero by default: nothing returns, nothing
extra is drawn).
"""
struct ReactionWheel <: AbstractActuator
    spin_axis::SVector{3, Float64}
    max_torque_nm::Float64
    max_momentum_nms::Float64
    spin_inertia_kgm2::Float64
    base_power_w::Float64
    motor_efficiency::Float64
    regenerative_efficiency::Float64

    function ReactionWheel(spin_axis, max_torque_nm::Real, max_momentum_nms::Real, spin_inertia_kgm2::Real,
                           base_power_w::Real, motor_efficiency::Real, regenerative_efficiency::Real)
        axis = SVector{3, Float64}(spin_axis)
        n = norm(axis)
        (isfinite(n) && n > 0.0) || throw(ArgumentError("ReactionWheel spin_axis must be finite and nonzero, got $(axis)."))
        (max_torque_nm > 0.0 && !isnan(max_torque_nm)) || throw(ArgumentError("ReactionWheel max_torque_nm must be positive, got $(max_torque_nm)."))
        (max_momentum_nms > 0.0 && !isnan(max_momentum_nms)) || throw(ArgumentError("ReactionWheel max_momentum_nms must be positive, got $(max_momentum_nms)."))
        (isnan(spin_inertia_kgm2) || (isfinite(spin_inertia_kgm2) && spin_inertia_kgm2 > 0.0)) ||
            throw(ArgumentError("ReactionWheel spin_inertia_kgm2 must be positive and finite (or NaN for none), got $(spin_inertia_kgm2)."))
        (isfinite(base_power_w) && base_power_w >= 0.0) || throw(ArgumentError("ReactionWheel base_power_w must be finite and nonnegative, got $(base_power_w)."))
        (0.0 < motor_efficiency <= 1.0) || throw(ArgumentError("ReactionWheel motor_efficiency must be in (0, 1], got $(motor_efficiency)."))
        (0.0 <= regenerative_efficiency <= 1.0) || throw(ArgumentError("ReactionWheel regenerative_efficiency must be in [0, 1], got $(regenerative_efficiency)."))
        return new(axis / n, Float64(max_torque_nm), Float64(max_momentum_nms), Float64(spin_inertia_kgm2),
                   Float64(base_power_w), Float64(motor_efficiency), Float64(regenerative_efficiency))
    end
end

ReactionWheel(spin_axis, max_torque_nm::Real, max_momentum_nms::Real) =
    ReactionWheel(spin_axis, max_torque_nm, max_momentum_nms, NaN, 0.0, 1.0, 0.0)

ReactionWheel(; spin_axis, max_torque_nm::Real=Inf, max_momentum_nms::Real=Inf, spin_inertia_kgm2::Real=NaN,
              base_power_w::Real=0.0, motor_efficiency::Real=1.0, regenerative_efficiency::Real=0.0) =
    ReactionWheel(spin_axis, max_torque_nm, max_momentum_nms, spin_inertia_kgm2, base_power_w, motor_efficiency, regenerative_efficiency)

# ---- Methods of the actuator contract --------------------------------------

actuator_limits(wheel::ReactionWheel) = (-wheel.max_torque_nm, wheel.max_torque_nm)

"""
    actuator_wrench(wheel::ReactionWheel, level, ctx) -> (force_body, torque_body)

A torque of `level` newton-metres on the vehicle about the wheel's spin axis,
in the body frame; no force.
"""
@inline function actuator_wrench(wheel::ReactionWheel, level::Real, ctx::ActuatorContext)
    return SVector{3, Float64}(0.0, 0.0, 0.0), ctx.link_to_body * (Float64(level) * wheel.spin_axis)
end

"""
    wheel_power(wheel::ReactionWheel, level, h_nms) -> Float64

The electrical power (W) a wheel storing `h_nms` along its spin axis draws
while it puts a torque of `level` N m on the vehicle: its base power, plus the
mechanical power into the wheel, `-level * h_nms / spin_inertia_kgm2`, over the
motor efficiency while the motor speeds the wheel up, minus the regenerated
share of it while the motor slows the wheel down. A NaN `h_nms`, or a wheel
without a spin inertia, gives the base power.
"""
function wheel_power(wheel::ReactionWheel, level::Real, h_nms::Real)::Float64
    power = wheel.base_power_w
    (isfinite(wheel.spin_inertia_kgm2) && isfinite(h_nms) && isfinite(level)) || return power
    mechanical = -Float64(level) * (Float64(h_nms) / wheel.spin_inertia_kgm2)
    return mechanical >= 0.0 ? power + mechanical / wheel.motor_efficiency :
        power + wheel.regenerative_efficiency * mechanical
end

"""
    actuator_consumption(wheel::ReactionWheel, level, ctx) -> (mass_rate_kg_s, power_w)

No propellant. The power is the base power: the context does not carry the
wheel's stored momentum, which the mechanical part needs. `mounted_power`, and
the `actuator_power` result column, use `wheel_power` with the integrated
momentum.
"""
@inline actuator_consumption(wheel::ReactionWheel, ::Real, ::ActuatorContext) =
    (mass_rate_kg_s = 0.0, power_w = wheel.base_power_w)

"""
    wheel_momentum_saturate(wheel::ReactionWheel, level, h_nms) -> Float64

The level a wheel storing `h_nms` along its spin axis can still produce: zero
when the wheel is at its momentum limit and the torque would store more
(the wheel's momentum changes at `-level`), `level` otherwise.
"""
@inline function wheel_momentum_saturate(wheel::ReactionWheel, level::Float64, h_nms::Float64)::Float64
    limit = wheel.max_momentum_nms
    (isfinite(limit) && abs(h_nms) >= limit && h_nms * (-level) > 0.0) && return 0.0
    return level
end

# ---- The assembly the engine integrates -------------------------------------

"""
    ReactionWheelAssembly

The reaction wheels of a link, whose stored momenta the engine integrates for
the root link. Wheels come either through the `Link` keywords `J_rw`,
`max_torque` and `max_h` (their limits are recorded but not applied) or as
mounted [`ReactionWheel`](@ref) actuators, whose limits are applied. The wheel
count is a runtime field, not a type parameter; see the note on the struct for
why.
"""
# Deliberately NOT parameterised on the wheel count. The parameter existed only
# to size J_rw/h_wheels statically, but it propagated outward as Link{N_RW},
# which made `Link` an abstract UnionAll and therefore made
# SpacecraftModel.links::Vector{Link} and .root::Link abstract fields. Every
# read of a link field in the RHS was then a dynamic lookup that boxed a
# Float64: measured at ~835 bytes allocated per satellite per RHS call, which
# accounted for essentially all of the 1.53 GB a 10000-satellite run allocated,
# and a 35.7x slowdown on field access in isolation. The wheel arrays are
# touched once per satellite per attitude step; the link fields are touched
# constantly, so the static sizing was being bought at far too high a price.
@kwdef mutable struct ReactionWheelAssembly
    # --- Parameters (Fixed) ---
    n_wheels::Int = 0 # Number of reaction wheels in the assembly

    "Jacobian mapping wheel angular velocity to body angular momentum"
    J_rw::Matrix{Float64} = zeros(Float64, 3, 0)

    "Precomputed pseudo-inverse of J_rw, mapping a desired body torque/momentum onto per-wheel values"
    J_rw_pinv::Matrix{Float64} = Matrix{Float64}(pinv(J_rw))

    "Maximum torque (Nm) *each* wheel can produce"
    max_wheel_torque::Float64 = 0.1
    
    "Maximum angular momentum (Nms) *each* wheel can store"
    max_wheel_h::Float64 = 0.1

    # --- State (Mutable) ---
    "Angular momentum (h) of each wheel"
    h_wheels::Vector{Float64} = zeros(Float64, n_wheels)

    "Derivative of wheel angular momentum (h_dot)"
    h_dot_wheels::Vector{Float64} = zeros(Float64, n_wheels)
    
    "Net torque (τ) applied to the body *by* the RW assembly"
    tau_body_net::MVector{3, Float64} = MVector{3, Float64}(zeros(3))

    "Wheels mounted with add_actuator!, the last columns of J_rw; wheels given through J_rw directly come first"
    wheels::Vector{ReactionWheel} = ReactionWheel[]

    "Momentum rate of each mounted wheel (N m) commanded in the current right-hand-side evaluation: minus its saturated level, written by the engine"
    mounted_rates::Vector{Float64} = Float64[]

    "Control moment gyros mounted with add_actuator! (their gimbal angles are the `cmg_angles` state)"
    cmgs::Vector{ControlMomentGyro} = ControlMomentGyro[]

    "Gimbal rate of each mounted CMG (rad/s) commanded in the current right-hand-side evaluation, written by the engine; the outer gimbals' first, then the inner gimbals', as the `cmg_angles` state"
    cmg_rates::Vector{Float64} = Float64[]

    "Inner gimbals of double-gimbal CMGs mounted with add_actuator!, after every CMG in the `cmg_angles` state and the mounting order"
    cmg_inner::Vector{CMGInnerGimbal} = CMGInnerGimbal[]

    "For each mounted CMG, the position of its inner gimbal in `cmg_inner` (0 for a single-gimbal CMG)"
    cmg_inner_of::Vector{Int} = Int[]
end
