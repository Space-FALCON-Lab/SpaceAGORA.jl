module SpacecraftModels

using StaticArrays
using LinearAlgebra

using ..Components
using ..AbstractTypes: AbstractActuator, AbstractSensor
using ..Geodesy: geodetic_altitude, ellipsoid_surface_radius, radius_for_geodetic_altitude
using ..OrbitalGeometry
using ..EphemeridesModels: SpiceEphemeridesModel, ephemerides_epoch_seconds, ephemerides_frame_lpi

export Link, Joint, CompliantAttachment, attachment_body_count, SpacecraftModel, DynamicsModel, InitialCondition, CartesianInitialCondition, AbstractInitialCondition, GuidanceModel, NavigationModel, ControlModel
export ConstellationSlot, OrbitalPlane, build_constellation, walker_planes, walker_constellation

const I3 = SMatrix{3, 3, Float64}(diagm(ones(3)))
const DEFAULT_INITIAL_CONDITION_Q = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
const DEFAULT_INITIAL_CONDITION_ANG_VEL = SVector{3, Float64}(0.0, 0.0, 0.0)

abstract type AbstractInitialCondition end

# An initial attitude quaternion (x, y, z, w) or body rate (rad/s) given as any vector
# or tuple of 4 or 3 numbers; an SVector passes through with the same values.
function _initial_attitude(q)::SVector{4, Float64}
    length(q) == 4 || throw(ArgumentError(
        "An initial condition's q is a quaternion (x, y, z, w), scalar last: give 4 numbers, not $(length(q))."))
    return SVector{4, Float64}(Tuple(q))
end

function _initial_rate(ang_vel)::SVector{3, Float64}
    length(ang_vel) == 3 || throw(ArgumentError(
        "An initial condition's ang_vel is a body rate in rad/s: give 3 numbers, not $(length(ang_vel))."))
    return SVector{3, Float64}(Tuple(ang_vel))
end

struct InitialCondition <: AbstractInitialCondition
    a::Float64 # Semimajor axis (m)
    e::Float64 # Eccentricity (nd)
    i::Float64 # Inclination (rad)
    ω::Float64 # Argument of periapsis (rad)
    Ω::Float64 # RAAN (rad)
    ν::Float64 # True anomaly (rad)
    q::SVector{4, Float64} # Initial orientation quaternion (x, y, z, w)
    ang_vel::SVector{3, Float64} # Initial angular velocity (rad/s)

    function InitialCondition(
        a::Float64,
        e::Float64,
        i::Float64,
        ω::Float64,
        Ω::Float64,
        ν::Float64,
        q::SVector{4, Float64},
        ang_vel::SVector{3, Float64},
        ::Val{:radians}
    )
        return new(a, e, i, ω, Ω, ν, q, ang_vel)
    end
end

function InitialCondition(
    a,
    e,
    i,
    ω,
    Ω,
    ν,
    q=DEFAULT_INITIAL_CONDITION_Q,
    ang_vel=DEFAULT_INITIAL_CONDITION_ANG_VEL
)
    return InitialCondition(
        Float64(a),
        Float64(e),
        deg2rad(Float64(i)),
        deg2rad(Float64(ω)),
        deg2rad(Float64(Ω)),
        deg2rad(Float64(ν)),
        _initial_attitude(q),
        _initial_rate(ang_vel),
        Val(:radians)
    )
end

"""
    InitialCondition(; a=0.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, ν=0.0, q=(0, 0, 0, 1), ang_vel=(0, 0, 0))
    InitialCondition(; ra, rp, i=0.0, ω=0.0, Ω=0.0, ν=180.0, q=(0, 0, 0, 1), ang_vel=(0, 0, 0))

Keplerian elements: the semimajor axis `a` (m) and eccentricity `e`, or the
apoapsis and periapsis radii `ra` and `rp` (m); the inclination `i`, argument of
periapsis `ω`, RAAN `Ω` and true anomaly `ν`, in degrees. Pass the mean anomaly
`M` (degrees, `0 <= e < 1`) in place of `ν` to place the spacecraft by time
since periapsis; it is converted to the true anomaly once, here. Without either,
`ν` is 0 (periapsis) for `a` and `e`, and 180 (apoapsis) for `ra` and `rp`. Any
real numbers are accepted. `q` is the initial attitude quaternion (x, y, z, w),
inertial to body, and `ang_vel` the initial body rate (rad/s).
"""
function InitialCondition(;
    a::Real=0.0,
    e::Real=0.0,
    ra::Union{Nothing, Real}=nothing,
    rp::Union{Nothing, Real}=nothing,
    i::Real=0.0,
    ω::Real=0.0,
    Ω::Real=0.0,
    ν::Union{Nothing, Real}=nothing,
    M::Union{Nothing, Real}=nothing,
    q=DEFAULT_INITIAL_CONDITION_Q,
    ang_vel=DEFAULT_INITIAL_CONDITION_ANG_VEL
)
    if (ra === nothing) != (rp === nothing)
        throw(ArgumentError("InitialCondition keyword construction requires both ra and rp when using apsides inputs."))
    end

    if ra !== nothing
        ra, rp = Float64(ra), Float64(rp)
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        ν_eff = _true_anomaly_deg(ν, M, e, 180.0)
        return InitialCondition(a, e, i, ω, Ω, ν_eff, q, ang_vel)
    end

    ν_eff = _true_anomaly_deg(ν, M, Float64(e), 0.0)
    return InitialCondition(a, e, i, ω, Ω, ν_eff, q, ang_vel)
end

# The true anomaly (degrees) from the keywords: `ν`, the mean anomaly `M`
# (degrees) converted on an ellipse of eccentricity `e`, or `default`.
function _true_anomaly_deg(ν, M, e::Float64, default::Float64)::Float64
    ν === nothing || M === nothing ||
        throw(ArgumentError("InitialCondition takes the true anomaly ν or the mean anomaly M, not both."))
    M === nothing && return ν === nothing ? default : Float64(ν)
    0.0 <= e < 1.0 || throw(ArgumentError(
        "InitialCondition: the mean anomaly M needs an elliptic orbit, 0 <= e < 1; got e = $e."))
    e == 0.0 && return Float64(M)   # on a circle the mean and true anomalies are equal
    return rad2deg(OrbitalGeometry._true_anomaly_from_mean(deg2rad(Float64(M)), e))
end

@inline function _initial_condition_lpi(planet, L_PI, initial_time, ephemerides_model)::SMatrix{3, 3, Float64, 9}
    if L_PI !== nothing
        return SMatrix{3, 3, Float64}(L_PI)
    end

    if initial_time !== nothing
        model = ephemerides_model === nothing ? SpiceEphemeridesModel() : ephemerides_model
        return ephemerides_frame_lpi(model, planet, ephemerides_epoch_seconds(model, initial_time))
    end

    if hasproperty(planet, :L_PI)
        planet_lpi = SMatrix{3, 3, Float64}(getproperty(planet, :L_PI))
        sum(abs2, planet_lpi) > 0.0 && return planet_lpi
    end

    return I3
end

@inline function _initial_condition_apsis_direction_ii(
    i::Float64,
    ω::Float64,
    Ω::Float64,
    ν::Float64
)::SVector{3, Float64}
    q = SVector{3, Float64}(cos(ν), sin(ν), 0.0)
    Q = OrbitalGeometry._perifocal_from_inertial(i, ω, Ω)
    u = Q' * q
    return u / norm(u)
end

@inline _initial_condition_oblate_altitude(radius::Float64, u_pp::SVector{3, Float64}, planet)::Float64 =
    geodetic_altitude(radius, u_pp, planet)

@inline _initial_condition_oblate_surface_radius(u_pp::SVector{3, Float64}, planet)::Float64 =
    ellipsoid_surface_radius(u_pp, planet)

function _initial_condition_radius_for_oblate_altitude(
    target_altitude::Float64,
    u_pp::SVector{3, Float64},
    planet
)::Float64
    target_altitude >= 0.0 ||
        throw(ArgumentError("Oblate InitialCondition altitudes must be nonnegative; got $target_altitude m."))
    return radius_for_geodetic_altitude(target_altitude, u_pp, planet)
end

"""
    InitialCondition(planet; ra, hp, i=0.0, ω=0.0, Ω=0.0, ν=180.0, q=(0, 0, 0, 1), ang_vel=(0, 0, 0),
                     initial_time=nothing, ephemerides_model=nothing, L_PI=nothing)

Construct orbital elements from apoapsis altitude `ra` and periapsis altitude
`hp`, measured above the planet's oblate reference ellipsoid along the apsis
directions. The angles are in degrees; pass the mean anomaly `M` in place of the
true anomaly `ν` to place the spacecraft by time since periapsis. Pass `initial_time` and `ephemerides_model` to compute the
inertial-to-planet-fixed frame at construction time. Pass `L_PI` to use a
specific frame directly; otherwise the constructor uses `planet.L_PI` when it
has been initialized. `q` is the initial attitude quaternion (x, y, z, w),
inertial to body, and `ang_vel` the initial body rate (rad/s), each any vector
or tuple of 4 or 3 numbers.
"""
function InitialCondition(
    planet;
    ra::Union{Nothing, Real}=nothing,
    hp::Union{Nothing, Real}=nothing,
    i::Real=0.0,
    ω::Real=0.0,
    Ω::Real=0.0,
    ν::Union{Nothing, Real}=nothing,
    M::Union{Nothing, Real}=nothing,
    q=DEFAULT_INITIAL_CONDITION_Q,
    ang_vel=DEFAULT_INITIAL_CONDITION_ANG_VEL,
    initial_time=nothing,
    ephemerides_model=nothing,
    L_PI::Union{Nothing, AbstractMatrix}=nothing
)
    (ra !== nothing && hp !== nothing) ||
        throw(ArgumentError("Oblate InitialCondition construction requires both ra and hp altitude inputs."))

    i_rad = deg2rad(Float64(i))
    ω_rad = deg2rad(Float64(ω))
    Ω_rad = deg2rad(Float64(Ω))
    l_pi = _initial_condition_lpi(planet, L_PI, initial_time, ephemerides_model)

    u_apo_pp = SVector{3, Float64}(l_pi * _initial_condition_apsis_direction_ii(i_rad, ω_rad, Ω_rad, Float64(pi)))
    u_peri_pp = SVector{3, Float64}(l_pi * _initial_condition_apsis_direction_ii(i_rad, ω_rad, Ω_rad, 0.0))
    u_apo_pp /= norm(u_apo_pp)
    u_peri_pp /= norm(u_peri_pp)

    ra_radius = _initial_condition_radius_for_oblate_altitude(Float64(ra), u_apo_pp, planet)
    rp_radius = _initial_condition_radius_for_oblate_altitude(Float64(hp), u_peri_pp, planet)
    ra_radius > rp_radius ||
        throw(ArgumentError("Oblate InitialCondition requires apoapsis radius > periapsis radius; got ra=$ra_radius m, rp=$rp_radius m."))
    return InitialCondition(
        ra=ra_radius,
        rp=rp_radius,
        i=Float64(i),
        ω=Float64(ω),
        Ω=Float64(Ω),
        ν=ν === nothing ? nothing : Float64(ν),
        M=M,
        q=q,
        ang_vel=ang_vel
    )
end

struct CartesianInitialCondition <: AbstractInitialCondition
    pos::SVector{3, Float64}     # ECI position (m)
    vel::SVector{3, Float64}     # ECI velocity (m/s)
    q::SVector{4, Float64}       # Initial orientation quaternion (x, y, z, w)
    ang_vel::SVector{3, Float64} # Initial angular velocity (rad/s)
end

function CartesianInitialCondition(
    pos,
    vel;
    q=DEFAULT_INITIAL_CONDITION_Q,
    ang_vel=DEFAULT_INITIAL_CONDITION_ANG_VEL
)
    return CartesianInitialCondition(
        SVector{3, Float64}(pos),
        SVector{3, Float64}(vel),
        _initial_attitude(q),
        _initial_rate(ang_vel)
    )
end

mutable struct Link
    root::Bool # Whether this link is a root link (i.e., the main bus or core body of the spacecraft).
    r::MVector{3, Float64} # Position of COM (Body frame for non-root, inertial frame for root)
    q::MVector{4, Float64} # Orientation (Body frame for non-root, inertial frame for root)
    ṙ::MVector{3, Float64} # velocity (body frame for non-root, inertial frame for root)
    ω::MVector{3, Float64} # angular velocity (body frame for non-root, inertial frame for root)
    dims::MVector{3, Float64} # Size[x,y,z] Box= x=thikness, y=z=width, height
    ref_area::Float64 # Reference area for aerodynamic calculations
    m::Float64 # Mass
    mass::SMatrix{3, 3, Float64} # Mass Matrix
    inertia::SMatrix{3, 3, Float64} # Inertial matrix
    aᵇ::SVector{3, Float64} # Left extent (Body frame)
    bᵇ::SVector{3, Float64} # Right extent (Body frame)
    α::Float64 # Angle of attack, rad
    β::Float64 # Sideslip angle, rad
    θ::Float64 # Flow angle, rad
    reflection_coefficient::Float64 # Reflection coefficient for aerodynamic calculations
    rw_assembly::ReactionWheelAssembly # Reaction wheel assembly
    net_force::MVector{3, Float64} # Net force acting on the link, to be updated at each simulation step
    net_torque::MVector{3, Float64} # Net torque acting on the link, to be updated at each simulation step
    attitude_control_rate::Float64 # Rate at which the attitude control function is called, in seconds
    SRP_facets::Vector{Facet}
    thrusters::Vector{Thruster}
    magnets::Vector{Magnet} # List of magnetic dipoles attached to the link
    cop_offset_b::MVector{3, Float64} # Aerodynamic center-of-pressure offset from the link COM (link frame, m); lever arm for the aero wrench torque, zeros = no intrinsic aero torque
    magnetorquers::Vector{Magnetorquer} # Torque rods mounted on the link (see add_actuator!)
    laser_terminals::Vector{LaserTerminal} # Laser terminals mounted on the link (see add_actuator!)
    hinges::Vector{PanelHinge} # The hinge that turns this (non-root) link, if one is mounted (see add_actuator!)
    gimbals::Vector{ThrusterGimbal} # Gimbal axes turning the thrust of thrusters on this link (see add_actuator!)
    # Actuators of user-defined kinds mounted on the link (see add_actuator!). Built-in
    # kinds keep their own concretely typed lists above; this one is evaluated through
    # a function barrier per actuator.
    custom_actuators::Vector{AbstractActuator}
    # Power-system parts mounted on the link (see PowerSystems.add_power!): solar
    # arrays, the battery (root link, at most one) and constant loads.
    solar_arrays::Vector{SolarArray}
    batteries::Vector{Battery}
    power_loads::Vector{PowerLoad}
    # Actuator faults of the spacecraft (kept on its root link; see ActuatorFaults.add_fault!).
    faults::FaultSet
    # Sensors mounted on the link (see SensorMounting.add_sensor!); the link's frame is theirs.
    sensors::Vector{AbstractSensor}

    function Link(; root=false,
        r=MVector{3, Float64}(0, 0, 0),
        q=MVector{4, Float64}(0, 0, 0, 1),
        ṙ=MVector{3, Float64}(0, 0, 0),
        ω=MVector{3, Float64}(0, 0, 0),
        dims=MVector{3, Float64}(0.5, 0.5, 0.1),
        ref_area=1.0,
        m=3.0,
        mass=SMatrix{3, 3, Float64}(m * I3),
        inertia=SMatrix{3, 3, Float64}(1 / 12 * m * diagm([dims[2]^2 + dims[3]^2; dims[1]^2 + dims[3]^2; dims[1]^2 + dims[2]^2])),
        a=SVector{3, Float64}(-0.5 * dims[1], 0, 0),
        b=SVector{3, Float64}(0.5 * dims[1], 0, 0),
        α=pi / 2.0,
        β=0.0,
        θ=0.0,
        reflection_coefficient=1.0,
        max_torque=0.25,
        max_h=70.0,
        rw=Float64[],
        J_rw=zeros(Float64, 3, 0),
        rw_τ=MVector{3, Float64}(zeros(3)),
        net_force=MVector{3, Float64}(zeros(3)),
        net_torque=MVector{3, Float64}(zeros(3)),
        attitude_control_rate=0.1,
        SRP_facets=Facet[],
        J_thruster=nothing,   # accepted for compatibility; ignored (see allocation_matrix)
        thrusters=Thruster[],
        magnets=Magnet[],
        cop_offset_b=MVector{3, Float64}(0, 0, 0),
        magnetorquers=Magnetorquer[],
        laser_terminals=LaserTerminal[],
        hinges=PanelHinge[],
        gimbals=ThrusterGimbal[],
        custom_actuators=AbstractActuator[],
        solar_arrays=SolarArray[],
        batteries=Battery[],
        power_loads=PowerLoad[],
        faults=FaultSet(),
        sensors=AbstractSensor[])

        n_rw = size(J_rw, 2)
        rw_assembly = ReactionWheelAssembly(
            n_wheels=n_rw,
            J_rw=Matrix{Float64}(J_rw),
            max_wheel_torque=max_torque,
            max_wheel_h=max_h,
            h_wheels=zeros(Float64, n_rw),
            h_dot_wheels=zeros(Float64, n_rw),
            tau_body_net=MVector{3, Float64}(zeros(3))      # tau_body_net
        )
        new(root, r, q, ṙ, ω, dims, ref_area, m, mass, inertia, a, b, α, β, θ, reflection_coefficient, rw_assembly, net_force, net_torque, attitude_control_rate, SRP_facets, thrusters, magnets, cop_offset_b, Vector{Magnetorquer}(magnetorquers), Vector{LaserTerminal}(laser_terminals), Vector{PanelHinge}(hinges), Vector{ThrusterGimbal}(gimbals), Vector{AbstractActuator}(custom_actuators),
            Vector{SolarArray}(solar_arrays), Vector{Battery}(batteries), Vector{PowerLoad}(power_loads), faults,
            Vector{AbstractSensor}(sensors))
    end
end


"""Joint types supported by `Joint(...; joint_type=...)`: `:fixed` (default, rigid), `:hinge` (one rotation about `axis`), `:slide` (one translation along `axis`), `:ball` (three rotations, quaternion coordinates)."""
const JOINT_TYPES = (:fixed, :hinge, :slide, :ball)

const _JOINT_Q_IDENTITY = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)

# Validate and normalize the joint-type parameters. Returns (type, axis, stiffness, damping, rest).
function _joint_spec(joint_type, axis, stiffness, damping, rest, initial_q, initial_qd)
    joint_type isa Symbol && joint_type in JOINT_TYPES ||
        throw(ArgumentError("Joint joint_type must be one of $(JOINT_TYPES), got $(repr(joint_type))."))
    has_axis = joint_type === :hinge || joint_type === :slide
    if has_axis && axis === nothing
        throw(ArgumentError("Joint of type $(repr(joint_type)) requires an `axis` (unit vector in the parent link frame)."))
    end
    ax = axis === nothing ? SVector{3, Float64}(0.0, 0.0, 1.0) : SVector{3, Float64}(axis)
    nax = norm(ax)
    (isfinite(nax) && nax > 0.0) || throw(ArgumentError("Joint axis must have finite, nonzero norm, got $(ax)."))
    ax = ax / nax

    gain(name, v) = begin
        if joint_type === :ball
            if v === nothing
                zeros(SMatrix{3, 3, Float64, 9})
            elseif v isa Real
                (isfinite(v) && v >= 0) || throw(ArgumentError("Joint $name must be finite and >= 0, got $v."))
                SMatrix{3, 3, Float64, 9}(Float64(v) * I)
            else
                m = SMatrix{3, 3, Float64, 9}(v)
                all(isfinite, m) || throw(ArgumentError("Joint $name must be finite."))
                s = (m + m') / 2
                lam = eigvals(Symmetric(Matrix(s)))
                minimum(lam) >= -1.0e-12 * max(1.0, maximum(abs, lam)) ||
                    throw(ArgumentError("Joint $name matrix must be positive semidefinite."))
                m
            end
        else
            v === nothing && return 0.0
            v isa Real || throw(ArgumentError("Joint $name must be a scalar for joint type $(repr(joint_type))."))
            (isfinite(v) && v >= 0) || throw(ArgumentError("Joint $name must be finite and >= 0, got $v."))
            Float64(v)
        end
    end
    k = gain("stiffness", stiffness)
    c = gain("damping", damping)

    r = if joint_type === :ball
        if rest === nothing
            _JOINT_Q_IDENTITY
        else
            qr = SVector{4, Float64}(rest)
            nq = norm(qr)
            (isfinite(nq) && nq > 0.0) || throw(ArgumentError("Joint rest quaternion must be finite and nonzero."))
            qr / nq
        end
    else
        rest === nothing && (rest = 0.0)
        rest isa Real && isfinite(rest) || throw(ArgumentError("Joint rest must be a finite scalar for joint type $(repr(joint_type))."))
        Float64(rest)
    end

    q0 = if joint_type === :ball
        if initial_q === nothing
            _JOINT_Q_IDENTITY
        else
            qi = SVector{4, Float64}(initial_q)
            nq = norm(qi)
            (isfinite(nq) && nq > 0.0) || throw(ArgumentError("Joint initial_q quaternion must be finite and nonzero."))
            qi / nq
        end
    else
        initial_q === nothing && (initial_q = 0.0)
        initial_q isa Real && isfinite(initial_q) ||
            throw(ArgumentError("Joint initial_q must be a finite scalar for joint type $(repr(joint_type))."))
        Float64(initial_q)
    end
    qd0 = if joint_type === :ball
        initial_qd === nothing ? zero(SVector{3, Float64}) : (all(isfinite, initial_qd) ? SVector{3, Float64}(initial_qd) :
            throw(ArgumentError("Joint initial_qd must be finite.")))
    else
        initial_qd === nothing && (initial_qd = 0.0)
        initial_qd isa Real && isfinite(initial_qd) ||
            throw(ArgumentError("Joint initial_qd must be a finite scalar for joint type $(repr(joint_type))."))
        Float64(initial_qd)
    end
    return joint_type, ax, k, c, r, q0, qd0
end

"""
    Joint(link1, link2; p1, p2, joint_type=:fixed, axis, stiffness, damping, rest, ...)
    Joint(link1, p1ᵇ, link2, p2ᵇ[, Kx, Kt, Cx, Ct]; joint_type=:fixed, ...)
    Joint(; link1, link2, p1, p2, joint_type=:fixed, ...)

Connection between a parent `link1` and a child `link2`. `p1ᵇ` is the joint point in the
`link1` body frame and `p2ᵇ` the same point in the `link2` body frame. Joint coordinate 0 is the
configured geometry (the links' configured `r`/`q`).

Joint-type keywords used by `build_articulated_tree` and the simulation engine:

- `joint_type`: `:fixed` (default, rigid, merged into the parent body), `:hinge` (rotation about
  `axis`), `:slide` (translation along `axis`) or `:ball` (three rotations).
- `axis`: unit vector in the parent (`link1`) frame; required for `:hinge`/`:slide`, normalized on
  construction.
- `stiffness`, `damping`: joint-space gains, `>= 0`. Scalars (N m/rad or N/m, N m s/rad or N s/m)
  for hinge/slide; a scalar or a 3x3 positive-semidefinite matrix for `:ball`.
- `rest`: joint coordinate at which the spring is unloaded. Scalar for hinge/slide (default 0),
  scalar-last quaternion for `:ball` (default identity).
- `initial_q`, `initial_qd`: initial joint coordinate and rate used by the simulation engine
  (scalars; for `:ball` a scalar-last quaternion and a 3-vector of parent-link angular velocity).
  Default 0 / identity and 0.
- `options`: reserved for future motor and limit settings; unused.

The legacy `Kx/Kt/Cx/Ct` fields stay unused by the articulated-body library.

Ball-joint coordinates, rates and joint-space gains use the parent `link1` frame, even when
`link1` is merged into another body through fixed joints. Output `joint_q` and `joint_qd`
retain that frame; merging does not re-express the public state in the composite body's frame.
"""
mutable struct Joint
    link1::Link
    link2::Link
    p1ᵇ::SVector{3, Float64} # Position of joint in body frame of link1
    p2ᵇ::SVector{3, Float64} # Position of joint in body frame of link2
    Kx::SMatrix{3, 3} # Stiffness Matrix
    Kt::SMatrix{3, 3} # Stiffness Matrix
    Cx::SMatrix{3, 3} # Damping Matrix
    Ct::SMatrix{3, 3} # Damping Matrix
    translational_displacement::SVector{3, Float64} # Translational displacement vector
    rotational_displacement::SVector{4, Float64} # Rotational displacement vector (quaternion)
    joint_type::Symbol # :fixed, :hinge, :slide or :ball
    axis::SVector{3, Float64} # Unit axis in the link1 frame (hinge/slide)
    stiffness::Union{Float64, SMatrix{3, 3, Float64, 9}} # Joint-space stiffness (scalar, or 3x3 for :ball)
    damping::Union{Float64, SMatrix{3, 3, Float64, 9}} # Joint-space damping (scalar, or 3x3 for :ball)
    rest::Union{Float64, SVector{4, Float64}} # Spring rest coordinate (scalar, or quaternion for :ball)
    options::Union{Nothing, NamedTuple} # Reserved for motors and limits; unused
    initial_q::Union{Float64, SVector{4, Float64}} # Initial joint coordinate (scalar, or quaternion for :ball)
    initial_qd::Union{Float64, SVector{3, Float64}} # Initial joint rate (scalar, or 3-vector for :ball)

    function Joint(link1::Link, p1ᵇ::SVector{3, Float64},
        link2::Link, p2ᵇ::SVector{3, Float64},
        Kx=SMatrix{3, 3, Float64}(0.0I),
        Kt=SMatrix{3, 3, Float64}(0.0I),
        Cx=zeros(SMatrix{3, 3, Float64}),
        Ct=zeros(SMatrix{3, 3, Float64});
        joint_type::Symbol=:fixed, axis=nothing, stiffness=nothing, damping=nothing, rest=nothing, options=nothing, initial_q=nothing, initial_qd=nothing)

        jt, ax, k, c, r, iq, iqd = _joint_spec(joint_type, axis, stiffness, damping, rest, initial_q, initial_qd)
        new(link1, link2, p1ᵇ, p2ᵇ, Kx, Kt, Cx, Ct,
            SVector{3, Float64}(0.0, 0.0, 0.0), SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
            jt, ax, k, c, r, options, iq, iqd)
    end

    function Joint(link1::Link, link2::Link; p1=link1.bᵇ, 
                                 p2=link2.aᵇ, 
                                 Kx=SMatrix{3,3, Float64}(0.0I), 
                                 Kt=SMatrix{3,3, Float64}(0.0I), 
                                 Cx=zeros(SMatrix{3,3, Float64}),
                                 Ct=zeros(SMatrix{3,3, Float64}),
                                 translational_displacement=SVector{3, Float64}(0.0, 0.0, 0.0),
                                 rotational_displacement=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
                                 joint_type::Symbol=:fixed, axis=nothing, stiffness=nothing, damping=nothing, rest=nothing, options=nothing, initial_q=nothing, initial_qd=nothing)
        
        jt, ax, k, c, r, iq, iqd = _joint_spec(joint_type, axis, stiffness, damping, rest, initial_q, initial_qd)
        new(link1, link2, p1, p2, Kx, Kt, Cx, Ct, 
            translational_displacement, rotational_displacement, jt, ax, k, c, r, options, iq, iqd)
    end

    function Joint(;link1=Link(), link2=Link(), p1=link1.bᵇ, 
        p2=link2.aᵇ, 
        Kx=SMatrix{3,3, Float64}(1.0I), 
        Kt=SMatrix{3,3, Float64}(1.0I), 
        Cx=zeros(SMatrix{3,3, Float64}),
        Ct=zeros(SMatrix{3,3, Float64}),
        translational_displacement=SVector{3, Float64}(0.0, 0.0, 0.0),
        rotational_displacement=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
        joint_type::Symbol=:fixed, axis=nothing, stiffness=nothing, damping=nothing, rest=nothing, options=nothing, initial_q=nothing, initial_qd=nothing)

        jt, ax, k, c, r, iq, iqd = _joint_spec(joint_type, axis, stiffness, damping, rest, initial_q, initial_qd)
        new(link1, link2, p1, p2, Kx, Kt, Cx, Ct, 
            translational_displacement, rotational_displacement, jt, ax, k, c, r, options, iq, iqd)
    end

    function Joint(joint::Joint)
        new(joint.link1, joint.link2, joint.p1ᵇ, joint.p2ᵇ, 
            joint.Kx, joint.Kt, joint.Cx, joint.Ct,
            joint.translational_displacement, joint.rotational_displacement,
            joint.joint_type, joint.axis, joint.stiffness, joint.damping, joint.rest, joint.options,
            joint.initial_q, joint.initial_qd)
    end
end

"""
    GuidanceModel(; guidance_effectors=(), guidance_rates=Float64[])

Tuple of guidance effectors and the rate (s) at which each is called. `GuidanceModel()` is the empty model (no effectors).
"""
@kwdef struct GuidanceModel{T_Effectors<:Tuple}
    guidance_effectors::T_Effectors = () # Tuple of guidance effector models (maneuver planning, etc.)
    guidance_rates::Vector{Float64} = Float64[] # Rates at which to call each guidance effector, in seconds
    function GuidanceModel(guidance_effectors::T_Effectors, guidance_rates::Vector{Float64}) where {T_Effectors<:Tuple}
        n_effectors = length(guidance_effectors)
        if length(guidance_rates) != n_effectors
            throw(ArgumentError("GuidanceModel guidance_rates length ($(length(guidance_rates))) must match guidance_effectors length ($n_effectors)."))
        end
        @inbounds for (idx, rate) in pairs(guidance_rates)
            if !isfinite(rate) || rate <= 0.0
                throw(ArgumentError("GuidanceModel guidance_rate at index $idx must be finite and > 0.0, got $rate."))
            end
        end
        new{T_Effectors}(guidance_effectors, guidance_rates)
    end
end

"""
    NavigationModel(; navigation_effectors=(), navigation_rates=Float64[])

Tuple of navigation effectors and the rate (s) at which each is called. `NavigationModel()` is the empty model (no effectors).
"""
@kwdef struct NavigationModel{T_Effectors<:Tuple}
    navigation_effectors::T_Effectors = () # Tuple of navigation effector models (sensors, etc.)
    navigation_rates::Vector{Float64} = Float64[] # Rates at which to call each navigation effector, in seconds
    function NavigationModel(navigation_effectors::T_Effectors, navigation_rates::Vector{Float64}) where {T_Effectors<:Tuple}
        n_effectors = length(navigation_effectors)
        if length(navigation_rates) != n_effectors
            throw(ArgumentError("NavigationModel navigation_rates length ($(length(navigation_rates))) must match navigation_effectors length ($n_effectors)."))
        end
        @inbounds for (idx, rate) in pairs(navigation_rates)
            if !isfinite(rate) || rate <= 0.0
                throw(ArgumentError("NavigationModel navigation_rate at index $idx must be finite and > 0.0, got $rate."))
            end
        end
        new{T_Effectors}(navigation_effectors, navigation_rates)
    end
end

"""
    ControlModel(; control_effectors=(), control_rates=Float64[])

Tuple of control effectors and the rate (s) at which each is called. `ControlModel()` is the empty model (no effectors).
"""
@kwdef struct ControlModel{T_Effectors<:Tuple}
    control_effectors::T_Effectors = () # Tuple of control effector models (reaction wheels, thrusters, etc.)
    control_rates::Vector{Float64} = Float64[] # Control rates for each effector, in seconds

     function ControlModel(control_effectors::T_Effectors, control_rates::Vector{Float64}) where {T_Effectors<:Tuple}
        n_effectors = length(control_effectors)
        if length(control_rates) != n_effectors
            throw(ArgumentError("ControlModel control_rates length ($(length(control_rates))) must match control_effectors length ($n_effectors)."))
        end
        @inbounds for (idx, rate) in pairs(control_rates)
            if !isfinite(rate) || rate <= 0.0
                throw(ArgumentError("ControlModel control_rate at index $idx must be finite and > 0.0, got $rate."))
            end
        end
        new{T_Effectors}(control_effectors, control_rates)
    end
end

# Compliant multibody models mounted on links (needs `Link`, needed by `SpacecraftModel`).
include(joinpath(@__DIR__, "..", "..", "dynamics", "multibody", "compliant", "compliant_attachment.jl"))

"""
    ConstellationSlot(plane, slot)

Where a spacecraft sits in a constellation: its orbital plane and its slot in
that plane, both numbered from 1. `build_constellation` and
`walker_constellation` set it on each spacecraft they build, as
`spacecraft.constellation_slot`; any other spacecraft has `nothing` there.
"""
struct ConstellationSlot
    plane::Int
    slot::Int

    function ConstellationSlot(plane::Integer, slot::Integer)
        plane >= 1 && slot >= 1 || throw(ArgumentError(
            "A constellation slot's plane and slot are numbered from 1; got plane $plane, slot $slot."))
        return new(Int(plane), Int(slot))
    end
end

mutable struct SpacecraftModel
    joints::Vector{Joint} # List of joints
    links::Vector{Link} # List of links (bodies)
    root::Link # Root link (main bus or core body)
    instant_actuation::Bool # Whether control inputs (e.g., solar panel angles) are applied instantly
    dry_mass::Float64 # Dry mass of the spacecraft
    prop_mass::Float64 # Fuel mass available for maneuvers
    inertia_tensor::SMatrix{3, 3, Float64} # Inertia tensor of the spacecraft in the body frame
    initial_condition::AbstractInitialCondition # Initial conditions for the simulation (orbit, attitude, etc.)
    id::Int64 # Unique identifier for the spacecraft (useful for multi-spacecraft simulations)
    # Optional per-spacecraft GNC. Empty by default. `flatten_spacecraft_gnc` folds these into the
    # configuration-level GNC tuples once at run setup, binding each effector to this spacecraft's index.
    guidance::GuidanceModel
    navigation::NavigationModel
    control::ControlModel
    # Compliant multibody models mounted on links (compliant panels, flexible appendages). Empty by default;
    # their bodies are not links and not part of `dry_mass`/`mass`.
    attachments::Vector{CompliantAttachment}
    n_terminal::Int64 # Number of interlink terminals
    battery_energy_index::Float64
    tempurature_index::Float64
    constellation_slot::Union{Nothing, ConstellationSlot} # Plane and slot in a constellation (see build_constellation), or nothing
end

# The full positional form before the constellation slot was added: no slot.
SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                initial_condition, id, guidance, navigation, control, attachments,
                n_terminal, battery_energy_index, tempurature_index) =
    SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                    initial_condition, id, guidance, navigation, control, attachments,
                    n_terminal, battery_energy_index, tempurature_index, nothing)

_empty_guidance() = GuidanceModel((), Float64[])
_empty_navigation() = NavigationModel((), Float64[])
_empty_control() = ControlModel((), Float64[])

# The spacecraft once stored its thruster and reaction-wheel counts, which
# nothing read. The mounted hardware is now the record (`mounted_thruster_count`,
# `mounted_wheel_count`, `root.rw_assembly.n_wheels`). The positional forms that
# take the two counts after the inertia tensor stay for compatibility, and
# ignore them.
SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                n_reaction_wheels::Integer, n_thrusters::Integer, initial_condition, id,
                guidance::GuidanceModel, navigation::NavigationModel, control::ControlModel,
                attachments::Vector{CompliantAttachment}, n_terminal, battery_energy_index, tempurature_index) =
    SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                    initial_condition, id, guidance, navigation, control, attachments,
                    n_terminal, battery_energy_index, tempurature_index)

SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                n_reaction_wheels::Integer, n_thrusters::Integer, initial_condition, id,
                guidance::GuidanceModel, navigation::NavigationModel, control::ControlModel,
                attachments::Vector{CompliantAttachment}) =
    SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                    initial_condition, id, guidance, navigation, control, attachments)

SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                n_reaction_wheels::Integer, n_thrusters::Integer, initial_condition, id,
                guidance::GuidanceModel, navigation::NavigationModel, control::ControlModel) =
    SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                    initial_condition, id, guidance, navigation, control)

SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                n_reaction_wheels::Integer, n_thrusters::Integer, initial_condition, id) =
    SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                    initial_condition, id)

# Positional form before interlink metadata was added.
SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                initial_condition, id,
                guidance::GuidanceModel, navigation::NavigationModel, control::ControlModel,
                attachments::Vector{CompliantAttachment}) =
    SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                    initial_condition, id,
                    guidance, navigation, control, attachments, 1, 100.0, 100.0)

# Positional form with GNC but no attachments.
SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                initial_condition, id,
                guidance::GuidanceModel, navigation::NavigationModel, control::ControlModel) =
    SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                    initial_condition, id,
                    guidance, navigation, control, CompliantAttachment[])

# Positional form without GNC: declares no per-spacecraft GNC.
SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                initial_condition, id) =
    SpacecraftModel(joints, links, root, instant_actuation, dry_mass, prop_mass, inertia_tensor,
                    initial_condition, id,
                    _empty_guidance(), _empty_navigation(), _empty_control(), CompliantAttachment[])

function SpacecraftModel(; joints::AbstractVector{<:Joint}=Joint[], links::AbstractVector{<:Link}=Link[], root::Link=Link(root=true),
                            instant_actuation::Bool=true,
                            prop_mass::Float64=0.0,
                            inertia_tensor::SMatrix{3,3,Float64}=SMatrix{3, 3, Float64}(zeros(3,3)),
                            n_reaction_wheels::Integer=0,   # accepted for compatibility; ignored
                            n_thrusters::Integer=0,         # accepted for compatibility; ignored
                            initial_condition::AbstractInitialCondition=InitialCondition(),
                            id::Int64=1,
                            guidance::GuidanceModel=_empty_guidance(),
                            navigation::NavigationModel=_empty_navigation(),
                            control::ControlModel=_empty_control(),
                            attachments::AbstractVector{CompliantAttachment}=CompliantAttachment[],
                            n_terminal::Int64=1,
                            battery_energy_index::Real=100.0,
                            tempurature_index::Real=100.0,
                            constellation_slot::Union{Nothing, ConstellationSlot}=nothing)
    n_terminal >= 0 || throw(ArgumentError("n_terminal must be nonnegative."))
    0 <= battery_energy_index <= 100 || throw(ArgumentError("battery_energy_index must be in [0, 100]."))
    0 <= tempurature_index <= 100 || throw(ArgumentError("tempurature_index must be in [0, 100]."))
    joints_vec = Vector{Joint}(joints)
    links_vec = Vector{Link}(links)

    dry_mass = 0.0
    if !any(link -> link === root, links_vec)
        # Keep the root in the link list so dry-mass aggregation sees the full assembly.
        push!(links_vec, root) # Include root in the links list for mass calculation
    end
    for link in links_vec
        dry_mass += link.m
    end
    for (k, a) in pairs(attachments)
        any(l -> l === a.link, links_vec) ||
            throw(ArgumentError("Attachment $k is mounted on a link that is not one of the spacecraft's links (or its root)."))
    end

    return SpacecraftModel(joints_vec, links_vec, root, instant_actuation, dry_mass, prop_mass, inertia_tensor, initial_condition, id,
                           guidance, navigation, control, Vector{CompliantAttachment}(attachments),
                           n_terminal, Float64(battery_energy_index), Float64(tempurature_index),
                           constellation_slot)
end

# Preserve the collection type in this module while giving configuration its own file.
include(joinpath(@__DIR__, "..", "..", "..", "orchestration", "simulation", "config", "constellation_configuration.jl"))

# Constellation layouts: orbital planes, the Walker patterns and the builder.
include(joinpath(@__DIR__, "..", "constellations", "builder.jl"))

end # module SpacecraftModels
