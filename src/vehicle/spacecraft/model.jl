module SpacecraftModels

using StaticArrays
using LinearAlgebra

using ..Components
using ..EphemeridesModels: SpiceEphemeridesModel, ephemerides_time_seconds, planet_frame_lpi

export Link, Joint, SpacecraftModel, DynamicsModel, InitialCondition, CartesianInitialCondition, AbstractInitialCondition, GuidanceModel, NavigationModel, ControlModel

const I3 = SMatrix{3, 3, Float64}(diagm(ones(3)))
const DEFAULT_INITIAL_CONDITION_Q = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
const DEFAULT_INITIAL_CONDITION_ANG_VEL = SVector{3, Float64}(0.0, 0.0, 0.0)

abstract type AbstractInitialCondition end

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
        q,
        ang_vel,
        Val(:radians)
    )
end

function InitialCondition(;
    a::Float64=0.0,
    e::Float64=0.0,
    ra::Union{Nothing, Float64}=nothing,
    rp::Union{Nothing, Float64}=nothing,
    i::Float64=0.0,
    ω::Float64=0.0,
    Ω::Float64=0.0,
    ν::Union{Nothing, Float64}=nothing,
    q::SVector{4, Float64}=DEFAULT_INITIAL_CONDITION_Q,
    ang_vel::SVector{3, Float64}=DEFAULT_INITIAL_CONDITION_ANG_VEL
)
    if (ra === nothing) != (rp === nothing)
        throw(ArgumentError("InitialCondition keyword construction requires both ra and rp when using apsides inputs."))
    end

    if ra !== nothing
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        ν_eff = ν === nothing ? 180.0 : ν
        return InitialCondition(a, e, i, ω, Ω, ν_eff, q, ang_vel)
    end

    ν_eff = ν === nothing ? 0.0 : ν
    return InitialCondition(a, e, i, ω, Ω, ν_eff, q, ang_vel)
end

@inline function _initial_condition_lpi(planet, L_PI, initial_time, ephemerides_model)::SMatrix{3, 3, Float64, 9}
    if L_PI !== nothing
        return SMatrix{3, 3, Float64}(L_PI)
    end

    if initial_time !== nothing
        model = ephemerides_model === nothing ? SpiceEphemeridesModel() : ephemerides_model
        return planet_frame_lpi(planet, ephemerides_time_seconds(initial_time, model), model)
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
    Q = @SMatrix [
        -sin(Ω)*cos(i)*sin(ω)+cos(Ω)*cos(ω)  cos(Ω)*cos(i)*sin(ω)+sin(Ω)*cos(ω)  sin(i)*sin(ω)
        -sin(Ω)*cos(i)*cos(ω)-cos(Ω)*sin(ω)  cos(Ω)*cos(i)*cos(ω)-sin(Ω)*sin(ω)  sin(i)*cos(ω)
         sin(Ω)*sin(i)                       -cos(Ω)*sin(i)                       cos(i)
    ]
    u = Q' * q
    return u / norm(u)
end

@inline function _initial_condition_oblate_altitude(radius::Float64, u_pp::SVector{3, Float64}, planet)::Float64
    x = radius * u_pp[1]
    y = radius * u_pp[2]
    z = radius * u_pp[3]

    f = (planet.Rp_e - planet.Rp_p) / planet.Rp_e
    e2 = 1.0 - (1.0 - f)^2
    ep2 = e2 / (1.0 - e2)
    p_xy = sqrt(x^2 + y^2)

    θ = atan(z * planet.Rp_e, p_xy * planet.Rp_p)
    lat = atan(z + ep2 * planet.Rp_p * sin(θ)^3, p_xy - e2 * planet.Rp_e * cos(θ)^3)
    N = planet.Rp_e / sqrt(1.0 - e2 * sin(lat)^2)
    return p_xy * cos(lat) + (z + e2 * N * sin(lat)^2) * sin(lat) - N
end

@inline function _initial_condition_oblate_surface_radius(u_pp::SVector{3, Float64}, planet)::Float64
    return inv(sqrt((u_pp[1]^2 + u_pp[2]^2) / planet.Rp_e^2 + u_pp[3]^2 / planet.Rp_p^2))
end

function _initial_condition_radius_for_oblate_altitude(
    target_altitude::Float64,
    u_pp::SVector{3, Float64},
    planet
)::Float64
    target_altitude >= 0.0 ||
        throw(ArgumentError("Oblate InitialCondition altitudes must be nonnegative; got $target_altitude m."))

    lo = _initial_condition_oblate_surface_radius(u_pp, planet)
    hi = lo + target_altitude + abs(planet.Rp_e - planet.Rp_p) + 1.0
    while _initial_condition_oblate_altitude(hi, u_pp, planet) < target_altitude
        hi += max(target_altitude, abs(planet.Rp_e - planet.Rp_p), 1.0)
    end

    for _ in 1:80
        mid = 0.5 * (lo + hi)
        if _initial_condition_oblate_altitude(mid, u_pp, planet) < target_altitude
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

"""
    InitialCondition(planet; ra, hp, i=0.0, ω=0.0, Ω=0.0, ν=180.0, initial_time=nothing, ephemerides_model=nothing, L_PI=nothing)

Construct orbital elements from apoapsis altitude `ra` and periapsis altitude
`hp`, measured above the planet's oblate reference ellipsoid along the apsis
directions. Pass `initial_time` and `ephemerides_model` to compute the
inertial-to-planet-fixed frame at construction time. Pass `L_PI` to use a
specific frame directly; otherwise the constructor uses `planet.L_PI` when it
has been initialized.
"""
function InitialCondition(
    planet;
    ra::Union{Nothing, Real}=nothing,
    hp::Union{Nothing, Real}=nothing,
    i::Real=0.0,
    ω::Real=0.0,
    Ω::Real=0.0,
    ν::Union{Nothing, Real}=nothing,
    q::SVector{4, Float64}=DEFAULT_INITIAL_CONDITION_Q,
    ang_vel::SVector{3, Float64}=DEFAULT_INITIAL_CONDITION_ANG_VEL,
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
    ν_eff = ν === nothing ? 180.0 : Float64(ν)

    return InitialCondition(
        ra=ra_radius,
        rp=rp_radius,
        i=Float64(i),
        ω=Float64(ω),
        Ω=Float64(Ω),
        ν=ν_eff,
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
    q::SVector{4, Float64}=DEFAULT_INITIAL_CONDITION_Q,
    ang_vel::SVector{3, Float64}=DEFAULT_INITIAL_CONDITION_ANG_VEL
)
    return CartesianInitialCondition(
        SVector{3, Float64}(pos),
        SVector{3, Float64}(vel),
        q,
        ang_vel
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
    J_thruster::Matrix{Float64} # Thruster Jacobian matrix
    thrusters::Vector{Thruster}
    magnets::Vector{Magnet} # List of magnetic dipoles attached to the link
    cop_offset_b::MVector{3, Float64} # Aerodynamic center-of-pressure offset from the link COM (link frame, m); lever arm for the aero wrench torque, zeros = no intrinsic aero torque

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
        J_thruster=Matrix{Float64}(zeros(3, 1)),
        thrusters=Thruster[],
        magnets=Magnet[],
        cop_offset_b=MVector{3, Float64}(0, 0, 0))

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
        new(root, r, q, ṙ, ω, dims, ref_area, m, mass, inertia, a, b, α, β, θ, reflection_coefficient, rw_assembly, net_force, net_torque, attitude_control_rate, SRP_facets, J_thruster, thrusters, magnets, cop_offset_b)
    end
end


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

    function Joint(link1::Link, p1ᵇ::SVector{3, Float64},
        link2::Link, p2ᵇ::SVector{3, Float64},
        Kx=SMatrix{3, 3, Float64}(0.0I),
        Kt=SMatrix{3, 3, Float64}(0.0I),
        Cx=zeros(SMatrix{3, 3, Float64}),
        Ct=zeros(SMatrix{3, 3, Float64}))

        new(link1, link2, p1ᵇ, p2ᵇ, Kx, Kt, Cx, Ct,
            SVector{3, Float64}(0.0, 0.0, 0.0), SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
    end

    function Joint(link1::Link, link2::Link; p1=link1.bᵇ, 
                                 p2=link2.aᵇ, 
                                 Kx=SMatrix{3,3, Float64}(0.0I), 
                                 Kt=SMatrix{3,3, Float64}(0.0I), 
                                 Cx=zeros(SMatrix{3,3, Float64}),
                                 Ct=zeros(SMatrix{3,3, Float64}),
                                 translational_displacement=SVector{3, Float64}(0.0, 0.0, 0.0),
                                 rotational_displacement=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
        
        new(link1, link2, p1, p2, Kx, Kt, Cx, Ct, 
            translational_displacement, rotational_displacement)
    end

    function Joint(;link1=Link(), link2=Link(), p1=link1.bᵇ, 
        p2=link2.aᵇ, 
        Kx=SMatrix{3,3, Float64}(1.0I), 
        Kt=SMatrix{3,3, Float64}(1.0I), 
        Cx=zeros(SMatrix{3,3, Float64}),
        Ct=zeros(SMatrix{3,3, Float64}),
        translational_displacement=SVector{3, Float64}(0.0, 0.0, 0.0),
        rotational_displacement=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))

        new(link1, link2, p1, p2, Kx, Kt, Cx, Ct, 
            translational_displacement, rotational_displacement)
    end

    function Joint(joint::Joint)
        new(joint.link1, joint.link2, joint.p1ᵇ, joint.p2ᵇ, 
            joint.Kx, joint.Kt, joint.Cx, joint.Ct,
            joint.translational_displacement, joint.rotational_displacement)
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
    n_reaction_wheels::Int64 # Number of reaction wheels in the spacecraft model
    n_thrusters::Int64 # Number of thrusters in the spacecraft model
    initial_condition::AbstractInitialCondition # Initial conditions for the simulation (orbit, attitude, etc.)
    id::Int64 # Unique identifier for the spacecraft (useful for multi-spacecraft simulations)
end

function SpacecraftModel(; joints::AbstractVector{<:Joint}=Joint[], links::AbstractVector{<:Link}=Link[], root::Link=Link(root=true),
                            instant_actuation::Bool=true,
                            prop_mass::Float64=0.0,
                            inertia_tensor::SMatrix{3,3,Float64}=SMatrix{3, 3, Float64}(zeros(3,3)),
                            n_reaction_wheels::Int64=0,
                            n_thrusters::Int64=0,
                            initial_condition::AbstractInitialCondition=InitialCondition(),
                            id::Int64=1)
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

    return SpacecraftModel(joints_vec, links_vec, root, instant_actuation, dry_mass, prop_mass, inertia_tensor, n_reaction_wheels, n_thrusters, initial_condition, id)
end

"""
Struct containing the information needed for the dynamics propagation of all spacecraft

Roots: Vector of root links (main bus or core bodies). This is the vector over which the threads will be parallelized, so each thread will handle the dynamics propagation of one root link and its associated sub-links.
DynamicEffectors: Tuple of dynamic effector models (gravity, drag, etc.) to be applied during the dynamics propagation.
"""
struct DynamicsModel{T_Effectors<:Tuple}
    spacecraft::Vector{SpacecraftModel} # Vector of root links (main bus or core bodies)
    dynamic_effectors::T_Effectors # Tuple of dynamic effector models (gravity, drag, etc.)

    function DynamicsModel(roots::Vector{SpacecraftModel}, dynamic_effectors::T_Effectors) where {T_Effectors<:Tuple}
        new{T_Effectors}(roots, dynamic_effectors)
    end
end

@kwdef struct GuidanceModel{T_Effectors<:Tuple}
    guidance_effectors::T_Effectors # Tuple of guidance effector models (maneuver planning, etc.)
    guidance_rates::Vector{Float64} # Rates at which to call each guidance effector, in seconds
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

@kwdef struct NavigationModel{T_Effectors<:Tuple}
    navigation_effectors::T_Effectors # Tuple of navigation effector models (sensors, etc.)
    navigation_rates::Vector{Float64} # Rates at which to call each navigation effector, in seconds
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

@kwdef struct ControlModel{T_Effectors<:Tuple}
    control_effectors::T_Effectors # Tuple of control effector models (reaction wheels, thrusters, etc.)
    control_rates::Vector{Float64} # Control rates for each effector, in seconds

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
end # module SpacecraftModels
