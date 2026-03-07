module PhysicalModel

using StaticArrays
using LinearAlgebra

# Import types defined in our other files
# using ..AbstractTypes
using ..Components

export Link, Joint, SpacecraftModel, DynamicsModel, InitialCondition, GuidanceModel, NavigationModel, ControlModel

const I3 = SMatrix{3, 3, Float64}(diagm(ones(3)))
const DEFAULT_INITIAL_CONDITION_Q = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
const DEFAULT_INITIAL_CONDITION_ANG_VEL = SVector{3, Float64}(0.0, 0.0, 0.0)

struct InitialCondition
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

# --- CRITICAL REFACTOR 1 ---
# Link is now parametric on N_RW (number of reaction wheels).
# This avoids the type-instability of using `gyro` as a value.
# I also removed the ..._function fields to decouple data from logic.
mutable struct Link{N_RW}
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
    rw_assembly::ReactionWheelAssembly{N_RW} # Reaction wheel assembly
    net_force::MVector{3, Float64} # Net force acting on the link, to be updated at each simulation step
    net_torque::MVector{3, Float64} # Net torque acting on the link, to be updated at each simulation step
    attitude_control_rate::Float64 # Rate at which the attitude control function is called, in seconds
    SRP_facets::Vector{Facet}
    J_thruster::Matrix{Float64} # Thruster Jacobian matrix
    thrusters::Vector{Thruster}
    magnets::Vector{Magnet} # List of magnetic dipoles attached to the link

    function Link{N_RW}(; root=false,
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
        rw=MVector{N_RW, Float64}(zeros(N_RW)),
        J_rw=MMatrix{3, N_RW, Float64}(zeros(3, N_RW)),
        rw_τ=MVector{3, Float64}(zeros(3)),
        net_force=MVector{3, Float64}(zeros(3)),
        net_torque=MVector{3, Float64}(zeros(3)),
        attitude_control_rate=0.1,
        SRP_facets=Facet[],
        J_thruster=Matrix{Float64}(zeros(3, 1)),
        thrusters=Thruster[],
        magnets=Magnet[]) where {N_RW}

        rw_assembly = ReactionWheelAssembly{N_RW}(
            J_rw=J_rw,
            max_wheel_torque=max_torque,
            max_wheel_h=max_h,
            h_wheels=MVector{N_RW, Float64}(zeros(N_RW)), # h_wheels
            h_dot_wheels=MVector{N_RW, Float64}(zeros(N_RW)), # h_dot_wheels
            tau_body_net=MVector{3, Float64}(zeros(3))      # tau_body_net
        )
        new{N_RW}(root, r, q, ṙ, ω, dims, ref_area, m, mass, inertia, a, b, α, β, θ, reflection_coefficient, rw_assembly, net_force, net_torque, attitude_control_rate, SRP_facets, J_thruster, thrusters, magnets)
    end
end

Link(; kwargs...) = Link{0}(; kwargs...)

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

    function Joint(;link1=Link{0}(), link2=Link{0}(), p1=link1.bᵇ, 
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

# --- CRITICAL REFACTOR 2 ---
# SpacecraftModel is now parametric on T_Effectors, a Tuple type.
# This ensures the `dynamic_effectors` field is type-stable.
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
    initial_condition::InitialCondition # Initial conditions for the simulation (orbit, attitude, etc.)
    id::Int64 # Unique identifier for the spacecraft (useful for multi-spacecraft simulations)
    # dynamic_effectors::T_Effectors # List of dynamic effector models (gravity, drag, etc.)

    # function SpacecraftModel(; joints=Joint[], children=Link[], root=Link{0}(),
    #     instant_actuation=true,
    #     prop_mass=0.0,
    #     inertia_tensor=SMatrix{3, 3, Float64}(zeros(3, 3)),
    #     n_reaction_wheels=0,
    #     n_thrusters=0,
    #     initial_condition) # Set default as empty tuple
        
    #     # Calculate dry mass from children
    #     dry_mass = root.m
    #     for child in children
    #         dry_mass += child.m
    #     end

    #     new(joints, children, root, instant_actuation, dry_mass, prop_mass, inertia_tensor, n_reaction_wheels, n_thrusters)
    # end
end

function SpacecraftModel(; joints::AbstractVector{<:Joint}=Joint[], links::AbstractVector{<:Link}=Link[], root::Link=Link{0}(root=true),
                            instant_actuation::Bool=true,
                            prop_mass::Float64=0.0,
                            inertia_tensor::SMatrix{3,3,Float64}=SMatrix{3, 3, Float64}(zeros(3,3)),
                            n_reaction_wheels::Int64=0,
                            n_thrusters::Int64=0,
                            initial_condition::InitialCondition=InitialCondition(),
                            id::Int64=1)
    joints_vec = Vector{Joint}(joints)
    links_vec = Vector{Link}(links)

    # Calculate dry mass from links
    dry_mass = 0.0
    if !any(link -> link === root, links_vec)
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
end # module Model
