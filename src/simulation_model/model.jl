module PhysicalModel

using StaticArrays
using LinearAlgebra

# Import types defined in our other files
# using ..AbstractTypes
using ..Components

export Link, Joint, SpacecraftModel

const I3 = SMatrix{3, 3, Float64}(diagm(ones(3)))

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
    rw_assembly::ReactionWheelAssembly{N_RW} # Reaction wheel assembly
    net_force::MVector{3, Float64} # Net force acting on the link, to be updated at each simulation step
    net_torque::MVector{3, Float64} # Net torque acting on the link, to be updated at each simulation step
    attitude_control_rate::Float64 # Rate at which the attitude control function is called, in seconds
    ω_wheel_derivatives::MVector{N_RW, Float64} # Angular momentum derivatives of the reaction wheels, to be updated at each simulation step
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
        max_torque=0.25,
        max_h=70.0,
        rw=MVector{N_RW, Float64}(zeros(N_RW)),
        J_rw=MMatrix{3, N_RW, Float64}(zeros(3, N_RW)),
        rw_τ=MVector{3, Float64}(zeros(3)),
        net_force=MVector{3, Float64}(zeros(3)),
        net_torque=MVector{3, Float64}(zeros(3)),
        attitude_control_rate=0.1,
        ω_wheel_derivatives=MVector{N_RW, Float64}(zeros(N_RW)),
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
        println(length(rw))
        new{N_RW}(root, r, q, ṙ, ω, dims, ref_area, m, mass, inertia, a, b, α, β, θ, rw_assembly, net_force, net_torque, attitude_control_rate, ω_wheel_derivatives, SRP_facets, J_thruster, thrusters, magnets)
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
end

# --- CRITICAL REFACTOR 2 ---
# SpacecraftModel is now parametric on T_Effectors, a Tuple type.
# This ensures the `dynamic_effectors` field is type-stable.
mutable struct SpacecraftModel{T_Effectors<:Tuple}
    joints::Vector{Joint} # List of joints
    links::Vector{Link} # List of links (bodies)
    roots::Vector{Link} # Vector of root links (main bus or core bodies)
    instant_actuation::Bool # Whether control inputs (e.g., solar panel angles) are applied instantly
    dry_mass::Float64 # Dry mass of the spacecraft
    prop_mass::Vector{Float64} # Fuel mass available for maneuvers
    inertia_tensors::Vector{SMatrix{3, 3, Float64}} # Inertia tensors of the spacecraft bodies in the body frame
    n_reaction_wheels::Int64 # Number of reaction wheels in the spacecraft model
    n_thrusters::Int64 # Number of thrusters in the spacecraft model
    dynamic_effectors::T_Effectors # List of dynamic effector models (gravity, drag, etc.)

    function SpacecraftModel(; joints=Joint[], links=Link[], roots=Link[],
        instant_actuation=true,
        prop_mass=Float64[],
        inertia_tensors=SMatrix{3, 3, Float64}[],
        n_reaction_wheels=0,
        n_thrusters=0,
        dynamic_effectors::T_Effectors=()) where {T_Effectors<:Tuple} # Set default as empty tuple
        
        # Calculate dry mass from links
        dry_mass = 0.0
        for link in links
            dry_mass += link.m
        end

        new{T_Effectors}(joints, links, roots, instant_actuation, dry_mass, prop_mass, inertia_tensors, n_reaction_wheels, n_thrusters, dynamic_effectors)
    end
end

end # module Model