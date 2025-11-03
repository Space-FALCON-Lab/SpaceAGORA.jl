module Components

using StaticArrays
using LinearAlgebra
export Facet, Thruster, Magnet, ReactionWheelAssembly, create_facet_list

const I3 = SMatrix{3, 3, Float64}(diagm(ones(3)))

@kwdef mutable struct Facet
    area::Float64 = 0.0 # Area of the facet
    attitude::MVector{4, Float64} = @MVector [0.0, 0.0, 0.0, 1.0] # Attitude of the facet relative to the body frame
    normal_vector::MVector{3, Float64} = @MVector [1.0, 0.0, 0.0] # Normal vector in the facet body frame, for flat plate this is in the x-direction
    cp::MVector{3, Float64} = @MVector [0.0, 0.0, 0.0] # Center of pressure with respect to the center of mass of the Link containing the facet, in the Link frame
    ρ::Float64 = 0.0 # Diffuse coefficient
    δ::Float64 = 0.0 # Specular coefficient
    name::String = "" # Identifying name
end

@kwdef mutable struct Thruster
    max_thrust::Float64 = 1.0 # Maximum thrust, N
    location::MVector{3, Float64} = MVector{3, Float64}(zeros(3)) # Location in the link frame, relative to the CoM of the link, m
    direction::MVector{3, Float64} = MVector{3, Float64}(zeros(3)) # Unit vector direction of thrust in the link frame, n/d
    Isp::Float64 = 0.0 # Specific impulse of the thruster, s
    cutoff_frequency::Float64 = 1.0 # Governing parameter for thruster ramp-up/ramp-down, rad/s.
    min_firing_time::Float64 = 0.0 # Minimum time for the thruster to be firing, s.
    level_on::Float64 = 0.75 # Upper threshold for the Schmitt trigger, nd
    level_off::Float64 = 0.25 # Lower threshold for the Schmitt trigger, nd
    κ::Float64 = 0.0 # Thrust factor, nd.
    thrust::Float64 = 0.0 # Current thrust magnitude, to be updated during simulation.
    stop_firing_time::Float64 = 0.0 # Time at which the thruster should stop firing, s.
end

@kwdef mutable struct Magnet
    m::MVector{3, Float64} = MVector{3, Float64}(zeros(3)) # Magnetic dipole moment in the body frame, nT
    location::MVector{3, Float64} = MVector{3, Float64}(zeros(3)) # Location in the link frame, relative to the CoM of the link, m
end

"""
    ReactionWheelAssembly{N}

A struct to hold the parameters and state of a
reaction wheel assembly with N wheels.
"""
@kwdef mutable struct ReactionWheelAssembly{N}
    # --- Parameters (Fixed) ---
    n_wheels::Int = N # Number of reaction wheels in the assembly

    "Jacobian mapping wheel angular velocity to body angular momentum"
    J_rw::SMatrix{3, N, Float64} = SMatrix{3, N, Float64}(eye(3, N)) # Default to identity mapping for 3 wheels
    
    "Maximum torque (Nm) *each* wheel can produce"
    max_wheel_torque::Float64 = 0.1
    
    "Maximum angular momentum (Nms) *each* wheel can store"
    max_wheel_h::Float64 = 0.1

    # --- State (Mutable) ---
    "Angular momentum (h) of each wheel"
    h_wheels::MVector{N, Float64} = MVector{N, Float64}(zeros(N))
    
    "Derivative of wheel angular momentum (h_dot)"
    h_dot_wheels::MVector{N, Float64} = MVector{N, Float64}(zeros(N))
    
    "Net torque (τ) applied to the body *by* the RW assembly"
    tau_body_net::MVector{3, Float64} = MVector{3, Float64}(zeros(3))
end

function create_facet_list(area_list::Vector{Float64}, attitude_list::Vector{SVector{4, Float64}}, normal_vector_list::Vector{SVector{3, Float64}},
    cp_loc_list::Vector{SVector{3, Float64}},
    diffuse_coeffs_list::Vector{Float64}, specular_coeffs_list::Vector{Float64}, facet_names_list::Vector{String})
    list_length = length(area_list)
    lists = Vector[area_list, attitude_list, normal_vector_list, cp_loc_list, specular_coeffs_list, diffuse_coeffs_list]
    @assert all(l -> length(l) == list_length, lists) "Not all lists are the same length"
    facet_vector = Vector{Facet}(undef, list_length)
    for i in eachindex(area_list)
        facet_vector[i] = Facet(area_list[i], attitude_list[i], normal_vector_list[i], cp_loc_list[i], diffuse_coeffs_list[i], specular_coeffs_list[i], facet_names_list[i])
    end
    return facet_vector
end

end # module Components