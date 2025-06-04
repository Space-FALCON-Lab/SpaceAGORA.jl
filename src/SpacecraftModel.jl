using StaticArrays
using LinearAlgebra

const I3 = SMatrix{3, 3, Float64}(diagm(ones(3)))

abstract type Joint end

struct FixedJoint <: Joint end
struct RevoluteJoint <: Joint
    axis::SVector{3, Float64}
end

abstract type RigidBody end

# Define various rigid body types with their properties
# Flat plate for solar panels, boxes for main bus
struct FlatPlate <: RigidBody
    name::String
    mass::Float64
    inertia_tensor::SMatrix{3, 3, Float64} # In principal axes
    dimensions::SVector{2, Float64} # Length/width dimensions
    area::Float64 # Reference area
    COM::SVector{3, Float64} # Center of mass in body frame, i.e., offset from geometric center
end

struct Box <: RigidBody
    name::String
    mass::Float64
    inertia_tensor::SMatrix{3, 3, Float64} # In principal axes
    dimensions::SVector{3, Float64} # Lengths in x, y, z directions
    area::Float64 # Aerodynamic reference area
    COM::SVector{3, Float64} # Center of mass in body frame, i.e., offset from geometric center
end

struct ReactionWheel <: RigidBody
    name::String
    mass::Float64
    inertia_tensor::SMatrix{3, 3, Float64} # In principal axes
    COM::SVector{3, Float64} # Center of mass in body frame, i.e., offset from geometric center
    max_speed::Float64 # Maximum speed of the wheel
end

mutable struct BodyNode
    body::RigidBody
    joint::Joint
    parent::Union{Nothing, Int} # Parent node in the tree structure
    children::Vector{Int} # Children nodes in the tree structure
    transform::SMatrix{4, 4, Float64} # Transformation matrix from parent to this node
end

mutable struct SpacecraftModel
    bodies::Vector{BodyNode} # List of body nodes
    root::Int # Index of the root node in the bodies vector
    q::Vector{Float64} # Generalized coordinates (joint angles, etc.)
    q_dot::Vector{Float64} # Generalized velocities
    joint_indices::Dict{Int, Int} # Body index => joint index in q, q_dot
    instant_actuation::Bool # Whether control inputs (e.g., solar panel angles) are applied instantly
end

# Function to add a body to the spacecraft model
function add_body!(model::SpacecraftModel, 
                   body::RigidBody, 
                   joint::Joint, 
                   parent_index::Union{Nothing, Int},
                   transform::SMatrix{4, 4, Float64})
    """
    Adds a new body to the spacecraft model.
    - `model`: The spacecraft model to which the body will be added.
    - `body`: The rigid body to be added.
    - `joint`: The joint connecting this body to its parent.
    - `parent_index`: The index of the parent body in the model's bodies vector, or `nothing` if this is the root body.
    - `transform`: The transformation matrix from the parent body to this body.
    """

    index = length(model.bodies) + 1
    new_node = BodyNode(body, joint, parent_index, Int[], transform)

    if isnothing(parent_index)
        push!(model.bodies, new_node)
        model.root = length(model.bodies)
    else
        parent_node = model.bodies[parent_index]
        new_node.parent = parent_index
        push!(parent_node.children, index)
        push!(model.bodies, new_node)
    end

    if joint isa RevoluteJoint
        push!(model.q, 0.0) # Initialize joint angle
        push!(model.q_dot, 0.0) # Initialize joint velocity
        model.joint_indices[index] = length(model.q) # Map body index to joint index
    end
end

function get_spacecraft_reference_area(model::SpacecraftModel)
    """
    Calculates the total reference area of the spacecraft model by summing the areas of all bodies,
    excluding actuators (reaction wheels, cmg, etc.) because they are internal.
    """
    total_area = 0.0
    for node in model.bodies
        if node.body isa FlatPlate || node.body isa Box
            # Only consider bodies with defined area
            total_area += node.body.area
        end
    end
    return total_area
end

# Function to create a transformation matrix for a translation
function translation(r::SVector{3, Float64})
    """
    Creates a transformation matrix for a translation by vector `r`.
    - `r`: The translation vector in the body frame.
    """
    T = @SMatrix [
        1.0 0.0 0.0 r[1];
        0.0 1.0 0.0 r[2];
        0.0 0.0 1.0 r[3];
        0.0 0.0 0.0 1.0
    ]
    return T
end