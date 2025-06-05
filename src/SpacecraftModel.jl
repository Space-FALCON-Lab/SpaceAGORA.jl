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
#TODO: Add blunted cone
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
    r::MVector{3, Float64} # Position of the COM in the parent frame
    q::MVector{4, Float64} # Quaternion representing the orientation in the parent frame
end

mutable struct SpacecraftModel
    bodies::Vector{BodyNode} # List of body nodes
    root::Int # Index of the root node in the bodies vector
    q::Vector{Float64} # Generalized coordinates (joint angles, etc.)
    q_dot::Vector{Float64} # Generalized velocities
    joint_indices::Dict{Int, Int} # Body index => joint index in q, q_dot
    instant_actuation::Bool # Whether control inputs (e.g., solar panel angles) are applied instantly
    dry_mass::Float64 # Dry mass of the spacecraft
    prop_mass::Float64 # Fuel mass available for maneuvers
end

# Function to add a body to the spacecraft model
function add_body!(model::SpacecraftModel, 
                body::RigidBody, 
                joint::Joint, 
                parent_index::Union{Nothing, Int},
                r::MVector{3, Float64},
                q::MVector{4, Float64})
    """
    Adds a new body to the spacecraft model.
    - `model`: The spacecraft model to which the body will be added.
    - `body`: The rigid body to be added.
    - `joint`: The joint connecting this body to its parent.
    - `parent_index`: The index of the parent body in the model's bodies vector, or `nothing` if this is the root body.
    - `r`: The position of the COM of the body relative to the COM of the parent.
    - `q`: The quaternion representing the orientation of the body in the parent frame. For the root body, this is the orientation in the inertial frame.
    """

    index = length(model.bodies) + 1
    new_node = BodyNode(body, joint, parent_index, Int[], r, q)
    model.dry_mass += body.mass # Update dry mass with the new body's mass

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

function get_spacecraft_length(model::SpacecraftModel)
    """
    Calculates the total length of the spacecraft model by summing the lengths of all bodies.
    This is useful for determining the overall size of the spacecraft.
    """
    total_length = 0.0
    for node in model.bodies
        if node.body isa Box
            # Only consider bodies with defined dimensions
            dims = node.body.dimensions
            total_length += dims[1] # Use the maximum dimension as the length
        end
    end
    return total_length
end

function get_SA_area(model::SpacecraftModel)
    """
    Calculates the total surface area of the solar array by summing the areas of all flat plates.
    This is useful for aerodynamic calculations.
    """
    total_area = 0.0
    for node in model.bodies
        if node.body isa FlatPlate
            # Only consider bodies with defined area
            total_area += node.body.area
        end
    end
    return total_area
end

function get_SC_area(model::SpacecraftModel)
    """
    Calculates the total surface area of the spacecraft bus by summing the areas of all boxes.
    This is useful for aerodynamic calculations.
    """
    total_area = 0.0
    for node in model.bodies
        if node.body isa Box
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
    r = MVector{3, Float64}(r)
    q = MVector{4, Float64}(0.0, 0.0, 0.0, 1.0) # Identity quaternion for no rotation
    return r, q
end
