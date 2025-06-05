include("utils/quaternion_utils.jl")
import .quaternion_utils
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

    function FlatPlate(name::String, mass::Float64, dimensions::SVector{2, Float64}, area::Float64, COM::SVector{3})
        thickness = 0.002 # Assume a small thickness for the plate (2mm, standard for solar panels)
        # Inertia tensor for a box in principal axes
        inertia_tensor = SMatrix{3, 3, Float64}(
            mass * (dimensions[1]^2 + dimensions[2]^2) / 12, 0, 0,
            0, mass * (thickness^2 + dimensions[2]^2) / 12, 0,
            0, 0, mass * (thickness^2 + dimensions[1]^2) / 12
        )
        return new(name, mass, inertia_tensor, dimensions, area, COM)
    end
end

struct Box <: RigidBody
    name::String
    mass::Float64
    inertia_tensor::SMatrix{3, 3, Float64} # In principal axes
    dimensions::SVector{3, Float64} # Lengths in x, y, z directions
    area::Float64 # Aerodynamic reference area
    COM::SVector{3, Float64} # Center of mass in body frame, i.e., offset from geometric center

    function Box(name::String, mass::Float64, inertia_tensor::SMatrix{3, 3, Float64}, dimensions::SVector{3, Float64})
        area = dimensions[1] * dimensions[2] # Assume area is the front face area
        COM = SVector{3}(0.0, 0.0, 0.0) # Center of mass at geometric center
        return new(name, mass, inertia_tensor, dimensions, area, COM)
    end

    function Box(name::String, mass::Float64, dimensions::SVector{3, Float64}, area::Float64, COM::SVector{3})
        # Inertia tensor for a box in principal axes
        inertia_tensor = SMatrix{3, 3, Float64}(
            mass * (dimensions[2]^2 + dimensions[3]^2) / 12, 0, 0,
            0, mass * (dimensions[1]^2 + dimensions[3]^2) / 12, 0,
            0, 0, mass * (dimensions[1]^2 + dimensions[2]^2) / 12
        )
        return new(name, mass, inertia_tensor, dimensions, area, COM)
    end
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
    r::SVector{3, Float64} # Position of the COM in the parent frame
    r_body::SVector{3, Float64} # Position of the COM in the body frame (relative to origin of the body)
    q::SVector{4, Float64} # Quaternion representing the orientation in the parent frame
    O_body::SMatrix{3, 3, Float64} # Orientation matrix representing the orientation in the body frame
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
    inertia_tensor::MMatrix{3, 3, Float64} # Inertia tensor of the spacecraft in the inertial frame
    COM::SVector{3, Float64} # Center of mass in the body frame (relative to origin of the root body)
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
    new_node = BodyNode(body, joint, parent_index, Int[], r, SVector{3}(zeros(3)), q, SMatrix{3, 3, Float64}(zeros(3, 3)))
    

    if isnothing(parent_index)
        push!(model.bodies, new_node)
        model.root = length(model.bodies)
        model.inertia_tensor = body.inertia_tensor # Set inertia tensor for the root body
        model.COM = body.COM # Set center of mass for the root body
    else
        parent_node = model.bodies[parent_index]
        new_node.parent = parent_index
        push!(parent_node.children, index)
        push!(model.bodies, new_node)
        # Update inertia tensor and COM for the spacecraft model
        # Update COM
        path_to_root = get_path_to_root(model, index) # Get path from body to root, starting from the new body
        r = [0;0;0]
        O_total = I3 # Initialize orientation matrix to identity
        # Calculate the position in the body frame
        for i in path_to_root
            body_node = model.bodies[i]
            R = rot(body_node.q) # Rotation matrix from quaternion ()
            O = R' # Orientation matrix from rotation matrix
            r = O * (r + body_node.r) # Update position into the current frame and add the new position
            O_total = O*O_total # Update the total orientation matrix
        end
        model.COM = (model.COM * model.dry_mass + r * body.mass) / (model.dry_mass + body.mass) # Update COM with the new body's COM
        new_node.r_body = r # Set the position of the COM in the body frame
        new_node.O_body = O_total # Update the orientation in the body frame
        # Update inertia tensor
        model.inertia_tensor = update_inertia_tensor(model) # Inertia tensor in the body frame
        println("Adding body: $(body.name) at index $index with parent $parent_index")
        println("Position in body frame: $(new_node.r_body)")
        println("Orientation in body frame: $(new_node.O_body)")
        println("Updated COM: $(model.COM)")
        println("Updated inertia tensor: $(model.inertia_tensor)")
    end

    model.dry_mass += body.mass # Update dry mass with the new body's mass

    if joint isa RevoluteJoint
        push!(model.q, 0.0) # Initialize joint angle
        push!(model.q_dot, 0.0) # Initialize joint velocity
        model.joint_indices[index] = length(model.q) # Map body index to joint index
    end
end

function get_path_to_root(model::SpacecraftModel, body_index::Int)
    """
    Returns the path from the given body index to the root of the spacecraft model.
    - `model`: The spacecraft model.
    - `body_index`: The index of the body for which to find the path.
    """
    path = Int[]
    current_index = body_index
    while current_index != model.root
        push!(path, current_index)
        current_index = model.bodies[current_index].parent
    end
    push!(path, model.root) # Include the root node
    return path # Return path from body to root
end

function update_inertia_tensor(model::SpacecraftModel)
    """
    Updates the inertia tensor of the spacecraft model based on the current configuration of bodies.
    This is useful for dynamics calculations.
    """
    inertia_tensor = SMatrix{3, 3, Float64}(zeros(3, 3))
    for node in model.bodies
        # Apply parallel axis theorem to update inertia tensor
        R = node.O_body'
        I_body = R*node.body.inertia_tensor*R' # Transform inertia tensor to the body frame
        inertia_tensor += I_body + node.body.mass * hat(node.r_body) * hat(node.r_body)' # Parallel axis theorem
    end
    return inertia_tensor
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
