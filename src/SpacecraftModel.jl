include("utils/quaternion_utils.jl")
import .quaternion_utils
using StaticArrays
using LinearAlgebra
const I3 = SMatrix{3, 3, Float64}(diagm(ones(3)))

mutable struct Link
    root::Bool # Whether this link is a root link (i.e., the main bus or core body of the spacecraft). 
    #^ This will have orientation expressed relative to the inertial frame.
    r::MVector{3} # Position of COM (Body frame)
    q::MVector{4} # Orientation (Body frame)
    ṙ::MVector{3} # velocity (body frame)
    ω::MVector{3} # angular velocity 
    dims::MVector{3} # Size[x,y,z] Box= x=thikness, y=z=width, height
    ref_area::Float64 # Reference area for aerodynamic calculations
    m::Float64 # Mass
    mass::SMatrix{3, 3} # Mass Matrix
    inertia::SMatrix{3,3} # Inertial matrix
    aᵇ::SVector{3} # Left extent (Body frame)
    bᵇ::SVector{3} # Right extent (Body frame)
    gyro::Int64 # Number of Gyroscope
    rw::MVector{} # angular velocity reaction wheels
    J_rw::Matrix{Float64}#J_rw::SMatrix{3,} # reaction wheel jacobian
    net_force::MVector{3, Float64} # Net force acting on the link, to be updated at each simulation step
    net_torque::MVector{3, Float64} # Net torque acting on the link, to be updated at each simulation step
    function Link(;root=false,
                    r=SVector{3, Float64}([0, 0, 0]), 
                    q=SVector{4, Float64}([0,0,0,1]), 
                    ṙ=SVector{3, Float64}([0,0,0]), 
                    ω=SVector{3, Float64}([0,0,0]), 
                    dims=SVector{3, Float64}([0.5, 0.5, 0.1]), 
                    ref_area=1.0,
                    m = 3.0,
                    mass=SMatrix{3,3}(m*I3), 
                    inertia=SMatrix{3,3}(1 / 12 * m * diagm([dims[2]^2 + dims[3]^2; dims[1]^2 + dims[3]^2; dims[1]^2 + dims[2]^2])),
                    a=SVector{3, Float64}([-0.5*dims[1], 0, 0]),
                    b=SVector{3, Float64}([0.5*dims[1], 0, 0]),
                    gyro = 3,
                    rw = zeros(SVector{Int(gyro), Float64}),
                    J_rw=zeros(3, gyro),
                    net_force=MVector{3, Float64}(zeros(3)),
                    net_torque=MVector{3, Float64}(zeros(3)))#SMatrix{3,Int(gyro)}(1.0I))
        new(root, r, q, ṙ, ω, dims, ref_area, m, mass, inertia, a, b, gyro, rw, J_rw, net_force, net_torque)
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
                    Kx=SMatrix{3,3}(0.0I), 
                    Kt=SMatrix{3,3}(0.0I), 
                    Cx=zeros(SMatrix{3,3}),
                    Ct=zeros(SMatrix{3,3}))
        new(link1, link2, p1ᵇ, p2ᵇ, Kx, Kt, Cx, Ct,
            SVector{3}(0.0, 0.0, 0.0), SVector{4}(0.0, 0.0, 0.0, 1.0))
    end

    function Joint(link1::Link, link2::Link; p1=link1.bᵇ, 
                                 p2=link2.aᵇ, 
                                 Kx=SMatrix{3,3}(0.0I), 
                                 Kt=SMatrix{3,3}(0.0I), 
                                 Cx=zeros(SMatrix{3,3}),
                                 Ct=zeros(SMatrix{3,3}),
                                 translational_displacement=SVector{3, Float64}(0.0, 0.0, 0.0),
                                 rotational_displacement=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
        
        new(link1, link2, p1, p2, Kx, Kt, Cx, Ct, 
            translational_displacement, rotational_displacement)
    end

    function Joint(;link1=Link(), link2=Link(), p1=link1.bᵇ, 
        p2=link2.aᵇ, 
        Kx=SMatrix{3,3}(1.0I), 
        Kt=SMatrix{3,3}(1.0I), 
        Cx=zeros(SMatrix{3,3}),
        Ct=zeros(SMatrix{3,3}),
        translational_displacement=SVector{3}(0.0, 0.0, 0.0),
        rotational_displacement=SVector{4}(0.0, 0.0, 0.0, 1.0))

        new(link1, link2, p1, p2, Kx, Kt, Cx, Ct, 
            translational_displacement, rotational_displacement)
    end
end
# abstract type RigidBody end

# # Define various rigid body types with their properties
# # Flat plate for solar panels, boxes for main bus
# #TODO: Add blunted cone
# struct FlatPlate <: RigidBody
#     name::String
#     mass::Float64
#     inertia_tensor::SMatrix{3, 3, Float64} # In principal axes
#     dimensions::SVector{2, Float64} # Length/width dimensions
#     area::Float64 # Reference area
#     COM::SVector{3, Float64} # Center of mass in body frame, i.e., offset from geometric center

#     function FlatPlate(name::String, mass::Float64, dimensions::SVector{2, Float64}, area::Float64, COM::SVector{3})
#         thickness = 0.002 # Assume a small thickness for the plate (2mm, standard for solar panels)
#         # Inertia tensor for a box in principal axes
#         inertia_tensor = SMatrix{3, 3, Float64}(
#             mass * (dimensions[1]^2 + dimensions[2]^2) / 12, 0, 0,
#             0, mass * (thickness^2 + dimensions[2]^2) / 12, 0,
#             0, 0, mass * (thickness^2 + dimensions[1]^2) / 12
#         )
#         return new(name, mass, inertia_tensor, dimensions, area, COM)
#     end
# end

# struct Box <: RigidBody
#     name::String
#     mass::Float64
#     inertia_tensor::SMatrix{3, 3, Float64} # In principal axes
#     dimensions::SVector{3, Float64} # Lengths in x, y, z directions
#     area::Float64 # Aerodynamic reference area
#     COM::SVector{3, Float64} # Center of mass in body frame, i.e., offset from geometric center

#     function Box(name::String, mass::Float64, inertia_tensor::SMatrix{3, 3, Float64}, dimensions::SVector{3, Float64})
#         area = dimensions[1] * dimensions[2] # Assume area is the front face area
#         COM = SVector{3}(0.0, 0.0, 0.0) # Center of mass at geometric center
#         return new(name, mass, inertia_tensor, dimensions, area, COM)
#     end

#     function Box(name::String, mass::Float64, dimensions::SVector{3, Float64}, area::Float64, COM::SVector{3})
#         # Inertia tensor for a box in principal axes
#         inertia_tensor = SMatrix{3, 3, Float64}(
#             mass * (dimensions[2]^2 + dimensions[3]^2) / 12, 0, 0,
#             0, mass * (dimensions[1]^2 + dimensions[3]^2) / 12, 0,
#             0, 0, mass * (dimensions[1]^2 + dimensions[2]^2) / 12
#         )
#         return new(name, mass, inertia_tensor, dimensions, area, COM)
#     end
# end

# mutable struct BodyNode
#     body::RigidBody
#     joint::Joint
#     parent::Union{Nothing, Int} # Parent node in the tree structure
#     children::Vector{Int} # Children nodes in the tree structure
#     r::SVector{3, Float64} # Position of the COM in the parent frame
#     r_body::SVector{3, Float64} # Position of the COM in the body frame (relative to origin of the body)
#     q::SVector{4, Float64} # Quaternion representing the orientation in the parent frame
#     O_body::SMatrix{3, 3, Float64} # Orientation matrix representing the orientation in the body frame
# end

mutable struct SpacecraftModel
    joints::Vector{Joint} # List of joints
    links::Vector{Link} # List of links (bodies)
    roots::Vector{Link} # Vector of root links (main bus or core bodies)
    instant_actuation::Bool # Whether control inputs (e.g., solar panel angles) are applied instantly
    # dry_mass::Float64 # Dry mass of the spacecraft
    prop_mass::Vector{Float64} # Fuel mass available for maneuvers
    inertia_tensors::Vector{SMatrix{3, 3, Float64}} # Inertia tensors of the spacecraft bodies
    # inertia_tensor::MMatrix{3, 3, Float64} # Inertia tensor of the spacecraft in the inertial frame
    # COM::SVector{3, Float64} # Center of mass in the body frame (relative to origin of the root body)

    function SpacecraftModel(;joints=Joint[], links=Link[], roots=Link[], 
                        instant_actuation=true, 
                        prop_mass=Float64[],
                        inertia_tensors=SMatrix{3, 3, Float64}[])
        new(joints, links, roots, instant_actuation, prop_mass, inertia_tensors)
    end
end

# Function to add a body to the spacecraft model
function add_body!(model::SpacecraftModel, 
                body::Link; prop_mass::Union{Nothing, Float64}=nothing)
    """
    Adds a new body to the spacecraft model.
    - `model`: The spacecraft model to which the body will be added.
    - `body`: The rigid body to be added.
    """

    push!(model.links, body) # Add the body to the links vector
    if body.root
        # If the body is a root body, it should be added to the roots vector
        @assert prop_mass !== nothing "Propellant mass must be provided for root body"
        push!(model.roots, body) # Add the body to the roots vector
        push!(model.prop_mass, prop_mass) # Add propellant mass for root body
        # Initialize inertia tensor for the root body
    end
end

# Function to add a joint to the spacecraft model
function add_joint!(model::SpacecraftModel, joint::Joint)
    """
    Adds a new joint to the spacecraft model.
    - `model`: The spacecraft model to which the joint will be added.
    - `joint`: The joint to be added.
    """
    push!(model.joints, joint) # Add the joint to the joints vector
    update_inertia_tensor!(model, joint.link1) # Calculate inertia tensor
end

function traverse_bodies(model::SpacecraftModel, body::Link)
    """
    Traverses the spacecraft model starting from the given body and returns all bodies connected to it.
    - `model`: The spacecraft model.
    - `body`: The body from which to start the traversal.
    """
    visited = Set{Link}() # Set to keep track of visited bodies
    root_index = 0 # Index of the root body, if needed
    queue = [body] # Initialize the queue with the starting body
    push!(visited, body) # Mark the initial body as visited
    while !isempty(queue)
        current_body = popfirst!(queue) # Get the next body in the queue
        if current_body.root
            root_index = findfirst(isequal(current_body), model.roots) # Find the index of the root body
        end
        # Add children bodies to the queue
        for joint in model.joints
            if joint.link1 == current_body && !in(joint.link2, visited)
                push!(queue, joint.link2) # Add the second link of the joint to the queue
                push!(visited, joint.link2) # Mark it as visited
            elseif joint.link2 == current_body && !in(joint.link1, visited)
                push!(queue, joint.link1) # Add the first link of the joint to the queue
                push!(visited, joint.link1) # Mark it as visited
            end
        end
    end
    @assert root_index != 0 "Root body not found in the model" # Ensure a root body was found
    return collect(visited), root_index # Return all visited bodies as a vector
end
function get_COM(model::SpacecraftModel, body::Link)
    """
    Returns the center of mass of the spacecraft that the body is a part of.
    - `model`: The spacecraft model.
    - `body`: The body for which to get the center of mass.
    """
    # Initialize COM to zero
    COM = MVector{3, Float64}(0.0, 0.0, 0.0)
    total_mass = 0.0
    # BFS starting from body to find all bodies attached to the current body
    bodies, root_index = traverse_bodies(model, body)
    println("Bodies found: ", length(bodies))
    COM = sum([b.r * b.m for b in bodies]) # Sum the position vectors weighted by mass
    total_mass = sum([b.m for b in bodies]) # Sum the mass of all bodies
    return COM / total_mass # Return the center of mass
end

function get_COM(bodies::Vector{Link})
    """
    Returns the center of mass of a collection of bodies.
    - `bodies`: A vector of bodies for which to calculate the center of mass.
    """
    COM = MVector{3, Float64}(0.0, 0.0, 0.0)
    total_mass = 0.0
    for body in bodies
        COM += body.r * body.m # Sum the position vectors weighted by mass
        total_mass += body.m # Sum the mass
    end
    return COM / total_mass # Return the center of mass
end

# Function to get the inertia tensor of a body in the spacecraft model
function update_inertia_tensor!(model::SpacecraftModel, body::Link)
    """
    Returns the inertia tensor of the spacecraft that the body is a part of.
    - `model`: The spacecraft model.
    - `body`: The body for which to get the inertia tensor.
    """
    # Initialize inertia tensor to zero
    inertia_tensor = SMatrix{3, 3, Float64}(zeros(3, 3))
    total_mass = 0.0
    # BFS starting from body to find all bodies attached to the current body
    bodies, root_index = traverse_bodies(model, body)
    for b in bodies
        # Apply parallel axis theorem to update inertia tensor
        R = rot(b.q) # Rotation matrix from quaternion
        I_body = R * b.inertia * R' # Transform inertia tensor to the body frame
        inertia_tensor += I_body + b.m * hat(b.r) * hat(b.r)' # Parallel axis theorem
    end
    if root_index > length(model.inertia_tensors)
        # If the root index is out of bounds, extend the inertia tensors vector
        push!(model.inertia_tensors, inertia_tensor) # Add the inertia tensor for the root body
    else
        model.inertia_tensors[root_index] = inertia_tensor # Update the inertia tensor for the root body
    end
    return inertia_tensor # Return the inertia tensor
end

function update_inertia_tensor(bodies::Vector{Link})
    """
    Returns the inertia tensor of a collection of bodies.
    - `bodies`: A vector of bodies for which to calculate the inertia tensor.
    """
    inertia_tensor = SMatrix{3, 3, Float64}(zeros(3, 3))
    total_mass = 0.0
    for body in bodies
        # Apply parallel axis theorem to update inertia tensor
        R = rot(body.q) # Rotation matrix from quaternion
        I_body = R * body.inertia * R' # Transform inertia tensor to the body frame
        inertia_tensor += I_body + body.m * hat(body.r) * hat(body.r)' # Parallel axis theorem
        total_mass += body.m # Sum the mass
    end
    return inertia_tensor # Return the inertia tensor
end

function get_inertia_tensor(model::SpacecraftModel, body::Link)
    """
    Returns the inertia tensor of the spacecraft that the body is a part of.
    - `model`: The spacecraft model.
    - `body`: The body for which to get the inertia tensor.
    """
    # Initialize inertia tensor to zero
    inertia_tensor = SMatrix{3, 3, Float64}(zeros(3, 3))
    total_mass = 0.0
    # BFS starting from body to find all bodies attached to the current body
    bodies, root_index = traverse_bodies(model, body)
    return model.inertia_tensors[root_index] # Return the inertia tensor of the root body    
end

function get_inertia_tensor(model::SpacecraftModel, root_index::Int)
    """
    Returns the inertia tensor of the spacecraft that the body is a part of.
    - `model`: The spacecraft model.
    - `root_index`: The index of the root body for which to get the inertia tensor.
    """
    @assert root_index <= length(model.inertia_tensors) "Root index out of bounds"
    return model.inertia_tensors[root_index] # Return the inertia tensor of the root body
end


function get_spacecraft_mass(model::SpacecraftModel, body::Link; dry=false)
    """
    Calculates the total mass of the spacecraft model by summing the masses of all bodies connected to body.
    This is useful for dynamics calculations.
    - `model`: The spacecraft model.
    - `body`: The body for which to calculate the total mass.
    """
    total_mass = 0.0
    # BFS starting from body to find all bodies attached to the current body
    bodies, root_index = traverse_bodies(model, body)
    for b in bodies
        total_mass += b.m # Sum the mass of each body
    end
    return dry ? total_mass : total_mass + model.prop_mass[root_index] # Return the total mass
end

function get_spacecraft_mass(model::SpacecraftModel, bodies::Vector{Link}, root_index::Int;dry=false)
    """
    Calculates the total mass of a collection of bodies by summing their masses.
    - `bodies`: A vector of bodies for which to calculate the total mass.
    """
    total_mass = sum([b.m for b in bodies]) # Sum the mass of each body
    return dry ? total_mass : total_mass + model.prop_mass[root_index] # Return the total mass
end

function get_spacecraft_mass(model::SpacecraftModel;dry=false)
    """
    Calculates the total mass of the spacecraft model by summing the masses of all bodies.
    This is useful for dynamics calculations.
    - `model`: The spacecraft model.
    """
    masses = Float64[]
    for root in model.roots
        # BFS starting from root to find all bodies attached to the root body
        bodies, root_index = traverse_bodies(model, root)
        total_mass = sum([b.m for b in bodies]) # Sum the mass of each body
        push!(masses, dry ? total_mass : total_mass + model.prop_mass[root_index]) # Return the total mass
    end
    return length(model.roots) == 1 ? masses[1] : masses # Return the total mass of all roots
end

function get_spacecraft_reference_area(model::SpacecraftModel)
    """
    Calculates the total reference area of the spacecraft model by summing the areas of all bodies.
    This is useful for aerodynamic calculations.
    - `model`: The spacecraft model.
    """
    total_area = Float64[]
    for root in model.roots
        area = 0.0
        # BFS starting from root to find all bodies attached to the root body
        bodies, root_index = traverse_bodies(model, root)
        for body in bodies
            if body.ref_area > 0.0 # Only consider bodies with defined area
                area += body.ref_area # Sum the reference areas of each body
            end
        end
        push!(total_area, area) # Add the total area for this root
    end
    return length(total_area) == 1 ? total_area[1] : total_area # Return the total area of all roots
end

function get_spacecraft_reference_area(model::SpacecraftModel, body::Link)
    """
    Calculates the total reference area of a collection of bodies by summing their areas.
    - `bodies`: A vector of bodies for which to calculate the total reference area.
    """
    total_area = 0.0
    # BFS starting from body to find all bodies attached to the current body
    bodies, root_index = traverse_bodies(model, body)
    for body in bodies
        if body.ref_area > 0.0 # Only consider bodies with defined area
            total_area += body.ref_area # Sum the reference areas of each body
        end
    end
    return total_area # Return the total area
end

function get_spacecraft_length(model::SpacecraftModel, body::Link)
    """
    Calculates the total length of the spacecraft model by summing the lengths of all bodies.
    This is useful for determining the overall size of the spacecraft.
    - `model`: The spacecraft model.
    - `body`: The body for which to calculate the total length.
    """
    max_length = 0.0
    # BFS starting from body to find all bodies attached to the current body
    bodies, root_index = traverse_bodies(model, body)
    for b in bodies
        if b.dims[1] > max_length # Use the maximum dimension as the length
            max_length = b.dims[1] # Update the maximum length
        end
    end
    return max_length # Return the maximum length found
end

function get_SA_area(model::SpacecraftModel, body::Link)
    """
    Calculates the total surface area of the solar array by summing the areas of all flat plates.
    This is useful for aerodynamic calculations.
    - `model`: The spacecraft model.
    - `body`: The body for which to calculate the total surface area of the solar array.
    """
    total_area = 0.0
    # BFS starting from body to find all bodies attached to the current body
    bodies, root_index = traverse_bodies(model, body)
    for b in bodies
        if !b.root # Only consider flat plates
            total_area += b.ref_area # Sum the areas of flat plates
        end
    end
    return total_area # Return the total area of the solar array
end

function get_SC_area(model::SpacecraftModel, body::Link)
    """
    Calculates the total surface area of the spacecraft bus by summing the areas of all boxes.
    This is useful for aerodynamic calculations.
    - `model`: The spacecraft model.
    - `body`: The body for which to calculate the total surface area of the spacecraft bus.
    """
    total_area = 0.0
    # BFS starting from body to find all bodies attached to the current body
    bodies, root_index = traverse_bodies(model, body)
    for b in bodies
        if b.root # Only consider boxes
            total_area += b.ref_area # Sum the areas of boxes
        end
    end
    return total_area # Return the total area of the spacecraft bus
end
# function get_path_to_root(model::SpacecraftModel, body_index::Int)
#     """
#     Returns the path from the given body index to the root of the spacecraft model.
#     - `model`: The spacecraft model.
#     - `body_index`: The index of the body for which to find the path.
#     """
#     path = Int[]
#     current_index = body_index
#     while current_index != model.root
#         push!(path, current_index)
#         current_index = model.bodies[current_index].parent
#     end
#     push!(path, model.root) # Include the root node
#     return path # Return path from body to root
# end

# function update_inertia_tensor(model::SpacecraftModel)
#     """
#     Updates the inertia tensor of the spacecraft model based on the current configuration of bodies.
#     This is useful for dynamics calculations.
#     """
#     inertia_tensor = SMatrix{3, 3, Float64}(zeros(3, 3))
#     for node in model.bodies
#         # Apply parallel axis theorem to update inertia tensor
#         R = node.O_body'
#         I_body = R*node.body.inertia_tensor*R' # Transform inertia tensor to the body frame
#         inertia_tensor += I_body + node.body.mass * hat(node.r_body) * hat(node.r_body)' # Parallel axis theorem
#     end
#     return inertia_tensor
# end
# function get_spacecraft_reference_area(model::SpacecraftModel)
#     """
#     Calculates the total reference area of the spacecraft model by summing the areas of all bodies,
#     excluding actuators (reaction wheels, cmg, etc.) because they are internal.
#     """
#     total_area = 0.0
#     for node in model.bodies
#         if node.body isa FlatPlate || node.body isa Box
#             # Only consider bodies with defined area
#             total_area += node.body.area
#         end
#     end 
#     return total_area
# end

# function get_spacecraft_length(model::SpacecraftModel)
#     """
#     Calculates the total length of the spacecraft model by summing the lengths of all bodies.
#     This is useful for determining the overall size of the spacecraft.
#     """
#     total_length = 0.0
#     for node in model.bodies
#         if node.body isa Box
#             # Only consider bodies with defined dimensions
#             dims = node.body.dimensions
#             total_length += dims[1] # Use the maximum dimension as the length
#         end
#     end
#     return total_length
# end

# function get_SA_area(model::SpacecraftModel)
#     """
#     Calculates the total surface area of the solar array by summing the areas of all flat plates.
#     This is useful for aerodynamic calculations.
#     """
#     total_area = 0.0
#     for node in model.bodies
#         if node.body isa FlatPlate
#             # Only consider bodies with defined area
#             total_area += node.body.area
#         end
#     end
#     return total_area
# end

# function get_SC_area(model::SpacecraftModel)
#     """
#     Calculates the total surface area of the spacecraft bus by summing the areas of all boxes.
#     This is useful for aerodynamic calculations.
#     """
#     total_area = 0.0
#     for node in model.bodies
#         if node.body isa Box
#             # Only consider bodies with defined area
#             total_area += node.body.area
#         end
#     end
#     return total_area
# end

# # Function to create a transformation matrix for a translation
# function translation(r::SVector{3, Float64})
#     """
#     Creates a transformation matrix for a translation by vector `r`.
#     - `r`: The translation vector in the body frame.
#     """
#     r = MVector{3, Float64}(r)
#     q = MVector{4, Float64}(0.0, 0.0, 0.0, 1.0) # Identity quaternion for no rotation
#     return r, q
# end
