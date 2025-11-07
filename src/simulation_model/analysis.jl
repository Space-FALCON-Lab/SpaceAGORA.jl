module Analysis

# --- External Dependencies ---
using StaticArrays
using LinearAlgebra

# --- Internal Dependencies ---
# Import the data structs defined in other modules
using ..PhysicalModel: SpacecraftModel, Link, Joint
using ..Kinematics: rotate_to_inertial

# Import the quaternion helpers
# (Assuming quaternion_utils is available in the main module)
# using ..quaternion_utils: rot, hat

# --- Public API ---
export traverse_bodies,
       get_COM,
       update_inertia_tensor!,
       update_inertia_tensor,
       get_inertia_tensor,
       set_inertia_tensor!,
       get_spacecraft_mass,
       get_spacecraft_reference_area,
       get_spacecraft_length,
       get_SA_area,
       get_SC_area,
       get_normal_vector,
       get_tangent_vector

"""
    traverse_bodies(model::SpacecraftModel, body::Link)

Traverses the spacecraft model starting from the given body and returns all
bodies connected to it via joints, as well as the index of the root body
for this assembly.
"""
function traverse_bodies(model::SpacecraftModel, body::Link)
    visited = Set{Link}() # Set to keep track of visited bodies
    root_index = 0 # Index of the root body, if needed
    queue = Link[body] # Initialize the queue with the starting body
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

"""
    get_COM(model::SpacecraftModel, body::Link)

Returns the center of mass of the spacecraft assembly that the body is a part of.
"""
function get_COM(model::SpacecraftModel, body::Link)
    # BFS starting from body to find all bodies attached to the current body
    bodies, root_index = traverse_bodies(model, body)
    return get_COM(bodies) # Call the version for a list of bodies
end

"""
    get_COM(bodies::Vector{Link})

Returns the center of mass of a collection of bodies.
"""
function get_COM(bodies::Vector{Link})
    COM = MVector{3, Float64}(0.0, 0.0, 0.0)
    total_mass = 0.0
    for body in bodies
        COM += body.r * body.m # Sum the position vectors weighted by mass
        total_mass += body.m # Sum the mass
    end
    return COM / total_mass # Return the center of mass
end

"""
    update_inertia_tensor!(model::SpacecraftModel, body::Link)

Calculates the inertia tensor of the entire assembly connected to `body`
and updates it in the `model.inertia_tensors` array.
"""
function update_inertia_tensor!(model::SpacecraftModel, body::Link)
    # BFS starting from body to find all bodies attached to the current body
    bodies, root_index = traverse_bodies(model, body)

    # Calculate the inertia tensor
    inertia_tensor = update_inertia_tensor(bodies, model.prop_mass[root_index])
    
    # Store it in the model
    if root_index > length(model.inertia_tensors)
        # If the root index is out of bounds, extend the inertia tensors vector
        push!(model.inertia_tensors, inertia_tensor) # Add the inertia tensor for the root body
    else
        model.inertia_tensors[root_index] = inertia_tensor # Update the inertia tensor for the root body
    end
    return inertia_tensor # Return the inertia tensor
end

"""
    update_inertia_tensor(bodies::Vector{Link}, prop_mass::Float64 = 0.0)

Returns the total inertia tensor of a collection of bodies using the parallel axis theorem.
"""
function update_inertia_tensor(bodies::Vector{Link}, prop_mass::Float64 = 0.0)
    inertia_tensor = SMatrix{3, 3, Float64}(zeros(3, 3))
    for b in bodies
        # Apply parallel axis theorem to update inertia tensor
        R = b.root ? SMatrix{3, 3, Float64}(I(3)) : rot(b.q) # Rotation matrix from quaternion
        I_body = R * b.inertia * R' # Transform inertia tensor to the body frame
        r = SVector{3, Float64}(b.r) # Position vector of the body
        fuel_mass = b.root ? prop_mass : 0.0 # Get propellant mass if root body
        inertia_tensor += I_body + (b.m + fuel_mass) * hat(r) * hat(r)' # Parallel axis theorem
    end
    return inertia_tensor # Return the inertia tensor
end

"""
    get_inertia_tensor(model::SpacecraftModel, body::Link)

Returns the pre-calculated inertia tensor for the assembly connected to `body`.
"""
function get_inertia_tensor(model::SpacecraftModel, body::Link)
    # BFS starting from body to find all bodies attached to the current body
    bodies, root_index = traverse_bodies(model, body)
    return model.inertia_tensors[root_index] # Return the inertia tensor of the root body
end

"""
    get_inertia_tensor(model::SpacecraftModel, root_index::Int)

Returns the pre-calculated inertia tensor for the assembly at `root_index`.
"""
function get_inertia_tensor(model::SpacecraftModel, root_index::Int)
    @assert root_index <= length(model.inertia_tensors) "Root index out of bounds"
    return model.inertia_tensors[root_index] # Return the inertia tensor of the root body
end

"""
    set_inertia_tensor!(model::SpacecraftModel, body::Link, inertia_tensor::SMatrix{3, 3, Float64})

Manually sets the inertia tensor for the assembly connected to `body`.
"""
function set_inertia_tensor!(model::SpacecraftModel, body::Link, inertia_tensor::SMatrix{3, 3, Float64})
    # Find the index of the body in the links vector
    bodies, root_index = traverse_bodies(model, body)
    # @assert index !== nothing "Body not found in the model"
    model.inertia_tensors[root_index] = inertia_tensor # Set the inertia tensor for the body
end

"""
    get_spacecraft_mass(model::SpacecraftModel; dry=false)

Calculates the total mass of all assemblies in the spacecraft model.
Returns a Vector{Float64} if multiple roots exist, or a single Float64 if only one.
"""
function get_spacecraft_mass(model::SpacecraftModel; dry=false)
    masses = Float64[]
    for root in model.roots
        # BFS starting from root to find all bodies attached to the root body
        bodies, root_index = traverse_bodies(model, root)
        push!(masses, get_spacecraft_mass(model, bodies, root_index; dry=dry))
    end
    return length(model.roots) == 1 ? masses[1] : masses # Return the total mass of all roots
end

"""
    get_spacecraft_mass(model::SpacecraftModel, body::Link; dry=false)

Calculates the total mass of the assembly connected to `body`.
"""
function get_spacecraft_mass(model::SpacecraftModel, body::Link; dry=false)
    bodies, root_index = traverse_bodies(model, body)
    return get_spacecraft_mass(model, bodies, root_index; dry=dry)
end

"""
    get_spacecraft_mass(model::SpacecraftModel, bodies::Vector{Link}, root_index::Int; dry=false)

Calculates the total mass of a collection of bodies by summing their masses.
"""
function get_spacecraft_mass(model::SpacecraftModel, bodies::Vector{Link}, root_index::Int; dry=false)
    total_mass = sum([b.m for b in bodies]) # Sum the mass of each body
    return dry ? total_mass : total_mass + model.prop_mass[root_index] # Return the total mass
end

"""
    get_spacecraft_reference_area(model::SpacecraftModel)

Calculates the total reference area of all assemblies in the model.
"""
function get_spacecraft_reference_area(model::SpacecraftModel)
    total_area = Float64[]
    for root in model.roots
        bodies, root_index = traverse_bodies(model, root)
        push!(total_area, get_spacecraft_reference_area(bodies))
    end
    return length(total_area) == 1 ? total_area[1] : total_area # Return the total area of all roots
end

"""
    get_spacecraft_reference_area(model::SpacecraftModel, body::Link)

Calculates the total reference area of the assembly connected to `body`.
"""
function get_spacecraft_reference_area(model::SpacecraftModel, body::Link)
    bodies, root_index = traverse_bodies(model, body)
    return get_spacecraft_reference_area(bodies)
end

"""
    get_spacecraft_reference_area(bodies::Vector{Link})

Calculates the total reference area of a collection of bodies.
"""
function get_spacecraft_reference_area(bodies::Vector{Link})
    total_area = 0.0
    for body in bodies
        if body.ref_area > 0.0 # Only consider bodies with defined area
            total_area += body.ref_area # Sum the reference areas of each body
        end
    end
    return total_area # Return the total area
end

"""
    get_spacecraft_length(model::SpacecraftModel, body::Link)

Calculates the maximum dimension of any body in the assembly connected to `body`.
"""
function get_spacecraft_length(model::SpacecraftModel, body::Link)
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

"""
    get_SA_area(model::SpacecraftModel, body::Link)

Calculates the total surface area of all non-root bodies (assumed to be solar arrays).
"""
function get_SA_area(model::SpacecraftModel, body::Link)
    bodies, root_index = traverse_bodies(model, body)
    return get_SA_area(bodies)
end

"""
    get_SA_area(bodies::Vector{Link})

Calculates the total surface area of all non-root bodies in a list.
"""
function get_SA_area(bodies::Vector{Link})
    total_area = 0.0
    for b in bodies
        if !b.root # Only consider non-root bodies (flat plates)
            total_area += b.ref_area # Sum the areas of flat plates
        end
    end
    return total_area # Return the total area of the solar array
end

"""
    get_SC_area(model::SpacecraftModel, body::Link)

Calculates the total surface area of all root bodies (assumed to be the s/c bus).
"""
function get_SC_area(model::SpacecraftModel, body::Link)
    bodies, root_index = traverse_bodies(model, body)
    return get_SC_area(bodies)
end

"""
    get_SC_area(bodies::Vector{Link})

Calculates the total surface area of all root bodies in a list.
"""
function get_SC_area(bodies::Vector{Link})
    total_area = 0.0
    for b in bodies
        if b.root # Only consider root bodies (boxes)
            total_area += b.ref_area # Sum the areas of boxes
        end
    end
    return total_area # Return the total area of the spacecraft bus
end

"""
    get_normal_vector(model::SpacecraftModel, body::Link, root_index::Int; normalized=false)

Returns the normal vector (body x-axis) of the body in the inertial frame.
"""
function get_normal_vector(model::SpacecraftModel, body::Link, root_index::Int; normalized=false)
    # Convert to inertial frame using the orientation quaternion
    R = rotate_to_inertial(model, body, root_index) # Get the rotation matrix to convert from body frame to inertial frame

    # Return the normal vector in the inertial frame
    normal = R * SVector{3, Float64}(1.0, 0.0, 0.0) # Normal vector in inertial frame
    
    return normalized ? normalize(normal) : normal
end

"""
    get_tangent_vector(model::SpacecraftModel, body::Link, root_index::Int; normalized=false)

Returns the tangent vector (body z-axis) of the body in the inertial frame.
"""
function get_tangent_vector(model::SpacecraftModel, body::Link, root_index::Int; normalized=false)
    # Convert to inertial frame using the orientation quaternion
    R = rotate_to_inertial(model, body, root_index) # Get the rotation matrix to convert from body frame to inertial frame

    # Return the normal vector in the inertial frame
    tangent = R * SVector{3, Float64}(0.0, 0.0, 1.0) # Tangent vector in inertial frame
    
    return normalized ? normalize(tangent) : tangent
end

end # module Analysis
