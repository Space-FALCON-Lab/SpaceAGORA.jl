"""
    get_COM(model::SpacecraftModel, body::Link)

Returns the center of mass of the spacecraft assembly that the body is a part of.
"""
function get_COM(model::SpacecraftModel)
    bodies = model.links # Get all bodies in the model
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
        R = b.root ? SMatrix{3, 3, Float64}(I(3)) : rot(b.q) # Rotation matrix from quaternion
        I_body = R * b.inertia * R' # Transform inertia tensor to the body frame
        r = SVector{3, Float64}(b.r) # Position vector of the body
        # Only the root assembly carries shared propellant mass in this legacy structure model.
        fuel_mass = b.root ? prop_mass : 0.0 # Get propellant mass if root body
        inertia_tensor += I_body + (b.m + fuel_mass) * hat(r) * hat(r)' # Parallel axis theorem
    end
    return inertia_tensor # Return the inertia tensor
end

"""
    get_inertia_tensor(model::SpacecraftModel, body::Link)

Returns the pre-calculated inertia tensor for the assembly connected to `body`.
"""
function get_inertia_tensor(model::SpacecraftModel)
    return model.inertia_tensor # Return the inertia tensor of the root body
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
    bodies, root_index = traverse_bodies(model, body)
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
