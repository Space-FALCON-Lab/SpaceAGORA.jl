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
    get_inertia_tensor(model::SpacecraftModel, body::Link)

Returns the pre-calculated inertia tensor for the assembly connected to `body`.
"""
function get_inertia_tensor(model::SpacecraftModel)
    return model.inertia_tensor # Return the inertia tensor of the root body
end
