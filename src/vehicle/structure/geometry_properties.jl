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
