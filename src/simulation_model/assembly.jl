module Assembly

using LinearAlgebra
using StaticArrays
using ..PhysicalModel
using ..Components
# using ..AbstractTypes

export add_body!, add_joint!, add_facet!, add_magnet!, add_thruster!

# NOTE: I have removed `add_dynamic_effector!`.
# This function is not type-stable. You cannot append to a tuple
# without changing the type of the parent struct.
# The effectors must be passed in as a tuple to the
# SpacecraftModel constructor (see src/model.jl).

# Function to add a body to the spacecraft model
function add_body!(model::SpacecraftModel,
    body::Link;
    prop_mass::Union{Nothing, Float64}=nothing)
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
    model.n_reaction_wheels += length(body.rw_assembly.h_wheels) # Increment the number of reaction wheels
    if body.root
        append!(model.inertia_tensors, [SMatrix{3, 3, Float64}(I(3))]) # Add inertia tensor for root body
        # update_inertia_tensor!(model, body) # This function is in Analysis.jl, creates circular dependency
        # We should call this from a higher-level script after assembly.
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
    # update_inertia_tensor!(model, joint.link1) # Calculate inertia tensor
end

function add_facet!(link::Link, facet::Facet)
    push!(link.SRP_facets, facet)
end

function add_facet!(link::Link, facets::Vector{Facet})
    push!(link.SRP_facets, facets...)
end

function add_magnet!(link::Link, magnet::Magnet)
    push!(link.magnets, magnet)
end

function add_magnet!(link::Link, magnets::Vector{Magnet})
    push!(link.magnets, magnets...)
end

function add_thruster!(model::SpacecraftModel, link::Link, thruster::Thruster)
    """
    Adds a thruster to the Link
    - `link` : The Link to which the thruster is added
    - `thruster` : The thruster to add
    """
    # Add the thruster torque to the Jacobian matrix
    normalize!(thruster.direction) # Ensure the thruster direction is a unit vector
    
    if isempty(link.thrusters)
        link.J_thruster = cross(thruster.location, thruster.direction) # If this is the first thruster to be added, simply set the Jacobian equal to the r x F vector
    else
        link.J_thruster = hcat(link.J_thruster, cross(thruster.location, thruster.direction)) # If there are already thrusters in the link, append the new r x F vector to the Jacobian
    end

    # Append the thruster to the list of thrusters for future reference
    push!(link.thrusters, thruster)
    model.n_thrusters += 1 # Increment the number of thrusters in the spacecraft model
end

end # module Assembly