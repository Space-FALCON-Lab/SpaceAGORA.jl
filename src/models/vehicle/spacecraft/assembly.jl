module Assembly

using LinearAlgebra
using StaticArrays
using ..SpacecraftModels
using ..Components
# using ..AbstractTypes

export add_body!, add_joint!, add_facet!, add_magnet!, add_thruster!

# NOTE: I have removed `add_dynamic_effector!`.
# This function is not type-stable. You cannot append to a tuple
# without changing the type of the parent struct.
# The effectors must be passed in as a tuple to the
# SpacecraftModels.SpacecraftModel constructor.

# Function to add a body to the spacecraft model
function add_body!(model::SpacecraftModel,
    body::Link;
    prop_mass::Union{Nothing, Float64}=nothing)
    """
    Adds a new body to the spacecraft model.
    - `model`: The spacecraft model to which the body will be added.
    - `body`: The rigid body to be added.
    """

    
    if body.root
        # If the body is a root body, it should be added to the roots vector
        @assert prop_mass !== nothing "Propellant mass must be provided for root body"
        # push!(model.roots, body) # Add the body to the roots vector
        model.root = body # Set the root field to the body
        model.prop_mass = prop_mass # Add propellant mass for root body
        push!(model.links, body) # Add the body to the links vector
        # Initialize inertia tensor for the root body
    else
        push!(model.links, body) # Add the body to the links vector
    end
    if body.root
        model.inertia_tensor = SMatrix{3, 3, Float64}(I(3)) # Add inertia tensor for root body
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
    # `direction` is the thrust direction; a zero direction is kept as is rather
    # than normalized to NaN. The thrusters' force and torque columns come from
    # `allocation_matrix(model)`.
    n = norm(thruster.direction)
    (isfinite(n) && n > 0.0) && (thruster.direction ./= n) # Ensure the thruster direction is a unit vector
    push!(link.thrusters, thruster)
    return thruster
end

end # module Assembly
