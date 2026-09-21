module Kinematics

using StaticArrays
using LinearAlgebra
using ..SpacecraftModels
include(joinpath(@__DIR__, "..", "..", "core", "numerics", "quaternion_utils.jl"))

export rotate_to_inertial, rotate_to_body, rotate_link

function rotate_to_inertial(model::SpacecraftModel, body::Link, root_index::Int)
    """
    Returns the rotation matrix to convert from the link frame to the inertial frame.
    """
    if body.root
        return rot(body.q)' # Rotation matrix from quaternion
    else
        # Child link attitudes are stored relative to the root body frame.
        return rot(model.root.q)' * rot(body.q)' # Rotation matrix from quaternion
    end
end

function rotate_to_body(body::Link)
    """
    Returns the rotation matrix to convert from the link frame to the body frame
    """
    if body.root
        return I(3)
    else
        return rot(body.q)'
    end
end
function rotate_link(body::Link, q::SVector{4, Float64})
    """
    Rotates the link to a new orientation defined by the quaternion `q`.
    """
    @assert !body.root "Cannot rotate a root body directly"
    # Update the orientation of the body
    body.q .= project_unit_quaternion(q)
end

function rotate_link(body::Link, dcm::SMatrix{3, 3, Float64})
    """
    Rotates the link to a new orientation defined by the direction cosine matrix `dcm`.
    """
    @assert !body.root "Cannot rotate a root body directly"
    # Update the orientation of the body
    body.q .= dcm_to_quaternion(dcm)
end

function rotate_link(body::Link, axis::SVector{3, Float64}, θ::Float64)
    """
    Rotates the link to a new orientation defined by the Euler angles `θ`.
    """
    @assert !body.root "Cannot rotate a root body directly"
    if norm(axis) <= 1e-6
        axis = SVector{3, Float64}(0.0, 1.0, 0.0) # Default axis if norm is too small
    end
    axis = axis / norm(axis) # Normalize the rotation axis
    # Update the orientation of the body
    body.q .= SVector{4, Float64}([axis .* sin(θ / 2); cos(θ / 2)]) # Convert Euler angles to quaternion
end

end # module Kinematics
