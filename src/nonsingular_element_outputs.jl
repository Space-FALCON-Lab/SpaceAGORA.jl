module NonsingularElementOutputs

using ForwardDiff
using LinearAlgebra

export eccentricity_vector_inertial,
       rtn_basis_from_rv,
       eccentricity_vector_orbital_frame,
       inclination_vector,
       semimajor_axis_from_rv,
       circular_tracking_elements_from_rv,
       inertial_eccentricity_tracking_elements_from_rv,
       circular_tracking_jacobian_inertial,
       inertial_eccentricity_tracking_jacobian_inertial,
       circular_tracking_sensitivity_rtn,
       inertial_eccentricity_tracking_sensitivity_rtn,
       circular_tracking_jacobian_inertial_ad,
       inertial_eccentricity_tracking_jacobian_inertial_ad,
       circular_tracking_sensitivity_rtn_ad,
       inertial_eccentricity_tracking_sensitivity_rtn_ad,
       embed_target_element_matrix

function eccentricity_vector_inertial(
    position::AbstractVector{T},
    velocity::AbstractVector{T};
    mu::Real = 1.0,
) where {T<:Real}
    length(position) == 3 && length(velocity) == 3 ||
        throw(ArgumentError("Position and velocity must each contain three components."))
    mu > 0.0 || throw(ArgumentError("The gravitational parameter must be positive."))
    mu_t = one(T) * mu
    radius = norm(position)
    return (
        (dot(velocity, velocity) - mu_t / radius) .* position .-
        dot(position, velocity) .* velocity
    ) ./ mu_t
end

function inertial_eccentricity_tracking_elements_from_rv(
    state::AbstractVector{T};
    mu::Real = 1.0,
) where {T<:Real}
    length(state) == 6 ||
        throw(ArgumentError("A Cartesian state must contain three positions and three velocities."))
    position = state[1:3]
    velocity = state[4:6]
    return vcat(
        [semimajor_axis_from_rv(position, velocity; mu = mu)],
        eccentricity_vector_inertial(position, velocity; mu = mu),
    )
end

function rtn_basis_from_rv(
    position::AbstractVector{T},
    velocity::AbstractVector{T},
) where {T<:Real}
    length(position) == 3 && length(velocity) == 3 ||
        throw(ArgumentError("Position and velocity must each contain three components."))
    radial = position / norm(position)
    angular_momentum = cross(position, velocity)
    normal = angular_momentum / norm(angular_momentum)
    tangential = cross(normal, radial)
    return hcat(radial, tangential, normal)
end

function eccentricity_vector_orbital_frame(
    position::AbstractVector{T},
    velocity::AbstractVector{T};
    mu::Real = 1.0,
) where {T<:Real}
    eccentricity_vector = eccentricity_vector_inertial(position, velocity; mu = mu)
    frame = rtn_basis_from_rv(position, velocity)
    return [
        dot(frame[:, 1], eccentricity_vector),
        dot(frame[:, 2], eccentricity_vector),
    ]
end

function inclination_vector(
    position::AbstractVector{T},
    velocity::AbstractVector{T},
) where {T<:Real}
    angular_momentum = cross(position, velocity)
    normal = angular_momentum / norm(angular_momentum)
    return [normal[1], normal[2]]
end

function semimajor_axis_from_rv(
    position::AbstractVector{T},
    velocity::AbstractVector{T};
    mu::Real = 1.0,
) where {T<:Real}
    length(position) == 3 && length(velocity) == 3 ||
        throw(ArgumentError("Position and velocity must each contain three components."))
    mu > 0.0 || throw(ArgumentError("The gravitational parameter must be positive."))
    mu_t = one(T) * mu
    energy = one(T) * 0.5 * dot(velocity, velocity) - mu_t / norm(position)
    return -mu_t / (one(T) * 2.0 * energy)
end

function circular_tracking_elements_from_rv(
    state::AbstractVector{T};
    mu::Real = 1.0,
    include_plane::Bool = false,
) where {T<:Real}
    length(state) == 6 ||
        throw(ArgumentError("A Cartesian state must contain three positions and three velocities."))
    position = state[1:3]
    velocity = state[4:6]
    semimajor_axis = semimajor_axis_from_rv(position, velocity; mu = mu)
    eccentricity_rtn = eccentricity_vector_orbital_frame(position, velocity; mu = mu)
    if include_plane
        return vcat([semimajor_axis], eccentricity_rtn, inclination_vector(position, velocity))
    end
    return vcat([semimajor_axis], eccentricity_rtn)
end

function eccentricity_differentials(r, v; mu)
    skew(x) = [0.0 -x[3] x[2]; x[3] 0.0 -x[1]; -x[2] x[1] 0.0]
    eye = Matrix{Float64}(I, 3, 3)
    dr, dv = hcat(eye, zeros(3, 3)), hcat(zeros(3, 3), eye)
    radius = norm(r)
    h = cross(r, v)
    dh = -skew(v)*dr + skew(r)*dv
    de = (-skew(h)*dv + skew(v)*dh)/mu - (eye/radius - r*r'/radius^3)*dr
    a = semimajor_axis_from_rv(r, v; mu=mu)
    da = 2a^2/mu * (mu*r'/radius^3*dr + v'*dv)
    radial, normal = r/radius, h/norm(h)
    dR = (eye-radial*radial')*dr/radius
    dN = (eye-normal*normal')*dh/norm(h)
    dT = -skew(radial)*dN + skew(normal)*dR
    return (; da, de, dR, dT, dN)
end

function inertial_eccentricity_tracking_jacobian_inertial(r, v; mu=1.0)
    d = eccentricity_differentials(r, v; mu=mu)
    return vcat(d.da, d.de)
end

function circular_tracking_jacobian_inertial(r, v; mu=1.0, include_plane=false)
    d = eccentricity_differentials(r, v; mu=mu)
    frame = rtn_basis_from_rv(r, v)
    ev = eccentricity_vector_inertial(r, v; mu=mu)
    output = vcat(d.da, frame[:, 1]'*d.de + ev'*d.dR,
        frame[:, 2]'*d.de + ev'*d.dT)
    return include_plane ? vcat(output, d.dN[1:2, :]) : output
end

function circular_tracking_sensitivity_rtn(r, v, frame; mu=1.0, include_plane=false)
    rotation = [frame zeros(3, 3); zeros(3, 3) frame]
    return circular_tracking_jacobian_inertial(r, v; mu=mu, include_plane=include_plane) * rotation
end

function inertial_eccentricity_tracking_sensitivity_rtn(r, v, frame; mu=1.0)
    rotation = [frame zeros(3, 3); zeros(3, 3) frame]
    return inertial_eccentricity_tracking_jacobian_inertial(r, v; mu=mu) * rotation
end

# AD variants below are reference checks, not production output matrices.
function circular_tracking_jacobian_inertial_ad(
    position::AbstractVector{<:Real},
    velocity::AbstractVector{<:Real};
    mu::Real = 1.0,
    include_plane::Bool = false,
)
    length(position) == 3 && length(velocity) == 3 ||
        throw(ArgumentError("Position and velocity must each contain three components."))
    state = vcat(Float64.(position), Float64.(velocity))
    return ForwardDiff.jacobian(state) do local_state
        circular_tracking_elements_from_rv(
            local_state;
            mu = mu,
            include_plane = include_plane,
        )
    end
end

function inertial_eccentricity_tracking_jacobian_inertial_ad(
    position::AbstractVector{<:Real},
    velocity::AbstractVector{<:Real};
    mu::Real = 1.0,
)
    length(position) == 3 && length(velocity) == 3 ||
        throw(ArgumentError("Position and velocity must each contain three components."))
    state = vcat(Float64.(position), Float64.(velocity))
    return ForwardDiff.jacobian(state) do local_state
        inertial_eccentricity_tracking_elements_from_rv(local_state; mu = mu)
    end
end

function circular_tracking_sensitivity_rtn_ad(
    position::AbstractVector{<:Real},
    velocity::AbstractVector{<:Real},
    rtn_frame::AbstractMatrix{<:Real};
    mu::Real = 1.0,
    include_plane::Bool = false,
)
    size(rtn_frame) == (3, 3) || throw(ArgumentError("The RTN frame must be 3 by 3."))
    inertial_jacobian = circular_tracking_jacobian_inertial_ad(
        position,
        velocity;
        mu = mu,
        include_plane = include_plane,
    )
    rtn_to_inertial = zeros(Float64, 6, 6)
    rtn_to_inertial[1:3, 1:3] .= rtn_frame
    rtn_to_inertial[4:6, 4:6] .= rtn_frame
    return inertial_jacobian * rtn_to_inertial
end

function inertial_eccentricity_tracking_sensitivity_rtn_ad(
    position::AbstractVector{<:Real},
    velocity::AbstractVector{<:Real},
    rtn_frame::AbstractMatrix{<:Real};
    mu::Real = 1.0,
)
    size(rtn_frame) == (3, 3) || throw(ArgumentError("The RTN frame must be 3 by 3."))
    inertial_jacobian = inertial_eccentricity_tracking_jacobian_inertial_ad(
        position,
        velocity;
        mu = mu,
    )
    rtn_to_inertial = zeros(Float64, 6, 6)
    rtn_to_inertial[1:3, 1:3] .= rtn_frame
    rtn_to_inertial[4:6, 4:6] .= rtn_frame
    return inertial_jacobian * rtn_to_inertial
end

function embed_target_element_matrix(
    satellite_output::AbstractMatrix{<:Real},
    target_satellite::Int,
    satellite_count::Int,
)
    satellite_count >= 1 || throw(ArgumentError("At least one satellite is required."))
    1 <= target_satellite <= satellite_count ||
        throw(ArgumentError("The target satellite index is invalid."))
    size(satellite_output, 2) == 6 ||
        throw(ArgumentError("A satellite output matrix must have six columns."))
    output = zeros(Float64, size(satellite_output, 1), 6 * satellite_count)
    columns = (6 * (target_satellite - 1) + 1):(6 * target_satellite)
    output[:, columns] .= satellite_output
    return output
end

end
