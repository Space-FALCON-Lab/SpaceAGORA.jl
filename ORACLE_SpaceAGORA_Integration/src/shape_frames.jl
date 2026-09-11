module ShapeFrames

using LinearAlgebra

export barycenter_free_basis,
       mass_weighted_shape_basis,
       rtn_basis,
       barycenter_state,
       common_reference_frame,
       rtn_pos_cols,
       rtn_vel_cols,
       position_ref_projection_from_rtn,
       shape_position_projection_ref,
       reconstruct_full_shape_from_Q,
       edge_incidence_matrix,
       edge_relative_position_projection_ref

function barycenter_free_basis(nsat::Int)::Matrix{Float64}
    nsat >= 2 || throw(ArgumentError("At least two satellites are required."))
    factors = svd(ones(Float64, 1, nsat); full = true)
    return Matrix(factors.V[:, 2:nsat])
end

function mass_weighted_shape_basis(
    masses::AbstractVector{<:Real},
)::Matrix{Float64}
    nsat = length(masses)
    nsat >= 2 || throw(ArgumentError("At least two satellites are required."))
    masses_f = Float64.(masses)
    all(>(0.0), masses_f) || throw(ArgumentError("All masses must be positive."))
    weights = masses_f ./ sum(masses_f)
    factors = svd(reshape(weights, 1, nsat); full = true)
    return Matrix(factors.V[:, 2:nsat])
end

function rtn_basis(
    r::AbstractVector{<:Real},
    v::AbstractVector{<:Real},
)::Matrix{Float64}
    r_f = Float64.(r)
    v_f = Float64.(v)
    radius = norm(r_f)
    radius > 0.0 || throw(ArgumentError("Cannot form RTN frame for zero position."))

    r_hat = r_f ./ radius
    angular_momentum = cross(r_f, v_f)
    momentum_norm = norm(angular_momentum)
    momentum_norm > 0.0 ||
        throw(ArgumentError("Cannot form RTN frame for zero angular momentum."))
    n_hat = angular_momentum ./ momentum_norm
    t_hat = cross(n_hat, r_hat)
    return hcat(r_hat, t_hat, n_hat)
end

function barycenter_state(
    r_list::AbstractVector{<:AbstractVector},
    v_list::AbstractVector{<:AbstractVector},
    masses::AbstractVector{<:Real},
)
    nsat = length(r_list)
    length(v_list) == nsat ||
        throw(ArgumentError("Position and velocity lists must have equal length."))
    length(masses) == nsat ||
        throw(ArgumentError("Mass vector must have one entry per satellite."))

    masses_f = Float64.(masses)
    all(>(0.0), masses_f) || throw(ArgumentError("All masses must be positive."))
    total_mass = sum(masses_f)
    r_barycenter = zeros(Float64, 3)
    v_barycenter = zeros(Float64, 3)
    for sat in 1:nsat
        length(r_list[sat]) == 3 || throw(ArgumentError("Positions must be 3-vectors."))
        length(v_list[sat]) == 3 || throw(ArgumentError("Velocities must be 3-vectors."))
        r_barycenter .+= masses_f[sat] .* r_list[sat]
        v_barycenter .+= masses_f[sat] .* v_list[sat]
    end
    return r_barycenter ./ total_mass, v_barycenter ./ total_mass
end

function common_reference_frame(
    r_list::AbstractVector{<:AbstractVector},
    v_list::AbstractVector{<:AbstractVector};
    frame_type::Symbol = :barycenter_rtn,
    ref_sat::Int = 1,
    masses = ones(length(r_list)),
)::Matrix{Float64}
    nsat = length(r_list)
    length(v_list) == nsat ||
        throw(ArgumentError("Position and velocity lists must have equal length."))

    if frame_type == :barycenter_rtn
        r_barycenter, v_barycenter = barycenter_state(r_list, v_list, masses)
        return rtn_basis(r_barycenter, v_barycenter)
    elseif frame_type == :chief_rtn
        1 <= ref_sat <= nsat || throw(ArgumentError("Reference satellite is out of bounds."))
        return rtn_basis(r_list[ref_sat], v_list[ref_sat])
    elseif frame_type == :inertial
        return Matrix{Float64}(I, 3, 3)
    end
    throw(
        ArgumentError(
            "Unknown frame type $(frame_type); use :barycenter_rtn, :chief_rtn, or :inertial.",
        ),
    )
end

rtn_pos_cols(sat::Int) = (6 * (sat - 1) + 1):(6 * (sat - 1) + 3)
rtn_vel_cols(sat::Int) = (6 * (sat - 1) + 4):(6 * sat)

function position_ref_projection_from_rtn(
    m_list::AbstractVector{<:AbstractMatrix{<:Real}},
    m_ref::AbstractMatrix{<:Real},
)::Matrix{Float64}
    nsat = length(m_list)
    size(m_ref) == (3, 3) || throw(ArgumentError("Reference frame must be 3 by 3."))
    projection = zeros(Float64, 3 * nsat, 6 * nsat)
    for sat in 1:nsat
        size(m_list[sat]) == (3, 3) ||
            throw(ArgumentError("Each satellite frame must be 3 by 3."))
        rows = (3 * (sat - 1) + 1):(3 * sat)
        projection[rows, rtn_pos_cols(sat)] .= transpose(m_ref) * m_list[sat]
    end
    return projection
end

function shape_position_projection_ref(
    m_list::AbstractVector{<:AbstractMatrix{<:Real}},
    m_ref::AbstractMatrix{<:Real},
    q::AbstractMatrix{<:Real},
)::Matrix{Float64}
    nsat = length(m_list)
    size(q) == (nsat, nsat - 1) ||
        throw(ArgumentError("Q must have size Ns by (Ns-1)."))
    position_projection = position_ref_projection_from_rtn(m_list, m_ref)
    return kron(transpose(Matrix{Float64}(q)), Matrix{Float64}(I, 3, 3)) *
        position_projection
end

function reconstruct_full_shape_from_Q(
    v_shape::AbstractVector{<:Real},
    q::AbstractMatrix{<:Real},
)
    nsat = size(q, 1)
    size(q, 2) == nsat - 1 || throw(ArgumentError("Q must have Ns-1 columns."))
    length(v_shape) == 3 * (nsat - 1) ||
        throw(ArgumentError("Shape vector must have length 3(Ns-1)."))
    full = kron(Matrix{Float64}(q), Matrix{Float64}(I, 3, 3)) * Float64.(v_shape)
    blocks = [full[(3 * (sat - 1) + 1):(3 * sat)] for sat in 1:nsat]
    return full, blocks
end

function edge_incidence_matrix(
    nsat::Int,
    edges::AbstractVector{<:Tuple{Int, Int}},
)::Matrix{Float64}
    nsat >= 2 || throw(ArgumentError("At least two satellites are required."))
    incidence = zeros(Float64, length(edges), nsat)
    for (edge_index, (start_sat, end_sat)) in enumerate(edges)
        1 <= start_sat <= nsat || throw(ArgumentError("Edge start is out of bounds."))
        1 <= end_sat <= nsat || throw(ArgumentError("Edge end is out of bounds."))
        start_sat != end_sat || throw(ArgumentError("Self-edges are not allowed."))
        incidence[edge_index, start_sat] = -1.0
        incidence[edge_index, end_sat] = 1.0
    end
    return incidence
end

function edge_relative_position_projection_ref(
    m_list::AbstractVector{<:AbstractMatrix{<:Real}},
    m_ref::AbstractMatrix{<:Real},
    edges::AbstractVector{<:Tuple{Int, Int}},
)::Matrix{Float64}
    incidence = edge_incidence_matrix(length(m_list), edges)
    position_projection = position_ref_projection_from_rtn(m_list, m_ref)
    return kron(incidence, Matrix{Float64}(I, 3, 3)) * position_projection
end

end
