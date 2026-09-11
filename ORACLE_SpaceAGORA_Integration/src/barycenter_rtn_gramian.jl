module BarycenterRTNGramian

using LinearAlgebra

export barycenter_state,
       rtn_basis,
       barycenter_rtn_frame,
       rtn_pos_cols,
       rtn_vel_cols,
       position_projection_barycenter_frame,
       equal_mass_centering_matrix,
       mass_weighted_centering_matrix,
       barycenter_relative_shape_projection,
       c2d_zoh,
       terminal_output_gramian,
       barycenter_relative_gramian_eigenstructure,
       split_barycenter_mode_blocks,
       nonzero_gramian_modes,
       mode_barycenter_residual,
       gramian_metrics_from_eigs

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

function barycenter_rtn_frame(
    r_list::AbstractVector{<:AbstractVector},
    v_list::AbstractVector{<:AbstractVector},
    masses::AbstractVector{<:Real},
)
    r_barycenter, v_barycenter = barycenter_state(r_list, v_list, masses)
    return rtn_basis(r_barycenter, v_barycenter), r_barycenter, v_barycenter
end

rtn_pos_cols(sat::Int) = (6 * (sat - 1) + 1):(6 * (sat - 1) + 3)
rtn_vel_cols(sat::Int) = (6 * (sat - 1) + 4):(6 * sat)

function position_projection_barycenter_frame(
    m_list::AbstractVector{<:AbstractMatrix{<:Real}},
    m_barycenter::AbstractMatrix{<:Real},
)::Matrix{Float64}
    nsat = length(m_list)
    size(m_barycenter) == (3, 3) ||
        throw(ArgumentError("Barycenter frame must be 3 by 3."))
    projection = zeros(Float64, 3 * nsat, 6 * nsat)
    for sat in 1:nsat
        size(m_list[sat]) == (3, 3) ||
            throw(ArgumentError("Each satellite frame must be 3 by 3."))
        rows = (3 * (sat - 1) + 1):(3 * sat)
        projection[rows, rtn_pos_cols(sat)] .=
            transpose(m_barycenter) * m_list[sat]
    end
    return projection
end

function equal_mass_centering_matrix(nsat::Int)::Matrix{Float64}
    nsat >= 2 || throw(ArgumentError("At least two satellites are required."))
    one_vector = ones(Float64, nsat)
    return Matrix{Float64}(I, nsat, nsat) -
        (one_vector * transpose(one_vector)) ./ nsat
end

function mass_weighted_centering_matrix(
    masses::AbstractVector{<:Real},
)::Matrix{Float64}
    nsat = length(masses)
    nsat >= 2 || throw(ArgumentError("At least two satellites are required."))
    masses_f = Float64.(masses)
    all(>(0.0), masses_f) || throw(ArgumentError("All masses must be positive."))
    weights = masses_f ./ sum(masses_f)
    return Matrix{Float64}(I, nsat, nsat) -
        ones(Float64, nsat) * transpose(weights)
end

function barycenter_relative_shape_projection(
    m_list::AbstractVector{<:AbstractMatrix{<:Real}},
    m_barycenter::AbstractMatrix{<:Real};
    masses = nothing,
)::Matrix{Float64}
    nsat = length(m_list)
    position_projection =
        position_projection_barycenter_frame(m_list, m_barycenter)
    centering = isnothing(masses) ?
        equal_mass_centering_matrix(nsat) :
        mass_weighted_centering_matrix(masses)
    return kron(centering, Matrix{Float64}(I, 3, 3)) * position_projection
end

function c2d_zoh(
    a::AbstractMatrix{<:Real},
    b::AbstractMatrix{<:Real},
    dt::Real,
)
    dt >= 0.0 || throw(ArgumentError("Discrete time step must be nonnegative."))
    a_f = Matrix{Float64}(a)
    b_f = Matrix{Float64}(b)
    state_dimension = size(a_f, 1)
    input_dimension = size(b_f, 2)
    size(a_f, 2) == state_dimension || throw(ArgumentError("A must be square."))
    size(b_f, 1) == state_dimension ||
        throw(ArgumentError("B row count must match A."))

    augmented = zeros(Float64, state_dimension + input_dimension, state_dimension + input_dimension)
    augmented[1:state_dimension, 1:state_dimension] .= a_f
    augmented[1:state_dimension, (state_dimension + 1):end] .= b_f
    exponential = exp(float(dt) .* augmented)
    return (
        exponential[1:state_dimension, 1:state_dimension],
        exponential[1:state_dimension, (state_dimension + 1):end],
    )
end

function terminal_output_gramian(
    ad_list::AbstractVector{<:AbstractMatrix{<:Real}},
    bd_list::AbstractVector{<:AbstractMatrix{<:Real}},
    p_terminal::AbstractMatrix{<:Real},
    start_index::Int,
    horizon::Int;
    ru_inv_list = nothing,
)
    horizon >= 1 || throw(ArgumentError("Horizon must be at least one interval."))
    terminal_index = start_index + horizon
    start_index >= 1 || throw(ArgumentError("Start index must be positive."))
    terminal_index - 1 <= length(ad_list) ||
        throw(ArgumentError("Gramian window exceeds A history."))
    terminal_index - 1 <= length(bd_list) ||
        throw(ArgumentError("Gramian window exceeds B history."))

    s_next = Matrix{Float64}(p_terminal)
    gramian = zeros(Float64, size(s_next, 1), size(s_next, 1))
    for interval in (terminal_index - 1):-1:start_index
        gamma = s_next * bd_list[interval]
        if isnothing(ru_inv_list)
            gramian .+= gamma * transpose(gamma)
        else
            gramian .+= gamma * ru_inv_list[interval] * transpose(gamma)
        end
        s_next = s_next * ad_list[interval]
    end
    return 0.5 .* (gramian .+ transpose(gramian))
end

function barycenter_relative_gramian_eigenstructure(
    gramian::AbstractMatrix{<:Real};
    regularization::Real = 0.0,
)
    symmetric_gramian =
        0.5 .* (Matrix{Float64}(gramian) .+ transpose(Matrix{Float64}(gramian)))
    if regularization > 0.0
        symmetric_gramian .+= float(regularization) .* Matrix{Float64}(
            I,
            size(symmetric_gramian, 1),
            size(symmetric_gramian, 1),
        )
    end
    decomposition = eigen(Symmetric(symmetric_gramian))
    order = sortperm(decomposition.values; rev = true)
    return decomposition.values[order], decomposition.vectors[:, order]
end

function nonzero_gramian_modes(
    gramian::AbstractMatrix{<:Real};
    rtol::Real = 1e-10,
)
    lambda_all, vectors_all =
        barycenter_relative_gramian_eigenstructure(gramian)
    isempty(lambda_all) &&
        return lambda_all, vectors_all, lambda_all, vectors_all
    threshold = float(rtol) * maximum(abs, lambda_all)
    keep = findall(value -> value > threshold, lambda_all)
    return (
        lambda_all[keep],
        vectors_all[:, keep],
        lambda_all,
        vectors_all,
    )
end

function split_barycenter_mode_blocks(
    mode::AbstractVector{<:Real},
    nsat::Int,
)
    length(mode) == 3 * nsat ||
        throw(ArgumentError("Mode length must equal 3Ns."))
    mode_f = Float64.(mode)
    return [mode_f[(3 * (sat - 1) + 1):(3 * sat)] for sat in 1:nsat]
end

function mode_barycenter_residual(
    blocks::AbstractVector{<:AbstractVector};
    masses = nothing,
)
    nsat = length(blocks)
    if isnothing(masses)
        return sum(Float64.(block) for block in blocks) ./ nsat
    end
    length(masses) == nsat ||
        throw(ArgumentError("Mass vector must have one entry per block."))
    masses_f = Float64.(masses)
    all(>(0.0), masses_f) || throw(ArgumentError("All masses must be positive."))
    residual = zeros(Float64, 3)
    for sat in 1:nsat
        residual .+= masses_f[sat] .* blocks[sat]
    end
    return residual ./ sum(masses_f)
end

function gramian_metrics_from_eigs(
    lambda_nonzero::AbstractVector{<:Real};
    eps::Real = 1e-30,
)
    isempty(lambda_nonzero) &&
        return (
            lambda_max = 0.0,
            lambda_min_nonzero = 0.0,
            logdet_nonzero = -Inf,
            condition_nonzero = Inf,
        )
    clipped = max.(Float64.(lambda_nonzero), float(eps))
    return (
        lambda_max = maximum(lambda_nonzero),
        lambda_min_nonzero = minimum(lambda_nonzero),
        logdet_nonzero = sum(log, clipped),
        condition_nonzero = maximum(clipped) / minimum(clipped),
    )
end

end
