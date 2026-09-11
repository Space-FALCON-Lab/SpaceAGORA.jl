module KSNormalizedBarycentricGramian

using LinearAlgebra

export normalized_physical_rtn_system,
       KSScales,
       make_ks_scales,
       lambda_ks,
       S_ks,
       normalized_ks_to_cartesian,
       normalized_ks_rhs,
       normalize_ks_state,
       denormalize_ks_state,
       stack_state_scale,
       stack_pq_scale,
       normalize_ks_linearization_s,
       normalize_pq_linearization_s,
       sbar_to_tbar_left_scaling_pq,
       normalized_H_pq,
       stacked_normalized_H_pq,
       physical_reduction_basis,
       reduce_pq_dynamics,
       rtn_basis,
       blockdiag_rtn_rotation,
       normalized_rtn_system,
       barycenter_state,
       barycenter_rtn_frame,
       centering_matrix,
       position_projection_barycenter_frame,
       barycentric_output_projection,
       c2d_zoh_normalized_time,
       terminal_output_gramian,
       gramian_shape_eigenstructure,
       split_mode_blocks

struct KSScales
    L::Float64
    mu::Float64
    T::Float64
    V::Float64
    a_scale::Float64
    p_scale::Float64
    q_scale::Float64
    h_scale::Float64
    t_scale::Float64
    s_scale::Float64
    u_scale::Float64
end

function make_ks_scales(; mu::Real, L::Real, u_scale::Real = 1.0)
    mu > 0.0 || throw(ArgumentError("Gravitational parameter must be positive."))
    L > 0.0 || throw(ArgumentError("Characteristic length must be positive."))
    u_scale > 0.0 || throw(ArgumentError("Input scale must be positive."))
    length_scale = float(L)
    mu_scale = float(mu)
    time_scale = sqrt(length_scale^3 / mu_scale)
    velocity_scale = length_scale / time_scale
    return KSScales(
        length_scale,
        mu_scale,
        time_scale,
        velocity_scale,
        length_scale / time_scale^2,
        sqrt(length_scale),
        length_scale^(3 / 2) / time_scale,
        velocity_scale^2,
        time_scale,
        time_scale / length_scale,
        float(u_scale),
    )
end

function blockdiag_dense(blocks::AbstractVector{<:AbstractMatrix})
    nrows = sum(size(block, 1) for block in blocks)
    ncols = sum(size(block, 2) for block in blocks)
    result = zeros(Float64, nrows, ncols)
    row_start = 1
    col_start = 1
    for block in blocks
        nrow, ncol = size(block)
        result[row_start:(row_start + nrow - 1), col_start:(col_start + ncol - 1)] .= block
        row_start += nrow
        col_start += ncol
    end
    return result
end

repeat_blockdiag(block, count) =
    blockdiag_dense([Matrix{Float64}(block) for _ in 1:count])

function lambda_ks(p::AbstractVector{<:Real})
    length(p) == 4 || throw(ArgumentError("KS position must have four components."))
    p1, p2, p3, p4 = p
    return Float64[
        p1 -p2 -p3 p4
        p2 p1 -p4 -p3
        p3 p4 p1 p2
    ]
end

function S_ks(p::AbstractVector{<:Real})
    length(p) == 4 || throw(ArgumentError("KS position must have four components."))
    p1, p2, p3, p4 = p
    l_matrix = Float64[
        p1 -p2 -p3 p4
        p2 p1 -p4 -p3
        p3 p4 p1 p2
        p4 -p3 p2 -p1
    ]
    return transpose(l_matrix)
end

function normalized_ks_to_cartesian(
    pbar::AbstractVector{<:Real},
    qbar::AbstractVector{<:Real},
)
    length(pbar) == 4 || throw(ArgumentError("Normalized p must have four components."))
    length(qbar) == 4 || throw(ArgumentError("Normalized q must have four components."))
    pbar_f = Float64.(pbar)
    qbar_f = Float64.(qbar)
    rhobar = dot(pbar_f, pbar_f)
    rhobar > 0.0 || throw(ArgumentError("Normalized KS radius must be positive."))
    lambda_p = lambda_ks(pbar_f)
    return lambda_p * pbar_f, (2.0 / rhobar) .* (lambda_p * qbar_f), rhobar
end

function normalized_ks_rhs(
    xbar::AbstractVector{<:Real},
    ubar::AbstractVector{<:Real},
    accel_bar_callback::Function,
)
    length(xbar) == 10 || throw(ArgumentError("Normalized KS state must have ten components."))
    pbar = Float64.(xbar[1:4])
    qbar = Float64.(xbar[5:8])
    hbar = float(xbar[9])
    _, vbar, rhobar = normalized_ks_to_cartesian(pbar, qbar)
    accel_bar = Float64.(accel_bar_callback(xbar, ubar))
    length(accel_bar) == 3 || throw(ArgumentError("Acceleration callback must return three components."))

    derivative = zeros(Float64, 10)
    derivative[1:4] .= qbar
    derivative[5:8] .=
        -(hbar / 4.0) .* pbar .+
        (rhobar / 2.0) .* (S_ks(pbar) * vcat(accel_bar, 0.0))
    derivative[9] = -2.0 * rhobar * dot(vbar, accel_bar)
    derivative[10] = rhobar
    return derivative
end

function single_state_scale(scales::KSScales)
    return Matrix(
        Diagonal(
            vcat(
                fill(scales.p_scale, 4),
                fill(scales.q_scale, 4),
                scales.h_scale,
                scales.t_scale,
            ),
        ),
    )
end

function single_pq_scale(scales::KSScales)
    return Matrix(
        Diagonal(
            vcat(
                fill(scales.p_scale, 4),
                fill(scales.q_scale, 4),
            ),
        ),
    )
end

stack_state_scale(scales::KSScales, nsat::Int) =
    repeat_blockdiag(single_state_scale(scales), nsat)
stack_pq_scale(scales::KSScales, nsat::Int) =
    repeat_blockdiag(single_pq_scale(scales), nsat)

function normalize_ks_state(x::AbstractVector{<:Real}, scales::KSScales, nsat::Int)
    length(x) == 10 * nsat || throw(ArgumentError("KS state length must equal 10Ns."))
    return stack_state_scale(scales, nsat) \ Float64.(x)
end

function denormalize_ks_state(xbar::AbstractVector{<:Real}, scales::KSScales, nsat::Int)
    length(xbar) == 10 * nsat || throw(ArgumentError("Normalized KS state length must equal 10Ns."))
    return stack_state_scale(scales, nsat) * Float64.(xbar)
end

function input_scale_matrix(input_dimension, scales, input_scale)
    if isnothing(input_scale)
        return scales.u_scale .* Matrix{Float64}(I, input_dimension, input_dimension)
    end
    matrix = Matrix{Float64}(input_scale)
    size(matrix) == (input_dimension, input_dimension) ||
        throw(ArgumentError("Input scale must be square in the input dimension."))
    return matrix
end

"""Normalize the full physical KS tangent before its analytical RTN push-forward."""
function normalized_physical_rtn_system(physical, scales, nsat; input_scale=nothing)
    Sy = stack_state_scale(scales, nsat)
    Sx = repeat_blockdiag(Matrix(Diagonal(vcat(fill(scales.L, 3), fill(scales.V, 3)))), nsat)
    Su = input_scale_matrix(size(physical.B_time, 2), scales, input_scale)
    At = scales.T * (Sy \ (physical.A_time * Sy))
    Bt = scales.T * (Sy \ (physical.B_time * Su))
    H = Sx \ (physical.H * Sy)
    G = Sy \ (physical.G * Sx)
    Hd = scales.T * (Sx \ (physical.Hdot * Sy))
    R, Rd = physical.R, scales.T * physical.Rdot
    Ac = (Hd + H*At)*G
    return (A_RTN=Rd*R' + R*Ac*R', B_RTN=R*H*Bt, H=H, G=G,
        state_scale=Sx, input_scale=Su)
end

function normalize_ks_linearization_s(
    a_ks_s::AbstractMatrix{<:Real},
    b_ks_s::AbstractMatrix{<:Real},
    scales::KSScales,
    nsat::Int;
    input_scale = nothing,
)
    state_scale = stack_state_scale(scales, nsat)
    size(a_ks_s) == (10 * nsat, 10 * nsat) ||
        throw(ArgumentError("Full KS A matrix must be 10Ns by 10Ns."))
    size(b_ks_s, 1) == 10 * nsat ||
        throw(ArgumentError("Full KS B row count must equal 10Ns."))
    u_scale = input_scale_matrix(size(b_ks_s, 2), scales, input_scale)
    return (
        scales.s_scale .* (state_scale \ (Matrix{Float64}(a_ks_s) * state_scale)),
        scales.s_scale .* (state_scale \ (Matrix{Float64}(b_ks_s) * u_scale)),
    )
end

function normalize_pq_linearization_s(
    a_pq_s::AbstractMatrix{<:Real},
    b_pq_s::AbstractMatrix{<:Real},
    scales::KSScales,
    nsat::Int;
    input_scale = nothing,
)
    pq_scale = stack_pq_scale(scales, nsat)
    size(a_pq_s) == (8 * nsat, 8 * nsat) ||
        throw(ArgumentError("KS p/q A matrix must be 8Ns by 8Ns."))
    size(b_pq_s, 1) == 8 * nsat ||
        throw(ArgumentError("KS p/q B row count must equal 8Ns."))
    u_scale = input_scale_matrix(size(b_pq_s, 2), scales, input_scale)
    return (
        scales.s_scale .* (pq_scale \ (Matrix{Float64}(a_pq_s) * pq_scale)),
        scales.s_scale .* (pq_scale \ (Matrix{Float64}(b_pq_s) * u_scale)),
    )
end

function sbar_to_tbar_left_scaling_pq(
    a_pq_sbar::AbstractMatrix{<:Real},
    b_pq_sbar::AbstractMatrix{<:Real},
    pbar_list::AbstractVector{<:AbstractVector},
)
    blocks = Matrix{Float64}[]
    for pbar in pbar_list
        rhobar = dot(pbar, pbar)
        rhobar > 0.0 || throw(ArgumentError("Normalized KS radius must be positive."))
        push!(blocks, (1.0 / rhobar) .* Matrix{Float64}(I, 8, 8))
    end
    time_scaling = blockdiag_dense(blocks)
    size(a_pq_sbar) == size(time_scaling) ||
        throw(ArgumentError("Normalized p/q A dimensions do not match p list."))
    return time_scaling * a_pq_sbar, time_scaling * b_pq_sbar
end

function normalized_H_pq(
    pbar::AbstractVector{<:Real},
    qbar::AbstractVector{<:Real},
)
    pbar_f = Float64.(pbar)
    qbar_f = Float64.(qbar)
    rhobar = dot(pbar_f, pbar_f)
    rhobar > 0.0 || throw(ArgumentError("Normalized KS radius must be positive."))
    lambda_p = lambda_ks(pbar_f)
    z = lambda_p * qbar_f
    return [
        2.0 .* lambda_p zeros(Float64, 3, 4)
        (2.0 / rhobar) .* lambda_ks(qbar_f) .-
        (4.0 / rhobar^2) .* (z * transpose(pbar_f)) (2.0 / rhobar) .* lambda_p
    ]
end

function stacked_normalized_H_pq(
    pbar_list::AbstractVector{<:AbstractVector},
    qbar_list::AbstractVector{<:AbstractVector},
)
    length(pbar_list) == length(qbar_list) ||
        throw(ArgumentError("Normalized p and q lists must have equal length."))
    return blockdiag_dense(
        [normalized_H_pq(pbar_list[i], qbar_list[i]) for i in eachindex(pbar_list)],
    )
end

function physical_reduction_basis(hbar::AbstractMatrix{<:Real}; rtol::Real = 1e-9)
    factorization = svd(Matrix{Float64}(hbar))
    isempty(factorization.S) && return zeros(Float64, size(hbar, 2), 0), 0, Float64[]
    maximum_singular_value = maximum(factorization.S)
    rank_physical = count(>(float(rtol) * maximum_singular_value), factorization.S)
    return factorization.V[:, 1:rank_physical], rank_physical, factorization.S
end

function reduce_pq_dynamics(a_pq_tbar, b_pq_tbar, basis)
    basis_f = Matrix{Float64}(basis)
    return (
        transpose(basis_f) * a_pq_tbar * basis_f,
        transpose(basis_f) * b_pq_tbar,
    )
end

function rtn_basis(r::AbstractVector{<:Real}, v::AbstractVector{<:Real})
    r_f = Float64.(r)
    v_f = Float64.(v)
    norm(r_f) > 0.0 || throw(ArgumentError("Position norm must be positive."))
    angular_momentum = cross(r_f, v_f)
    norm(angular_momentum) > 0.0 || throw(ArgumentError("Angular momentum norm must be positive."))
    r_hat = r_f ./ norm(r_f)
    n_hat = angular_momentum ./ norm(angular_momentum)
    return hcat(r_hat, cross(n_hat, r_hat), n_hat)
end

function blockdiag_rtn_rotation(frames::AbstractVector{<:AbstractMatrix})
    return blockdiag_dense(
        [
            [
                transpose(Matrix{Float64}(frame)) zeros(Float64, 3, 3)
                zeros(Float64, 3, 3) transpose(Matrix{Float64}(frame))
            ] for frame in frames
        ],
    )
end

function normalized_rtn_system(a_eta, b_eta, tbar, tbar_dot)
    tbar_f = Matrix{Float64}(tbar)
    size(tbar_f, 1) == size(tbar_f, 2) ||
        throw(ArgumentError("Normalized RTN transform must be square."))
    inverse_tbar = inv(tbar_f)
    return (
        Matrix{Float64}(tbar_dot) * inverse_tbar +
        tbar_f * Matrix{Float64}(a_eta) * inverse_tbar,
        tbar_f * Matrix{Float64}(b_eta),
    )
end

function barycenter_state(r_list, v_list, masses)
    nsat = length(r_list)
    length(v_list) == nsat == length(masses) ||
        throw(ArgumentError("Positions, velocities, and masses must have equal length."))
    masses_f = Float64.(masses)
    all(>(0.0), masses_f) || throw(ArgumentError("All masses must be positive."))
    total_mass = sum(masses_f)
    r_barycenter = sum(masses_f[i] .* Float64.(r_list[i]) for i in 1:nsat) ./ total_mass
    v_barycenter = sum(masses_f[i] .* Float64.(v_list[i]) for i in 1:nsat) ./ total_mass
    return r_barycenter, v_barycenter
end

function barycenter_rtn_frame(r_list, v_list, masses)
    r_barycenter, v_barycenter = barycenter_state(r_list, v_list, masses)
    return rtn_basis(r_barycenter, v_barycenter), r_barycenter, v_barycenter
end

function centering_matrix(nsat::Int)
    nsat >= 2 || throw(ArgumentError("At least two satellites are required."))
    one_vector = ones(Float64, nsat)
    return Matrix{Float64}(I, nsat, nsat) -
        one_vector * transpose(one_vector) ./ nsat
end

function centering_matrix(masses::AbstractVector{<:Real})
    length(masses) >= 2 || throw(ArgumentError("At least two satellites are required."))
    masses_f = Float64.(masses)
    all(>(0.0), masses_f) || throw(ArgumentError("All masses must be positive."))
    weights = masses_f ./ sum(masses_f)
    return Matrix{Float64}(I, length(masses), length(masses)) -
        ones(Float64, length(masses)) * transpose(weights)
end

rtn_pos_cols(sat::Int) = (6 * (sat - 1) + 1):(6 * (sat - 1) + 3)

function position_projection_barycenter_frame(frames, barycenter_frame)
    nsat = length(frames)
    projection = zeros(Float64, 3 * nsat, 6 * nsat)
    for sat in 1:nsat
        rows = (3 * (sat - 1) + 1):(3 * sat)
        projection[rows, rtn_pos_cols(sat)] .=
            transpose(Matrix{Float64}(barycenter_frame)) * frames[sat]
    end
    return projection
end

function barycentric_output_projection(frames, barycenter_frame; masses = nothing)
    nsat = length(frames)
    centering = isnothing(masses) ? centering_matrix(nsat) : centering_matrix(masses)
    return kron(centering, Matrix{Float64}(I, 3, 3)) *
        position_projection_barycenter_frame(frames, barycenter_frame)
end

function c2d_zoh_normalized_time(a_bar, b_bar, dt_physical, scales::KSScales)
    dt_physical >= 0.0 || throw(ArgumentError("Physical time step must be nonnegative."))
    a_f = Matrix{Float64}(a_bar)
    b_f = Matrix{Float64}(b_bar)
    nstate = size(a_f, 1)
    ninput = size(b_f, 2)
    size(a_f, 2) == nstate || throw(ArgumentError("A must be square."))
    size(b_f, 1) == nstate || throw(ArgumentError("B row count must match A."))
    augmented = zeros(Float64, nstate + ninput, nstate + ninput)
    augmented[1:nstate, 1:nstate] .= a_f
    augmented[1:nstate, (nstate + 1):end] .= b_f
    exponential = exp((float(dt_physical) / scales.T) .* augmented)
    return (
        exponential[1:nstate, 1:nstate],
        exponential[1:nstate, (nstate + 1):end],
    )
end

function terminal_output_gramian(
    ad_list,
    bd_list,
    terminal_projection,
    start_index::Int,
    horizon::Int;
    energy_weighted::Bool = false,
    dtbar_list = nothing,
)
    horizon >= 1 || throw(ArgumentError("Horizon must be positive."))
    terminal_index = start_index + horizon
    terminal_index - 1 <= length(ad_list) || throw(ArgumentError("Gramian window exceeds A history."))
    terminal_index - 1 <= length(bd_list) || throw(ArgumentError("Gramian window exceeds B history."))
    energy_weighted && isnothing(dtbar_list) &&
        throw(ArgumentError("Normalized time steps are required for energy weighting."))

    s_next = Matrix{Float64}(terminal_projection)
    gramian = zeros(Float64, size(s_next, 1), size(s_next, 1))
    for interval in (terminal_index - 1):-1:start_index
        gamma = s_next * bd_list[interval]
        weight = energy_weighted ? 1.0 / dtbar_list[interval] : 1.0
        gramian .+= weight .* (gamma * transpose(gamma))
        s_next = s_next * ad_list[interval]
    end
    return 0.5 .* (gramian .+ transpose(gramian))
end

function gramian_shape_eigenstructure(
    gramian::AbstractMatrix{<:Real};
    nsat::Int,
    eig_tol::Real = 1e-10,
)
    gramian_f = Matrix{Float64}(gramian)
    size(gramian_f) == (3 * nsat, 3 * nsat) ||
        throw(ArgumentError("Barycentric Gramian must be 3Ns by 3Ns."))
    decomposition = eigen(Symmetric(0.5 .* (gramian_f .+ transpose(gramian_f))))
    order = sortperm(decomposition.values; rev = true)
    lambda_all = max.(decomposition.values[order], 0.0)
    vectors_all = decomposition.vectors[:, order]
    expected_rank = 3 * (nsat - 1)
    if isempty(lambda_all) || lambda_all[1] == 0.0
        return (
            lambda_shape = Float64[],
            vectors_shape = zeros(Float64, 3 * nsat, 0),
            lambda_all = lambda_all,
            vectors_all = vectors_all,
            sqrt_lambda = Float64[],
            lambda_relative = Float64[],
            participation = Float64[],
            condition_shape = Inf,
            dominance_fraction = 0.0,
            effective_rank = 0.0,
            rank_shape = 0,
            expected_rank = expected_rank,
        )
    end
    keep = findall(>(float(eig_tol) * lambda_all[1]), lambda_all)
    length(keep) > expected_rank && (keep = keep[1:expected_rank])
    lambda_shape = lambda_all[keep]
    vectors_shape = vectors_all[:, keep]
    participation = lambda_shape ./ sum(lambda_shape)
    positive_participation = participation[participation .> 0.0]
    return (
        lambda_shape = lambda_shape,
        vectors_shape = vectors_shape,
        lambda_all = lambda_all,
        vectors_all = vectors_all,
        sqrt_lambda = sqrt.(lambda_shape),
        lambda_relative = lambda_shape ./ lambda_shape[1],
        participation = participation,
        condition_shape = lambda_shape[1] / lambda_shape[end],
        dominance_fraction = participation[1],
        effective_rank = exp(-sum(positive_participation .* log.(positive_participation))),
        rank_shape = length(lambda_shape),
        expected_rank = expected_rank,
    )
end

function split_mode_blocks(mode::AbstractVector{<:Real})
    length(mode) % 3 == 0 || throw(ArgumentError("Mode length must be divisible by three."))
    return [Float64.(mode[(3 * (sat - 1) + 1):(3 * sat)]) for sat in 1:(length(mode) ÷ 3)]
end

end
