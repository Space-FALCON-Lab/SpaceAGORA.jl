using LinearAlgebra

function barycenter_free_basis(nsat::Int)
    nsat >= 2 || throw(ArgumentError("At least two satellites are required."))
    onevec = ones(Float64, 1, nsat)
    factors = svd(onevec; full = true)
    q = factors.V[:, 2:nsat]
    return Matrix{Float64}(q)
end

rtn_pos_cols(sat::Int) = (6 * (sat - 1) + 1):(6 * (sat - 1) + 3)
rtn_vel_cols(sat::Int) = (6 * (sat - 1) + 4):(6 * sat)

function position_inertial_projection_from_rtn(m_list::Vector{Matrix{Float64}})
    nsat = length(m_list)
    p = zeros(Float64, 3 * nsat, 6 * nsat)

    for sat in 1:nsat
        rows = (3 * (sat - 1) + 1):(3 * sat)
        p[rows, rtn_pos_cols(sat)] .= m_list[sat]
    end

    return p
end

function shape_position_projection_inertial(m_list::Vector{Matrix{Float64}})
    nsat = length(m_list)
    q = barycenter_free_basis(nsat)
    ppos = position_inertial_projection_from_rtn(m_list)
    return kron(transpose(q), Matrix{Float64}(I, 3, 3)) * ppos
end

function position_ref_projection_from_rtn(
    m_list::Vector{Matrix{Float64}},
    m_ref::Matrix{Float64},
)
    nsat = length(m_list)
    p = zeros(Float64, 3 * nsat, 6 * nsat)

    for sat in 1:nsat
        rows = (3 * (sat - 1) + 1):(3 * sat)
        p[rows, rtn_pos_cols(sat)] .= transpose(m_ref) * m_list[sat]
    end

    return p
end

function shape_position_projection_ref(
    m_list::Vector{Matrix{Float64}},
    m_ref::Matrix{Float64},
)
    nsat = length(m_list)
    q = barycenter_free_basis(nsat)
    ppos = position_ref_projection_from_rtn(m_list, m_ref)
    return kron(transpose(q), Matrix{Float64}(I, 3, 3)) * ppos
end

function c2d_zoh(a::Matrix{Float64}, b::Matrix{Float64}, dt::Real)
    dt > 0.0 || throw(ArgumentError("Discrete time step must be positive."))
    n = size(a, 1)
    m = size(b, 2)
    size(a, 2) == n || throw(ArgumentError("A must be square."))
    size(b, 1) == n || throw(ArgumentError("B row count must match A."))

    aug = zeros(Float64, n + m, n + m)
    aug[1:n, 1:n] .= a
    aug[1:n, (n + 1):(n + m)] .= b
    exp_aug = exp(float(dt) .* aug)

    return exp_aug[1:n, 1:n], exp_aug[1:n, (n + 1):(n + m)]
end

c2d_zoh_noinput(a::Matrix{Float64}, dt::Real) = exp(float(dt) .* a)

function terminal_output_gramian(
    ad_list::Vector{Matrix{Float64}},
    bd_list::Vector{Matrix{Float64}},
    p_terminal::Matrix{Float64},
    k0::Int,
    horizon::Int;
    ru_inv_list::Union{Nothing, Vector{Matrix{Float64}}} = nothing,
)
    horizon >= 1 || throw(ArgumentError("Horizon must be at least one interval."))
    kf = k0 + horizon
    k0 >= 1 && kf - 1 <= length(ad_list) ||
        throw(ArgumentError("Gramian window is outside the discrete system history."))

    ny = size(p_terminal, 1)
    w = zeros(Float64, ny, ny)
    s_next = copy(p_terminal)

    for step in (kf - 1):-1:k0
        gamma = s_next * bd_list[step]
        if isnothing(ru_inv_list)
            w .+= gamma * transpose(gamma)
        else
            w .+= gamma * ru_inv_list[step] * transpose(gamma)
        end
        s_next = s_next * ad_list[step]
    end

    return 0.5 .* (w .+ transpose(w))
end

function gramian_eigenstructure(w::Matrix{Float64}; regularization::Real = 0.0)
    wsym = 0.5 .* (w .+ transpose(w))
    if regularization > 0.0
        wsym .+= float(regularization) .* Matrix{Float64}(I, size(wsym, 1), size(wsym, 1))
    end
    eig = eigen(Symmetric(wsym))
    order = sortperm(eig.values; rev = true)
    return eig.values[order], eig.vectors[:, order]
end

function gramian_metrics(lambda::Vector{Float64}; eps::Real = 1e-30)
    clipped = max.(lambda, float(eps))
    return (
        lambda_max = maximum(lambda),
        lambda_min = minimum(lambda),
        logdet = sum(log.(clipped)),
        condition = maximum(clipped) / minimum(clipped),
    )
end

function reconstruct_full_shape_from_Q(v_shape::Vector{Float64}, q::Matrix{Float64})
    nsat = size(q, 1)
    full = kron(q, Matrix{Float64}(I, 3, 3)) * v_shape
    blocks = [full[(3 * (sat - 1) + 1):(3 * sat)] for sat in 1:nsat]
    return full, blocks
end
