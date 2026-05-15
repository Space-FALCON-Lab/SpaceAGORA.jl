function rpo_hcw_continuous_mats(n::Real)
    n = Float64(n)
    A = [
        0.0 0.0 0.0 1.0 0.0 0.0
        0.0 0.0 0.0 0.0 1.0 0.0
        0.0 0.0 0.0 0.0 0.0 1.0
        3*n^2 0.0 0.0 0.0 2*n 0.0
        0.0 0.0 0.0 -2*n 0.0 0.0
        0.0 0.0 -n^2 0.0 0.0 0.0
    ]
    B = [
        0.0 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]
    return A, B
end

function rpo_discretize_zoh(A, B, dt::Real)
    nx = size(A, 1)
    nu = size(B, 2)
    M = zeros(nx + nu, nx + nu)
    M[1:nx, 1:nx] .= A
    M[1:nx, nx+1:end] .= B
    Md = exp(M * Float64(dt))
    return Md[1:nx, 1:nx], Md[1:nx, nx+1:end]
end

function rpo_prediction_mats(Ad, Bd, horizon::Int)
    nx = size(Ad, 1)
    nu = size(Bd, 2)
    Abar = zeros(nx * horizon, nx)
    Bbar = zeros(nx * horizon, nu * horizon)
    powers = [Matrix{Float64}(I, nx, nx)]
    for _ in 1:horizon
        push!(powers, Ad * powers[end])
    end
    for i in 1:horizon
        rows = (nx * (i - 1) + 1):(nx * i)
        Abar[rows, :] .= powers[i + 1]
        for j in 1:i
            cols = (nu * (j - 1) + 1):(nu * j)
            Bbar[rows, cols] .= powers[i - j + 1] * Bd
        end
    end
    return Abar, Bbar
end

function rpo_block_diag(blocks)
    rows = sum(size(b, 1) for b in blocks)
    cols = sum(size(b, 2) for b in blocks)
    out = zeros(rows, cols)
    r = 1
    c = 1
    for b in blocks
        nr, nc = size(b)
        out[r:r+nr-1, c:c+nc-1] .= b
        r += nr
        c += nc
    end
    return out
end

mutable struct RpoLQMPCController
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
    horizon::Int
    H::SparseMatrixCSC{Float64, Int}
    E::Matrix{Float64}
    F::Matrix{Float64}
    G::SparseMatrixCSC{Float64, Int}
    W::Vector{Float64}
    qp_model::OSQP.Model
    qp_results::OSQP.Results
    U_prev::Vector{Float64}
    u_min::Vector{Float64}
    u_max::Vector{Float64}
end

function init_rpo_lqmpc(n, dt, Q, R, Qf, horizon; u_min=nothing, u_max=nothing)
    A, B = rpo_hcw_continuous_mats(n)
    Ad, Bd = rpo_discretize_zoh(A, B, dt)
    nx = size(Ad, 1)
    nu = size(Bd, 2)
    u_min = u_min === nothing ? fill(-Inf, nu) : Vector{Float64}(u_min)
    u_max = u_max === nothing ? fill(Inf, nu) : Vector{Float64}(u_max)
    Abar, Bbar = rpo_prediction_mats(Ad, Bd, horizon)
    Qbar = rpo_block_diag([fill(Matrix{Float64}(Q), max(horizon - 1, 0)); [Matrix{Float64}(Qf)]])
    Rbar = rpo_block_diag(fill(Matrix{Float64}(R), horizon))
    H = 0.5 .* (Bbar' * Qbar * Bbar + Rbar)
    H .= H .+ H'
    E = Bbar' * Qbar * Abar
    F = Bbar' * Qbar
    G = sparse(vcat(Matrix{Float64}(I, nu * horizon, nu * horizon), -Matrix{Float64}(I, nu * horizon, nu * horizon)))
    W = vcat(repeat(u_max, horizon), -repeat(u_min, horizon))
    model = OSQP.Model()
    nU = nu * horizon
    OSQP.setup!(
        model;
        P=triu(2.0 .* sparse(H)),
        q=zeros(nU),
        A=G,
        l=fill(-Inf, size(G, 1)),
        u=W,
        verbose=false,
        polish=false,
        warm_start=true,
        eps_abs=1.0e-4,
        eps_rel=1.0e-4,
        max_iter=1000,
    )
    return RpoLQMPCController(Ad, Bd, horizon, sparse(H), E, F, G, W, model, OSQP.Results(), zeros(nU), u_min, u_max)
end

function rpo_ref_preview(plan::RPOPlan, t_elapsed_s::Real, dt::Real, horizon::Int)
    nx = 6
    out = zeros(nx, horizon + 1)
    n_ref = size(plan.r_ref_rtn, 2)
    n_ref == 0 && return out
    start_idx = clamp(Int(floor(Float64(t_elapsed_s) / max(Float64(dt), 1.0e-9))) + 1, 1, n_ref)
    for j in 0:horizon
        idx = min(start_idx + j, n_ref)
        out[1:3, j + 1] .= plan.r_ref_rtn[:, idx]
        out[4:6, j + 1] .= plan.v_ref_rtn[:, idx]
    end
    return out
end

function rpo_lqmpc_control(ctrl::RpoLQMPCController, x, x_ref)
    nx = size(ctrl.Ad, 1)
    nu = size(ctrl.Bd, 2)
    ref = Matrix{Float64}(x_ref)
    size(ref, 1) == nx || throw(ArgumentError("RPO MPC reference has $(size(ref, 1)) rows, expected $nx."))
    if size(ref, 2) < ctrl.horizon + 1
        padded = repeat(ref[:, end], 1, ctrl.horizon + 1)
        padded[:, 1:size(ref, 2)] .= ref
        ref = padded
    end
    ref_stack = vec(ref[:, 2:(ctrl.horizon + 1)])
    q = ctrl.E * Vector{Float64}(x) - ctrl.F * ref_stack
    OSQP.update_q!(ctrl.qp_model, 2.0 .* q)
    OSQP.warm_start_x!(ctrl.qp_model, ctrl.U_prev)
    results = OSQP.solve!(ctrl.qp_model, ctrl.qp_results)
    if results.info.status == :Solved || results.info.status == :Solved_inaccurate
        U = results.x
        ctrl.U_prev[1:end-nu] .= U[nu+1:end]
        ctrl.U_prev[end-nu+1:end] .= U[end-nu+1:end]
        return SVector{3, Float64}(U[1], U[2], U[3])
    end
    return SVector{3, Float64}(0.0, 0.0, 0.0)
end
