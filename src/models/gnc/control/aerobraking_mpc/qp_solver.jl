# Condensed OSQP formulation for energy targeting and maximum depletion.
const _WCM2_TO_WM2 = 1.0e4
const _JCM2_TO_JM2 = 1.0e4

function _append_constraint_rows!(
    rows::Vector{Vector{Float64}},
    lower::Vector{Float64},
    upper::Vector{Float64},
    A,
    l,
    u,
)
    for i in axes(A, 1)
        push!(rows, vec(Float64.(A[i, :])))
        push!(lower, Float64(l[i]))
        push!(upper, Float64(u[i]))
    end
    return nothing
end

function _osqp_status_symbol(status)
    return Symbol(replace(lowercase(String(status)), " " => "_"))
end

function _slack_offsets(config::AerobrakingMPCConfig, N::Int)
    ns_qdot = (config.use_constraints && config.use_qdot_constraint) ? N : 0
    ns_Q = (config.use_constraints && config.use_heat_load_constraint) ? 1 : 0
    ns_drag = (config.use_constraints && config.use_drag_constraint) ? N : 0
    qdot_range = (N + 1):(N + ns_qdot)
    heat_load_idx = N + ns_qdot + 1
    drag_range = (N + ns_qdot + ns_Q + 1):(N + ns_qdot + ns_Q + ns_drag)
    return ns_qdot, ns_Q, ns_drag, qdot_range, heat_load_idx, drag_range
end

"""
    solve_mpc_qp(problem, config) -> AerobrakingMPCSolution

Solve the condensed MPC quadratic program with OSQP for the objective of
`config.mode` under the area bounds and the active limits of `config`. The
heat-rate rows also bound the predicted heat rate below by zero. Slacks
(weighted by `slack_weight`) soften the heat-rate upper bound, heat load and
drag in [`TargetEnergyMode`](@ref), and only the drag limit in
[`MaxEnergyDepletionMode`](@ref). `solution.slacks` always holds `N` heat-rate,
one heat-load and `N` drag entries for the active limits, in that order.
"""
function solve_mpc_qp(problem::AerobrakingMPCProblem, config::AerobrakingMPCConfig)
    N = problem.N
    ny = problem.ny
    nU = N
    H = problem.H
    Mx = problem.Mx
    δX0 = problem.δX0
    Ybar = problem.Ybar
    t = problem.t
    all(isfinite, H) || throw(ArgumentError("MPC prediction matrix H contains non-finite values."))
    all(isfinite, Mx) || throw(ArgumentError("MPC initial-state map Mx contains non-finite values."))
    all(isfinite, Ybar) || throw(ArgumentError("MPC nominal outputs contain non-finite values."))
    all(isfinite, t) || throw(ArgumentError("MPC node times contain non-finite values."))

    H_Erow, E_off = terminal_energy_affine_terms(problem)

    ns_qdot, ns_Q, ns_drag, qdot_range, heat_load_idx, drag_range = _slack_offsets(config, N)
    ns = ns_qdot + ns_Q + ns_drag
    nz = nU + ns

    Ddu = finite_difference_matrix(N)
    area_min = config.bus_reference_area_m2
    area_max = config.bus_reference_area_m2 + config.controllable_area_m2
    area_span = max(area_max - area_min, eps(Float64))
    Abar = clamp(problem.Abar_m2, area_min, area_max)
    Abar_vec = fill(Abar, N)

    # Cost ½ zᵀQz + qᵀz: area and area-change penalties, the objective of the
    # mode on the predicted end-of-pass energy E_N = H_Erow ⋅ ΔU + E_off, and the
    # slack penalty.
    Q_uu = 2.0 * config.area_weight * Matrix(I, nU, nU) +
        2.0 * config.area_slew_weight * (Ddu' * Ddu)
    q_u_vec = zeros(nU)
    if config.mode isa TargetEnergyMode
        # target_energy_weight * (1e6 (E_N - E*))², with E in MJ/kg.
        energy_cost_scale = 1.0e6
        hE = energy_cost_scale .* H_Erow
        cE = energy_cost_scale * (E_off - config.target_energy_mj_kg)
        Q_uu .+= 2.0 * config.target_energy_weight * (hE * hE')
        q_u_vec .+= 2.0 * config.target_energy_weight * cE .* hE
    else
        # Minimise E_N (J/kg), scaled by the energy and area magnitudes.
        energy_scale = max(1.0, maximum(abs.(Ybar[:, 4])))
        q_u_vec .+= (config.max_depletion_energy_weight / (energy_scale * area_span)) .* H_Erow
    end

    Q = zeros(nz, nz)
    Q[1:nU, 1:nU] .= Q_uu
    Q[nU+1:nz, nU+1:nz] .= 2.0 * config.slack_weight * Matrix(I, ns, ns)
    # A small diagonal keeps Q positive definite when the area weights are zero.
    Q = Matrix(Symmetric(Q + 1.0e-10 * Matrix(I, nz, nz)))
    q = zeros(nz)
    q[1:nU] .= q_u_vec

    rows = Vector{Vector{Float64}}()
    lower = Float64[]
    upper = Float64[]

    # Decision vector z = [ΔU (N); heat-rate (N), heat-load (1) and drag (N)
    # slacks of the active limits]. The slacks are nonnegative.
    if ns > 0
        A = zeros(ns, nz)
        A[:, nU+1:nz] .= Matrix(I, ns, ns)
        _append_constraint_rows!(rows, lower, upper, A, zeros(ns), fill(Inf, ns))
    end

    # Area bounds: A_bus <= Abar + ΔU_k <= A_bus + A_panel.
    A = zeros(N, nz)
    A[:, 1:nU] .= Matrix(I, N, N)
    _append_constraint_rows!(rows, lower, upper, A, fill(area_min - Abar, N), fill(area_max - Abar, N))

    Δt_nodes = length(t) > 1 ? vcat(diff(t), t[end] - t[end - 1]) : fill(0.0, N)
    # Slew: the first increment is zero and |A_{k+1} - A_k| <= Ȧ_max Δt_k.
    if config.use_constraints && config.use_slew_constraint && N > 1
        A = zeros(1, nz)
        A[1, 1] = 1.0
        _append_constraint_rows!(rows, lower, upper, A, [0.0], [0.0])

        Δt_du = Δt_nodes[1:end-1]
        Dbar = Ddu * Abar_vec
        A = zeros(N - 1, nz)
        A[:, 1:nU] .= Ddu
        _append_constraint_rows!(
            rows,
            lower,
            upper,
            A,
            -config.area_slew_max_m2_s .* Δt_du .- Dbar,
            config.area_slew_max_m2_s .* Δt_du .- Dbar,
        )
    end

    # Slack coupling follows the two benchmark formulations: energy targeting
    # softens the heat-rate upper bound, heat load and drag; maximum depletion
    # keeps heat rate and heat load hard and softens drag only. The unused slack
    # variables stay in the decision vector and are driven to zero by their cost.
    soft_path_limits = config.mode isa TargetEnergyMode
    yoff = Mx * δX0 .+ vec(Ybar')
    S_qdot = build_output_selector(N, ny; yidx=3)
    # Heat rate, normalised by the limit: q̇_k <= q̇_max (+ slack) and q̇_k >= 0.
    if config.use_constraints && config.use_qdot_constraint
        Hq = S_qdot * H
        bq = S_qdot * (Mx * δX0)
        qdot_bar = Ybar[:, 3]
        qdot_limit = fill(config.qdot_max_w_cm2 * _WCM2_TO_WM2, N)

        A = zeros(N, nz)
        A[:, 1:nU] .= Hq ./ qdot_limit
        soft_path_limits && (A[:, qdot_range] .= -Matrix(I, N, N))
        _append_constraint_rows!(
            rows,
            lower,
            upper,
            A,
            fill(-Inf, N),
            (qdot_limit .- qdot_bar .- bq) ./ qdot_limit,
        )

        A = zeros(N, nz)
        A[:, 1:nU] .= Hq ./ qdot_limit
        _append_constraint_rows!(
            rows,
            lower,
            upper,
            A,
            (-qdot_bar .- bq) ./ qdot_limit,
            fill(Inf, N),
        )
    end

    # Heat load: Σ_k q̇_k Δt_k <= Q_max (+ slack), a rectangle rule over the nodes.
    if config.use_constraints && config.use_heat_load_constraint
        Wrow = Δt_nodes' * S_qdot
        H_Q = vec(Wrow * H)
        heat_load_limit = config.heat_load_max_j_cm2 * _JCM2_TO_JM2
        rhs_Q = heat_load_limit - (Wrow * yoff)[1]
        A = zeros(1, nz)
        # Normalize the accumulated-load row so OSQP's absolute feasibility
        # tolerance remains meaningful when heat-load and slew limits coexist.
        A[1, 1:nU] .= H_Q ./ heat_load_limit
        soft_path_limits && (A[1, heat_load_idx] = -1.0)
        _append_constraint_rows!(
            rows, lower, upper, A, [-Inf], [rhs_Q / heat_load_limit])
    end

    # Drag, normalised by the limit: D_k <= D_max + slack.
    if config.use_constraints && config.use_drag_constraint
        S_drag = build_output_selector(N, ny; yidx=2)
        H_drag = S_drag * H
        b_drag = S_drag * (Mx * δX0)
        drag_bar = Ybar[:, 2]
        drag_limit = fill(config.drag_max_n, N)
        A = zeros(N, nz)
        A[:, 1:nU] .= H_drag ./ drag_limit
        A[:, drag_range] .= -Matrix(I, N, N)
        _append_constraint_rows!(
            rows,
            lower,
            upper,
            A,
            fill(-Inf, N),
            (drag_limit .- drag_bar .- b_drag) ./ drag_limit,
        )
    end

    A_mat = isempty(rows) ? zeros(0, nz) : reduce(vcat, transpose.(rows))
    model = OSQP.Model()
    OSQP.setup!(
        model;
        P=triu(sparse(Q)),
        q=q,
        A=sparse(A_mat),
        l=lower,
        u=upper,
        verbose=false,
        polish=false,
        warm_start=true,
        eps_abs=config.osqp_eps_abs,
        eps_rel=config.osqp_eps_rel,
        max_iter=config.osqp_max_iter,
    )
    results = OSQP.solve!(model)
    status = _osqp_status_symbol(results.info.status)
    ok = status in (:solved, :solved_inaccurate)
    if !ok
        return solution_from_area_delta(
            config.mode,
            problem,
            zeros(nU);
            ok=false,
            objective=NaN,
            slacks=zeros(ns),
            solver_status=status,
        )
    end
    z = Vector{Float64}(results.x)
    if !all(isfinite, z)
        return solution_from_area_delta(
            config.mode,
            problem,
            zeros(nU);
            ok=false,
            objective=NaN,
            slacks=zeros(ns),
            solver_status=:nonfinite_solution,
        )
    end
    return solution_from_area_delta(
        config.mode,
        problem,
        z[1:nU];
        ok=true,
        objective=Float64(results.info.obj_val),
        slacks=ns > 0 ? z[nU+1:end] : Float64[],
        solver_status=status,
    )
end
