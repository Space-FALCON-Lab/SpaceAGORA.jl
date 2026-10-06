# Aerobraking MPC prediction matrices and nominal outputs.
function aerobraking_mpc_step_matrices(
    state::AbstractVector,
    params::AerobrakingMPCParams,
    config::AerobrakingMPCConfig;
    Δs::Real=1.7e-7,
    area_m2::Real=config.bus_reference_area_m2 + config.controllable_area_m2,
    density::Function=(altitude_m, elapsed_time_s) -> 0.0,
)
    X0 = Float64.(collect(state))
    p = SVector{4, Float64}(X0[1:4])
    p_prime = SVector{4, Float64}(X0[5:8])
    h_ks = X0[9]
    elapsed_time_s = length(X0) >= 10 ? X0[10] : 0.0

    r_vec = ks_position(p)
    r_norm = norm(r_vec)
    v_vec = ks_velocity(p, p_prime)
    area = Float64(area_m2)

    Ωx = ks_rotation_cross_matrix(params)
    Vrel = v_vec - Ωx * r_vec
    Vn = sqrt(dot(Vrel, Vrel)) + eps(Float64)
    ρ, dρ_dr = ks_density_value_gradient(
        density, r_vec, elapsed_time_s, params; gradient_step_m=10.0)

    Jr_u, Jv_u, Jv_w, _, _, _ = ks_kinematics_jacobians(p, p_prime)

    step_s = Float64(Δs)
    tangent = ks_first_order_tangent_map(
        X0,
        params,
        area,
        step_s;
        density_kg_m3=density,
        drag_coefficient=config.drag_coefficient,
        mass_kg=config.mass_kg,
        use_drag=true,
    )
    A_k = tangent.transition[1:9, 1:9]
    B_k = tangent.input_transition[1:9, :]

    rhat = r_vec / (r_norm + eps(Float64))
    Ch_u = transpose(rhat) * Jr_u
    Ch_w = zeros(1, 4)

    Cq_r = (0.5 * Vn^2) * transpose(dρ_dr) + ρ * transpose(Vrel) * (-Ωx)
    Cq_v = ρ * transpose(Vrel)
    Cdrag_u = (config.drag_coefficient * area) .* (Cq_r * Jr_u + Cq_v * Jv_u)
    Cdrag_w = (config.drag_coefficient * area) .* (Cq_v * Jv_w)
    Ddrag = 0.5 * ρ * Vn^2 * config.drag_coefficient

    panel_area = max(config.controllable_area_m2, eps(Float64))
    Γ = 0.5 * ρ * Vn^3 / panel_area
    Vn2 = max(Vn * Vn, eps(Float64))
    Γ_over_ρ = 0.5 * Vn^3 / panel_area
    dΓ_dr = Γ_over_ρ * transpose(dρ_dr) + (3.0 * Γ / Vn2) * transpose(Vrel) * (-Ωx)
    dΓ_dv = (3.0 * Γ / Vn2) * transpose(Vrel)
    Cqd_u = (area - config.bus_reference_area_m2) * (dΓ_dr * Jr_u + dΓ_dv * Jv_u)
    Cqd_w = (area - config.bus_reference_area_m2) * (dΓ_dv * Jv_w)
    Dqd = Γ

    energy_unit_scale = config.mode isa TargetEnergyMode ? 1.0e6 : 1.0
    C_k = vcat(
        hcat(Ch_u, Ch_w, 0.0),
        hcat(Cdrag_u, Cdrag_w, 0.0),
        hcat(Cqd_u, Cqd_w, 0.0),
        hcat(zeros(1, 8), -1.0 / energy_unit_scale),
    )
    D_k = reshape([0.0, Ddrag, Dqd, 0.0], 4, 1)

    h_out = r_norm - params.Re
    drag_out = 0.5 * ρ * Vn^2 * config.drag_coefficient * area
    qdot_out = Γ * (area - config.bus_reference_area_m2)
    E_out = -h_ks / energy_unit_scale
    y_out = @SVector [h_out, drag_out, qdot_out, E_out]

    return Matrix(A_k), B_k, Matrix(C_k), D_k, y_out
end

function build_mpc_problem(
    reference,
    params::AerobrakingMPCParams,
    config::AerobrakingMPCConfig;
    density::Function=(altitude_m, elapsed_time_s) -> 0.0,
    max_nodes::Union{Nothing, Int}=nothing,
    zero_initial_deviation::Bool=true,
)
    states = reference.states[:, 1:9]
    times = collect(reference.time_s)
    idx_full = collect(1:size(states, 1))
    idx = if max_nodes !== nothing && length(idx_full) > max_nodes
        unique(round.(Int, LinRange(first(idx_full), last(idx_full), max_nodes)))
    else
        idx_full
    end
    N = length(idx)
    ny = 4
    Xbar = states[idx, 1:9]
    t = times[idx]
    X0 = vec(Xbar[1, :])
    δX0 = zero_initial_deviation ? zeros(length(X0)) : X0 .- vec(Xbar[min(2, N), :])

    A_hist = Vector{Matrix{Float64}}(undef, max(N - 1, 1))
    B_hist = Vector{Matrix{Float64}}(undef, N)
    C_hist = Vector{Matrix{Float64}}(undef, N)
    D_hist = Vector{Matrix{Float64}}(undef, N)
    Ybar = zeros(Float64, N, ny)
    area = Float64(reference.nominal_area_m2)

    for k in 1:N
        Δt_k = k < N ? (t[k + 1] - t[k]) : (N > 1 ? (t[k] - t[k - 1]) : 0.0)
        r_sund_k = max(sum(abs2, Xbar[k, 1:4]), eps(Float64))
        Δs_k = Δt_k / r_sund_k
        state_k = vcat(Xbar[k, :], t[k])
        Ak, Bk, Ck, Dk, yk = aerobraking_mpc_step_matrices(
            state_k,
            params,
            config;
            Δs=Δs_k,
            area_m2=area,
            density=density,
        )
        k < N && (A_hist[k] = Ak)
        B_hist[k] = Bk
        C_hist[k] = Ck
        D_hist[k] = Dk
        Ybar[k, :] .= yk
    end
    N == 1 && (A_hist[1] = Matrix(I, 9, 9))
    _, Mx, H = condensed_form(A_hist, B_hist, C_hist, D_hist)
    return AerobrakingMPCProblem(
        params=params,
        H=H,
        Mx=Mx,
        δX0=δX0,
        N=N,
        ny=ny,
        t=t,
        Ybar=Ybar,
        Xbar=Xbar,
        Abar_m2=area,
    )
end
