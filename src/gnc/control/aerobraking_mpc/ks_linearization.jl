#=
"""
    Linearized KS model for MPC prediction.

    The output rows are altitude, drag, heat rate, and specific energy. The
    spacecraft mass, area, planet constants, and density model come from the
    user-selected SpaceAGORA case.
"""
=#
@inline function _density_with_gradient(density::Function, altitude_m::Real, elapsed_time_s::Real; dh_m::Real=10.0)
    h = Float64(altitude_m)
    t = Float64(elapsed_time_s)
    step = max(abs(Float64(dh_m)), 1.0)
    ρ = max(0.0, Float64(density(h, t)))
    ρp = max(0.0, Float64(density(h + step, t)))
    ρm = max(0.0, Float64(density(h - step, t)))
    return ρ, (ρp - ρm) / (2.0 * step)
end

@inline function _dΛ_du(i::Int)
    if i == 1
        return @SMatrix [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0]
    elseif i == 2
        return @SMatrix [0.0 -1.0 0.0 0.0; 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0]
    elseif i == 3
        return @SMatrix [0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 -1.0; 1.0 0.0 0.0 0.0]
    end
    return @SMatrix [0.0 0.0 0.0 1.0; 0.0 0.0 -1.0 0.0; 0.0 1.0 0.0 0.0]
end

@inline function _ks_kinematics_jacobians(u, w)
    L3 = Λ(u)
    r = dot(u, u)
    v = (2.0 / r) * (L3 * w)
    Jr_u = 2.0 * L3
    Jv_u = (2.0 / r) * Λ(w) - (4.0 / (r * r)) * ((L3 * w) * transpose(u))
    Jv_w = (2.0 / r) * L3
    return Jr_u, Jv_u, Jv_w, r, v, L3
end

function _j2_acceleration_jacobian_si(rvec, params::AerobrakingMPCParams)
    x, y, z = rvec
    r2 = dot(rvec, rvec)
    r1 = sqrt(r2)
    r1 <= eps(Float64) && return @SMatrix zeros(3, 3)
    r4 = r2 * r2
    r5 = r1^5
    r7 = r1^7

    k = 1.5 * params.J2 * params.μ * params.Re^2
    z2_r2 = z * z / r2
    U = 5.0 * z2_r2 - 1.0
    W = 5.0 * z2_r2 - 3.0
    b = @SVector [x * U, y * U, z * W]

    dUdx = -10.0 * z * z * x / r4
    dUdy = -10.0 * z * z * y / r4
    dUdz = 10.0 * z * (x * x + y * y) / r4

    colx = @SVector [U + x * dUdx, y * dUdx, z * dUdx]
    coly = @SVector [x * dUdy, U + y * dUdy, z * dUdy]
    colz = @SVector [x * dUdz, y * dUdz, W + z * dUdz]
    db = SMatrix{3, 3, Float64}(hcat(colx, coly, colz))

    df = (-5.0 * k / r7) * rvec
    return (k / r5) * db + b * transpose(df)
end

function _drag_jacobians_si(
    rvec,
    vvec,
    area_m2::Real,
    params::AerobrakingMPCParams,
    config::AerobrakingMPCConfig,
    density::Function,
    elapsed_time_s::Real,
)
    altitude_m = norm(rvec) - params.Re
    ρ, dρ_dh = _density_with_gradient(density, altitude_m, elapsed_time_s)
    Ωx = _skew_rotation(params)
    V = vvec - Ωx * rvec
    Vn = sqrt(dot(V, V)) + eps(Float64)

    K = 0.5 * ρ * config.drag_coefficient * Float64(area_m2) / config.mass_kg
    I3 = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    Mv = Vn * I3 + (V * transpose(V)) / Vn
    Jv = -K * Mv

    r1 = sqrt(dot(rvec, rvec)) + eps(Float64)
    rhat = (1.0 / r1) * rvec
    dKdr = (0.5 * config.drag_coefficient * Float64(area_m2) / config.mass_kg) * dρ_dh * rhat
    c = Vn * V
    Jr = -(c * transpose(dKdr)) + K * Mv * Ωx
    JA = -(0.5 * ρ * config.drag_coefficient / config.mass_kg) * c

    return Jr, Jv, JA, ρ
end

function linearized_ks_dynamics(
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
    h = X0[9]
    elapsed_time_s = length(X0) >= 10 ? X0[10] : 0.0

    r_vec = Λ(p) * p
    r_norm = norm(r_vec)
    v_vec = (2.0 / r_norm) * Λ(p) * p_prime
    ω = sqrt(max(h, eps(Float64)) / 4.0)
    area = Float64(area_m2)

    Ωx = _skew_rotation(params)
    Vrel = v_vec - Ωx * r_vec
    Vn = sqrt(dot(Vrel, Vrel)) + eps(Float64)
    a_J2_vec = j2_acceleration_si(r_vec, params)
    JrD, JvD, JAD, ρ = _drag_jacobians_si(r_vec, v_vec, area, params, config, density, elapsed_time_s)
    a_D_vec = -0.5 * ρ * config.drag_coefficient * area / config.mass_kg * Vn * Vrel
    a_nom = a_J2_vec + a_D_vec

    Jr_u, Jv_u, Jv_w, r_sund, _, L3 = _ks_kinematics_jacobians(p, p_prime)
    J2r = _j2_acceleration_jacobian_si(r_vec, params)
    Ja_r = J2r + JrD
    Ja_v = JvD
    Ja_A = JAD

    G = 0.5 * r_sund * transpose(L3)
    Ga_u = zeros(Float64, 4, 4)
    for i in 1:4
        dG = 0.5 * (2.0 * p[i] * transpose(L3) + r_sund * transpose(_dΛ_du(i)))
        Ga_u[:, i] = dG * a_nom + G * (Ja_r * Jr_u[:, i] + Ja_v * Jv_u[:, i])
    end

    Ga_w = zeros(Float64, 4, 4)
    for i in 1:4
        Ga_w[:, i] = G * (Ja_v * Jv_w[:, i])
    end

    dGadx = zeros(Float64, 8, 8)
    dGadx[5:8, 1:4] .= Ga_u
    dGadx[5:8, 5:8] .= Ga_w
    A_k = Φ(ω, Float64(Δs)) + Float64(Δs) * dGadx

    B_lower = G * Ja_A
    B_k = reshape(Float64(Δs) .* vcat(zeros(4), B_lower), 8, 1)

    rhat = r_vec / (r_norm + eps(Float64))
    Ch_u = transpose(rhat) * Jr_u
    Ch_w = zeros(1, 4)

    _, dρ_dh = _density_with_gradient(density, r_norm - params.Re, elapsed_time_s)
    dρ_dr = dρ_dh * rhat
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
    CE_r = transpose((params.μ / (r_norm^3)) * r_vec) ./ energy_unit_scale
    CE_v = transpose(v_vec) ./ energy_unit_scale
    CE_u = CE_r * Jr_u + CE_v * Jv_u
    CE_w = CE_v * Jv_w

    C_k = vcat(
        hcat(Ch_u, Ch_w),
        hcat(Cdrag_u, Cdrag_w),
        hcat(Cqd_u, Cqd_w),
        hcat(CE_u, CE_w),
    )
    D_k = reshape([0.0, Ddrag, Dqd, 0.0], 4, 1)

    h_out = r_norm - params.Re
    drag_out = 0.5 * ρ * Vn^2 * config.drag_coefficient * area
    qdot_out = Γ * (area - config.bus_reference_area_m2)
    E_out = (0.5 * dot(v_vec, v_vec) - params.μ / r_norm) / energy_unit_scale
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
    Xbar = states[idx, 1:8]
    hbar = states[idx, 9]
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
        state_k = vcat(Xbar[k, :], hbar[k], t[k])
        Ak, Bk, Ck, Dk, yk = linearized_ks_dynamics(
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
    N == 1 && (A_hist[1] = Matrix(I, 8, 8))
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
