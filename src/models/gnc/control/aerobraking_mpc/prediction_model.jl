# Aerobraking MPC prediction matrices and nominal outputs.
#
# Outputs per node: y = [altitude (m), drag (N), panel heat rate (W/m²), specific
# energy]. Drag is ½ρV²C_D A and the heat rate is the kinetic-energy flux on the
# exposed panel, ½ρV³(A - A_bus)/A_panel, with V the speed relative to the
# co-rotating atmosphere. The energy is -h_KS, in MJ/kg for energy targeting and
# J/kg for maximum depletion, the units of the two benchmark formulations.

# Discrete maps of one node for deviations from the nominal pass, on the 9-state
# [u, u′, h_KS] (time is not a state of the prediction):
# δx_{k+1} ≈ A δx_k + B ΔA_k and y_k ≈ ȳ_k + C δx_k + D ΔA_k.
function aerobraking_mpc_step_matrices(
    state::AbstractVector,
    params::KSPropagationParams,
    config::AerobrakingMPCConfig;
    Δs::Real=1.7e-7,
    area_m2::Real=config.bus_reference_area_m2 + config.controllable_area_m2,
    density::Function=(altitude_m, elapsed_time_s) -> 0.0,
)
    length(state) == 10 || throw(ArgumentError("KS state must have 10 components."))
    area = Float64(area_m2)
    tangent = ks_first_order_tangent_map(
        state,
        params,
        area,
        Float64(Δs);
        density_kg_m3=density,
        drag_coefficient=config.drag_coefficient,
        mass_kg=config.mass_kg,
        use_drag=true,
    )
    A_k = tangent.transition[1:9, 1:9]
    B_k = tangent.input_transition[1:9, :]

    partials = tangent.evaluation.partials
    ρ, dρ_dr = partials.density, partials.density_gradient
    Jr_u, Jv_u, Jv_w = partials.position_u, partials.velocity_u, partials.velocity_u_prime
    r_vec, Vrel = partials.position, partials.relative_velocity
    r_norm = norm(r_vec)
    Vn = norm(Vrel) + eps(Float64)
    dV_dr = -ks_rotation_cross_matrix(params)    # V = v - Ω × r

    # Altitude row: ∂|r|/∂u.
    Ch_u = transpose(r_vec / r_norm) * Jr_u
    Ch_w = zeros(1, 4)

    # Drag row: ∂(½ρV²)/∂r and ∂(½ρV²)/∂v, times C_D A.
    Cq_r = (0.5 * Vn^2) * transpose(dρ_dr) + ρ * transpose(Vrel) * dV_dr
    Cq_v = ρ * transpose(Vrel)
    Cdrag_u = (config.drag_coefficient * area) .* (Cq_r * Jr_u + Cq_v * Jv_u)
    Cdrag_w = (config.drag_coefficient * area) .* (Cq_v * Jv_w)
    Ddrag = 0.5 * ρ * Vn^2 * config.drag_coefficient

    # Heat-rate row: Γ = ½ρV³/A_panel per m² of exposed panel area.
    panel_area = max(config.controllable_area_m2, eps(Float64))
    exposed_panel_area = area - config.bus_reference_area_m2
    Γ = 0.5 * ρ * Vn^3 / panel_area
    dΓ_dr = (0.5 * Vn^3 / panel_area) * transpose(dρ_dr) +
        (3.0 * Γ / Vn^2) * transpose(Vrel) * dV_dr
    dΓ_dv = (3.0 * Γ / Vn^2) * transpose(Vrel)
    Cqd_u = exposed_panel_area * (dΓ_dr * Jr_u + dΓ_dv * Jv_u)
    Cqd_w = exposed_panel_area * (dΓ_dv * Jv_w)

    # Energy row: E = -h_KS, so it reads the ninth state directly.
    energy_unit_scale = config.mode isa TargetEnergyMode ? 1.0e6 : 1.0
    C_k = vcat(
        hcat(Ch_u, Ch_w, 0.0),
        hcat(Cdrag_u, Cdrag_w, 0.0),
        hcat(Cqd_u, Cqd_w, 0.0),
        hcat(zeros(1, 8), -1.0 / energy_unit_scale),
    )
    D_k = reshape([0.0, Ddrag, Γ, 0.0], 4, 1)

    y_out = @SVector [
        r_norm - params.Re,
        Ddrag * area,
        Γ * exposed_panel_area,
        -Float64(state[9]) / energy_unit_scale,
    ]
    return Matrix(A_k), B_k, Matrix(C_k), D_k, y_out
end

"""
    build_mpc_problem(reference, params, config; density, max_nodes, zero_initial_deviation)

Linearise the KS dynamics along a [`build_reference_drag_pass`](@ref) result and
condense the prediction into an [`AerobrakingMPCProblem`](@ref). `max_nodes`
subsamples the reference nodes.
"""
function build_mpc_problem(
    reference,
    params::KSPropagationParams,
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
    # As in the benchmark scripts, the optional initial deviation is the first
    # node's offset from the second.
    δX0 = zero_initial_deviation ? zeros(length(X0)) : X0 .- vec(Xbar[min(2, N), :])

    A_hist = Vector{Matrix{Float64}}(undef, max(N - 1, 1))
    B_hist = Vector{Matrix{Float64}}(undef, N)
    C_hist = Vector{Matrix{Float64}}(undef, N)
    D_hist = Vector{Matrix{Float64}}(undef, N)
    Ybar = zeros(Float64, N, ny)
    area = Float64(reference.nominal_area_m2)

    for k in 1:N
        # Fictitious-time step to the next node, from dt = |r| ds at this node
        # (the last node reuses the previous interval).
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
