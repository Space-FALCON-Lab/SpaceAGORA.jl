using Test
using LinearAlgebra
using SpaceAGORA

function _test_mpc_config(mode; constraints=mpc_constraints(
        heat_rate=false, heat_load=false, drag=false, slew=false))
    return AerobrakingMPCConfig(
        mode=mode,
        bus_reference_area_m2=2.0,
        controllable_area_m2=4.0,
        mass_kg=500.0,
        drag_coefficient=2.2,
        qdot_max_w_cm2=1.0,
        heat_load_max_j_cm2=10.0,
        drag_max_n=100.0,
        area_slew_max_m2_s=0.2,
        use_constraints=!isempty(constraint_names(constraints)),
        use_slew_constraint=constraint_active(constraints, :slew),
        use_qdot_constraint=constraint_active(constraints, :heat_rate),
        use_heat_load_constraint=constraint_active(constraints, :heat_load),
        use_drag_constraint=constraint_active(constraints, :drag),
        target_energy_mj_kg=-20.0,
        area_weight=1.0e-5,
        area_slew_weight=0.0,
        slack_weight=1.0e3,
        target_energy_weight=1.0e-8,
        max_depletion_energy_weight=1.0,
        osqp_eps_abs=1.0e-7,
        osqp_eps_rel=1.0e-7,
        osqp_max_iter=10_000,
    )
end


@testset "paper-convention nine-state KS MPC linearization" begin
    hooks = SpaceAGORA.SimulationModel.ControlHooks
    params = AerobrakingMPCParams(
        Re=6_378_137.0,
        μ=3.986004418e14,
        J2=1.08262668e-3,
        Ω=7.292115e-5,
    )
    config = _test_mpc_config(TargetEnergyMode())
    position = [params.Re + 140.0e3, 0.0, 50.0e3]
    velocity = [-250.0, 7.75e3, 300.0]
    state = cartesian_to_ks_state(position, velocity, params; elapsed_time_s=30.0)
    density_value = 2.0e-9
    density = (_, _) -> density_value
    area_m2 = 5.0
    delta_s = 1.0e-9

    A, B, C, D, output = hooks.linearized_ks_dynamics(
        state,
        params,
        config;
        Δs=delta_s,
        area_m2=area_m2,
        density=density,
    )
    @test size(A) == (9, 9)
    @test size(B) == (9, 1)
    @test size(C) == (4, 9)
    @test size(D) == (4, 1)
    @test length(output) == 4
    first_order_h_column = -0.5 .* delta_s .* state[1:4]
    @test norm(A[5:8, 9] - first_order_h_column) /
        norm(first_order_h_column) < 1.0e-5
    @test B[9, 1] > 0.0
    @test C[4, 1:8] == zeros(8)
    @test C[4, 9] == -1.0e-6
    @test output[4] == -state[9] / 1.0e6

    function step_at(candidate_state, candidate_area)
        return ks_rk4_step(
            candidate_state,
            params,
            candidate_area,
            delta_s;
            config=config,
            density_kg_m3=density_value,
            use_drag=true,
        )[1:9]
    end
    numerical_step_jacobian = zeros(9, 9)
    for column in 1:9
        state_step = max(abs(state[column]), 1.0) * 2.0e-4
        state_plus = copy(state)
        state_minus = copy(state)
        state_plus[column] += state_step
        state_minus[column] -= state_step
        coarse = (
            step_at(state_plus, area_m2) - step_at(state_minus, area_m2)
        ) ./ (2.0 * state_step)
        half_step = 0.5 * state_step
        state_plus[column] = state[column] + half_step
        state_minus[column] = state[column] - half_step
        fine = (
            step_at(state_plus, area_m2) - step_at(state_minus, area_m2)
        ) ./ (2.0 * half_step)
        numerical_step_jacobian[:, column] .= (4.0 .* fine .- coarse) ./ 3.0
    end
    @test norm(A - numerical_step_jacobian) / norm(numerical_step_jacobian) < 5.0e-6

    area_step = 2.0e-3
    numerical_B = (
        -step_at(state, area_m2 + 2.0 * area_step) +
        8.0 * step_at(state, area_m2 + area_step) -
        8.0 * step_at(state, area_m2 - area_step) +
        step_at(state, area_m2 - 2.0 * area_step)
    ) ./ (12.0 * area_step)
    @test norm(vec(B) - numerical_B) / max(norm(numerical_B), eps(Float64)) < 1.0e-3

    state_perturbation = 1.0e-6 .* max.(abs.(state[1:9]), 1.0) .*
        [1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0]
    area_perturbation = 1.0e-3
    nominal_next = step_at(state, area_m2)
    nonlinear_next = step_at(
        vcat(state[1:9] .+ state_perturbation, state[10]),
        area_m2 + area_perturbation,
    )
    linearized_next = nominal_next + A * state_perturbation + vec(B) * area_perturbation
    @test norm(nonlinear_next - linearized_next) /
        max(norm(nonlinear_next - nominal_next), eps(Float64)) < 2.0e-5

    function output_at(candidate_state, candidate_area)
        return hooks.linearized_ks_dynamics(
            candidate_state,
            params,
            config;
            Δs=delta_s,
            area_m2=candidate_area,
            density=density,
        )[5]
    end
    numerical_C = zeros(4, 9)
    for column in 1:9
        state_step = max(abs(state[column]), 1.0) * 1.0e-7
        state_plus = copy(state)
        state_minus = copy(state)
        state_plus[column] += state_step
        state_minus[column] -= state_step
        numerical_C[:, column] .= (
            output_at(state_plus, area_m2) - output_at(state_minus, area_m2)
        ) ./ (2.0 * state_step)
    end
    numerical_D = (
        output_at(state, area_m2 + area_step) -
        output_at(state, area_m2 - area_step)
    ) ./ (2.0 * area_step)
    for row in 1:4
        @test norm(C[row, :] - numerical_C[row, :]) /
            max(norm(numerical_C[row, :]), 1.0) < 2.0e-6
    end
    @test norm(vec(D) - numerical_D) / max(norm(numerical_D), 1.0) < 2.0e-6
end

@testset "aerobraking MPC campaign foundations" begin
    cfg = _test_mpc_config(TargetEnergyMode())
    for area in range(cfg.bus_reference_area_m2,
            cfg.bus_reference_area_m2 + cfg.controllable_area_m2; length=9)
        alpha = alpha_from_commanded_area(cfg, area; min_alpha_rad=0.0, max_alpha_rad=pi / 2)
        @test commanded_area_from_alpha(cfg, alpha; min_alpha_rad=0.0, max_alpha_rad=pi / 2) ≈ area
    end

end
