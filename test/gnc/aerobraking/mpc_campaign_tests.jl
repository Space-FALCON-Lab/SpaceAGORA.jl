using Test
using LinearAlgebra
using SpaceAGORA

@testset "selected planet polyfit feeds the MPC density adapter" begin
    planet = SpaceAGORA.SimulationModel.Earth()
    model = PolynomialFitAtmosphereModel(planet)
    context = (args=(environment_model=(planet=planet, density_model=model),),)
    density = density_function_from_spaceagora(context)
    for altitude_m in (130.0e3, 160.0e3, 300.0e3)
        expected = getDensity(model, altitude_m, 0.0, 0.0, 0.0, false,
            context)[1]
        @test density(altitude_m, 0.0) == expected
    end
    position = [0.0, 0.0, planet.Rp_p + 130.0e3]
    spherical_altitude = norm(position) - planet.Rp_e
    expected_geodetic = getDensity(model, 130.0e3, 0.0, 0.0, 0.0,
        false, context)[1]
    @test isapprox(density(spherical_altitude, 0.0, position),
        expected_geodetic; rtol=1.0e-12)

    gradient_position = [planet.Rp_e + 160.0e3, 0.0, 0.0]
    analytical_gradient = density(Val(:gradient), gradient_position, 0.0)
    step_m = 1.0
    numerical_gradient = [
        (density(0.0, 0.0, gradient_position + step_m .* unit_vector) -
         density(0.0, 0.0, gradient_position - step_m .* unit_vector)) /
        (2.0 * step_m)
        for unit_vector in ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0])
    ]
    @test isapprox(analytical_gradient, numerical_gradient; rtol=1.0e-7,
        atol=1.0e-20)
end

@testset "shared orbital elements to Cartesian initialization" begin
    planet = SpaceAGORA.SimulationModel.Earth(
        Rp_e=6_378.13e3, Rp_p=6_378.13e3, Rp_m=6_378.13e3,
        μ=3.986e14, J2=1.08263e-3)
    apoapsis_radius = 56_378.0e3
    periapsis_radius = planet.Rp_e + 130.0e3
    for anomaly_deg in (0.0, 45.0, 180.0)
        initial_condition = SpaceAGORA.SimulationModel.InitialCondition(
            ra=apoapsis_radius,
            rp=periapsis_radius,
            i=89.876,
            Ω=104.115,
            ω=75.505,
            ν=anomaly_deg,
        )
        position, velocity = orbital_elements_to_cartesian(
            initial_condition, planet)
        engine_position, engine_velocity =
            SpaceAGORA.SimulationEngine.orbitalelemtorv(
                initial_condition, planet)
        @test position == engine_position
        @test velocity == engine_velocity
        @test isapprox(
            norm(cross(position, velocity))^2 / planet.μ,
            initial_condition.a * (1 - initial_condition.e^2);
            rtol=1.0e-12)
        @test isapprox(
            0.5 * dot(velocity, velocity) - planet.μ / norm(position),
            -planet.μ / (2 * initial_condition.a);
            rtol=1.0e-12)
        if anomaly_deg == 0.0
            @test isapprox(norm(position), periapsis_radius; rtol=1.0e-12)
        elseif anomaly_deg == 180.0
            @test isapprox(norm(position), apoapsis_radius; rtol=1.0e-12)
        end
    end
end

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


@testset "nine-state KS MPC linearization" begin
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

    A, B, C, D, output = hooks.aerobraking_mpc_step_matrices(
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
        return ks_implicit_midpoint_step(
            candidate_state,
            params,
            candidate_area,
            delta_s;
            density_kg_m3=density_value,
            drag_coefficient=config.drag_coefficient,
            mass_kg=config.mass_kg,
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
    @test norm(A - numerical_step_jacobian) / norm(numerical_step_jacobian) < 1.0e-5

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
        return hooks.aerobraking_mpc_step_matrices(
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

    params = AerobrakingMPCParams(
        Re=6.378e6, μ=3.986e14, J2=1.08263e-3, Ω=1.0e-4)
    position = [params.Re + 100.0e3, 0.0, 0.0]
    velocity = [0.0, 7.6e3, 0.0]
    area = 4.0
    density_value = 2.0e-9
    density = (_, _) -> density_value
    output = evaluate_cartesian_mpc_outputs(
        position, velocity, 12.0, area, params, cfg; density=density)
    relative_speed = velocity[2] - params.Ω * position[1]
    @test output.altitude_m ≈ 100.0e3
    @test output.density_kg_m3 == density_value
    @test output.relative_speed_m_s ≈ relative_speed
    @test output.drag_n ≈
        0.5 * density_value * relative_speed^2 * cfg.drag_coefficient * area
    @test output.heat_rate_w_cm2 ≈
        0.5 * density_value * relative_speed^3 *
        commanded_area_fraction(cfg, area) / 1.0e4

    ks_state = cartesian_to_ks_state(
        position, velocity, params; elapsed_time_s=12.0)
    ks_output = evaluate_ks_mpc_outputs(
        ks_state, area, params, cfg; density=density)
    @test ks_output.drag_n ≈ output.drag_n
    @test ks_output.heat_rate_w_cm2 ≈ output.heat_rate_w_cm2
    @test cumulative_mpc_heat_load([1.0, 3.0, 2.0], [0.0, 2.0, 5.0]) ≈
        [0.0, 4.0, 11.5]

end

@testset "combined MPC constraints retain heat-load feasibility" begin
    N = 3
    ny = 4
    H = zeros(N * ny, N)
    H[3, 1] = 1.0e3
    H[7, 2] = 1.0e3
    H[11, 3] = 1.0e3
    H[12, :] .= -1.0
    problem = AerobrakingMPCProblem(
        params=AerobrakingMPCParams(
            Re=6.378e6, μ=3.986e14, J2=0.0, Ω=0.0),
        H=H,
        Mx=zeros(N * ny, 9),
        δX0=zeros(9),
        N=N,
        ny=ny,
        t=[0.0, 100.0, 200.0],
        Ybar=zeros(N, ny),
        Xbar=zeros(N, 9),
        Abar_m2=2.0,
    )
    config = AerobrakingMPCConfig(
        mode=MaxEnergyDepletionMode(),
        bus_reference_area_m2=2.0,
        controllable_area_m2=4.0,
        mass_kg=500.0,
        drag_coefficient=2.2,
        qdot_max_w_cm2=1.0,
        heat_load_max_j_cm2=10.0,
        drag_max_n=100.0,
        area_slew_max_m2_s=0.2,
        use_constraints=true,
        use_slew_constraint=true,
        use_qdot_constraint=true,
        use_heat_load_constraint=true,
        use_drag_constraint=true,
        target_energy_mj_kg=-20.0,
        area_weight=1.0e-8,
        area_slew_weight=0.0,
        slack_weight=1.0e3,
        target_energy_weight=0.0,
        max_depletion_energy_weight=1.0,
        osqp_eps_abs=1.0e-7,
        osqp_eps_rel=1.0e-7,
        osqp_max_iter=10_000,
    )
    solution = solve_mpc_qp(problem, config)
    rectangle_heat_load_j_m2 = sum(
        solution.predicted_outputs[:, 3] .* [100.0, 100.0, 100.0])
    @test solution.ok
    @test rectangle_heat_load_j_m2 <=
        config.heat_load_max_j_cm2 * 1.0e4 + 1.0e-2
end

@testset "MPC slew constraint changes a rate-limited optimum" begin
    N = 3
    ny = 4
    H = zeros(N * ny, N)
    H[N * ny, N] = -1.0
    problem = AerobrakingMPCProblem(
        params=AerobrakingMPCParams(
            Re=6.378e6, μ=3.986e14, J2=0.0, Ω=0.0),
        H=H,
        Mx=zeros(N * ny, 9),
        δX0=zeros(9),
        N=N,
        ny=ny,
        t=[0.0, 1.0, 2.0],
        Ybar=zeros(N, ny),
        Xbar=zeros(N, 9),
        Abar_m2=2.0,
    )
    unconstrained = _test_mpc_config(
        MaxEnergyDepletionMode(); constraints=mpc_constraints(
            heat_rate=false, heat_load=false, drag=false, slew=false))
    rate_limited = _test_mpc_config(
        MaxEnergyDepletionMode(); constraints=mpc_constraints(:slew))
    unconstrained_solution = solve_mpc_qp(problem, unconstrained)
    rate_limited_solution = solve_mpc_qp(problem, rate_limited)
    rate_limited_slew = diff(rate_limited_solution.commanded_area_m2) ./
        diff(problem.t)

    @test unconstrained_solution.ok
    @test rate_limited_solution.ok
    @test maximum(abs.(rate_limited_slew)) <=
        rate_limited.area_slew_max_m2_s + 1.0e-6
    @test rate_limited_solution.commanded_area_m2[1] ≈ problem.Abar_m2 atol=1.0e-6
    @test maximum(abs.(unconstrained_solution.commanded_area_m2 .-
        rate_limited_solution.commanded_area_m2)) > 1.0
end
