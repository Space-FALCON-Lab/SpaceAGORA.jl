using Test
using LinearAlgebra
using SpaceAGORA
using SpaceAGORA.SimulationModel.DynamicsKS
using SpaceAGORA.SimulationModel.ControlHooks: commanded_area_fraction, interpolate_mpc_plan
using StaticArrays

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

@testset "planet rotation sign is preserved" begin
    venus = SpaceAGORA.SimulationModel.Venus()
    context = (environment_model=(planet=venus,),)
    @test venus.ω[3] < 0.0
    @test mpc_params_from_spaceagora(context).Ω == venus.ω[3]
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
    params = KSPropagationParams(
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

    function first_order_flow(candidate_state, candidate_area)
        derivative = nonlinear_ks_rhs(
            candidate_state, params, candidate_area;
            density_kg_m3=density_value,
            drag_coefficient=config.drag_coefficient,
            mass_kg=config.mass_kg,
            use_drag=true,
        )
        return candidate_state[1:9] + delta_s .* derivative[1:9]
    end
    numerical_step_jacobian = zeros(9, 9)
    for column in 1:9
        state_step = max(abs(state[column]), 1.0) * 2.0e-4
        state_plus = copy(state)
        state_minus = copy(state)
        state_plus[column] += state_step
        state_minus[column] -= state_step
        coarse = (
            first_order_flow(state_plus, area_m2) -
            first_order_flow(state_minus, area_m2)
        ) ./ (2.0 * state_step)
        half_step = 0.5 * state_step
        state_plus[column] = state[column] + half_step
        state_minus[column] = state[column] - half_step
        fine = (
            first_order_flow(state_plus, area_m2) -
            first_order_flow(state_minus, area_m2)
        ) ./ (2.0 * half_step)
        numerical_step_jacobian[:, column] .= (4.0 .* fine .- coarse) ./ 3.0
    end
    @test norm(A - numerical_step_jacobian) / norm(numerical_step_jacobian) < 1.0e-5

    area_step = 2.0e-3
    numerical_B = (
        -first_order_flow(state, area_m2 + 2.0 * area_step) +
        8.0 * first_order_flow(state, area_m2 + area_step) -
        8.0 * first_order_flow(state, area_m2 - area_step) +
        first_order_flow(state, area_m2 - 2.0 * area_step)
    ) ./ (12.0 * area_step)
    @test norm(vec(B) - numerical_B) / max(norm(numerical_B), eps(Float64)) < 1.0e-3

    state_perturbation = 1.0e-6 .* max.(abs.(state[1:9]), 1.0) .*
        [1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0]
    area_perturbation = 1.0e-3
    nominal_next = first_order_flow(state, area_m2)
    nonlinear_next = first_order_flow(
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

    params = KSPropagationParams(
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
        params=KSPropagationParams(
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
        params=KSPropagationParams(
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

# A pass whose nominal heat rate (2 W/cm²) and drag (200 N) exceed their limits
# with no control authority (H = 0): only slacks can make it feasible.
function _mpc_slack_probe(mode)
    N = 2
    ny = 4
    Ybar = zeros(N, ny)
    Ybar[:, 2] .= 200.0
    Ybar[:, 3] .= 2.0e4
    problem = AerobrakingMPCProblem(
        params=KSPropagationParams(Re=6.378e6, μ=3.986e14, J2=0.0, Ω=0.0),
        H=zeros(N * ny, N), Mx=zeros(N * ny, 9), δX0=zeros(9), N=N, ny=ny,
        t=[0.0, 100.0], Ybar=Ybar, Xbar=zeros(N, 9), Abar_m2=2.0)
    config = AerobrakingMPCConfig(
        mode=mode, bus_reference_area_m2=2.0, controllable_area_m2=4.0,
        mass_kg=500.0, drag_coefficient=2.2, qdot_max_w_cm2=1.0,
        heat_load_max_j_cm2=10.0, drag_max_n=100.0, area_slew_max_m2_s=0.2,
        use_constraints=true, use_slew_constraint=false, use_qdot_constraint=true,
        use_heat_load_constraint=true, use_drag_constraint=true,
        target_energy_mj_kg=-20.0, area_weight=1.0e-5, area_slew_weight=0.0,
        slack_weight=1.0e3, target_energy_weight=0.0,
        max_depletion_energy_weight=0.0, osqp_eps_abs=1.0e-7, osqp_eps_rel=1.0e-7,
        osqp_max_iter=10_000)
    return solve_mpc_qp(problem, config), N
end

@testset "MPC slack coupling follows the benchmark formulations" begin
    # Energy targeting: heat-rate, heat-load and drag slacks restore feasibility.
    solution, N = _mpc_slack_probe(TargetEnergyMode())
    @test solution.ok
    @test length(solution.slacks) == 2N + 1
    @test all(solution.slacks .>= -1.0e-8)
    @test maximum(solution.slacks[1:N]) > 0.5           # heat rate
    @test solution.slacks[N + 1] > 0.5                  # heat load
    @test maximum(solution.slacks[N + 2:end]) > 0.5     # drag

    # Maximum depletion keeps heat rate and heat load hard: the same pass is infeasible.
    solution, _ = _mpc_slack_probe(MaxEnergyDepletionMode())
    @test !solution.ok
    @test solution.solver_status in (:primal_infeasible, :primal_infeasible_inaccurate)
    @test all(iszero, solution.ΔU)
end

@testset "KS density gradient is one-sided at a density cutoff" begin
    params = KSPropagationParams(Re=6.378137e6, μ=3.986004418e14)
    top_m = 300.0e3
    density = (altitude_m, t) -> altitude_m < top_m ? 1.0e-11 * exp(-(altitude_m - top_m) / 5.0e4) : 0.0
    position = SVector(params.Re + top_m - 1.0, 0.0, 0.0)
    _, gradient = ks_density_value_gradient(density, position, 0.0, params; gradient_step_m=10.0)
    expected = -1.0e-11 * exp(1.0 / 5.0e4) / 5.0e4
    @test gradient[1] ≈ expected rtol=1.0e-3
    @test abs(gradient[2]) < 1.0e-20 && abs(gradient[3]) < 1.0e-20
end

# A density model that reports the time it was sampled at, to check the epoch.
struct _MPCTimeProbeDensity <: SpaceAGORA.AbstractDensityModel end
SpaceAGORA.getDensity(::_MPCTimeProbeDensity, h, lat, lon, t, wind, p) = (Float64(t), 200.0, zeros(3))

@testset "MPC density adapter: run-time epoch and planet-fixed sampling" begin
    planet = SpaceAGORA.SimulationModel.Earth()
    args = (environment_model=(planet=planet, density_model=_MPCTimeProbeDensity()),)
    hooks = SpaceAGORA.SimulationModel.ControlHooks
    source = density_function_from_spaceagora(args; epoch_s=1234.0)
    @test source(150.0e3, 6.0) == 1240.0
    @test source(150.0e3, 6.0, SVector(planet.Rp_e + 150.0e3, 0.0, 0.0)) == 1240.0

    # After a quarter turn of the planet, a fixed prediction-frame point on +x
    # lies at longitude -90° in the planet-fixed frame.
    Ω = source.rotation_rate_rad_s
    @test Ω ≈ 7.2921159e-5 rtol=1.0e-6
    τ = (pi / 2) / Ω
    position_pf = hooks._mpc_prediction_to_planet_fixed(source, SVector(7.0e6, 0.0, 0.0), τ)
    @test position_pf ≈ SVector(0.0, -7.0e6, 0.0) atol=1.0e-6
    back = hooks._mpc_planet_fixed_to_prediction(source, position_pf, τ)
    @test back ≈ SVector(7.0e6, 0.0, 0.0) atol=1.0e-6

    # The analytic polynomial-fit gradient points along the vertical in the
    # prediction frame, also after the planet has turned.
    polyfit = (environment_model=(planet=planet, density_model=PolynomialFitAtmosphereModel(planet)),)
    adapter = density_function_from_spaceagora(polyfit; epoch_s=0.0)
    position = SVector(planet.Rp_e + 150.0e3, 0.0, 0.0)
    gradient = adapter(Val(:gradient), position, τ)
    @test gradient[1] < 0.0
    @test norm(gradient[2:3]) <= 1.0e-9 * abs(gradient[1])
end

@testset "MPC control effectors bind to a spacecraft" begin
    lifecycle = SpaceAGORA.SimulationLifecycle
    runtime = SpaceAGORA.GNCRuntime
    config = AerobrakingMPCConfig(
        mode=MaxEnergyDepletionMode(), bus_reference_area_m2=5.74,
        controllable_area_m2=5.7, mass_kg=650.0, drag_coefficient=2.2,
        qdot_max_w_cm2=0.29, heat_load_max_j_cm2=40.0, drag_max_n=7.55,
        area_slew_max_m2_s=0.2, use_constraints=true, use_slew_constraint=false,
        use_qdot_constraint=true, use_heat_load_constraint=false,
        use_drag_constraint=true, target_energy_mj_kg=-6.41, area_weight=1.0e-5,
        area_slew_weight=0.0, slack_weight=1.0e-2, target_energy_weight=0.0,
        max_depletion_energy_weight=1.0, osqp_eps_abs=1.0e-5, osqp_eps_rel=1.0e-5,
        osqp_max_iter=10_000)
    controller = AerobrakingMPCControlModel(
        config=config, spacecraft_index=1, controlled_panel_links=(2, 3),
        min_alpha_rad=1.0e-4, max_alpha_rad=pi / 2, solve_interval_s=10.0,
        build_reference_on_tick=true, qp_max_nodes=40,
        reference=AerobrakingMPCReferenceConfig(h_cut_m=300.0e3, max_coast_steps=10, max_pass_steps=10),
        fallback_area_m2=11.44, prediction_latitude_rad=0.0,
        prediction_longitude_rad=0.0, prediction_wind=false)
    @test runtime.spacecraft_scope(controller) isa runtime.BoundSpacecraft
    bound = bind_spacecraft(controller, 3)
    @test bound.spacecraft_index == 3 && controller.spacecraft_index == 1
    @test bound.state !== controller.state
    @test runtime.held_command(bound, 3) == (commanded_area_m2=NaN, alpha_rad=NaN) ||
        all(isnan, values(runtime.held_command(bound, 3)))
    @test runtime.held_command(bound, 1) === nothing

    # Each atmosphere entry starts a new heat-load budget.
    controller.state.estimated_heat_load_j_cm2 = 12.0
    lifecycle.atmosphere_transition!(controller, nothing, nothing, 50.0, 1, false)
    @test controller.state.estimated_heat_load_j_cm2 == 12.0
    lifecycle.atmosphere_transition!(controller, nothing, nothing, 60.0, 1, true)
    @test controller.state.estimated_heat_load_j_cm2 == 0.0
    @test lifecycle.requires_atmosphere_events(controller)
end

isdefined(@__MODULE__, :ConformanceInstances) ||
    include(joinpath(@__DIR__, "..", "..", "fixtures", "conformance_instances.jl"))

@testset "MPC save fields read the run's controllers; campaign passes and preflight" begin
    template, campaign = ConformanceInstances.aerobraking_mpc_controls()
    # A run works on its own copy of the controller (here bound to spacecraft 2).
    running = bind_spacecraft(template, 2)
    running.state.held_commanded_area_m2 = 7.5
    running.state.estimated_heat_load_j_cm2 = 3.0
    integrator = (p=(args=(control_model=(control_effectors=(running,),),
        dynamics_model=(spacecraft=[nothing, nothing],)),),)
    fields = Dict(field.name => field for field in mpc_control_save_fields(template))
    @test isequal(fields[:mpc_commanded_area_m2].getter(nothing, 0.0, integrator), [NaN, 7.5])
    @test isequal(fields[:mpc_heat_load_j_cm2].getter(nothing, 0.0, integrator), [NaN, 3.0])
    @test all(isnan, fields[:mpc_drag_n].getter(nothing, 0.0, integrator))

    # Passes are counted at the run's atmosphere exits; entries restart heat accounting.
    lifecycle = SpaceAGORA.SimulationLifecycle
    lifecycle.atmosphere_transition!(campaign, nothing, nothing, 10.0, 1, true)
    @test campaign.state.completed_passes == 0
    lifecycle.atmosphere_transition!(campaign, nothing, nothing, 400.0, 1, false)
    lifecycle.atmosphere_transition!(campaign, nothing, nothing, 500.0, 2, false)
    @test campaign.state.completed_passes == 1
    campaign_fields = Dict(field.name => field for field in mpc_campaign_save_fields(campaign))
    run = (p=(args=(control_model=(control_effectors=(campaign,),), dynamics_model=(spacecraft=[nothing],)),),)
    @test campaign_fields[:mpc_campaign_completed_passes].getter(nothing, 0.0, run) == [1.0]
    @test campaign_fields[:mpc_campaign_phase_code].getter(nothing, 0.0, run) == [1.0]

    # The pair is checked before the run.
    swapped = AerobrakingMPCCampaignControlModel(
        maximum_depletion_control=campaign.targeting_control,
        targeting_control=campaign.maximum_depletion_control)
    @test_throws ArgumentError lifecycle.preflight_control(swapped, nothing; isolate_state=true)
    @test lifecycle.preflight_control(campaign, nothing; isolate_state=true) === nothing
end
