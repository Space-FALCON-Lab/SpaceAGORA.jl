@testset "Rigid-Body Angular Dynamics Uses Inertia Tensor" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.02, -0.03, 0.015)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    inertia_tensor = SMatrix{3, 3, Float64}(2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0)
    sc.inertia_tensor = inertia_tensor
    applied_torque = SVector{3, Float64}(0.12, -0.08, 0.05)

    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(ConstantTorqueModel(applied_torque),),
        keplerian=true
    )

    u0 = build_initial_conditions(args)
    du0 = copy(u0)
    du0 .= 0.0
    p = ODEParams(n_sats=1, args=args)
    spacecraft_dynamics!(du0, u0, p, 0.0)

    ω = SVector{3, Float64}(u0.sc[1].ω)
    qdot_expected = SimulationModel.DynamicsRotational.quaternion_derivative(ω, SVector{4, Float64}(u0.sc[1].q))
    ωdot_expected = SimulationModel.DynamicsRotational.angular_acceleration(
        ω,
        inertia_tensor,
        applied_torque;
        include_gyroscopic=true
    )
    ωdot_legacy = inertia_tensor \ (applied_torque - cross(ω, inertia_tensor * ω))
    @test isapprox(SVector{4, Float64}(du0.sc[1].q), qdot_expected; atol=1e-12, rtol=1e-10)
    @test isapprox(ωdot_expected, ωdot_legacy; atol=1e-12, rtol=1e-10)
    @test isapprox(SVector{3, Float64}(du0.sc[1].ω), ωdot_expected; atol=1e-12, rtol=1e-10)
end

@testset "Gravity-Gradient Torque Couples Into Rotational RHS" begin
    θ = deg2rad(45.0)
    q0 = normalize(SVector{4, Float64}(0.0, 0.0, sin(θ / 2), cos(θ / 2)))
    ω0 = SVector{3, Float64}(0.0, 0.0, 0.0)
    pos_ii = SVector{3, Float64}(EARTH.Rp_e + 500e3, 0.0, 0.0)
    inertia_tensor = SMatrix{3, 3, Float64}(2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0)

    r_hat_body = normalize(SimulationModel.rot(q0) * pos_ii)
    τ_expected = 3.0 * EARTH.μ / norm(pos_ii)^3 * cross(r_hat_body, inertia_tensor * r_hat_body)
    ωdot_expected = inertia_tensor \ τ_expected
    @test norm(τ_expected) > 0.0

    function gravity_gradient_ωdot(effector)
        sc = make_spacecraft(
            ra_alt_m=500e3,
            rp_alt_m=500e3,
            i_deg=35.0,
            ω_deg=40.0,
            Ω_deg=10.0,
            ν_deg=175.0,
            orientation_state=(q0, ω0)
        )
        sc.inertia_tensor = inertia_tensor

        args = build_config(
            spacecraft=sc,
            density_model=NoAtmosphereModel(),
            orientation_sim=true,
            mission_time=10.0,
            EI_km=120.0,
            dynamic_effectors=(effector,),
            keplerian=true
        )

        u0 = build_initial_conditions(args)
        u0.sc[1].pos .= pos_ii
        u0.sc[1].vel .= 0.0
        u0.sc[1].q .= q0
        u0.sc[1].ω .= ω0

        du0 = copy(u0)
        du0 .= 0.0
        p = ODEParams(n_sats=1, args=args)
        spacecraft_dynamics!(du0, u0, p, 0.0)
        return SVector{3, Float64}(du0.sc[1].ω)
    end

    @test isapprox(
        gravity_gradient_ωdot(InverseSquaredGravityModel(gravity_gradient=false)),
        SVector{3, Float64}(0.0, 0.0, 0.0);
        atol=1e-12,
        rtol=1e-10
    )

    for effector in (
        ConstantGravityModel(gravity_gradient=true),
        InverseSquaredGravityModel(gravity_gradient=true),
        InverseSquaredJ2GravityModel(gravity_gradient=true),
    )
        @test isapprox(gravity_gradient_ωdot(effector), ωdot_expected; atol=1e-12, rtol=1e-10)
    end
end

@testset "Orientation Simulation Rejects Invalid Inertia Tensor" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.01, -0.015, 0.02)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    sc.inertia_tensor = SMatrix{3, 3, Float64}(zeros(3, 3))

    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )

    @test_throws ArgumentError run_simulation(args)
end

@testset "Heat Loads Are Not Coupled To Force Magnitude" begin
    sc = make_single_link_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=175.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(ConstantForceModel(SVector{3, Float64}(1.0e5, -2.0e5, 3.0e5)),),
        keplerian=true
    )
    args.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))

    u0 = build_initial_conditions(args)
    du0 = copy(u0)
    du0 .= 0.0
    p = ODEParams(n_sats=1, args=args)
    spacecraft_dynamics!(du0, u0, p, 0.0)

    @test norm(SVector{3, Float64}(du0.sc[1].vel)) > 0.0
    @test all(==(0.0), du0.sc[1].heat_loads)
end

@testset "Atmosphere-Implicit Split RHS Recombines To Full RHS" begin
    q0_split = normalize(SVector{4, Float64}(0.15, -0.1, 0.25, 0.95))
    ω0_split = SVector{3, Float64}(0.01, -0.015, 0.02)
    sc_split = make_spacecraft(
        ra_alt_m=220e3,
        rp_alt_m=100e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=0.0,
        orientation_state=(q0_split, ω0_split)
    )
    explicit_torque = SVector{3, Float64}(0.12, -0.08, 0.05)
    args_split = build_config(
        spacecraft=sc_split,
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=true,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            AerodynamicCoefficientfM(),
            ConstantTorqueModel(explicit_torque),
        ),
        keplerian=true
    )
    args_split.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))

    u0_split = build_initial_conditions(args_split)
    du_full = copy(u0_split)
    du_implicit = copy(u0_split)
    du_explicit = copy(u0_split)
    du_sum = copy(u0_split)
    du_full .= 0.0
    du_implicit .= 0.0
    du_explicit .= 0.0
    du_sum .= 0.0
    p_split = ODEParams(n_sats=1, args=args_split)
    _initialize_heat_rate_buffers!(p_split)
    p_split.shared_buffers.in_atmosphere[1] = false
    p_split.shared_buffers.in_atmosphere_sample_t[1] = -1.0
    @test !SimulationEngine._drag_state_buffer_current(p_split, 1, 0.0)
    @test !SimulationEngine._all_active_spacecraft_outside_atmosphere(u0_split.sc, p_split, 0.0)

    spacecraft_dynamics!(du_full, u0_split, p_split, 0.0)
    spacecraft_dynamics_implicit_atmosphere!(du_implicit, u0_split, p_split, 0.0)
    spacecraft_dynamics_explicit_remainder!(du_explicit, u0_split, p_split, 0.0)
    du_sum .= du_implicit .+ du_explicit

    @test isapprox(du_sum.sc[1].pos, du_full.sc[1].pos; atol=1e-12, rtol=1e-10)
    @test isapprox(du_sum.sc[1].vel, du_full.sc[1].vel; atol=1e-12, rtol=1e-10)
    @test isapprox(du_sum.sc[1].mass, du_full.sc[1].mass; atol=1e-12, rtol=1e-10)
    @test isapprox(du_sum.sc[1].q, du_full.sc[1].q; atol=1e-12, rtol=1e-10)
    @test isapprox(du_sum.sc[1].ω, du_full.sc[1].ω; atol=1e-12, rtol=1e-10)
    @test isapprox(du_sum.sc[1].heat_loads, du_full.sc[1].heat_loads; atol=1e-12, rtol=1e-10)

    @test all(==(0.0), du_implicit.sc[1].pos)
    @test du_implicit.sc[1].mass == 0.0
    @test all(==(0.0), du_implicit.sc[1].q)
    @test all(==(0.0), du_implicit.sc[1].heat_loads)
    @test norm(SVector{3, Float64}(du_implicit.sc[1].vel)) > 0.0
    @test norm(SVector{3, Float64}(du_implicit.sc[1].ω)) > 0.0

    @test isapprox(
        SVector{3, Float64}(du_explicit.sc[1].pos),
        SVector{3, Float64}(u0_split.sc[1].vel);
        atol=1e-12,
        rtol=1e-10
    )
    @test du_explicit.sc[1].mass == 0.0
    @test any(>(0.0), du_explicit.sc[1].heat_loads)

    args_split_no_aero = build_config(
        spacecraft=sc_split,
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=true,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            ConstantTorqueModel(explicit_torque),
        ),
        keplerian=true
    )
    args_split_no_aero.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))
    u0_split_no_aero = build_initial_conditions(args_split_no_aero)
    du_explicit_no_aero = copy(u0_split_no_aero)
    du_explicit_no_aero .= 0.0
    p_split_no_aero = ODEParams(n_sats=1, args=args_split_no_aero)
    _initialize_heat_rate_buffers!(p_split_no_aero)
    spacecraft_dynamics_explicit_remainder!(du_explicit_no_aero, u0_split_no_aero, p_split_no_aero, 0.0)
    @test isapprox(du_explicit.sc[1].ω, du_explicit_no_aero.sc[1].ω; atol=1e-12, rtol=1e-10)
end

@testset "Atmosphere-Implicit Fast Path Requires Current Outside Decision" begin
    q0_outside = normalize(SVector{4, Float64}(0.1, 0.2, -0.05, 0.97))
    ω0_outside = SVector{3, Float64}(0.01, 0.0, -0.015)
    sc_outside = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        orientation_state=(q0_outside, ω0_outside)
    )
    args_outside = build_config(
        spacecraft=sc_outside,
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=true,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            AerodynamicCoefficientfM(),
        ),
        keplerian=true
    )
    args_outside.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))

    u0_outside = build_initial_conditions(args_outside)
    du_outside = copy(u0_outside)
    du_outside .= 0.0
    p_outside = ODEParams(n_sats=1, args=args_outside)
    _initialize_heat_rate_buffers!(p_outside)

    p_outside.shared_buffers.in_atmosphere[1] = true
    p_outside.shared_buffers.in_atmosphere_sample_t[1] = -1.0
    @test !SimulationEngine._drag_state_buffer_current(p_outside, 1, 0.0)
    @test SimulationEngine._all_active_spacecraft_outside_atmosphere(u0_outside.sc, p_outside, 0.0)

    spacecraft_dynamics_implicit_atmosphere!(du_outside, u0_outside, p_outside, 0.0)
    @test all(==(0.0), du_outside.sc[1].pos)
    @test all(==(0.0), du_outside.sc[1].vel)
    @test du_outside.sc[1].mass == 0.0
    @test all(==(0.0), du_outside.sc[1].q)
    @test all(==(0.0), du_outside.sc[1].ω)
    @test all(==(0.0), du_outside.sc[1].heat_loads)

    p_outside.shared_buffers.in_atmosphere[1] = false
    p_outside.shared_buffers.in_atmosphere_sample_t[1] = 0.0
    @test SimulationEngine._drag_state_buffer_current(p_outside, 1, 0.0)
    @test SimulationEngine._all_active_spacecraft_outside_atmosphere(u0_outside.sc, p_outside, 0.0)

    p_outside.shared_buffers.in_atmosphere[1] = true
    p_outside.shared_buffers.in_atmosphere_sample_t[1] = 0.0
    @test SimulationEngine._drag_state_buffer_current(p_outside, 1, 0.0)
    @test !SimulationEngine._all_active_spacecraft_outside_atmosphere(u0_outside.sc, p_outside, 0.0)

    sc_mixed = [
        make_spacecraft(ra_alt_m=220e3, rp_alt_m=100e3, ν_deg=0.0),
        make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=0.0),
    ]
    args_mixed = build_config_multi(
        spacecraft=sc_mixed,
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(AerodynamicCoefficientfM(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    args_mixed.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))
    u0_mixed = build_initial_conditions(args_mixed)
    du_mixed = copy(u0_mixed)
    du_mixed .= 0.0
    p_mixed = ODEParams(n_sats=2, args=args_mixed)
    _initialize_heat_rate_buffers!(p_mixed)
    p_mixed.shared_buffers.in_atmosphere[1] = false
    p_mixed.shared_buffers.in_atmosphere[2] = true
    p_mixed.shared_buffers.in_atmosphere_sample_t .= -1.0
    @test !SimulationEngine._all_active_spacecraft_outside_atmosphere(u0_mixed.sc, p_mixed, 0.0)

    withenv(
        "SPACEAGORA_RHS_EXECUTION_MODE" => "flat",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "off",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(max(1, Threads.nthreads()))
    ) do
        spacecraft_dynamics_implicit_atmosphere!(du_mixed, u0_mixed, p_mixed, 0.0)
    end
    @test norm(SVector{3, Float64}(du_mixed.sc[1].vel)) > 0.0
    @test all(==(0.0), du_mixed.sc[2].pos)
    @test all(==(0.0), du_mixed.sc[2].vel)
    @test du_mixed.sc[2].mass == 0.0
    @test all(==(0.0), du_mixed.sc[2].heat_loads)
end

@testset "Split IMEX Allows Zero Implicit Partition" begin
    args_zero_implicit = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    u0_zero_implicit = build_initial_conditions(args_zero_implicit)
    p_zero_implicit = ODEParams(n_sats=1, args=args_zero_implicit)
    _initialize_heat_rate_buffers!(p_zero_implicit)
    split_prob_zero_implicit = withenv("SPACEAGORA_SOLVER_MODE" => "split_imex") do
        _build_typed_solver_problem(u0_zero_implicit, (0.0, 30.0), p_zero_implicit, CallbackSet(), _solver_policy_mode())
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "split_imex", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol_zero_implicit, meta_zero_implicit = _solve_with_solver_policy(split_prob_zero_implicit, _active_solver_config(), args_zero_implicit, 1e-8, 1e-8)
        @test string(sol_zero_implicit.retcode) == "Success"
        @test meta_zero_implicit.solver == "KenCarp4(IMEX)"
    end
end

@testset "Heat Load Derivative Uses Thermal Model" begin
    sc = make_single_link_spacecraft(
        ra_alt_m=220e3,
        rp_alt_m=100e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=0.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=true
    )
    args.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))

    u0 = build_initial_conditions(args)
    du0 = copy(u0)
    du0 .= 0.0
    p = ODEParams(n_sats=1, args=args)
    _initialize_heat_rate_buffers!(p)
    spacecraft_dynamics!(du0, u0, p, 0.0)
    drag_force, drag_torque = SimulationModel.calcForceTorque(AerodynamicCoefficientfM(), u0.sc[1], p, 1)

    @test p.shared_buffers.densities[1] > 0.0
    @test isfinite(p.shared_buffers.temperatures[1])
    @test all(isfinite, p.shared_buffers.winds[1])
    @test norm(drag_force) > 0.0
    @test all(isfinite, drag_torque)
    @test all(isfinite, p.shared_buffers.heat_rates[1])
    @test any(>(0.0), p.shared_buffers.heat_rates[1])
    @test all(isfinite, du0.sc[1].heat_loads)
    @test any(>(0.0), du0.sc[1].heat_loads)
end

