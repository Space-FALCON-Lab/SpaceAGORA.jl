@testset "Orbital Elements Round-Trip Invariant" begin
    ic = InitialCondition(
        ra=EARTH.Rp_e + 550e3,
        rp=EARTH.Rp_e + 300e3,
        i=37.0,
        ω=55.0,
        Ω=25.0,
        ν=130.0
    )
    r, v = orbitalelemtorv(ic, EARTH)
    oe = rvtoorbitalelement(SVector{3, Float64}(r), SVector{3, Float64}(v), EARTH)

    @test isapprox(oe[1], ic.a; rtol=1e-10, atol=1e-6)
    @test isapprox(oe[2], ic.e; rtol=1e-10, atol=1e-10)
    @test angle_distance(oe[3], ic.i) < 1e-8
    @test angle_distance(oe[4], ic.Ω) < 1e-8
    @test angle_distance(oe[5], ic.ω) < 1e-8
    @test angle_distance(oe[6], ic.ν) < 1e-8
end

@testset "Reference System Branch Coverage" begin
    ic_arg_wrap = InitialCondition(
        ra=EARTH.Rp_e + 700e3,
        rp=EARTH.Rp_e + 300e3,
        i=55.0,
        ω=300.0,
        Ω=40.0,
        ν=25.0
    )
    r_wrap, v_wrap = orbitalelemtorv(ic_arg_wrap, EARTH)
    oe_wrap = rvtoorbitalelement(SVector{3, Float64}(r_wrap), SVector{3, Float64}(v_wrap), EARTH)
    @test oe_wrap[5] > pi

    ic_circ_incl = InitialCondition(
        ra=EARTH.Rp_e + 500e3,
        rp=EARTH.Rp_e + 500e3,
        i=40.0,
        ω=0.0,
        Ω=15.0,
        ν=250.0
    )
    r_circ, v_circ = orbitalelemtorv(ic_circ_incl, EARTH)
    oe_circ = rvtoorbitalelement(SVector{3, Float64}(r_circ), SVector{3, Float64}(v_circ), EARTH)
    @test abs(oe_circ[2]) < 1e-10
    @test oe_circ[6] > pi

    planet_topo = (
        Rp_e=10.0,
        Rp_p=9.0,
        topography_function=(a, Clm, Slm, lat, lon, A) -> 9.5,
        Clm_topo=zeros(1, 1),
        Slm_topo=zeros(1, 1),
        A_topo=0.0
    )
    rp_topo = SVector{3, Float64}(0.1, 0.5, -8.944723618090451)
    lla = rtolatlong(rp_topo, planet_topo)
    @test all(isfinite, lla)

    topo_module_name = gensym(:ReferenceSystemTopoSandbox)
    Core.eval(Main, :(module $topo_module_name
        using LinearAlgebra
        using StaticArrays
        args = :legacy_topography_args
        include(joinpath(Main.REPO_ROOT, "src", "core", "interfaces", "reference_system.jl"))
    end))
    topo_sandbox = getfield(Main, topo_module_name)
    lla_topo = topo_sandbox.rtolatlong(rp_topo, planet_topo, true)
    @test isapprox(lla_topo[1], norm(rp_topo) - 9.5; atol=1e-12, rtol=0.0)

    q_branch2 = orbital_elements_to_lvlh_quaternion(
        1.186823891356144,
        0.13305160284932352,
        0.0,
        2.0245819323134224
    )
    q_branch3 = orbital_elements_to_lvlh_quaternion(
        0.0,
        0.1,
        0.0,
        0.10471975511965977
    )
    @test isapprox(norm(q_branch2), 1.0; atol=1e-12, rtol=0.0)
    @test isapprox(norm(q_branch3), 1.0; atol=1e-12, rtol=0.0)
end

@testset "Circular Orbit Robustness" begin
    ic_eq = InitialCondition(
        ra=EARTH.Rp_e + 500e3,
        rp=EARTH.Rp_e + 500e3,
        i=0.0,
        ω=0.0,
        Ω=0.0,
        ν=30.0
    )
    r, v = orbitalelemtorv(ic_eq, EARTH)
    oe = rvtoorbitalelement(SVector{3, Float64}(r), SVector{3, Float64}(v), EARTH)
    @test all(isfinite, oe)
    @test abs(oe[2]) < 1e-10

    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=0.0,
        ω_deg=0.0,
        Ω_deg=0.0,
        ν_deg=30.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=300.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    df = run_case(args)
    @test nrow(df) > 10
end

@testset "Quaternion Norm Invariant" begin
    q0 = 2.3 * normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.01, -0.015, 0.02)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=420e3,
        i_deg=40.0,
        ω_deg=10.0,
        Ω_deg=20.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-9,
            abstol_orbit=1e-9,
            reltol_quaternion=1e-9,
            abstol_quaternion=1e-9,
            dt_max_orbit=5.0
        )
    )

    u0 = build_initial_conditions(args)
    @test isapprox(norm(SVector{4, Float64}(u0.sc[1].q)), 1.0; atol=1e-12, rtol=0.0)

    df = run_case(args)
    qnorm = sqrt.(df.sc1_q_1.^2 .+ df.sc1_q_2.^2 .+ df.sc1_q_3.^2 .+ df.sc1_q_4.^2)
    @test maximum(abs.(qnorm .- 1.0)) < 1e-6
end

@testset "DynamicsRotational Core Contracts" begin
    rotational = SimulationModel.DynamicsRotational

    q = normalize(SVector{4, Float64}(0.15, -0.25, 0.35, 0.85))
    ω = SVector{3, Float64}(0.02, -0.03, 0.01)

    qdot_new = rotational.quaternion_derivative(ω, q)
    qdot_legacy = 0.5 * SVector{4, Float64}(SimulationModel.quat_mult(SVector{4, Float64}(ω..., 0.0), q))
    @test isapprox(qdot_new, qdot_legacy; atol=1e-12, rtol=1e-12)
    @test abs(dot(q, qdot_new)) < 1e-12

    J = SMatrix{3, 3, Float64}(2.0, 0.1, 0.0, 0.1, 3.0, 0.0, 0.0, 0.0, 4.0)
    τ = SVector{3, Float64}(0.12, -0.08, 0.05)
    ωdot_new = rotational.angular_acceleration(ω, J, τ; include_gyroscopic=true)
    ωdot_legacy = J \ (τ - cross(ω, J * ω))
    @test isapprox(ωdot_new, ωdot_legacy; atol=1e-12, rtol=1e-12)

    ωdot_no_gyro = rotational.angular_acceleration(ω, J, τ; include_gyroscopic=false)
    @test isapprox(ωdot_no_gyro, J \ τ; atol=1e-12, rtol=1e-12)

    @test rotational.body_torque(MVector{3, Float64}(τ)) == τ
    @test rotational.body_angular_velocity(MVector{3, Float64}(ω)) == ω
    @test rotational.combine_torques(τ, -τ) == SVector{3, Float64}(0.0, 0.0, 0.0)

    J_spherical = SMatrix{3, 3, Float64}(2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0)
    ω_any = SVector{3, Float64}(0.3, -0.25, 0.11)
    ωdot_torque_free = rotational.angular_acceleration(
        ω_any,
        J_spherical,
        SVector{3, Float64}(0.0, 0.0, 0.0);
        include_gyroscopic=true
    )
    @test isapprox(ωdot_torque_free, SVector{3, Float64}(0.0, 0.0, 0.0); atol=1e-12, rtol=0.0)
end

@testset "Torque-Free Rigid-Body Regression Cases" begin
    function make_torque_free_spacecraft(
        inertia_tensor::SMatrix{3, 3, Float64},
        q0::SVector{4, Float64},
        w0::SVector{3, Float64};
        id::Int64=1
    )
        root = Link{0}(root=true, m=500.0, ref_area=12.0, inertia=inertia_tensor)
        ic = InitialCondition(
            ra=EARTH.Rp_e + 500e3,
            rp=EARTH.Rp_e + 500e3,
            i=0.0,
            ω=0.0,
            Ω=0.0,
            ν=0.0,
            q=q0,
            ang_vel=w0
        )
        return SpacecraftModel(
            joints=Joint[],
            links=[root],
            root=root,
            prop_mass=0.0,
            inertia_tensor=inertia_tensor,
            initial_condition=ic,
            id=id
        )
    end

    function run_torque_free_case(
        spacecraft::SpacecraftModel;
        mission_time::Float64,
        dt_max_orbit::Float64=0.5
    )
        args = build_config(
            spacecraft=spacecraft,
            density_model=NoAtmosphereModel(),
            orientation_sim=true,
            mission_time=mission_time,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true,
            simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
            tolerances=IntegrationTolerances(
                reltol_orbit=1e-10,
                abstol_orbit=1e-10,
                reltol_quaternion=1e-10,
                abstol_quaternion=1e-10,
                reltol_angular_rate=1e-10,
                abstol_angular_rate=1e-10,
                dt_max_orbit=dt_max_orbit
            ),
            ephemerides_model=SimpleEphemeridesModel()
        )
        return run_simulation(args; return_solution=true)
    end

    withenv("SPACEAGORA_SOLVER_MODE" => "tsit5") do
        J_axis = SMatrix{3, 3, Float64}(8.0, 0.0, 0.0, 0.0, 12.0, 0.0, 0.0, 0.0, 20.0)
        q_axis0 = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
        ω_axis0 = SVector{3, Float64}(0.0, 0.0, 0.2)
        sol_axis = run_torque_free_case(
            make_torque_free_spacecraft(J_axis, q_axis0, ω_axis0);
            mission_time=40.0,
            dt_max_orbit=0.2
        )

        @test length(sol_axis.t) > 10
        for (t, state) in zip(sol_axis.t, sol_axis.u)
            q_state = SimulationModel.project_unit_quaternion(state.sc[1].q)
            ω_state = SVector{3, Float64}(state.sc[1].ω)
            q_expected = SVector{4, Float64}(0.0, 0.0, sin(0.5 * ω_axis0[3] * Float64(t)), cos(0.5 * ω_axis0[3] * Float64(t)))

            @test isapprox(ω_state, ω_axis0; atol=1e-9, rtol=1e-9)
            @test abs(dot(q_state, q_expected)) > 1.0 - 5e-8
        end

        J_free = SMatrix{3, 3, Float64}(8.0, 0.0, 0.0, 0.0, 11.0, 0.0, 0.0, 0.0, 15.0)
        q_free0 = normalize(SVector{4, Float64}(0.2, -0.3, 0.1, 0.92))
        ω_free0 = SVector{3, Float64}(0.19, -0.27, 0.31)
        sol_free = run_torque_free_case(
            make_torque_free_spacecraft(J_free, q_free0, ω_free0);
            mission_time=180.0
        )

        energies = Float64[]
        angular_momentum_norms = Float64[]
        quaternion_norm_errors = Float64[]
        for state in sol_free.u
            ω_state = SVector{3, Float64}(state.sc[1].ω)
            push!(energies, 0.5 * dot(ω_state, J_free * ω_state))
            push!(angular_momentum_norms, norm(J_free * ω_state))
            push!(quaternion_norm_errors, abs(norm(SVector{4, Float64}(state.sc[1].q)) - 1.0))
        end

        energy_drift = (maximum(energies) - minimum(energies)) / abs(first(energies))
        angular_momentum_drift = (maximum(angular_momentum_norms) - minimum(angular_momentum_norms)) / first(angular_momentum_norms)

        @test energy_drift < 1e-7
        @test angular_momentum_drift < 1e-7
        @test maximum(quaternion_norm_errors) < 1e-10
    end
end

@testset "Solver Tolerances Apply Quaternion Overrides" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.01, -0.015, 0.02)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=420e3,
        i_deg=40.0,
        ω_deg=10.0,
        Ω_deg=20.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    tols = IntegrationTolerances(
        reltol_orbit=1e-5,
        abstol_orbit=2e-6,
        reltol_quaternion=3e-7,
        abstol_quaternion=4e-8,
        reltol_mass=5e-7,
        abstol_mass=6e-8,
        reltol_heat_load=7e-7,
        abstol_heat_load=8e-8,
        reltol_angular_rate=9e-7,
        abstol_angular_rate=1e-7,
        dt_max_orbit=5.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=tols
    )

    u0 = build_initial_conditions(args)
    reltol_vec, abstol_vec = _build_solver_tolerances(u0, args)
    @test all(isapprox.(reltol_vec.sc[1].pos, tols.reltol_orbit; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec.sc[1].pos, tols.abstol_orbit; atol=0.0, rtol=0.0))
    @test all(isapprox.(reltol_vec.sc[1].q, tols.reltol_quaternion; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec.sc[1].q, tols.abstol_quaternion; atol=0.0, rtol=0.0))
    @test all(isapprox.(reltol_vec.sc[1].ω, tols.reltol_angular_rate; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec.sc[1].ω, tols.abstol_angular_rate; atol=0.0, rtol=0.0))
    @test isapprox(reltol_vec.sc[1].mass, tols.reltol_mass; atol=0.0, rtol=0.0)
    @test isapprox(abstol_vec.sc[1].mass, tols.abstol_mass; atol=0.0, rtol=0.0)
    @test all(isapprox.(reltol_vec.sc[1].heat_loads, tols.reltol_heat_load; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec.sc[1].heat_loads, tols.abstol_heat_load; atol=0.0, rtol=0.0))

    tols_no_orient = IntegrationTolerances(
        reltol_orbit=1e-5,
        abstol_orbit=2e-6,
        reltol_mass=0.0,
        abstol_mass=0.0,
        reltol_heat_load=0.0,
        abstol_heat_load=0.0,
        reltol_angular_rate=0.0,
        abstol_angular_rate=0.0,
        dt_max_orbit=5.0
    )
    args_no_orient = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=tols_no_orient
    )
    u0_no_orient = build_initial_conditions(args_no_orient)
    reltol_scalar, abstol_scalar = _build_solver_tolerances(u0_no_orient, args_no_orient)
    @test reltol_scalar == tols_no_orient.reltol_orbit
    @test abstol_scalar == tols_no_orient.abstol_orbit

    tols_no_orient_component = IntegrationTolerances(
        reltol_orbit=1e-5,
        abstol_orbit=2e-6,
        reltol_mass=3e-7,
        abstol_mass=4e-8,
        reltol_heat_load=5e-7,
        abstol_heat_load=6e-8,
        dt_max_orbit=5.0
    )
    args_no_orient_component = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=tols_no_orient_component
    )
    u0_no_orient_component = build_initial_conditions(args_no_orient_component)
    reltol_vec_no_orient, abstol_vec_no_orient = _build_solver_tolerances(u0_no_orient_component, args_no_orient_component)
    @test isapprox(reltol_vec_no_orient.sc[1].mass, tols_no_orient_component.reltol_mass; atol=0.0, rtol=0.0)
    @test isapprox(abstol_vec_no_orient.sc[1].mass, tols_no_orient_component.abstol_mass; atol=0.0, rtol=0.0)
    @test all(isapprox.(reltol_vec_no_orient.sc[1].heat_loads, tols_no_orient_component.reltol_heat_load; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec_no_orient.sc[1].heat_loads, tols_no_orient_component.abstol_heat_load; atol=0.0, rtol=0.0))

    tols_ω_no_orient = IntegrationTolerances(
        reltol_orbit=1e-5,
        abstol_orbit=2e-6,
        reltol_mass=0.0,
        abstol_mass=0.0,
        reltol_heat_load=0.0,
        abstol_heat_load=0.0,
        reltol_angular_rate=3e-7,
        dt_max_orbit=5.0
    )
    args_ω_no_orient = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=tols_ω_no_orient
    )
    u0_ω_no_orient = build_initial_conditions(args_ω_no_orient)
    reltol_ω_no_orient, abstol_ω_no_orient = _build_solver_tolerances(u0_ω_no_orient, args_ω_no_orient)
    @test reltol_ω_no_orient == tols_ω_no_orient.reltol_orbit
    @test abstol_ω_no_orient == tols_ω_no_orient.abstol_orbit
end
