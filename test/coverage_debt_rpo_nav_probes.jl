using Test
using LinearAlgebra
using Random
using StaticArrays

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel

if !isdefined(@__MODULE__, :RPOStationAssets)
    include(joinpath(REPO_ROOT, "src", "assets", "rpo_station_assets.jl"))
end

const SM = SimulationModel
const NH = SM.NavigationHooks
const CH = SM.ControlHooks
const CRA = SM.ClothRobotArmDynamics
const ASSETS = RPOStationAssets

const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH = Earth("", SPICE_PATH)

"""Minimal simulation configuration used only to type ODEParams arguments."""
function build_min_config()
    root = Link{0}(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra=EARTH.Rp_e + 500e3,
        rp=EARTH.Rp_e + 450e3,
        i=35.0,
        ω=40.0,
        Ω=10.0,
        ν=175.0
    )
    spacecraft = SpacecraftModel(
        Joint[],
        [root],
        root,
        true,
        root.m,
        0.0,
        root.inertia,
        0,
        0,
        ic,
        1
    )
    environment_model = EnvironmentModel(
        planet=EARTH,
        EI=120.0,
        density_model=NoAtmosphereModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
        topography=false,
        wind=false
    )
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=60.0,
            orientation_sim=false,
            num_steps_to_save=10
        ),
        environment_model=environment_model,
        dynamics_model=DynamicsModel([spacecraft], (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances()
    )
end

"""Brute-force nearest neighbor over a 3 x N point cloud."""
function brute_force_nearest(p::SVector{3, Float64}, points::Matrix{Float64})
    best_idx = 0
    best_d2 = Inf
    for j in 1:size(points, 2)
        q = SVector{3, Float64}(points[1, j], points[2, j], points[3, j])
        d2 = sum(abs2, p - q)
        if d2 < best_d2
            best_d2 = d2
            best_idx = j
        end
    end
    return best_idx, best_d2
end

@testset "RPO Nav Coverage Probes" begin

    # ---------------------------------------------------------------
    # station_geometry.jl
    # ---------------------------------------------------------------
    @testset "RPOStationGeometry construction" begin
        @test_throws ArgumentError SM.RPOStationGeometry(zeros(2, 4))
        @test_throws ArgumentError SM.RPOStationGeometry(zeros(3, 0))
        @test_throws ArgumentError SM.RPOStationGeometry(zeros(3, 2); keepout_radius_m=-1.0)

        single = SM.RPOStationGeometry(reshape([1.0, 2.0, 3.0], 3, 1); keepout_radius_m=0.5, name="single")
        @test single.name == "single"
        @test single.keepout_radius_m == 0.5
        @test single.kd_root !== nothing
        @test single.kd_root.idx == 1
        @test single.kd_root.axis == 1
        @test single.kd_root.left === nothing
        @test single.kd_root.right === nothing

        # Three collinear-in-x points: root should be the x-median (idx of x == 0.0).
        three = SM.RPOStationGeometry([-1.0 0.0 2.0; 0.0 0.0 0.0; 0.0 0.0 0.0])
        @test three.name == "station"
        @test three.keepout_radius_m == 0.0
        root = three.kd_root
        @test root.axis == 1
        @test three.points_body[1, root.idx] == 0.0
        @test root.left !== nothing
        @test root.right !== nothing
        @test three.points_body[1, root.left.idx] == -1.0
        @test three.points_body[1, root.right.idx] == 2.0
        # Children alternate to axis 2.
        @test root.left.axis == 2
        @test root.right.axis == 2

        # Integer matrix input is converted to Float64.
        int_geom = SM.RPOStationGeometry([0 1; 0 0; 0 0]; keepout_radius_m=1)
        @test int_geom.points_body isa Matrix{Float64}
        @test int_geom.keepout_radius_m === 1.0

        # Two-point cloud: mid = 1, so no left child and one right child.
        two = SM.RPOStationGeometry([0.0 5.0; 0.0 0.0; 0.0 0.0])
        @test two.kd_root.left === nothing
        @test two.kd_root.right !== nothing
    end

    # ---------------------------------------------------------------
    # mesh_distance.jl
    # ---------------------------------------------------------------
    @testset "nearest station point queries" begin
        pts = [0.0 2.0 -3.0; 0.0 0.0 4.0; 0.0 0.0 0.0]
        station = SM.RPOStationGeometry(pts)

        @test NH._rpo_station_point(station.points_body, 2) == SVector{3, Float64}(2.0, 0.0, 0.0)

        nearest, distance, idx = SM.nearest_station_point(SVector{3, Float64}(1.6, 0.2, 0.0), station)
        @test idx == 2
        @test nearest ≈ SVector{3, Float64}(2.0, 0.0, 0.0)
        @test distance ≈ sqrt(0.4^2 + 0.2^2)

        # Plain-Vector query path and squared-distance query.
        d2 = SM.nearest_station_distance_sq([1.6, 0.2, 0.0], station)
        @test d2 ≈ 0.4^2 + 0.2^2

        # Query exactly on a station point.
        nearest0, distance0, idx0 = SM.nearest_station_point([-3.0, 4.0, 0.0], station)
        @test idx0 == 3
        @test distance0 == 0.0
        @test SM.nearest_station_distance_sq([-3.0, 4.0, 0.0], station) == 0.0

        # Empty KD-tree error paths (bypass the validating constructor).
        empty_station = SM.RPOStationGeometry(Matrix{Float64}(undef, 3, 0), nothing, 0.0, "empty")
        @test_throws ArgumentError SM.nearest_station_point(SVector{3, Float64}(0.0, 0.0, 0.0), empty_station)
        @test_throws ArgumentError SM.nearest_station_distance_sq(SVector{3, Float64}(0.0, 0.0, 0.0), empty_station)

        # Randomized brute-force sweep exercises both near/far recursion branches.
        rng = MersenneTwister(2026)
        cloud = 4.0 .* (rand(rng, 3, 64) .- 0.5)
        rand_station = SM.RPOStationGeometry(cloud)
        for _ in 1:48
            p = SVector{3, Float64}(6.0 .* (rand(rng, 3) .- 0.5))
            bf_idx, bf_d2 = brute_force_nearest(p, cloud)
            q, dist, idx_kd = SM.nearest_station_point(p, rand_station)
            @test dist ≈ sqrt(bf_d2) atol=1e-12
            @test q == SVector{3, Float64}(cloud[1, idx_kd], cloud[2, idx_kd], cloud[3, idx_kd])
            @test sum(abs2, p - q) ≈ bf_d2 atol=1e-12
            @test SM.nearest_station_distance_sq(p, rand_station) ≈ bf_d2 atol=1e-12
        end
    end

    # ---------------------------------------------------------------
    # rpo_reference_geometry.jl + cubesat_geometry error path
    # ---------------------------------------------------------------
    origin_station = SM.RPOStationGeometry(reshape(zeros(3), 3, 1); keepout_radius_m=2.0)
    geom = SM.RPOReferenceGeometry(origin_station; chaser=SM.RPOCubeSatGeometry(dims_m=(0.2, 0.2, 0.3)))

    @testset "RPOReferenceGeometry" begin
        @test geom.station === origin_station
        @test geom.chaser.half_extents_body ≈ SVector{3, Float64}(0.1, 0.1, 0.15)

        default_geom = SM.RPOReferenceGeometry(origin_station)
        @test default_geom.chaser.half_extents_body ≈ SVector{3, Float64}(0.05, 0.05, 0.15)

        direct = SM.RPOReferenceGeometry(origin_station, SM.RPOCubeSatGeometry(dims_m=(0.4, 0.4, 0.4)))
        @test direct.chaser.half_extents_body ≈ SVector{3, Float64}(0.2, 0.2, 0.2)

        @test_throws ArgumentError SM.RPOCubeSatGeometry(dims_m=(0.0, 0.1, 0.1))
        @test_throws ArgumentError SM.RPOCubeSatGeometry(dims_m=(Inf, 0.1, 0.1))
    end

    # ---------------------------------------------------------------
    # clearance.jl
    # ---------------------------------------------------------------
    @testset "clearance queries" begin
        # Station point at origin, keepout 2.0, max chaser half extent 0.15.
        result = SM.rpo_clearance_to_station(SVector{3, Float64}(10.0, 0.0, 0.0), geom)
        @test result.distance ≈ 10.0
        @test result.clearance ≈ 10.0 - 2.0 - 0.15
        @test result.nearest_point == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test result.nearest_index == 1

        @test SM.rpo_clearance_distance_to_station(SVector{3, Float64}(10.0, 0.0, 0.0), geom) ≈ result.clearance
        @test SM.rpo_clearance_distance_to_station([0.0, 1.0, 0.0], geom) ≈ 1.0 - 2.0 - 0.15

        inside = SM.rpo_clearance_to_station([0.5, 0.0, 0.0], geom)
        @test inside.clearance < 0.0

        # Path stats: one violating sample out of four.
        path = [
            10.0 5.0 0.5 -8.0;
             0.0 0.0 0.0  0.0;
             0.0 0.0 0.0  0.0
        ]
        stats = SM.rpo_path_clearance_stats(path, geom)
        @test stats.violation_count == 1
        @test stats.violation_fraction ≈ 0.25
        @test stats.min_clearance ≈ 0.5 - 2.0 - 0.15

        # All-clear path.
        clear_stats = SM.rpo_path_clearance_stats([9.0 8.0; 0.0 0.0; 0.0 0.0], geom; safe_distance_m=0.5)
        @test clear_stats.violation_count == 0
        @test clear_stats.violation_fraction == 0.0
        @test clear_stats.min_clearance ≈ 8.0 - 2.0 - 0.15

        # BUG (tracked upstream): safe_distance_m is accepted but never used by
        # rpo_path_clearance_stats — violations are counted against clearance < 0
        # regardless of the margin.  The closest sample here has ~0.35 m clearance,
        # inside the requested 1.0 m margin, so a margin-aware check must count it
        # as a violation.  @test_broken asserts the CORRECT expectation: it records
        # Broken while the bug exists and flips to an unexpected-pass failure once
        # the source fix lands, so CI never rejects the fix.
        margin_stats = SM.rpo_path_clearance_stats([2.5 9.0; 0.0 0.0; 0.0 0.0], geom; safe_distance_m=1.0)
        @test margin_stats.min_clearance ≈ 2.5 - 2.0 - 0.15
        @test_broken margin_stats.violation_count == 1

        @test_throws ArgumentError SM.rpo_path_clearance_stats(zeros(4, 3), geom)
    end

    # ---------------------------------------------------------------
    # Station geometry assets (gated on files being present)
    # ---------------------------------------------------------------
    @testset "station geometry assets" begin
        demo_path = ASSETS.station_geometry_path(:demo)
        @test endswith(demo_path, joinpath("demo", "station_pointcloud.csv"))
        @test ASSETS.station_geometry_path(:gateway) == ASSETS.station_cad_path(:gateway)
        @test_throws ArgumentError ASSETS.station_geometry_path(:bogus)
        @test_throws ArgumentError ASSETS.station_cad_path(:bogus)
        @test_throws ArgumentError ASSETS.load_rpo_station_pointcloud(:gateway)

        if isfile(demo_path)
            demo_points = ASSETS.load_rpo_station_pointcloud(:demo)
            @test size(demo_points, 1) == 3
            @test size(demo_points, 2) > 0
        else
            @test_throws ArgumentError ASSETS.load_rpo_station_pointcloud(:demo)
        end

        cad_path = ASSETS.station_cad_path(:gateway)
        if isfile(cad_path)
            triangles = ASSETS.load_rpo_station_cad_triangles(:gateway)
            @test size(triangles, 1) == 3
            @test size(triangles, 2) % 3 == 0
            @test size(triangles, 2) > 0

            gateway_points = ASSETS.load_rpo_station_cad_pointcloud(:gateway; n_points=48, rng=MersenneTwister(740))
            @test size(gateway_points) == (3, 48)
            gateway_station = SM.RPOStationGeometry(gateway_points; keepout_radius_m=0.5, name="gateway")
            gateway_geom = SM.RPOReferenceGeometry(gateway_station)
            rng = MersenneTwister(7)
            for _ in 1:8
                p = SVector{3, Float64}(30.0 .* (rand(rng, 3) .- 0.5))
                bf_idx, bf_d2 = brute_force_nearest(p, gateway_points)
                _, dist, _ = SM.nearest_station_point(p, gateway_station)
                @test dist ≈ sqrt(bf_d2) atol=1e-9
            end
            far_clearance = SM.rpo_clearance_distance_to_station(SVector{3, Float64}(1.0e3, 0.0, 0.0), gateway_geom)
            @test far_clearance > 0.0
        end
    end

    # ---------------------------------------------------------------
    # Shared fixtures for control / cloth-arm tests
    # ---------------------------------------------------------------
    min_config = build_min_config()
    P1 = SM.ODEParams{1}(args=min_config)
    P2 = SM.ODEParams{2}(args=min_config)

    arm = SM.default_cloth_arm_model()
    base_pose = SM.ClothArmBasePose(SVector{3, Float64}(0.0, 0.0, 0.0), SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
    q_goal_true = [0.3, -0.2, 0.4]
    target = SM.cloth_fk(arm, base_pose, q_goal_true).end_effector_position
    plan = SM.plan_robot_arm_motion(
        arm,
        base_pose,
        zeros(3),
        target;
        config=SM.RobotArmPlannerConfig(dt_s=0.05, duration_s=1.0),
    )

    # ---------------------------------------------------------------
    # robot_arm_control.jl
    # ---------------------------------------------------------------
    @testset "robot arm joint MPC init" begin
        inertia = CH._robot_arm_default_joint_inertia(plan)
        @test length(inertia) == 3
        @test all(>(0.0), inertia)
        # Joint 1 carries every downstream link, so it dominates the wrist joint.
        @test inertia[1] > inertia[3]

        ctrl = SM.init_robot_arm_joint_mpc(plan; dt_s=0.05, horizon=4)
        @test ctrl.horizon == 4
        @test size(ctrl.H) == (12, 12)
        @test size(ctrl.E) == (12, 6)
        @test ctrl.H ≈ ctrl.H'  # symmetrized Hessian
        @test ctrl.U_prev == zeros(12)
        @test ctrl.joint_inertia_kg_m2 == inertia

        explicit = SM.init_robot_arm_joint_mpc(plan; joint_inertia_kg_m2=[0.5, 0.4, 0.3])
        @test explicit.joint_inertia_kg_m2 == [0.5, 0.4, 0.3]
        @test explicit.Bd[4, 1] ≈ 0.1 / 0.5  # dt / J1 with default dt 0.1

        @test_throws ArgumentError SM.init_robot_arm_joint_mpc(plan; horizon=0)
        @test_throws ArgumentError SM.init_robot_arm_joint_mpc(plan; dt_s=0.0)
        @test_throws ArgumentError SM.init_robot_arm_joint_mpc(plan; joint_inertia_kg_m2=[1.0, 2.0])
        @test_throws ArgumentError SM.init_robot_arm_joint_mpc(plan; joint_inertia_kg_m2=[1.0, 0.0, 1.0])
    end

    @testset "robot arm joint MPC preview and control" begin
        ctrl = SM.init_robot_arm_joint_mpc(plan; dt_s=0.05, horizon=4)

        preview = SM.robot_arm_joint_mpc_reference_preview(plan, 0.0, 0.05, 4)
        @test size(preview) == (6, 5)
        s0 = SM.robot_arm_plan_sample(plan, 0.0)
        @test preview[1:3, 1] ≈ s0.q
        @test preview[4:6, 1] ≈ s0.dq

        # Steady state at the plan goal is an equilibrium: command must vanish.
        x_eq = vcat(plan.q_goal, zeros(3))
        x_ref_eq = repeat(x_eq, 1, ctrl.horizon + 1)
        u_eq = SM.robot_arm_joint_mpc_control(ctrl, x_eq, x_ref_eq)
        @test length(u_eq) == 3
        @test norm(u_eq) < 1.0e-8
        @test norm(ctrl.U_prev) < 1.0e-8

        # Short reference triggers the padding branch and must agree at steady state.
        u_pad = SM.robot_arm_joint_mpc_control(ctrl, x_eq, reshape(x_eq, 6, 1))
        @test norm(u_pad) < 1.0e-8

        # Positive position error must produce a restoring (negative) acceleration.
        x_off = vcat(plan.q_goal .+ 0.1, zeros(3))
        u_off = SM.robot_arm_joint_mpc_control(ctrl, x_off, x_ref_eq)
        @test all(<(0.0), u_off)

        @test_throws ArgumentError SM.robot_arm_joint_mpc_control(ctrl, x_eq, zeros(5, 5))
    end

    @testset "robot arm measured state helpers" begin
        qc = CH._robot_arm_control_quat_conj(SVector{4, Float64}(0.1, 0.2, 0.3, 0.9))
        @test qc[4] > 0.0
        @test norm(qc) ≈ 1.0
        @test qc[1:3] ≈ -SVector{3, Float64}(0.1, 0.2, 0.3) / norm([0.1, 0.2, 0.3, 0.9]) atol=1e-12

        θ = 0.3
        q_rel = SVector{4, Float64}(0.0, 0.0, sin(θ / 2), cos(θ / 2))
        @test CH._robot_arm_control_axis_angle_about(q_rel, SVector{3, Float64}(0.0, 0.0, 1.0)) ≈ θ atol=1e-12
        # Negative-scalar quaternion is flipped before extraction.
        @test CH._robot_arm_control_axis_angle_about(-q_rel, SVector{3, Float64}(0.0, 0.0, 1.0)) ≈ θ atol=1e-12

        q_true = [0.2, -0.1, 0.3]
        pose_meas = SM.cloth_fk(arm, base_pose, q_true)
        arm_q_mat = zeros(4, 3)
        for i in 1:3
            arm_q_mat[:, i] .= pose_meas.link_quaternions[i]
        end
        sc_view = (
            q=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
            ω=SVector{3, Float64}(0.0, 0.0, 0.0),
            arm_q=arm_q_mat,
            arm_ω=zeros(3, 3),
        )
        x_meas = SM.robot_arm_measured_joint_state(plan, sc_view)
        @test x_meas !== nothing
        @test x_meas[1:3] ≈ q_true atol=1e-8
        @test x_meas[4:6] ≈ zeros(3) atol=1e-10

        # Missing arm state properties return nothing.
        @test SM.robot_arm_measured_joint_state(plan, (q=sc_view.q,)) === nothing
        @test SM.robot_arm_measured_joint_state(plan, (arm_q=arm_q_mat,)) === nothing

        # Without a base quaternion the plan base pose is used (identity here).
        x_no_base = SM.robot_arm_measured_joint_state(plan, (arm_q=arm_q_mat, arm_ω=zeros(3, 3)))
        @test x_no_base[1:3] ≈ q_true atol=1e-8

        ref_state = CH._robot_arm_control_reference_state(plan, 0.0)
        @test ref_state ≈ vcat(SM.robot_arm_plan_sample(plan, 0.0).q, SM.robot_arm_plan_sample(plan, 0.0).dq)

        u_mock = (sc=[sc_view],)
        @test CH._robot_arm_control_spacecraft_state(u_mock, 1) === sc_view
        @test CH._robot_arm_control_spacecraft_state(u_mock, 2) === nothing
        @test CH._robot_arm_control_spacecraft_state((foo=1,), 1) === nothing
    end

    @testset "robot arm control effector" begin
        held_default = CH.RobotArmHeldActuation()
        @test held_default.joint_torque_nm == Float64[]
        @test held_default.base_force_ii == SVector{3, Float64}(0.0, 0.0, 0.0)

        # Measured state at the goal, plan sampled far past its end: equilibrium.
        pose_goal = SM.cloth_fk(arm, base_pose, plan.q_goal)
        arm_q_goal = zeros(4, 3)
        for i in 1:3
            arm_q_goal[:, i] .= pose_goal.link_quaternions[i]
        end
        sc_goal = (
            q=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
            ω=SVector{3, Float64}(0.0, 0.0, 0.0),
            arm_q=arm_q_goal,
            arm_ω=zeros(3, 3),
        )
        eff = SM.RobotArmControlEffector(
            plan=plan,
            controller=SM.init_robot_arm_joint_mpc(plan; dt_s=0.05, horizon=4),
            control_dt_s=0.05,
        )
        @test SM.calcControlEffect!(eff, (sc=[sc_goal],), P1, 100.0, 1) === nothing
        @test length(eff.held.joint_torque_nm) == 3
        @test norm(eff.held.joint_torque_nm) < 1.0e-8
        @test eff.held.base_force_ii == SVector{3, Float64}(0.0, 0.0, 0.0)

        # Wrong satellite index leaves held untouched.
        eff.held = CH.RobotArmHeldActuation(joint_torque_nm=[9.0, 9.0, 9.0])
        SM.calcControlEffect!(eff, (sc=[sc_goal],), P1, 0.0, 2)
        @test eff.held.joint_torque_nm == [9.0, 9.0, 9.0]

        # No plan: early return.
        eff_no_plan = SM.RobotArmControlEffector()
        @test SM.calcControlEffect!(eff_no_plan, (sc=[sc_goal],), P1, 0.0, 1) === nothing
        @test eff_no_plan.held.joint_torque_nm == Float64[]

        # No controller: zero torque held command of plan dimension.
        eff_no_ctrl = SM.RobotArmControlEffector(plan=plan)
        SM.calcControlEffect!(eff_no_ctrl, (sc=[sc_goal],), P1, 0.0, 1)
        @test eff_no_ctrl.held.joint_torque_nm == zeros(3)

        # Measured state unavailable: falls back to the reference state (still equilibrium).
        eff_ref = SM.RobotArmControlEffector(
            plan=plan,
            controller=SM.init_robot_arm_joint_mpc(plan; dt_s=0.05, horizon=4),
            control_dt_s=0.05,
        )
        SM.calcControlEffect!(eff_ref, (sc=[(pos=zeros(3),)],), P1, 100.0, 1)
        @test norm(eff_ref.held.joint_torque_nm) < 1.0e-8

        # Held wrench accessors.
        f, τ = SM.calcControlForceTorque(eff, Float64[], P1, 1, 0.0)
        @test f == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test τ == SVector{3, Float64}(0.0, 0.0, 0.0)
        f2, τ2 = SM.calcControlForceTorque(eff, Float64[], P1, 2, 0.0)
        @test f2 == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test τ2 == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test SM.calcControlMassFlowRate(eff, Float64[], P1, 1, 0.0) == 0.0
    end

    # ---------------------------------------------------------------
    # rpo_mpc_control_model.jl
    # ---------------------------------------------------------------
    @testset "RPO MPC control model" begin
        r_t = SVector{3, Float64}(7000.0e3, 0.0, 0.0)
        v_t = SVector{3, Float64}(0.0, 7500.0, 0.0)
        r_c = r_t + SVector{3, Float64}(100.0, 0.0, 0.0)  # 100 m radial offset
        v_c = v_t

        chaser_view = (pos=r_c, vel=v_c, mass=50.0, q=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), ω=SVector{3, Float64}(0.0, 0.0, 0.0))
        target_view = (pos=r_t, vel=v_t, mass=1000.0, q=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), ω=SVector{3, Float64}(0.0, 0.0, 0.0))
        u_rpo = (sc=[chaser_view, target_view],)

        rp, vp = CH._rpo_control_state_pos_vel(u_rpo, 1)
        @test rp == r_c
        @test vp == v_c
        rp2, vp2 = CH._rpo_control_state_pos_vel(u_rpo, 2)
        @test rp2 == r_t

        # Early-return guards.
        model = SM.RPOMPCControlModel()
        @test SM.calcControlEffect!(model, u_rpo, P2, 0.0, 2) === nothing  # wrong sat
        @test SM.calcControlEffect!(model, u_rpo, P2, 0.0, 1) === nothing  # invalid plan
        model.plan_buffer.valid = true
        model.plan_buffer.plan.valid = true
        model.plan_buffer.updated_at_s = 0.0
        @test SM.calcControlEffect!(model, u_rpo, P2, 0.0, 1) === nothing  # no controller
        @test model.held.force_ii == SVector{3, Float64}(0.0, 0.0, 0.0)

        # Full control path: hold-at-origin plan with a radial offset.
        n_mm = sqrt(3.986004418e14 / (7000.0e3)^3)
        Q = Matrix{Float64}(I, 6, 6)
        R = 0.1 .* Matrix{Float64}(I, 3, 3)
        ctrl = CH.init_rpo_lqmpc(n_mm, 1.0, Q, R, 10.0 .* Q, 5)
        nsteps = 30
        rpo_plan = CH.RPOPlan(
            valid=true,
            t_ref_s=collect(0.0:1.0:(nsteps - 1)),
            r_ref_rtn=zeros(3, nsteps),
            v_ref_rtn=zeros(3, nsteps),
        )
        model.controller = ctrl
        model.plan_buffer.plan = rpo_plan
        model.plan_buffer.valid = true
        model.plan_buffer.updated_at_s = 0.0

        @test SM.calcControlEffect!(model, u_rpo, P2, 0.0, 1) === nothing
        @test all(isfinite, model.held.force_ii)
        @test norm(model.held.force_ii) > 0.0
        @test all(>=(0.0), model.held.thruster_forces_n)
        @test any(>(0.0), model.held.thruster_forces_n)
        # The commanded force must push the chaser back toward the target
        # (negative radial in RTN ~ negative inertial x here).
        @test model.held.force_ii[1] < 0.0
        # Default thrusters act through the COM: no thruster torque, no RW command.
        @test model.held.torque_body ≈ SVector{3, Float64}(0.0, 0.0, 0.0) atol=1e-12

        f, τ = SM.calcControlForceTorque(model, Float64[], P2, 1, 0.0)
        @test f == model.held.force_ii
        @test τ == model.held.torque_body
        f0, τ0 = SM.calcControlForceTorque(model, Float64[], P2, 2, 0.0)
        @test f0 == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test τ0 == SVector{3, Float64}(0.0, 0.0, 0.0)

        # Mass flow: negative when thrusters fire, zero for the other satellite.
        mdot = SM.calcControlMassFlowRate(model, Float64[], P2, 1, 0.0)
        @test mdot < 0.0
        expected_mdot = -sum(
            model.held.thruster_forces_n[j] > 0.0 ?
                model.held.thruster_forces_n[j] / (model.thrusters.isp_s[j] * 9.80665) : 0.0
            for j in 1:6
        )
        @test mdot ≈ expected_mdot
        @test SM.calcControlMassFlowRate(model, Float64[], P2, 2, 0.0) == 0.0

        idle = SM.RPOMPCControlModel()
        @test SM.calcControlMassFlowRate(idle, Float64[], P2, 1, 0.0) == 0.0
    end

    # ---------------------------------------------------------------
    # cloth_robot_arm_dynamics.jl
    # ---------------------------------------------------------------
    cloth_plan = SM.plan_robot_arm_motion(
        arm,
        base_pose,
        zeros(3),
        target;
        config=SM.RobotArmPlannerConfig(dt_s=0.1, duration_s=0.4),
    )

    @testset "cloth arm internal helpers" begin
        @test CRA._unit_quat(zeros(4)) == SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
        @test CRA._unit_quat([0.0, 0.0, 0.0, 2.0]) == SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
        @test CRA._quat_conj([0.1, 0.2, 0.3, 0.9])[4] ≈ 0.9 / norm([0.1, 0.2, 0.3, 0.9])
        qz = SVector{4, Float64}(0.0, 0.0, sin(0.2), cos(0.2))
        @test CRA._quat_mul(SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), qz) ≈ qz
        @test CRA._rot(SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)) ≈ Matrix{Float64}(I, 3, 3)

        # Axis-angle error: identity, finite rotation, sign flip, small-angle branch.
        @test CRA._axis_angle_error(SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)) == SVector{3, Float64}(0.0, 0.0, 0.0)
        θ = 0.4
        qx = SVector{4, Float64}(sin(θ / 2), 0.0, 0.0, cos(θ / 2))
        @test CRA._axis_angle_error(qx) ≈ SVector{3, Float64}(θ, 0.0, 0.0) atol=1e-12
        @test CRA._axis_angle_error(-qx) ≈ SVector{3, Float64}(θ, 0.0, 0.0) atol=1e-12

        # Compliance matrix normalization.
        @test CRA._compliance_matrix(2.0, :k) ≈ 2.0 .* Matrix{Float64}(I, 3, 3)
        @test CRA._compliance_matrix([1.0, 2.0, 3.0], :k) ≈ diagm([1.0, 2.0, 3.0])
        @test CRA._compliance_matrix((1.0, 2.0, 3.0), :k) ≈ diagm([1.0, 2.0, 3.0])
        m33 = [4.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 6.0]
        @test CRA._compliance_matrix(m33, :k) ≈ m33
        @test_throws ArgumentError CRA._compliance_matrix(zeros(2, 2), :k)
        @test_throws ArgumentError CRA._compliance_matrix([1.0, 2.0], :k)
        @test_throws ArgumentError CRA._compliance_matrix((1.0, 2.0), :k)
        @test_throws ArgumentError CRA._compliance_matrix("bad", :k)

        fills = CRA._joint_compliance_matrices(3.0, 3, :k)
        @test length(fills) == 3
        @test all(M -> M ≈ 3.0 .* Matrix{Float64}(I, 3, 3), fills)
        per_joint = CRA._joint_compliance_matrices([1.0, 2.0, 3.0], 3, :k)
        @test per_joint[2] ≈ 2.0 .* Matrix{Float64}(I, 3, 3)
        @test_throws ArgumentError CRA._joint_compliance_matrices([1.0, 2.0], 3, :k)

        # Cylindrical link inertia.
        link = arm.links[1]
        Jl = CRA._link_inertia(link)
        len = norm(link.vector_parent)
        @test Jl[1, 1] ≈ 0.5 * link.mass_kg * link.radius_m^2
        @test Jl[2, 2] ≈ (1 / 12) * link.mass_kg * (3 * link.radius_m^2 + len^2)
        @test Jl[3, 3] ≈ Jl[2, 2]

        v_off = CRA._body_offset_velocity_world(
            SVector{3, Float64}(1.0, 0.0, 0.0),
            SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
            SVector{3, Float64}(0.0, 0.0, 1.0),
            SVector{3, Float64}(0.0, 1.0, 0.0),
        )
        @test v_off ≈ SVector{3, Float64}(1.0 - 1.0, 0.0, 0.0) atol=1e-14
    end

    @testset "cloth arm reference and state packing" begin
        ref = SM.cloth_reference_state(cloth_plan, 0.0)
        @test ref.t_s == 0.0
        @test ref.end_effector_position ≈ SVector{3, Float64}(cloth_plan.ee_ref[:, 1])
        @test ref.state.pose.end_effector_position ≈ ref.end_effector_position
        @test ref.state.q ≈ cloth_plan.q_ref[:, 1]

        # Rest quaternions recompose the FK chain: q_child = q_parent ∘ rest.
        rests = SM.cloth_robot_arm_rest_quaternions(cloth_plan, 0.0)
        pose0 = SM.cloth_fk(arm, base_pose, cloth_plan.q_ref[:, 1])
        parent_q = base_pose.quaternion
        for i in 1:3
            recomposed = CRA._quat_mul(parent_q, rests[i])
            @test recomposed ≈ pose0.link_quaternions[i] atol=1e-12
            parent_q = pose0.link_quaternions[i]
        end

        x0 = SM.cloth_robot_arm_initial_state(cloth_plan)
        @test length(x0) == 13 * 3
        @test SM.cloth_robot_arm_end_effector(cloth_plan, x0) ≈ SVector{3, Float64}(cloth_plan.ee_ref[:, 1]) atol=1e-10

        shape = SM.coupled_cloth_robot_arm_state_shape(cloth_plan)
        @test size(shape.arm_r) == (3, 3)
        @test size(shape.arm_q) == (4, 3)
        @test size(shape.arm_v) == (3, 3)
        @test size(shape.arm_ω) == (3, 3)

        model_mb = SM.cloth_robot_arm_multibody(cloth_plan)
        @test length(model_mb.bodies) == 3
        @test length(model_mb.joints) == 3
        @test [j.parent for j in model_mb.joints] == [0, 1, 2]
        @test model_mb.joints[1].rest_child_parent_quat ≈ rests[1] atol=1e-12

        actuators = SM.cloth_robot_arm_actuators(cloth_plan; kp_n_m_rad=2.0, kd_n_m_s_rad=0.1)
        @test length(actuators) == 3
        @test actuators[2].joint == 2
        @test actuators[1].kp_n_m_rad ≈ 2.0 .* Matrix{Float64}(I, 3, 3)
    end

    @testset "coupled cloth arm state init and RHS" begin
        n = 3
        sc_view = (
            pos=SVector{3, Float64}(0.0, 0.0, 0.0),
            vel=SVector{3, Float64}(1.0, 2.0, 3.0),
            q=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
            ω=SVector{3, Float64}(0.0, 0.0, 0.0),
            arm_r=zeros(3, n),
            arm_q=zeros(4, n),
            arm_v=zeros(3, n),
            arm_ω=zeros(3, n),
        )
        @test SM.initialize_coupled_cloth_robot_arm_state!(sc_view, cloth_plan) === nothing
        pose0 = SM.cloth_fk(arm, base_pose, cloth_plan.q_ref[:, 1])
        for i in 1:n
            @test sc_view.arm_r[:, i] ≈ Vector(pose0.link_com_positions[i]) atol=1e-12
            @test sc_view.arm_q[:, i] ≈ Vector(pose0.link_quaternions[i]) atol=1e-12
            # Zero base rotation: every link inherits the base translational velocity.
            @test sc_view.arm_v[:, i] ≈ [1.0, 2.0, 3.0] atol=1e-12
            @test sc_view.arm_ω[:, i] ≈ zeros(3) atol=1e-12
        end
        # Missing arm state: no-op.
        @test SM.initialize_coupled_cloth_robot_arm_state!((pos=zeros(3),), cloth_plan) === nothing

        # RHS at the exact rest configuration with zero velocities is an equilibrium.
        rest_view = (
            pos=SVector{3, Float64}(0.0, 0.0, 0.0),
            vel=SVector{3, Float64}(0.0, 0.0, 0.0),
            q=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
            ω=SVector{3, Float64}(0.0, 0.0, 0.0),
            arm_r=zeros(3, n),
            arm_q=zeros(4, n),
            arm_v=zeros(3, n),
            arm_ω=zeros(3, n),
        )
        SM.initialize_coupled_cloth_robot_arm_state!(rest_view, cloth_plan)
        du_view = (arm_r=fill(NaN, 3, n), arm_q=fill(NaN, 4, n), arm_v=fill(NaN, 3, n), arm_ω=fill(NaN, 3, n))
        base_force = zeros(3)
        base_torque = zeros(3)
        @test SM.assign_coupled_cloth_robot_arm_rhs!(du_view, rest_view, cloth_plan, 0.0, base_force, base_torque) === nothing
        @test norm(du_view.arm_r) < 1.0e-10
        @test norm(du_view.arm_q) < 1.0e-10
        @test norm(du_view.arm_v) < 1.0e-8
        @test norm(du_view.arm_ω) < 1.0e-8
        @test norm(base_force) < 1.0e-8
        @test norm(base_torque) < 1.0e-8

        # Sampling the plan at its end makes the rest orientation differ from the
        # held configuration: compliance plus actuators produce angular accelerations.
        du_late = (arm_r=zeros(3, n), arm_q=zeros(4, n), arm_v=zeros(3, n), arm_ω=zeros(3, n))
        base_force_late = zeros(3)
        base_torque_late = zeros(3)
        actuators = SM.cloth_robot_arm_actuators(
            cloth_plan;
            kp_n_m_rad=5.0,
            kd_n_m_s_rad=0.2,
            torque_limit_n_m=1.0e-3,
            feedforward_torque_child_body=(1.0e-4, 0.0, 0.0),
        )
        SM.assign_coupled_cloth_robot_arm_rhs!(
            du_late,
            rest_view,
            cloth_plan,
            cloth_plan.t_ref_s[end],
            base_force_late,
            base_torque_late;
            joint_actuators=actuators,
        )
        @test norm(du_late.arm_ω) > 0.0
        @test all(isfinite, du_late.arm_ω)
        @test all(isfinite, base_torque_late)

        # Missing arm properties or an empty arm: early no-op returns.
        @test SM.assign_coupled_cloth_robot_arm_rhs!(du_view, (pos=zeros(3),), cloth_plan, 0.0, zeros(3), zeros(3)) === nothing
        empty_model = SM.ClothArmModel(SM.ClothArmLink[], SM.ClothArmJoint[], SVector{3, Float64}(0.0, 0.0, 0.0))
        empty_plan = SM.RobotArmPlan(
            empty_model,
            base_pose,
            [0.0],
            zeros(0, 1),
            zeros(0, 1),
            zeros(0, 1),
            zeros(3, 1),
            Float64[],
            Float64[],
            SVector{3, Float64}(0.0, 0.0, 0.0),
            0.0,
            :manual,
        )
        @test SM.assign_coupled_cloth_robot_arm_rhs!(du_view, (arm_r=zeros(3, 0), pos=zeros(3), vel=zeros(3)), empty_plan, 0.0, zeros(3), zeros(3)) === nothing
        @test SM.cloth_robot_arm_end_effector(empty_plan, Float64[]) == base_pose.position
    end

    @testset "cloth arm plan simulation" begin
        sim = SM.simulate_cloth_robot_arm_plan(cloth_plan; dt_s=0.1, duration_s=0.25)
        times = sim.trajectory.t_s
        @test times[1] == 0.0
        @test times[end] ≈ 0.25  # duration is appended when the grid falls short
        nt = length(times)
        @test size(sim.end_effector_positions) == (3, nt)
        @test size(sim.reference_end_effector_positions) == (3, nt)
        @test length(sim.tracking_error_m) == nt
        @test all(isfinite, sim.tracking_error_m)
        @test sim.tracking_error_m[1] ≈ 0.0 atol=1e-9
        @test size(sim.joint_compliance_torques_body) == (3, 3, nt)
        @test size(sim.joint_actuator_torques_body) == (3, 3, nt)
        @test size(sim.joint_total_torques_body) == (3, 3, nt)
        @test sim.joint_total_torques_body ≈ sim.joint_compliance_torques_body .+ sim.joint_actuator_torques_body atol=1e-9

        # RK4 stepper with explicit actuators.
        actuators = SM.cloth_robot_arm_actuators(cloth_plan; kp_n_m_rad=1.0, kd_n_m_s_rad=0.05)
        sim_rk4 = SM.simulate_cloth_robot_arm_plan(
            cloth_plan;
            dt_s=0.1,
            duration_s=0.2,
            integrator=:rk4,
            joint_actuators=actuators,
        )
        @test sim_rk4.trajectory.t_s[end] ≈ 0.2
        @test all(isfinite, sim_rk4.tracking_error_m)

        @test_throws ArgumentError SM.simulate_cloth_robot_arm_plan(cloth_plan; dt_s=0.1, duration_s=0.2, integrator=:euler)
    end
end
