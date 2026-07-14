using Test
using LinearAlgebra
using StaticArrays
# SpaceAGORA is already bound at Main scope by test/helpers/bootstrap.jl's
# raw-include of src/SpaceAGORA.jl; `using SpaceAGORA` here would conflict
# with that binding (raw-include and package-load can't coexist in one
# process -- see test/TEST_RESTRUCTURE_PLAN.md Phase 2 notes).
using Random
using ComponentArrays

const ROBOTICS_PLANNING_CASES = (
    (
        name=:cross_body_sweep,
        q_start=[0.0, 0.0, 0.0],
        q_goal=[0.25, -0.35, 0.20],
    ),
    (
        name=:folded_reach,
        q_start=[-0.30, 0.20, -0.15],
        q_goal=[0.32, -0.28, 0.18],
    ),
    (
        name=:reverse_elbow,
        q_start=[0.42, -0.18, 0.22],
        q_goal=[-0.18, 0.34, -0.26],
    ),
    (
        name=:x_reach_change,
        q_start=[0.08, 0.95, -0.85],
        q_goal=[-0.12, -0.08, 0.06],
    ),
)

const ROBOTICS_HYPR_CASES = (
    (
        name=:mid_path_block,
        q_start=[0.0, -0.25, 0.35],
        q_goal=[0.35, -0.65, 0.30],
        obstacles=(
            (fraction=0.34, offset=SVector{3, Float64}(-0.035, 0.015, 0.010), radius=0.020),
            (fraction=0.68, offset=SVector{3, Float64}(0.040, -0.006, 0.004), radius=0.014),
        ),
    ),
    (
        name=:offset_corridor,
        q_start=[-0.28, 0.72, -0.58],
        q_goal=[0.22, -0.10, 0.08],
        obstacles=(
            (fraction=0.24, offset=SVector{3, Float64}(-0.050, -0.012, 0.008), radius=0.016),
            (fraction=0.58, offset=SVector{3, Float64}(0.055, 0.014, -0.004), radius=0.018),
            (fraction=0.82, offset=SVector{3, Float64}(-0.030, 0.006, 0.012), radius=0.012),
        ),
    ),
    (
        name=:x_reach_corridor,
        q_start=[0.18, 1.05, -0.92],
        q_goal=[-0.16, -0.14, 0.10],
        obstacles=(
            (fraction=0.20, offset=SVector{3, Float64}(0.060, 0.010, 0.006), radius=0.014),
            (fraction=0.52, offset=SVector{3, Float64}(-0.065, -0.010, 0.010), radius=0.016),
            (fraction=0.86, offset=SVector{3, Float64}(0.045, 0.012, -0.006), radius=0.012),
        ),
    ),
)

function _robotics_path_point(path::AbstractMatrix, fraction::Real)
    n = size(path, 2)
    n == 1 && return SVector{3, Float64}(path[:, 1])
    u = clamp(Float64(fraction), 0.0, 1.0) * (n - 1)
    lo = clamp(floor(Int, u) + 1, 1, n)
    hi = min(lo + 1, n)
    α = clamp(u - (lo - 1), 0.0, 1.0)
    return SVector{3, Float64}((1.0 - α) .* path[:, lo] .+ α .* path[:, hi])
end

function _robotics_obstacles_from_nominal(plan, specs)
    return SpaceAGORA.RobotArmSphereObstacle[
        SpaceAGORA.RobotArmSphereObstacle(
            _robotics_path_point(plan.ee_ref, spec.fraction) + spec.offset,
            spec.radius,
        ) for spec in specs
    ]
end

@testset "Cloth-arm FK/IK planning bridge" begin
    model = SpaceAGORA.default_cloth_arm_model()
    base = SpaceAGORA.ClothArmBasePose([0.0, 0.0, 0.0])
    q0 = [0.0, 0.0, 0.0]

    pose0 = SpaceAGORA.cloth_fk(model, base, q0)
    @test pose0.end_effector_position ≈ [0.46, 0.0, 0.0] atol=1e-12
    @test SpaceAGORA.cloth_total_reach(model) ≈ 0.46

    q = ROBOTICS_PLANNING_CASES[1].q_goal
    target = SpaceAGORA.cloth_fk(model, base, q).end_effector_position
    q_sol = SpaceAGORA.cloth_ik(model, base, target; q_seed=q0, position_tol_m=5.0e-4)
    solved = SpaceAGORA.cloth_fk(model, base, q_sol)
    @test norm(solved.end_effector_position - target) < 5.0e-3

    cfg = SpaceAGORA.RobotArmPlannerConfig(dt_s=0.25, duration_s=1.0, ik_tol_m=5.0e-4)
    plan = SpaceAGORA.plan_robot_arm_motion(model, base, q0, target; config=cfg)
    @test size(plan.q_ref) == (3, 5)
    @test plan.t_ref_s[end] == 1.0
    @test plan.final_error_m < 5.0e-3

    sample = SpaceAGORA.robot_arm_plan_sample(plan, 0.5)
    @test length(sample.q) == 3
    @test length(sample.dq) == 3
    @test length(sample.ddq) == 3

    ref_state = SpaceAGORA.cloth_reference_state(plan, 0.5)
    ee, _ = SpaceAGORA.cloth_end_effector_pose(ref_state.state)
    @test ee ≈ sample.ee

    mpc = SpaceAGORA.init_robot_arm_joint_mpc(plan; dt_s=0.1, horizon=4)
    preview = SpaceAGORA.robot_arm_joint_mpc_reference_preview(plan, 0.0, 0.1, mpc.horizon)
    τ = SpaceAGORA.robot_arm_joint_mpc_control(mpc, preview[:, 1], preview)
    @test length(τ) == length(q0)
    @test all(isfinite, τ)

    shape = merge((
        pos=zeros(3),
        vel=zeros(3),
        mass=10.0,
        heat_loads=zeros(1),
        q=Float64[0.0, 0.0, 0.0, 1.0],
        ω=zeros(3),
    ), SpaceAGORA.coupled_cloth_robot_arm_state_shape(plan))
    sc = ComponentVector(shape)
    SpaceAGORA.initialize_coupled_cloth_robot_arm_state!(sc, plan)
    x_meas = SpaceAGORA.robot_arm_measured_joint_state(plan, sc)
    @test x_meas[1:length(q0)] ≈ plan.q_ref[:, 1] atol=1.0e-10
    @test x_meas[(length(q0) + 1):end] ≈ plan.dq_ref[:, 1] atol=1.0e-10

    for case in ROBOTICS_PLANNING_CASES[2:end]
        target_i = SpaceAGORA.cloth_fk(model, base, case.q_goal).end_effector_position
        plan_i = SpaceAGORA.plan_robot_arm_motion(
            model,
            base,
            case.q_start,
            target_i;
            config=cfg,
        )
        @test plan_i.q_start ≈ case.q_start
        @test plan_i.final_error_m < 5.0e-3
        @test all(isfinite, plan_i.q_ref)
        @test all(isfinite, plan_i.ee_ref)
    end
end

@testset "Robot arm HYPR planner" begin
    model = SpaceAGORA.default_cloth_arm_model()
    base = SpaceAGORA.ClothArmBasePose([0.0, 0.0, 0.0])
    primary_case = ROBOTICS_HYPR_CASES[1]
    q0 = primary_case.q_start
    q_goal_truth = primary_case.q_goal
    planner_cfg = SpaceAGORA.RobotArmPlannerConfig(dt_s=0.1, duration_s=0.5, ik_tol_m=1.0e-3)
    hypr_cfg = SpaceAGORA.RobotArmHYPRConfig(
        n_waypoints=2,
        n_particles=18,
        n_iters=6,
        n_samples=14,
        safe_distance_m=0.004,
        w_obs=250.0,
        early_stopping_enable=false,
    )
    target = SpaceAGORA.cloth_fk(model, base, q_goal_truth).end_effector_position
    obstacle = nothing
    for (case_idx, case) in pairs(ROBOTICS_HYPR_CASES)
        case_target = SpaceAGORA.cloth_fk(model, base, case.q_goal).end_effector_position
        nominal = SpaceAGORA.plan_robot_arm_motion(
            model,
            base,
            case.q_start,
            case_target;
            config=planner_cfg,
        )
        obstacles = _robotics_obstacles_from_nominal(nominal, case.obstacles)
        case_idx == 1 && (obstacle = obstacles[1])
        result = SpaceAGORA.plan_robot_arm_motion_hypr(
            model,
            base,
            case.q_start,
            case_target;
            planner_config=planner_cfg,
            hypr_config=hypr_cfg,
            obstacles=obstacles,
            rng=MersenneTwister(41 + case_idx),
        )
        @test result.plan.planner == :hypr
        @test size(result.control_points) == (3, 4)
        @test size(result.sampled_path) == (3, hypr_cfg.n_samples)
        @test length(result.cost_history) == hypr_cfg.n_iters
        @test isfinite(result.cost)
        @test isfinite(result.components.J_len_norm)
        @test result.config.refinement_enable === true
        @test result.components.refinement_enabled === true
        @test result.components.refinement_rounds >= 0
        @test result.plan.final_error_m < 5.0e-3
    end

    warmstart_cfg = SpaceAGORA.RobotArmHYPRConfig(
        n_waypoints=2,
        n_particles=8,
        n_iters=2,
        n_samples=10,
        rrt_warmstart_enable=true,
        rrt_warmstart_iters=6,
        rrt_warmstart_shortcut_iters=2,
        early_stopping_enable=false,
    )
    warmstarted = SpaceAGORA.plan_robot_arm_motion_hypr(
        model,
        base,
        q0,
        target;
        planner_config=SpaceAGORA.RobotArmPlannerConfig(dt_s=0.1, duration_s=0.5, ik_tol_m=1.0e-3),
        hypr_config=warmstart_cfg,
        obstacles=SpaceAGORA.RobotArmSphereObstacle[],
        rng=MersenneTwister(46),
    )
    @test warmstarted.components.rrt_warmstart_enabled === true
    @test warmstarted.components.rrt_warmstart_attempted === true
    @test warmstarted.components.rrt_warmstart_path_found === true
    @test warmstarted.components.rrt_warmstart_n_points == 2
    @test warmstarted.control_points[:, 1] ≈ q0
    @test warmstarted.control_points[:, end] ≈ warmstarted.plan.q_goal

    retimed_cfg = SpaceAGORA.RobotArmHYPRConfig(
        n_waypoints=0,
        n_particles=4,
        n_iters=1,
        n_samples=8,
        retime_enable=true,
        retime_max_joint_velocity_rad_s=0.2,
        retime_max_joint_acceleration_rad_s2=0.8,
        early_stopping_enable=false,
    )
    retimed = SpaceAGORA.plan_robot_arm_motion_hypr(
        model,
        base,
        q0,
        target;
        planner_config=SpaceAGORA.RobotArmPlannerConfig(dt_s=0.1, duration_s=0.5, ik_tol_m=1.0e-3),
        hypr_config=retimed_cfg,
        obstacles=SpaceAGORA.RobotArmSphereObstacle[],
        rng=MersenneTwister(44),
    )
    @test retimed.plan.t_ref_s[end] > 0.5
    @test retimed.components.retime_enabled === true
    @test retimed.components.retime_duration_s == retimed.plan.t_ref_s[end]
    @test maximum(abs.(retimed.plan.dq_ref)) <= retimed_cfg.retime_max_joint_velocity_rad_s + 1.0e-8

    cloth_retimed_target = SpaceAGORA.cloth_fk(model, base, [0.08, -0.12, 0.08]).end_effector_position
    cloth_retimed = SpaceAGORA.plan_robot_arm_motion_hypr(
        model,
        base,
        [0.0, -0.05, 0.05],
        cloth_retimed_target;
        planner_config=SpaceAGORA.RobotArmPlannerConfig(dt_s=0.1, duration_s=0.5, ik_tol_m=1.0e-3),
        hypr_config=SpaceAGORA.RobotArmHYPRConfig(
            n_waypoints=0,
            n_particles=4,
            n_iters=1,
            n_samples=8,
            retime_enable=true,
            retime_max_joint_velocity_rad_s=Inf,
            retime_max_joint_acceleration_rad_s2=Inf,
            retime_max_base_force_n=0.05,
            retime_max_base_torque_nm=0.25,
            retime_base_wrench_iters=4,
            retime_cloth_dt_s=0.05,
            early_stopping_enable=false,
        ),
        obstacles=SpaceAGORA.RobotArmSphereObstacle[],
        rng=MersenneTwister(45),
    )
    @test cloth_retimed.plan.t_ref_s[end] > 0.5
    @test cloth_retimed.components.retime_base_force_ratio <= 1.0
    @test isfinite(cloth_retimed.components.retime_base_torque_ratio)

    dispatch_plan = SpaceAGORA.plan_robot_arm_motion(
        model,
        base,
        q0,
        target;
        config=SpaceAGORA.RobotArmPlannerConfig(dt_s=0.1, duration_s=0.5, ik_tol_m=1.0e-3),
        planner=:hypr,
        hypr_config=SpaceAGORA.RobotArmHYPRConfig(n_waypoints=1, n_particles=8, n_iters=2, n_samples=8, early_stopping_enable=false),
        obstacles=[obstacle],
        rng=MersenneTwister(43),
    )
    @test dispatch_plan.planner == :hypr
end

@testset "Closest surface target" begin
    points = [
        0.0 2.0 3.0;
        0.0 0.0 0.0;
        0.0 0.0 0.0
    ]
    target = SpaceAGORA.closest_surface_target(points, [2.1, 1.0, 0.0]; standoff_m=0.5)
    @test target.index == 2
    @test target.surface_point ≈ [2.0, 0.0, 0.0]
    @test norm(target.target - target.surface_point) ≈ 0.5
end

@testset "Compliant multibody Cloth dynamics" begin
    body = SpaceAGORA.CompliantBody(
        :body,
        2.0,
        SMatrix{3, 3, Float64}(1.0, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0, 1.4),
    )
    model = SpaceAGORA.CompliantMultibodyModel(
        [body],
        SpaceAGORA.CompliantJoint[],
        SVector{3, Float64}(0.0, 0.0, 0.0),
        SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
    )
    x0 = SpaceAGORA.compliant_state_vector(
        [SVector{3, Float64}(0.0, 0.0, 0.0)],
        [SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)];
        velocities=[SVector{3, Float64}(1.0, 0.0, 0.0)],
    )
    dx = SpaceAGORA.compliant_multibody_dynamics(model, x0)
    @test dx[1:3] ≈ [1.0, 0.0, 0.0]
    @test dx[8:10] ≈ [0.0, 0.0, 0.0]

    joint = SpaceAGORA.CompliantJoint(
        :anchor,
        0,
        1,
        SVector{3, Float64}(0.0, 0.0, 0.0),
        SVector{3, Float64}(0.0, 0.0, 0.0),
        SMatrix{3, 3, Float64}(100.0I),
        SMatrix{3, 3, Float64}(1.0I),
        SMatrix{3, 3, Float64}(0.0I),
        SMatrix{3, 3, Float64}(0.0I),
        SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
    )
    anchored = SpaceAGORA.CompliantMultibodyModel(
        [body],
        [joint],
        SVector{3, Float64}(0.0, 0.0, 0.0),
        SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
    )
    x_offset = SpaceAGORA.compliant_state_vector(
        [SVector{3, Float64}(0.1, 0.0, 0.0)],
        [SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)],
    )
    dx_offset = SpaceAGORA.compliant_multibody_dynamics(anchored, x_offset)
    @test dx_offset[8] < 0.0
end

@testset "Compliant topology builders and actuator loads" begin
    panel = SpaceAGORA.build_rectangular_compliant_grid(
        2,
        2;
        spacing_m=(0.5, 0.25),
        tile_size_m=(0.45, 0.20),
        mass_kg=0.2,
        thickness_m=0.005,
        anchor_index=1,
        k_translation_n_m=100.0,
        c_translation_n_s_m=1.0,
        k_rotation_n_m_rad=10.0,
        c_rotation_n_m_s_rad=0.1,
    )
    @test length(panel.model.bodies) == 4
    @test length(panel.model.joints) == 5
    @test length(panel.initial_state) == 52
    @test panel.model.joints[end].parent == 0

    actuator = SpaceAGORA.CompliantJointActuator(
        :limited_motor,
        1;
        torque_limit_n_m=0.05,
        kp_n_m_rad=10.0,
        kd_n_m_s_rad=0.0,
    )
    x = copy(panel.initial_state)
    qz = SVector{4, Float64}(0.0, 0.0, sin(0.5), cos(0.5))
    x[17:20] .= qz
    loads = SpaceAGORA.compliant_joint_loads(panel.model, x; joint_actuators=[actuator])
    @test length(loads) == length(panel.model.joints)
    @test norm(loads[1].actuator_torque_child_body) <= sqrt(3) * 0.05 + 1.0e-12
    @test norm(loads[1].compliance_torque_child_body) > 0.0
end

@testset "Cloth robot arm dynamics realization" begin
    model = SpaceAGORA.default_cloth_arm_model()
    base = SpaceAGORA.ClothArmBasePose([0.0, 0.0, 0.0])
    q0 = [0.0, 0.0, 0.0]
    target = SpaceAGORA.cloth_fk(model, base, [0.2, -0.2, 0.15]).end_effector_position
    cfg = SpaceAGORA.RobotArmPlannerConfig(dt_s=0.05, duration_s=0.15, ik_tol_m=1.0e-3)
    plan = SpaceAGORA.plan_robot_arm_motion(model, base, q0, target; config=cfg)
    per_joint_model = SpaceAGORA.cloth_robot_arm_multibody(
        plan;
        k_translation_n_m=[1.0e4, 2.0e4, 3.0e4],
        c_translation_n_s_m=[40.0, 50.0, 60.0],
        k_rotation_n_m_rad=[
            SMatrix{3, 3, Float64}(10.0I),
            SMatrix{3, 3, Float64}(20.0I),
            SMatrix{3, 3, Float64}(30.0I),
        ],
        c_rotation_n_m_s_rad=[(0.5, 0.6, 0.7), (0.8, 0.9, 1.0), (1.1, 1.2, 1.3)],
    )
    @test per_joint_model.joints[1].k_translation_n_m[1, 1] ≈ 1.0e4
    @test per_joint_model.joints[2].k_translation_n_m[1, 1] ≈ 2.0e4
    @test per_joint_model.joints[3].k_rotation_n_m_rad[1, 1] ≈ 30.0
    @test per_joint_model.joints[2].c_rotation_n_m_s_rad[2, 2] ≈ 0.9

    sim = SpaceAGORA.simulate_cloth_robot_arm_plan(
        plan;
        dt_s=0.05,
        duration_s=0.15,
        k_translation_n_m=[1.0e4, 1.1e4, 1.2e4],
        c_translation_n_s_m=[45.0, 50.0, 55.0],
        k_rotation_n_m_rad=[25.0, 30.0, 35.0],
        c_rotation_n_m_s_rad=[0.8, 1.0, 1.2],
        actuator_torque_limit_n_m=0.2,
        actuator_kp_n_m_rad=2.0,
        actuator_kd_n_m_s_rad=0.1,
    )
    @test length(sim.trajectory.t_s) == 4
    @test size(sim.end_effector_positions) == (3, 4)
    @test all(isfinite, sim.tracking_error_m)
    @test sim.tracking_error_m[1] < 1.0e-10
    @test size(sim.joint_compliance_torques_body) == (3, length(model.joints), 4)
    @test size(sim.joint_actuator_torques_body) == (3, length(model.joints), 4)
    @test sim.joint_total_torques_body ≈ sim.joint_compliance_torques_body .+ sim.joint_actuator_torques_body
end

@testset "Coupled cloth arm base reaction" begin
    model = SpaceAGORA.default_cloth_arm_model()
    base = SpaceAGORA.ClothArmBasePose([0.0, 0.0, 0.0])
    q0 = [0.0, 0.0, 0.0]
    target = SpaceAGORA.cloth_fk(model, base, [0.1, -0.1, 0.05]).end_effector_position
    cfg = SpaceAGORA.RobotArmPlannerConfig(dt_s=0.05, duration_s=0.1, ik_tol_m=1.0e-3)
    plan = SpaceAGORA.plan_robot_arm_motion(model, base, q0, target; config=cfg)
    shape = merge((
        pos=zeros(3),
        vel=zeros(3),
        mass=10.0,
        heat_loads=zeros(1),
        q=Float64[0.0, 0.0, 0.0, 1.0],
        ω=zeros(3),
    ), SpaceAGORA.coupled_cloth_robot_arm_state_shape(plan))
    sc = ComponentVector(shape)
    du = zero(sc)
    SpaceAGORA.initialize_coupled_cloth_robot_arm_state!(sc, plan)
    sc.arm_r[1, 1] += 0.01
    forces = MVector{3, Float64}(0.0, 0.0, 0.0)
    torques = MVector{3, Float64}(0.0, 0.0, 0.0)
    SpaceAGORA.assign_coupled_cloth_robot_arm_rhs!(
        du,
        sc,
        plan,
        0.0,
        forces,
        torques;
        k_translation_n_m=100.0,
        c_translation_n_s_m=0.0,
        k_rotation_n_m_rad=0.0,
        c_rotation_n_m_s_rad=0.0,
    )
    @test forces[1] > 0.0
    @test du.arm_v[1, 1] < 0.0
end
