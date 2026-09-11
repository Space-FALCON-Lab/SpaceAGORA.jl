using Test
using Dates
using DiffEqBase
using DiffEqCallbacks
using OrdinaryDiffEq
using Quaternions
using Serialization
using StaticArrays
using ComponentArrays
using TOML
using LinearAlgebra
using Random

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel

const SM = SimulationModel
const RA = SimulationModel.RobotArmPlanning

"Assert every column of a joint matrix respects the model's joint limits."
function _probe_within_joint_limits(model, q_matrix; tol=1.0e-9)
    lo = [joint.lower_rad for joint in model.joints]
    hi = [joint.upper_rad for joint in model.joints]
    for k in axes(q_matrix, 2)
        all(q_matrix[:, k] .>= lo .- tol) || return false
        all(q_matrix[:, k] .<= hi .+ tol) || return false
    end
    return true
end

"Maximum consecutive joint-space step along a sampled path."
function _probe_max_step(samples)
    m = 0.0
    for k in 2:size(samples, 2)
        m = max(m, norm(samples[:, k] .- samples[:, k - 1]))
    end
    return m
end

@testset "coverage debt: robot arm HYPR probes" begin
    model = SM.default_cloth_arm_model()
    base = SM.ClothArmBasePose([0.0, 0.0, 0.0])
    planner_cfg = SM.RobotArmPlannerConfig(dt_s=0.1, duration_s=0.5, ik_tol_m=1.0e-3)
    q0 = [0.0, 0.3, -0.3]
    q_goal_truth = [0.6, -0.4, 0.4]
    target = SM.cloth_fk(model, base, q_goal_truth).end_effector_position
    q_goal_ik = SM.cloth_ik(
        model,
        base,
        SVector{3, Float64}(target);
        q_seed=q0,
        max_iters=planner_cfg.ik_max_iters,
        position_tol_m=planner_cfg.ik_tol_m,
        damping=planner_cfg.ik_damping,
    )
    no_obstacles = SM.RobotArmSphereObstacle[]

    @testset "config construction and validation" begin
        obs_tuple = SM.RobotArmSphereObstacle((0.1, 0.2, 0.3), 0.05)
        @test obs_tuple.center == SVector{3, Float64}(0.1, 0.2, 0.3)
        @test obs_tuple.radius_m == 0.05

        cfg_ok = SM.RobotArmHYPRConfig()
        @test RA._validate_robot_arm_hypr_config(cfg_ok) === cfg_ok
        cfg_poly = SM.RobotArmHYPRConfig(curve_type=:polyline)
        @test RA._validate_robot_arm_hypr_config(cfg_poly) === cfg_poly

        bad_cfg_kwargs = [
            (n_waypoints=-1,),
            (n_particles=0,),
            (n_iters=0,),
            (n_samples=1,),
            (curve_type=:zigzag,),
            (safe_distance_m=-0.1,),
            (rrt_warmstart_iters=-1,),
            (rrt_warmstart_step_size_rad=0.0,),
            (rrt_warmstart_goal_sample_rate=1.5,),
            (rrt_warmstart_collision_sample_ds_rad=0.0,),
            (rrt_warmstart_connect_max_steps=0,),
            (rrt_warmstart_shortcut_iters=-1,),
            (rrt_warmstart_runtime_limit_s=-1.0,),
            (retime_max_joint_velocity_rad_s=0.0,),
            (retime_max_joint_acceleration_rad_s2=0.0,),
            (retime_reaction_time_s=-1.0,),
            (retime_reaction_gain=-1.0,),
            (retime_min_scale=0.5,),
            (retime_min_scale=2.0, retime_max_scale=1.5),
            (retime_max_base_force_n=0.0,),
            (retime_max_base_torque_nm=0.0,),
            (retime_base_wrench_model_gain=0.0,),
            (retime_base_wrench_margin=0.9,),
            (retime_base_wrench_iters=-1,),
            (refinement_rounds=-1,),
            (refinement_step_fraction=0.0,),
            (refinement_shrink=1.0,),
            (refinement_min_abs_cost_improvement=-1.0,),
            (refinement_min_rel_cost_improvement=-1.0,),
            (retime_cloth_dt_s=0.0,),
            (retime_cloth_k_translation_n_m=-1.0,),
            (retime_cloth_c_translation_n_s_m=-1.0,),
            (retime_cloth_k_rotation_n_m_rad=-1.0,),
            (retime_cloth_c_rotation_n_m_s_rad=-1.0,),
        ]
        for kw in bad_cfg_kwargs
            @test_throws ArgumentError RA._validate_robot_arm_hypr_config(
                SM.RobotArmHYPRConfig(; kw...),
            )
        end
    end

    @testset "clearance stats" begin
        points = hcat(q0, q_goal_truth)
        samples = SM.robot_arm_sample_hypr_path(points, 16; curve_type=:polyline)

        empty_stats = SM.robot_arm_clearance_stats_from_samples(
            model, base, samples, no_obstacles, 0.01,
        )
        @test empty_stats.min_clearance == Inf
        @test empty_stats.violation_count == 0
        @test empty_stats.violation_fraction == 0.0
        @test empty_stats.clearance_penalty == 0.0

        far_obs = [SM.RobotArmSphereObstacle([50.0, 50.0, 50.0], 0.1)]
        far_stats = SM.robot_arm_clearance_stats_from_samples(
            model, base, samples, far_obs, 0.01,
        )
        @test far_stats.min_clearance > 0.0
        @test isfinite(far_stats.min_clearance)
        @test far_stats.violation_count == 0
        @test far_stats.clearance_penalty == 0.0

        # Sphere swallowing the arm base guarantees negative clearance.
        near_obs = [SM.RobotArmSphereObstacle([0.0, 0.0, 0.0], 0.06)]
        near_stats = SM.robot_arm_clearance_stats_from_samples(
            model, base, samples, near_obs, 0.005,
        )
        @test near_stats.min_clearance < 0.0
        @test near_stats.violation_count > 0
        @test near_stats.clearance_penalty > 0.0
        n_checks = size(samples, 2) * length(model.links) * length(near_obs)
        @test near_stats.violation_fraction ≈ near_stats.violation_count / n_checks
        @test 0.0 < near_stats.violation_fraction <= 1.0
    end

    @testset "cost components and segment distance" begin
        cfg = SM.RobotArmHYPRConfig(n_samples=12, curve_type=:polyline)
        points = hcat(q0, q_goal_truth)
        comps = SM.robot_arm_hypr_path_cost_components(points, model, base, no_obstacles, cfg)
        @test comps.total ≈ cfg.w_len * comps.J_len_norm^2 +
            cfg.w_smooth * comps.J_smooth + cfg.w_obs * comps.J_obs
        @test comps.J_len ≈ norm(q_goal_truth .- q0) atol=1.0e-9
        @test comps.J_len_norm ≈ 1.0 atol=1.0e-9
        @test comps.J_obs == 0.0
        @test comps.min_clearance == Inf
        @test comps.len_ref ≈ norm(q_goal_truth .- q0)

        cut = SM.robot_arm_hypr_path_cost_components(
            points, model, base, no_obstacles, cfg; cost_cutoff=-1.0,
        )
        @test cut.total == Inf
        @test cut.J_len ≈ comps.J_len
        @test cut.violation_count == comps.violation_count

        p = SVector{3, Float64}(1.0, 0.0, 0.0)
        a = SVector{3, Float64}(0.0, 0.0, 0.0)
        @test RA._robot_arm_segment_distance(p, a, a) ≈ 1.0
        b = SVector{3, Float64}(2.0, 0.0, 0.0)
        @test RA._robot_arm_segment_distance(p, a, b) ≈ 0.0 atol=1.0e-12
        c = SVector{3, Float64}(0.0, 1.0, 0.0)
        @test RA._robot_arm_segment_distance(p, a, c) ≈ 1.0

        better = RA._robot_arm_hypr_refinement_better
        cfg_ref = SM.RobotArmHYPRConfig()
        lower_obs = (J_obs=0.0, min_clearance=1.0, total=10.0)
        higher_obs = (J_obs=1.0, min_clearance=0.0, total=1.0)
        @test better(lower_obs, higher_obs, cfg_ref) === true
        @test better(higher_obs, lower_obs, cfg_ref) === false
        clr_a = (J_obs=1.0, min_clearance=0.5, total=3.0)
        clr_b = (J_obs=1.0, min_clearance=0.1, total=3.0)
        @test better(clr_a, clr_b, cfg_ref) === true
        cost_a = (J_obs=0.0, min_clearance=1.0, total=1.0)
        cost_b = (J_obs=0.0, min_clearance=1.0, total=2.0)
        @test better(cost_a, cost_b, cfg_ref) === true
        @test better(cost_b, cost_b, cfg_ref) === false

        # Post-refinement early return when refinement is disabled.
        cfg_noref = SM.RobotArmHYPRConfig(refinement_enable=false, n_samples=10)
        pts4 = hcat(q0, 0.7 .* q0 .+ 0.3 .* q_goal_truth, 0.3 .* q0 .+ 0.7 .* q_goal_truth, q_goal_truth)
        refined, refined_comps, improved, rounds = RA._robot_arm_hypr_post_refine_points(
            pts4, model, base, no_obstacles, cfg_noref,
        )
        @test refined == pts4
        @test improved === false
        @test rounds == 0
        @test refined_comps.total ≈ SM.robot_arm_hypr_path_cost_components(
            pts4, model, base, no_obstacles, cfg_noref,
        ).total

        # A strongly kinked midpoint gives refinement room to actually improve.
        cfg_refine = SM.RobotArmHYPRConfig(
            n_samples=12, curve_type=:polyline, refinement_rounds=4,
        )
        kink = 0.5 .* (q0 .+ q_goal_truth) .+ [1.5, -1.5, 1.5]
        pts_kinked = hcat(q0, kink, q_goal_truth)
        base_comps = SM.robot_arm_hypr_path_cost_components(
            pts_kinked, model, base, no_obstacles, cfg_refine,
        )
        better_pts, better_comps, did_improve, ran_rounds = RA._robot_arm_hypr_post_refine_points(
            pts_kinked, model, base, no_obstacles, cfg_refine,
        )
        @test did_improve === true
        @test ran_rounds >= 1
        @test better_comps.total < base_comps.total
        @test better_pts[:, 1] == q0 && better_pts[:, end] == q_goal_truth
        @test _probe_within_joint_limits(model, better_pts)
    end

    @testset "planner fast path (n_waypoints == 0)" begin
        cfg0 = SM.RobotArmHYPRConfig(n_waypoints=0, n_samples=10)
        result = SM.plan_robot_arm_motion_hypr(
            model, base, q0, target;
            planner_config=planner_cfg,
            hypr_config=cfg0,
            obstacles=no_obstacles,
            rng=MersenneTwister(11),
        )
        @test result.plan.planner == :hypr
        @test result.control_points ≈ hcat(q0, q_goal_ik)
        @test size(result.sampled_path) == (3, cfg0.n_samples)
        @test result.sampled_path[:, 1] ≈ q0
        @test result.sampled_path[:, end] ≈ q_goal_ik atol=1.0e-9
        @test result.cost_history == [result.cost]
        @test result.cost == result.components.total
        @test result.early_stopped === false
        @test result.early_stop_iter == 0
        @test result.cull_replacements == 0
        @test result.components.refinement_rounds == 0
        @test result.components.rrt_warmstart_enabled === false
        @test result.components.rrt_warmstart_attempted === false
        @test result.components.rrt_warmstart_n_points == 0
        @test result.components.retime_enabled === false
        @test result.components.retime_nominal_duration_s == planner_cfg.duration_s
        @test result.components.retime_scale_max ≈ 1.0
        @test result.plan.t_ref_s[1] == 0.0
        @test issorted(result.plan.t_ref_s) && allunique(result.plan.t_ref_s)
        @test result.plan.t_ref_s[end] ≈ planner_cfg.duration_s
        @test _probe_within_joint_limits(model, result.plan.q_ref)
        @test all(isfinite, result.plan.dq_ref)
        @test all(isfinite, result.plan.ddq_ref)
        @test result.plan.final_error_m < 5.0e-3
        @test result.plan.q_start ≈ q0
        @test result.plan.q_goal ≈ q_goal_ik
    end

    @testset "planner full PSO path, determinism, and culling" begin
        cfg = SM.RobotArmHYPRConfig(
            n_waypoints=2,
            n_particles=10,
            n_iters=6,
            n_samples=12,
            cull_enable=true,
            cull_start_iter=2,
            early_stopping_enable=false,
            refinement_rounds=2,
        )
        run_once(seed) = SM.plan_robot_arm_motion_hypr(
            model, base, q0, target;
            planner_config=planner_cfg,
            hypr_config=cfg,
            obstacles=no_obstacles,
            rng=MersenneTwister(seed),
        )
        result = run_once(21)
        @test size(result.control_points) == (3, cfg.n_waypoints + 2)
        @test result.control_points[:, 1] ≈ q0
        @test result.control_points[:, end] ≈ q_goal_ik
        @test _probe_within_joint_limits(model, result.control_points)
        @test _probe_within_joint_limits(model, result.sampled_path)
        @test size(result.sampled_path) == (3, cfg.n_samples)
        @test all(isfinite, result.sampled_path)
        @test _probe_max_step(result.sampled_path) < 1.0
        @test length(result.cost_history) == cfg.n_iters
        @test all(diff(result.cost_history) .<= 1.0e-12)
        @test result.cost == result.components.total
        @test result.cost <= result.cost_history[end] + 1.0e-9
        @test result.components.total ≈ cfg.w_len * result.components.J_len_norm^2 +
            cfg.w_smooth * result.components.J_smooth + cfg.w_obs * result.components.J_obs
        @test result.cull_replacements > 0
        @test result.components.refinement_enabled === true
        @test result.components.refinement_rounds >= 0
        @test result.plan.final_error_m < 5.0e-3
        @test issorted(result.plan.t_ref_s)

        repeat_result = run_once(21)
        @test repeat_result.control_points == result.control_points
        @test repeat_result.cost == result.cost

        # Schedule disabled variant exercises the static-weights branch.
        cfg_static = SM.RobotArmHYPRConfig(
            n_waypoints=1,
            n_particles=6,
            n_iters=2,
            n_samples=10,
            schedule_enable=false,
            cull_enable=false,
            early_stopping_enable=false,
            refinement_enable=false,
        )
        static_result = SM.plan_robot_arm_motion_hypr(
            model, base, q0, target;
            planner_config=planner_cfg,
            hypr_config=cfg_static,
            obstacles=no_obstacles,
            rng=MersenneTwister(5),
        )
        @test static_result.cull_replacements == 0
        @test static_result.components.refinement_improved === false
        @test isfinite(static_result.cost)

        # The q_start length check must fire before IK runs, with the planner's own message.
        @test_throws ArgumentError("q_start does not match the arm model.") SM.plan_robot_arm_motion_hypr(
            model, base, [0.0, 0.0], target;
            planner_config=planner_cfg,
            hypr_config=cfg_static,
            obstacles=no_obstacles,
            rng=MersenneTwister(5),
        )
    end

    @testset "planner early stopping" begin
        cfg_es = SM.RobotArmHYPRConfig(
            n_waypoints=1,
            n_particles=6,
            n_iters=40,
            n_samples=10,
            cull_enable=false,
            refinement_enable=false,
            early_stopping_enable=true,
            early_stopping_min_iters=2,
            early_stopping_patience=2,
            early_stopping_require_feasible=true,
        )
        result = SM.plan_robot_arm_motion_hypr(
            model, base, q0, target;
            planner_config=planner_cfg,
            hypr_config=cfg_es,
            obstacles=no_obstacles,
            rng=MersenneTwister(31),
        )
        @test result.early_stopped === true
        @test 0 < result.early_stop_iter < cfg_es.n_iters
        @test length(result.cost_history) == result.early_stop_iter

        # A harder swarm keeps improving past min_iters before stalling, so the
        # stale-counter reset path is exercised as well.
        cfg_es_hard = SM.RobotArmHYPRConfig(
            n_waypoints=2,
            n_particles=12,
            n_iters=60,
            n_samples=12,
            cull_enable=false,
            refinement_enable=false,
            early_stopping_enable=true,
            early_stopping_min_iters=1,
            early_stopping_patience=5,
        )
        q_mid_es = 0.5 .* (q0 .+ q_goal_ik)
        blocker_es = SM.RobotArmSphereObstacle(
            SM.cloth_fk(model, base, q_mid_es).end_effector_position, 0.015,
        )
        hard = SM.plan_robot_arm_motion_hypr(
            model, base, q0, target;
            planner_config=planner_cfg,
            hypr_config=cfg_es_hard,
            obstacles=[blocker_es],
            rng=MersenneTwister(32),
        )
        @test length(hard.cost_history) <= cfg_es_hard.n_iters
        @test all(diff(hard.cost_history) .<= 1.0e-12)
        if hard.early_stopped
            @test hard.early_stop_iter == length(hard.cost_history)
        end
        # The swarm improved after the early-stopping window opened.
        @test hard.cost_history[end] < hard.cost_history[1]
    end

    @testset "planner with obstacles" begin
        q_mid = 0.5 .* (q0 .+ q_goal_ik)
        blocker_center = SM.cloth_fk(model, base, q_mid).end_effector_position
        blocker = SM.RobotArmSphereObstacle(blocker_center, 0.02)
        cfg_obs = SM.RobotArmHYPRConfig(
            n_waypoints=2,
            n_particles=12,
            n_iters=8,
            n_samples=14,
            safe_distance_m=0.004,
            early_stopping_enable=false,
            refinement_rounds=2,
        )
        result = SM.plan_robot_arm_motion_hypr(
            model, base, q0, target;
            planner_config=planner_cfg,
            hypr_config=cfg_obs,
            obstacles=[blocker],
            rng=MersenneTwister(77),
        )
        @test isfinite(result.components.min_clearance)
        recomputed = SM.robot_arm_clearance_stats_from_samples(
            model,
            base,
            SM.robot_arm_sample_hypr_path(result.control_points, cfg_obs.n_samples;
                curve_type=cfg_obs.curve_type),
            [blocker],
            cfg_obs.safe_distance_m,
        )
        @test recomputed.min_clearance ≈ result.components.min_clearance
        @test recomputed.violation_count == result.components.violation_count
        @test recomputed.clearance_penalty ≈ result.components.clearance_penalty
        @test result.components.J_obs ≈ recomputed.violation_count +
            recomputed.clearance_penalty / max(cfg_obs.safe_distance_m^2, 1.0e-8)
        @test _probe_within_joint_limits(model, result.control_points)
    end

    @testset "rrt warmstart" begin
        # Direct-path branch: no obstacles means the straight segment is safe.
        cfg_direct = SM.RobotArmHYPRConfig(
            n_waypoints=1,
            n_particles=6,
            n_iters=2,
            n_samples=10,
            rrt_warmstart_enable=true,
            rrt_warmstart_iters=8,
            early_stopping_enable=false,
            refinement_enable=false,
        )
        direct = SM.plan_robot_arm_motion_hypr(
            model, base, q0, target;
            planner_config=planner_cfg,
            hypr_config=cfg_direct,
            obstacles=no_obstacles,
            rng=MersenneTwister(46),
        )
        @test direct.components.rrt_warmstart_enabled === true
        @test direct.components.rrt_warmstart_attempted === true
        @test direct.components.rrt_warmstart_path_found === true
        @test direct.components.rrt_warmstart_iterations == 0
        @test direct.components.rrt_warmstart_n_points == 2
        @test isfinite(direct.components.rrt_warmstart_cost)
        @test direct.components.rrt_warmstart_cost == direct.components.rrt_warmstart_raw_cost
        @test direct.components.rrt_warmstart.path_found === true

        # Blocked direct path: the trees must grow and connect around the sphere.
        q_mid = 0.5 .* (q0 .+ q_goal_ik)
        blocker = SM.RobotArmSphereObstacle(
            SM.cloth_fk(model, base, q_mid).end_effector_position, 0.02,
        )
        cfg_rrt = SM.RobotArmHYPRConfig(
            n_waypoints=2,
            n_particles=6,
            n_iters=2,
            n_samples=12,
            safe_distance_m=0.004,
            rrt_warmstart_enable=true,
            rrt_warmstart_iters=600,
            rrt_warmstart_step_size_rad=0.35,
            rrt_warmstart_goal_sample_rate=0.2,
            rrt_warmstart_collision_sample_ds_rad=0.15,
            rrt_warmstart_shortcut_iters=15,
            early_stopping_enable=false,
            refinement_enable=false,
        )
        warm_path, warm_diag = RA._robot_arm_rrt_connect_warmstart_path(
            q0, q_goal_ik, model, base, [blocker], cfg_rrt, MersenneTwister(101),
        )
        @test warm_diag.attempted === true
        @test warm_diag.iterations >= 1
        @test warm_diag.path_found === true
        @test warm_path !== nothing
        @test size(warm_path, 1) == 3
        @test warm_diag.n_points == size(warm_path, 2) >= 2
        @test warm_path[:, 1] ≈ q0
        @test warm_path[:, end] ≈ q_goal_ik
        @test isfinite(warm_diag.cost) && isfinite(warm_diag.raw_cost)
        # The warm-started path must itself be collision free at the checked tolerance.
        warm_stats = SM.robot_arm_clearance_stats_from_samples(
            model,
            base,
            SM.robot_arm_sample_hypr_path(warm_path, 40; curve_type=:polyline),
            [blocker],
            cfg_rrt.safe_distance_m,
        )
        @test warm_stats.min_clearance > 0.0

        warmstarted = SM.plan_robot_arm_motion_hypr(
            model, base, q0, target;
            planner_config=planner_cfg,
            hypr_config=cfg_rrt,
            obstacles=[blocker],
            rng=MersenneTwister(101),
        )
        @test warmstarted.components.rrt_warmstart_path_found === true
        @test warmstarted.control_points[:, 1] ≈ q0
        @test warmstarted.control_points[:, end] ≈ warmstarted.plan.q_goal

        # Runtime-limited warm start bails out after the first iteration.
        cfg_budget = SM.RobotArmHYPRConfig(
            rrt_warmstart_enable=true,
            rrt_warmstart_iters=500,
            rrt_warmstart_step_size_rad=0.35,
            rrt_warmstart_goal_sample_rate=0.0,
            rrt_warmstart_collision_sample_ds_rad=0.15,
            rrt_warmstart_runtime_limit_s=0.0,
            rrt_warmstart_connect_max_steps=1,
            safe_distance_m=0.004,
        )
        budget_path, budget_diag = RA._robot_arm_rrt_connect_warmstart_path(
            q0, q_goal_ik, model, base, [blocker], cfg_budget, MersenneTwister(102),
        )
        @test budget_diag.attempted === true
        @test budget_diag.iterations <= 1
        @test budget_diag.path_found === false
        @test budget_path === nothing
        @test budget_diag.n_points == 2
        @test budget_diag.cost == budget_diag.raw_cost

        # Disabled warm start returns empty diagnostics.
        cfg_off = SM.RobotArmHYPRConfig(rrt_warmstart_enable=false)
        off_path, off_diag = RA._robot_arm_rrt_connect_warmstart_path(
            q0, q_goal_ik, model, base, [blocker], cfg_off, MersenneTwister(1),
        )
        @test off_path === nothing
        @test off_diag.enabled === false
        @test off_diag.attempted === false
        @test off_diag.cost == Inf

        # Segment sampling covers coarse and fine discretizations.
        seg = RA._robot_arm_rrt_segment_samples([0.0, 0.0, 0.0], [0.1, 0.0, 0.0], 1.0)
        @test size(seg) == (3, 2)
        seg_fine = RA._robot_arm_rrt_segment_samples([0.0, 0.0, 0.0], [0.1, 0.0, 0.0], 0.01)
        @test size(seg_fine, 2) == 11
        @test seg_fine[1, 1] == 0.0 && seg_fine[1, end] ≈ 0.1

        # Tree primitives.
        tree = RA.RobotArmRRTConnectTree([0.0, 0.0, 0.0])
        @test tree.nodes == [[0.0, 0.0, 0.0]]
        @test RA._robot_arm_rrt_nearest_index(tree, [1.0, 0.0, 0.0]) == 1
        status, idx = RA._robot_arm_rrt_extend!(
            tree, [1.0, 0.0, 0.0], model, base, no_obstacles,
            SM.RobotArmHYPRConfig(rrt_warmstart_step_size_rad=0.25),
        )
        @test status == :advanced
        @test idx == 2
        @test norm(tree.nodes[2] .- [0.25, 0.0, 0.0]) < 1.0e-9
        @test tree.costs[2] ≈ 0.25

        # Shortcutting a zigzag polyline with no obstacles shrinks or preserves length.
        zig = hcat(
            [0.0, 0.0, 0.0],
            [0.2, 0.4, 0.0],
            [0.4, -0.4, 0.0],
            [0.6, 0.0, 0.0],
        )
        cut = RA._robot_arm_rrt_shortcut_path(
            zig, model, base, no_obstacles, cfg_rrt, MersenneTwister(3),
        )
        @test RA._robot_arm_path_length(cut) <= RA._robot_arm_path_length(zig) + 1.0e-9
        @test cut[:, 1] == zig[:, 1] && cut[:, end] == zig[:, end]
        two = hcat([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        @test RA._robot_arm_rrt_shortcut_path(
            two, model, base, no_obstacles, cfg_rrt, MersenneTwister(3),
        ) == two

        # Polyline resampling branches.
        @test_throws ArgumentError RA._robot_arm_resample_polyline_points(zig, 1)
        same = RA._robot_arm_resample_polyline_points(zig, 4)
        @test same == zig && same !== zig
        expanded = RA._robot_arm_resample_polyline_points(zig, 7)
        @test size(expanded) == (3, 7)
        @test expanded[:, 1] == zig[:, 1] && expanded[:, end] == zig[:, end]
        single = reshape([1.0, 2.0, 3.0], 3, 1)
        widened = RA._robot_arm_resample_polyline_points(single, 3)
        @test widened == hcat(single, zeros(3), single)
        degenerate = hcat([1.0, 1.0, 1.0], [1.0, 1.0, 1.0])
        deg = RA._robot_arm_resample_polyline_points(degenerate, 4)
        @test all(deg .== 1.0)
    end

    @testset "retiming and base wrench" begin
        # Velocity plus acceleration limited retime stretches the schedule.
        cfg_vel = SM.RobotArmHYPRConfig(
            n_waypoints=0,
            n_samples=8,
            retime_enable=true,
            retime_max_joint_velocity_rad_s=0.2,
            retime_max_joint_acceleration_rad_s2=0.8,
        )
        retimed = SM.plan_robot_arm_motion_hypr(
            model, base, q0, target;
            planner_config=planner_cfg,
            hypr_config=cfg_vel,
            obstacles=no_obstacles,
            rng=MersenneTwister(44),
        )
        @test retimed.plan.t_ref_s[end] > planner_cfg.duration_s
        @test issorted(retimed.plan.t_ref_s) && allunique(retimed.plan.t_ref_s)
        @test retimed.components.retime_enabled === true
        @test retimed.components.retime_duration_s == retimed.plan.t_ref_s[end]
        @test retimed.components.retime_scale_max ≈
            retimed.plan.t_ref_s[end] / planner_cfg.duration_s
        @test maximum(abs.(retimed.plan.dq_ref)) <=
            cfg_vel.retime_max_joint_velocity_rad_s + 1.0e-8

        # Reaction-time tapered retime (no wrench limits) stays within scale bounds.
        cfg_react = SM.RobotArmHYPRConfig(
            n_waypoints=0,
            n_samples=8,
            retime_enable=true,
            retime_max_joint_velocity_rad_s=0.2,
            retime_reaction_time_s=0.5,
            retime_reaction_gain=2.0,
            retime_max_scale=50.0,
        )
        reacted = SM.plan_robot_arm_motion_hypr(
            model, base, q0, target;
            planner_config=planner_cfg,
            hypr_config=cfg_react,
            obstacles=no_obstacles,
            rng=MersenneTwister(44),
        )
        @test reacted.plan.t_ref_s[end] >= retimed.plan.t_ref_s[end] - 1.0e-9
        @test reacted.plan.t_ref_s[end] <=
            cfg_react.retime_max_scale * planner_cfg.duration_s + 1.0e-9

        @test RA._robot_arm_hypr_reaction_scale(SM.RobotArmHYPRConfig(), 3.0, 2, 5) == 1.0
        cfg_scale = SM.RobotArmHYPRConfig(retime_reaction_time_s=1.0, retime_reaction_gain=1.0)
        @test RA._robot_arm_hypr_reaction_scale(cfg_scale, 1.0, 1, 3) == 1.0
        @test RA._robot_arm_hypr_reaction_scale(cfg_scale, 1.0, 2, 3) ≈ 2.0

        t_from_scales = RA._robot_arm_hypr_reference_times_from_scales(0.1, [1.0, 2.0, 1.0])
        @test t_from_scales ≈ [0.0, 0.1, 0.3, 0.4]

        # Rigid-model wrench ratios: disabled limits short-circuit to zero.
        q_ref = SM.robot_arm_sample_hypr_path(hcat(q0, q_goal_ik), 6; curve_type=:polyline)
        t_ref = collect(0.0:0.1:0.5)
        cfg_plain = SM.RobotArmHYPRConfig()
        zero_ratios = RA._robot_arm_hypr_rigid_base_wrench_ratios(model, base, q_ref, t_ref, cfg_plain)
        @test zero_ratios.force_ratio == 0.0
        @test zero_ratios.torque_ratio == 0.0

        cfg_rigid = SM.RobotArmHYPRConfig(
            retime_cloth_physics_enable=false,
            retime_max_base_force_n=0.5,
            retime_max_base_torque_nm=0.5,
        )
        rigid_ratios = RA._robot_arm_hypr_rigid_base_wrench_ratios(model, base, q_ref, t_ref, cfg_rigid)
        @test rigid_ratios.force_ratio > 0.0
        @test rigid_ratios.torque_ratio > 0.0
        @test length(rigid_ratios.node_ratio) == size(q_ref, 2)
        @test maximum(rigid_ratios.node_ratio) >= max(rigid_ratios.force_ratio, rigid_ratios.torque_ratio) - 1.0e-12

        single_ratios = RA._robot_arm_hypr_rigid_base_wrench_ratios(
            model, base, reshape(q0, 3, 1), [0.0], cfg_rigid,
        )
        @test single_ratios.force_ratio == 0.0
        @test single_ratios.node_ratio == [0.0]

        # Dispatch: cloth physics disabled routes to the rigid estimator.
        dispatch_ratios = RA._robot_arm_hypr_base_wrench_ratios(model, base, q_ref, t_ref, cfg_rigid)
        @test dispatch_ratios.force_ratio ≈ rigid_ratios.force_ratio
        @test dispatch_ratios.torque_ratio ≈ rigid_ratios.torque_ratio

        # Rigid wrench-limited retime iterates the scale loop.
        cfg_wrench = SM.RobotArmHYPRConfig(
            n_waypoints=0,
            n_samples=8,
            retime_enable=true,
            retime_cloth_physics_enable=false,
            retime_max_base_force_n=0.02,
            retime_max_base_torque_nm=0.05,
            retime_base_wrench_iters=4,
        )
        wrenched = SM.plan_robot_arm_motion_hypr(
            model, base, q0, target;
            planner_config=planner_cfg,
            hypr_config=cfg_wrench,
            obstacles=no_obstacles,
            rng=MersenneTwister(45),
        )
        @test wrenched.plan.t_ref_s[end] > planner_cfg.duration_s
        @test wrenched.components.retime_base_force_ratio <= 1.0 + 1.0e-6
        @test wrenched.components.retime_base_force_ratio >= 0.0
        @test isfinite(wrenched.components.retime_base_torque_ratio)

        # Cloth-coupled wrench estimation drives the coupled simulation path.
        cfg_cloth = SM.RobotArmHYPRConfig(
            n_waypoints=0,
            n_samples=8,
            retime_enable=true,
            retime_max_base_force_n=0.05,
            retime_max_base_torque_nm=0.25,
            retime_base_wrench_iters=2,
            retime_cloth_dt_s=0.05,
        )
        cloth_retimed = SM.plan_robot_arm_motion_hypr(
            model, base, [0.0, -0.05, 0.05],
            SM.cloth_fk(model, base, [0.08, -0.12, 0.08]).end_effector_position;
            planner_config=planner_cfg,
            hypr_config=cfg_cloth,
            obstacles=no_obstacles,
            rng=MersenneTwister(45),
        )
        @test cloth_retimed.plan.t_ref_s[end] > planner_cfg.duration_s
        @test isfinite(cloth_retimed.components.retime_base_force_ratio)
        @test isfinite(cloth_retimed.components.retime_base_torque_ratio)
        @test issorted(cloth_retimed.plan.t_ref_s)

        # Direct retime-reference branches.
        nominal = RA._robot_arm_hypr_retime_reference(
            model, base, reshape(q0, 3, 1), planner_cfg, cfg_vel,
        )
        @test nominal[1] == 0.0 && nominal[end] ≈ planner_cfg.duration_s
        off = RA._robot_arm_hypr_retime_reference(
            model, base, q_ref, planner_cfg, SM.RobotArmHYPRConfig(retime_enable=false),
        )
        @test off[end] ≈ planner_cfg.duration_s
        @test_throws ArgumentError RA._robot_arm_hypr_retime_reference(
            model, base, q_ref, SM.RobotArmPlannerConfig(duration_s=0.0), cfg_vel,
        )
    end

    @testset "swarm helpers" begin
        seeded = RA._robot_arm_seed_control_points(q0, q_goal_ik, 3)
        @test size(seeded) == (3, 5)
        @test seeded[:, 1] == q0
        @test seeded[:, end] == q_goal_ik
        @test seeded[:, 3] ≈ 0.5 .* (q0 .+ q_goal_ik)

        flat = RA._robot_arm_flatten_internal_points(seeded)
        @test length(flat) == 9
        rebuilt = RA._robot_arm_control_points(q0, q_goal_ik, flat, 3)
        @test rebuilt ≈ seeded

        pts = hcat(q0, q_goal_ik)
        poly = SM.robot_arm_sample_hypr_path(pts, 5; curve_type=:polyline)
        @test size(poly) == (3, 5)
        @test poly[:, 1] ≈ q0 && poly[:, end] ≈ q_goal_ik
        bez = SM.robot_arm_sample_hypr_path(seeded, 9; curve_type=:bezier)
        @test size(bez) == (3, 9)
        @test bez[:, 1] ≈ q0 && bez[:, end] ≈ q_goal_ik

        @test RA._robot_arm_path_length(poly) ≈ norm(q_goal_ik .- q0)
        @test RA._robot_arm_path_smoothness(hcat(q0, q_goal_ik)) == 0.0
        @test RA._robot_arm_path_smoothness(poly) ≈ 0.0 atol=1.0e-12
        kinked = hcat([0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0])
        @test RA._robot_arm_path_smoothness(kinked) > 0.0

        cfg_default = SM.RobotArmHYPRConfig()
        @test RA._robot_arm_hypr_material_improvement(1.0, 2.0, cfg_default) === true
        @test RA._robot_arm_hypr_material_improvement(2.0, 1.0, cfg_default) === false
        @test RA._robot_arm_hypr_early_stopping_feasible((J_obs=5.0,), cfg_default) === true
        cfg_feas = SM.RobotArmHYPRConfig(early_stopping_require_feasible=true)
        @test RA._robot_arm_hypr_early_stopping_feasible((J_obs=0.0,), cfg_feas) === true
        @test RA._robot_arm_hypr_early_stopping_feasible((J_obs=1.0,), cfg_feas) === false

        weights_sched = RA._robot_arm_hypr_iteration_weights(SM.RobotArmHYPRConfig(n_iters=10), 5)
        @test all(isfinite, (weights_sched.w_inertia, weights_sched.c1, weights_sched.c2))
        weights_flat = RA._robot_arm_hypr_iteration_weights(
            SM.RobotArmHYPRConfig(schedule_enable=false), 3,
        )
        @test weights_flat.w_inertia == SM.RobotArmHYPRConfig().w_inertia

        # Culling: replaces the worst particles and resets their bests.
        cull_cfg = SM.RobotArmHYPRConfig(cull_enable=true, cull_start_iter=1, cull_fraction_max=0.5)
        positions = [0.1 0.2 0.3 0.4; 0.1 0.2 0.3 0.4; 0.1 0.2 0.3 0.4]
        velocities = ones(3, 4)
        pbest = copy(positions)
        pbest_cost = [1.0, 2.0, 3.0, 4.0]
        lo = fill(-1.0, 3)
        hi = fill(1.0, 3)
        gbest = zeros(3)
        n_replaced = RA._robot_arm_hypr_cull_swarm!(
            positions, velocities, pbest, pbest_cost, lo, hi, gbest, cull_cfg, 5, MersenneTwister(9),
        )
        @test n_replaced == 2
        @test count(isinf, pbest_cost) == 2
        @test pbest_cost[1] == 1.0 && pbest_cost[2] == 2.0
        @test all(velocities[:, 3] .== 0.0) && all(velocities[:, 4] .== 0.0)
        @test all(-1.0 .<= positions .<= 1.0)
        @test pbest[:, 3] == positions[:, 3]

        off_cfg = SM.RobotArmHYPRConfig(cull_enable=false)
        @test RA._robot_arm_hypr_cull_swarm!(
            positions, velocities, pbest, pbest_cost, lo, hi, gbest, off_cfg, 50, MersenneTwister(9),
        ) == 0
        early_cfg = SM.RobotArmHYPRConfig(cull_enable=true, cull_start_iter=100)
        @test RA._robot_arm_hypr_cull_swarm!(
            positions, velocities, pbest, pbest_cost, lo, hi, gbest, early_cfg, 5, MersenneTwister(9),
        ) == 0
        all_inf_cost = fill(Inf, 4)
        @test RA._robot_arm_hypr_cull_swarm!(
            positions, velocities, pbest, all_inf_cost, lo, hi, gbest, cull_cfg, 5, MersenneTwister(9),
        ) == 0
        one_particle = reshape([0.0, 0.0, 0.0], 3, 1)
        @test RA._robot_arm_hypr_cull_swarm!(
            one_particle, zeros(3, 1), copy(one_particle), [1.0], lo, hi, gbest, cull_cfg, 5,
            MersenneTwister(9),
        ) == 0
    end
end
