using Test
using Dates
using DiffEqBase
using DiffEqCallbacks
using OrdinaryDiffEq
using Serialization
using StaticArrays
using ComponentArrays
using TOML
using LinearAlgebra
using Random

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using SpaceAGORA
const SimulationModel = SpaceAGORA.SimulationModel
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

end
