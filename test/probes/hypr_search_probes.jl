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
const GH = SimulationModel.GuidanceHooks

"Station geometry far away from the local RTN frame so paths near the origin are free."
function _hypr_far_geometry(; keepout_radius_m::Real=0.0)
    station = SM.RPOStationGeometry(
        reshape([100.0, 0.0, 0.0], 3, 1);
        keepout_radius_m=keepout_radius_m,
        name="hypr_far",
    )
    return SM.RPOReferenceGeometry(station; chaser=SM.RPOCubeSatGeometry(dims_m=(0.02, 0.02, 0.02)))
end

"Single-point station at the origin acting as one sphere obstacle of given keepout radius."
function _hypr_origin_geometry(; keepout_radius_m::Real=0.0)
    station = SM.RPOStationGeometry(
        reshape([0.0, 0.0, 0.0], 3, 1);
        keepout_radius_m=keepout_radius_m,
        name="hypr_origin",
    )
    return SM.RPOReferenceGeometry(station; chaser=SM.RPOCubeSatGeometry(dims_m=(0.02, 0.02, 0.02)))
end

"Tiny deterministic-ish PSO config with every optional subsystem off unless requested."
function _hypr_cfg(; kwargs...)
    return SM.RPOPSOConfig(;
        n_waypoints=1,
        n_particles=6,
        n_iters=3,
        sample_ds_m=0.25,
        refinement_sample_ds_m=0.1,
        curve_type=:polyline,
        search_margin_m=1.0,
        spread_scale=0.2,
        adaptive_enable=false,
        adaptive_sampling_enable=false,
        probe_enable=false,
        refinement_enable=false,
        refinement_rounds=2,
        refinement_waypoint_passes=4,
        rrt_warmstart_enable=false,
        cull_enable=false,
        reexplore_enable=false,
        stagnation_learning_enable=false,
        early_stopping_enable=false,
        schedule_enable=false,
        kwargs...,
    )
end

const _SV3 = SVector{3, Float64}

@testset "coverage debt hypr search probes" begin

    @testset "rrt primitives" begin
        far = _hypr_far_geometry()
        blocked = _hypr_origin_geometry(; keepout_radius_m=0.5)

        tree = GH.RPORRTConnectTree(_SV3(-1.5, 0.0, 0.0))
        @test length(tree.nodes) == 1
        @test tree.parents == [0]
        @test tree.costs == [0.0]

        @test GH.rpo_rrt_nearest_index(tree, _SV3(0.0, 0.0, 0.0)) == 1
        @test GH.rpo_rrt_near_indices(tree, _SV3(-1.5, 0.0, 0.0), -1.0) == [1]
        @test GH.rpo_rrt_near_indices(tree, _SV3(5.0, 0.0, 0.0), 0.5) == Int[]

        rng = MersenneTwister(3)
        lo = _SV3(-1.0, 0.0, 2.0)
        hi = _SV3(1.0, 2.0, 3.0)
        for _ in 1:5
            q = GH.rpo_rrt_random_state(rng, lo, hi)
            @test all(lo .<= q .<= hi)
        end
        @test GH.rpo_rrt_random_state(rng, lo, lo) == lo

        q, status = GH.rpo_rrt_steer(_SV3(0.0, 0.0, 0.0), _SV3(0.0, 0.0, 0.0), 0.5)
        @test status == :trapped
        q, status = GH.rpo_rrt_steer(_SV3(0.0, 0.0, 0.0), _SV3(0.2, 0.0, 0.0), 0.5)
        @test status == :reached && q ≈ _SV3(0.2, 0.0, 0.0)
        q, status = GH.rpo_rrt_steer(_SV3(0.0, 0.0, 0.0), _SV3(2.0, 0.0, 0.0), 0.5)
        @test status == :advanced && q ≈ _SV3(0.5, 0.0, 0.0)

        # Collision spacing: safety fraction and obstacle guard both bind.
        @test GH.rpo_rrt_collision_min_ds_m(0.2, far, 0.0, 0.5, 0.5) ≈ 0.005  # guard on 0.01 radius
        @test GH.rpo_rrt_collision_min_ds_m(0.2, blocked, 0.1, 0.5, 0.5) ≈ 0.05
        zero_station = SM.RPOStationGeometry(reshape([100.0, 0.0, 0.0], 3, 1); keepout_radius_m=0.0)
        zero_geom = SM.RPOReferenceGeometry(
            zero_station;
            chaser=SM.RPOCubeSatGeometry(_SV3(0.0, 0.0, 0.0)),
        )
        @test GH.rpo_rrt_collision_min_ds_m(0.2, zero_geom, 0.0, 0.5, 0.5) ≈ 0.2  # no inflation

        # Segment safety in kwarg form, uniform and adaptive sampling.
        @test GH.rpo_rrt_segment_is_safe(
            (-1.0, 2.0, 0.0), (1.0, 2.0, 0.0), blocked;
            safe_distance_m=0.05, sample_ds_m=0.05, adaptive_enable=false)
        @test !GH.rpo_rrt_segment_is_safe(
            (-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), blocked;
            safe_distance_m=0.05, sample_ds_m=0.05, adaptive_enable=false)
        @test GH.rpo_rrt_segment_is_safe(
            (-1.0, 2.0, 0.0), (1.0, 2.0, 0.0), blocked;
            safe_distance_m=0.05, sample_ds_m=0.05)
        @test !GH.rpo_rrt_segment_is_safe(
            (-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), blocked;
            safe_distance_m=0.05, sample_ds_m=0.05)

        settings = GH.RPORRTConnectSettings(; step_size_m=0.5, collision_sample_ds_m=0.05)
        star_settings = GH.RPORRTStarSettings(; step_size_m=0.5, collision_sample_ds_m=0.05)
        @test GH.rpo_rrt_segment_is_safe((-1.0, 2.0, 0.0), (1.0, 2.0, 0.0), blocked, settings; safe_distance_m=0.05)
        @test !GH.rpo_rrt_segment_is_safe((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), blocked, star_settings; safe_distance_m=0.05)

        # Extend: trapped by steer, advanced, trapped by collision.
        st, idx = GH.rpo_rrt_extend!(tree, _SV3(-1.5, 0.0, 0.0), blocked, settings; safe_distance_m=0.05)
        @test st == :trapped && idx == 1
        @test length(tree.nodes) == 1
        st, idx = GH.rpo_rrt_extend!(tree, _SV3(0.0, 0.0, 0.0), blocked, settings; safe_distance_m=0.05)
        @test st == :advanced && idx == 2
        @test tree.nodes[2] ≈ _SV3(-1.0, 0.0, 0.0)
        @test tree.costs[2] ≈ 0.5
        @test tree.parents[2] == 1
        st, idx = GH.rpo_rrt_extend!(tree, _SV3(0.0, 0.0, 0.0), blocked, settings; safe_distance_m=0.05)
        @test st == :trapped  # next step would enter the keepout sphere
        @test length(tree.nodes) == 2

        # Connect reaches a target in free space and respects the step budget.
        free_tree = GH.RPORRTConnectTree(_SV3(-1.5, 0.0, 0.0))
        st, idx = GH.rpo_rrt_connect!(free_tree, _SV3(-1.5, 1.3, 0.0), blocked, settings; safe_distance_m=0.05)
        @test st == :reached
        @test free_tree.nodes[idx] ≈ _SV3(-1.5, 1.3, 0.0)
        limited = GH.RPORRTConnectSettings(; step_size_m=0.5, collision_sample_ds_m=0.05, connect_max_steps=1)
        budget_tree = GH.RPORRTConnectTree(_SV3(-1.5, 0.0, 0.0))
        st, _ = GH.rpo_rrt_connect!(budget_tree, _SV3(-1.5, 5.0, 0.0), blocked, limited; safe_distance_m=0.05)
        @test st == :advanced
        @test length(budget_tree.nodes) == 2

        # Path reconstruction and joining through the shared HYPR helpers.
        chain = GH.RPORRTConnectTree(
            [_SV3(0.0, 0.0, 0.0), _SV3(1.0, 0.0, 0.0), _SV3(2.0, 0.0, 0.0)],
            [0, 1, 2],
            [0.0, 1.0, 2.0],
        )
        @test GH.rpo_rrt_tree_path(chain, 3) == [0.0 1.0 2.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        goal_chain = GH.RPORRTConnectTree(
            [_SV3(4.0, 0.0, 0.0), _SV3(3.0, 0.0, 0.0), _SV3(2.0, 0.0, 0.0)],
            [0, 1, 2],
            [0.0, 1.0, 2.0],
        )
        joined = GH.rpo_rrt_connect_join_paths(chain, 3, goal_chain, 3)
        @test joined == [0.0 1.0 2.0 3.0 4.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0]

        stale = GH.RPORRTConnectTree(
            [_SV3(0.0, 0.0, 0.0), _SV3(2.0, 0.0, 0.0), _SV3(2.0, 1.0, 0.0)],
            [0, 1, 2],
            [0.0, 99.0, 99.0],
        )
        GH.rpo_rrt_refresh_subtree_costs!(stale, 1)
        @test stale.costs ≈ [0.0, 2.0, 3.0]

        # Shortcutting keeps endpoints and never lengthens the path.
        zigzag = [0.0 1.0 2.0 3.0; 0.0 1.0 -1.0 0.0; 0.0 0.0 0.0 0.0]
        cut = GH.rpo_rrt_shortcut_path(zigzag, far, settings; safe_distance_m=0.0, rng=MersenneTwister(9))
        @test cut[:, 1] ≈ zigzag[:, 1]
        @test cut[:, end] ≈ zigzag[:, end]
        @test GH.rpo_path_length(cut) <= GH.rpo_path_length(zigzag) + 1.0e-9
        @test size(cut, 2) < size(zigzag, 2)
        @test GH.rpo_rrt_shortcut_path(zigzag[:, [1, 4]], far, settings; safe_distance_m=0.0, rng=MersenneTwister(9)) == zigzag[:, [1, 4]]
    end

    @testset "rrt star add node" begin
        far = _hypr_far_geometry()
        blocked = _hypr_origin_geometry(; keepout_radius_m=0.5)
        settings = GH.RPORRTStarSettings(; step_size_m=0.6, collision_sample_ds_m=0.05, neighbor_radius_m=1.5)

        # Rewire: the new node offers a cheaper route to an overpriced sibling.
        tree = GH.RPORRTConnectTree(
            [_SV3(0.0, 0.0, 0.0), _SV3(2.0, 0.0, 0.0)],
            [0, 1],
            [0.0, 10.0],
        )
        new_idx = GH.rpo_rrt_star_add_node!(tree, _SV3(1.0, 0.0, 0.0), far, settings; safe_distance_m=0.0)
        @test new_idx == 3
        @test tree.parents[3] == 1
        @test tree.costs[3] ≈ 1.0
        @test tree.parents[2] == 3  # rewired through the new node
        @test tree.costs[2] ≈ 2.0

        # Unsafe candidate is rejected outright.
        b_tree = GH.RPORRTConnectTree(_SV3(-1.5, 0.0, 0.0))
        @test GH.rpo_rrt_star_add_node!(b_tree, _SV3(0.0, 0.0, 0.0), blocked, settings; safe_distance_m=0.05) == 0
        @test length(b_tree.nodes) == 1

        # Empty neighborhood falls back to the nearest node.
        tight = GH.RPORRTStarSettings(; step_size_m=0.6, collision_sample_ds_m=0.05, neighbor_radius_m=0.01)
        f_tree = GH.RPORRTConnectTree(_SV3(0.0, 0.0, 0.0))
        idx = GH.rpo_rrt_star_add_node!(f_tree, _SV3(0.5, 0.0, 0.0), far, tight; safe_distance_m=0.0)
        @test idx == 2
        @test f_tree.parents[2] == 1
        @test f_tree.costs[2] ≈ 0.5
    end

    @testset "rrt connect planning" begin
        far = _hypr_far_geometry()
        blocked = _hypr_origin_geometry(; keepout_radius_m=0.4)
        cfg = _hypr_cfg(; refinement_enable=true, refinement_rounds=1, refinement_waypoint_passes=2)
        start = _SV3(-1.5, 0.0, 0.0)
        goal = _SV3(1.5, 0.0, 0.0)

        # Free space returns the direct two-point path without any tree work.
        direct = GH.rpo_rrt_connect_plan_path(start, goal, far, cfg; safe_distance_m=0.1, rng=MersenneTwister(1))
        @test direct.path_found
        @test direct.iterations == 0
        @test direct.path == hcat(collect(start), collect(goal))
        @test direct.raw_path == direct.path
        @test direct.cost == direct.raw_cost
        @test direct.cost_history == [direct.cost]
        @test direct.components.violation_count == 0
        @test !direct.refinement_improved

        # Blocked corridor forces the bidirectional search around the sphere.
        settings = GH.RPORRTConnectSettings(;
            n_iters=200, step_size_m=0.5, collision_sample_ds_m=0.05, shortcut_iters=30)
        plan = GH.rpo_rrt_connect_plan_path(
            start, goal, blocked, cfg;
            safe_distance_m=0.1, settings=settings, rng=MersenneTwister(11))
        @test plan.path_found
        @test plan.iterations >= 1
        @test size(plan.raw_path, 2) >= 3
        @test plan.raw_path[:, 1] ≈ collect(start)
        @test plan.raw_path[:, end] ≈ collect(goal)
        @test plan.path[:, 1] ≈ collect(start) atol = 1.0e-9
        @test plan.path[:, end] ≈ collect(goal) atol = 1.0e-9
        @test plan.raw_components.violation_count == 0
        @test isfinite(plan.cost)
        @test plan.config.curve_type == :polyline

        # Goal buried inside the keepout is unreachable: planner reports failure.
        buried_goal = _SV3(0.05, 0.0, 0.0)
        fail = GH.rpo_rrt_connect_plan_path(
            start, buried_goal, blocked, cfg;
            safe_distance_m=0.1,
            settings=GH.RPORRTConnectSettings(; n_iters=15, step_size_m=0.5, collision_sample_ds_m=0.05),
            rng=MersenneTwister(2))
        @test !fail.path_found
        @test fail.iterations >= 1
        @test fail.raw_path == hcat(collect(start), collect(buried_goal))

        # Zero runtime budget stops after the first iteration.
        timed = GH.rpo_rrt_connect_plan_path(
            start, goal, blocked, cfg;
            safe_distance_m=0.1, max_runtime_s=0.0,
            settings=GH.RPORRTConnectSettings(; n_iters=50, step_size_m=0.5, collision_sample_ds_m=0.05),
            rng=MersenneTwister(4))
        @test timed.iterations == 1
        @test !timed.path_found

        # Bezier variant refits the polyline plan to control points.
        bez = GH.rpo_rrt_connect_bezier_plan_path(
            start, goal, blocked, cfg;
            safe_distance_m=0.1, settings=settings, rng=MersenneTwister(11))
        @test bez.path_found
        @test bez.config.curve_type == :bezier
        @test bez.path[:, 1] ≈ collect(start) atol = 1.0e-9
        @test bez.path[:, end] ≈ collect(goal) atol = 1.0e-9
        @test size(bez.smoothed_path, 2) >= 2
        @test bez.smoothed_path[:, 1] ≈ collect(start) atol = 1.0e-9
        @test isfinite(bez.cost)
        @test bez.raw_components.violation_count == 0
    end

    @testset "rrt star planning" begin
        far = _hypr_far_geometry()
        blocked = _hypr_origin_geometry(; keepout_radius_m=0.4)
        cfg = _hypr_cfg(; refinement_enable=true, refinement_rounds=1, refinement_waypoint_passes=2)
        start = _SV3(-1.5, 0.0, 0.0)
        goal = _SV3(1.5, 0.0, 0.0)

        direct = GH.rpo_rrt_star_plan_path(start, goal, far, cfg; safe_distance_m=0.1, rng=MersenneTwister(1))
        @test direct.path_found
        @test direct.iterations == 0
        @test direct.path == hcat(collect(start), collect(goal))
        @test direct.cost == direct.raw_cost

        settings = GH.RPORRTStarSettings(;
            n_iters=400, step_size_m=0.6, collision_sample_ds_m=0.05,
            neighbor_radius_m=1.2, goal_sample_rate=0.15, shortcut_iters=30)
        plan = GH.rpo_rrt_star_plan_path(
            start, goal, blocked, cfg;
            safe_distance_m=0.1, settings=settings, rng=MersenneTwister(5))
        @test plan.path_found
        @test plan.iterations >= 1
        @test plan.raw_path[:, 1] ≈ collect(start)
        @test plan.raw_path[:, end] ≈ collect(goal)
        @test plan.raw_components.violation_count == 0
        @test length(plan.cost_history) == length(plan.history)
        @test length(plan.cost_history) >= 1
        @test isfinite(plan.cost)
        @test plan.path[:, 1] ≈ collect(start) atol = 1.0e-9
        @test plan.path[:, end] ≈ collect(goal) atol = 1.0e-9

        buried_goal = _SV3(0.05, 0.0, 0.0)
        fail = GH.rpo_rrt_star_plan_path(
            start, buried_goal, blocked, cfg;
            safe_distance_m=0.1,
            settings=GH.RPORRTStarSettings(; n_iters=15, step_size_m=0.5, collision_sample_ds_m=0.05),
            rng=MersenneTwister(2))
        @test !fail.path_found
        @test fail.raw_path == hcat(collect(start), collect(buried_goal))
        @test length(fail.cost_history) == 1
        @test length(fail.history) == 1

        timed = GH.rpo_rrt_star_plan_path(
            start, goal, blocked, cfg;
            safe_distance_m=0.1, max_runtime_s=0.0,
            settings=GH.RPORRTStarSettings(; n_iters=50, step_size_m=0.5, collision_sample_ds_m=0.05),
            rng=MersenneTwister(4))
        @test timed.iterations == 1
        @test !timed.path_found
    end

    @testset "pso plan fast paths and diagnostics" begin
        far = _hypr_far_geometry()
        start = _SV3(0.0, 0.0, 0.0)
        goal = _SV3(2.0, 0.0, 0.0)

        # Zero waypoints: straight chord, no swarm, callback fired once at iter 0.
        callback_log = Tuple{Int, Float64}[]
        cfg0 = _hypr_cfg(; n_waypoints=0)
        res0 = GH.rpo_pso_plan_path(
            start, goal, far, cfg0;
            safe_distance_m=0.0, rng=MersenneTwister(1),
            iteration_callback=(iter, cost, comps) -> push!(callback_log, (iter, cost)))
        @test res0.path == hcat(start, goal)
        @test isfinite(res0.cost)
        @test res0.cost ≈ res0.components.total
        @test !res0.warmstart.attempted
        @test res0.cost_history == [res0.cost]
        @test callback_log == [(0, res0.cost)]
        @test !res0.early_stopped && !res0.iteration_timed_out
        @test res0.iteration_timeout_phase == :none

        # Early stopping: impossible improvement thresholds trip patience at iter 2.
        cfg_es = _hypr_cfg(;
            n_iters=10, early_stopping_enable=true, early_stopping_patience=1,
            early_stopping_min_iters=1, early_stopping_min_abs_improvement=1.0e12,
            early_stopping_min_rel_improvement=1.0e12)
        res_es = GH.rpo_pso_plan_path(start, goal, far, cfg_es; safe_distance_m=0.0, rng=MersenneTwister(2))
        @test res_es.early_stopped
        @test res_es.early_stop_iter == 2
        @test length(res_es.cost_history) == 2
        @test !res_es.iteration_timed_out

        # Iteration runtime budget of ~zero times out in the first phase checkpoint;
        # with reexploration disabled the phase reports the iteration start, not :reexplore.
        cfg_to = _hypr_cfg(; n_iters=3, iteration_runtime_limit_s=1.0e-9)
        res_to = GH.rpo_pso_plan_path(start, goal, far, cfg_to; safe_distance_m=0.0, rng=MersenneTwister(3))
        @test res_to.iteration_timed_out
        @test res_to.iteration_timeout_iter == 1
        @test res_to.iteration_timeout_phase == :iteration_start
        @test length(res_to.iteration_timeout_events) == 1
        @test res_to.iteration_timeout_events[1].iter == 1
        @test res_to.iteration_timeout_events[1].elapsed_s >= 0.0
        @test length(res_to.cost_history) == 1
        @test !res_to.early_stopped
    end

    @testset "pso plan full runs" begin
        far = _hypr_far_geometry()
        blocked = _hypr_origin_geometry(; keepout_radius_m=0.6)
        start = _SV3(-1.5, 0.0, 0.0)
        goal = _SV3(1.5, 0.0, 0.0)

        # Free-space polyline run with schedule, culling, stagnation learning,
        # refinement, and a finite (generous) per-iteration budget.
        callback_iters = Int[]
        cfg_main = _hypr_cfg(;
            n_waypoints=2, n_particles=6, n_iters=5,
            schedule_enable=true,
            cull_enable=true, cull_start_iter=1, cull_fraction_max=0.4,
            stagnation_learning_enable=true, stagnation_learning_threshold=1,
            stagnation_learning_elite_fraction=0.2,
            refinement_enable=true, refinement_rounds=1, refinement_waypoint_passes=2,
            iteration_runtime_limit_s=60.0)
        res = GH.rpo_pso_plan_path(
            start, goal, far, cfg_main;
            safe_distance_m=0.0, rng=MersenneTwister(7),
            iteration_callback=(iter, cost, comps) -> push!(callback_iters, iter))
        @test res.path[:, 1] ≈ collect(start) atol = 1.0e-9
        @test res.path[:, end] ≈ collect(goal) atol = 1.0e-9
        @test isfinite(res.cost)
        @test res.components.violation_count == 0
        @test length(res.cost_history) == cfg_main.n_iters
        @test all(diff(res.cost_history) .<= 1.0e-9)  # global best never worsens
        @test callback_iters == collect(1:cfg_main.n_iters)
        @test !res.adaptive.enabled
        @test !res.iteration_timed_out && !res.early_stopped
        # Free-space optimum is nearly the straight chord.
        @test GH.rpo_path_length(res.path) < 1.25 * norm(goal - start)

        # Blocked corridor triggers reexploration growth of waypoints and margin.
        cfg_re = _hypr_cfg(;
            n_waypoints=1, n_particles=6, n_iters=4,
            search_margin_m=0.3, spread_scale=0.1,
            reexplore_enable=true, reexplore_trigger_iter=1,
            reexplore_search_margin_scale=2.0, reexplore_waypoint_scale=1.5,
            reexplore_waypoint_increment=2, reexplore_max_waypoints=3)
        res_re = GH.rpo_pso_plan_path(
            start, goal, blocked, cfg_re;
            safe_distance_m=0.05, rng=MersenneTwister(21))
        @test res_re.config.n_waypoints == 3
        @test res_re.config.search_margin_m > cfg_re.search_margin_m
        @test length(res_re.cost_history) == cfg_re.n_iters
        @test res_re.path[:, 1] ≈ collect(start) atol = 1.0e-9
        @test res_re.path[:, end] ≈ collect(goal) atol = 1.0e-9

        # Bezier planning seeded by an RRT-Connect warm start.
        cfg_ws = _hypr_cfg(;
            n_waypoints=2, n_particles=8, n_iters=3,
            curve_type=:bezier,
            cull_enable=true, cull_start_iter=1, cull_fraction_max=0.4,
            rrt_warmstart_enable=true, rrt_warmstart_iters=200,
            rrt_warmstart_step_size_m=0.5, rrt_warmstart_collision_sample_ds_m=0.05,
            rrt_warmstart_shortcut_iters=20,
            reexplore_max_waypoints=6)
        res_ws = GH.rpo_pso_plan_path(
            _SV3(-1.5, 0.0, 0.0), _SV3(1.5, 0.0, 0.0),
            _hypr_origin_geometry(; keepout_radius_m=0.4), cfg_ws;
            safe_distance_m=0.1, rng=MersenneTwister(13))
        @test res_ws.warmstart.enabled && res_ws.warmstart.attempted
        @test res_ws.warmstart.path_found
        @test res_ws.warmstart.n_points >= 3  # had to route around the sphere
        @test res_ws.config.curve_type == :bezier
        @test res_ws.config.n_waypoints >= cfg_ws.n_waypoints
        @test res_ws.config.n_waypoints <= 6
        @test res_ws.path[:, 1] ≈ [-1.5, 0.0, 0.0] atol = 1.0e-9
        @test res_ws.path[:, end] ≈ [1.5, 0.0, 0.0] atol = 1.0e-9
        @test isfinite(res_ws.cost)

        # Adaptive policy integrated into the planner entry point.
        cfg_adapt = _hypr_cfg(;
            n_waypoints=1, n_particles=6, n_iters=2,
            adaptive_enable=true, adaptive_allow_downscale=true,
            adaptive_n_waypoints_min=1, adaptive_n_waypoints_max=2,
            adaptive_n_particles_min=4, adaptive_n_particles_max=8,
            adaptive_n_iters_min=1, adaptive_n_iters_max=2)
        res_ad = GH.rpo_pso_plan_path(start, goal, far, cfg_adapt; rng=MersenneTwister(31))
        @test res_ad.adaptive.enabled
        @test 0.0 <= res_ad.adaptive.explore <= 1.0
        @test 4 <= res_ad.config.n_particles <= 8
        @test isfinite(res_ad.cost)
        @test res_ad.path[:, 1] ≈ collect(start) atol = 1.0e-9
        @test res_ad.path[:, end] ≈ collect(goal) atol = 1.0e-9
    end
end
