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

    @testset "adaptive policy: complexity and probe metrics" begin
        far = _hypr_far_geometry()
        blocked = _hypr_origin_geometry(; keepout_radius_m=0.5)

        # Free corridor: only the smooth clearance term contributes.
        c_free = GH.rpo_estimate_geometry_complexity((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), far)
        @test c_free ≈ 0.3 / (1.0 + 98.99) atol = 1.0e-9
        # Corridor through the keepout sphere saturates to full complexity.
        c_hit = GH.rpo_estimate_geometry_complexity((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), blocked)
        @test c_hit == 1.0
        # Positive clearance below the safety buffer adds the buffer term.
        c_buf = GH.rpo_estimate_geometry_complexity(
            (-1.0, 0.6, 0.0), (1.0, 0.6, 0.0), blocked; safe_distance_m=0.2)
        @test c_buf ≈ 0.7 + 0.3 / (1.0 + (0.6 - 0.51)) atol = 1.0e-9
        @test 0.0 <= c_buf <= 1.0

        ok = GH.rpo_probe_geometry_metrics((-1.0, 0.6, 0.0), (1.0, 0.6, 0.0), blocked)
        @test ok.detour_ratio == 1.0
        @test ok.min_clearance ≈ 0.6 - 0.51 atol = 1.0e-9
        @test ok.violation_fraction == 0.0
        @test ok.success  # safe_distance defaults to zero
        tight = GH.rpo_probe_geometry_metrics((-1.0, 0.6, 0.0), (1.0, 0.6, 0.0), blocked; safe_distance_m=0.2)
        @test !tight.success
        crossing = GH.rpo_probe_geometry_metrics((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), blocked)
        @test crossing.min_clearance < 0.0
        @test crossing.violation_fraction > 0.0
        @test !crossing.success
    end

    @testset "adaptive policy: config adaptation" begin
        far = _hypr_far_geometry()
        start = (0.0, 0.0, 0.0)
        goal = (1.0, 0.0, 0.0)

        base_off = _hypr_cfg()
        cfg_off, info_off = GH.rpo_adaptive_pso_config(base_off, start, goal, far)
        @test !info_off.enabled
        @test info_off.complexity == 0.0
        @test info_off.explore == 0.0
        @test info_off.distance_m ≈ 1.0
        @test cfg_off.n_particles == base_off.n_particles
        @test cfg_off.n_iters == base_off.n_iters

        # Downscale allowed with zero maxima triggers the raw fallback sizing.
        base_dn = _hypr_cfg(;
            adaptive_enable=true,
            adaptive_allow_downscale=true,
            adaptive_n_particles_max=0,
            adaptive_n_iters_max=0,
            adaptive_n_particles_min=2,
            adaptive_n_iters_min=1,
            adaptive_n_waypoints_min=1,
            adaptive_n_waypoints_max=2,
        )
        cfg_dn, info_dn = GH.rpo_adaptive_pso_config(base_dn, start, goal, far)
        @test info_dn.enabled
        @test 0.0 <= info_dn.complexity <= 1.0
        @test 0.0 <= info_dn.explore <= 1.0
        @test info_dn.distance_m ≈ 1.0
        @test 1 <= cfg_dn.n_waypoints <= 2
        @test 2 <= cfg_dn.n_particles <= max(20, 3 * base_dn.n_particles)
        @test 1 <= cfg_dn.n_iters <= max(5, 3 * base_dn.n_iters)
        @test base_dn.adaptive_c1_min <= cfg_dn.c1 <= base_dn.adaptive_c1_max
        @test base_dn.adaptive_c2_min <= cfg_dn.c2 <= base_dn.adaptive_c2_max
        @test base_dn.adaptive_w_inertia_min <= cfg_dn.w_inertia <= base_dn.adaptive_w_inertia_max
        @test base_dn.adaptive_spread_scale_min <= cfg_dn.spread_scale <= base_dn.adaptive_spread_scale_max
        @test base_dn.adaptive_w_len_min <= cfg_dn.w_len <= base_dn.adaptive_w_len_max

        # Downscale disallowed keeps at least the base sizing even with tiny maxima.
        base_up = _hypr_cfg(;
            adaptive_enable=true,
            adaptive_allow_downscale=false,
            adaptive_n_particles_max=4,
            adaptive_n_iters_max=1,
            adaptive_n_particles_min=2,
            adaptive_n_iters_min=1,
            adaptive_n_waypoints_min=1,
            adaptive_n_waypoints_max=1,
        )
        cfg_up, info_up = GH.rpo_adaptive_pso_config(base_up, start, goal, far)
        @test info_up.enabled
        @test cfg_up.n_particles >= base_up.n_particles
        @test cfg_up.n_iters >= base_up.n_iters
        @test cfg_up.n_waypoints >= base_up.n_waypoints
    end

    @testset "refinement primitives" begin
        far = _hypr_far_geometry()
        blocked = _hypr_origin_geometry(; keepout_radius_m=0.5)
        cfg = _hypr_cfg(; sample_ds_m=0.25, refinement_sample_ds_m=0.1)

        refine_cfg = GH.rpo_refinement_config(cfg)
        @test refine_cfg.sample_ds_m == 0.1
        @test refine_cfg.refinement_sample_ds_m == 0.1

        better = (J_obs=0.0, total=1.0)
        worse = (J_obs=0.0, total=2.0)
        @test GH.rpo_refinement_better(better, worse, cfg)
        @test !GH.rpo_refinement_better(worse, better, cfg)
        @test !GH.rpo_refinement_better(better, better, cfg)  # zero improvement fails thresholds
        obs_worse = (J_obs=5.0, total=0.1)
        @test !GH.rpo_refinement_better(obs_worse, worse, cfg)  # obstacle regression vetoes

        wild = [0.0 10.0 2.0; 0.0 -50.0 0.0; 0.0 0.0 0.0]
        clamped = GH.rpo_refinement_clamp_path(wild, cfg)
        lo, hi = GH.rpo_pso_bounds(wild[:, 1], wild[:, end], cfg)
        @test clamped[:, 1] == wild[:, 1]
        @test clamped[:, end] == wild[:, end]
        @test all(lo .<= clamped[:, 2] .<= hi)
        @test clamped[2, 2] == lo[2]

        seg = GH.rpo_refinement_segment_samples((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.4)
        @test size(seg, 2) == 4  # ceil(1/0.4)+1
        @test seg[:, 1] ≈ [0.0, 0.0, 0.0]
        @test seg[:, end] ≈ [1.0, 0.0, 0.0]
        @test issorted(seg[1, :])
        seg_one = GH.rpo_refinement_segment_samples((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 10.0)
        @test size(seg_one, 2) == 2
        seg_degen = GH.rpo_refinement_segment_samples((1.0, 1.0, 1.0), (1.0, 1.0, 1.0), 0.1)
        @test size(seg_degen, 2) == 2
        @test seg_degen[:, 1] ≈ seg_degen[:, 2]

        # Uniform-sampling safety check.
        @test GH.rpo_refinement_segment_is_safe((-1.0, 2.0, 0.0), (1.0, 2.0, 0.0), blocked, cfg)
        @test !GH.rpo_refinement_segment_is_safe((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), blocked, cfg)
        # Adaptive-sampling safety check hits the other branch.
        cfg_ad = _hypr_cfg(; adaptive_sampling_enable=true, sample_ds_m=0.1)
        @test GH.rpo_refinement_segment_is_safe((-1.0, 2.0, 0.0), (1.0, 2.0, 0.0), blocked, cfg_ad)
        @test !GH.rpo_refinement_segment_is_safe(
            (-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), blocked, cfg_ad; safe_distance_m=0.05)

        @test GH.rpo_refinement_bernstein(3, 0, 0.0) == 1.0
        @test GH.rpo_refinement_bernstein(3, 3, 1.0) == 1.0
        @test sum(GH.rpo_refinement_bernstein(3, j, 0.3) for j in 0:3) ≈ 1.0

        @test GH.rpo_refinement_sample_params(reshape([1.0, 2.0, 3.0], 3, 1)) == [0.0]
        same = [1.0 1.0 1.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        @test GH.rpo_refinement_sample_params(same) == [0.0, 0.5, 1.0]  # degenerate arc length
        line3 = [0.0 1.0 3.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        u = GH.rpo_refinement_sample_params(line3)
        @test u ≈ [0.0, 1.0 / 3.0, 1.0]

        dense_line = GH.rpo_refinement_segment_samples((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), 0.1)
        controls2 = GH.rpo_fit_bezier_fixed_endpoints(dense_line, 2, cfg)
        @test size(controls2) == (3, 2)
        @test controls2[:, 1] ≈ [0.0, 0.0, 0.0]
        @test controls2[:, 2] ≈ [2.0, 0.0, 0.0]
        controls4 = GH.rpo_fit_bezier_fixed_endpoints(dense_line, 4, cfg)
        @test size(controls4) == (3, 4)
        @test controls4[:, 1] ≈ dense_line[:, 1]
        @test controls4[:, end] ≈ dense_line[:, end]
        @test all(abs.(controls4[2, :]) .< 1.0e-6)  # collinear controls for straight data
        @test all(abs.(controls4[3, :]) .< 1.0e-6)
        # n_control below 2 clamps up to endpoints only.
        @test size(GH.rpo_fit_bezier_fixed_endpoints(dense_line, 1, cfg), 2) == 2

        a = _SV3(0.0, 0.0, 0.0)
        b = _SV3(1.0, 0.0, 0.0)
        @test GH.rpo_refinement_project_to_segment((0.5, 1.0, 0.0), a, b) ≈ _SV3(0.5, 0.0, 0.0)
        @test GH.rpo_refinement_project_to_segment((-2.0, 1.0, 0.0), a, b) ≈ a
        @test GH.rpo_refinement_project_to_segment((5.0, -1.0, 0.0), a, b) ≈ b
        @test GH.rpo_refinement_project_to_segment((0.3, 0.4, 0.0), a, a) ≈ a  # degenerate segment
    end

    @testset "refinement passes" begin
        far = _hypr_far_geometry()
        cfg = _hypr_cfg(; refinement_enable=true, sample_ds_m=0.1, refinement_sample_ds_m=0.1)

        zigzag = [0.0 1.0 2.0 3.0; 0.0 1.0 -1.0 0.0; 0.0 0.0 0.0 0.0]
        straightened = GH.rpo_refinement_shortcut_samples(zigzag, far, cfg)
        @test size(straightened, 2) == 2  # whole zigzag collapses in free space
        @test straightened[:, 1] ≈ zigzag[:, 1]
        @test straightened[:, end] ≈ zigzag[:, end]
        two = [0.0 1.0; 0.0 0.0; 0.0 0.0]
        @test GH.rpo_refinement_shortcut_samples(two, far, cfg) == two

        zig_comps = GH.rpo_normalized_path_cost_components(zigzag, far, cfg)
        straight = [0.0 1.0 2.0 3.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
        accepted, comps, ok = GH.rpo_try_accept_refinement(straight, far, cfg, zig_comps)
        @test ok
        @test comps.total < zig_comps.total
        @test accepted[:, 1] ≈ zigzag[:, 1]
        rejected, kept, no = GH.rpo_try_accept_refinement(zigzag, far, cfg, zig_comps)
        @test !no
        @test rejected === nothing
        @test kept === zig_comps

        # Shortcut+refit accepts a straightened candidate in free space.
        cand, cand_comps, did = GH.rpo_refine_shortcut_refit(zigzag, far, cfg, zig_comps)
        @test did
        @test cand_comps.total < zig_comps.total
        @test cand[:, 1] ≈ zigzag[:, 1]
        @test cand[:, end] ≈ zigzag[:, end]

        # A blocked elbow with coarse refinement sampling cannot shortcut at all.
        blocked = _hypr_origin_geometry(; keepout_radius_m=0.9)
        elbow = [-1.0 0.0 1.0; 0.0 1.0 0.0; 0.0 0.0 0.0]
        cfg_coarse = _hypr_cfg(; sample_ds_m=0.05, refinement_sample_ds_m=5.0)
        elbow_comps = GH.rpo_normalized_path_cost_components(elbow, blocked, cfg_coarse)
        same_path, same_comps, changed = GH.rpo_refine_shortcut_refit(elbow, blocked, cfg_coarse, elbow_comps)
        @test !changed
        @test same_path == elbow
        @test same_comps === elbow_comps

        # Handle tightening pulls an off-chord handle toward the chord.
        bent = [0.0 1.5 3.0; 0.0 2.0 0.0; 0.0 0.0 0.0]
        bent_comps = GH.rpo_normalized_path_cost_components(bent, far, cfg)
        tight, tight_comps, timproved = GH.rpo_refine_tighten_handles(bent, far, cfg, bent_comps)
        @test timproved
        @test tight_comps.total < bent_comps.total
        @test tight[:, 1] ≈ bent[:, 1]
        @test tight[:, end] ≈ bent[:, end]
        @test abs(tight[2, 2]) < 2.0  # handle moved toward the chord
        short = [0.0 1.0; 0.0 0.0; 0.0 0.0]
        s_path, s_comps, s_improved = GH.rpo_refine_tighten_handles(short, far, cfg, bent_comps)
        @test !s_improved
        @test s_path == short

        # Degree lowering finds the straight 2-control representation.
        detour = [0.0 1.0 2.0 3.0; 0.0 2.0 2.0 0.0; 0.0 0.0 0.0 0.0]
        detour_comps = GH.rpo_normalized_path_cost_components(detour, far, cfg)
        lower, lower_comps, limproved = GH.rpo_refine_lower_degree(detour, far, cfg, detour_comps)
        @test limproved
        @test lower_comps.total < detour_comps.total
        @test size(lower, 2) < size(detour, 2)
        @test lower[:, 1] ≈ detour[:, 1]
        @test lower[:, end] ≈ detour[:, end]
        tri = [0.0 1.0 2.0; 0.0 1.0 0.0; 0.0 0.0 0.0]
        t_path, t_comps, t_improved = GH.rpo_refine_lower_degree(tri, far, cfg, detour_comps)
        @test !t_improved
        @test t_path == tri
    end

    @testset "post refine path" begin
        far = _hypr_far_geometry()
        zigzag = [0.0 1.0 2.0 3.0; 0.0 1.0 -1.0 0.0; 0.0 0.0 0.0 0.0]

        cfg_off = _hypr_cfg(; refinement_enable=false)
        path_off, cost_off, improved_off = GH.rpo_post_refine_path(zigzag, far, cfg_off)
        @test !improved_off
        @test path_off == zigzag
        @test cost_off ≈ GH.rpo_path_cost(zigzag, far, cfg_off)

        cfg_zero = _hypr_cfg(; refinement_enable=true, refinement_rounds=0)
        _, cost_zero, improved_zero = GH.rpo_post_refine_path(zigzag, far, cfg_zero)
        @test !improved_zero
        @test cost_zero ≈ GH.rpo_path_cost(zigzag, far, cfg_zero)

        cfg_on = _hypr_cfg(; refinement_enable=true, sample_ds_m=0.1)
        initial_cost = GH.rpo_path_cost(zigzag, far, cfg_on)
        refined, cost_on, improved_on = GH.rpo_post_refine_path(zigzag, far, cfg_on)
        @test improved_on
        @test cost_on < initial_cost
        @test refined[:, 1] ≈ zigzag[:, 1]
        @test refined[:, end] ≈ zigzag[:, end]
        # Refined free-space path is nearly straight.
        @test GH.rpo_path_length(refined) < GH.rpo_path_length(zigzag)
    end

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
        @test GH.rpo_rrt_shortcut_path(zigzag[:, [1, 4]], far, settings; safe_distance_m=0.0) == zigzag[:, [1, 4]]
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
        direct = GH.rpo_rrt_connect_plan_path(start, goal, far, cfg; safe_distance_m=0.1)
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

        direct = GH.rpo_rrt_star_plan_path(start, goal, far, cfg; safe_distance_m=0.1)
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

    @testset "pso swarm helpers" begin
        far = _hypr_far_geometry()

        cfg_flat = _hypr_cfg(; schedule_enable=false)
        w = GH.rpo_pso_iteration_weights(cfg_flat, 2)
        @test w == (w_inertia=cfg_flat.w_inertia, c1=cfg_flat.c1, c2=cfg_flat.c2)
        cfg_sched = _hypr_cfg(; schedule_enable=true, n_iters=10)
        w_last = GH.rpo_pso_iteration_weights(cfg_sched, 10)
        @test w_last.w_inertia < cfg_sched.w_inertia

        a = _SV3(0.0, 0.0, 0.0)
        b = _SV3(2.0, 0.0, 0.0)
        @test GH.rpo_pso_project_to_segment((1.0, 3.0, 0.0), a, b) ≈ _SV3(1.0, 0.0, 0.0)
        @test GH.rpo_pso_project_to_segment((-9.0, 0.0, 0.0), a, b) ≈ a
        @test GH.rpo_pso_project_to_segment((1.0, 1.0, 1.0), a, a) ≈ a

        @test GH.rpo_pso_tapered_noise_scale(1, 1) == 1.0
        @test GH.rpo_pso_tapered_noise_scale(1, 5) == clamp(1.0 / 3.0, 0.25, 1.0)
        @test GH.rpo_pso_tapered_noise_scale(3, 5) == 1.0  # middle waypoint keeps full noise

        # Culling: disabled, warm-up, and unevaluated swarms are all no-ops.
        start = _SV3(0.0, 0.0, 0.0)
        goal = _SV3(1.0, 0.0, 0.0)
        mk_state() = (
            zeros(3, 4), ones(3, 4), zeros(3, 4), [1.0, 2.0, 3.0, 4.0],
            -ones(3), ones(3), zeros(3),
        )
        cfg_cull = _hypr_cfg(;
            n_particles=4, cull_enable=true, cull_start_iter=1, cull_fraction_max=0.5,
            cull_noise_scale=0.1, cull_arc_velocity_scale=0.2)
        pos, vel, pb, pbc, lo, hi, gb = mk_state()
        cfg_no = _hypr_cfg(; n_particles=4, cull_enable=false)
        @test GH.rpo_pso_cull_swarm!(pos, vel, pb, pbc, lo, hi, gb, start, goal, 1, cfg_no, 5, MersenneTwister(1)) == 0
        @test GH.rpo_pso_cull_swarm!(pos, vel, pb, pbc, lo, hi, gb, start, goal, 1, cfg_cull, 0, MersenneTwister(1)) == 0
        inf_costs = fill(Inf, 4)
        @test GH.rpo_pso_cull_swarm!(pos, vel, pb, inf_costs, lo, hi, gb, start, goal, 1, cfg_cull, 5, MersenneTwister(1)) == 0

        # Polyline branch replaces the two worst particles with jittered gbest copies.
        pos, vel, pb, pbc, lo, hi, gb = mk_state()
        n = GH.rpo_pso_cull_swarm!(pos, vel, pb, pbc, lo, hi, gb, start, goal, 1, cfg_cull, 3, MersenneTwister(1))
        @test n == 2
        @test pbc[1] == 1.0 && pbc[2] == 2.0
        @test pbc[3] == Inf && pbc[4] == Inf
        @test all(vel[:, 3] .== 0.0) && all(vel[:, 4] .== 0.0)
        @test all(lo .<= pos[:, 3] .<= hi) && all(lo .<= pos[:, 4] .<= hi)
        @test pb[:, 3] == pos[:, 3] && pb[:, 4] == pos[:, 4]

        # Bezier branch reseeds along corridor-biased targets with arc velocities.
        cfg_bez = _hypr_cfg(;
            n_particles=4, curve_type=:bezier, cull_enable=true, cull_start_iter=1,
            cull_fraction_max=0.5, cull_noise_scale=0.1, cull_arc_velocity_scale=0.2)
        pos, vel, pb, pbc, lo, hi, gb = mk_state()
        n = GH.rpo_pso_cull_swarm!(pos, vel, pb, pbc, lo, hi, gb, start, goal, 1, cfg_bez, 3, MersenneTwister(1))
        @test n == 2
        @test pbc[3] == Inf && pbc[4] == Inf
        @test all(lo .<= pos[:, 3] .<= hi) && all(lo .<= pos[:, 4] .<= hi)
        @test pb[:, 3] == pos[:, 3] && pb[:, 4] == pos[:, 4]

        # Effective safety distance resolution.
        cfg_sd = _hypr_cfg(; safe_distance_m=0.3, sample_ds_m=0.3)
        @test GH.rpo_pso_effective_safe_distance(cfg_sd, nothing) == 0.3
        @test GH.rpo_pso_effective_safe_distance(cfg_sd, 0.5) == 0.5
        @test GH.rpo_pso_effective_safe_distance(cfg_sd, 0.0) == 0.3  # zero defers to config
        cfg_sd0 = _hypr_cfg()
        @test GH.rpo_pso_effective_safe_distance(cfg_sd0, 0.0) == 0.0

        diag = GH.rpo_pso_empty_warmstart_diagnostics(_hypr_cfg(; rrt_warmstart_enable=true))
        @test diag.enabled && !diag.attempted && !diag.path_found
        @test diag.iterations == 0 && diag.cost == Inf && diag.n_points == 0

        # Warm start generation: disabled and free-space direct hit.
        p_off, d_off = GH.rpo_pso_rrt_warmstart_path(start, goal, far, _hypr_cfg(), 0.0, MersenneTwister(1))
        @test p_off === nothing
        @test !d_off.enabled && !d_off.attempted
        cfg_ws = _hypr_cfg(; rrt_warmstart_enable=true, rrt_warmstart_iters=50)
        p_on, d_on = GH.rpo_pso_rrt_warmstart_path(start, goal, far, cfg_ws, 0.0, MersenneTwister(1))
        @test d_on.attempted && d_on.path_found
        @test d_on.n_points == 2
        @test p_on == hcat(collect(start), collect(goal))
        @test isfinite(d_on.cost)

        mask = GH.rpo_pso_protected_particle_mask([3.0, 1.0, 2.0], 0.34)
        @test mask == [false, true, true]

        # Failed stagnation-learning attempts keep the particle eligible for another
        # attempt next iteration; only an accepted attempt clears the counter.
        @test GH.rpo_pso_stagnation_count_after_learning(3, false) == 3
        @test GH.rpo_pso_stagnation_count_after_learning(3, true) == 0
        @test GH.rpo_pso_stagnation_count_after_learning(0, false) == 0

        cfg_imp = _hypr_cfg()
        @test GH.rpo_pso_material_improvement(0.5, 1.0, cfg_imp)
        @test !GH.rpo_pso_material_improvement(1.0, 1.0, cfg_imp)
        @test GH.rpo_pso_material_improvement(1.0, Inf, cfg_imp)

        cfg_feas = _hypr_cfg(; early_stopping_require_feasible=true)
        @test !GH.rpo_pso_early_stopping_feasible(nothing, cfg_feas)
        @test GH.rpo_pso_early_stopping_feasible((violation_count=0,), cfg_feas)
        @test !GH.rpo_pso_early_stopping_feasible((violation_count=3,), cfg_feas)
        cfg_lax = _hypr_cfg(; early_stopping_require_feasible=false)
        @test GH.rpo_pso_early_stopping_feasible(nothing, cfg_lax)
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
