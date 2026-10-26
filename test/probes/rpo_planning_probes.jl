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
const HU = SimulationModel.HYPRUtils

"Small station geometry far away from the RTN origin so it does not disturb local paths."
function _probe_far_geometry(; keepout_radius_m::Real=0.0)
    station = SM.RPOStationGeometry(
        reshape([100.0, 0.0, 0.0], 3, 1);
        keepout_radius_m=keepout_radius_m,
        name="probe_far",
    )
    return SM.RPOReferenceGeometry(station; chaser=SM.RPOCubeSatGeometry(dims_m=(0.02, 0.02, 0.02)))
end

"Station geometry with a single point at the origin."
function _probe_origin_geometry(; keepout_radius_m::Real=0.0)
    station = SM.RPOStationGeometry(
        reshape([0.0, 0.0, 0.0], 3, 1);
        keepout_radius_m=keepout_radius_m,
        name="probe_origin",
    )
    return SM.RPOReferenceGeometry(station; chaser=SM.RPOCubeSatGeometry(dims_m=(0.02, 0.02, 0.02)))
end

"Cheap deterministic-ish PSO config for fast plan builds."
function _probe_tiny_cfg(; kwargs...)
    return SM.RPOPSOConfig(;
        n_waypoints=1,
        n_particles=8,
        n_iters=2,
        sample_ds_m=0.25,
        curve_type=:polyline,
        adaptive_enable=false,
        adaptive_sampling_enable=false,
        probe_enable=false,
        refinement_enable=false,
        retime_a_max_mps2=0.1,
        retime_dt_s=1.0,
        rrt_warmstart_iters=120,
        rrt_warmstart_step_size_m=0.5,
        kwargs...,
    )
end

"Straight-line reference plan along the RTN x axis used by replanning decision tests."
function _probe_straight_plan(x0::Float64, x1::Float64; ds::Float64=0.25)
    xs = collect(x0:ds:x1)
    r_ref = zeros(3, length(xs))
    r_ref[1, :] .= xs
    return GH.RPOPlan(
        valid=true,
        t_ref_s=collect(0.0:1.0:(length(xs) - 1)),
        r_ref_rtn=r_ref,
        v_ref_rtn=zeros(3, length(xs)),
        path_rtn=[x0 x1; 0.0 0.0; 0.0 0.0],
        cost=1.0,
    )
end

@testset "coverage debt rpo planning probes" begin

    @testset "HYPRUtils path helpers" begin
        pts = [0.0 3.0; 0.0 4.0; 0.0 0.0]
        @test HU.hypr_path_length(pts) ≈ 5.0
        @test HU.hypr_path_length(reshape([1.0, 2.0, 3.0], 3, 1)) == 0.0

        # Linear Bezier midpoint is the chord midpoint.
        mid = HU.hypr_bezier_point(pts, 0.5)
        @test mid ≈ [1.5, 2.0, 0.0]
        @test HU.hypr_bezier_point(pts, 0.0) ≈ pts[:, 1]
        @test HU.hypr_bezier_point(pts, 1.0) ≈ pts[:, 2]

        # Quadratic Bezier at t=0.5 equals 0.25 p0 + 0.5 p1 + 0.25 p2.
        quad = [0.0 1.0 2.0; 0.0 2.0 0.0; 0.0 0.0 0.0]
        @test HU.hypr_bezier_point(quad, 0.5) ≈ [1.0, 1.0, 0.0]

        out = zeros(3)
        work = similar(quad)
        HU.hypr_bezier_point!(out, work, quad, 0.25)
        @test out ≈ HU.hypr_bezier_point(quad, 0.25)

        samples_bez = HU.hypr_sample_count_path(pts, 5; curve_type=:bezier)
        @test size(samples_bez) == (3, 5)
        @test samples_bez[:, 1] ≈ pts[:, 1]
        @test samples_bez[:, end] ≈ pts[:, 2]
        @test samples_bez[:, 3] ≈ [1.5, 2.0, 0.0]

        tri = [0.0 1.0 1.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
        samples_poly = HU.hypr_sample_count_path(tri, 5; curve_type=:polyline)
        @test samples_poly[:, 3] ≈ tri[:, 2]  # midpoint sample lands on the interior knot
        @test samples_poly[:, 2] ≈ [0.5, 0.0, 0.0]
        @test_throws ArgumentError HU.hypr_sample_count_path(pts, 1)
        @test_throws ArgumentError HU.hypr_sample_count_path(pts, 4; curve_type=:spline)
    end

    @testset "HYPRUtils RRT helpers" begin
        nodes = [SVector{2, Float64}(0.0, 0.0), SVector{2, Float64}(1.0, 0.0), SVector{2, Float64}(2.0, 0.0)]
        tree = (nodes=nodes, parents=[0, 1, 2], costs=[0.0, 1.0, 2.0])

        @test HU.hypr_rrt_nearest_index(tree, SVector{2, Float64}(1.9, 0.0)) == 3
        @test HU.hypr_rrt_near_indices(tree, SVector{2, Float64}(1.0, 0.0), 1.05) == [1, 2, 3]
        @test HU.hypr_rrt_near_indices(tree, SVector{2, Float64}(10.0, 0.0), 0.5) == Int[]

        q, status = HU.hypr_rrt_steer(nodes[1], nodes[1], 0.5)
        @test status == :trapped && q == nodes[1]
        q, status = HU.hypr_rrt_steer(nodes[1], nodes[2], 5.0)
        @test status == :reached && q == nodes[2]
        q, status = HU.hypr_rrt_steer(nodes[1], nodes[3], 0.5)
        @test status == :advanced
        @test q ≈ SVector{2, Float64}(0.5, 0.0)

        path = HU.hypr_rrt_tree_path(tree, 3)
        @test path == [0.0 1.0 2.0; 0.0 0.0 0.0]

        goal_tree = (nodes=[SVector{2, Float64}(4.0, 0.0), SVector{2, Float64}(3.0, 0.0), SVector{2, Float64}(2.0, 0.0)],
                     parent=[0, 1, 2], cost=[0.0, 1.0, 2.0])
        joined = HU.hypr_rrt_join_paths(tree, 3, goal_tree, 3)
        @test joined == [0.0 1.0 2.0 3.0 4.0; 0.0 0.0 0.0 0.0 0.0]

        # Rewire bookkeeping via robot-arm-style :parent/:cost fields.
        rew = (nodes=[SVector{2, Float64}(0.0, 0.0), SVector{2, Float64}(2.0, 0.0), SVector{2, Float64}(2.0, 1.0)],
               parent=[0, 1, 2], cost=[0.0, 99.0, 99.0])
        HU.hypr_rrt_refresh_subtree_costs!(rew, 1)
        @test rew.cost ≈ [0.0, 2.0, 3.0]
    end

    @testset "RPO plan buffer" begin
        plan = GH.RPOPlan()
        @test !plan.valid
        @test plan.cost == Inf
        @test size(plan.r_ref_rtn) == (3, 0)

        buffer = GH.RPOPlanBuffer()
        @test !buffer.valid
        @test isnan(buffer.updated_at_s)

        valid_plan = _probe_straight_plan(0.0, 1.0)
        out = GH.update_rpo_plan_buffer!(buffer, valid_plan, 12.5)
        @test out === buffer
        @test buffer.valid
        @test buffer.plan === valid_plan
        @test buffer.updated_at_s == 12.5

        GH.update_rpo_plan_buffer!(buffer, GH.RPOPlan(), 13.0)
        @test !buffer.valid
        @test buffer.updated_at_s == 13.0
    end

    @testset "RPO PSO parameters" begin
        cfg = SM.rpo_pso_config()
        @test cfg.n_particles == 200
        @test cfg.curve_type == :bezier

        aliased = SM.rpo_pso_config(cfg;
            pso_particles=16,
            pso_iters=3,
            pso_curve_type=:polyline,
            pso_safe_distance=0.2,
        )
        @test aliased.n_particles == 16
        @test aliased.n_iters == 3
        @test aliased.curve_type == :polyline
        @test aliased.safe_distance_m == 0.2
        # Sampling density syncs with the keepout distance.
        @test aliased.sample_ds_m == 0.2

        configurator = SM.RPOPSOConfigurator(
            swarm=SM.RPOPSOSwarmSettings(n_particles=10, n_iters=4),
            retiming=SM.RPOPSORetimingSettings(dt_s=0.5),
        )
        built = SM.RPOPSOConfig(configurator)
        @test built.n_particles == 10
        @test built.n_iters == 4
        @test built.retime_dt_s == 0.5
        via_helper = SM.rpo_pso_config(configurator; pso_particles=11)
        @test via_helper.n_particles == 11

        @test GH.rpo_hypr_sampling_density_m(SM.RPOPSOConfig(sample_ds_m=0.5), 0.2) == 0.2
        @test GH.rpo_hypr_sampling_density_m(SM.rpo_pso_config(; safe_distance_m=0.3)) == 0.3
        @test GH.rpo_hypr_sampling_density_m(SM.RPOPSOConfig(sample_ds_m=0.5), 0.0) == 0.5
        @test GH.rpo_hypr_refinement_sampling_density_m(SM.RPOPSOConfig(), 0.2) == 0.2
        @test GH.rpo_hypr_refinement_sampling_density_m(SM.rpo_pso_config(; safe_distance_m=0.3)) == 0.3
        @test GH.rpo_hypr_refinement_sampling_density_m(SM.RPOPSOConfig(refinement_sample_ds_m=0.04), 0.0) == 0.04

        @test_throws ArgumentError SM.rpo_pso_config(; n_particles=0)
        @test_throws ArgumentError SM.rpo_pso_config(; curve_type=:spline)
        @test_throws ArgumentError SM.rpo_pso_config(; retime_dt_s=0.0)
        @test_throws ArgumentError SM.rpo_pso_config(; sample_ds_m=-1.0)
        @test_throws ArgumentError SM.rpo_pso_config(; retime_max_speed_mps=0.1, retime_min_speed_mps=0.2)
    end

    @testset "RPO path sampling" begin
        line = [0.0 2.0; 0.0 0.0; 0.0 0.0]
        geom = _probe_origin_geometry(; keepout_radius_m=0.2)

        @test GH.rpo_path_length(line) ≈ 2.0
        @test GH.rpo_bezier_point(line, 0.5) ≈ [1.0, 0.0, 0.0]

        bez = GH.rpo_sample_path_bezier(line, 0.5)
        @test size(bez, 2) == 5
        @test bez[:, 1] ≈ line[:, 1]
        @test bez[:, end] ≈ line[:, 2]
        @test all(abs.(bez[2, :]) .< 1.0e-12)
        @test issorted(bez[1, :])
        # Degenerate one-point input passes through unchanged.
        single = reshape([1.0, 2.0, 3.0], 3, 1)
        @test GH.rpo_sample_path_bezier(single, 0.5) == single

        elbow = [0.0 1.0 1.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
        poly = GH.rpo_sample_path_polyline(elbow, 0.5)
        @test poly[:, 1] ≈ elbow[:, 1]
        @test poly[:, end] ≈ elbow[:, end]
        # Half arc length of the elbow is the interior knot.
        resampled = GH.rpo_resample_polyline_points(elbow, 3)
        @test resampled[:, 2] ≈ elbow[:, 2]
        @test size(GH.rpo_resample_polyline_points(elbow, 1), 2) == 2  # clamps to at least 2 samples

        @test GH.rpo_sample_path(line, 0.5; curve_type=:bezier) ≈ bez
        @test GH.rpo_sample_path(elbow, 0.5; curve_type=:polyline) ≈ poly
        @test_throws ArgumentError GH.rpo_sample_path(line, 0.5; curve_type=:spline)

        @test GH.rpo_inflated_obstacle_radius_m(geom, 0.1) ≈ 0.2 + 0.01 + 0.1
        @test GH.rpo_inflated_obstacle_radius_m(geom, -0.5) ≈ 0.2 + 0.01

        cfg_off = SM.RPOPSOConfig(adaptive_sampling_enable=false, sample_ds_m=0.5, curve_type=:polyline)
        @test GH.rpo_adaptive_sampling_min_ds_m(0.5, geom, cfg_off) == 0.5
        cfg_on = SM.RPOPSOConfig(
            adaptive_sampling_enable=true,
            sample_ds_m=0.5,
            adaptive_sampling_max_ds_m=0.5,
            adaptive_sampling_safe_distance_fraction=0.5,
            adaptive_sampling_obstacle_guard_fraction=0.5,
            curve_type=:polyline,
        )
        min_ds = GH.rpo_adaptive_sampling_min_ds_m(0.5, geom, cfg_on; safe_distance_m=0.1)
        @test min_ds ≈ min(0.5, 0.5 * 0.1, 0.5 * GH.rpo_inflated_obstacle_radius_m(geom, 0.1))

        @test GH.rpo_adaptive_sampling_step_m(0.0, 0.05, 0.5, 1.0, 1.0) ≈ 0.05
        @test GH.rpo_adaptive_sampling_step_m(10.0, 0.05, 0.5, 1.0, 1.0) ≈ 0.5
        near = GH.rpo_adaptive_sampling_step_m(0.02, 0.05, 0.5, 1.0, 1.0)
        @test near <= 0.02 + 0.05 + 1.0e-12

        seg = GH.rpo_adaptive_segment_samples(
            SVector{3, Float64}(1.0, 1.0, 0.0),
            SVector{3, Float64}(1.0, 1.0, 0.0),
            geom;
            min_ds_m=0.05,
            max_ds_m=0.5,
            far_clearance_m=1.0,
        )
        @test size(seg) == (3, 1)

        # cfg dispatcher: uniform when disabled, adaptive when enabled.
        uniform = GH.rpo_sample_path(elbow, cfg_off, geom)
        @test uniform ≈ GH.rpo_sample_path_polyline(elbow, 0.5)
        adaptive_poly = GH.rpo_sample_path(elbow, cfg_on, geom; safe_distance_m=0.05)
        @test adaptive_poly[:, 1] ≈ elbow[:, 1]
        @test adaptive_poly[:, end] ≈ elbow[:, end]
        @test size(adaptive_poly, 2) >= size(uniform, 2)
        adaptive_bez = GH.rpo_sample_path(line, cfg_on, geom; safe_distance_m=0.05, curve_type=:bezier)
        @test adaptive_bez[:, 1] ≈ line[:, 1]
        @test adaptive_bez[:, end] ≈ line[:, 2]
        @test_throws ArgumentError GH.rpo_sample_path(line, cfg_on, geom; curve_type=:spline)

        work = similar(line)
        point = zeros(3)
        GH.rpo_bezier_point!(point, work, line, 0.0)
        speed = GH.rpo_bezier_speed_estimate(line, work, point, 0.0)
        @test speed ≈ 2.0 rtol = 1.0e-3
        GH.rpo_bezier_point!(point, work, line, 1.0)
        @test GH.rpo_bezier_speed_estimate(line, work, point, 1.0) ≈ 2.0 rtol = 1.0e-3
    end

    @testset "RPO path retiming" begin
        pts = [0.0 3.0 3.0; 0.0 0.0 4.0; 0.0 0.0 0.0]
        s = GH.rpo_arc_length_params(pts)
        @test s ≈ [0.0, 3.0, 7.0]

        clean = GH.rpo_remove_near_duplicate_samples(pts)
        @test clean == pts
        dup = [0.0 0.0 1.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        deduped = @test_logs (:warn, r"near-duplicate") GH.rpo_remove_near_duplicate_samples(dup)
        @test size(deduped, 2) == 2
        @test GH.rpo_remove_near_duplicate_samples(reshape([1.0, 0.0, 0.0], 3, 1)) == reshape([1.0, 0.0, 0.0], 3, 1)

        @test GH.rpo_interpolate_along_path(reshape([1.0, 2.0, 3.0], 3, 1), [0.0], 0.5) == [1.0, 2.0, 3.0]
        @test GH.rpo_interpolate_along_path(pts, s, -1.0) == pts[:, 1]
        @test GH.rpo_interpolate_along_path(pts, s, 99.0) == pts[:, end]
        @test GH.rpo_interpolate_along_path(pts, s, 1.5) ≈ [1.5, 0.0, 0.0]
        @test GH.rpo_interpolate_along_path(pts, s, 5.0) ≈ [3.0, 2.0, 0.0]

        @test GH.rpo_curvature_from_samples(pts[:, 1:2], s[1:2]) == [0.0, 0.0]
        straight = [0.0 1.0 2.0 3.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
        κ_line = GH.rpo_curvature_from_samples(straight, GH.rpo_arc_length_params(straight))
        @test all(κ_line .< 1.0e-9)
        θ = 0.0:0.05:(pi / 2)
        circle = vcat((2.0 .* cos.(θ))', (2.0 .* sin.(θ))', zeros(length(θ))')
        κ_circle = GH.rpo_curvature_from_samples(circle, GH.rpo_arc_length_params(circle))
        @test all(abs.(κ_circle[5:(end - 5)] .- 0.5) .< 0.02)
        # Zero-length inner segment zeroes local curvature estimates.
        κ_degenerate = GH.rpo_curvature_from_samples(pts, [0.0, 0.0, 1.0])
        @test κ_degenerate == zeros(3)
        # Out-and-back path makes the central-difference tangent vanish (rpn fallback).
        outback = [0.0 1.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        @test GH.rpo_curvature_from_samples(outback, [0.0, 1.0, 2.0]) == zeros(3)

        # Reference construction from retimed samples.
        cfg = _probe_tiny_cfg(; sample_ds_m=0.5)
        far_geom = _probe_far_geometry()
        path = [0.0 4.0; 0.0 0.0; 0.0 0.0]
        t_ref, r_ref, v_ref = GH.rpo_reference_from_path(path, far_geom, cfg)
        @test length(t_ref) == size(r_ref, 2) == size(v_ref, 2)
        @test t_ref[1] == 0.0
        @test all(diff(t_ref) .≈ cfg.retime_dt_s)
        @test v_ref[:, end] ≈ v_ref[:, end - 1]
    end

    @testset "RPO replanning spheres and config" begin
        sphere = GH.RPOReplanningSphere((1.0, 2.0, 3.0), 0.5)
        @test sphere.center_rtn == SVector{3, Float64}(1.0, 2.0, 3.0)
        @test sphere.radius_m == 0.5
        @test sphere.appear_time_s == 0.0
        @test sphere.disappear_time_s == Inf
        @test sphere.label == "dynamic_sphere"
        @test_throws ArgumentError GH.RPOReplanningSphere((0.0, 0.0, 0.0), 0.0)
        @test_throws ArgumentError GH.RPOReplanningSphere((0.0, 0.0, 0.0), 1.0; appear_time_s=5.0, disappear_time_s=1.0)

        # NamedTuple normalization with primary and alias field names.
        s_primary = GH._rpo_replanning_sphere((center_rtn=(1.0, 0.0, 0.0), radius_m=0.2, label="p"))
        @test s_primary.center_rtn == SVector{3, Float64}(1.0, 0.0, 0.0)
        @test s_primary.label == "p"
        s_alias = GH._rpo_replanning_sphere((center=(2.0, 0.0, 0.0), radius=0.3, appear_time=1.0, disappear_time=9.0, velocity=(0.1, 0.0, 0.0)))
        @test s_alias.radius_m == 0.3
        @test s_alias.appear_time_s == 1.0
        @test s_alias.disappear_time_s == 9.0
        @test s_alias.velocity_rtn_mps == SVector{3, Float64}(0.1, 0.0, 0.0)
        @test GH._rpo_replanning_sphere(s_alias) === s_alias
        @test_throws ArgumentError GH._rpo_replanning_sphere((radius_m=1.0,))
        @test_throws ArgumentError GH._rpo_replanning_sphere((center=(0.0, 0.0, 0.0),))

        config = GH.RPOReplanningConfig(; spheres=[(center=(1.0, 0.0, 0.0), radius=0.5)], desired_goal_rtn=(2.0, 0.0, 0.0))
        @test config.enabled
        @test length(config.spheres) == 1
        @test config.desired_goal_rtn == SVector{3, Float64}(2.0, 0.0, 0.0)
        @test config.tracking_error_retime_m == Inf
        # Negative tracking thresholds clamp to zero without erroring.
        @test GH.RPOReplanningConfig(; tracking_error_retime_m=-3.0).tracking_error_retime_m == 0.0

        @test_throws ArgumentError GH.RPOReplanningConfig(; goal_change_tolerance_m=-1.0)
        @test_throws ArgumentError GH.RPOReplanningConfig(; safe_distance_m=-1.0)
        @test_throws ArgumentError GH.RPOReplanningConfig(; safe_distance_m=0.5, retime_clearance_m=0.1)
        @test_throws ArgumentError GH.RPOReplanningConfig(; min_replan_interval_s=-1.0)
        @test_throws ArgumentError GH.RPOReplanningConfig(; hysteresis_samples=0)
        @test_throws ArgumentError GH.RPOReplanningConfig(; sphere_surface_samples=6)
        @test_throws ArgumentError GH.RPOReplanningConfig(; remaining_sample_ds_m=0.0)
        @test_throws ArgumentError GH.RPOReplanningConfig(; tracking_error_retime_m=NaN)
        @test_throws ArgumentError GH.RPOReplanningConfig(; tracking_error_replan_m=NaN)

        moving = GH.RPOReplanningSphere((0.0, 0.0, 0.0), 1.0; appear_time_s=5.0, disappear_time_s=10.0, velocity_rtn_mps=(1.0, 0.0, 0.0))
        @test GH.rpo_replanning_sphere_center(moving, 2.0) == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test GH.rpo_replanning_sphere_center(moving, 7.0) == SVector{3, Float64}(2.0, 0.0, 0.0)

        window_cfg = GH.RPOReplanningConfig(; spheres=[moving])
        @test isempty(GH.rpo_active_replanning_spheres(window_cfg, 4.0))
        active = GH.rpo_active_replanning_spheres(window_cfg, 6.0)
        @test length(active) == 1
        @test active[1].center_rtn == SVector{3, Float64}(1.0, 0.0, 0.0)
        @test isempty(GH.rpo_active_replanning_spheres(window_cfg, 11.0))

        s1 = GH.RPOReplanningSphere((1.0, 0.0, 0.0), 0.5; label="a")
        s2 = GH.RPOReplanningSphere((0.0, 1.0, 0.0), 0.4; label="b")
        @test GH.rpo_replanning_signature(GH.RPOReplanningSphere[]) == UInt(0)
        @test GH.rpo_replanning_signature([s1, s2]) == GH.rpo_replanning_signature([s2, s1])
        @test GH.rpo_replanning_signature([s1]) != GH.rpo_replanning_signature([s2])

        pts = GH.rpo_sphere_surface_points(s1; n_points=24)
        @test size(pts) == (3, 24)
        radii = [norm(pts[:, j] - s1.center_rtn) for j in 1:24]
        @test all(abs.(radii .- 0.5) .< 1.0e-9)
    end

    @testset "RPO replanning geometry and decisions" begin
        geom = _probe_far_geometry()
        s1 = GH.RPOReplanningSphere((0.0, 1.0, 0.0), 0.3)
        @test GH.rpo_geometry_with_replanning_spheres(geom, GH.RPOReplanningSphere[]) === geom
        augmented = GH.rpo_geometry_with_replanning_spheres(geom, [s1]; sphere_surface_samples=16)
        @test size(augmented.station.points_body, 2) == 1 + 16
        @test endswith(augmented.station.name, "+replanning")
        @test augmented.station.keepout_radius_m == geom.station.keepout_radius_m
        @test augmented.chaser === geom.chaser

        plan = _probe_straight_plan(-2.0, 2.0)
        empty_plan = GH.RPOPlan(valid=true)

        remaining_empty = GH.rpo_remaining_reference_path(empty_plan, (1.0, 2.0, 3.0))
        @test remaining_empty == reshape([1.0, 2.0, 3.0], 3, 1)
        remaining = GH.rpo_remaining_reference_path(plan, (0.05, 0.0, 0.0); sample_ds_m=0.25)
        @test remaining[:, 1] ≈ [0.05, 0.0, 0.0]
        @test remaining[:, end] ≈ [2.0, 0.0, 0.0]
        @test size(remaining, 2) >= 2

        @test GH.rpo_reference_tracking_error(empty_plan, (0.0, 0.0, 0.0)) == Inf
        @test GH.rpo_reference_tracking_error(plan, (0.0, 0.7, 0.0)) ≈ 0.7
        @test GH.rpo_reference_tracking_error(plan, (0.0, 0.0, 0.0)) ≈ 0.0 atol = 1.0e-12

        current = SVector{3, Float64}(-2.0, 0.0, 0.0)

        disabled = GH.rpo_replanning_decision(plan, current, geom, GH.RPOReplanningConfig(; enabled=false), 0.0)
        @test disabled.action == :none && disabled.reason == :disabled
        invalid = GH.rpo_replanning_decision(GH.RPOPlan(), current, geom, GH.RPOReplanningConfig(), 0.0)
        @test invalid.action == :none && invalid.reason == :disabled

        goal_cfg = GH.RPOReplanningConfig(; desired_goal_rtn=(5.0, 0.0, 0.0))
        goal_dec = GH.rpo_replanning_decision(plan, current, geom, goal_cfg, 0.0)
        @test goal_dec.action == :replan && goal_dec.reason == :goal_changed
        # An empty reference treats the current position as the active goal.
        empty_goal_dec = GH.rpo_replanning_decision(empty_plan, current, geom, goal_cfg, 0.0)
        @test empty_goal_dec.action == :replan && empty_goal_dec.reason == :goal_changed

        none_dec = GH.rpo_replanning_decision(plan, current, geom, GH.RPOReplanningConfig(), 0.0)
        @test none_dec.action == :none && none_dec.reason == :clearance_ok

        track_cfg = GH.RPOReplanningConfig(; tracking_error_retime_m=1.0, tracking_error_replan_m=10.0)
        offside = SVector{3, Float64}(0.0, 3.0, 0.0)
        track_retime = GH.rpo_replanning_decision(plan, offside, geom, track_cfg, 0.0)
        @test track_retime.action == :retime && track_retime.reason == :tracking_error
        far_off = SVector{3, Float64}(0.0, 30.0, 0.0)
        track_replan = GH.rpo_replanning_decision(plan, far_off, geom, track_cfg, 0.0)
        @test track_replan.action == :replan && track_replan.reason == :tracking_error

        base_sphere_kwargs = (safe_distance_m=0.5, retime_clearance_m=1.0, sphere_surface_samples=96, remaining_sample_ds_m=0.25)
        unsafe_cfg = GH.RPOReplanningConfig(; spheres=[(center=(0.0, 0.5, 0.0), radius=0.3)], base_sphere_kwargs...)
        unsafe_dec = GH.rpo_replanning_decision(plan, current, geom, unsafe_cfg, 0.0)
        @test unsafe_dec.action == :replan && unsafe_dec.reason == :unsafe_remaining_path
        @test unsafe_dec.min_clearance < 0.5
        @test length(unsafe_dec.spheres) == 1
        @test unsafe_dec.geometry !== geom

        retime_cfg = GH.RPOReplanningConfig(; spheres=[(center=(0.0, 1.0, 0.0), radius=0.3)], base_sphere_kwargs...)
        retime_dec = GH.rpo_replanning_decision(plan, current, geom, retime_cfg, 0.0)
        @test retime_dec.action == :retime && retime_dec.reason == :low_clearance_remaining_path
        @test 0.5 <= retime_dec.min_clearance < 1.0

        ok_cfg = GH.RPOReplanningConfig(; spheres=[(center=(0.0, 5.0, 0.0), radius=0.3)], base_sphere_kwargs...)
        ok_dec = GH.rpo_replanning_decision(plan, current, geom, ok_cfg, 0.0)
        @test ok_dec.action == :none && ok_dec.reason == :clearance_ok
        @test ok_dec.min_clearance > 1.0
    end

    @testset "RPO plan-from-path and retime-existing" begin
        geom = _probe_far_geometry()
        cfg = _probe_tiny_cfg(; sample_ds_m=0.5)
        path = [-2.0 2.0; 0.0 0.0; 0.0 0.0]

        plan = GH.rpo_plan_from_path(path, geom, cfg, 0.05, 7.0; cost=1.5, diagnostics=(tag=:probe,))
        @test plan.valid
        @test plan.cost == 1.5
        @test plan.diagnostics.planned_at_s == 7.0
        @test plan.diagnostics.tag == :probe
        @test plan.path_rtn == path
        @test plan.r_ref_rtn[:, 1] ≈ path[:, 1] atol = 1.0e-9
        @test plan.r_ref_rtn[:, end] ≈ path[:, 2] atol = 1.0e-9
        @test plan.t_ref_s[1] == 0.0

        current = SVector{3, Float64}(-1.4, 0.1, 0.0)
        retimed = GH.rpo_retime_existing_plan(plan, current, geom, cfg, 0.05, 9.0)
        @test retimed.valid
        @test retimed.cost == plan.cost
        @test retimed.diagnostics.replanning_action == :retime
        # Retiming refreshes planned_at_s; inherited diagnostics cannot override it.
        @test retimed.diagnostics.planned_at_s == 9.0
        @test retimed.r_ref_rtn[:, 1] ≈ collect(current) atol = 1.0e-9
        @test retimed.r_ref_rtn[:, end] ≈ path[:, 2] atol = 1.0e-9
    end
end
