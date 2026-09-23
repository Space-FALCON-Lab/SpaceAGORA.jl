using Test
using LinearAlgebra
using Random
using StaticArrays
using SpaceAGORA

const HM_SM = SpaceAGORA.SimulationModel
const HM_GH = HM_SM.GuidanceHooks

# A station reduced to one point far from the test paths, so clearance never limits them.
function _hm_far_geometry(point=(0.0, 0.0, 500.0); keepout=0.0)
    station = HM_SM.RPOStationGeometry(reshape(collect(Float64, point), 3, 1); keepout_radius_m=keepout)
    return HM_SM.RPOReferenceGeometry(station; chaser=HM_SM.RPOCubeSatGeometry(dims_m=(0.02, 0.02, 0.02)))
end

function _hm_manuscript_cfg(; kwargs...)
    return HM_SM.rpo_pso_config(HM_SM.RPOPSOConfig(;
        hypr_mode=:manuscript,
        retime_accel_limit_enable=true,
        retime_dt_s=0.1,
        retime_a_max_mps2=0.00625,
        retime_speed_scale=0.5,
        sample_ds_m=1.0,
        adaptive_sampling_max_ds_m=4.0,
        adaptive_sampling_far_clearance_m=20.0,
        obstacle_sigmoid_k=1.0e5,
        mean_motion_radps=0.00113,
        mass_kg=5.2,
        kwargs...,
    ))
end

@testset "RPO HyPR manuscript HCW fuel proxy" begin
    n = 0.00113
    dt = 0.1
    mass, isp, g0 = 5.2, 60.0, 9.80665
    # Constant velocity along T at a fixed radial and cross-track offset: every
    # finite-difference acceleration is zero, so each u_k is the HCW term
    # [-3n^2 x - 2n v_y, 0, n^2 z].
    x0, z0, v = 3.0, -2.0, 0.2
    M = 101
    r = zeros(3, M)
    for k in 1:M
        r[:, k] .= (x0, 5.0 + v * (k - 1) * dt, z0)
    end
    proxy = HM_GH.rpo_hcw_fuel_proxy(r, dt, n, mass, isp, g0)
    u_norm = sqrt((3.0 * n^2 * x0 + 2.0 * n * v)^2 + (n^2 * z0)^2)
    @test proxy.delta_v_eq_mps ≈ (M - 2) * u_norm * dt rtol = 1.0e-12
    @test proxy.J_fuel ≈ mass * proxy.delta_v_eq_mps / (isp * g0) rtol = 1.0e-14
    # Without orbital motion the same reference needs no control.
    @test HM_GH.rpo_hcw_fuel_proxy(r, dt, 0.0, mass, isp, g0).delta_v_eq_mps ≈ 0.0 atol = 1.0e-9
    # Uniform acceleration along R without orbital motion: u_k = a on every step.
    a = 0.004
    ra = zeros(3, M)
    for k in 1:M
        ra[1, k] = 0.5 * a * ((k - 1) * dt)^2
    end
    @test HM_GH.rpo_hcw_fuel_proxy(ra, dt, 0.0, mass, isp, g0).delta_v_eq_mps ≈ (M - 2) * a * dt rtol = 1.0e-9
end

@testset "RPO HyPR manuscript Eq. 6 threshold" begin
    @test HM_GH.rpo_obstacle_sigmoid_threshold(0.10, 0.01, :manuscript) ≈ 0.11
    @test HM_GH.rpo_obstacle_sigmoid_threshold(0.10, 0.01, :legacy) ≈ 0.09
    geom = _hm_far_geometry((0.0, 0.0, 0.0))
    # Clearance 0.105 m lies between d_safe = 0.10 and d_safe + tol = 0.11: Eq. 6
    # penalizes it fully, the legacy threshold (d_safe - tol) not at all.
    half = maximum(geom.chaser.half_extents_body)
    band = reshape([0.105 + half, 0.0, 0.0], 3, 1)
    manuscript = HM_GH.rpo_clearance_stats_from_samples(band, geom, 0.10;
        obstacle_sigmoid_k=1.0e5, obstacle_sigmoid_tol_m=0.01, threshold_mode=:manuscript)
    legacy = HM_GH.rpo_clearance_stats_from_samples(band, geom, 0.10;
        obstacle_sigmoid_k=1.0e5, obstacle_sigmoid_tol_m=0.01)
    @test manuscript.obstacle_score > 1.0 - 1.0e-9
    @test legacy.obstacle_score < 1.0e-9
    clear = reshape([0.115 + half, 0.0, 0.0], 3, 1)
    @test HM_GH.rpo_clearance_stats_from_samples(clear, geom, 0.10;
        obstacle_sigmoid_k=1.0e5, obstacle_sigmoid_tol_m=0.01, threshold_mode=:manuscript).obstacle_score < 1.0e-9

    # Eq. 5 has no length term: the total is w_obs J_obs + w_fuel J_fuel.
    cfg = _hm_manuscript_cfg(w_obs=1.0e6, w_fuel=2.0, obstacle_sigmoid_tol_m=0.5, safe_distance_m=1.0)
    path = [0.0 20.0 40.0; -60.0 0.0 60.0; 0.0 0.0 0.0]
    comps = HM_GH.rpo_normalized_path_cost_components(path, _hm_far_geometry(), cfg; safe_distance_m=1.0)
    @test comps.J_obs < 1.0e-12
    @test comps.J_fuel > 0.0
    @test comps.total ≈ cfg.w_obs * comps.J_obs + cfg.w_fuel * comps.J_fuel rtol = 1.0e-14
    @test comps.violation_count == 0
    # The same path through a station point is penalized sample by sample.
    blocked = _hm_far_geometry((20.0, 0.0, 0.0); keepout=1.0)
    blocked_comps = HM_GH.rpo_normalized_path_cost_components(path, blocked, cfg; safe_distance_m=1.0)
    @test blocked_comps.violation_count > 0
    @test blocked_comps.J_obs > 1.0
end

@testset "RPO HyPR retimer boundary speeds and acceleration limit" begin
    a_max = 0.00625
    cfg = _hm_manuscript_cfg(safe_distance_m=1.0)
    geom = _hm_far_geometry()
    path = [0.0 10.0 20.0 30.0; -100.0 -40.0 40.0 100.0; 0.0 0.0 0.0 0.0]
    ref = HM_GH.rpo_retimed_reference(path, geom, cfg; safe_distance_m=1.0)
    @test ref.speed_mps[1] == 0.0
    @test ref.speed_mps[end] == 0.0
    @test ref.v_rtn[:, 1] == zeros(3)
    @test ref.v_rtn[:, end] == zeros(3)
    @test ref.r_rtn[:, 1] ≈ path[:, 1] atol = 1.0e-9
    @test ref.r_rtn[:, end] == path[:, end]
    @test ref.t_s[2] - ref.t_s[1] ≈ cfg.retime_dt_s
    @test maximum(abs.(ref.profile.a_seg)) <= a_max * (1.0 + 1.0e-9)
    @test all(ref.profile.v[2:(end - 1)] .<= ref.profile.v_point[2:(end - 1)] .+ 1.0e-12)
    # Speed changes between fixed steps stay within the limit.
    dv = diff(ref.speed_mps) ./ diff(ref.t_s)
    @test maximum(abs.(dv)) <= a_max * (1.0 + 1.0e-6)
    # The discrete reference itself: second differences of the positions stay
    # near the tangential limit (the path is almost straight here).
    dt = cfg.retime_dt_s
    acc = [norm(ref.r_rtn[:, k + 1] - 2.0 .* ref.r_rtn[:, k] + ref.r_rtn[:, k - 1]) / dt^2 for k in 2:(size(ref.r_rtn, 2) - 1)]
    @test maximum(acc) <= 1.05 * a_max
    # The velocity reference is the speed along the curve.
    @test maximum(abs.([norm(ref.v_rtn[:, k]) for k in eachindex(ref.speed_mps)] .- ref.speed_mps)) < 1.0e-12

    # A nonzero initial speed is honoured and the reference still ends at rest.
    moving = HM_SM.rpo_pso_config(cfg; retime_initial_speed_mps=0.05)
    ref_moving = HM_GH.rpo_retimed_reference(path, geom, moving; safe_distance_m=1.0)
    @test ref_moving.speed_mps[1] ≈ 0.05
    @test ref_moving.speed_mps[end] == 0.0

    # The reference builder returns the same profile and rests at the goal.
    t_ref, r_ref, v_ref = HM_SM.rpo_reference_from_path(path, geom, cfg; safe_distance_m=1.0)
    @test t_ref == ref.t_s
    @test r_ref == ref.r_rtn
    @test v_ref[:, end] == zeros(3)

    # The streamed fuel proxy equals the proxy of the stored reference held at the goal.
    streamed = HM_GH.rpo_profile_hcw_fuel_proxy(ref.profile, dt, cfg.mean_motion_radps, cfg.mass_kg, cfg.isp_s, cfg.g0_mps2)
    stored = HM_GH.rpo_hcw_fuel_proxy(hcat(ref.r_rtn, ref.r_rtn[:, end]), dt, cfg.mean_motion_radps, cfg.mass_kg, cfg.isp_s, cfg.g0_mps2)
    @test streamed.delta_v_eq_mps ≈ stored.delta_v_eq_mps rtol = 1.0e-12

    # Acceleration demand: tangential part within the limit, HCW part as computed from the state.
    demand = HM_GH.rpo_reference_accel_demand(ref, cfg.mean_motion_radps)
    @test maximum(demand.tangential_mps2) <= a_max * (1.0 + 1.0e-9)
    @test size(demand.u_rtn) == size(ref.r_rtn)

    # Pointwise limits: near a station point the manuscript's d_avail = c - d_safe slows the reference.
    near = _hm_far_geometry((15.0, 0.0, 6.0); keepout=1.0)
    ref_near = HM_GH.rpo_retimed_reference(path, near, cfg; safe_distance_m=1.0)
    @test maximum(ref_near.speed_mps) <= maximum(ref.speed_mps)
    @test ref_near.t_s[end] > ref.t_s[end]
end

@testset "RPO HyPR manuscript exploration score and coefficients" begin
    s = HM_GH.rpo_manuscript_exploration_score(200.0, 100.0, 100)
    @test s.detour_score == 1.0
    @test s.search_effort_score ≈ 1.0 - exp(-1.0)
    @test s.eta ≈ (2.0 + (1.0 - exp(-1.0))) / 3.0
    @test HM_GH.rpo_manuscript_exploration_score(100.0, 100.0, 0).eta == 0.0
    @test HM_GH.rpo_manuscript_exploration_score(500.0, 100.0, 0).detour_score == 1.0

    base = HM_SM.rpo_pso_config(HM_SM.RPOPSOConfig(
        hypr_mode=:manuscript,
        adaptive_n_waypoints_min=3, adaptive_n_waypoints_max=8,
        adaptive_n_particles_min=60, adaptive_n_particles_max=160,
        adaptive_n_iters_min=10, adaptive_n_iters_max=60,
        adaptive_w_inertia_min=0.4, adaptive_w_inertia_max=0.75,
        adaptive_c1_min=1.2, adaptive_c1_max=1.8,
        adaptive_c2_min=1.2, adaptive_c2_max=2.2,
    ))
    start = SVector(0.0, 0.0, 0.0)
    goal = SVector(0.0, 100.0, 0.0)
    warm(length_m, iterations) = (path_found=true, path_length_m=length_m, iterations=iterations)
    easy, easy_diag = HM_GH.rpo_manuscript_adaptive_pso_config(base, start, goal, warm(100.0, 0))
    hard, hard_diag = HM_GH.rpo_manuscript_adaptive_pso_config(base, start, goal, warm(250.0, 10_000))
    @test easy_diag.eta == 0.0
    @test hard_diag.eta ≈ 1.0
    @test (easy.w_inertia, easy.c1, easy.c2) == (0.4, 1.2, 2.2)
    @test hard.w_inertia ≈ 0.75
    @test hard.c1 ≈ 1.8          # equations: c1 grows with eta
    @test hard.c2 ≈ 1.2          # equations: c2 shrinks with eta
    @test (easy.n_waypoints, easy.n_particles, easy.n_iters) == (3, 60, 10)
    @test (hard.n_waypoints, hard.n_particles, hard.n_iters) == (8, 160, 60)
    mid, mid_diag = HM_GH.rpo_manuscript_adaptive_pso_config(base, start, goal, warm(125.0, 0))
    @test mid_diag.detour_score ≈ 0.25
    @test mid_diag.eta ≈ 0.5 / 3.0
    @test mid.n_particles == round(Int, (1 - mid_diag.eta) * 60 + mid_diag.eta * 160)
    # No warm-start path counts as the largest detour.
    _, failed = HM_GH.rpo_manuscript_adaptive_pso_config(base, start, goal, (path_found=false, path_length_m=0.0, iterations=50))
    @test failed.detour_score == 1.0

    @test_throws ArgumentError HM_SM.rpo_pso_config(HM_SM.RPOPSOConfig(); hypr_mode=:paper)
    @test_throws ArgumentError HM_SM.rpo_pso_config(HM_SM.RPOPSOConfig(); mean_motion_radps=-1.0)
    @test HM_SM.RPOPSOConfig().hypr_mode === :legacy
    @test HM_SM.RPOPSOConfig().retime_accel_limit_enable === false
    configured = HM_SM.RPOPSOConfig(HM_SM.RPOPSOConfigurator(
        hypr_mode=:manuscript,
        objective=HM_SM.RPOPSOObjectiveSettings(mean_motion_radps=0.001, fuel_proxy_dt_s=0.5),
        retiming=HM_SM.RPOPSORetimingSettings(accel_limit_enable=true, initial_speed_mps=0.01),
        adaptive=HM_SM.RPOPSOAdaptiveSettings(search_effort_scale=50.0),
    ))
    @test configured.hypr_mode === :manuscript
    @test configured.mean_motion_radps == 0.001
    @test HM_GH.rpo_fuel_proxy_dt_s(configured) == 0.5
    @test configured.retime_accel_limit_enable === true
    @test configured.retime_initial_speed_mps == 0.01
    @test configured.adaptive_search_effort_scale == 50.0
end

@testset "RPO HyPR manuscript planner" begin
    rng = MersenneTwister(3)
    points = zeros(3, 1500)
    for j in 1:size(points, 2)
        d = randn(rng, 3)
        points[:, j] .= 8.0 .* d ./ norm(d)
    end
    geom = HM_SM.RPOReferenceGeometry(HM_SM.RPOStationGeometry(points; keepout_radius_m=1.0);
        chaser=HM_SM.RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.3)))
    cfg = _hm_manuscript_cfg(
        w_obs=1.0e6,
        w_fuel=1.0,
        obstacle_sigmoid_tol_m=0.5,
        adaptive_sampling_max_ds_m=3.0,
        adaptive_sampling_far_clearance_m=10.0,
        rrt_warmstart_enable=true,
        rrt_warmstart_iters=1000,
        rrt_warmstart_step_size_m=2.0,
        rrt_warmstart_collision_sample_ds_m=0.5,
        rrt_warmstart_shortcut_iters=40,
        rrt_warmstart_box_margin_m=10.0,
        search_margin_m=25.0,
        adaptive_n_waypoints_min=3, adaptive_n_waypoints_max=5,
        adaptive_n_particles_min=8, adaptive_n_particles_max=16,
        adaptive_n_iters_min=3, adaptive_n_iters_max=6,
        schedule_transition_fraction=0.5,
        early_stopping_enable=false,
    )
    start = SVector(0.0, -40.0, 0.0)
    goal = SVector(0.0, 40.0, 0.0)
    plan = HM_SM.rpo_pso_plan_path(start, goal, geom, cfg; safe_distance_m=1.0, rng=MersenneTwister(741))
    @test plan.warmstart.attempted
    @test plan.warmstart.path_found
    @test plan.warmstart.iterations > 0          # the straight segment is blocked
    @test plan.warmstart.path_length_m >= norm(goal - start)
    @test plan.adaptive.mode === :manuscript
    @test 0.0 <= plan.adaptive.eta <= 1.0
    @test plan.config.n_waypoints == plan.adaptive.n_waypoints
    @test plan.config.n_particles == plan.adaptive.n_particles
    @test length(plan.cost_history) == plan.config.n_iters
    @test issorted(plan.cost_history; rev=true)
    @test plan.components.violation_count == 0
    @test plan.cost ≈ cfg.w_obs * plan.components.J_obs + cfg.w_fuel * plan.components.J_fuel rtol = 1.0e-12
    @test plan.path[:, 1] == collect(start)
    @test plan.path[:, end] == collect(goal)
    again = HM_SM.rpo_pso_plan_path(start, goal, geom, cfg; safe_distance_m=1.0, rng=MersenneTwister(741))
    @test again.path == plan.path
    @test again.cost == plan.cost
end

@testset "RPO MPC command log default" begin
    @test HM_SM.RPOMPCControlModel().command_log === nothing
    log = HM_SM.RPOControlCommandLog()
    @test isempty(log.t_s) && isempty(log.thruster_forces_n)
end

@testset "RPO MPC reference preview on update times" begin
    # A reference sampled every 0.1 s whose x coordinate is the sample index.
    n = 400
    plan = HM_SM.RPOPlan(valid=true, t_ref_s=collect(0:(n - 1)) .* 0.1,
        r_ref_rtn=vcat(reshape(collect(0.0:(n - 1)), 1, n), zeros(2, n)), v_ref_rtn=zeros(3, n))
    # Update times as the simulation produces them, including ones where t/dt
    # falls just below an integer in floating point.
    for k in (3, 14, 92, 162, 164, 250)
        t = k * 0.1
        preview = HM_SM.ControlHooks.rpo_ref_preview(plan, t, 0.1, 3)
        @test preview[1, 1] == k
        @test preview[1, 2] == k + 1
    end
    @test 16.2 / 0.1 < 162      # the floating-point case the preview must absorb
    @test HM_SM.ControlHooks.rpo_ref_preview(plan, 1.0e6, 0.1, 2)[1, end] == n - 1   # held at the end
end
