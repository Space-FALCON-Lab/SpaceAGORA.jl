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
    # Constant velocity along T at a fixed radial and cross-track offset. The
    # reference is at rest before its first sample and after its last, so the
    # first term departs from rest (v/Δt along T), the last brakes (-v/Δt), and
    # every term between is the HCW term [-3n^2 x - 2n v_y, 0, n^2 z].
    x0, z0, v = 3.0, -2.0, 0.2
    M = 101
    r = zeros(3, M)
    for k in 1:M
        r[:, k] .= (x0, 5.0 + v * (k - 1) * dt, z0)
    end
    proxy = HM_GH.rpo_hcw_fuel_proxy(r, dt, n, mass, isp, g0)
    u_norm = sqrt((3.0 * n^2 * x0 + 2.0 * n * v)^2 + (n^2 * z0)^2)
    u_departure = norm([-3.0 * n^2 * x0, v / dt, n^2 * z0])
    u_arrival = norm([-3.0 * n^2 * x0 - 2.0 * n * v, -v / dt, n^2 * z0])
    @test proxy.delta_v_eq_mps ≈ dt * (u_departure + (M - 2) * u_norm + u_arrival) rtol = 1.0e-12
    @test proxy.J_fuel ≈ mass * proxy.delta_v_eq_mps / (isp * g0) rtol = 1.0e-14
    # Without orbital motion the same reference costs its departure and its braking.
    @test HM_GH.rpo_hcw_fuel_proxy(r, dt, 0.0, mass, isp, g0).delta_v_eq_mps ≈ 2.0 * v rtol = 1.0e-9
    # Uniform acceleration from rest along R without orbital motion: a/2 to leave
    # rest, a on every interior step, and braking from the last velocity.
    a = 0.004
    ra = zeros(3, M)
    for k in 1:M
        ra[1, k] = 0.5 * a * ((k - 1) * dt)^2
    end
    @test HM_GH.rpo_hcw_fuel_proxy(ra, dt, 0.0, mass, isp, g0).delta_v_eq_mps ≈ (2M - 3) * a * dt rtol = 1.0e-9
    # A single sample needs no control.
    @test HM_GH.rpo_hcw_fuel_proxy(r[:, 1:1], dt, n, mass, isp, g0).delta_v_eq_mps == 0.0
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

@testset "RPO HyPR reference builder, fuel proxy and acceleration demand" begin
    a_max = 0.00625
    cfg = _hm_manuscript_cfg(safe_distance_m=1.0)
    geom = _hm_far_geometry()
    path = [0.0 10.0 20.0 30.0; -100.0 -40.0 40.0 100.0; 0.0 0.0 0.0 0.0]
    # The same reference HYPR's retimer builds, assembled from the services it uses.
    samples, params, clearances = HM_GH.rpo_sample_path_with_params(path, cfg, geom; safe_distance_m=1.0,
        base_ds_m=HM_GH.rpo_retime_sampling_ds_m(cfg, 1.0), curve_type=cfg.curve_type)
    profile = HM_GH.rpo_retime_profile(HM_GH.RPORetimeCurve(path, cfg.curve_type), samples, params, clearances,
        geom, cfg; safe_distance_m=1.0)
    ref = HM_GH.rpo_retimed_reference_from_profile(profile, cfg.retime_dt_s)
    dt = cfg.retime_dt_s

    # The reference builder returns the same profile and rests at the goal.
    t_ref, r_ref, v_ref = HM_SM.rpo_reference_from_path(path, geom, cfg; safe_distance_m=1.0)
    @test t_ref == ref.t_s
    @test r_ref == ref.r_rtn
    @test v_ref[:, end] == zeros(3)

    # The streamed fuel proxy equals the proxy of the stored reference.
    streamed = HM_GH.rpo_profile_hcw_fuel_proxy(ref.profile, dt, cfg.mean_motion_radps, cfg.mass_kg, cfg.isp_s, cfg.g0_mps2)
    stored = HM_GH.rpo_hcw_fuel_proxy(ref.r_rtn, dt, cfg.mean_motion_radps, cfg.mass_kg, cfg.isp_s, cfg.g0_mps2)
    @test streamed.delta_v_eq_mps ≈ stored.delta_v_eq_mps rtol = 1.0e-12

    # Acceleration demand: tangential part within the limit, HCW part as computed from the state.
    demand = HM_GH.rpo_reference_accel_demand(ref, cfg.mean_motion_radps)
    @test maximum(demand.tangential_mps2) <= a_max * (1.0 + 1.0e-9)
    @test size(demand.u_rtn) == size(ref.r_rtn)
end

@testset "RPO HyPR manuscript coefficients" begin
    @test_throws ArgumentError HM_SM.rpo_pso_config(HM_SM.RPOPSOConfig(); hypr_mode=:paper)
    @test_throws ArgumentError HM_SM.rpo_pso_config(HM_SM.RPOPSOConfig(); mean_motion_radps=-1.0)
    @test HM_SM.RPOPSOConfig().hypr_mode === :legacy
    @test HM_SM.RPOPSOConfig().retime_accel_limit_enable === false
    configured = HM_SM.RPOPSOConfig(HM_SM.RPOPSOConfigurator(
        hypr_mode=:manuscript,
        swarm=HM_SM.RPOPSOSwarmSettings(station_box_margin_m=(1.0, 2.0, 3.0)),
        objective=HM_SM.RPOPSOObjectiveSettings(mean_motion_radps=0.001, fuel_proxy_dt_s=0.5),
        retiming=HM_SM.RPOPSORetimingSettings(accel_limit_enable=true, initial_speed_mps=0.01),
        adaptive=HM_SM.RPOPSOAdaptiveSettings(search_effort_scale=50.0),
        cull=HM_SM.RPOPSOCullSettings(noise_abs_m=0.4),
    ))
    @test configured.station_box_margin_m == (1.0, 2.0, 3.0)
    @test configured.cull_noise_abs_m == 0.4
    @test HM_SM.rpo_pso_config(configured; pso_station_box_margin=(4, 5, 6)).station_box_margin_m == (4.0, 5.0, 6.0)
    @test HM_SM.rpo_pso_config(configured; pso_cull_noise_abs=0.1).cull_noise_abs_m == 0.1
    @test_throws ArgumentError HM_SM.rpo_pso_config(configured; station_box_margin_m=(1.0, -1.0, 1.0))
    @test_throws ArgumentError HM_SM.rpo_pso_config(configured; station_box_margin_m=(1.0, Inf, 1.0))
    @test_throws ArgumentError HM_SM.rpo_pso_config(configured; cull_noise_abs_m=-0.1)
    @test configured.hypr_mode === :manuscript
    @test configured.mean_motion_radps == 0.001
    @test configured.fuel_proxy_dt_s == 0.5
    @test configured.retime_accel_limit_enable === true
    @test configured.retime_initial_speed_mps == 0.01
    @test configured.adaptive_search_effort_scale == 50.0
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
