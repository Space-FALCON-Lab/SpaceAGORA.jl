using Test
using StaticArrays
using LinearAlgebra
using Random
using SpaceAGORA
const SM = SpaceAGORA.SimulationModel

@testset "RPO RTN frame helpers" begin
    r_t = SVector{3, Float64}(7000e3, 0.0, 0.0)
    v_t = SVector{3, Float64}(0.0, 7500.0, 0.0)
    r_rel = SVector{3, Float64}(10.0, -2.0, 1.0)
    v_rel = SVector{3, Float64}(0.01, 0.02, -0.01)
    r_c, v_c = SM.GuidanceHooks.rtn_to_inertial_relative_state(r_rel, v_rel, r_t, v_t)
    x = SM.GuidanceHooks.inertial_to_rtn_relative_state(r_c, v_c, r_t, v_t)
    @test x[1:3] ≈ r_rel atol=1e-9
    @test x[4:6] ≈ v_rel atol=1e-9
end

@testset "RPO geometry and PSO basics" begin
    @test isfile(SpaceAGORA.station_cad_path(:gateway))
    gateway_points = SpaceAGORA.load_rpo_station_cad_pointcloud(:gateway; n_points=16, rng=MersenneTwister(740))
    @test size(gateway_points) == (3, 16)
    points = SpaceAGORA.load_rpo_station_pointcloud(:demo)
    station = SM.RPOStationGeometry(points; keepout_radius_m=0.25)
    geom = SM.RPOReferenceGeometry(station; chaser=SM.RPOCubeSatGeometry(dims_m=(0.2, 0.2, 0.3)))
    stats = SM.rpo_clearance_to_station(SVector{3, Float64}(8.0, 0.0, 0.0), geom)
    @test stats.clearance > 0.0
    simple_station = SM.RPOStationGeometry([0.0 2.0 -3.0; 0.0 0.0 4.0; 0.0 0.0 0.0])
    nearest, distance, idx = SM.nearest_station_point(SVector{3, Float64}(1.6, 0.2, 0.0), simple_station)
    @test idx == 2
    @test nearest ≈ SVector{3, Float64}(2.0, 0.0, 0.0)
    @test distance ≈ sqrt(0.2)

    cfg = SM.RPOPSOConfig(
        n_waypoints=1,
        n_particles=8,
        n_iters=2,
        sample_ds_m=0.5,
        retime_a_max_mps2=0.1,
        tf_s=20.0,
        cost_ref_distance_m=10.0,
        curve_type=:bezier,
        adaptive_enable=false,
    )
    plan = SM.rpo_pso_plan_path(
        SVector{3, Float64}(-6.0, 0.0, 0.0),
        SVector{3, Float64}(6.0, 0.0, 0.0),
        geom,
        cfg;
        safe_distance_m=0.1,
        rng=MersenneTwister(1),
    )
    @test size(plan.path, 1) == 3
    @test isfinite(plan.cost)

    modular_cfg = SM.rpo_pso_config(SM.RPOPSOConfigurator(
        swarm=SM.RPOPSOSwarmSettings(
            n_waypoints=2,
            n_particles=9,
            n_iters=3,
            sample_ds_m=0.5,
            curve_type=:bezier,
        ),
        objective=SM.RPOPSOObjectiveSettings(
            w_len=1.2,
            w_obs=1.0e5,
            cost_ref_distance_m=10.0,
            tf_s=20.0,
        ),
        adaptive=SM.RPOPSOAdaptiveSettings(
            enabled=false,
            n_particles_min=6,
            n_particles_max=12,
            n_iters_min=2,
            n_iters_max=5,
        ),
        cull=SM.RPOPSOCullSettings(
            enabled=true,
            fraction_max=0.25,
            start_iter=1,
        ),
        schedule=SM.RPOPSOScheduleSettings(enabled=true),
        refinement=SM.RPOPSORefinementSettings(enabled=false),
        retiming=SM.RPOPSORetimingSettings(a_max_mps2=0.1),
    ))
    @test modular_cfg.n_waypoints == 2
    @test modular_cfg.n_particles == 9
    @test modular_cfg.adaptive_enable === false
    @test modular_cfg.cull_enable === true
    @test modular_cfg.schedule_enable === true
    @test modular_cfg.refinement_enable === false
    updated_cfg = SM.rpo_pso_config(modular_cfg; pso_particles=11, pso_cull_start_iter=2)
    @test updated_cfg.n_particles == 11
    @test updated_cfg.cull_start_iter == 2

    fixed_cfg, adaptive_stats = SM.GuidanceHooks.rpo_adaptive_pso_config(
        modular_cfg,
        SVector{3, Float64}(-6.0, 0.0, 0.0),
        SVector{3, Float64}(6.0, 0.0, 0.0),
        geom;
        safe_distance_m=0.1,
    )
    @test fixed_cfg.n_particles == modular_cfg.n_particles
    @test adaptive_stats.enabled === false
end

@testset "RPO six-axis thrust allocation" begin
    thrusters = SM.SixAxisThrusterModel(max_thrust_n=SVector{6, Float64}(1, 1, 1, 1, 1, 1))
    forces = SM.rpo_allocate_six_axis_thrusters(SVector{3, Float64}(0.5, -0.25, 0.1), thrusters)
    @test forces ≈ SVector{6, Float64}(0.5, 0.0, 0.0, 0.25, 0.1, 0.0)
    body_force, body_torque = SM.rpo_thruster_wrench_body(forces, thrusters)
    @test body_force ≈ SVector{3, Float64}(0.5, -0.25, 0.1)
    @test body_torque ≈ SVector{3, Float64}(0.0, 0.0, 0.0)
end
