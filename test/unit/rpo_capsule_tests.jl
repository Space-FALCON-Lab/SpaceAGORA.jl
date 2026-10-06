using Test, Random, LinearAlgebra, StaticArrays, SpaceAGORA
const CapsuleNav = SpaceAGORA.SimulationModel.NavigationHooks
const CapsuleGuidance = SpaceAGORA.SimulationModel.GuidanceHooks

@testset "Continuous capsule clearance" begin
    N = CapsuleNav
    chaser = N.RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.3))
    radius = norm(chaser.half_extents_body)
    geometry = N.RPOReferenceGeometry(N.RPOStationGeometry(zeros(3,1); keepout_radius_m=0.25); chaser=chaser)
    a, b = SVector(-1.,0.,0.), SVector(1.,0.,0.)
    @test N.rpo_clearance_distance_to_station(a,geometry) > 0.25
    @test N.rpo_clearance_distance_to_station(b,geometry) > 0.25
    @test N.rpo_capsule_clearance_to_station(a,b,geometry) ≈ -0.25-radius
    @test N.rpo_capsule_clearance_to_station(a,a,geometry) ≈ N.rpo_clearance_distance_to_station(a,geometry)
    @test N.rpo_capsule_clearance_to_station(b,a,geometry) ≈ N.rpo_capsule_clearance_to_station(a,b,geometry)
    y = 0.25+radius+0.245
    @test N.rpo_capsule_clearance_to_station(SVector(-1.,y,0.),SVector(1.,y,0.),geometry) ≈ 0.245
    @test radius ≈ sqrt(0.0275)
    rng = MersenneTwister(213)
    points = randn(rng,3,100)
    geometry = N.RPOReferenceGeometry(N.RPOStationGeometry(points;keepout_radius_m=0.25);chaser=chaser)
    for _ in 1:100
        a,b = SVector{3}(randn(rng,3)),SVector{3}(randn(rng,3))
        ab = b-a
        expected = minimum(norm(points[:,i]-(a+clamp(dot(points[:,i]-a,ab)/dot(ab,ab),0.,1.)*ab)) for i in axes(points,2))-0.25-radius
        @test N.rpo_capsule_clearance_to_station(a,b,geometry) ≈ expected atol=1e-12
    end
end

@testset "Retiming rejects an interior capsule intersection" begin
    N,G = CapsuleNav,CapsuleGuidance
    geometry = N.RPOReferenceGeometry(N.RPOStationGeometry(zeros(3,1);keepout_radius_m=0.25))
    cfg = G.RPOPSOConfig(curve_type=:polyline,sample_ds_m=10.0,adaptive_sampling_enable=false,clearance_feasibility_tol_m=0.005)
    path = [-1.0 1.0; 0.0 0.0; 0.0 0.0]
    result = G.rpo_prepare_retimed_candidate(path,geometry,cfg;safe_distance_m=0.25)
    @test result.profile === nothing
    @test result.components.violation_count > 0
    @test result.components.min_clearance < 0.0
    @test result.components.J_obs > 0.99
    @test G.RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M == 0.25
    path[2,:] .= 2.0
    result = G.rpo_prepare_retimed_candidate(path,geometry,cfg;safe_distance_m=0.25)
    @test result.profile !== nothing
    @test result.components.violation_count == 0
end

@testset "Shared planner collision model" begin
    N,G = CapsuleNav,CapsuleGuidance
    geometry = N.RPOReferenceGeometry(N.RPOStationGeometry(zeros(3,1);keepout_radius_m=0.25))
    path = [-1.0 1.0; 0.0 0.0; 0.0 0.0]
    cfg = G.RPOPSOConfig(curve_type=:polyline,sample_ds_m=10.0,adaptive_sampling_enable=false,clearance_feasibility_tol_m=0.005)
    components = G.rpo_normalized_path_cost_components(path,geometry,cfg;safe_distance_m=0.25)
    @test components.violation_count > 0
    @test components.min_clearance ≈ N.rpo_capsule_clearance_to_station(path[:,1],path[:,2],geometry)
    for settings in (G.RPORRTConnectSettings(clearance_feasibility_tol_m=0.005),G.RPORRTStarSettings(clearance_feasibility_tol_m=0.005))
        @test !G.rpo_rrt_segment_is_safe(path[:,1],path[:,2],geometry,settings;safe_distance_m=0.25)
        y = 0.25 + norm(geometry.chaser.half_extents_body) + 0.246
        @test G.rpo_rrt_segment_is_safe([-1.,y,0.],[1.,y,0.],geometry,settings;safe_distance_m=0.25)
        @test !G.rpo_rrt_segment_is_safe([-1.,y-0.002,0.],[1.,y-0.002,0.],geometry,settings;safe_distance_m=0.25)
    end
    @test G.rpo_soft_obstacle_cost_from_samples(path,geometry;safe_distance_m=0.25,obstacle_margin_m=0.5) > 0
    waypoints = [-1.0 1.0 2.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    @test G.rpo_stomp_waypoint_state_cost(waypoints,1,geometry;safe_distance_m=0.25,w_obs=1.0,w_len=0.0,w_smooth=0.0,obstacle_margin_m=0.5) > 0
end
