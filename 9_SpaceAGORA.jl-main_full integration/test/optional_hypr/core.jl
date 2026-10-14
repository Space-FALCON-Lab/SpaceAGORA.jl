using Test, Pkg, Random, LinearAlgebra, StaticArrays, SpaceAGORA
const S = SpaceAGORA.SimulationModel
const G = S.GuidanceHooks
const P = SpaceAGORA.RPOPlannerInterfaces
@testset "Independent core installation has no HYPR execution" begin
    @test Base.find_package("SpaceAGORAHYPR") === nothing
    @test all(p.name != "SpaceAGORAHYPR" for p in values(Pkg.dependencies()))
    @test !hypr_available()
    @test all(id.name != "SpaceAGORAHYPR" for id in keys(Base.loaded_modules))
    @test RPOPSOConfig === G.RPOPSOConfig
    @test RobotArmHYPRConfig === S.RobotArmPlanning.RobotArmHYPRConfig
    @test parentmodule(HYPRRPOPlanner) === SpaceAGORA.HYPRRPOPlanning
    @test all(m.module === G for m in methods(G.rpo_pso_plan_path))
    planner = HYPRRPOPlanner(RPOPSOConfig())
    @test_throws HYPRUnavailableError P.initialize_planner(planner, nothing)
    @test_throws HYPRUnavailableError make_rpo_configuration(planner=planner)
    @test_throws HYPRUnavailableError S.rpo_pso_plan_path(zeros(3),ones(3),nothing,RPOPSOConfig())
    @test_throws HYPRUnavailableError S.rpo_reference_from_path(zeros(3,2),nothing,RPOPSOConfig())
    @test_throws HYPRUnavailableError robot_arm_hypr_path_cost_components(zeros(3,2),nothing,nothing,RobotArmSphereObstacle[],RobotArmHYPRConfig())
    error_text = try P.initialize_planner(planner,nothing); "" catch e; sprint(showerror,e) end
    @test occursin("using SpaceAGORAHYPR",error_text)
    @test_throws HYPRUnavailableError SpaceAGORA.SimulationLifecycle.preflight_guidance(
        S.RPOGuidanceModel(), nothing; isolate_state=true)
    mktempdir() do folder
        destination=joinpath(folder,"refused-run")
        args=make_rpo_configuration(planner=DirectRPOPlanner(),mission_time_s=.2,
            simulation_settings=SimulationSettings(results=true,results_directory=destination,
                verbose=false,generate_plots=false))
        owner=only(args.guidance_model.guidance_effectors)
        owner.planner=planner
        @test_throws HYPRUnavailableError run_simulation(args)
        @test !ispath(destination)
        @test owner.runtime===nothing
    end
end
@testset "Configured comparison boundaries explain missing HYPR" begin
    a=SVector(-2.,0.,0.); b=SVector(2.,0.,0.); path=hcat(a,b)
    cfg=RPOPSOConfig(n_particles=4,n_iters=1)
    comparison=G.RPOPlannerComparisonConfig(pso_config=cfg)
    rng=MersenneTwister(20261003); untouched=copy(rng)
    for planner in (:hypr,:pso_unrefined,:rrt_connect,:rrt_connect_bezier,:rrt_star,:chomp,:stomp)
        @test_throws HYPRUnavailableError G.rpo_plan_comparison_path(planner,a,b,nothing,comparison;rng)
    end
    @test rand(rng,4)==rand(untouched,4)
    # Guard before collecting cases, consuming an iterator, or creating a controller.
    @test_throws HYPRUnavailableError G.rpo_run_planner_comparison_batch(
        (error("cases must not be consumed") for _ in 1:1),nothing,comparison)
    @test_throws HYPRUnavailableError G.rpo_track_retimed_path_lqmpc(path,b,nothing,cfg,G.RPOLQMPCTrackingSettings())
    for f in (G.rpo_retime_path,G.rpo_retimed_reference,G.rpo_path_cost)
        @test_throws HYPRUnavailableError f(path,nothing,cfg)
    end
    for f in (G.rpo_rrt_connect_plan_path,G.rpo_rrt_connect_bezier_plan_path,G.rpo_rrt_star_plan_path)
        @test_throws HYPRUnavailableError f(a,b,nothing,cfg)
    end
    for f in (G.rpo_sample_path,G.rpo_sample_path_with_params)
        @test_throws HYPRUnavailableError f(path,cfg,nothing)
    end
    for f in (G.rpo_sample_path_polyline_adaptive,G.rpo_sample_path_bezier_adaptive,G.rpo_sample_path_bezier_adaptive_with_params)
        @test_throws HYPRUnavailableError f(path,nothing,cfg)
    end
    # Those broad compatibility fallbacks must not hide the shared scalar overload.
    @test G.rpo_sample_path(path,.1;curve_type=:polyline)==G.rpo_sample_path_polyline(path,.1)
    @test !hypr_available()
end
@testset "Shared search, metrics and retiming stay usable" begin
    a=SVector(-2.,0.,0.); b=SVector(2.,0.,0.)
    for (planner,settings) in ((G.rpo_rrt_connect_plan_path,G.RPORRTConnectSettings(n_iters=5)),
            (G.rpo_rrt_star_plan_path,G.RPORRTStarSettings(n_iters=5)))
        result=planner(a,b,nothing; bounds=(SVector(-3.,-3.,-3.),SVector(3.,3.,3.)),settings=settings,
            evaluate_components=p->(total=G.rpo_path_length(p),),evaluate_cost=G.rpo_path_length,
            edge_is_safe=(a,b)->true,rng=MersenneTwister(741))
        @test result.path==hcat(a,b)
        @test result.cost==4.
        @test result.path_found
    end
    geo=S.RPOReferenceGeometry(S.RPOStationGeometry(reshape([0.,0.,50.],3,1);keepout_radius_m=.25);
        chaser=S.RPOCubeSatGeometry(dims_m=(.1,.1,.3)))
    path=hcat(a,b)
    samples=G.rpo_sample_path_polyline(path,.1)
    @test G.rpo_path_length(samples)≈4.
    ref,_,_=G.rpo_retime_samples(samples,geo;max_speed_mps=.5,min_speed_mps=0.,dt_s=.1,max_steps=1000,
        available_distance=(c,d,s)->d,pointwise_speed=(d,k)->.5)
    @test isapprox(ref[:,1],a;atol=1e-12,rtol=0) && isapprox(ref[:,end],b;atol=1e-12,rtol=0)
    @test !hypr_available()
end
@testset "Baseline simulation and quintic arm run without HYPR" begin
    args=make_rpo_configuration(planner=DirectRPOPlanner(),mission_time_s=.2)
    sol=run_simulation(args;return_solution=true)
    @test sol.t[end]≈.2
    @test rpo_run_report(sol).planners[1].request_id==1
    model=default_cloth_arm_model();base=ClothArmBasePose([0.,0.,0.]);q0=[0.,.3,-.3]
    target=cloth_fk(model,base,[.6,-.4,.4]).end_effector_position
    out=plan_robot_arm_motion(model,base,q0,target;config=RobotArmPlannerConfig(duration_s=.5))
    @test robot_arm_sample_hypr_path(hcat(q0,q0),3)==repeat(q0,1,3)
    @test robot_arm_clearance_stats_from_samples(model,base,hcat(q0,q0),RobotArmSphereObstacle[],0.).violation_count==0
    @test out.planner===:cloth_quintic
    @test all(isfinite,out.q_ref)
    @test_throws HYPRUnavailableError plan_robot_arm_motion(model,base,q0,target;planner=:hypr)
    @test !hypr_available()
end
include(joinpath(@__DIR__,"..","unit","gnc","rpo_planner_contract_tests.jl"))
include(joinpath(@__DIR__,"..","unit","gnc","direct_rpo_planner_tests.jl"))
include(joinpath(@__DIR__,"..","unit","gnc","rrt_boundary_tests.jl"))
include(joinpath(@__DIR__,"..","unit","gnc","shared_retiming_tests.jl"))
