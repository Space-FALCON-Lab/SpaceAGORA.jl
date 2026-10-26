# Execute with the optional HYPR integration environment.
module FinalISpaRoReconstructionChecks
using Test, LinearAlgebra, StaticArrays, Random, SpaceAGORA, HYPR
const SM=SpaceAGORA.SimulationModel
const G=SM.GuidanceHooks
const C=SM.ControlHooks
# Retained numerical regression fixture from the paper reconstruction at bfa32f6a.
# It supports generic actuator/replanning/retiming assertions; no swarm policy is tested here.
function paper_config()
    cfg = RPOPSOConfig(hypr_mode=:manuscript, n_waypoints=5,n_particles=100,n_iters=60,
        mass_kg=5.,isp_s=60.,w_obs=1e6,w_fuel=1.,safe_distance_m=.5,
        obstacle_sigmoid_k=1e5,obstacle_sigmoid_tol_m=.001,
        adaptive_allow_downscale=true,adaptive_n_particles_min=60,adaptive_n_particles_max=160,
        adaptive_n_iters_min=10,adaptive_n_iters_max=60,
        adaptive_n_waypoints_min=3,adaptive_n_waypoints_max=8,
        retime_dt_s=.1,fuel_proxy_dt_s=.1,retime_reaction_time_s=.25,
        retime_a_max_mps2=.00625,retime_speed_scale=.5,retime_max_speed_mps=.25,
        retime_accel_limit_enable=false,rrt_warmstart_enable=true,
        rrt_warmstart_iters=1000,rrt_warmstart_step_size_m=.75,
        rrt_warmstart_shortcut_iters=80,cull_start_iter=10,cull_fraction_max=.25,
        stagnation_learning_threshold=10,schedule_transition_fraction=.5,
        refinement_start_iter=60,refinement_period=10,refinement_sample_ds_m=.025,
        early_stopping_enable=false,iteration_runtime_limit_s=Inf)
    n=sqrt(SM.Earth().μ/(SM.Earth().Rp_e+420000.0)^3)
    return G.rpo_pso_config(cfg;mean_motion_radps=n, reexplore_max_waypoints=8)
end

cfg=paper_config()
@testset "Actuator limits" begin
    jets=SM.SixAxisThrusterModel(max_thrust_n=fill(.05,6))
    @test sum(C.rpo_allocate_six_axis_thrusters([.005,0.,0.],jets)) == .005
    @test sum(C.rpo_allocate_six_axis_thrusters([.005,0.,0.],jets;minimum_throttle_fraction=.2)) == 0
    @test sum(C.rpo_allocate_six_axis_thrusters([.01,0.,0.],jets;minimum_throttle_fraction=.2)) == 0
    @test C.rpo_allocate_six_axis_thrusters([.011,-.1,0.],jets;minimum_throttle_fraction=.2) ≈ [.011,0,0,.05,0,0]
    @test_throws ArgumentError C.rpo_allocate_six_axis_thrusters(zeros(3),jets;minimum_throttle_fraction=NaN)
end
@testset "Replanning independent triggers and reference clock" begin
    far=SM.RPOReferenceGeometry(SM.RPOStationGeometry(reshape([100.,100.,100.],3,1);keepout_radius_m=0.))
    plan=SM.RPOPlan(valid=true,t_ref_s=[0.,1.,2.],r_ref_rtn=[0. 1. 2.;0. 0. 0.;0. 0. 0.],v_ref_rtn=zeros(3,3),path_rtn=[0. 2.;0. 0.;0. 0.],diagnostics=(planned_at_s=10.,))
    pc=G.RPOReplanningConfig(safe_distance_m=.5,retime_clearance_m=.6,tracking_error_retime_m=.75,tracking_error_replan_m=1.1)
    @test G.rpo_reference_tracking_error(plan,[0.,0.,0.];t=12.) ≈ 2.
    @test G.rpo_reference_tracking_error(plan,[.5,0.,0.];t=10.5) ≈ 0.
    @test G.rpo_replanning_decision(plan,[0.,2.,0.],far,pc,10.).action==:replan
    ps=G.RPOReplanningConfig(spheres=[G.RPOReplanningSphere([1000.,1000.,1000.],1.)],safe_distance_m=.5,retime_clearance_m=.6,tracking_error_retime_m=.75,tracking_error_replan_m=1.1)
    @test G.rpo_replanning_decision(plan,[0.,2.,0.],far,ps,10.).action==:replan
    obstacle=SM.RPOReferenceGeometry(SM.RPOStationGeometry(reshape([1.,0.,0.],3,1);keepout_radius_m=0.))
    @test G.rpo_replanning_decision(plan,[0.,0.,0.],obstacle,pc,10.).action==:replan
    @test G.rpo_replanning_decision(plan,[0.,.9,0.],far,pc,10.).action==:retime
    @test G.rpo_replanning_decision(plan,[0.,0.,0.],far,pc,10.).action==:none
end
@testset "Completed reference does not reset the arrival clock" begin
    geo=SM.RPOReferenceGeometry(SM.RPOStationGeometry(reshape([.7,0.,0.],3,1);keepout_radius_m=0.))
    plan=SM.RPOPlan(valid=true,t_ref_s=[0.],r_ref_rtn=zeros(3,1),v_ref_rtn=zeros(3,1),path_rtn=zeros(3,1),diagnostics=(planned_at_s=10.,))
    pc=G.RPOReplanningConfig(safe_distance_m=.5,retime_clearance_m=.6,tracking_error_retime_m=.75,tracking_error_replan_m=1.1)
    d=G.rpo_replanning_decision(plan,zeros(3),geo,pc,11.)
    @test d.action==:none && d.reason==:reference_complete
    @test .5 <= d.min_clearance < .6
    @test G.rpo_replanning_decision(plan,[0.,.9,0.],geo,pc,11.).action==:retime
    @test G.rpo_replanning_decision(plan,[0.,2.,0.],geo,pc,11.).action==:replan
    blocked=SM.RPOReferenceGeometry(SM.RPOStationGeometry(zeros(3,1);keepout_radius_m=0.))
    @test G.rpo_replanning_decision(plan,zeros(3),blocked,pc,11.).action==:replan
    changed=G.RPOReplanningConfig(safe_distance_m=.5,retime_clearance_m=.6,desired_goal_rtn=[2.,0.,0.])
    @test G.rpo_replanning_decision(plan,zeros(3),geo,changed,11.).reason==:goal_changed
end
@testset "Remaining route retains corners" begin
    points=[0. .04 .04;0. 0. .04;0. 0. 0.]
    plan=SM.RPOPlan(valid=true,t_ref_s=[0.,1.,2.],r_ref_rtn=points,v_ref_rtn=zeros(3,3),path_rtn=points)
    remaining=G.rpo_remaining_reference_path(plan,points[:,1];sample_ds_m=.1,t=0.)
    @test any(remaining[:,j]==points[:,2] for j in axes(remaining,2))
    @test G.rpo_path_length(remaining) ≈ .08
    @test all(p[2]==0 || p[1]==.04 for p in eachcol(remaining))
    @test_throws ArgumentError G.rpo_remaining_reference_path(plan,points[:,1];sample_ds_m=0.)
end
@testset "Retiming uses current progress without recursive connectors" begin
    geo=SM.RPOReferenceGeometry(SM.RPOStationGeometry(reshape([100.,100.,100.],3,1);keepout_radius_m=0.))
    route=[0. 1. 2.;0. 0. 0.;0. 0. 0.]
    plan=SM.RPOPlan(valid=true,t_ref_s=[0.,10.,20.],r_ref_rtn=route,v_ref_rtn=zeros(3,3),path_rtn=route)
    cc=G.rpo_pso_config(paper_config();curve_type=:polyline)
    p1=G.rpo_retime_existing_plan(plan,[.5,.001,0.],geo,cc,.5,15.;current_velocity_rtn=[.04,0.,0.])
    @test p1.r_ref_rtn[:,1]==[.5,.001,0.]
    @test p1.v_ref_rtn[:,1]==[.04,0.,0.]
    @test p1.diagnostics.retiming_geometry_rtn[:,1]≈[.5,0.,0.]
    @test p1.r_ref_rtn[:,end] ≈ [2.,0.,0.] atol=1e-12 rtol=0.
    p2=G.rpo_retime_existing_plan(p1,[.6,.001,0.],geo,cc,.5,15.1;current_velocity_rtn=[.04,0.,0.])
    @test p2.diagnostics.retiming_geometry_rtn[:,1]≈[.6,0.,0.]
    @test all(iszero,p2.diagnostics.retiming_geometry_rtn[2:3,:])
    @test last(p2.t_ref_s)<last(p1.t_ref_s)
end
@testset "Retiming gives the measured-state connector its own traversal time" begin
    geo=SM.RPOReferenceGeometry(SM.RPOStationGeometry(reshape([100.,100.,100.],3,1);keepout_radius_m=0.))
    route=[0. 1. 2.;0. 0. 0.;0. 0. 0.]
    plan=SM.RPOPlan(valid=true,t_ref_s=[0.,10.,20.],r_ref_rtn=route,v_ref_rtn=zeros(3,3),path_rtn=route)
    cc=G.rpo_pso_config(paper_config();curve_type=:polyline)
    limit=cc.retime_speed_scale*cc.retime_max_speed_mps
    for offset in (0.,.001,.9)
        current=[.5,offset,0.]
        p=G.rpo_retime_existing_plan(plan,current,geo,cc,.5,15.;current_velocity_rtn=[.04,0.,0.])
        chords=diff(p.r_ref_rtn;dims=2)./reshape(diff(p.t_ref_s),1,:)
        @test maximum(norm.(eachcol(chords)))<=limit+1e-9
        @test last(p.t_ref_s)>=G.rpo_path_length(p.path_rtn)/limit-1e-9
        @test p.r_ref_rtn[:,1]≈current
        @test p.r_ref_rtn[:,end]≈route[:,end]
        @test p.v_ref_rtn[:,1]==[.04,0.,0.]
        @test p.v_ref_rtn[:,end]==zeros(3)
        @test p.v_ref_rtn[:,2:end-1]≈chords[:,2:end]
        # A new retime refers to the original projected route, never to the
        # displaced connector generated by the preceding update.
        next_current=[.6,offset/2,0.]
        q=G.rpo_retime_existing_plan(p,next_current,geo,cc,.5,15.1;current_velocity_rtn=[.04,0.,0.])
        @test q.diagnostics.retiming_geometry_rtn[:,1]≈[.6,0.,0.]
        @test all(iszero,q.diagnostics.retiming_geometry_rtn[2:3,:])
        @test maximum(norm.(eachcol(diff(q.r_ref_rtn;dims=2)./reshape(diff(q.t_ref_s),1,:))))<=limit+1e-9
    end
    corner_route=[0. 1. 1. 2.;0. 0. 1. 1.;0. 0. 0. 0.]
    corner_plan=SM.RPOPlan(valid=true,t_ref_s=[0.,10.,20.,30.],r_ref_rtn=corner_route,v_ref_rtn=zeros(3,4),path_rtn=corner_route)
    q=G.rpo_retime_existing_plan(corner_plan,[.5,-.9,0.],geo,cc,.5,1.;current_velocity_rtn=zeros(3))
    @test q.path_rtn[:,2:end]==corner_route[:,2:end]
    # The endpoint-only route also needs a timed connector after displacement.
    point=SM.RPOPlan(valid=true,t_ref_s=[0.],r_ref_rtn=reshape([2.,0.,0.],3,1),v_ref_rtn=zeros(3,1),path_rtn=reshape([2.,0.,0.],3,1))
    q=G.rpo_retime_existing_plan(point,[2.,.9,0.],geo,cc,.5,20.;current_velocity_rtn=zeros(3))
    @test q.r_ref_rtn[:,1]≈[2.,.9,0.]
    @test q.r_ref_rtn[:,end]≈[2.,0.,0.]
    @test last(q.t_ref_s)>=.9/limit-1e-9
end
@testset "Retiming is invariant to collinear time sampling" begin
    thin=[0. 2.;0. 0.;0. 0.]
    dense=vcat(reshape(collect(range(0.,2.;length=1001)),1,:),zeros(2,1001))
    @test G._rpo_retime_route_vertices(dense)==thin
    bend=[0. 1. 1. 2.;0. 0. 1. 1.;0. 0. 0. 0.]
    @test G._rpo_retime_route_vertices(bend)==bend
    reverse=[0. 1. .5 2.;0. 0. 0. 0.;0. 0. 0. 0.]
    @test G._rpo_retime_route_vertices(reverse)==reverse
    excursion=[0. 1. 0.;0. 0. 0.;0. 0. 0.]
    @test G._rpo_retime_route_vertices(excursion)==excursion
    smallbend=[0. 1. 2.;0. 1e-8 0.;0. 0. 0.]
    @test G._rpo_retime_route_vertices(smallbend)==smallbend
    geo=SM.RPOReferenceGeometry(SM.RPOStationGeometry(reshape([100.,100.,100.],3,1);keepout_radius_m=0.))
    cc=G.rpo_pso_config(paper_config();curve_type=:polyline)
    makeplan(r)=SM.RPOPlan(valid=true,t_ref_s=collect(range(0.,20.;length=size(r,2))),r_ref_rtn=r,v_ref_rtn=zeros(size(r)),path_rtn=r)
    for offset in (.001,.9)
        a=G.rpo_retime_existing_plan(makeplan(thin),[.5,offset,0.],geo,cc,.5,15.;current_velocity_rtn=[.04,0.,0.])
        b=G.rpo_retime_existing_plan(makeplan(dense),[.5,offset,0.],geo,cc,.5,15.;current_velocity_rtn=[.04,0.,0.])
        @test a.t_ref_s==b.t_ref_s
        @test a.r_ref_rtn≈b.r_ref_rtn
        @test a.v_ref_rtn≈b.v_ref_rtn
    end
end
@testset "Nonlinear propagation and swept clearance" begin
    p=G.rpo_init_two_body_plant(zeros(6),cfg.mean_motion_radps)
    @test norm(G.rpo_step_two_body!(p,zeros(3),10.)) < 1e-8
    p1=G.rpo_init_two_body_plant([10.,5.,2.,0.,0.,0.],cfg.mean_motion_radps)
    p2=deepcopy(p1)
    x1=G.rpo_step_two_body!(p1,[.001,0.,0.],10.;max_step_s=.05)
    x2=G.rpo_step_two_body!(p2,[.001,0.,0.],10.;max_step_s=.025)
    @test norm(x1-x2)<1e-5
    geo=SM.RPOReferenceGeometry(SM.RPOStationGeometry(reshape([0.,0.,0.],3,1);keepout_radius_m=0.))
    @test G.rpo_comparison_segment_clearance([-1.,0.,0.],[1.,0.,0.],geo)<0
end
@testset "Executed reference and mass accounting" begin
    geo=SM.RPOReferenceGeometry(SM.RPOStationGeometry(reshape([20.,20.,20.],3,1);keepout_radius_m=0.))
    path=[0. .5;0. 0.;0. 0.]
    tr=G.RPOLQMPCTrackingSettings(mean_motion_radps=cfg.mean_motion_radps,horizon=60,plant=:two_body,
        minimum_throttle_fraction=.2,settle_time_s=2.,q_pos=20.,q_vel=2.,qf_pos=400.,qf_vel=40.)
    run=G.rpo_track_retimed_path_lqmpc(path,path[:,end],geo,cfg,tr;safe_distance_m=.5)
    @test run.v_ref_rtn[:,1] == zeros(3)
    @test run.v_ref_rtn[:,end] == zeros(3)
    @test run.fuel_used ≈ first(run.mass_history_kg)-last(run.mass_history_kg) atol=1e-12
    @test all(f==0 || f>.01 for f in run.thruster_forces_n)
    @test all(f<=.05 for f in run.thruster_forces_n)
    @test sum(values(run.qp_status_counts))==size(run.u_hist,2)-1
end
println("CHECKS_COMPLETE");flush(stdout)

@testset "Detection retention and moving snapshot" begin
    sphere=G.RPOReplanningSphere([0.,0.,0.],1.;appear_time_s=1.,disappear_time_s=2.,velocity_rtn_mps=[1.,0.,0.])
    pc=G.RPOReplanningConfig(spheres=[sphere],observation_window_s=2.)
    @test isempty(G.rpo_active_replanning_spheres(pc,.9))
    active=only(G.rpo_active_replanning_spheres(pc,3.))
    @test active.center_rtn≈[2.,0.,0.]
    @test G.rpo_replanning_sphere_center(active,3.)≈active.center_rtn
    @test length(G.rpo_active_replanning_spheres(pc,4.))==1
    @test isempty(G.rpo_active_replanning_spheres(pc,4.0001))
    @test_throws ArgumentError G.RPOReplanningConfig(observation_window_s=-1.)
end
@testset "Canonical curved geometry is independent of the reference time grid" begin
    geo=SM.RPOReferenceGeometry(SM.RPOStationGeometry(reshape([100.,100.,100.],3,1);keepout_radius_m=0.))
    cc=G.rpo_pso_config(paper_config();curve_type=:bezier)
    controls=[0. 1. 2. 3.;0. 1. -1. 0.;0. 0. 0. 0.]
    p=G.rpo_plan_from_path(controls,geo,cc,.5,0.)
    @test haskey(p.diagnostics,:retiming_geometry_rtn)
    nominal=copy(p.diagnostics.retiming_geometry_rtn)
    @test size(nominal,2)<size(p.r_ref_rtn,2)
    @test nominal[:,1]==controls[:,1] && nominal[:,end]==controls[:,end]
    # Replace only the time grid, retaining the same nominal geometric route.
    thinidx=unique(vcat(1:7:length(p.t_ref_s),length(p.t_ref_s)))
    thin=SM.RPOPlan(valid=true,t_ref_s=p.t_ref_s[thinidx],r_ref_rtn=p.r_ref_rtn[:,thinidx],
        v_ref_rtn=p.v_ref_rtn[:,thinidx],path_rtn=p.path_rtn,diagnostics=p.diagnostics)
    x=nominal[:,3]+[0.,0.,.0002]
    a=G.rpo_retime_existing_plan(p,x,geo,cc,.5,1.;current_velocity_rtn=[.01,0.,0.])
    b=G.rpo_retime_existing_plan(thin,x,geo,cc,.5,1.;current_velocity_rtn=[.01,0.,0.])
    @test a.t_ref_s==b.t_ref_s
    @test a.r_ref_rtn==b.r_ref_rtn
    @test a.v_ref_rtn==b.v_ref_rtn
    @test p.diagnostics.retiming_geometry_rtn==nominal
    @test a.r_ref_rtn[:,1]==x
    @test a.r_ref_rtn[:,end]≈controls[:,end]
    @test all(iszero,a.diagnostics.retiming_geometry_rtn[3,:])
    # Repeated measured states remain off-curve; nominal samples must not drift.
    for k in 1:12
        x=a.r_ref_rtn[:,min(3,end)]+[0.,0.,.00001]
        a=G.rpo_retime_existing_plan(a,x,geo,cc,.5,1.0 + .1k;current_velocity_rtn=[.01,0.,0.])
        @test all(iszero,a.diagnostics.retiming_geometry_rtn[3,:])
        @test a.r_ref_rtn[:,1]==x
        @test norm(a.r_ref_rtn[:,2]-a.r_ref_rtn[:,1])/.1>.001
        @test maximum(norm.(eachcol(diff(a.r_ref_rtn;dims=2))))/.1<=.125+1e-9
    end
end
@testset "Canonical polyline retains corners shorter than the sampling step" begin
    geo=SM.RPOReferenceGeometry(SM.RPOStationGeometry(reshape([100.,100.,100.],3,1);keepout_radius_m=0.))
    corners=[0. .001 .001 .002;0. 0. .001 .001;0. 0. 0. 0.]
    for adaptive in (false,true)
        cc=G.rpo_pso_config(paper_config();curve_type=:polyline,adaptive_sampling_enable=adaptive)
        @test G._rpo_canonical_retiming_geometry(corners,geo,cc,.5)==corners
        @test G._rpo_canonical_retiming_geometry(corners[:,end:end],geo,cc,.5)==corners[:,end:end]
    end
end

@testset "Expired multipoint reference preserves callback completion clock" begin
    geo=SM.RPOReferenceGeometry(SM.RPOStationGeometry(reshape([1.,.7,0.],3,1);keepout_radius_m=0.))
    plan=SM.RPOPlan(valid=true,t_ref_s=[0.,.1],r_ref_rtn=[0. 1.;0. 0.;0. 0.],v_ref_rtn=zeros(3,2),path_rtn=[0. 1.;0. 0.;0. 0.],diagnostics=(planned_at_s=10.,))
    pc=G.RPOReplanningConfig(safe_distance_m=.5,retime_clearance_m=.6,tracking_error_retime_m=.75,tracking_error_replan_m=1.1)
    x=[1.,1e-6,0.,0.,0.,0.]
    @test G.rpo_replanning_decision(plan,x[1:3],geo,pc,10.05).action==:retime
    d=G.rpo_replanning_decision(plan,x[1:3],geo,pc,10.2)
    @test d.action==:none && d.reason==:reference_complete
    @test .5<=d.min_clearance<.6
    @test G.rpo_replanning_decision(plan,[1.,0.,.9],geo,pc,10.2).action==:retime
    @test G.rpo_replanning_decision(plan,[1.,0.,2.],geo,pc,10.2).action==:replan
    @test G.rpo_replanning_decision(plan,[1.,.1,0.],geo,pc,10.2).reason==:unsafe_remaining_path
    changed=G.RPOReplanningConfig(safe_distance_m=.5,retime_clearance_m=.6,desired_goal_rtn=[2.,0.,0.])
    @test G.rpo_replanning_decision(plan,x[1:3],geo,changed,10.2).reason==:goal_changed
    for epoch in (0.,10.,1e6)
        shifted=SM.RPOPlan(valid=true,t_ref_s=plan.t_ref_s,r_ref_rtn=plan.r_ref_rtn,v_ref_rtn=plan.v_ref_rtn,path_rtn=plan.path_rtn,diagnostics=(planned_at_s=epoch,))
        deadline=epoch+.1
        @test G.rpo_replanning_decision(shifted,x[1:3],geo,pc,prevfloat(deadline)).action==:retime
        for now in (deadline,nextfloat(deadline))
            @test G.rpo_replanning_decision(shifted,x[1:3],geo,pc,now).reason==:reference_complete
            @test size(G.rpo_future_reference_path(shifted,now),2)==1
        end
    end
end
@testset "Malformed or stale retained geometry cannot be released" begin
    geo=SM.RPOReferenceGeometry(SM.RPOStationGeometry(reshape([100.,100.,100.],3,1);keepout_radius_m=0.))
    route=[0. 1.;0. 0.;0. 0.]
    cc=G.rpo_pso_config(paper_config();curve_type=:polyline)
    pc=G.RPOReplanningConfig(safe_distance_m=.5,retime_clearance_m=.6)
    for cached in (zeros(2,2),zeros(3,0),[0. NaN;0. 0.;0. 0.],zeros(3,2),[0.,1.])
        p=SM.RPOPlan(valid=true,t_ref_s=[0.,1.],r_ref_rtn=route,v_ref_rtn=zeros(3,2),path_rtn=route,diagnostics=(retiming_geometry_rtn=cached,planned_at_s=0.))
        @test !G._rpo_retiming_geometry_valid(p)
        @test_throws ArgumentError G.rpo_retime_existing_plan(p,zeros(3),geo,cc,.5,.1)
    end
end

end
