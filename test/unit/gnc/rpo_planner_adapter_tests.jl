module RPOPlannerAdapterTests
using Test, SpaceAGORA, Random, LinearAlgebra
const P=SpaceAGORA.RPOPlannerInterfaces
const H=SpaceAGORA.HYPRRPOPlanning
const D=SpaceAGORA.DirectRPOPlanning
const S=SpaceAGORA.SimulationModel
const G=S.GuidanceHooks
function geometry(;center=[0.,0.,50.])
    S.RPOReferenceGeometry(S.RPOStationGeometry(reshape(center,3,1);keepout_radius_m=.25);
        chaser=S.RPOCubeSatGeometry(dims_m=(.1,.1,.3)))
end
function request(;kwargs...)
    defaults=(request_id=1,chaser_id=101,target_id=201,epoch="fixture",time_s=0.,
        x_rtn=(3.,0.,0.,0.,0.,0.),target_state_ii=(7e6,0.,0.,0.,7500.,0.),goal_rtn_m=(5.,1.,0.),
        geometry=geometry(),geometry_revision="fixture-v1",
        constraints=P.RPOPlanningConstraints(clearance_m=.2,max_speed_mps=.5,max_acceleration_mps2=.1),
        reference_dt_s=.1,preview_horizon_steps=2,valid_until_s=100.)
    P.RPOPlanningRequest(;merge(defaults,(;kwargs...))...)
end
function config(mode=:legacy; kwargs...)
    defaults=(hypr_mode=mode,n_waypoints=2,n_particles=8,n_iters=3,
        adaptive_enable=false,adaptive_n_waypoints_min=2,adaptive_n_waypoints_max=2,
        adaptive_n_particles_min=8,adaptive_n_particles_max=8,adaptive_n_iters_min=3,adaptive_n_iters_max=3,
        sample_ds_m=.2,safe_distance_m=.2,retime_dt_s=.1,retime_a_max_mps2=.1,
        retime_max_speed_mps=.5,retime_accel_limit_enable=true,iteration_runtime_limit_s=Inf,
        rrt_warmstart_enable=true)
    S.RPOPSOConfig(;merge(defaults,(;kwargs...))...)
end
clearance(p,g)=S.rpo_clearance_distance_to_station(p,g)
function analytic_segment(a,b,g)
    x,y=collect(a),collect(b);d=y-x;den=dot(d,d)
    minimum(norm(x+(den==0 ? 0. : clamp(dot(q-x,d)/den,0.,1.))*d-q) for q in eachcol(g.station.points_body)) -
        g.station.keepout_radius_m-maximum(g.chaser.half_extents_body)
end
validate(req,out)=P.validate_rpo_result(req,out;clearance_at=clearance)

@testset "Matched-input HYPR adapter parity and unchanged validator" begin
    for (mode,accel,seed,adaptive,reason) in ((:legacy,false,741,false,:initial),(:legacy,true,742,false,:initial),
            (:manuscript,true,743,false,:forced),(:legacy,true,744,true,:replan))
        req=request(reason=reason,constraints=P.RPOPlanningConstraints(clearance_m=.2,max_speed_mps=.5,
            max_acceleration_mps2=accel ? .1 : nothing))
        cfg=config(mode;retime_accel_limit_enable=accel,adaptive_enable=adaptive)
        planner=H.HYPRRPOPlanner(cfg;rrt_on_replan=true)
        rng=MersenneTwister(seed);control=MersenneTwister(seed)
        out=P.plan_rpo!(nothing,planner,req,rng)
        @test out.status===:candidate
        @test validate(req,out).accepted
        d=out.diagnostics
        @test d.planning_budget.fraction==.01
        @test d.input_config.retime_max_speed_mps==.495
        @test planner.config.retime_max_speed_mps==.5
        @test req.constraints.max_acceleration_mps2==(accel ? .1 : nothing)
        raw=S.rpo_pso_plan_path(req.x_rtn[1:3],req.goal_rtn_m,req.geometry,d.input_config;safe_distance_m=.2,rng=control)
        t,r,v=G.rpo_reference_from_path(raw.path,req.geometry,raw.config;safe_distance_m=.2)
        @test out.reference.t_ref_s==t
        @test out.reference.r_ref_rtn_m==r
        @test out.reference.v_ref_rtn_mps==v
        @test d.path==raw.path && d.cost==raw.cost && d.cost_history==raw.cost_history
        @test d.effective_config==raw.config
        @test rand(rng)==rand(control)
        @test !P.planner_capabilities(planner).restart
        @test !P.planner_capabilities(planner).retiming
    end
    # An output from the unchanged retimer at the unchanged physical limit
    # still fails for real turning excess. Planning headroom is not validator slack.
    req=request();cfg=config();raw=S.rpo_pso_plan_path(req.x_rtn[1:3],req.goal_rtn_m,req.geometry,cfg;safe_distance_m=.2,rng=MersenneTwister(742))
    t,r,v=G.rpo_reference_from_path(raw.path,req.geometry,raw.config;safe_distance_m=.2)
    ref=P.RPOReference(t_ref_s=t,r_ref_rtn_m=r,v_ref_rtn_mps=v,origin_time_s=0.,valid_until_s=100.,chaser_id=101,target_id=201,geometry_revision="fixture-v1")
    bad=P.RPOPlanningResult(request_id=1,status=:candidate,termination=:completed,reference=ref)
    @test validate(req,bad).reason===:acceleration_limit
    insufficient=P.plan_rpo!(nothing,H.HYPRRPOPlanner(cfg;headroom=P.RPOPlanningHeadroom(fraction=1e-8)),req,MersenneTwister(742))
    @test insufficient.status===:failed
    @test insufficient.termination===:reference_rejected
    @test insufficient.diagnostics.validation.reason===:acceleration_limit
    @test insufficient.diagnostics.planning_budget.fraction==1e-8
end

@testset "Two planners, common request and validation" begin
    req=request();baseline=D.DirectRPOPlanner(segment_clearance_at=analytic_segment)
    for planner in (baseline,H.HYPRRPOPlanner(config()))
        @test P.validate_rpo_capabilities(planner,req).accepted
        out=P.plan_rpo!(P.initialize_planner(planner,nothing),planner,deepcopy(req),MersenneTwister(742))
        @test validate(req,out).accepted
        @test out.reference.chaser_id==101 && out.reference.target_id==201
        @test out.reference.frame===:target_rtn
        @test out.reference.r_ref_rtn_m[:,1]==[3.,0.,0.]
        @test out.reference.r_ref_rtn_m[:,end]==[5.,1.,0.]
    end
    near=request(geometry=geometry(center=[4.,.5,0.]))
    @test P.plan_rpo!(nothing,baseline,near,MersenneTwister(1)).termination===:direct_segment_blocked
end

@testset "Unsupported modes, budgets and predictable failure" begin
    req=request();rng=MersenneTwister(5);before=copy(rng)
    plan(r;cfg=config(),kwargs...)=P.plan_rpo!(nothing,H.HYPRRPOPlanner(cfg;kwargs...),r,rng)
    @test plan(request(reason=:retime)).termination===:retiming_not_implemented
    @test plan(request(geometry=nothing)).termination===:unsupported_geometry
    @test plan(req;cfg=config(retime_min_speed_mps=.499)).termination===:configured_speed_exceeds_planning_budget
    @test plan(req;cfg=config(retime_accel_limit_enable=false)).termination===:acceleration_limited_retimer_required
    @test plan(request(x_rtn=(1e16,0.,0.,0.,0.,0.),goal_rtn_m=(nextfloat(1e16),0.,0.))).termination===:insufficient_reference_precision
    @test rand(rng)==rand(before) # these refusals happen before the optimizer
    @test plan(req;max_reference_samples=3).termination===:reference_work_limit
    @test plan(request(valid_until_s=.5)).termination===:insufficient_reference_lifetime
    zero=plan(req;cfg=config(iteration_runtime_limit_s=0.))
    @test zero.status===:candidate && !zero.diagnostics.iteration_timed_out
    timed=plan(req;cfg=config(iteration_runtime_limit_s=1e-20))
    @test timed.status===:failed && timed.diagnostics.validation.reason===:time_budget_not_allowed
    @test_throws ArgumentError H.HYPRRPOPlanner(config();max_reference_samples=1)
end
end
