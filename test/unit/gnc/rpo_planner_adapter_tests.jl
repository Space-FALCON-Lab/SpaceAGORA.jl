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
        @test P.planner_capabilities(planner).retiming # retime_rpo! (PR 19b)
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

@testset "Lifetime refusals agree across retiming policies" begin
    for mode in (:legacy, :manuscript), accel in (false, true)
        planner=H.HYPRRPOPlanner(config(mode;retime_accel_limit_enable=accel))
        constraints=P.RPOPlanningConstraints(clearance_m=.2,max_speed_mps=.5,
            max_acceleration_mps2=accel ? .1 : nothing)
        req=request(time_s=17.,constraints=constraints)
        seed=742
        accepted=P.plan_rpo!(nothing,planner,req,MersenneTwister(seed))
        @test accepted.status===:candidate
        endpoint=req.time_s+last(accepted.reference.t_ref_s)
        # An exact lifetime fit is legal in either mode, at a nonzero epoch.
        exact=request(time_s=17.,constraints=constraints,valid_until_s=endpoint)
        fits=P.plan_rpo!(nothing,planner,exact,MersenneTwister(seed))
        @test fits.status===:candidate
        @test validate(exact,fits).accepted
        @test fits.reference.t_ref_s==accepted.reference.t_ref_s
        # The legacy retimer used to report these as failed/reference_rejected.
        # Cover an ordinary shortage and a one-ULP shortage, with no tolerance
        # expansion in the adapter's prospective lifetime budget.
        for expiry in (endpoint-req.reference_dt_s, prevfloat(endpoint))
            short=request(time_s=17.,constraints=constraints,valid_until_s=expiry)
            result=P.plan_rpo!(nothing,planner,short,MersenneTwister(seed))
            @test result.status===:infeasible
            @test result.termination===:insufficient_reference_lifetime
            @test result.reference===nothing
        end
    end
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

@testset "Retiming an active reference" begin
    # One second in, the chaser is on the reference or a few centimetres off it.
    for (mode,accel,offset) in ((:legacy,false,.03),(:legacy,true,0.),(:manuscript,true,.03))
        constraints=P.RPOPlanningConstraints(clearance_m=.2,max_speed_mps=.5,
            max_acceleration_mps2=accel ? .1 : nothing)
        planner=H.HYPRRPOPlanner(config(mode;retime_accel_limit_enable=accel))
        @test P.validate_rpo_capabilities(planner,request(reason=:retime)).accepted
        active=P.plan_rpo!(nothing,planner,request(constraints=constraints),MersenneTwister(742)).reference
        x=active.r_ref_rtn_m[:,min(11,end)].+[0.,offset,0.]
        req=request(request_id=2,reason=:retime,time_s=1.,x_rtn=(x...,0.,0.,0.),
            constraints=constraints,valid_until_s=101.)
        rng=MersenneTwister(5);untouched=copy(rng)
        out=P.retime_rpo!(nothing,planner,deepcopy(req),deepcopy(active),rng)
        @test out.status===:candidate && out.termination===:completed
        @test validate(req,out).accepted
        @test rand(rng)==rand(untouched) # retiming draws no random numbers
        ref=out.reference
        @test ref.r_ref_rtn_m[:,1]==x && ref.r_ref_rtn_m[:,end]==[5.,1.,0.]
        @test (ref.origin_time_s,ref.valid_until_s,ref.geometry_revision)==(1.,101.,"fixture-v1")
        d=out.diagnostics
        path=H.retiming_route(active.r_ref_rtn_m,x)
        @test d.path==path && d.revision===:retime && d.path_kind===:polyline
        @test d.input_config.curve_type===:polyline && d.input_config.retime_max_speed_mps==.495
        # Matched input: the shared retiming kernels give the same reference.
        t,r,v=G.rpo_reference_from_path(path,req.geometry,d.input_config;safe_distance_m=.2)
        @test (ref.t_ref_s,ref.r_ref_rtn_m,ref.v_ref_rtn_mps)==(t,r,v)
        # The general name reaches it through the family.
        again=SpaceAGORA.propose_revision!(nothing,planner,deepcopy(req),deepcopy(active),MersenneTwister(5))
        @test again.reference.r_ref_rtn_m==ref.r_ref_rtn_m
    end
    planner=H.HYPRRPOPlanner(config())
    active=P.plan_rpo!(nothing,planner,request(),MersenneTwister(742)).reference
    # The legacy tangential retimer can turn the corner at an off-reference start
    # faster than the common acceleration limit; the unchanged validator refuses it.
    x=active.r_ref_rtn_m[:,11].+[0.,.01,0.]
    req=request(request_id=2,reason=:retime,time_s=1.,x_rtn=(x...,0.,0.,0.),valid_until_s=101.)
    turned=P.retime_rpo!(nothing,planner,req,active,MersenneTwister(1))
    @test (turned.status,turned.termination)==(:failed,:reference_rejected)
    @test turned.diagnostics.validation.reason===:acceleration_limit
    retime(;kwargs...)=P.retime_rpo!(nothing,planner,request(;reason=:retime,kwargs...),active,MersenneTwister(1))
    @test retime().status===:candidate
    changed=retime(goal_rtn_m=(5.,2.,0.))
    @test (changed.status,changed.termination)==(:infeasible,:goal_not_on_reference)
    @test retime(chaser_id=102).termination===:reference_mismatch
    @test retime(geometry=nothing).termination===:unsupported_geometry
    short=retime(valid_until_s=1.)
    @test (short.status,short.termination)==(:infeasible,:insufficient_reference_lifetime)
    @test retime(constraints=P.RPOPlanningConstraints(clearance_m=.2,max_speed_mps=.5,
        max_acceleration_mps2=.1)).termination===:completed
end

@testset "Replaying a retained path" begin
    for (mode,accel) in ((:legacy,false),(:legacy,true),(:manuscript,true))
        constraints=P.RPOPlanningConstraints(clearance_m=.2,max_speed_mps=.5,
            max_acceleration_mps2=accel ? .1 : nothing)
        req=request(constraints=constraints)
        planned=P.plan_rpo!(nothing,H.HYPRRPOPlanner(config(mode;retime_accel_limit_enable=accel)),req,MersenneTwister(742))
        @test planned.status===:candidate
        d=planned.diagnostics
        # The optimizer's path and configuration, replayed for the same request, give
        # the same reference without searching or drawing.
        retained=H.RetainedPathRPOPlanner(d.path,d.effective_config)
        rng=MersenneTwister(5);untouched=copy(rng)
        out=P.plan_rpo!(nothing,retained,deepcopy(req),rng)
        @test out.status===:candidate && out.termination===:completed && validate(req,out).accepted
        @test rand(rng)==rand(untouched)
        @test (out.reference.t_ref_s,out.reference.r_ref_rtn_m,out.reference.v_ref_rtn_mps)==
            (planned.reference.t_ref_s,planned.reference.r_ref_rtn_m,planned.reference.v_ref_rtn_mps)
        @test out.diagnostics.planner===:retained_path && out.diagnostics.path==d.path
        @test SpaceAGORA.propose!(nothing,retained,deepcopy(req),MersenneTwister(5)).reference.r_ref_rtn_m==
            out.reference.r_ref_rtn_m
    end
    cfg=config();path=[3. 4. 5.;0. 0. 1.;0. .5 0.]
    retained=H.RetainedPathRPOPlanner(path,cfg)
    @test !SpaceAGORA.consumes_rng(retained)
    caps=P.planner_capabilities(retained)
    @test caps.state_sources==(:truth,) && caps.frames==(:target_rtn,) && !caps.retiming && !caps.restart
    @test P.validate_rpo_capabilities(retained,request(reason=:retime)).reason===:unsupported_retiming
    # Copies at construction: later edits to the caller's arrays change nothing.
    path[1,2]=99.;@test retained.path_rtn[1,2]==4.
    replay(r;planner=retained)=P.plan_rpo!(nothing,planner,r,MersenneTwister(1))
    # End points within the tolerance are set to the request's; beyond it, refused.
    near=replay(request(x_rtn=(3.0+5e-7,0.,0.,0.,0.,0.)))
    @test near.status===:candidate && near.reference.r_ref_rtn_m[:,1]==[3.0+5e-7,0.,0.]
    @test near.reference.r_ref_rtn_m[:,end]==[5.,1.,0.]
    @test replay(request(x_rtn=(3.1,0.,0.,0.,0.,0.))).termination===:start_not_on_path
    @test replay(request(goal_rtn_m=(5.,1.1,0.))).termination===:goal_not_on_path
    @test replay(request(x_rtn=(3.1,0.,0.,0.,0.,0.));planner=H.RetainedPathRPOPlanner(retained.path_rtn,cfg;
        endpoint_tolerance_m=.2)).status===:candidate
    @test replay(request(reason=:retime)).termination===:retiming_not_implemented
    @test replay(request(geometry=nothing)).termination===:unsupported_geometry
    @test replay(request(valid_until_s=1.)).termination===:insufficient_reference_lifetime
    # The same request mapping as HYPRRPOPlanner: a requested acceleration limit needs
    # the acceleration-limited retimer.
    @test replay(request();planner=H.RetainedPathRPOPlanner(retained.path_rtn,config(retime_accel_limit_enable=false))).termination===
        :acceleration_limited_retimer_required
    for bad in (zeros(2,3),zeros(3,1),[0. NaN;0. 0.;0. 0.])
        @test_throws ArgumentError H.RetainedPathRPOPlanner(bad,cfg)
    end
    @test_throws ArgumentError H.RetainedPathRPOPlanner(retained.path_rtn,cfg;max_reference_samples=1)
    @test_throws ArgumentError H.RetainedPathRPOPlanner(retained.path_rtn,cfg;endpoint_tolerance_m=-1.)
    @test SpaceAGORA.RetainedPathRPOPlanner===H.RetainedPathRPOPlanner
end
end
