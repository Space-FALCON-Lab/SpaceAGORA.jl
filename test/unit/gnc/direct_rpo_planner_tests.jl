# The baseline and planning-budget policy must work without the simulator or HYPR.
module DirectRPOPlannerTests
using Test, Random, LinearAlgebra
include(joinpath(@__DIR__, "../../../src/gnc/interfaces/rpo_planner.jl"))
include(joinpath(@__DIR__, "../../../src/gnc/direct_rpo/direct_rpo_planner.jl"))
const P=RPOPlannerInterfaces
const D=DirectRPOPlanning
function request(; kwargs...)
    defaults=(request_id=1,chaser_id=101,target_id=201,epoch="fixture",time_s=0.,
        x_rtn=(3.,0.,0.,0.,0.,0.),target_state_ii=(7e6,0.,0.,0.,7500.,0.),goal_rtn_m=(5.,0.,0.),
        geometry=(center=[0.,0.,0.],radius=0.4),geometry_revision="sphere-v1",
        constraints=P.RPOPlanningConstraints(clearance_m=.2,max_speed_mps=.5,max_acceleration_mps2=.1),
        reference_dt_s=.1,preview_horizon_steps=2,valid_until_s=100.)
    P.RPOPlanningRequest(; merge(defaults,(;kwargs...))...)
end
function segment_clearance(a,b,g)
    x,y=collect(a),collect(b);d=y-x
    u=norm(d)==0 ? 0. : clamp(dot(g.center-x,d)/dot(d,d),0.,1.)
    norm(x+u*d-g.center)-g.radius
end
clearance(p,g)=norm(collect(p)-g.center)-g.radius
plan(req;kwargs...)=P.plan_rpo!(nothing,D.DirectRPOPlanner(;segment_clearance_at=segment_clearance,kwargs...),req,MersenneTwister(1))
validate(req,res)=P.validate_rpo_result(req,res;clearance_at=clearance)

@testset "Prospective reserve and precision refusal" begin
    @test !isdefined(@__MODULE__,:SpaceAGORA)
    @test !isdefined(P,:RPOPSOConfig)
    for x in (0.,-1.,1.,Inf,NaN)
        @test_throws ArgumentError P.RPOPlanningHeadroom(fraction=x)
    end
    for R in (3.,30.,300.,3000.), dt in (.1,.01,.001)
        req=request(x_rtn=(R,0.,0.,0.,0.,0.),goal_rtn_m=(R+2.,0.,0.),reference_dt_s=dt)
        out=plan(req); @test out.status===:candidate
        checked=validate(req,out); @test checked.accepted
        @test out.diagnostics.planning_budget.fraction==.01
        @test checked.metrics.max_speed_mps <= .5
        @test checked.metrics.max_acceleration_mps2 <= .1
        @test out.diagnostics.analytic_max_speed_mps <= .495
        @test out.diagnostics.analytic_max_acceleration_mps2 <= .099
        @test out.reference.v_ref_rtn_mps[:,end]==zeros(3)
    end
    huge=request(x_rtn=(1e16,0.,0.,0.,0.,0.),goal_rtn_m=(nextfloat(1e16),0.,0.))
    @test plan(huge).termination===:insufficient_reference_precision
    @test plan(request();headroom=P.RPOPlanningHeadroom(fraction=eps(Float64))).termination===:insufficient_reference_precision
    # Large velocity scale with a tiny acceleration budget must also refuse.
    b=P.rpo_planning_budget(request(),P.RPOPlanningHeadroom();position_scale_m=5.,velocity_scale_mps=1e16)
    @test !b.supported
    @test request().constraints.max_speed_mps==.5
end

@testset "Headroom for finite differences at different initial speeds" begin
    for dt in (.1,.01,.001), v0 in (0.,.5,5.)
        times=[i*dt for i in 0:199]; acc=.099
        pos=hcat(([3+v0*t+acc*t^2/2,0.,0.] for t in times)...)
        vel=hcat(([v0+acc*t,0.,0.] for t in times)...)
        req=request(reference_dt_s=dt,goal_rtn_m=Tuple(pos[:,end]),
            constraints=P.RPOPlanningConstraints(clearance_m=.2,max_acceleration_mps2=.1),
            validation=P.RPOValidationSettings(limit_roundoff_rtol=0.,clearance_sample_ds_m=1.))
        budget=P.rpo_planning_budget(req,P.RPOPlanningHeadroom();
            position_scale_m=maximum(norm,eachcol(pos)),velocity_scale_mps=maximum(norm,eachcol(vel)))
        @test budget.supported
        ref=P.RPOReference(t_ref_s=times,r_ref_rtn_m=pos,v_ref_rtn_mps=vel,origin_time_s=0.,valid_until_s=100.,
            chaser_id=101,target_id=201,geometry_revision="sphere-v1")
        result=P.RPOPlanningResult(request_id=1,status=:candidate,termination=:completed,reference=ref)
        @test validate(req,result).accepted
    end
end

@testset "Direct segment, lifetime, allocation and callback failures" begin
    @test_throws ArgumentError D.DirectRPOPlanner(max_reference_samples=1)
    req=request();planner=D.DirectRPOPlanner(segment_clearance_at=segment_clearance)
    @test P.initialize_planner(planner,nothing)===nothing
    @test !P.planner_capabilities(planner).restart
    @test !P.planner_capabilities(planner).retiming
    @test P.plan_rpo!(nothing,D.DirectRPOPlanner(),req,MersenneTwister(1)).termination===:segment_clearance_not_supplied
    @test plan(request(reason=:retime)).status===:unsupported
    @test plan(request(constraints=P.RPOPlanningConstraints(clearance_m=.2))).termination===:speed_limit_required
    blocked=request(x_rtn=(-3.,0.,0.,0.,0.,0.),goal_rtn_m=(3.,0.,0.))
    @test plan(blocked).termination===:direct_segment_blocked # endpoints clear, segment crosses sphere
    @test plan(request(valid_until_s=.5)).termination===:insufficient_reference_lifetime
    @test plan(req;max_reference_samples=2).termination===:reference_work_limit
    stationary=request(goal_rtn_m=(3.,0.,0.))
    @test validate(stationary,plan(stationary)).accepted
    nanplanner=D.DirectRPOPlanner(segment_clearance_at=(a,b,g)->NaN)
    @test P.plan_rpo!(nothing,nanplanner,req,MersenneTwister(1)).termination===:nonfinite_segment_clearance
    bad=D.DirectRPOPlanner(segment_clearance_at=(a,b,g)->error("query failure"))
    @test_throws ErrorException P.plan_rpo!(nothing,bad,req,MersenneTwister(1))
    rng=MersenneTwister(123);control=copy(rng)
    out=P.plan_rpo!(nothing,planner,req,rng)
    @test rand(rng)==rand(control)
    @test req.geometry.center==zeros(3)
    req.geometry.center[1]=10
    @test out.reference.r_ref_rtn_m[:,1]==[3.,0.,0.]
end

module ExternalPlanner
using ..RPOPlannerInterfaces
const P=RPOPlannerInterfaces
struct Hold <: P.AbstractRPOPlanner end
P.planner_capabilities(::Hold)=P.RPOPlannerCapabilities(state_sources=(:truth,),frames=(:target_rtn,))
function P.plan_rpo!(::Nothing,::Hold,req::P.RPOPlanningRequest,rng::P.AbstractRNG)
    ref=P.RPOReference(t_ref_s=[0.,req.reference_dt_s],r_ref_rtn_m=hcat(collect(req.goal_rtn_m),collect(req.goal_rtn_m)),
        v_ref_rtn_mps=zeros(3,2),origin_time_s=req.time_s,valid_until_s=req.valid_until_s,
        chaser_id=req.chaser_id,target_id=req.target_id,geometry_revision=req.geometry_revision)
    P.RPOPlanningResult(request_id=req.request_id,status=:candidate,termination=:completed,reference=ref)
end
end
@testset "External implementation shares the contract" begin
    req=request(goal_rtn_m=(3.,0.,0.));other=ExternalPlanner.Hold()
    @test P.validate_rpo_capabilities(other,req).accepted
    @test validate(req,P.plan_rpo!(nothing,other,req,MersenneTwister(1))).accepted
    @test !P.validate_rpo_capabilities(other,req;checkpoint=true).accepted
end
end
