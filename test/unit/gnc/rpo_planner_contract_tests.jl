# Deliberately standalone: this contract must load without SpaceAGORA/HYPR.
module RPOPlannerContractTests
using Test, Random, LinearAlgebra
include(joinpath(@__DIR__, "../../../src/gnc/interfaces/rpo_planner.jl"))
const P = RPOPlannerInterfaces

function request(; kwargs...)
    defaults = (request_id=1, chaser_id=101, target_id=201, epoch="2026-01-01T00:00:00Z",
        time_s=0.0, x_rtn=(3.,0.,0.,0.,0.,0.), target_state_ii=(7e6,0.,0.,0.,7500.,0.),
        goal_rtn_m=(5.,0.,0.), geometry=zeros(3,1), geometry_revision="synthetic-origin-v1",
        constraints=P.RPOPlanningConstraints(clearance_m=0.1), reference_dt_s=1.,
        preview_horizon_steps=2, valid_until_s=10.)
    return P.RPOPlanningRequest(; merge(defaults, (;kwargs...))...)
end
function reference(; kwargs...)
    defaults = (t_ref_s=[0.,1.,2.], r_ref_rtn_m=[3. 4. 5.; 0. 0. 0.; 0. 0. 0.],
        v_ref_rtn_mps=[0. 1. 0.; 0. 0. 0.; 0. 0. 0.], origin_time_s=0.,
        valid_until_s=10., chaser_id=101, target_id=201, geometry_revision="synthetic-origin-v1")
    return P.RPOReference(; merge(defaults, (;kwargs...))...)
end
candidate(ref=reference(); kwargs...) = P.RPOPlanningResult(;
    merge((request_id=1,status=:candidate,termination=:completed,reference=ref),(;kwargs...))...)
clearance(p, geometry) = norm(collect(p)-geometry[:,1]) - 0.4
validate(req, res; kwargs...) = P.validate_rpo_result(req,res;clearance_at=clearance,kwargs...)

@testset "Independent contract and snapshot ownership" begin
    @test !isdefined(P, :RPOPSOConfig)
    @test !isdefined(@__MODULE__, :SpaceAGORA)
    geom=zeros(3,1); state=[3.,0.,0.,0.,0.,0.]; epoch=[2026,1,1]
    r=request(geometry=geom,x_rtn=state,epoch=epoch)
    geom[1]=50;state[1]=99;epoch[1]=2000
    @test r.geometry==zeros(3,1)
    @test r.x_rtn[1]==3
    @test r.epoch[1]==2026
    planner_copy=deepcopy(r);planner_copy.geometry[1]=99
    @test r.geometry[1]==0
    times=[0.,1.,2.]; pos=[3. 4. 5.;0. 0. 0.;0. 0. 0.]; vel=zeros(3,3)
    ref=reference(t_ref_s=times,r_ref_rtn_m=pos,v_ref_rtn_mps=vel)
    times[1]=99;pos[1]=99;vel[1]=99
    @test ref.t_ref_s[1]==0 && ref.r_ref_rtn_m[1]==3 && ref.v_ref_rtn_mps[1]==0
    diag=(history=[1.,2.],);out=candidate(ref;diagnostics=diag)
    ref.r_ref_rtn_m[1]=99;diag.history[1]=99
    @test out.reference.r_ref_rtn_m[1]==3 && out.diagnostics.history[1]==1
end

@testset "Request preconditions" begin
    for kwargs in ((request_id=0,), (chaser_id=201,), (target_id=-1,),
        (state_source=:estimate,), (frame=:body,), (reason=:unknown,),
        (reference_dt_s=0.,), (reference_dt_s=Inf,), (time_s=NaN,), (time_s=1e30,valid_until_s=2e30,),
        (preview_horizon_steps=0,), (valid_until_s=1.,), (observation_time_s=1.,),
        (x_rtn=(NaN,0.,0.,0.,0.,0.),), (x_rtn=(1,2,3),),
        (goal_rtn_m=(1.,Inf,0.),), (geometry_revision="",),
        (target_state_ii=(0.,0.,0.,0.,1.,0.),),
        (target_state_ii=(1.,0.,0.,1.,0.,0.),),
        (validation=P.RPOValidationSettings(time_atol_s=1.),))
        @test_throws ArgumentError request(;kwargs...)
    end
    @test_throws ArgumentError P.RPOPlanningConstraints(clearance_m=-1.)
    @test_throws ArgumentError P.RPOPlanningConstraints(clearance_m=0.,max_speed_mps=Inf)
    @test_throws ArgumentError P.RPOPlanningConstraints(clearance_m=0.,max_acceleration_mps2=0.)
    @test_throws ArgumentError P.RPOValidationSettings(max_clearance_samples=1)
    @test_throws ArgumentError P.RPOValidationSettings(clearance_sample_ds_m=NaN)
    @test_throws ArgumentError candidate(;status=:failed)
    @test_throws ArgumentError candidate(;status=:success)
end

@testset "Reference acceptance, identity, time and shape" begin
    r=request(); out=candidate(); before=deepcopy(out.reference.r_ref_rtn_m)
    result=validate(r,out)
    @test result.accepted && result.reason===:validated
    @test result.metrics.min_clearance_m≈2.6
    @test result.metrics.clearance_samples==41
    @test result.metrics.max_speed_mps==1
    @test result.metrics.max_acceleration_mps2==1
    @test out.reference.r_ref_rtn_m==before
    @test validate(r,candidate(;termination=:error)).reason===:unusable_termination
    @test validate(r,candidate(;termination=:time_budget)).reason===:time_budget_not_allowed
    @test validate(r,candidate(;termination=:iteration_limit)).accepted
    permitted=request(validation=P.RPOValidationSettings(allow_time_budget_candidate=true))
    @test validate(permitted,candidate(;termination=:time_budget)).accepted
    @test validate(r,candidate(;request_id=2)).reason===:request_mismatch
    for (kwargs,reason) in [
        ((chaser_id=2,),:spacecraft_mismatch), ((target_id=2,),:spacecraft_mismatch),
        ((frame=:body,),:frame_mismatch), ((geometry_revision="old",),:geometry_mismatch),
        ((terminal_policy=:hold,),:unsupported_terminal_policy),
        ((t_ref_s=Float64[],),:reference_shape), ((v_ref_rtn_mps=zeros(2,3),),:reference_shape),
        ((origin_time_s=1.,),:origin_mismatch), ((valid_until_s=Inf,),:nonfinite_reference),
        ((valid_until_s=11.,),:invalid_validity), ((valid_until_s=0.,),:invalid_validity),
        ((t_ref_s=[0.1,1.,2.],),:time_origin), ((t_ref_s=[0.,1.,1.],),:nonuniform_time_grid),
        ((t_ref_s=[0.,0.9,2.],),:nonuniform_time_grid), ((valid_until_s=1.,),:reference_exceeds_validity),
        ((r_ref_rtn_m=zeros(3,3),),:endpoint_mismatch)]
        @test validate(r,candidate(reference(;kwargs...))).reason===reason
    end
    for bad in (NaN,Inf,-Inf), field in (:t_ref_s,:r_ref_rtn_m,:v_ref_rtn_mps)
        ref=reference();getfield(ref,field)[2]=bad
        @test validate(r,candidate(ref)).reason===:nonfinite_reference
    end
    @test P.rpo_reference_is_current(r,reference(),8.)
    @test !P.rpo_reference_is_current(r,reference(),8.01)
    @test !P.rpo_reference_is_current(r,reference(),-1.)
    @test !P.rpo_reference_is_current(r,reference(),NaN)
    @test validate(r,out;time_s=8.01).reason===:expired_or_not_yet_valid
    # Origin is simulation time; relative grid does not become absolute time.
    shifted=request(time_s=100.,valid_until_s=110.)
    @test validate(shifted,candidate(reference(origin_time_s=100.,valid_until_s=110.))).accepted
end

@testset "Independent feasibility checks and bounded work" begin
    @test P.validate_rpo_result(request(),candidate()).reason===:clearance_not_checked
    @test validate(request(constraints=P.RPOPlanningConstraints(clearance_m=0.1,max_speed_mps=0.9)),candidate()).reason===:speed_limit
    @test validate(request(constraints=P.RPOPlanningConstraints(clearance_m=0.1,max_acceleration_mps2=0.9)),candidate()).reason===:acceleration_limit
    @test validate(request(constraints=P.RPOPlanningConstraints(clearance_m=3.)),candidate()).reason===:clearance_limit
    calls=Ref(0)
    query(p,g)=(calls[]+=1;clearance(p,g))
    cap=request(validation=P.RPOValidationSettings(max_clearance_samples=3))
    @test P.validate_rpo_result(cap,candidate();clearance_at=query).reason===:clearance_work_limit
    @test calls[]==0
    for bad in (NaN,Inf,-Inf)
        @test P.validate_rpo_result(request(),candidate();clearance_at=(p,g)->bad).reason===:nonfinite_clearance
    end
    @test_throws ErrorException P.validate_rpo_result(request(),candidate();clearance_at=(p,g)->error("broken query"))
    # Safe endpoint knots cannot conceal an obstacle in the segment interior.
    r=request(x_rtn=(-1.,0.,0.,0.,0.,0.),goal_rtn_m=(1.,0.,0.))
    ref=reference(t_ref_s=[0.,1.],r_ref_rtn_m=[-1. 1.;0. 0.;0. 0.],v_ref_rtn_mps=zeros(3,2))
    @test validate(r,candidate(ref)).reason===:clearance_limit
    # Diagnostic Bezier control vertices are not the executable reference.
    @test validate(request(),candidate(;diagnostics=(path=(representation=:bezier_control_polygon,points=zeros(3,3)),))).accepted
    for status in (:infeasible,:unsupported,:failed)
        @test validate(request(),P.RPOPlanningResult(request_id=1,status=status,termination=:reported)).reason===status
    end
end

struct MissingPlanner <: P.AbstractRPOPlanner end
struct DeclaredPlanner <: P.AbstractRPOPlanner end
P.planner_capabilities(::DeclaredPlanner)=P.RPOPlannerCapabilities(state_sources=(:truth,),frames=(:target_rtn,))
@testset "Explicit extension and refusal behavior" begin
    r=request(); rng=MersenneTwister(7); untouched=copy(rng)
    @test P.initialize_planner(MissingPlanner(),nothing)===nothing
    @test !P.validate_rpo_capabilities(MissingPlanner(),r).accepted
    @test P.validate_rpo_capabilities(DeclaredPlanner(),r).accepted
    @test P.validate_rpo_capabilities(DeclaredPlanner(),r;checkpoint=true).reason===:unsupported_restart
    @test P.validate_rpo_capabilities(DeclaredPlanner(),request(reason=:retime)).reason===:unsupported_retiming
    @test P.plan_rpo!(nothing,MissingPlanner(),r,rng).status===:unsupported
    @test P.retime_rpo!(nothing,DeclaredPlanner(),r,reference(),rng).termination===:retiming_not_implemented
    @test rand(rng)==rand(untouched)
    @test_throws MethodError P.plan_rpo!(nothing,DeclaredPlanner(),r)
end
end
