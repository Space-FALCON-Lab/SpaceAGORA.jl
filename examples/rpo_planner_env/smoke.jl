module RPOPublicPilot
using SpaceAGORA, Test, Random, StaticArrays

# Every package symbol in this file comes from the supported root surface.
function hypr()
    HYPRRPOPlanner(RPOPSOConfig(n_waypoints=2,n_particles=8,n_iters=3,
        adaptive_enable=false,adaptive_n_waypoints_min=2,adaptive_n_waypoints_max=2,
        adaptive_n_particles_min=8,adaptive_n_particles_max=8,
        adaptive_n_iters_min=3,adaptive_n_iters_max=3,sample_ds_m=.2,
        retime_accel_limit_enable=true,retime_max_speed_mps=.5,
        retime_a_max_mps2=.5*.05/5.2,iteration_runtime_limit_s=Inf,rrt_warmstart_enable=true))
end

# A small user-owned planner. Valid only for an unchanged goal, by design.
struct HoldPlanner <: AbstractRPOPlanner
    fail_reason::Symbol
end
HoldPlanner()=HoldPlanner(:none)
SpaceAGORA.planner_capabilities(::HoldPlanner)=RPOPlannerCapabilities(
    state_sources=(:truth,),frames=(:target_rtn,),retiming=true)
SpaceAGORA.initialize_planner(::HoldPlanner,context)=Ref(0)
function SpaceAGORA.plan_rpo!(state,p::HoldPlanner,r::RPOPlanningRequest,rng::AbstractRNG)
    state[]+=1
    p.fail_reason===r.reason && return RPOPlanningResult(request_id=r.request_id,
        status=:infeasible,termination=:test_infeasibility)
    ref=RPOReference(t_ref_s=[0.,r.reference_dt_s],
        r_ref_rtn_m=hcat(collect(r.x_rtn[1:3]),collect(r.goal_rtn_m)),v_ref_rtn_mps=zeros(3,2),
        origin_time_s=r.time_s,valid_until_s=r.valid_until_s,chaser_id=r.chaser_id,
        target_id=r.target_id,geometry_revision=r.geometry_revision)
    RPOPlanningResult(request_id=r.request_id,status=:candidate,termination=:completed,
        reference=ref,diagnostics=(calls=state[],draw=rand(rng)))
end
SpaceAGORA.retime_rpo!(state,p::HoldPlanner,r::RPOPlanningRequest,active::RPOReference,rng::AbstractRNG)=
    plan_rpo!(state,p,r,rng)

struct SmallAcceleration <: AbstractForceTorqueModel end
SpaceAGORA.wrench(::SmallAcceleration,x::StateSample,env::EnvironmentSample,t::Float64)=
    (SVector(0.,0.,.01*x.mass_kg),SVector(0.,0.,0.))

@testset "Public planner interchangeability, external planner and force" begin
    for planner in (DirectRPOPlanner(),hypr())
        args=make_rpo_configuration(planner=planner)
        sol=run_simulation(args;return_solution=true)
        report=rpo_run_report(sol)
        @test last(sol.t)==2.
        @test all(isfinite,sol.u[end])
        @test report.planners[1].request_id==1
        installed=only(filter(x->x.event===:installed,report.planners[1].records))
        @test installed.validation.accepted
        @test installed.reference.chaser_id==101
        @test !isempty(report.commands[1].log.t_s)
        @test all(s->s in (:Solved,:Solved_inaccurate),report.commands[1].log.qp_status)
        @test all(f->maximum(f)<=.05,report.commands[1].log.thruster_forces_n)
        println("PILOT ",nameof(typeof(planner))," samples=",length(sol.t),
            " control_updates=",length(report.commands[1].log.t_s),
            " reference_duration=",last(installed.reference.t_ref_s))
    end
    hold=make_rpo_configuration(planner=HoldPlanner(),goal_rtn_m=(3.,0.,0.),mission_time_s=.05)
    base=run_simulation(hold;return_solution=true)
    forced=run_simulation(make_rpo_configuration(planner=HoldPlanner(),goal_rtn_m=(3.,0.,0.),
        mission_time_s=.05,extra_effectors=(SmallAcceleration(),));return_solution=true)
    @test rpo_run_report(base).planners[1].request_id==1
    @test forced.u[end].sc[1].pos[3]-base.u[end].sc[1].pos[3] ≈ .5*.01*.05^2 atol=1e-9 rtol=1e-8
    @test forced.u[end].sc[1].vel[3]-base.u[end].sc[1].vel[3] ≈ .01*.05 atol=1e-9 rtol=1e-8
end
end
