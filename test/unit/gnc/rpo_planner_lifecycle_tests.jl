module RPOPlannerLifecycleTests
using Test, SpaceAGORA, Random, LinearAlgebra
const S=SpaceAGORA.SimulationModel
const E=SpaceAGORA.SimulationEngine
const L=SpaceAGORA.RPOPlannerLifecycle
const Hooks=SpaceAGORA.SimulationLifecycle
include(joinpath(@__DIR__,"..","..","..","examples","rpo_planner_env","smoke.jl"))
const Hold=RPOPublicPilot.HoldPlanner
function holdargs(;kwargs...)
    make_rpo_configuration(;planner=Hold(),goal_rtn_m=(3.,0.,0.),kwargs...)
end
function prepare(args)
    owned=deepcopy(args)
    u=E.build_initial_conditions(owned)
    p=S.ODEParams(n_sats=length(owned.dynamics_model.spacecraft),args=owned)
    for g in owned.guidance_model.guidance_effectors
        Hooks.initialize_guidance!(g,u,p,0.)
    end
    return owned,u,p,only(owned.guidance_model.guidance_effectors)
end
function outcome(f)
    try f(); nothing catch e; e end
end

@testset "Preparation, callback order, duplicates, held control and RHS" begin
    a=holdargs(planning_events=(RPOPlanningEvent(.1;token=1),RPOPlanningEvent(.1;token=1),
        RPOPlanningEvent(.1;token=2),RPOPlanningEvent(.2;reason=:retime,token=3)),mission_time_s=.3)
    owned,u,p,g=prepare(a)
    @test g.runtime.request_id==1
    @test !a.guidance_model.guidance_effectors[1].plan_buffer.valid
    @test g.plan_buffer===only(owned.control_model.control_effectors).plan_buffer
    g2=deepcopy(g); @test g2.runtime.rng!==g.runtime.rng
    c=only(owned.control_model.control_effectors)
    S.calcControlEffect!(c,u,p,.1,1)
    @test last(g.runtime.records).event===:control && last(g.runtime.records).request_id==1
    S.GuidanceHooks.calcGuidanceEffect!(g,u,p,.1,1)
    @test g.runtime.request_id==3
    rng=copy(g.runtime.rng)
    S.GuidanceHooks.calcGuidanceEffect!(g,u,p,.1,1)
    @test g.runtime.request_id==3 && rand(copy(rng))==rand(copy(g.runtime.rng))
    S.GuidanceHooks.calcGuidanceEffect!(g,u,p,.2,1)
    @test g.runtime.request_id==4
    @test last(g.runtime.records).reason===:retime
    push!(g.events,RPOPlanningEvent(.2;token=4))
    S.GuidanceHooks.calcGuidanceEffect!(g,u,p,.2,1)
    @test g.runtime.request_id==5 # A new token at the same accepted time is not a duplicate.
    S.GuidanceHooks.calcGuidanceEffect!(g,u,p,.2,1)
    @test g.runtime.request_id==5
    sol=run_simulation(a;return_solution=true)
    report=rpo_run_report(sol)
    rec=report.planners[1].records
    @test report.planners[1].request_id==4
    at=filter(x->hasproperty(x,:time_s)&&x.time_s==.1,rec)
    @test findfirst(x->x.event===:control,at)<findfirst(x->x.event===:guidance,at)
    @test first(filter(x->x.event===:control,at)).request_id==1
    # Full RHS calls at arbitrary trial times cannot call planner or consume RNG.
    params=sol.prob.p;run_g=only(params.args.guidance_model.guidance_effectors)
    saved=deepcopy((run_g.runtime.request_id,run_g.runtime.rng,repr(run_g.runtime.records)))
    ctrl=only(params.args.control_model.control_effectors);held=ctrl.held
    du=similar(sol.u[end])
    for t in (.301,.299,.301)
        E.spacecraft_dynamics!(du,sol.u[end],params,t)
    end
    @test run_g.runtime.request_id==saved[1]
    @test rand(copy(run_g.runtime.rng))==rand(saved[2])
    @test repr(run_g.runtime.records)==saved[3]
    @test ctrl.held===held
    # Accessors cannot modify the run-owned reference or command arrays.
    report.planners[1].active_reference.r_ref_rtn_m .= 99
    empty!(report.commands[1].log.t_s)
    fresh=rpo_run_report(sol)
    @test maximum(fresh.planners[1].active_reference.r_ref_rtn_m)<99
    @test !isempty(fresh.commands[1].log.t_s)
end

@testset "Explicit failures and consumption expiry" begin
    for reason in (:initial,:forced,:replan,:retime)
        a=make_rpo_configuration(planner=Hold(reason),goal_rtn_m=(3.,0.,0.),mission_time_s=.2,
            planning_events=reason===:initial ? () : (RPOPlanningEvent(.1;reason=reason),))
        e=outcome(()->run_simulation(a;return_solution=true))
        @test e isa RPOPlanningError
        @test e.chaser_id==101 && e.reason===:infeasible
        @test last(e.report.records).details.termination===:test_infeasibility
    end
    a=make_rpo_configuration(planner=DirectRPOPlanner(),goal_rtn_m=(3.,0.,0.),mission_time_s=.2,
        planning_events=(RPOPlanningEvent(.1;reason=:retime),))
    e=outcome(()->run_simulation(a;return_solution=true))
    @test e isa RPOPlanningError && e.reason===:unsupported_retiming
    a=holdargs(plan_validity_s=12*.1,mission_time_s=.2)
    owned,u,p,g=prepare(a);c=only(owned.control_model.control_effectors)
    @test isnothing(Hooks.before_reference_control!(g,c,u,p,0.))
    e=outcome(()->Hooks.before_reference_control!(g,c,u,p,.1))
    @test e isa RPOPlanningError && e.reason===:expired_reference
    @test e.report.request_id==1
    @test g.plan_buffer.valid # preserved for forensics, never consumed by control
    e=outcome(()->run_simulation(a;return_solution=true))
    @test e isa RPOPlanningError && e.reason===:expired_reference
    a=holdargs(replan_interval_s=.1,mission_time_s=.2)
    owned,u,p,g=prepare(a)
    S.GuidanceHooks.calcGuidanceEffect!(g,u,p,.1,1)
    @test last(g.runtime.records).reason===:replan
end

@testset "Periodic cadence and controller-preview lifetime" begin
    dt = 0.1
    for stride in (1, 2, 3, 5)
        interval = stride * dt
        sol = run_simulation(holdargs(replan_interval_s=interval, mission_time_s=2.0);
            return_solution=true)
        records = only(rpo_run_report(sol).planners).records
        installed = filter(x -> x.event === :installed && x.reason === :replan, records)
        expected = [k * dt for k in stride:stride:19]
        actual = [x.time_s for x in installed]
        @test length(actual) == length(expected)
        @test length(actual) == length(expected) &&
            all(isapprox(a, b; atol=1e-10, rtol=0) for (a, b) in zip(actual, expected))
        @test [x.request_id for x in installed] == collect(2:length(expected)+1)
        println("PERIODIC_CADENCE interval=", interval, " times=", actual)
    end
    # Control uses the old plan before guidance replaces it at the same tick.
    # Each lifetime covers that interval plus the full 12-step preview.
    for (interval, lifetime) in ((0.1, 1.3), (0.2, 1.4))
        sol = nothing
        err = outcome(() -> (sol = run_simulation(holdargs(replan_interval_s=interval,
            plan_validity_s=lifetime, mission_time_s=2.0); return_solution=true)))
        @test isnothing(err)
        if isnothing(err)
            records = only(rpo_run_report(sol).planners).records
            @test last(sol.t) == 2.0
            @test count(x -> x.event === :control, records) == 19
            @test count(x -> x.event === :failure, records) == 0
        end
        println("PERIODIC_LIFETIME interval=", interval, " lifetime=", lifetime,
            " outcome=", isnothing(err) ? :completed : err.reason)
    end
    # The existing absolute time tolerance admits a near-boundary tick only.
    owned, u, p, g = prepare(holdargs(replan_interval_s=0.2))
    atol = g.validation.time_atol_s
    S.GuidanceHooks.calcGuidanceEffect!(g, u, p, 0.2 - 2atol, 1)
    @test g.runtime.request_id == 1
    S.GuidanceHooks.calcGuidanceEffect!(g, u, p, 0.2 - atol/2, 1)
    @test g.runtime.request_id == 2
    rng = copy(g.runtime.rng)
    S.GuidanceHooks.calcGuidanceEffect!(g, u, p, 0.2 - atol/2, 1)
    @test g.runtime.request_id == 2
    @test rand(copy(rng)) == rand(copy(g.runtime.rng))
    owned, u, p, g = prepare(holdargs())
    for tick in 1:19
        S.GuidanceHooks.calcGuidanceEffect!(g, u, p, tick * dt, 1)
    end
    @test g.runtime.request_id == 1 # Infinite interval keeps background planning disabled.
end

@testset "Independent runs and ID-based streams" begin
    args=holdargs(mission_time_s=.2)
    s1=run_simulation(args;return_solution=true);s2=run_simulation(args;return_solution=true)
    @test s1.t==s2.t && s1.u==s2.u
    @test args.guidance_model.guidance_effectors[1].runtime===nothing
    tasks=[Threads.@spawn run_simulation(args;return_solution=true) for _ in 1:2]
    concurrent=fetch.(tasks)
    @test all(s->s.t==s1.t&&s.u==s1.u,concurrent)
    specs=((id=101,start_rtn_m=(3.,0.,0.),goal_rtn_m=(3.,0.,0.)),
           (id=102,start_rtn_m=(5.,0.,0.),goal_rtn_m=(5.,0.,0.)))
    reports=[rpo_run_report(run_simulation(make_rpo_configuration(planner=Hold(),chasers=x,
        mission_time_s=.2);return_solution=true)) for x in (specs,reverse(specs))]
    for id in (101,102)
        firsts=[only(filter(g->g.chaser_id==id,r.planners)) for r in reports]
        @test firsts[1].stream_words==firsts[2].stream_words
        indices=[findfirst(g->g.chaser_id==id,r.planners) for r in reports]
        @test reports[1].states[end].sc[indices[1]].pos ≈ reports[2].states[end].sc[indices[2]].pos atol=1e-8 rtol=1e-8
        @test reports[1].states[end].sc[indices[1]].vel ≈ reports[2].states[end].sc[indices[2]].vel atol=1e-8 rtol=1e-8
        @test only(filter(x->x.event===:installed,firsts[1].records)).diagnostics.draw==
              only(filter(x->x.event===:installed,firsts[2].records)).diagnostics.draw
    end
    @test reports[1].planners[1].stream_words!=reports[1].planners[2].stream_words
end

@testset "Preflight refusal before outputs" begin
    mktempdir() do folder
        output=joinpath(folder,"must-not-exist")
        for settings in (SimulationSettings(results=true,results_directory=output,checkpoint_enabled=true),
                SimulationSettings(results=true,results_directory=output,resume_from_checkpoint=true))
            a=holdargs(simulation_settings=settings)
            @test_throws ArgumentError run_simulation(a)
            @test !ispath(output)
        end
        @test_throws ArgumentError run_simulation(holdargs();isolate_state=false)
    end
end

struct InjectPlanner <: AbstractRPOPlanner
    defect::Symbol
end
SpaceAGORA.planner_capabilities(::InjectPlanner)=RPOPlannerCapabilities(state_sources=(:truth,),frames=(:target_rtn,))
SpaceAGORA.initialize_planner(::InjectPlanner,ctx)=Ref{Any}(nothing)
function SpaceAGORA.plan_rpo!(state,planner::InjectPlanner,r::RPOPlanningRequest,rng::AbstractRNG)
    d=planner.defect
    d===:throw && error("preserve original diagnostic")
    if d===:mutate_geometry
        r.geometry.station.points_body .= 100
    end
    times=[0.,r.reference_dt_s];positions=hcat(collect(r.x_rtn[1:3]),collect(r.goal_rtn_m));velocities=zeros(3,2)
    d===:nan && (positions[1,1]=NaN)
    d===:inf && (velocities[1,1]=Inf)
    d===:shape && (positions=zeros(2,2))
    d===:empty && (times=Float64[];positions=zeros(3,0);velocities=zeros(3,0))
    d===:dt && (times[2]*=2)
    d===:endpoint && (positions[1,end]+=1)
    d===:speed && (velocities[1,2]=100)
    d===:acceleration && (velocities[1,2]=.01)
    ref=RPOReference(t_ref_s=times,r_ref_rtn_m=positions,v_ref_rtn_mps=velocities,
        origin_time_s=r.time_s,valid_until_s=r.valid_until_s,
        chaser_id=d===:id ? 999 : r.chaser_id,target_id=r.target_id,
        frame=d===:frame ? :inertial : :target_rtn,
        geometry_revision=d===:revision ? "stale" : r.geometry_revision)
    result=RPOPlanningResult(request_id=r.request_id,status=:candidate,termination=:completed,
        reference=ref,diagnostics=(mutable_values=[1,2],))
    state[]=result
    return result
end
@testset "Untrusted request and result boundaries" begin
    expected=(:nan=>:nonfinite_reference,:inf=>:nonfinite_reference,:shape=>:reference_shape,
        :empty=>:reference_shape,:dt=>:nonuniform_time_grid,:endpoint=>:endpoint_mismatch,
        :speed=>:speed_limit,:acceleration=>:acceleration_limit,:id=>:spacecraft_mismatch,:frame=>:frame_mismatch,:revision=>:geometry_mismatch)
    for (defect,reason) in expected
        a=make_rpo_configuration(planner=InjectPlanner(defect),goal_rtn_m=(3.,0.,0.))
        e=outcome(()->prepare(a))
        @test e isa RPOPlanningError && e.reason===reason
        @test count(x->x.event===:installed,e.report.records)==0
    end
    a=make_rpo_configuration(planner=InjectPlanner(:mutate_geometry),goal_rtn_m=(3.,0.,0.),
        station_points_m=reshape([3.,0.,0.],3,1))
    e=outcome(()->prepare(a))
    @test e isa RPOPlanningError && e.reason===:clearance_limit
    @test a.guidance_model.guidance_effectors[1].geometry.station.points_body==reshape([3.,0.,0.],3,1)
    a=make_rpo_configuration(planner=InjectPlanner(:throw),goal_rtn_m=(3.,0.,0.))
    e=outcome(()->prepare(a))
    @test e isa RPOPlanningError && e.reason===:planner_exception
    @test e.cause isa ErrorException && occursin("original diagnostic",sprint(showerror,e.cause))
    @test !isempty(e.backtrace)
    a=make_rpo_configuration(planner=InjectPlanner(:none),goal_rtn_m=(3.,0.,0.))
    owned,u,p,g=prepare(a)
    retained=g.runtime.state[]
    retained.reference.r_ref_rtn_m .= 99
    retained.diagnostics.mutable_values .= 99
    @test g.plan_buffer.plan.r_ref_rtn[1,1]==3.
    @test g.plan_buffer.plan.diagnostics.planner.mutable_values==[1,2]
    @test g.runtime.reference.r_ref_rtn_m[1,1]==3.
end

using DiffEqCallbacks: PresetTimeCallback
@testset "Matched legacy HYPR references and propagated control" begin
    for (mode,adaptive,forced) in ((:legacy,false,false),(:manuscript,false,false),(:legacy,true,true))
        cfg=rpo_pso_config(RPOPublicPilot.hypr().config;hypr_mode=mode,adaptive_enable=adaptive)
        args=make_rpo_configuration(planner=HYPRRPOPlanner(cfg),mission_time_s=.4,
            planning_events=forced ? (RPOPlanningEvent(.1),) : ())
        sol=run_simulation(args;return_solution=true)
        report=rpo_run_report(sol);g=report.planners[1]
        requests=filter(x->x.event===:request,g.records)
        installed=filter(x->x.event===:installed,g.records)
        rng=MersenneTwister(g.stream_words)
        plans=S.RPOPlan[]
        for (rq,rec) in zip(requests,installed)
            req=rq.request;cfg=rec.diagnostics.input_config
            raw=S.rpo_pso_plan_path(req.x_rtn[1:3],req.goal_rtn_m,req.geometry,cfg;
                safe_distance_m=req.constraints.clearance_m,rng=rng)
            times,positions,velocities=S.rpo_reference_from_path(raw.path,req.geometry,raw.config;
                safe_distance_m=req.constraints.clearance_m)
            @test rec.reference.t_ref_s==times
            @test rec.reference.r_ref_rtn_m==positions
            @test rec.reference.v_ref_rtn_mps==velocities
            @test rec.diagnostics.effective_config==raw.config
            push!(plans,S.RPOPlan(valid=true,t_ref_s=times,r_ref_rtn=positions,v_ref_rtn=velocities,
                path_rtn=raw.path,cost=raw.cost))
        end
        legacy=deepcopy(args)
        control=only(legacy.control_model.control_effectors)
        S.update_rpo_plan_buffer!(control.plan_buffer,plans[1],0.)
        guide=S.RPOGuidanceModel(chaser_idx=1,target_idx=2,plan_buffer=control.plan_buffer)
        legacy=S.SimConfig._with_configuration(legacy;guidance_model=S.GuidanceModel(
            guidance_effectors=(guide,),guidance_rates=[.1]))
        callback=PresetTimeCallback([.1],integ->begin
            ctrl=only(integ.p.args.control_model.control_effectors)
            S.update_rpo_plan_buffer!(ctrl.plan_buffer,deepcopy(plans[2]),.1)
        end)
        old=run_simulation(legacy;return_solution=true,extra_callbacks=forced ? (callback,) : ())
        @test all(isapprox(sol(t),old(t);atol=1e-8,rtol=1e-8) for t in 0.:.1:.4)
        oldlog=only(old.prob.p.args.control_model.control_effectors).command_log
        newlog=report.commands[1].log
        @test newlog.t_s==oldlog.t_s
        @test all(isapprox(a,b;atol=1e-8,rtol=1e-8) for (a,b) in zip(newlog.accel_cmd_rtn,oldlog.accel_cmd_rtn))
        println("LEGACY_PARITY ",mode," adaptive=",adaptive," forced=",forced,
            " exact_states=",all(sol(t)==old(t) for t in 0.:.1:.4),
            " exact_commands=",newlog.accel_cmd_rtn==oldlog.accel_cmd_rtn)
    end
end
end
