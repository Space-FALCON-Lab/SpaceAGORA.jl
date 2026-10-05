module EnergyDepletionLifecycleTests
using Test, SpaceAGORA, OrdinaryDiffEq, DiffEqCallbacks, ComponentArrays
import .._edg_test_context
const S = SpaceAGORA.SimulationModel
const E = SpaceAGORA.SimulationEngine
const L = SpaceAGORA.SimulationLifecycle
const CB = S.SimulationCallbacks
const SC = S.SimConfig

function context(; n=1, targeting=false, control_only=false)
    base = _edg_test_context(guidance_modes=targeting ? (:targeting,) : (:max_energy_depletion,),
        max_energy_submodes=(:heat_rate,), heat_rate_limit_w_cm2=Inf,
        planning_horizon_s=1.0, target_apoapsis_radius_m=1e12)
    state = S.AerobrakingEnergyDepletionState(num_sats=n)
    g = S.AerobrakingEnergyDepletionGuidanceModel(base.config, state)
    c = S.AerobrakingEnergyDepletionControlModel(base.config, state)
    members = [deepcopy(base.spacecraft) for _ in 1:n]
    for (i,sc) in enumerate(members); sc.id = i; end
    a = SC._with_configuration(base.args;
        dynamics_model=S.DynamicsModel(members, base.args.dynamics_model.dynamic_effectors),
        guidance_model=control_only ? S.GuidanceModel((), Float64[]) : S.GuidanceModel((g,), [1.0]),
        control_model=S.ControlModel((c,), [0.5]),
        mission_configuration=S.MissionConfiguration(mission_type=S.MissionTime,
            mission_time=2.0, orientation_sim=false),
        simulation_settings=S.SimulationSettings(results=false, verbose=false),
        solver_config=S.SolverConfig(solver_mode=:tsit5))
    p = S.ODEParams(n_sats=n, args=a)
    p.shared_buffers.et_start[] = base.p.shared_buffers.et_start[]
    u = ComponentVector(sc=[(pos=collect(base.u.sc[1].pos), vel=collect(base.u.sc[1].vel),
        mass=461.0, heat_loads=fill(7.0,3)) for _ in 1:n])
    return (; state, g, c, a, p, u)
end
function seed!(s, i=1)
    s.selected_mode[i] = :targeting
    s.targeting_active[i] = true
    s.safe_low_drag[i] = true
    s.energy_bracketing_evaluated[i] = true
    s.energy_bracketing_count[i] = 7
    s.target_energy_jkg[i] = -4e6
    s.bracket_min_energy_jkg[i] = -5e6
    s.bracket_max_energy_jkg[i] = -3e6
    s.heat_load_switches_s[i] = (11.0,12.0)
    s.heat_load_switch_solved[i] = true
    s.heat_load_drag_passage_active[i] = true
    s.targeting_switch_s[i] = 42.0
    s.last_switch_solve_t[i] = 4.0
    s.last_alpha_rad[i] = 0.3
    s.last_alpha_heat_rate_rad[i] = 0.4
    s.last_alpha_structural_rad[i] = 0.5
    s.last_heat_rate_w_cm2[i] = 0.6
    s.last_heat_load_j_cm2[i] = 0.7
    s.last_dynamic_pressure_pa[i] = 0.8
    return s
end
snapshot(s, i) = Tuple(getfield(s,f)[i] for f in fieldnames(typeof(s)))
function initialize!(x)
    for g in x.a.guidance_model.guidance_effectors; L.initialize_guidance!(g,x.u,x.p,0.0); end
    for c in x.a.control_model.control_effectors; L.initialize_control!(c,x.u,x.p,0.0); end
end
mutable struct MaskOpts
    dtmax::Float64
    reltol::Any
    abstol::Any
    adaptive::Bool
end
function mask_integrator(x; t=3e7)
    rt,at = E._build_solver_tolerances(x.u,x.a)
    (; p=x.p, u=x.u, t, opts=MaskOpts(20.,rt,at,true))
end
function apply_mask!(x, events)
    boundary=x.a.environment_model.planet.Rp_e + x.a.environment_model.EI*1e3
    for i in eachindex(events)
        x.u.sc[i].pos .= (boundary + (events[i]>0 ? -1e-6 : events[i]<0 ? 1e-6 : -100.0),0,0)
        x.u.sc[i].vel .= (events[i]<0 ? -100.0 : 100.0,0,0)
    end
    CB.get_drag_state_callback(length(events)).affect!(mask_integrator(x),Int8.(events))
end

@testset "EDG fresh initialization and spacecraft-owned pass resets" begin
    for control_only in (false,true), offset in (-100.0,0.0,100.0)
        x=context(;control_only)
        seed!(x.state)
        boundary=x.a.environment_model.planet.Rp_e+x.a.environment_model.EI*1e3
        x.u.sc[1].pos .= (boundary+offset,0,0)
        before=copy(x.u)
        E._initialize_in_atmosphere_flags!(x.p,x.u)
        initialize!(x)
        @test x.p.shared_buffers.in_atmosphere == [offset<=0]
        @test isequal(snapshot(x.state,1),snapshot(S.AerobrakingEnergyDepletionState(num_sats=1),1))
        @test x.u == before
        @test x.g.state === x.c.state
    end
    for control_only in (false,true)
        x=context(n=3;control_only)
        for i in 1:3; seed!(x.state,i); end
        x.p.shared_buffers.in_atmosphere .= [true,false,true]
        untouched=snapshot(x.state,3)
        physical=copy(x.u.sc[1].heat_loads)
        apply_mask!(x,[1,-1,0])
        @test x.p.shared_buffers.in_atmosphere == [false,true,true]
        for i in 1:2
            @test !x.state.energy_bracketing_evaluated[i] && !x.state.targeting_active[i]
            @test !x.state.safe_low_drag[i] && x.state.selected_mode[i] == :inactive
            @test x.state.targeting_switch_s[i] == Inf
            @test all(isnan,(x.state.target_energy_jkg[i],x.state.bracket_min_energy_jkg[i],x.state.bracket_max_energy_jkg[i]))
            @test !x.state.heat_load_switch_solved[i] && !x.state.heat_load_drag_passage_active[i]
            @test x.state.heat_load_switches_s[i] == (Inf,Inf)
            @test x.state.last_switch_solve_t[i] == -Inf
            @test x.state.energy_bracketing_count[i] == 7
            @test x.state.last_alpha_rad[i] == 0.3 && x.state.last_heat_load_j_cm2[i] == 0.7
        end
        @test snapshot(x.state,3) == untouched
        @test x.u.sc[1].heat_loads == physical
        # A newly solved plan survives duplicate delivery of the same mask.
        seed!(x.state,1); seed!(x.state,2)
        before=[snapshot(x.state,i) for i in 1:3]
        apply_mask!(x,[1,-1,0])
        @test [snapshot(x.state,i) for i in 1:3] == before
        # Reverse the simultaneous directions; both changed members reset.
        apply_mask!(x,[-1,1,0])
        @test x.state.targeting_switch_s == [Inf,Inf,42.0]
        # An omitted simultaneous member is still reconciled by geometry.
        x.p.shared_buffers.in_atmosphere[3]=false
        apply_mask!(x,[-1,1,0])
        @test x.state.targeting_switch_s[3] == Inf
    end
end

@testset "Atmosphere startup and reconciliation share boundary direction" begin
    x=context()
    boundary=x.a.environment_model.planet.Rp_e+x.a.environment_model.EI*1e3
    for (offset, inbound, tangent, outbound) in (
        (-100.0, true, true, true),
        (-32eps(boundary), true, true, false),
        (0.0, true, true, false),
        (32eps(boundary), true, true, false),
        (100.0, false, false, false))
        for (speed, expected) in ((-100.0,inbound),(0.0,tangent),(100.0,outbound))
            x.u.sc[1].pos .= (boundary+offset,0,0)
            x.u.sc[1].vel .= (speed,10,0)
            E._initialize_in_atmosphere_flags!(x.p,x.u)
            @test x.p.shared_buffers.in_atmosphere == [expected]
            x.p.shared_buffers.in_atmosphere[1]=!expected
            CB._refresh_crossing_atmosphere_flags!(mask_integrator(x),Int8[0])
            @test x.p.shared_buffers.in_atmosphere == [expected]
        end
    end
end

@testset "EDG crossings remain scheduled with identical solver phases" begin
    for control_only in (false,true)
        x=context(;control_only)
        tol=S.IntegrationTolerances(reltol_orbit=1e-8,reltol_atmosphere=1e-8,
            abstol_orbit=1e-10,abstol_atmosphere=1e-10,dt_max_orbit=1.0,dt_max_atmosphere=1.0)
        a=SC._with_configuration(x.a;integration_tolerances=tol)
        @test CB._requires_drag_state_callback(a.dynamics_model.dynamic_effectors,a)
        # Removing EDG retains the old inert-default selection rule.
        legacy=SC._with_configuration(a;guidance_model=S.GuidanceModel((),Float64[]),
            control_model=S.ControlModel((),Float64[]))
        @test !CB._requires_drag_state_callback(legacy.dynamics_model.dynamic_effectors,legacy)
    end
end

@testset "EDG preflight rejects incomplete ownership and unsupported continuation" begin
    x=context()
    @test isnothing(L.preflight_guidance(x.g,x.a;isolate_state=true))
    @test isnothing(L.preflight_control(x.c,x.a;isolate_state=false))
    alone=context(control_only=true)
    @test isnothing(L.preflight_control(alone.c,alone.a;isolate_state=true))
    for field in fieldnames(typeof(x.state))
        y=deepcopy(x); pop!(getfield(y.state,field))
        @test_throws ArgumentError L.preflight_guidance(y.g,y.a;isolate_state=true)
        @test_throws ArgumentError L.preflight_control(y.c,y.a;isolate_state=true)
    end
    different=S.AerobrakingEnergyDepletionControlModel(x.c.config,deepcopy(x.state))
    a=SC._with_configuration(x.a;control_model=S.ControlModel((different,),[1.0]))
    @test_throws ArgumentError L.preflight_guidance(x.g,a;isolate_state=true)
    other_config=S.AerobrakingEnergyDepletionConfig(max_alpha_rad=0.1)
    different=S.AerobrakingEnergyDepletionControlModel(other_config,x.state)
    a=SC._with_configuration(x.a;control_model=S.ControlModel((different,),[1.0]))
    @test_throws ArgumentError L.preflight_control(different,a;isolate_state=true)
    a=SC._with_configuration(x.a;guidance_model=S.GuidanceModel((x.g,x.g),[1.0,2.0]))
    @test_throws ArgumentError L.preflight_guidance(x.g,a;isolate_state=true)
    a=SC._with_configuration(x.a;control_model=S.ControlModel((x.c,x.c),[1.0,2.0]))
    @test_throws ArgumentError L.preflight_control(x.c,a;isolate_state=true)
    # Single-sided duplicates must fail for the explicit EDG ownership rule,
    # without relying on only(...) in the paired-model path to throw.
    for guidance_only in (false,true)
        a=SC._with_configuration(x.a;
            guidance_model=guidance_only ? S.GuidanceModel((x.g,x.g),[1.0,2.0]) : S.GuidanceModel((),Float64[]),
            control_model=guidance_only ? S.ControlModel((),Float64[]) : S.ControlModel((x.c,x.c),[1.0,2.0]))
        err=try
            guidance_only ? L.preflight_guidance(x.g,a;isolate_state=true) : L.preflight_control(x.c,a;isolate_state=true)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("EDG requires at most one guidance model and one control model",sprint(showerror,err))
    end
    y=context(targeting=true,control_only=true)
    @test_throws ArgumentError L.preflight_control(y.c,y.a;isolate_state=true)
    # Both write and resume fail through the public engine before creating output.
    for control_only in (false,true), resume in (false,true)
        mktempdir() do dir
            x=context(;control_only)
            a=SC._with_configuration(x.a;simulation_settings=S.SimulationSettings(
                results=true,verbose=false,results_directory=joinpath(dir,"output"),
                checkpoint_enabled=!resume,resume_from_checkpoint=resume,
                checkpoint_directory=joinpath(dir,"checkpoints")))
            err=try run_simulation(a); nothing catch e; e end
            @test err isa ArgumentError
            @test occursin("EDG does not support checkpoint",sprint(showerror,err))
            @test isempty(readdir(dir))
        end
    end
end

@testset "EDG fresh and reused configurations through simulation initialization" begin
    for control_only in (false,true)
        x=context(;control_only)
        seed!(x.state)
        original=snapshot(x.state,1)
        first=run_simulation(x.a;return_solution=true)
        second=run_simulation(x.a;return_solution=true)
        fresh=run_simulation(context(;control_only).a;return_solution=true)
        @test first.u == second.u == fresh.u
        @test first.t == second.t == fresh.t
        @test snapshot(x.state,1) == original # Default isolation leaves caller state intact.
        owned=only(first.prob.p.args.control_model.control_effectors).state
        @test owned !== x.state && owned.energy_bracketing_count == [0]
        @test owned.targeting_switch_s == [Inf]
        @test owned.last_alpha_rad[1] == x.c.config.max_alpha_rad
        if !control_only
            @test only(first.prob.p.args.guidance_model.guidance_effectors).state === owned
        end
        # Explicit caller ownership still starts EDG afresh, rather than resuming.
        reused=run_simulation(x.a;isolate_state=false,return_solution=true)
        @test only(reused.prob.p.args.control_model.control_effectors).state === x.state
        @test x.state.energy_bracketing_count == [0] && x.state.targeting_switch_s == [Inf]
    end
end

@testset "EDG real scheduled callbacks across a manufactured exit and reentry" begin
    x=context(targeting=true)
    seed!(x.state); x.state.energy_bracketing_count[1]=1
    boundary=x.a.environment_model.planet.Rp_e+x.a.environment_model.EI*1e3
    x.u.sc[1].pos .= (boundary-100,0,0)
    x.u.sc[1].vel .= (0,0,0)
    E._initialize_in_atmosphere_flags!(x.p,x.u)
    function radial!(du,u,p,t)
        du .= 0
        du.sc[1].pos[1]=10pi*sin(pi*t/10)
        du.sc[1].vel[1]=pi^2*cos(pi*t/10)
    end
    records=NamedTuple[]
    record=PeriodicCallback(i->push!(records,(t=i.t,inside=i.p.shared_buffers.in_atmosphere[1],
        count=x.state.energy_bracketing_count[1],switch=x.state.targeting_switch_s[1])),0.5)
    callbacks=CallbackSet(CB.get_drag_state_callback(1),CB.get_control_callbacks(1,x.a)...,
        CB.get_guidance_callbacks(1,x.a)...,record)
    sol=solve(ODEProblem(radial!,x.u,(0.,19.),x.p),Tsit5();callback=callbacks,
        dtmax=0.25,abstol=1e-8,reltol=1e-10)
    @test string(sol.retcode)=="Success" && last(sol.t)==19.0
    @test any(!r.inside for r in records) && any(r.t>15 && r.inside for r in records)
    @test all(r->r.switch==42.0 && r.count==1,filter(r->r.t<4.5,records))
    # Counts measure evaluations, not passes. A guidance tick exactly at exit
    # can still satisfy the existing <= EI predictor predicate. Entry must
    # invalidate even that boundary solve, then cache one new bracket.
    before_reentry=last(filter(r->r.t<15,records)).count
    after_reentry=filter(r->r.t>=16,records)
    @test x.state.energy_bracketing_count==[before_reentry+1]
    @test all(r->r.count==before_reentry+1,after_reentry)
    @test x.state.targeting_switch_s==[Inf]
    @test x.u.sc[1].heat_loads==fill(7.0,3)
    println("EDG_SCHEDULED_REENTRY ",records)
end

@testset "EDG outbound boundary start invalidates t0 plan on reentry" begin
    x=context(targeting=true)
    boundary=x.a.environment_model.planet.Rp_e+x.a.environment_model.EI*1e3
    x.u.sc[1].pos .= (boundary,0,0)
    x.u.sc[1].vel .= (100pi/10.6,0,0)
    # A real thruster initialization invokes guidance at t0. No burn is requested.
    thruster=S.BaseThrusterModel(thrust=[0.0],direction=[0.0],Δv=[0.0],
        start_burn_time=[Inf],stop_burn_time=[Inf],Isp=[300.0])
    a=SC._with_configuration(x.a;control_model=S.ControlModel((x.c,thruster),[0.5,1.0]))
    p=S.ODEParams(n_sats=1,args=a)
    p.shared_buffers.et_start[]=x.p.shared_buffers.et_start[]
    E._initialize_in_atmosphere_flags!(p,x.u)
    @test p.shared_buffers.in_atmosphere == [false]
    function radial!(du,u,p,t)
        du .= 0
        du.sc[1].pos[1]=u.sc[1].vel[1]
        du.sc[1].vel[1]=-100*(pi/10.6)^2*sin(pi*t/10.6)
    end
    records=NamedTuple[]
    record=PeriodicCallback(i->push!(records,(t=i.t,inside=i.p.shared_buffers.in_atmosphere[1],
        count=x.state.energy_bracketing_count[1])),0.5)
    callbacks=CallbackSet(CB.get_drag_state_callback(1),CB.get_control_callbacks(1,a)...,
        CB.get_guidance_callbacks(1,a)...,record)
    integrator=init(ODEProblem(radial!,x.u,(0.,19.),p),Tsit5();callback=callbacks,
        dtmax=0.25,abstol=1e-8,reltol=1e-10)
    @test x.state.energy_bracketing_count == [1] # The production t0 trigger actually ran.
    solve!(integrator)
    @test string(integrator.sol.retcode)=="Success" && integrator.t==19.0
    coast=filter(r->r.t<10.5,records)
    nextpass=filter(r->r.t>=11,records)
    @test !isempty(coast) && all(r->!r.inside && r.count==1,coast)
    @test !isempty(nextpass) && all(r->r.inside && r.count==2,nextpass)
    @test x.state.energy_bracketing_count == [2]
    @test x.u.sc[1].heat_loads==fill(7.0,3)
    println("EDG_BOUNDARY_REENTRY ",records)
end
end
