module EnergyDepletionLifecycleTests
using Test, SpaceAGORA, OrdinaryDiffEq, DiffEqCallbacks, ComponentArrays, StaticArrays, LinearAlgebra
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
snapshot(s, i) = Tuple(deepcopy(getfield(s,f) isa AbstractVector ? getfield(s,f)[i] : getfield(s,f)) for f in fieldnames(typeof(s)))
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
# Unit-level geometric reconciliation, separate from the original spherical
# phase mask. Actual observer execution is exercised by the solver tests below.
function observe_members!(x; t=3e7)
    i=mask_integrator(x;t=t)
    for idx in eachindex(x.p.is_active)
        side=CB._edg_observed_boundary_side(x.u,x.p,t,idx)
        CB._edg_reconcile_observed_heat!(x.u,i,t,t,idx,side.inside)
    end
end
function initialize_observer_side!(x; t=0.)
    cb=CB.get_edg_heat_callback(length(x.p.is_active))
    d=first(cb.discrete_callbacks)
    d.initialize(d,x.u,t,mask_integrator(x;t=t))
end
function apply_mask!(x, events)
    boundary=x.a.environment_model.planet.Rp_e + x.a.environment_model.EI*1e3
    for i in eachindex(events)
        x.u.sc[i].pos .= (boundary + (events[i]==0 ? -100. : 0.),0,0)
        x.u.sc[i].vel .= (events[i]<0 ? -100.0 : 100.0,0,0)
    end
    observe_members!(x)
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
        x.state.heat_boundary_inside .= [true,false,true]
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
        @test isequal(snapshot(x.state,3),untouched)
        @test x.u.sc[1].heat_loads == physical
        # A newly solved plan survives duplicate delivery of the same mask.
        seed!(x.state,1); seed!(x.state,2)
        before=[snapshot(x.state,i) for i in 1:3]
        apply_mask!(x,[1,-1,0])
        @test isequal([snapshot(x.state,i) for i in 1:3],before)
        # Reverse the simultaneous directions; both changed members reset.
        apply_mask!(x,[-1,1,0])
        @test x.state.targeting_switch_s == [Inf,Inf,42.0]
        # An omitted simultaneous member is still reconciled by geometry.
        x.p.shared_buffers.in_atmosphere[3]=false
        x.state.heat_boundary_inside[3]=false
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
        getfield(x.state,field) isa AbstractVector || continue
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
        @test isequal(snapshot(x.state,1),original) # Default isolation leaves caller state intact.
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

# These manufactured trajectories test callback/lifecycle ownership only. Their
# radial velocities are not physical targeting inputs. Physical targeting and
# repeated guidance-owned replanning are covered by edg_replan_callback_tests.jl.
@testset "EDG scheduled cache invalidation across manufactured exit and reentry" begin
    x=context()
    seed!(x.state); x.state.targeting_active[1]=false
    x.state.energy_bracketing_count[1]=1
    boundary=x.a.environment_model.planet.Rp_e+x.a.environment_model.EI*1e3
    x.u.sc[1].pos .= (boundary-100,0,0)
    x.u.sc[1].vel .= (0,0,0)
    x.state.heat_boundary_inside[1]=true
    E._initialize_in_atmosphere_flags!(x.p,x.u)
    function radial!(du,u,p,t)
        du .= 0
        du.sc[1].pos[1]=10pi*sin(pi*t/10)
        du.sc[1].vel[1]=pi^2*cos(pi*t/10)
    end
    records=NamedTuple[]
    record=PeriodicCallback(i->begin
        # Insert a distinguishable stale coast cache. Reentry must clear it.
        if 9.9 < i.t < 10.1
            x.state.targeting_switch_s[1]=84.0
        end
        push!(records,(t=i.t,inside=i.p.shared_buffers.in_atmosphere[1],
            count=x.state.energy_bracketing_count[1],switch=x.state.targeting_switch_s[1]))
    end,0.5)
    crossing=CB.get_drag_state_callback(1)
    callbacks=CB._edg_observe_heat_callbacks(CallbackSet(crossing,CB.get_control_callbacks(1,x.a)...,
        CB.get_guidance_callbacks(1,x.a)...,record),1)
    sol=solve(ODEProblem(radial!,x.u,(0.,19.),x.p),Tsit5();callback=callbacks,
        dtmax=0.25,abstol=1e-8,reltol=1e-10)
    @test string(sol.retcode)=="Success" && last(sol.t)==19.0
    @test any(!r.inside for r in records) && any(r.t>15 && r.inside for r in records)
    @test all(r->r.switch==42.0 && r.count==1,filter(r->r.t<4.5,records))
    @test all(r->r.switch==Inf,filter(r->6<=r.t<=9,records))
    @test all(r->r.switch==84.0,filter(r->10.5<=r.t<=14.5,records))
    @test all(r->r.switch==Inf,filter(r->r.t>=16,records))
    @test x.state.energy_bracketing_count==[1]
    @test isempty(only(x.state.heat_load_exit_j_cm2)) # Cleared on the new entry.
    @test x.state.heat_load_entry_j_cm2==[fill(7.0,3)]
    @test x.u.sc[1].heat_loads==fill(7.0,3)
end

@testset "EDG outbound start and real t0 guidance preserve reentry invalidation" begin
    x=context()
    seed!(x.state); x.state.targeting_active[1]=false
    boundary=x.a.environment_model.planet.Rp_e+x.a.environment_model.EI*1e3
    x.u.sc[1].pos .= (boundary,0,0)
    x.u.sc[1].vel .= (100pi/10.6,0,0)
    thruster=S.BaseThrusterModel(thrust=[0.0],direction=[0.0],Δv=[0.0],
        start_burn_time=[Inf],stop_burn_time=[Inf],Isp=[300.0])
    a=SC._with_configuration(x.a;control_model=S.ControlModel((x.c,thruster),[0.5,1.0]))
    p=S.ODEParams(n_sats=1,args=a)
    p.shared_buffers.et_start[]=x.p.shared_buffers.et_start[]
    E._initialize_in_atmosphere_flags!(p,x.u)
    @test p.shared_buffers.in_atmosphere == [false]
    x.state.heat_boundary_inside[1]=false
    function radial!(du,u,p,t)
        du .= 0
        du.sc[1].pos[1]=u.sc[1].vel[1]
        du.sc[1].vel[1]=-100*(pi/10.6)^2*sin(pi*t/10.6)
    end
    records=NamedTuple[]
    record=PeriodicCallback(i->push!(records,(t=i.t,inside=i.p.shared_buffers.in_atmosphere[1],
        switch=x.state.targeting_switch_s[1])),0.5)
    crossing=CB.get_drag_state_callback(1)
    callbacks=CB._edg_observe_heat_callbacks(CallbackSet(crossing,CB.get_control_callbacks(1,a)...,
        CB.get_guidance_callbacks(1,a)...,record),1)
    integrator=init(ODEProblem(radial!,x.u,(0.,19.),p),Tsit5();callback=callbacks,
        dtmax=0.25,abstol=1e-8,reltol=1e-10)
    @test x.state.selected_mode == [:max_energy_depletion] # Real thruster t0 guidance ran.
    x.state.targeting_switch_s[1]=42.0
    solve!(integrator)
    @test string(integrator.sol.retcode)=="Success" && integrator.t==19.0
    coast=filter(r->r.t<10.5,records)
    nextpass=filter(r->r.t>=11,records)
    @test !isempty(coast) && all(r->!r.inside && r.switch==42.0,coast)
    @test !isempty(nextpass) && all(r->r.inside && r.switch==Inf,nextpass)
    @test x.state.heat_load_entry_j_cm2==[fill(7.0,3)]
    @test x.u.sc[1].heat_loads==fill(7.0,3)
end

@testset "Geodetic observer direction, deduplication and first-entry baseline" begin
    function place_geodetic!(x, i, latitude, height, vertical, northward, t)
        planet=x.a.environment_model.planet
        e2=1-(planet.Rp_p/planet.Rp_e)^2
        n=planet.Rp_e/sqrt(1-e2*sin(latitude)^2)
        position=SVector((n+height)*cos(latitude),0.,(n*(1-e2)+height)*sin(latitude))
        up=SVector(cos(latitude),0.,sin(latitude))
        north=SVector(-sin(latitude),0.,cos(latitude))
        rotation=S.EDGServices._edg_planet_frame_lpi(x.p,t)
        x.u.sc[i].pos .= rotation'*position
        x.u.sc[i].vel .= rotation'*(vertical*up+northward*north+LinearAlgebra.cross(planet.ω,position))
    end
    # At mid-latitude, ellipsoid-normal and radial motion have opposite signs.
    # Raw-side reconciliation is sampled after a resolved crossing, outside
    # transform roundoff. Exact registered roots are covered separately.
    # Reset only the crossing members.
    for latitude in (-pi/4,pi/4), direction in (-1,1), t in (0.,1e4)
        x=context(n=3)
        for i in 1:3;seed!(x.state,i);end
        x.state.heat_boundary_inside .= direction>0
        ei=1e3*x.a.environment_model.EI
        place_geodetic!(x,1,0.,ei+direction*0.01,direction*5.,0.,t)
        place_geodetic!(x,2,latitude,ei+direction*0.01,direction*5.,sign(latitude)*direction*5000.,t)
        place_geodetic!(x,3,latitude,ei+(direction>0 ? -100. : 100.),0.,0.,t)
        second=x.u.sc[2]
        @test sign(LinearAlgebra.dot(second.pos,second.vel))==-direction
        # Central difference of the actual root independently verifies direction.
        before=deepcopy(x.u);after=deepcopy(x.u);dt=0.001
        before.sc[2].pos .-= dt*second.vel;after.sc[2].pos .+= dt*second.vel
        rate=(S.ControlHooks._edg_heat_boundary_distance(after,x.p,t+dt,2)-S.ControlHooks._edg_heat_boundary_distance(before,x.p,t-dt,2))/(2dt)
        @test isapprox(rate,direction*5.;atol=1e-5,rtol=0.)
        physical=copy(x.u);third=snapshot(x.state,3)
        @test all(idx -> direction * CB._edg_observer_distance(x.u,x.p,t,idx) > 0.001, 1:2)
        observe_members!(x;t=t)
        @test x.state.targeting_switch_s[1:2]==[Inf,Inf]
        @test x.state.heat_boundary_inside[1:2]==fill(direction<0,2)
        @test isequal(snapshot(x.state,3),third)
        @test x.u==physical
        seed!(x.state,2)
        second_state=snapshot(x.state,2)
        observe_members!(x;t=t)
        @test isequal(snapshot(x.state,2),second_state) # Late explicit duplicate.
    end
    # Body-fixed stationarity can leave tiny positive transform residue.
    # Distinguish that tangent from physically resolved very slow crossings.
    for latitude in (-pi/4,pi/4), t in (0.,1e4), vertical in (-1e-8,0.,1e-8)
        x=context();seed!(x.state)
        S.ControlHooks._edg_initialize_heat_accounting!(x.a,x.u,0.)
        place_geodetic!(x,1,latitude,1e3*x.a.environment_model.EI,vertical,0.,t)
        physical=copy(x.u)
        rate=S.ControlHooks._edg_heat_boundary_vertical_speed(x.u,x.p,t,1)
        @test vertical==0. ? rate==0. : sign(rate)==sign(vertical)
        initialize_observer_side!(x;t=t)
        @test x.state.heat_boundary_inside[1]==(vertical<=0.)
        @test x.state.targeting_switch_s[1]==42.
        @test x.state.heat_load_entry_j_cm2==[zeros(3)]
        @test x.u==physical
    end
    # Initialize outside the geodetic boundary while inside the spherical one.
    # The first omitted entry must capture each panel's current physical heat.
    x=context(n=2);initialize!(x)
    S.ControlHooks._edg_initialize_heat_accounting!(x.a,x.u,0.)
    ei=1e3*x.a.environment_model.EI
    for (i,lat) in enumerate((0.,pi/4));place_geodetic!(x,i,lat,ei+1.,-5.,0.,0.);end
    E._initialize_in_atmosphere_flags!(x.p,x.u)
    @test x.p.shared_buffers.in_atmosphere[2]
    physical=copy(x.u)
    initialize_observer_side!(x;t=0.)
    @test x.state.heat_boundary_inside==[false,false]
    @test x.state.heat_load_entry_j_cm2==[zeros(3),zeros(3)]
    @test x.u==physical
    for (i,lat) in enumerate((0.,pi/4))
        place_geodetic!(x,i,lat,ei-0.01,-5.,0.,0.2)
        x.u.sc[i].heat_loads .= [7.2,8.2,9.2]
    end
    physical=copy(x.u)
    @test all(idx -> CB._edg_observer_distance(x.u,x.p,0.2,idx) < -0.001, 1:2)
    observe_members!(x;t=0.2)
    @test x.state.heat_load_entry_j_cm2==[[7.2,8.2,9.2],[7.2,8.2,9.2]]
    @test x.state.heat_boundary_inside==[true,true]
    @test x.u==physical
    observe_members!(x;t=0.2)
    @test x.state.heat_load_entry_j_cm2==[[7.2,8.2,9.2],[7.2,8.2,9.2]]
    # Initialization records side without treating a fresh in-pass start as entry
    # or discarding an already computed t0 guidance plan.
    for (height_offset,vertical,inside) in ((-100.,-5.,true),(0.,5.,false),(100.,-5.,false))
        x=context();seed!(x.state)
        S.ControlHooks._edg_initialize_heat_accounting!(x.a,x.u,0.)
        place_geodetic!(x,1,pi/4,1e3*x.a.environment_model.EI+height_offset,vertical,0.,0.)
        physical=copy(x.u)
        initialize_observer_side!(x;t=0.)
        @test x.state.heat_boundary_inside[1]===inside
        @test x.state.targeting_switch_s[1]==42.
        @test x.state.heat_load_entry_j_cm2==[zeros(3)]
        @test x.u==physical
    end
end
@testset "Dense observer preserves per-spacecraft ownership" begin
    x=context(n=3)
    boundary=x.a.environment_model.planet.Rp_e + 1e3*x.a.environment_model.EI
    for (idx,offset,speed) in ((1,-20000.,1000.),(2,20000.,-1000.),(3,-20000.,0.))
        x.u.sc[idx].pos .= [boundary+offset,0.,0.]
        x.u.sc[idx].vel .= [speed,0.,0.]
    end
    S.ControlHooks._edg_initialize_heat_accounting!(x.a,x.u,0.)
    for idx in 1:3;seed!(x.state,idx);end
    function rhs!(du,u,p,t)
        fill!(du,0.)
        for idx in 1:3
            du.sc[idx].pos .=u.sc[idx].vel
            du.sc[idx].heat_loads .=[0.,Float64(idx),2.0idx]
        end
    end
    sol=solve(ODEProblem(rhs!,copy(x.u),(0.,30.),x.p),Tsit5();
        callback=CB._edg_observe_heat_callbacks(CallbackSet(),3),dtmax=.17,reltol=1e-10,abstol=1e-10)
    @test sol.retcode==CB.DiffEqBase.SciMLBase.ReturnCode.Success
    @test x.state.heat_boundary_inside==[false,true,true]
    @test isapprox(x.state.heat_boundary_exit_time_s[1],20.;atol=1e-8)
    @test isapprox(x.state.heat_boundary_entry_time_s[2],20.;atol=1e-8)
    @test isapprox(x.state.heat_load_exit_j_cm2[1],[7.,27.,47.];atol=1e-8)
    @test isapprox(x.state.heat_load_entry_j_cm2[2],[7.,47.,87.];atol=1e-8)
    @test isnan(x.state.heat_boundary_entry_time_s[1])
    @test isnan(x.state.heat_boundary_exit_time_s[2])
    @test isnan(x.state.heat_boundary_entry_time_s[3]) && isnan(x.state.heat_boundary_exit_time_s[3])
    @test x.state.targeting_switch_s==[Inf,Inf,42.]
    @test x.state.heat_load_entry_j_cm2[3]==zeros(3)
    @test isempty(x.state.heat_load_exit_j_cm2[3])
end

end
