module ExternalPropagationLifecycleTests

using Test, LinearAlgebra, StaticArrays, SpaceAGORA, TOML, SHA, Serialization
using SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const EP = SpaceAGORA.ExternalPropagationServices
const LIFE = SpaceAGORA.SimulationLifecycle
const ZERO3 = SVector(0.0,0.0,0.0)
const R0 = SVector(7.0e6,0.0,0.0)
const V0 = SVector(0.0,7500.0,0.0)

# Deliberately ordinary configured control hooks, shared by the normal route and the external
# service. The owner never implements this force itself. gain>0 exercises state feedback.
mutable struct HeldAcceleration <: AbstractControlEffectorModel
    held::Float64
    gain::Float64
    slope::Float64
    ticks::Vector{Float64}
    commands::Vector{Float64}
    observed::Vector{SVector{6,Float64}}
end
HeldAcceleration(;gain=0.0,slope=0.002) = HeldAcceleration(0.01,gain,slope,Float64[],Float64[],SVector{6,Float64}[])
function SpaceAGORA.update_control!(c::HeldAcceleration,u,p,t::Float64,i::Int64)
    x = u.sc[i]
    c.held = 0.01 + c.slope*t - c.gain*x.vel[1]
    push!(c.ticks,t); push!(c.commands,c.held)
    push!(c.observed,SVector{6,Float64}(vcat(x.pos,x.vel)))
    return nothing
end
SpaceAGORA.control_wrench(c::HeldAcceleration,x,p,t) = (SVector(x.mass_kg*c.held,0.0,0.0),ZERO3)
EP.external_stage_safe(::HeldAcceleration) = true
EP.external_stage_key(c::HeldAcceleration) = (kind=:held_acceleration_v1,gain=c.gain,slope=c.slope)
LIFE.checkpoint_state(c::HeldAcceleration) = (held=c.held,ticks=copy(c.ticks),commands=copy(c.commands),observed=copy(c.observed))
function LIFE.restore_checkpoint_state!(c::HeldAcceleration,d)
    c.held=d.held; c.ticks=copy(d.ticks); c.commands=copy(d.commands); c.observed=copy(d.observed)
    return nothing
end

# A control effector that declares no stage safety: configured-force owners refuse it.
struct UnknownControl <: AbstractControlEffectorModel end

struct Owner <: EP.AbstractExternalPropagator
    dt::Float64
    key::Symbol
end
Owner(;dt=1/128,key=:rk4_translation_v1) = Owner(dt,key)
mutable struct Runtime
    dt::Float64
    n::Int
    y::Vector{Float64}
    service::Any
    changes::Vector{Float64}
    evaluations::Int
end
EP.external_spacecraft(::Owner) = [1]
EP.external_step(o::Owner) = o.dt
EP.external_capabilities(::Owner) = (configured_forces=true,checkpoint=true)
EP.external_checkpoint_key(o::Owner) = (kind=o.key,algorithm=:rk4,dt=o.dt)
function EP.external_prepare(o::Owner,args,u0)
    x = EP.external_initial_state(u0,1)
    Runtime(o.dt,0,collect(vcat(x.pos,x.vel)),nothing,Float64[],0)
end
EP.external_bind_force_service!(r::Runtime,service) = (r.service=service; nothing)
function EP.external_control_changed!(r::Runtime,t::Float64)
    @assert abs(EP.external_time(r)-t) <= 16eps(max(abs(t),1.0))
    push!(r.changes,t)
    return nothing
end
EP.external_time(r::Runtime) = r.n*r.dt
function rhs(r,y,t)
    r.evaluations+=1
    vcat(y[4:6],r.service(view(y,1:3),view(y,4:6),t))
end
function step(r,y,t,h)
    k1=rhs(r,y,t); k2=rhs(r,y+h*k1/2,t+h/2)
    k3=rhs(r,y+h*k2/2,t+h/2); k4=rhs(r,y+h*k3,t+h)
    y+(h/6)*(k1+2k2+2k3+k4)
end
function EP.external_sync!(r::Runtime,t::Float64)
    target=floor(Int,t/r.dt+1e-10)
    before=r.n
    while r.n<target
        r.y=step(r,r.y,EP.external_time(r),r.dt);r.n+=1
    end
    return r.n-before
end
function preview(r,t)
    h=t-EP.external_time(r)
    h>=-16eps(max(abs(t),1.0)) || error("owner passed physical query")
    h<=r.dt+16eps(max(abs(t),1.0)) || error("owner lag exceeds step")
    iszero(h) ? copy(r.y) : step(r,r.y,EP.external_time(r),h)
end
function EP.external_state(r::Runtime,t::Float64,slot::Int)
    y=preview(r,t)
    (pos=SVector{3,Float64}(y[1:3]),vel=SVector{3,Float64}(y[4:6]),q=SVector(0.0,0.0,0.0,1.0),ω=ZERO3)
end
function EP.external_acceleration(r::Runtime,t::Float64,slot::Int)
    # Trial engine stages can be ahead of the last accepted owner synchronization.
    # Advance a private preview, never commit a future owner state or tick control.
    y=copy(r.y);tt=EP.external_time(r)
    while t-tt>r.dt
        y=step(r,y,tt,r.dt);tt+=r.dt
    end
    t==tt || (y=step(r,y,tt,t-tt))
    r.service(view(y,1:3),view(y,4:6),t)
end
EP.external_checkpoint_state(r::Runtime,t::Float64) = (version=1,n=r.n,y=copy(r.y),changes=copy(r.changes),evaluations=r.evaluations)
function EP.external_restore_state!(r::Runtime,d,t::Float64)
    d.version==1 || throw(ArgumentError("unsupported test owner snapshot"))
    d.n*r.dt==t || throw(ArgumentError("test owner snapshot clock mismatch"))
    r.n=d.n;r.y=copy(d.y);r.changes=copy(d.changes);r.evaluations=d.evaluations
    return nothing
end

test_directory() = mktempdir()

# The checkpoint a run left on disk, read back: its time, its manifest and its
# external-owner entry (owner state and endpoint handoff at the checkpoint time).
function checkpoint_payload(cfg,expected_time)
    paths=SE._checkpoint_paths(cfg)
    payload=open(deserialize,paths.data)
    parsed=TOML.parsefile(paths.manifest)
    @test payload.t==expected_time
    @test parsed["time_s"]==expected_time
    @test parsed["data_sha256"]==bytes2hex(open(sha256,paths.data))
    owners=payload.runtime_state.external_owners
    @test owners.time_s==expected_time
    @test owners.handoff.time_s==expected_time
    return payload
end

function setup(;owned=true,owner=Owner(),controller=HeldAcceleration(),controllers=nothing,rates=nothing,
               gravity=false,interval=nothing,resume=false,dir=test_directory(),tend=2.0,period=0.25,
               data_rate=0.25,ephemeris=SM.SimpleEphemeridesModel(),initial_pos=R0,initial_vel=V0)
    root=SM.Link(root=true,m=100.0,dims=MVector(1.0,1.0,1.0))
    sc=SM.SpacecraftModel(links=[root],root=root,prop_mass=0.0,inertia_tensor=root.inertia,
        initial_condition=SM.CartesianInitialCondition(initial_pos,initial_vel),id=1)
    effects=gravity===:j2 ? (SM.InverseSquaredJ2GravityModel(),) : gravity===true ? (SM.InverseSquaredGravityModel(),) : ()
    base=make_example_config(planet=SM.make_no_gram_planet(:earth),spacecraft=sc,mission_time=tend,
        initial_time=SM.InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.0),
        dynamic_effectors=effects,density_model=SM.NoAtmosphereModel(),ephemerides_model=ephemeris,
        orientation_sim=false,keplerian=true,verbose=false,results=false,
        solver_config=SM.SolverConfig(solver_mode=:dp8))
    cfg=SM.SimConfig._with_configuration(base;
        simulation_settings=SM.SimulationSettings(results=false,verbose=false,generate_plots=false,normalize=false,
            save_csv=false,results_directory=dir,checkpoint_enabled=interval!==nothing,
            checkpoint_interval_s=something(interval,300.0),checkpoint_directory=joinpath(dir,"checkpoints"),resume_from_checkpoint=resume),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime,keplerian=true,
            mission_time=tend,orientation_sim=false,data_rate=data_rate),
        control_model=SM.ControlModel(control_effectors=controllers===nothing ? (controller,) : controllers,
            control_rates=rates===nothing ? [period] : rates))
    EP.with_external_propagators(cfg,owned ? (owner,) : ())
end

function run(cfg;isolate_state=false)
    last_runtime=Ref{Any}(nothing)
    function factory(args)
        SM.SimulationCallbacks.DiscreteCallback((u,t,int)->true,
            int->(isempty(EP.external_runtimes(int)) || (last_runtime[]=first(EP.external_runtimes(int)));nothing);
            save_positions=(false,false))
    end
    out=run_simulation(cfg;return_results=true,return_solution=true,isolate_state=isolate_state,extra_callback_factories=(factory,))
    (;out,rt=last_runtime[])
end
state_row(tbl,i) = vcat([tbl[i,"sc1_pos_$k"] for k in 1:3],[tbl[i,"sc1_vel_$k"] for k in 1:3])

function held_reference(c,t)
    position=R0[1];velocity=0.0;left=0.0;held=0.01
    for (tick, command) in zip(c.ticks,c.commands)
        tick>t && break
        h=tick-left;position+=velocity*h+held*h^2/2;velocity+=held*h
        left=tick;held=command
    end
    h=t-left
    (position=position+velocity*h+held*h^2/2,velocity=velocity+held*h)
end

@testset "experimental lifecycle contract and force guards" begin
    @test EP.require_version(v"2.0.0")===nothing
    @test EP.require_lifecycle_version(v"0.1.0")===nothing
    @test_throws EP.CompatibilityError EP.require_lifecycle_version(v"9.0.0")
    @test EP.validate_external_propagators(setup(),:dp8)===nothing
    @test_throws ArgumentError EP.validate_external_propagators(setup(interval=0.3),:dp8)
    @test_throws ArgumentError EP.validate_external_propagators(setup(period=0.3),:dp8)
    @test_throws ArgumentError EP.validate_external_propagators(setup(controller=UnknownControl()),:dp8)
    @test_throws ArgumentError EP.validate_external_propagators(setup(ephemeris=SM.SpiceEphemeridesModel()),:dp8)
    hardware=setup();sc=only(hardware.dynamics_model.spacecraft)
    add_actuator!(sc,sc.root,ReactionWheel(spin_axis=SVector(1.0,0.0,0.0),max_torque_nm=1.0,max_momentum_nms=10.0))
    @test_throws ArgumentError EP.validate_external_propagators(hardware,:dp8)
end

@testset "actual configured forces and held control through run_simulation" begin
    for gravity in (false,true), gain in (0.0,0.02)
        c1=HeldAcceleration(;gain);c2=HeldAcceleration(;gain)
        owned=run(setup(controller=c1,gravity=gravity))
        ordinary=run(setup(owned=false,controller=c2,gravity=gravity))
        @test owned.out.table.time==ordinary.out.table.time
        @test c1.ticks==c2.ticks
        @test c1.ticks==owned.rt.changes
        @test maximum(abs.(c1.commands.-c2.commands))<1e-9
        @test owned.rt.evaluations>0
        for i in eachindex(owned.out.table.time)
            d=state_row(owned.out.table,i)-state_row(ordinary.out.table,i)
            @test norm(d[1:3])<1e-4
            @test norm(d[4:6])<1e-7
            if !gravity && gain==0.0
                ref=held_reference(c1,owned.out.table.time[i]);y=state_row(owned.out.table,i)
                @test abs(y[1]-ref.position)<1e-5
                @test abs(y[4]-ref.velocity)<1e-10
            end
        end
        @test owned.rt.service(R0,V0,2.0)[1]≈(gravity ? -SM.make_no_gram_planet(:earth).μ/norm(R0)^2 : 0.0)+c1.held atol=1e-12
        @test all(owned.rt.service(R0,V0,2.0)[2:3].==0)
    end
    # A constant force has an independent analytic answer, not merely route agreement.
    c=HeldAcceleration(slope=0.0)
    r=run(setup(controller=c))
    for i in eachindex(r.out.table.time)
        t=r.out.table.time[i];y=state_row(r.out.table,i)
        @test abs(y[1]-(R0[1]+0.005t^2))<1e-5
        @test abs(y[4]-0.01t)<1e-10
    end
end

@testset "actual owner checkpoint and cold resume" begin
    full_control=HeldAcceleration(gain=0.02)
    full=run(setup(controller=full_control))
    dir=test_directory();first_control=HeldAcceleration(gain=0.02)
    first_cfg=setup(controller=first_control,interval=0.5,dir=dir,tend=1.0)
    first_run=run(first_cfg)
    checkpoint_payload(first_cfg,1.0)
    ck=SE._load_checkpoint(first_cfg)
    owner_y=ck.runtime_state.external_owners.owners[1].state.y
    shadow_y=collect(vcat(ck.u.sc[1].pos,ck.u.sc[1].vel))
    @test ck.t==1.0
    @test ck.runtime_state.external_owners.version==1
    @test ck.runtime_state.external_owners.owners[1].state.n==128
    @test owner_y==shadow_y
    handoff=ck.runtime_state.external_owners.handoff
    @test handoff.time_s==ck.t
    @test handoff.updates[1].pos==collect(ck.u.sc[1].pos)
    @test handoff.updates[1].vel==collect(ck.u.sc[1].vel)
    @test collect(first_run.out.solution.u[end].sc[1].pos)==collect(ck.u.sc[1].pos)
    @test collect(first_run.out.solution.u[end].sc[1].vel)==collect(ck.u.sc[1].vel)
    @test state_row(first_run.out.table,length(first_run.out.table.time))==shadow_y
    # The handoff updates only owned kinematics; it is also the state passed to continuation.
    mock_owned=(pos=collect(ck.u.sc[1].pos).+1.0,vel=collect(ck.u.sc[1].vel).+2.0,mass=123.0,heat=[456.0])
    mock_other=(pos=[1.0,2.0,3.0],vel=[4.0,5.0,6.0],mass=789.0,heat=[987.0])
    mock_u=(sc=[mock_owned,mock_other],)
    mock_before=deepcopy(mock_u)
    SE._capture_external_checkpoint_endpoint!(first_run.rt.service.parameters,mock_u,ck.t)
    @test mock_u.sc[2]==mock_before.sc[2]
    @test mock_u.sc[1].mass==mock_before.sc[1].mass && mock_u.sc[1].heat==mock_before.sc[1].heat
    @test mock_u.sc[1].pos==owner_y[1:3] && mock_u.sc[1].vel==owner_y[4:6]
    untouched=deepcopy(mock_u)
    @test_throws ArgumentError SE._capture_external_checkpoint_endpoint!(first_run.rt.service.parameters,mock_u,ck.t+0.25)
    @test mock_u==untouched
    resumed_control=HeldAcceleration(gain=0.02)
    resumed_cfg=setup(controller=resumed_control,interval=0.5,dir=dir,resume=true)
    resumed=run(resumed_cfg)
    checkpoint_payload(resumed_cfg,2.0)
    @test resumed_control.ticks==full_control.ticks
    @test length(resumed_control.ticks)==length(unique(resumed_control.ticks))
    @test resumed.rt.changes==full.rt.changes
    @test norm(resumed.rt.y-full.rt.y)<1e-8
    @test resumed.rt!==first_run.rt
    @test resumed.rt.y!==first_run.rt.y
    @test resumed.rt.service!==first_run.rt.service
    for i in eachindex(resumed.out.table.time)
        t=resumed.out.table.time[i]
        j=findfirst(==(t),full.out.table.time)
        @test j!==nothing
        j===nothing || @test norm(state_row(resumed.out.table,i)-state_row(full.out.table,j))<1e-8
    end
    @test_throws ArgumentError run(setup(owner=Owner(key=:different_algorithm),controller=HeldAcceleration(gain=0.02),interval=0.5,dir=dir,resume=true))
    @test_throws ArgumentError run(setup(controller=HeldAcceleration(gain=0.03),interval=0.5,dir=dir,resume=true))
    @test_throws ArgumentError run(setup(controller=HeldAcceleration(gain=0.02),interval=0.5,dir=dir,resume=true,
        ephemeris=SM.SimpleEphemeridesModel(prime_meridian_at_reference_rad=0.3)))
    @test_throws ArgumentError SE._check_external_checkpoint_state(resumed_cfg,(t=1.0,u=ck.u,runtime_state=nothing))
    corrupt=merge(ck,(runtime_state=merge(ck.runtime_state,(external_owners=merge(ck.runtime_state.external_owners,(time_s=2.0,)),)),))
    @test_throws ArgumentError SE._check_external_checkpoint_state(resumed_cfg,corrupt)
    entry=ck.runtime_state.external_owners; saved=entry.owners[1]
    wrong_version=merge(entry,(owners=[merge(saved,(state=merge(saved.state,(version=2,)),))],))
    @test_throws ArgumentError SE._restore_external_checkpoint_state!(resumed.rt.service.parameters,wrong_version,ck.u,ck.t)
    altered_y=copy(saved.state.y);altered_y[1]+=1.0
    wrong_state=merge(entry,(owners=[merge(saved,(state=merge(saved.state,(y=altered_y,)),))],))
    @test_throws ArgumentError SE._restore_external_checkpoint_state!(resumed.rt.service.parameters,wrong_state,ck.u,ck.t)
end

@testset "J2 configured force and non-control output epochs" begin
    pos=SVector(6.5e6,1.0e6,2.0e6);vel=SVector(-1000.0,7000.0,1200.0)
    ephem=SM.SimpleEphemeridesModel(prime_meridian_at_reference_rad=0.37)
    cfg=setup(gravity=:j2,initial_pos=pos,initial_vel=vel,ephemeris=ephem,data_rate=0.125)
    ext=run(cfg)
    normal=run(setup(owned=false,gravity=:j2,initial_pos=pos,initial_vel=vel,ephemeris=ephem,data_rate=0.125))
    for i in eachindex(ext.out.table.time)
        @test norm(state_row(ext.out.table,i)-state_row(normal.out.table,i))<1e-4
    end
    p=cfg.environment_model.planet;rr=norm(pos);u=pos[3]^2/rr^2
    central=-p.μ*pos/rr^3
    j2=(1.5*p.J2*p.μ*p.Rp_e^2/rr^5)*SVector(pos[1]*(5u-1),pos[2]*(5u-1),pos[3]*(5u-3))
    c=only(cfg.control_model.control_effectors)
    before=(y=copy(ext.rt.y),n=ext.rt.n,held=c.held,ticks=copy(c.ticks),commands=copy(c.commands),
        mounted=copy(only(cfg.dynamics_model.spacecraft).root.rw_assembly.mounted_rates),
        cmg=copy(only(cfg.dynamics_model.spacecraft).root.rw_assembly.cmg_rates))
    @test ext.rt.service(pos,vel,2.0)≈central+j2+SVector(c.held,0.0,0.0) atol=1e-11 rtol=1e-13
    @test ext.rt.y==before.y && ext.rt.n==before.n
    @test c.held==before.held && c.ticks==before.ticks && c.commands==before.commands
    @test only(cfg.dynamics_model.spacecraft).root.rw_assembly.mounted_rates==before.mounted
    @test only(cfg.dynamics_model.spacecraft).root.rw_assembly.cmg_rates==before.cmg
end

@testset "between-tick restart, multiple controllers and shared-template isolation" begin
    full=run(setup(controller=HeldAcceleration(gain=0.02),data_rate=0.125))
    dir=test_directory()
    stopped_cfg=setup(controller=HeldAcceleration(gain=0.02),interval=0.125,tend=1.125,dir=dir,data_rate=0.125)
    run(stopped_cfg)
    checkpoint_payload(stopped_cfg,1.125)
    resume_cfg=setup(controller=HeldAcceleration(gain=0.02),interval=0.125,resume=true,dir=dir,data_rate=0.125)
    resumed=run(resume_cfg)
    checkpoint_payload(resume_cfg,2.0)
    @test norm(resumed.rt.y-full.rt.y)<1e-8
    @test resumed.rt.changes==full.rt.changes
    cs=(HeldAcceleration(slope=0.001),HeldAcceleration(slope=-0.001))
    ns=(HeldAcceleration(slope=0.001),HeldAcceleration(slope=-0.001))
    multi=run(setup(controllers=cs,rates=[0.25,0.5],data_rate=0.125))
    ordinary=run(setup(owned=false,controllers=ns,rates=[0.25,0.5],data_rate=0.125))
    @test cs[1].ticks==ns[1].ticks && cs[2].ticks==ns[2].ticks
    @test sort(multi.rt.changes)==sort(vcat(cs[1].ticks,cs[2].ticks))
    for i in eachindex(multi.out.table.time)
        @test norm(state_row(multi.out.table,i)-state_row(ordinary.out.table,i))<1e-4
    end
    shared=setup(owner=Owner(),tend=0.5)
    first=run(shared;isolate_state=true);second=run(shared;isolate_state=true)
    @test first.rt!==second.rt && first.rt.y!==second.rt.y && first.rt.service!==second.rt.service
    @test first.rt.y==second.rt.y
    @test isempty(only(shared.control_model.control_effectors).ticks)
end

# Saved output times of configured-force owners are solver stops, so each saved row is the
# synchronized owner state rather than the outer solver's interpolant between stops.
const CB = SM.SimulationCallbacks

# Records integrator.t and the six owned shadow entries at every accepted step, after the sync
# and saving callbacks (extra callbacks are appended last), and the owner's runtime.
function recording_run(cfg)
    times=Float64[];states=Vector{Vector{Float64}}();runtime=Ref{Any}(nothing)
    function record(int)
        push!(times,Float64(int.t));push!(states,vcat(collect(int.u.sc[1].pos),collect(int.u.sc[1].vel)))
        isempty(EP.external_runtimes(int)) || (runtime[]=first(EP.external_runtimes(int)))
        return nothing
    end
    factory(args)=CB.DiscreteCallback((u,t,int)->true,record;save_positions=(false,false))
    out=run_simulation(cfg;return_results=true,return_solution=true,isolate_state=false,extra_callback_factories=(factory,))
    (;out,times,states,rt=runtime[])
end

@testset "configured-force owner output times are solver stops" begin
    @test CB._external_output_stops_required(setup())
    @test !CB._external_output_stops_required(setup(owned=false))
    @test CB.get_external_output_stop_callback(setup(owned=false),nothing)===nothing
    # Outputs every 0.125 s, controls every 0.25 s: half of the output times are not control stops.
    r=recording_run(setup(controller=HeldAcceleration(gain=0.02),data_rate=0.125))
    table_times=Float64.(r.out.table.time)
    @test table_times==collect(range(0.0,2.0;step=0.125))
    for (i,t) in enumerate(table_times)
        t==0.0 && continue
        j=findlast(==(t),r.times)
        @test j!==nothing
        # The saved row is the state the solver holds at that stop, after the owner sync.
        j===nothing || @test state_row(r.out.table,i)==r.states[j]
    end
    # Checkpointed segments end on the same output grid; rows still come from stops.
    c=recording_run(setup(controller=HeldAcceleration(gain=0.02),interval=0.5,tend=2.0,dir=test_directory(),data_rate=0.125))
    @test Float64.(c.out.table.time)==table_times
    for (i,t) in enumerate(table_times)
        t==0.0 && continue
        @test findlast(==(t),c.times)!==nothing
        @test norm(state_row(c.out.table,i)-state_row(r.out.table,i))<1e-8
    end
end

# Every column of two saved tables, compared with `isequal` (bit for bit, NaN equal to NaN).
same_rows(a,b) = propertynames(a)==propertynames(b) &&
    all(isequal(getproperty(a,name),getproperty(b,name)) for name in propertynames(a))

# A run stopped at a checkpoint and resumed from it in a new run saves the rows of the same run
# left uninterrupted, bit for bit, with a configured-force owner: restarted at a control tick
# (1 s) and between two ticks (1.125 s). Outputs every 0.125 s fall between control ticks too.
@testset "checkpointed and resumed run saves the uninterrupted run's rows" begin
    full_control=HeldAcceleration(gain=0.02)
    full=run(setup(controller=full_control,data_rate=0.125))
    @test Float64.(full.out.table.time)==collect(range(0.0,2.0;step=0.125))
    for (stop,interval) in ((1.0,0.5),(1.125,0.375))
        dir=test_directory()
        stopped_cfg=setup(controller=HeldAcceleration(gain=0.02),interval=interval,tend=stop,dir=dir,data_rate=0.125)
        run(stopped_cfg)
        @test SE._load_checkpoint(stopped_cfg).t==stop
        resumed_control=HeldAcceleration(gain=0.02)
        resumed=recording_run(setup(controller=resumed_control,interval=interval,dir=dir,resume=true,data_rate=0.125))
        # The new run resumed at the checkpoint: no step before it.
        @test minimum(resumed.times)>stop
        @test same_rows(resumed.out.table,full.out.table)
        @test resumed_control.ticks==full_control.ticks
        @test resumed_control.commands==full_control.commands
        @test resumed.rt.n==full.rt.n && resumed.rt.y==full.rt.y
    end
end

end # module
