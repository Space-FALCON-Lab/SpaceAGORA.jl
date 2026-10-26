module EDGHeatObserverRegression
using SpaceAGORA, Test, ComponentArrays, StaticArrays, LinearAlgebra, OrdinaryDiffEq, DiffEqCallbacks
const REPO_ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
const SimulationModel=SpaceAGORA.SimulationModel
const SimulationEngine=SpaceAGORA.SimulationEngine
const CB=SimulationModel.SimulationCallbacks
const CH=SimulationModel.ControlHooks
include_string(@__MODULE__,first(split(read(joinpath(REPO_ROOT,"test/gnc/aerobraking/energy_depletion_gnc_tests.jl"),String),"\n@testset")))
const _edg_observe_heat_callbacks=CB._edg_observe_heat_callbacks
function linear_context(;inbound=false)
 c=_edg_test_context(heat_load_limit_j_cm2=30.)
 c.u.sc[1].pos .= [c.args.environment_model.planet.Rp_e+(inbound ? 180000. : 140000.),0,0]
 c.u.sc[1].vel .= [inbound ? -1000. : 1000.,0,0]
 c.u.sc[1].heat_loads .= 0.
 CH._edg_initialize_heat_accounting!(c.args,c.u,0.)
 c
end
function linear_rhs!(du,u,p,t)
 fill!(du,0.);du.sc[1].pos .=u.sc[1].vel;du.sc[1].heat_loads .=[0.,1.,2.]
end
@testset "EDG observer preserves original callback behavior" begin
 @testset "Existing external and spherical callback winner preserved" begin
  for drag in (false,true), inbound in (false,true), scaled in (false,true)
   outputs=[]
   for enabled in (false,true)
    c=linear_context(;inbound);calls=Ref(0)
    cond=(u,t,i)->(norm(u.sc[1].pos)-i.p.args.environment_model.planet.Rp_e)/(scaled ? 1000. : 1.)-(scaled ? 160. : 160000.)
    effect=i->begin calls[]+=1;terminate!(i);end
    extra=ContinuousCallback(cond,inbound ? nothing : effect,inbound ? effect : nothing;rootfind=CB.DiffEqBase.SciMLBase.RightRootFind)
    raw=drag ? CallbackSet(CB.get_drag_state_callback(1),extra) : CallbackSet(extra)
    actual=enabled ? _edg_observe_heat_callbacks(raw,1) : raw
    sol=solve(ODEProblem(linear_rhs!,copy(c.u),(0.,40.),c.p),Tsit5();callback=actual,dtmax=.17,reltol=1e-10,abstol=1e-10)
    push!(outputs,(time=sol.t[end],u=Array(sol.u[end]),retcode=sol.retcode,calls=calls[]))
    if enabled
     @test c.state.heat_boundary_inside[1]==inbound
     @test isapprox(inbound ? c.state.heat_load_entry_j_cm2[1] : c.state.heat_load_exit_j_cm2[1],[0.,20.,40.];atol=1e-8)
    end
   end
   @test outputs[1]==outputs[2]
  end
 end
 @testset "Two crossings in one accepted interval and repeated passages" begin
  for cubic in (false,true), terminate_at_exit in (false,true)
   c=linear_context();a=c.args.environment_model.planet.Rp_e
   h(t)=cubic ? 100*(t-1)*(t-2)*(t-3) : 100*(t-1)*(t-2)
   v(t)=cubic ? 100*(3t^2-12t+11) : 100*(2t-3)
   accel(t)=cubic ? 100*(6t-12) : 200.
   c.u.sc[1].pos .=[a+160000+h(0.),0,0];c.u.sc[1].vel .=[v(0.),0,0]
   function rhs!(du,u,p,t)
    fill!(du,0.);du.sc[1].pos .=u.sc[1].vel;du.sc[1].vel[1]=accel(t);du.sc[1].heat_loads .=[0.,1.,2.]
   end
   finish=cubic ? 3. : 2.;captured=Ref(Float64[])
   condition=CB.EDGGeodeticBoundaryCondition((u,t,i)->CH._edg_heat_boundary_distance(u,i.p,Float64(t),1),1,160000.)
   ext=ContinuousCallback(condition,i->begin
    i.t < finish-1e-8 && return nothing
    captured[]=copy(c.state.heat_load_exit_j_cm2[1]);terminate!(i)
   end,nothing)
   raw=terminate_at_exit ? CallbackSet(ext) : CallbackSet()
   sol=solve(ODEProblem(rhs!,copy(c.u),(0.,4.),c.p),Tsit5();callback=_edg_observe_heat_callbacks(raw,1),adaptive=false,dt=terminate_at_exit ? .25 : 4.,abstol=1e-10,reltol=1e-10)
   @test isapprox(c.state.heat_boundary_entry_time_s[1],cubic ? 2. : 1.;atol=1e-8)
   @test isapprox(c.state.heat_boundary_exit_time_s[1],finish;atol=1e-8)
   @test isapprox(c.state.heat_load_entry_j_cm2[1],[0.,cubic ? 2. : 1.,cubic ? 4. : 2.];atol=1e-8)
   @test isapprox(c.state.heat_load_exit_j_cm2[1],[0.,finish,2finish];atol=1e-8)
   @test c.state.heat_boundary_observed_at_s[1] >= c.state.heat_boundary_exit_time_s[1]
   terminate_at_exit && (@test captured[]==c.state.heat_load_exit_j_cm2[1])
  end
 end
 @testset "Vector masks, nested callback sets and field preservation" begin
  for inbound in (false,true)
   c=linear_context(;inbound);masks=[];init_count=Ref(0);final_count=Ref(0)
   cond! = (out,u,t,i)->(out .= (inbound ? 20-t : t-20))
   effect! = (i,mask)->begin push!(masks,copy(mask));terminate!(i);end
   ext=VectorContinuousCallback(cond!,effect!,2;rootfind=CB.DiffEqBase.SciMLBase.RightRootFind,interp_points=7,save_positions=(false,true),idxs=1:3,dtrelax=.75,abstol=1e-9,reltol=0.,repeat_nudge=1//50,maybe_discontinuity=false,initialize_save_discretes=false,
    initialize=(cb,u,t,i)->(init_count[]+=1),finalize=(cb,u,t,i)->(final_count[]+=1))
   recorder=DiscreteCallback((u,t,i)->true,i->nothing;save_positions=(false,false))
   raw=CallbackSet(recorder,CallbackSet(ext))
   wrapped=_edg_observe_heat_callbacks(raw,1)
   for name in fieldnames(typeof(ext))
    name in (:affect!,:initialize,:finalize) && continue
    @test isequal(getfield(ext,name),getfield(wrapped.continuous_callbacks[1],name))
   end
   @test wrapped.discrete_callbacks[2].condition===recorder.condition
   @test wrapped.discrete_callbacks[2].save_positions==recorder.save_positions
   sol=solve(ODEProblem(linear_rhs!,copy(c.u),(0.,40.),c.p),Tsit5();callback=wrapped,dtmax=.17)
   @test sol.retcode==CB.DiffEqBase.SciMLBase.ReturnCode.Terminated
   @test masks==[[inbound ? -1 : 1,inbound ? -1 : 1]]
   @test init_count[]==final_count[]==1
   @test isapprox(inbound ? c.state.heat_load_entry_j_cm2[1] : c.state.heat_load_exit_j_cm2[1],[0.,20.,40.];atol=1e-8)
  end
 end
 @testset "SavingCallback and lifecycle hooks retain their original callback" begin
  outputs=[]
  for enabled in (false,true)
   c=linear_context()
   saved=SavedValues(Float64,Float64)
   saving=SavingCallback((u,t,i)->u.sc[1].heat_loads[2],saved;saveat=0.:2.:40.,save_start=true)
   periodic_times=Float64[]
   periodic=PeriodicCallback(i->push!(periodic_times,i.t),1/3;initial_affect=true)
   lifecycle_calls=Symbol[]
   owner=Ref{Any}()
   init_hook=(cb,u,t,i)->begin
    @test cb===owner[]
    push!(lifecycle_calls,:initialize)
   end
   final_hook=(cb,u,t,i)->begin
    @test cb===owner[]
    push!(lifecycle_calls,:finalize)
   end
   owner[]=DiscreteCallback((u,t,i)->false,i->nothing;initialize=init_hook,finalize=final_hook,save_positions=(false,false))
   raw=CallbackSet(periodic,saving,owner[])
   callbacks=enabled ? _edg_observe_heat_callbacks(raw,1) : raw
   integrator=init(ODEProblem(linear_rhs!,copy(c.u),(0.,40.),c.p),Tsit5();callback=callbacks,dtmax=.17,reltol=1e-10,abstol=1e-10)
   attempt_results=[]
   for attempt in 1:2
    if attempt==2
     empty!(periodic_times)
     CH._edg_initialize_heat_accounting!(c.args,c.u,0.)
     reinit!(integrator,copy(c.u);t0=0.,tf=40.,reinit_callbacks=true,reset_dt=true)
    end
    sol=solve!(integrator)
    @test saved.t==collect(0.:2.:40.)
    @test saved.saveval≈saved.t
    @test periodic_times[1]==0.
    @test all(isapprox.(diff(periodic_times),1/3;atol=1e-12))
    push!(attempt_results,(times=copy(saved.t),values=copy(saved.saveval),periodic=copy(periodic_times),u=Array(sol.u[end]),retcode=sol.retcode))
   end
   @test lifecycle_calls==[:initialize,:finalize,:initialize,:finalize]
   push!(outputs,attempt_results)
  end
  @test outputs[1]==outputs[2]
 end
 @testset "Reinitialization resets observer cursor" begin
  for has_continuous in (false,true)
   c=linear_context(); ext=ContinuousCallback((u,t,i)->t-20.,i->terminate!(i),nothing)
   raw=has_continuous ? CallbackSet(ext) : CallbackSet()
   wrapped=_edg_observe_heat_callbacks(raw,1)
   integrator=init(ODEProblem(linear_rhs!,copy(c.u),(0.,40.),c.p),Tsit5();callback=wrapped,dtmax=.17)
   for attempt in 1:2
    if attempt==2
     CH._edg_initialize_heat_accounting!(c.args,c.u,0.)
     reinit!(integrator,copy(c.u);t0=0.,tf=40.,reinit_callbacks=true,reset_dt=true)
    end
    solve!(integrator)
    @test isnan(c.state.heat_boundary_entry_time_s[1])
    @test isapprox(c.state.heat_boundary_exit_time_s[1],20.;atol=1e-8)
    @test isapprox(c.state.heat_load_exit_j_cm2[1],[0.,20.,40.];atol=1e-8)
   end
  end
 end
 @testset "Observer preserves RHS work and adaptive partition" begin
  runs=[]
  for enabled in (false,true)
   c=linear_context()
   raw=CallbackSet()
   sol=solve(ODEProblem(linear_rhs!,copy(c.u),(0.,40.),c.p),Tsit5();callback=enabled ? _edg_observe_heat_callbacks(raw,1) : raw,dtmax=.17,reltol=1e-10,abstol=1e-10)
   push!(runs,(t=sol.t,u=Array(sol.u[end]),nf=sol.stats.nf,accepted=sol.stats.naccept,rejected=sol.stats.nreject))
  end
  @test runs[1]==runs[2]
 end
 @testset "Discrete mutation is visible to same-tick consumer" begin
  c=linear_context(;inbound=true);seen=Ref(false)
  mutate=DiscreteCallback((u,t,i)->t==1.,i->begin i.u.sc[1].pos[1]-=40000.;i.u.sc[1].heat_loads .+=50.;end;save_positions=(false,false))
  consume=DiscreteCallback((u,t,i)->t==1.,i->begin
   @test c.state.heat_boundary_inside[1]===true
   @test c.state.heat_load_entry_j_cm2[1]==i.u.sc[1].heat_loads
   @test isapprox(c.state.heat_load_entry_j_cm2[1],[50.,51.,52.];atol=1e-10)
   seen[]=true;terminate!(i)
  end;save_positions=(false,false))
  solve(ODEProblem(linear_rhs!,copy(c.u),(0.,2.),c.p),Tsit5();callback=_edg_observe_heat_callbacks(CallbackSet(mutate,consume),1),tstops=[1.],dtmax=.17)
  @test seen[]
 end
 @testset "Exact-root reversal precedes same-time consumer" begin
  c=linear_context();seen=Ref(false);exit=Ref(Float64[])
  reverse=ContinuousCallback((u,t,i)->t-20.,i->begin exit[]=copy(c.state.heat_load_exit_j_cm2[1]);i.u.sc[1].vel[1]*=-1.;end,nothing)
  consume=DiscreteCallback((u,t,i)->t==20.,i->begin
   @test c.state.heat_boundary_inside[1]===true
   @test isapprox(exit[],[0.,20.,40.];atol=1e-8)
   @test c.state.heat_load_entry_j_cm2[1]==i.u.sc[1].heat_loads
   seen[]=true;terminate!(i)
  end;save_positions=(false,false))
  solve(ODEProblem(linear_rhs!,copy(c.u),(0.,40.),c.p),Tsit5();callback=_edg_observe_heat_callbacks(CallbackSet(reverse,consume),1),tstops=[20.],dtmax=.17)
  @test seen[]
 end
 @testset "Initializer changes define initial side without a passage" begin
  for discrete in (false,true)
   c=linear_context();c.state.targeting_switch_s[1]=42.
   first=ContinuousCallback((u,t,i)->t-100.,i->nothing,nothing;
    initialize=(cb,u,t,i)->begin i.u.sc[1].pos[1]+=40000.;cb.affect!(i);end)
   later_init=(cb,u,t,i)->begin
    @test c.state.targeting_switch_s[1]==42.
    @test c.state.heat_boundary_inside[1]===false
    @test isnan(c.state.heat_boundary_entry_time_s[1]) && isnan(c.state.heat_boundary_exit_time_s[1])
    @test c.state.heat_load_entry_j_cm2[1]==zeros(3)
   end
   later=discrete ? DiscreteCallback((u,t,i)->false,i->nothing;initialize=later_init) : ContinuousCallback((u,t,i)->t-101.,i->nothing,nothing;initialize=later_init)
   integrator=init(ODEProblem(linear_rhs!,copy(c.u),(0.,1.),c.p),Tsit5();callback=_edg_observe_heat_callbacks(CallbackSet(first,later),1),dtmax=.17)
   @test c.state.targeting_switch_s[1]==42.
   @test isnan(c.state.heat_boundary_entry_time_s[1]) && isnan(c.state.heat_boundary_exit_time_s[1])
   solve!(integrator)
   @test isnan(c.state.heat_boundary_entry_time_s[1]) && isnan(c.state.heat_boundary_exit_time_s[1])
  end
 end

end

@testset "Registered stop and pre-effect heat at equivalent roots" begin
for lat in (0.,45.,60.), form in (:heat,:frame,:scaled), right in (false,true), mutate in (false,true)
 c=_edg_test_context(heat_load_limit_j_cm2=30.)
 planet=c.args.environment_model.planet
 a=planet.Rp_e;b=planet.Rp_p;e2=1-b*b/(a*a);n=a/sqrt(1-e2*sind(lat)^2)
 c.u.sc[1].pos .= [(n+140000)*cosd(lat),0,(n*(1-e2)+140000)*sind(lat)]
 v=[1000*cosd(lat),0,1000*sind(lat)];c.u.sc[1].vel .=v
 c.u.sc[1].heat_loads .= [0.,0.,0.]
 CH._edg_initialize_heat_accounting!(c.args,c.u,0.)
 count=Ref(0);captured=Ref(Float64[])
 condition = form==:heat ? ((u,t,i)->CH._edg_heat_boundary_distance(u,i.p,Float64(t),1)) :
   form==:frame ? ((u,t,i)->SimulationEngine.sample_planet_frame(u.sc[1],i.p,1,Float64(t)).alt_m-160000.) :
   ((u,t,i)->SimulationEngine.sample_planet_frame(u.sc[1],i.p,1,Float64(t)).alt_m/1000-160.)
 extra=ContinuousCallback(CB.EDGGeodeticBoundaryCondition(condition,1,160000.),i->begin
  count[]+=1;captured[]=copy(c.state.heat_load_exit_j_cm2[1]);mutate && (i.u.sc[1].heat_loads .+= 100.);terminate!(i)
 end,nothing;rootfind=(right ? CB.DiffEqBase.SciMLBase.RightRootFind : CB.DiffEqBase.SciMLBase.LeftRootFind),abstol=1e-6,reltol=0.)
 function rhs!(du,u,p,t)
  fill!(du,0.);du.sc[1].pos .=v;du.sc[1].heat_loads .=[0.,1.,2.]
 end
 sol=solve(ODEProblem(rhs!,copy(c.u),(0.,40.),c.p),Tsit5();callback=_edg_observe_heat_callbacks(CallbackSet(extra),1),dtmax=0.17,reltol=1e-10,abstol=1e-10)
 ok=sol.retcode==CB.DiffEqBase.SciMLBase.ReturnCode.Terminated && count[]==1 && length(captured[])==3 && isapprox(captured[][2],20.;atol=1e-8) && isapprox(captured[][3],40.;atol=1e-8)
 @test ok
 @test abs(sol.t[end]-20.)<1e-8
 @test isapprox(sol.u[end].sc[1].heat_loads,[0.,20.,40.].+(mutate ? 100. : 0.);atol=1e-8)
end
end

# Exact manufactured interpolation isolates observer geometry from the tiny
# excursions that an ODE solver can introduce in an analytic tangent orbit.
mutable struct ObserverWitness{P,F}
 p::P
 trajectory::F
 u::Any
 tprev::Float64
 t::Float64
 derivative_changed::Bool
end
(w::ObserverWitness)(t)=w.trajectory(t)
CB.DiffEqBase.derivative_discontinuity!(w::ObserverWitness, value::Bool)=(w.derivative_changed=value)
function witness(c,trajectory,raw)
 w=ObserverWitness(c.p,trajectory,trajectory(0.),0.,0.,false)
 wrapped=_edg_observe_heat_callbacks(raw,1)
 for cb in (wrapped.continuous_callbacks...,wrapped.discrete_callbacks...)
  cb.initialize(cb,w.u,0.,w)
 end
 return w,wrapped
end
function witness_step!(w,wrapped,t)
 w.tprev=w.t;w.t=Float64(t);w.u=w.trajectory(t)
 for cb in wrapped.discrete_callbacks
  cb.condition(w.u,w.t,w) && cb.affect!(w)
 end
end
@testset "Exactly represented outside tangent does not invent a passage" begin
 for noop in (false,true)
  c=linear_context(;inbound=true);a=c.args.environment_model.planet.Rp_e
  c.state.targeting_switch_s[1]=42.
  trajectory=t->begin
   u=copy(c.u);u.sc[1].pos .=[a+160000. + 100*(t-1)^2,0.,0.]
   u.sc[1].vel .=[200*(t-1),0.,0.];u.sc[1].heat_loads .=[0.,t,2t]
   # Make the evaluated geometry nonnegative, including frame/norm rounding.
   # This fixture promises an actual computed tangent, not an analytic ideal
   # whose rounded frame transform can dip below the surface.
   for _ in 1:8
    CB._edg_observer_distance(u,c.p,Float64(t),1) >= 0.0 && break
    u.sc[1].pos[1]=nextfloat(u.sc[1].pos[1])
   end
   @assert CB._edg_observer_distance(u,c.p,Float64(t),1) >= 0.0
   u
  end
  raw=noop ? CallbackSet(DiscreteCallback((u,t,i)->t==1.,i->nothing;save_positions=(false,false))) : CallbackSet()
  w,wrapped=witness(c,trajectory,raw)
  @test CB._edg_observer_distance(trajectory(1.),c.p,1.,1)==0.
  witness_step!(w,wrapped,1.);witness_step!(w,wrapped,2.)
  @test isnan(c.state.heat_boundary_entry_time_s[1])
  @test isnan(c.state.heat_boundary_exit_time_s[1])
  @test c.state.heat_load_entry_j_cm2[1]==zeros(3)
  @test isempty(c.state.heat_load_exit_j_cm2[1])
  @test c.state.targeting_switch_s[1]==42.
  @test c.state.heat_boundary_inside[1]===false
 end
end
@testset "Near-surface generic events do not create an early crossing" begin
 for reverse in (false,true)
  c=linear_context(;inbound=true);a=c.args.environment_model.planet.Rp_e
  c.state.targeting_switch_s[1]=42.;unit=eps(a)
  trajectory=t->begin
   distance=t < 1. ? 1.0-(1.0-7unit)*t : (8-t)*unit
   speed=t < 1. ? -(1.0-7unit) : -unit
   u=copy(c.u);u.sc[1].pos .=[a+160000. + distance,0.,0.]
   u.sc[1].vel .=[speed,0.,0.];u.sc[1].heat_loads .=[0.,t,2t];u
  end
  effect=i->(reverse && (i.u.sc[1].vel[1]*=-1.);nothing)
  raw=CallbackSet(ContinuousCallback((u,t,i)->t-1.,effect,nothing))
  w,wrapped=witness(c,trajectory,raw)
  w.tprev=0.;w.t=1.;w.u=trajectory(1.)
  @test 0 < CB._edg_observer_distance(w.u,c.p,1.,1) < 64eps(Float64)*norm(w.u.sc[1].pos)
  wrapped.continuous_callbacks[1].affect!(w)
  @test c.state.heat_boundary_inside[1]===false
  @test isnan(c.state.heat_boundary_entry_time_s[1])
  @test c.state.targeting_switch_s[1]==42.
  # Continue the prescribed inward trajectory; an actual raw root is captured
  # near eight seconds, never at the unrelated one-second callback.
  witness_step!(w,wrapped,10.)
  @test 7. < c.state.heat_boundary_entry_time_s[1] < 9.
  @test isapprox(c.state.heat_load_entry_j_cm2[1][2],c.state.heat_boundary_entry_time_s[1];atol=1e-12)
  @test c.state.heat_boundary_observed_at_s[1]==10.
 end
end
@testset "Registered continuation preserves event snapshots" begin
 for inbound in (false,true), reverse in (false,true), mutate_heat in (false,true)
  c=linear_context(;inbound);a=c.args.environment_model.planet.Rp_e
  delta=8eps(a);direction=inbound ? -1. : 1.;after_event=Ref(false)
  trajectory=t->begin
   u=copy(c.u)
   travel=after_event[] && reverse ? -(t-1)-delta : t-1-delta
   u.sc[1].pos .=[a+160000. + direction*travel,0.,0.]
   u.sc[1].vel .=[after_event[] && reverse ? -direction : direction,0.,0.]
   u.sc[1].heat_loads .=[0.,t,2t].+(after_event[] && mutate_heat ? 100. : 0.)
   u
  end
  marker=CB.EDGGeodeticBoundaryCondition((u,t,i)->CH._edg_heat_boundary_distance(u,i.p,Float64(t),1),1,160000.)
  effect=i->begin
   reverse && (i.u.sc[1].vel[1]*=-1.)
   mutate_heat && (i.u.sc[1].heat_loads .+=100.)
   after_event[]=true
  end
  w,wrapped=witness(c,trajectory,CallbackSet(ContinuousCallback(marker,effect,effect)))
  w.t=1.;w.u=trajectory(1.)
  (inbound ? wrapped.continuous_callbacks[1].affect_neg! : wrapped.continuous_callbacks[1].affect!)(w)
  snapshot=()->(entry=c.state.heat_boundary_entry_time_s[1],exit=c.state.heat_boundary_exit_time_s[1],
   entry_heat=copy(c.state.heat_load_entry_j_cm2[1]),exit_heat=copy(c.state.heat_load_exit_j_cm2[1]),inside=c.state.heat_boundary_inside[1])
  before=snapshot()
  @test before.inside==(reverse ? !inbound : inbound)
  @test inbound ? before.entry==1. : (reverse ? before.entry==1. : before.exit==1.)
  witness_step!(w,wrapped,2.)
  @test isequal(before,snapshot())
 end
end

@testset "Later same-time effects consume precise registered event identity" begin
 for inbound in (false,true), jump in (false,true)
  c=linear_context(;inbound);a=c.args.environment_model.planet.Rp_e
  delta=8eps(a);direction=inbound ? -1. : 1.;later=Ref(false)
  trajectory=t->begin
   u=copy(c.u)
   travel=later[] ? -(t-1)-delta+(jump ? 1000. : 0.) : t-1-delta
   u.sc[1].pos .=[a+160000. + direction*travel,0.,0.]
   u.sc[1].vel .=[later[] ? -direction : direction,0.,0.]
   u.sc[1].heat_loads .=[0.,t,2t].+(later[] ? 100. : 0.);u
  end
  marker=CB.EDGGeodeticBoundaryCondition((u,t,i)->CH._edg_heat_boundary_distance(u,i.p,Float64(t),1),1,160000.)
  boundary=ContinuousCallback(marker,i->nothing,i->nothing)
  mutate=DiscreteCallback((u,t,i)->t==1.,i->begin
   i.u.sc[1].vel[1]*=-1.;i.u.sc[1].heat_loads .+=100.
   jump && (i.u.sc[1].pos[1]+=1000direction)
   later[]=true
  end;save_positions=(false,false))
  w,wrapped=witness(c,trajectory,CallbackSet(boundary,mutate))
  w.t=1.;w.u=trajectory(1.)
  (inbound ? wrapped.continuous_callbacks[1].affect_neg! : wrapped.continuous_callbacks[1].affect!)(w)
  wrapped.discrete_callbacks[2].affect!(w)
  expected=jump ? inbound : !inbound
  @test c.state.heat_boundary_inside[1]===expected
  heat=expected ? c.state.heat_load_entry_j_cm2[1] : c.state.heat_load_exit_j_cm2[1]
  if jump
   # Position moved to the already declared side: no new passage snapshot.
   @test heat==[0.,1.,2.]
  else
   @test heat==[100.,101.,102.]
  end
  snapshot=copy(heat)
  witness_step!(w,wrapped,2.)
  @test c.state.heat_boundary_inside[1]===expected
  @test (expected ? c.state.heat_load_entry_j_cm2[1] : c.state.heat_load_exit_j_cm2[1])==snapshot
 end
end

@testset "Registered surface validates identity and tangent direction" begin
 c=linear_context()
 for (idx,alt) in ((0,160000.),(2,160000.),(1,159999.),(1,NaN))
  marker=CB.EDGGeodeticBoundaryCondition((u,t,i)->t-20.,idx,alt)
  @test_throws ArgumentError CB._edg_validate_boundary_condition(marker,c.p,1)
 end
 marker=CB.EDGGeodeticBoundaryCondition((u,t,i)->CH._edg_heat_boundary_distance(u,i.p,Float64(t),1),1,160000.)
 @test CB._edg_validate_boundary_condition(marker,c.p,1)===nothing
 @test_throws ArgumentError CB._edg_registered_boundary_side(marker,c.u,c.p,0.)
 c.u.sc[1].pos[1]=c.args.environment_model.planet.Rp_e+160000.
 c.u.sc[1].vel .=0.
 @test CB._edg_registered_boundary_side(marker,c.u,c.p,0.)===nothing
end

end
