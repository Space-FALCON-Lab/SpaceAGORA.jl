# The planner path's contract extensions: geometry providers with revisions, the
# hold installed under `on_rejection=:install_fallback`, and retiming triggers.
module RPOPlannerExtensionsTests
using Test, SpaceAGORA, Random, LinearAlgebra, StaticArrays
const S=SpaceAGORA.SimulationModel
const E=SpaceAGORA.SimulationEngine
const G=S.GuidanceHooks
const Hooks=SpaceAGORA.SimulationLifecycle

# Plans a two-sample reference from the chaser to the goal and refuses the listed
# reasons. Declares retiming, which plans again. Valid only for a goal at the start.
struct StepPlanner <: AbstractRPOPlanner
    refuse::Tuple{Vararg{Symbol}}
end
SpaceAGORA.planner_capabilities(::StepPlanner)=RPOPlannerCapabilities(
    state_sources=(:truth,),frames=(:target_rtn,),retiming=true)
SpaceAGORA.consumes_rng(::StepPlanner)=false
function SpaceAGORA.plan_rpo!(state,p::StepPlanner,r::RPOPlanningRequest,rng::AbstractRNG)
    r.reason in p.refuse && return RPOPlanningResult(request_id=r.request_id,
        status=:infeasible,termination=:test_refusal)
    ref=RPOReference(t_ref_s=[0.,r.reference_dt_s],
        r_ref_rtn_m=hcat(collect(r.x_rtn[1:3]),collect(r.goal_rtn_m)),v_ref_rtn_mps=zeros(3,2),
        origin_time_s=r.time_s,valid_until_s=r.valid_until_s,chaser_id=r.chaser_id,
        target_id=r.target_id,geometry_revision=r.geometry_revision)
    RPOPlanningResult(request_id=r.request_id,status=:candidate,termination=:completed,reference=ref)
end
SpaceAGORA.retime_rpo!(state,p::StepPlanner,r::RPOPlanningRequest,active::RPOReference,rng::AbstractRNG)=
    plan_rpo!(state,p,r,rng)

# The small HYPR configuration of examples/rpo_planner_env/smoke.jl, in either
# retiming mode. Whether validation accepts a retimed reference under the common
# acceleration limit depends on the mode and the path (see the adapter tests).
hypr(mode=:legacy)=HYPRRPOPlanner(RPOPSOConfig(n_waypoints=2,n_particles=8,n_iters=3,
    adaptive_enable=false,adaptive_n_waypoints_min=2,adaptive_n_waypoints_max=2,
    adaptive_n_particles_min=8,adaptive_n_particles_max=8,
    adaptive_n_iters_min=3,adaptive_n_iters_max=3,sample_ds_m=.2,
    retime_accel_limit_enable=true,retime_max_speed_mps=.5,
    retime_a_max_mps2=.5*.05/5.2,iteration_runtime_limit_s=Inf,rrt_warmstart_enable=true,
    hypr_mode=mode))

function prepare(args)
    owned=deepcopy(args)
    u=E.build_initial_conditions(owned)
    p=S.ODEParams(n_sats=length(owned.dynamics_model.spacecraft),args=owned)
    SpaceAGORA.SimulationModel.SimulationCallbacks.initialize_gnc_links!(owned,u,p,0.)   # links first, as the engine does
    for g in owned.guidance_model.guidance_effectors
        Hooks.initialize_guidance!(g,u,p,0.)
    end
    return owned,u,p,only(owned.guidance_model.guidance_effectors)
end
# Put the chaser at an RTN position, at rest in the rotating frame.
function move_chaser!(u,rtn)
    b=u.sc[2]
    r,v=S.FrameTransforms.rtn_to_inertial_relative_state(SVector{3}(rtn),SVector(0.,0.,0.),
        SVector{3}(b.pos),SVector{3}(b.vel))
    u.sc[1].pos.=r;u.sc[1].vel.=v
    return u
end
tick!(g,u,p,t)=G.update_guidance!(g,u,p,Float64(t),1)
outcome(f)=try f(); nothing catch e; e end
lifecycle(records)=[x.event for x in records if x.event in (:installed,:rejected,:fallback_installed,:trigger)]
triggers(records)=[x for x in records if x.event===:trigger]
station()=S.RPOReferenceGeometry(S.RPOStationGeometry(zeros(3,1);keepout_radius_m=.25);
    chaser=S.RPOCubeSatGeometry(dims_m=(.1,.1,.3)))

@testset "Replanning spheres: revisions and snapshots" begin
    spheres=[(center_rtn=(4.,0.,0.),radius_m=.1,appear_time_s=.2,disappear_time_s=.5,label="a"),
        (center_rtn=(0.,4.,0.),radius_m=.2,appear_time_s=.3),
        (center_rtn=(0.,0.,4.),radius_m=.2,appear_time_s=1.,velocity_rtn_mps=(.1,0.,0.))]
    provider=RPOReplanningSpheres(spheres;sphere_surface_samples=24,observation_window_s=.1)
    @test provider isa AbstractRPOGeometryProvider
    revision(t)=rpo_geometry_revision(provider,"base",t)
    @test revision(0.)=="base" && revision(.2)=="base+spheres:1"
    @test revision(.3)=="base+spheres:1,2" && revision(.55)=="base+spheres:1,2" # observation window
    @test revision(.65)=="base+spheres:2"
    @test revision(1.)=="base+spheres:2,3@1.0" && revision(1.5)=="base+spheres:2,3@1.5" # moving
    base=station()
    @test rpo_geometry_snapshot(provider,base,0.)===base
    # The snapshot is the legacy model's geometry at the same time.
    legacy=G.RPOReplanningConfig(spheres=spheres,sphere_surface_samples=24,observation_window_s=.1)
    for t in (.25,.45,.62,1.2)
        snapshot=rpo_geometry_snapshot(provider,base,t)
        expected=G.rpo_geometry_with_replanning_spheres(base,G.rpo_active_replanning_spheres(legacy,t);
            sphere_surface_samples=24)
        @test snapshot.station.points_body==expected.station.points_body
        @test snapshot.station.name==expected.station.name && snapshot.station.keepout_radius_m==.25
    end
    # Equal revisions give equal snapshots; the base is never modified.
    @test revision(.31)==revision(.49)
    @test rpo_geometry_snapshot(provider,base,.31).station.points_body==
        rpo_geometry_snapshot(provider,base,.49).station.points_body
    @test base.station.points_body==zeros(3,1)
    # Both constructors: a sphere list builds the configuration, and an existing
    # configuration is used as it is, not read as a sphere list.
    @test length(provider.config.spheres)==3 && provider.config.sphere_surface_samples==24
    @test provider.config.observation_window_s==.1
    wrapped=RPOReplanningSpheres(legacy)
    @test wrapped isa RPOReplanningSpheres && wrapped.config===legacy
    @test rpo_geometry_revision(wrapped,"base",1.)==revision(1.)
    @test length(methods(RPOReplanningSpheres))==2
    @test_throws ArgumentError RPOReplanningSpheres([(center_rtn=(1.,0.,0.),radius_m=0.)])
    @test_throws ArgumentError RPOReplanningSpheres(spheres;sphere_surface_samples=4)
    @test_throws ArgumentError RPOReplanningSpheres(spheres;observation_window_s=-1.)
end

@testset "Configuration refusals" begin
    step=StepPlanner(())
    @test_throws ArgumentError make_rpo_configuration(planner=step,geometry_provider=:spheres)
    @test_throws ArgumentError make_rpo_configuration(planner=step,on_rejection=:keep_previous)
    @test_throws ArgumentError make_rpo_configuration(planner=DirectRPOPlanner(),retime_tracking_error_m=.5)
    for value in (0.,-1.,NaN)
        @test_throws ArgumentError make_rpo_configuration(planner=step,retime_tracking_error_m=value)
    end
    @test_throws ArgumentError make_rpo_configuration(planner=step,hold_position_tolerance_m=0.)
    @test_throws ArgumentError make_rpo_configuration(planner=step,hold_speed_tolerance_mps=Inf)
    g=only(make_rpo_configuration(planner=step,on_rejection=:install_fallback,retime_tracking_error_m=.5,
        hold_position_tolerance_m=.1,hold_speed_tolerance_mps=.01).guidance_model.guidance_effectors)
    @test bound_family(g)==RPOReferenceFamily(on_rejection=:install_fallback)
    @test (g.retime_tracking_error_m,g.hold_position_tolerance_m,g.hold_speed_tolerance_mps)==(.5,.1,.01)
    @test bound_family(only(make_rpo_configuration(planner=step).guidance_model.guidance_effectors))==RPOReferenceFamily()
end

@testset "Inactive provider and default options leave a run unchanged" begin
    later=RPOReplanningSpheres([(center_rtn=(0.,0.,9.),radius_m=.1,appear_time_s=100.)])
    for planner in (DirectRPOPlanner(),hypr())
        base=run_simulation(make_rpo_configuration(planner=planner,mission_time_s=.4,replan_interval_s=.2);
            return_solution=true)
        same=run_simulation(make_rpo_configuration(planner=planner,mission_time_s=.4,replan_interval_s=.2,
            geometry_provider=later,on_rejection=:fail,retime_tracking_error_m=Inf);return_solution=true)
        a,b=rpo_run_report(base),rpo_run_report(same)
        @test repr(a)==repr(b)
        @test base.t==same.t && all(x==y for (x,y) in zip(base.u,same.u))
        @test isempty(triggers(only(a.planners).records))
        @test keys(only(a.planners))==(:chaser_id,:target_id,:planner,:settings,:seed,:generator,
            :stream_mapping,:stream_words,:state_source,:frame,:geometry_revision,:constraints,
            :validation_settings,:request_id,:active_reference,:records)
    end
end

@testset "Geometry revisions: kept, replanned and refused" begin
    # A sphere far from the path: the reference is checked and kept.
    far=RPOReplanningSpheres([(center_rtn=(3.,10.,0.),radius_m=.25,appear_time_s=.15)])
    owned,u,p,g=prepare(make_rpo_configuration(planner=DirectRPOPlanner(),geometry_provider=far))
    tick!(g,u,p,.1)
    @test g.runtime.request_id==1 && isempty(triggers(g.runtime.records))
    tick!(g,u,p,.2)
    kept=only(triggers(g.runtime.records))
    @test kept.trigger===:geometry_revision && kept.time_s==.2
    @test kept.details.action===:kept && kept.details.revision=="rpo-pilot-v1+spheres:1"
    @test kept.details.min_clearance_m>=.2
    tick!(g,u,p,.3)
    @test g.runtime.request_id==1 && length(triggers(g.runtime.records))==1 # checked once per revision
    # A sphere on the remaining path: a replan from the displaced chaser clears it.
    block=RPOReplanningSpheres([(center_rtn=(4.,0.,0.),radius_m=.1,appear_time_s=.15)])
    owned,u,p,g=prepare(make_rpo_configuration(planner=DirectRPOPlanner(),geometry_provider=block))
    move_chaser!(u,(3.,3.,0.))
    tick!(g,u,p,.2)
    blocked=only(triggers(g.runtime.records))
    @test blocked.details.action===:replan && blocked.details.min_clearance_m<.2
    request=last(filter(x->x.event===:request,g.runtime.records)).request
    @test request.reason===:replan && request.geometry_revision=="rpo-pilot-v1+spheres:1"
    @test size(request.geometry.station.points_body,2)==1+96
    installed=last(filter(x->x.event===:installed,g.runtime.records))
    @test installed.time_s==.2 && installed.validation.accepted
    @test installed.reference.geometry_revision=="rpo-pilot-v1+spheres:1"
    @test installed.validation.metrics.min_clearance_m>=.2
    @test g.triggers.revision=="rpo-pilot-v1+spheres:1"
    # Without a fallback, a blocked replan stops the run as any rejection does.
    e=outcome(()->run_simulation(make_rpo_configuration(planner=DirectRPOPlanner(),
        geometry_provider=block,mission_time_s=.3)))
    @test e isa RPOPlanningError && e.reason===:infeasible
    @test last(e.report.records).details.termination===:direct_segment_blocked
end

# A provider that fails after a time, by throwing or with an empty revision.
struct BrokenProvider <: AbstractRPOGeometryProvider
    after::Float64
    empty::Bool
end
SpaceAGORA.rpo_geometry_revision(p::BrokenProvider,revision,t)=
    t<=p.after ? revision : p.empty ? "" : error("provider unavailable")
SpaceAGORA.rpo_geometry_snapshot(::BrokenProvider,geometry,t)=geometry

@testset "Provider failures are planning failures" begin
    for (empty,reason) in ((false,:geometry_provider_exception),(true,:invalid_geometry_revision))
        owned,u,p,g=prepare(make_rpo_configuration(planner=DirectRPOPlanner(),
            geometry_provider=BrokenProvider(.15,empty)))
        e=outcome(()->tick!(g,u,p,.2))
        @test e isa RPOPlanningError && e.reason===reason
        @test (e.cause isa ErrorException)==!empty
    end
    # Inside a request, the driver reports it like any other exception there.
    e=outcome(()->run_simulation(make_rpo_configuration(planner=DirectRPOPlanner(),
        geometry_provider=BrokenProvider(-1.,false),mission_time_s=.2)))
    @test e isa RPOPlanningError && e.reason===:planner_exception && e.cause isa ErrorException
end

@testset "Hold fallback: install, settled retry, renewal and recovery" begin
    still(;kwargs...)=make_rpo_configuration(;goal_rtn_m=(3.,0.,0.),on_rejection=:install_fallback,kwargs...)
    # Refused replans: a hold, one settled retry, then renewals before expiry.
    sol=run_simulation(still(planner=StepPlanner((:replan,)),mission_time_s=3.,plan_validity_s=2.,
        planning_events=(RPOPlanningEvent(.1;reason=:replan),));return_solution=true)
    records=only(rpo_run_report(sol).planners).records
    held=filter(x->x.event===:fallback_installed,records)
    # Each hold lasts plan_validity_s; it is renewed while the next tick's preview
    # still fits (12 steps of 0.1 s): every 0.8 s after the settled retry.
    @test isapprox([x.time_s for x in held],[.1,.2,1.,1.8,2.6];atol=1e-9)
    marks=triggers(records)
    @test [x.trigger for x in marks]==[:hold_settled,:hold_renewal,:hold_renewal,:hold_renewal]
    @test isapprox([x.time_s for x in marks],[.2,1.,1.8,2.6];atol=1e-9)
    hold=last(held).product
    @test hold.r_ref_rtn_m[:,1]==hold.r_ref_rtn_m[:,2] && hold.v_ref_rtn_mps==zeros(3,2)
    @test all(x->x.event!==:failure,records) # renewed, never expired
    buffer=only(sol.prob.p.args.control_model.control_effectors).plan_buffer
    @test buffer.plan.r_ref_rtn==hold.r_ref_rtn_m && buffer.plan.diagnostics.validation.reason===:infeasible
    # A refused initial request holds at the start; the settled retry recovers.
    sol=run_simulation(still(planner=StepPlanner((:initial,)),mission_time_s=.4);return_solution=true)
    records=only(rpo_run_report(sol).planners).records
    @test lifecycle(records)==[:rejected,:fallback_installed,:trigger,:installed]
    @test last(filter(x->x.event===:installed,records)).reason===:replan
    # A new geometry revision replaces a hold at once, without a clearance check.
    late=RPOReplanningSpheres([(center_rtn=(0.,0.,9.),radius_m=.1,appear_time_s=.35)])
    sol=run_simulation(still(planner=StepPlanner((:replan,)),mission_time_s=.6,geometry_provider=late,
        planning_events=(RPOPlanningEvent(.1;reason=:replan),));return_solution=true)
    marks=triggers(only(rpo_run_report(sol).planners).records)
    @test [x.trigger for x in marks]==[:hold_settled,:geometry_revision,:hold_settled]
    @test marks[2].details.action===:replan && isnan(marks[2].details.min_clearance_m)
    # Tracking-error replans wait while the chaser brakes onto a hold.
    owned,u,p,g=prepare(still(planner=StepPlanner((:replan,)),tracking_error_limit_m=.05,
        planning_events=(RPOPlanningEvent(.1;reason=:replan),)))
    tick!(g,u,p,.1)
    @test g.triggers.hold===:holding && g.runtime.request_id==2
    move_chaser!(u,(3.5,0.,0.))
    tick!(g,u,p,.2)
    @test g.runtime.request_id==2 && g.triggers.hold===:holding
    # Interval replans wait too: with an interval of one guidance tick, a planner that
    # keeps refusing is asked once more when the chaser settles, not at every tick.
    sol=run_simulation(still(planner=StepPlanner((:replan,)),mission_time_s=1.,replan_interval_s=.1,
        planning_events=(RPOPlanningEvent(.1;reason=:replan),));return_solution=true)
    report=only(rpo_run_report(sol).planners)
    @test report.request_id==3 # the initial plan, the refused event and one settled retry
    @test [x.trigger for x in triggers(report.records)]==[:hold_settled]
    @test count(x->x.event===:fallback_installed,report.records)==2
end

@testset "Retiming: events and tracking error" begin
    # A :retime event installs a retimed reference.
    sol=run_simulation(make_rpo_configuration(planner=hypr(),mission_time_s=.5,
        planning_events=(RPOPlanningEvent(.3;reason=:retime),));return_solution=true)
    installed=filter(x->x.event===:installed,only(rpo_run_report(sol).planners).records)
    @test [x.reason for x in installed]==[:initial,:retime]
    @test installed[2].diagnostics.revision===:retime && installed[2].validation.accepted
    @test installed[2].reference.r_ref_rtn_m[:,end]==installed[1].reference.r_ref_rtn_m[:,end]
    # Here the manuscript retimer turns faster than the common acceleration limit; the
    # refused retime becomes a hold under the fallback.
    sol=run_simulation(make_rpo_configuration(planner=hypr(:manuscript),mission_time_s=.5,
        on_rejection=:install_fallback,planning_events=(RPOPlanningEvent(.3;reason=:retime),));return_solution=true)
    records=only(rpo_run_report(sol).planners).records
    @test lifecycle(records)[1:3]==[:installed,:rejected,:fallback_installed]
    rejected=only(filter(x->x.event===:rejected,records))
    @test rejected.reason===:retime && rejected.validation.reason===:failed
    # A tracking error above the threshold asks for a retime from the chaser's
    # position (here with speed and clearance limits only).
    unlimited=RPOPlanningConstraints(clearance_m=.2,max_speed_mps=.5)
    owned,u,p,g=prepare(make_rpo_configuration(planner=hypr(),constraints=unlimited,retime_tracking_error_m=.1))
    move_chaser!(u,(3.,.05,0.))
    tick!(g,u,p,.1)
    @test g.runtime.request_id==1 && isempty(triggers(g.runtime.records))
    move_chaser!(u,(3.,.15,0.))
    tick!(g,u,p,.2)
    @test only(triggers(g.runtime.records)).trigger===:retime_tracking_error
    retimed=last(g.runtime.records)
    @test retimed.event===:installed && retimed.reason===:retime && retimed.validation.accepted
    request=last(filter(x->x.event===:request,g.runtime.records)).request
    @test retimed.reference.r_ref_rtn_m[:,1]==collect(request.x_rtn[1:3])
    @test norm(collect(request.x_rtn[1:3])-[3.,.15,0.])<1e-9
end
end
