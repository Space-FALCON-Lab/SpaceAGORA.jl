using StaticArrays: SVector, MVector
using LinearAlgebra: Diagonal, I

"""
    make_rpo_configuration(; planner, seed=741, start_rtn_m=(3,0,0),
        goal_rtn_m=(5,0,0), station_points_m=zeros(3,1), control_dt_s=0.1,
        mission_time_s=2.0, plan_validity_s=120.0, extra_effectors=(), ...)

Build a bounded Earth RPO pilot using public planner contracts: circular 420 km
orbit, truth observations, zero initial rotating-RTN relative velocity, static
RTN point-cloud geometry, inverse-square gravity, no atmosphere and LQ-MPC.
The run returns the usual simulation solution; inspect it with `rpo_run_report`.

`chasers` may be a tuple of `(id, start_rtn_m, goal_rtn_m)` named tuples sharing
one passive target (`target_id=201`). Each gets a copied planner and private RNG.
This geometry checks the station only, not collisions between chasers.
`planning_events` is a sequence of `RPOPlanningEvent`; optional positive
`replan_interval_s` and `tracking_error_limit_m` enable background replanning.
A failed request stops the run. Built-in planners do not implement retiming.

Defaults: station keepout 0.25 m plus largest chaser half-extent 0.15 m;
clearance 0.2 m; speed 0.5 m/s; reference acceleration 0.5*0.05/5.2 m/s².
The controller uses a 12-step preview. `extra_effectors` extends the existing
force interface. Checkpoints and `isolate_state=false` are refused at startup.
`simulation_settings`/`solver_config` may use their existing public types. Only
the default `:tsit5` route is validated for this lifecycle; other solver modes are
outside the pilot. The first control update needs 13 steps of validity, although
the constructor accepts 12; periodic replanning also needs to cover the interval
to the next guidance tick plus the 12-step preview. Events with no eligible tick
before the mission end remain undelivered; inspect the run report.
This is a reference/tracking demonstration, not a certified rendezvous model.
"""
function make_rpo_configuration(; planner::P.AbstractRPOPlanner, seed::Integer=741,
        start_rtn_m=(3.,0.,0.), goal_rtn_m=(5.,0.,0.), chaser_id::Integer=101,
        target_id::Integer=201, chasers=nothing, station_points_m=zeros(3,1),
        station_keepout_radius_m=.25, geometry_revision="rpo-pilot-v1",
        control_dt_s=.1, mission_time_s=2., plan_validity_s=120.,
        constraints=P.RPOPlanningConstraints(clearance_m=.2,max_speed_mps=.5,
            max_acceleration_mps2=.5*.05/5.2),
        validation=P.RPOValidationSettings(), planning_events=(),
        replan_interval_s=Inf, tracking_error_limit_m=Inf, extra_effectors=(),
        simulation_settings=S.SimulationSettings(results=false,verbose=false,generate_plots=false),
        solver_config=S.SolverConfig(solver_mode=:tsit5))
    dt=P._positive(control_dt_s,"control_dt_s")
    mission=P._positive(mission_time_s,"mission_time_s")
    validity=P._positive(plan_validity_s,"plan_validity_s")
    validity >= 12dt || throw(ArgumentError("Plan lifetime must cover the controller preview."))
    for (name,value) in (("replan_interval_s",replan_interval_s),("tracking_error_limit_m",tracking_error_limit_m))
        (isfinite(value) && value>0) || value==Inf || throw(ArgumentError("$name must be positive or Inf."))
    end
    specs=chasers===nothing ? ((id=Int(chaser_id),start_rtn_m=start_rtn_m,goal_rtn_m=goal_rtn_m),) : Tuple(chasers)
    isempty(specs) && throw(ArgumentError("At least one chaser is required."))
    ids=[Int(x.id) for x in specs]
    all(>(0),ids) && target_id>0 && !(target_id in ids) && length(unique(ids))==length(ids) ||
        throw(ArgumentError("Distinct positive spacecraft IDs are required."))
    points=copy(Matrix{Float64}(station_points_m))
    all(isfinite,points) || throw(ArgumentError("Station points must be finite."))
    keepout=P._nonnegative(station_keepout_radius_m,"station_keepout_radius_m")
    geometry=S.RPOReferenceGeometry(S.RPOStationGeometry(points;keepout_radius_m=keepout);
        chaser=S.RPOCubeSatGeometry(dims_m=(.1,.1,.3)))
    revision=P._nonempty(geometry_revision,"geometry_revision")
    events=collect(RPOPlanningEvent,planning_events)
    # Identical duplicate tokens are repeat delivery; conflicting reuse is ambiguous.
    for e in events, other in events
        e.token==other.token && (e.time_s!=other.time_s || e.reason!=other.reason) &&
            throw(ArgumentError("A planning token cannot name different events."))
    end
    sort!(events;by=e->e.time_s,alg=Base.Sort.MergeSort)
    selected = planner isa DirectRPOPlanner && planner.segment_clearance_at===nothing ?
        DirectRPOPlanner(segment_clearance_at=station_segment_clearance,
            headroom=planner.headroom,max_reference_samples=planner.max_reference_samples) : planner
    planet=S.make_no_gram_planet(:earth)
    radius=planet.Rp_e+420e3; n=sqrt(planet.μ/radius^3)
    rt=SVector(radius,0.,0.);vt=SVector(0.,sqrt(planet.μ/radius),0.)
    q=SVector(0.,0.,0.,1.)
    root=S.Link(root=true,m=500.,dims=MVector(4.,2.,2.),ref_area=8.)
    target=S.SpacecraftModel(joints=S.Joint[],links=[root],root=root,prop_mass=0.,
        inertia_tensor=root.inertia,id=Int(target_id),
        initial_condition=S.CartesianInitialCondition(rt,vt;q=q,ang_vel=SVector(0.,0.,n)))
    spacecraft=S.SpacecraftModel[];guidance=PlannerGuidance[];controls=S.RPOMPCControlModel[]
    for (i,spec) in enumerate(specs)
        start=P._state_tuple(spec.start_rtn_m,Val(3),"start_rtn_m")
        goal=P._state_tuple(spec.goal_rtn_m,Val(3),"goal_rtn_m")
        r,v=S.FrameTransforms.rtn_to_inertial_relative_state(SVector{3}(start),SVector(0.,0.,0.),rt,vt)
        body=S.Link(root=true,m=5.,dims=MVector(.1,.1,.3),ref_area=.03)
        sc=S.SpacecraftModel(joints=S.Joint[],links=[body],root=body,prop_mass=.2,
            inertia_tensor=body.inertia,id=Int(spec.id),
            initial_condition=S.CartesianInitialCondition(r,v;q=q))
        push!(spacecraft,sc)
        buffer=S.RPOPlanBuffer()
        push!(guidance,PlannerGuidance(planner=deepcopy(selected),chaser_id=Int(spec.id),
            target_id=Int(target_id),chaser_idx=i,target_idx=length(specs)+1,goal_rtn_m=goal,
            geometry=deepcopy(geometry),geometry_revision=revision,constraints=constraints,
            validation=validation,reference_dt_s=dt,preview_horizon_steps=12,validity_s=validity,
            seed=Int(seed),events=copy(events),replan_interval_s=Float64(replan_interval_s),
            tracking_error_limit_m=Float64(tracking_error_limit_m),plan_buffer=buffer))
        Q=Matrix(Diagonal([20.,20.,20.,2.,2.,2.]));R=.1*Matrix{Float64}(I,3,3)
        controller=S.init_rpo_lqmpc(n,dt,Q,R,10Q,12;
            u_min=fill(-.05/5.2,3),u_max=fill(.05/5.2,3))
        thrusters=S.SixAxisThrusterModel(max_thrust_n=SVector{6}(fill(.05,6)),isp_s=SVector{6}(fill(60.,6)))
        push!(controls,S.RPOMPCControlModel(chaser_idx=i,target_idx=length(specs)+1,
            thrusters=thrusters,controller=controller,plan_buffer=buffer,control_dt_s=dt,
            attitude_kp=.02,rate_kd=.08,max_rw_torque_nm=.002,command_log=S.RPOControlCommandLog()))
    end
    push!(spacecraft,target)
    return S.SimulationConfiguration(simulation_settings=simulation_settings,
        mission_configuration=S.MissionConfiguration(mission_type=S.MissionTime,keplerian=false,
            number_of_orbits=1,mission_time=mission,orientation_sim=true,num_steps_to_save=400,data_rate=dt),
        environment_model=S.EnvironmentModel(planet=planet,EI=120.,density_model=S.NoAtmosphereModel(),
            ephemerides_model=S.SimpleEphemeridesModel(),
            thermal_model=S.MaxwellianHeat(thermal_accomodation_factor=1.,planet=planet),topography=false,wind=false),
        dynamics_model=S.DynamicsModel(spacecraft,(S.InverseSquaredGravityModel(),extra_effectors...)),
        guidance_model=S.GuidanceModel(guidance_effectors=Tuple(guidance),guidance_rates=fill(dt,length(specs))),
        navigation_model=S.NavigationModel(navigation_effectors=(),navigation_rates=Float64[]),
        control_model=S.ControlModel(control_effectors=Tuple(controls),control_rates=fill(dt,length(specs))),
        initial_time=S.InitialTime(year=2026,month=1,day=1),
        integration_tolerances=S.IntegrationTolerances(reltol_orbit=1e-8,abstol_orbit=1e-8,
            reltol_quaternion=1e-8,abstol_quaternion=1e-8,dt_max_orbit=.05),solver_config=solver_config)
end
