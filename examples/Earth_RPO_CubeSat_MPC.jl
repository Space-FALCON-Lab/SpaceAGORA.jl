# Load from the companion environment before common.jl continues in the examples environment.
# Setup: julia --startup-file=no examples/rpo_planner_env/setup.jl
#        julia --startup-file=no examples/setup.jl
using SpaceAGORA, HYPR
include(joinpath(@__DIR__, "common.jl"))
using CSV
using DataFrames
using PlotlyJS
using Random
using StaticArrays
using LinearAlgebra
using DiffEqBase: DiscreteCallback, terminate!

const DEFAULT_STATION_DIMS_M = (4.0, 2.0, 2.0)
const DEFAULT_STATION_MASS_KG = 500.0
const DEFAULT_STATION_REF_AREA_M2 = 8.0
const DEFAULT_RPO_EPOCH = InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0)

"""
    _rpo_station_settings(; points, n_points, seed, keepout_radius_m, name, dims_m, mass_kg, ref_area_m2)

Resolve the station for `build_rpo_cubesat_mpc_demo`. With `points === nothing`
the Gateway core CAD is sampled as before (`source = :gateway`); otherwise any
3 x N body-frame point cloud in meters stands in (`source = :custom`, and the
sampling settings are unused). `ref_area_m2 === nothing` keeps the reviewed
8 m² of the default box and takes the y-z face area `dims_m[2] * dims_m[3]`
for other dimensions. The scenario has no atmosphere, so the area only
matters once a caller adds aerodynamic models.
"""
function _rpo_station_settings(; points, n_points, seed, keepout_radius_m, name, dims_m, mass_kg, ref_area_m2)
    length(dims_m) == 3 || throw(ArgumentError("station_dims_m must hold three lengths in meters, got $(dims_m)."))
    dims = (Float64(dims_m[1]), Float64(dims_m[2]), Float64(dims_m[3]))
    all(d -> isfinite(d) && d > 0.0, dims) || throw(ArgumentError("station_dims_m must be positive and finite, got $(dims_m)."))
    mass = Float64(mass_kg)
    (isfinite(mass) && mass > 0.0) || throw(ArgumentError("station_mass_kg must be positive and finite, got $(mass_kg)."))
    area = if ref_area_m2 === nothing
        dims == DEFAULT_STATION_DIMS_M ? DEFAULT_STATION_REF_AREA_M2 : dims[2] * dims[3]
    else
        Float64(ref_area_m2)
    end
    (isfinite(area) && area > 0.0) || throw(ArgumentError("station_ref_area_m2 must be positive and finite, got $(ref_area_m2)."))
    if points === nothing
        source = :gateway
        cloud = SpaceAGORA.load_rpo_station_cad_pointcloud(:gateway; n_points=n_points, rng=MersenneTwister(seed))
    else
        source = :custom
        cloud = Matrix{Float64}(points)
        size(cloud, 1) == 3 || throw(ArgumentError("station_points must be a 3 x N body-frame point cloud in meters, got size $(size(cloud))."))
    end
    return (
        name=String(name),
        source=source,
        points=cloud,
        n_points=size(cloud, 2),
        keepout_radius_m=Float64(keepout_radius_m),
        dims_m=dims,
        mass_kg=mass,
        ref_area_m2=area,
    )
end

"""
    rpo_reference_stop_callback(min_time_s, hold_s)

Stop the run `hold_s` after the reference the planner installed ends, but not before
`min_time_s`. The planner plans inside the run, so the run's length cannot be set
from the plan beforehand. Pass the callback to `run_simulation` through
`extra_callbacks`; without it the run lasts `max_mission_time_s`.
"""
function rpo_reference_stop_callback(min_time_s::Float64, hold_s::Float64)
    function reference_done(u, t, integrator)
        for c in integrator.p.args.control_model.control_effectors
            c isa SimulationModel.RPOMPCControlModel || continue
            plan, buffer = c.plan_buffer.plan, c.plan_buffer
            buffer.valid && !isempty(plan.t_ref_s) || return false
            return t >= max(min_time_s, buffer.updated_at_s + plan.t_ref_s[end] + hold_s)
        end
        return false
    end
    return DiscreteCallback(reference_done, terminate!; save_positions=(false, false))
end

"""
    RPOMissionCapPlanner(planner, end_by_s)

Plan with `planner`, but refuse a candidate whose reference ends after `end_by_s`
(`:infeasible`, `:reference_ends_after_mission_cap`). The request's lifetime also
covers the controller's preview past the end of the run, so `planner` alone would
accept a reference that ends up to one preview after it, and the run would end
before the reference does.
"""
struct RPOMissionCapPlanner{T<:AbstractRPOPlanner} <: AbstractRPOPlanner
    planner::T
    end_by_s::Float64
end
SpaceAGORA.planner_capabilities(p::RPOMissionCapPlanner) = planner_capabilities(p.planner)
SpaceAGORA.consumes_rng(p::RPOMissionCapPlanner) = consumes_rng(p.planner)
SpaceAGORA.RPOPlannerInterfaces.require_planner_support(p::RPOMissionCapPlanner) =
    SpaceAGORA.RPOPlannerInterfaces.require_planner_support(p.planner)
SpaceAGORA.initialize_planner(p::RPOMissionCapPlanner, run_context) = initialize_planner(p.planner, run_context)
SpaceAGORA.plan_rpo!(state, p::RPOMissionCapPlanner, request::RPOPlanningRequest, rng::AbstractRNG) =
    _rpo_within_mission_cap(plan_rpo!(state, p.planner, request, rng), request, p.end_by_s)
SpaceAGORA.retime_rpo!(state, p::RPOMissionCapPlanner, request::RPOPlanningRequest, active::RPOReference, rng::AbstractRNG) =
    _rpo_within_mission_cap(retime_rpo!(state, p.planner, request, active, rng), request, p.end_by_s)

function _rpo_within_mission_cap(result, request, end_by_s)
    result.status === :candidate || return result
    ends_s = request.time_s + last(result.reference.t_ref_s)
    ends_s <= end_by_s && return result
    return RPOPlanningResult(request_id=result.request_id, status=:infeasible,
        termination=:reference_ends_after_mission_cap,
        diagnostics=merge(result.diagnostics, (reference_end_s=ends_s, max_mission_time_s=end_by_s)))
end

"""
    rpo_installed_plan(solution)

The reference the planner installed in a run, read from `rpo_run_report`, with the
planner's own record of it: `t_ref_s`, `r_ref_rtn` and `v_ref_rtn` (the validated
reference), `path_rtn` and `config` (the geometric path and the configuration it was
timed with), `cost` and `cost_history` (`NaN` and empty for a replayed path),
`validation`, `planner` and
`time_s` (when it was installed).
"""
function rpo_installed_plan(solution)
    record = last(filter(x -> x.event === :installed, only(rpo_run_report(solution).planners).records))
    ref, d = record.reference, record.diagnostics
    return (t_ref_s=ref.t_ref_s, r_ref_rtn=ref.r_ref_rtn_m, v_ref_rtn=ref.v_ref_rtn_mps, path_rtn=d.path,
        config=get(d, :effective_config, d.input_config), cost=get(d, :cost, NaN),
        cost_history=collect(get(d, :cost_history, Float64[])),
        validation=record.validation, planner=d.planner, time_s=record.time_s)
end

"""
    run_rpo_case(demo; kwargs...) -> solution

Run a built scenario with its stop callback and return the solution; `kwargs` go to
`run_simulation`. In the example smoke mode the run is shortened as for other examples.
"""
function run_rpo_case(demo; kwargs...)
    args = SpaceAGORA.TelemetryVerification._example_smoke_args(demo.args)
    return run_simulation(args; return_solution=true, extra_callbacks=(demo.stop_callback,), kwargs...)
end

"""
    build_rpo_cubesat_mpc_demo(; mission_time=180.0, kwargs...)

Build a complete two-spacecraft RPO scenario: a passive station target and a
six-axis CubeSat chaser using HYPR/PSO guidance and LQ-MPC tracking.

Guidance is the planner contract: a `PlannerGuidance` with `HYPRRPOPlanner(pso_config)`
plans inside the run, at its first guidance tick, from the run's own random stream
(`seed` and the spacecraft IDs 101 and 201). The run validates the reference against
the request's clearance (`safe_distance_m`) and the retimer's limits before the
controller uses it, and refuses `isolate_state=false`. Run the
scenario with `run_rpo_case(demo)` (or pass `demo.stop_callback` to `run_simulation`)
and read the installed reference with `rpo_installed_plan(solution)`.

The station defaults to the Gateway core: a sampled CAD point cloud with a
0.25 m keep-out radius on a 500 kg, 4 x 2 x 2 m box. Pass `station_points`
(3 x N body-frame meters, for example from `sample_model_pointcloud`) with
`station_keepout_radius_m`, `station_name`, `station_dims_m`,
`station_mass_kg` and optionally `station_ref_area_m2` for another station,
and scale the planner with `safe_distance_m`, `cost_ref_distance_m`,
`search_margin_m` and `sample_ds_m`. `reference_max_speed_mps` optionally
caps the reference speed without changing the geometric plan. Defaults keep
the Gateway geometry, mass and reference area; the returned `station` record
states what was used.

`precomputed_plan` replays a retained, endpoint-matched geometric plan without
running another optimizer: `RetainedPathRPOPlanner(plan.path, plan.config)` times it
for the run's request, and the run validates it. `chaser_dry_mass_kg` defaults to 5.0 kg, plus the
existing 0.2 kg propellant; pass 4.8 kg to reconstruct a 5 kg wet-mass experiment.
`minimum_throttle_fraction` defaults to zero. A value of 0.2 suppresses commands
at or below 20% of the rated thrust, then holds the remaining average thrust.
This does not model individual valve pulses or establish the paper's exact PWM.

The planner always receives the orbit's mean motion (`mean_motion_radps`),
which the `:manuscript` HyPR mode uses in its HCW fuel proxy.
`retime_accel_mps2` sets the retiming acceleration (default half the per-axis
thrust acceleration), `mpc_terminal_weight_scale` sets Qf = scale * Q
(default 10), `mpc_horizon` is the LQ-MPC horizon in control steps (default
12). The stop callback ends the run `post_reference_hold_s` after the installed
reference ends (default 20 s), but not before `mission_time`; `max_mission_time_s`
(default 1800 s) is the longest run. The guidance planner is `planner` wrapped in
`RPOMissionCapPlanner`: a reference that ends after the cap is refused
(`:reference_ends_after_mission_cap`), and the run stops with an `RPOPlanningError`.
The plan's lifetime is the cap plus the controller's preview. `record_control_commands`
attaches an `RPOControlCommandLog` to the controller; read the run's copy from
`rpo_run_report(solution).commands`. `initial_time` is the epoch
(default 2026-01-01 00:00 UTC). The station starts at a fixed inertial point,
so the date sets the lighting; the relative motion does not depend on it.
"""
function build_rpo_cubesat_mpc_demo(;
    mission_time=180.0,
    results_directory=joinpath(REPO_ROOT, "output", "rpo_single_case"),
    seed::Integer=741,
    start_rtn=SVector{3, Float64}(-8.0, -4.0, 2.0),
    goal_rtn=SVector{3, Float64}(5.0, 0.0, 0.0),
    pso_n_particles::Integer=120,
    pso_n_iters::Integer=35,
    pso_config=nothing,
    precomputed_plan=nothing,
    pso_configurator=nothing,
    n_station_points::Integer=10000,
    station_geometry_seed::Integer=seed,
    station_points=nothing,
    station_keepout_radius_m::Real=0.25,
    station_name::AbstractString="gateway_core",
    station_dims_m=DEFAULT_STATION_DIMS_M,
    station_mass_kg::Real=DEFAULT_STATION_MASS_KG,
    station_ref_area_m2=nothing,
    safe_distance_m::Real=0.1,
    cost_ref_distance_m::Real=20.0,
    search_margin_m=nothing,
    sample_ds_m::Real=0.05,
    reference_max_speed_mps=nothing,
    data_rate_s::Real=10.0,
    pso_iteration_runtime_limit_s=nothing,
    pso_iteration_callback=nothing,
    retime_accel_mps2=nothing,
    chaser_dry_mass_kg::Real=5.0,
    minimum_throttle_fraction::Real=0.0,
    mpc_terminal_weight_scale::Real=10.0,
    mpc_horizon::Integer=12,
    post_reference_hold_s::Real=20.0,
    max_mission_time_s::Real=1800.0,
    record_control_commands::Bool=false,
    initial_time::InitialTime=DEFAULT_RPO_EPOCH,
    verbose::Bool=true,
)
    mpc_horizon >= 1 || throw(ArgumentError("mpc_horizon must be at least 1."))
    start_rtn = SVector{3, Float64}(start_rtn)
    goal_rtn = SVector{3, Float64}(goal_rtn)
    planet = Earth("", SPICE_PATH)
    orbit_radius = planet.Rp_e + 420e3
    μ = planet.μ
    n = sqrt(μ / orbit_radius^3)
    r_station_ii = SVector{3, Float64}(orbit_radius, 0.0, 0.0)
    v_station_ii = SVector{3, Float64}(0.0, sqrt(μ / orbit_radius), 0.0)
    r0_rel_rtn = start_rtn
    v0_rel_rtn = SVector{3, Float64}(0.0, 0.0, 0.0)
    r_chaser_ii, v_chaser_ii = SimulationModel.GuidanceHooks.rtn_to_inertial_relative_state(
        r0_rel_rtn,
        v0_rel_rtn,
        r_station_ii,
        v_station_ii,
    )

    station_spec = _rpo_station_settings(;
        points=station_points,
        n_points=n_station_points,
        seed=station_geometry_seed,
        keepout_radius_m=station_keepout_radius_m,
        name=station_name,
        dims_m=station_dims_m,
        mass_kg=station_mass_kg,
        ref_area_m2=station_ref_area_m2,
    )
    station_geometry = RPOStationGeometry(station_spec.points; keepout_radius_m=station_spec.keepout_radius_m, name=station_spec.name)
    chaser_geometry = RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.3))
    geometry = RPOReferenceGeometry(station_geometry; chaser=chaser_geometry)

    q_identity = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
    station_root = Link(
        root=true,
        m=station_spec.mass_kg,
        dims=MVector{3, Float64}(station_spec.dims_m...),
        ref_area=station_spec.ref_area_m2,
    )
    station = SpacecraftModel(
        joints=Joint[],
        links=[station_root],
        root=station_root,
        prop_mass=0.0,
        inertia_tensor=station_root.inertia,
        # The planner's body-frame station cloud is fixed in RTN. For this
        # circular, equatorial, torque-free orbit, spin about the normal at n
        # keeps the simulated body and its displayed model in that frame.
        initial_condition=CartesianInitialCondition(r_station_ii, v_station_ii;
            q=q_identity, ang_vel=SVector{3, Float64}(0.0, 0.0, n)),
        id=201,
    )

    chaser_root = Link(
        root=true,
        m=Float64(chaser_dry_mass_kg),
        dims=MVector{3, Float64}(0.1, 0.1, 0.3),
        ref_area=0.03,
    )
    chaser = SpacecraftModel(
        joints=Joint[],
        links=[chaser_root],
        root=chaser_root,
        prop_mass=0.2,
        inertia_tensor=chaser_root.inertia,
        initial_condition=CartesianInitialCondition(r_chaser_ii, v_chaser_ii; q=q_identity),
        id=101,
    )

    thrusters = SixAxisThrusterModel(
        max_thrust_n=SVector{6, Float64}(fill(0.05, 6)),
        isp_s=SVector{6, Float64}(fill(60.0, 6)),
    )
    chaser_initial_mass_kg = chaser_root.m + chaser.prop_mass
    max_axis_accel_mps2 = minimum(thrusters.max_thrust_n ./ chaser_initial_mass_kg)
    retime_accel_mps2 = retime_accel_mps2 === nothing ? 0.5 * max_axis_accel_mps2 : Float64(retime_accel_mps2)
    isfinite(retime_accel_mps2) && retime_accel_mps2 > 0.0 ||
        throw(ArgumentError("retime_accel_mps2 must be finite and positive."))

    if pso_config !== nothing && pso_configurator !== nothing
        throw(ArgumentError("Pass either pso_config or pso_configurator, not both."))
    end
    default_pso_configurator = RPOPSOConfigurator(
        swarm=RPOPSOSwarmSettings(
            n_waypoints=5,
            n_particles=Int(pso_n_particles),
            n_iters=Int(pso_n_iters),
            curve_type=:bezier,
            sample_ds_m=Float64(sample_ds_m),
        ),
        adaptive=RPOPSOAdaptiveSettings(
            allow_downscale=true,
            n_particles_min=80,
            n_particles_max=160,
            n_iters_min=8,
            n_iters_max=35,
        ),
        early_stopping=RPOPSOEarlyStoppingSettings(
            enabled=true,
            patience=8,
            min_iters=12,
            min_rel_improvement=1.0e-4,
        ),
        objective=RPOPSOObjectiveSettings(
            cost_ref_distance_m=Float64(cost_ref_distance_m),
            mass_kg=chaser_initial_mass_kg,
            tf_s=Float64(mission_time),
        ),
        retiming=RPOPSORetimingSettings(
            dt_s=0.1,
            a_max_mps2=retime_accel_mps2,
        ),
    )
    pso_cfg = if pso_config !== nothing
        rpo_pso_config(pso_config)
    elseif pso_configurator !== nothing
        rpo_pso_config(
            pso_configurator;
            mass_kg=chaser_initial_mass_kg,
            tf_s=Float64(mission_time),
            retime_a_max_mps2=retime_accel_mps2,
        )
    else
        rpo_pso_config(default_pso_configurator)
    end
    if pso_iteration_runtime_limit_s !== nothing
        pso_cfg = rpo_pso_config(pso_cfg; iteration_runtime_limit_s=Float64(pso_iteration_runtime_limit_s))
    end
    if search_margin_m !== nothing
        pso_cfg = rpo_pso_config(pso_cfg; search_margin_m=Float64(search_margin_m))
    end
    if reference_max_speed_mps !== nothing
        speed_limit = Float64(reference_max_speed_mps)
        isfinite(speed_limit) && speed_limit > 0.0 ||
            throw(ArgumentError("reference_max_speed_mps must be finite and positive."))
        pso_cfg = rpo_pso_config(pso_cfg; retime_max_speed_mps=speed_limit)
    end
    pso_cfg = rpo_pso_config(pso_cfg; mean_motion_radps=n)

    pso_iteration_callback === nothing || throw(ArgumentError(
        "pso_iteration_callback is not supported: the planner now plans inside the run through HYPRRPOPlanner."))
    # The planner plans inside the run, at its first guidance tick, from the run's own
    # random stream (`seed` and the two spacecraft IDs). The run validates the reference
    # before the controller uses it; read it afterwards with `rpo_installed_plan`.
    if precomputed_plan === nothing
        planner = HYPRRPOPlanner(pso_cfg)
        planning_cfg = pso_cfg
    else
        size(precomputed_plan.path,1) == 3 && size(precomputed_plan.path,2) >= 2 ||
            throw(ArgumentError("Precomputed RPO path must contain at least two 3D points."))
        norm(precomputed_plan.path[:,1]-r0_rel_rtn) <= 1e-8 && norm(precomputed_plan.path[:,end]-goal_rtn) <= 1e-8 ||
            throw(ArgumentError("Precomputed RPO path endpoints must match the scenario."))
        # A saved plan is replayed, timed for the run's request, and validated like any other.
        planner = RetainedPathRPOPlanner(precomputed_plan.path, precomputed_plan.config)
        planning_cfg = rpo_pso_config(precomputed_plan.config)
    end
    # Limits the run validates: the clearance always; `reference_max_speed_mps` when it
    # is given (the planner then keeps 1% headroom below it); the thrusters' per-axis
    # acceleration when the retimer limits acceleration (it plans at half of it by
    # default). A retimer setting alone is a planner choice, not a validated limit.
    constraints = RPOPlanningConstraints(
        clearance_m=Float64(safe_distance_m),
        max_speed_mps=reference_max_speed_mps === nothing ? nothing : Float64(reference_max_speed_mps),
        max_acceleration_mps2=planning_cfg.retime_accel_limit_enable ? max_axis_accel_mps2 : nothing,
    )
    # A wall-clock-limited PSO result is accepted only when the caller set that limit.
    validation = RPOValidationSettings(allow_time_budget_candidate=pso_iteration_runtime_limit_s !== nothing)
    max_mission_time_s = Float64(max_mission_time_s)
    isfinite(max_mission_time_s) && max_mission_time_s >= Float64(mission_time) ||
        throw(ArgumentError("max_mission_time_s must be finite and at least mission_time."))
    plan_buffer = RPOPlanBuffer()
    guidance = SpaceAGORA.RPOPlannerLifecycle.PlannerGuidance(
        # The reference itself must end by the cap; the lifetime below adds the preview.
        planner=RPOMissionCapPlanner(planner, max_mission_time_s),
        chaser_id=101,
        target_id=201,
        chaser_index=1,
        target_index=2,
        goal_rtn_m=Tuple(goal_rtn),
        geometry=geometry,
        geometry_revision=station_spec.name,
        constraints=constraints,
        validation=validation,
        reference_dt_s=pso_cfg.retime_dt_s,
        preview_horizon_steps=Int(mpc_horizon),
        # The controller previews the reference up to the end of the longest run allowed.
        validity_s=max_mission_time_s + (Int(mpc_horizon) + 1) * pso_cfg.retime_dt_s,
        seed=Int(seed),
        events=RPOPlanningEvent[],
        replan_interval_s=Inf,
        tracking_error_limit_m=Inf,
        plan_buffer=plan_buffer,
    )

    Q = Diagonal([20.0, 20.0, 20.0, 2.0, 2.0, 2.0])
    R = 0.1I(3)
    Qf = Float64(mpc_terminal_weight_scale) .* Matrix(Q)
    controller = init_rpo_lqmpc(
        n,
        pso_cfg.retime_dt_s,
        Matrix(Q),
        Matrix(R),
        Qf,
        Int(mpc_horizon);
        u_min=fill(-max_axis_accel_mps2, 3),
        u_max=fill(max_axis_accel_mps2, 3),
    )
    control = RPOMPCControlModel(
        chaser_index=1,
        target_index=2,
        thrusters=thrusters,
        minimum_throttle_fraction=Float64(minimum_throttle_fraction),
        controller=controller,
        plan_buffer=plan_buffer,
        control_dt_s=pso_cfg.retime_dt_s,
        attitude_kp=0.02,
        rate_kd=0.08,
        max_rw_torque_nm=0.002,
        command_log=record_control_commands ? RPOControlCommandLog() : nothing,
    )

    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=true,
            verbose=verbose,
            generate_plots=false,
            results_directory=results_directory,
            normalize=false,
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=false,
            number_of_orbits=1,
            mission_time=max_mission_time_s,
            orientation_sim=true,
            num_steps_to_save=400,
            data_rate=Float64(data_rate_s),
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false,
        ),
        dynamics_model=DynamicsModel([chaser, station], (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[pso_cfg.retime_dt_s]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(control,), control_rates=[pso_cfg.retime_dt_s]),
        initial_time=initial_time,
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_quaternion=1e-8,
            abstol_quaternion=1e-8,
            dt_max_orbit=0.05,
        ),
    )

    return (
        args=args,
        geometry=geometry,
        guidance=guidance,
        control=control,
        pso_config=pso_cfg,
        plan_buffer=plan_buffer,
        initial_relative_state_rtn=vcat(r0_rel_rtn, v0_rel_rtn),
        goal_rtn=goal_rtn,
        planner=planner,
        stop_callback=rpo_reference_stop_callback(Float64(mission_time), Float64(post_reference_hold_s)),
        mission=(min_time_s=Float64(mission_time), max_time_s=max_mission_time_s,
            post_reference_hold_s=Float64(post_reference_hold_s)),
        seed=Int(seed),
        orbit=(
            radius_m=orbit_radius,
            altitude_m=orbit_radius - planet.Rp_e,
            mu_m3ps2=μ,
            mean_motion_radps=n,
        ),
        chaser=(
            dry_mass_kg=chaser_root.m,
            prop_mass_kg=chaser.prop_mass,
            initial_mass_kg=chaser_initial_mass_kg,
            max_thrust_n=collect(thrusters.max_thrust_n),
            isp_s=collect(thrusters.isp_s),
            max_axis_accel_mps2=max_axis_accel_mps2,
            dims_m=(0.1, 0.1, 0.3),
        ),
        mpc=(
            Q_diag=collect(diag(Q)),
            R_diag=collect(diag(Matrix(R))),
            Qf_scale=Float64(mpc_terminal_weight_scale),
            horizon=Int(mpc_horizon),
            dt_s=pso_cfg.retime_dt_s,
            u_limit_mps2=max_axis_accel_mps2,
        ),
        station=(
            name=station_spec.name,
            source=station_spec.source,
            n_points=station_spec.n_points,
            keepout_radius_m=station_spec.keepout_radius_m,
            dims_m=station_spec.dims_m,
            mass_kg=station_spec.mass_kg,
            ref_area_m2=station_spec.ref_area_m2,
            safe_distance_m=Float64(safe_distance_m),
            cost_ref_distance_m=pso_cfg.cost_ref_distance_m,
            search_margin_m=pso_cfg.search_margin_m,
            sample_ds_m=pso_cfg.sample_ds_m,
        ),
    )
end

function _interp_ref_columns(values::Matrix{Float64}, t_ref::Vector{Float64}, t::Real)
    size(values, 2) == 0 && return zeros(3)
    tq = Float64(t)
    tq <= t_ref[1] && return copy(values[:, 1])
    tq >= t_ref[end] && return copy(values[:, end])
    hi = searchsortedfirst(t_ref, tq)
    lo = max(hi - 1, 1)
    denom = max(t_ref[hi] - t_ref[lo], 1.0e-9)
    α = (tq - t_ref[lo]) / denom
    return (1.0 - α) .* values[:, lo] .+ α .* values[:, hi]
end

function _earth_surface_trace(radius_m; n_lon=48, n_lat=24)
    lon = range(0.0, 2π; length=n_lon)
    lat = range(-π / 2, π / 2; length=n_lat)
    x = [radius_m * cos(φ) * cos(λ) for φ in lat, λ in lon]
    y = [radius_m * cos(φ) * sin(λ) for φ in lat, λ in lon]
    z = [radius_m * sin(φ) for φ in lat, λ in lon]
    return surface(
        x=x ./ 1000.0,
        y=y ./ 1000.0,
        z=z ./ 1000.0,
        opacity=0.35,
        colorscale=[[0.0, "rgb(40,90,170)"], [1.0, "rgb(70,130,190)"]],
        showscale=false,
        name="Earth",
    )
end

function _rpo_postprocess(csv_path::AbstractString, plan)
    df = CSV.read(csv_path, DataFrame)
    n_samples = nrow(df)
    actual_rtn = zeros(3, n_samples)
    ref_rtn = zeros(3, n_samples)
    tracking_error = zeros(n_samples)
    for k in 1:n_samples
        r_chaser = SVector{3, Float64}(df.sc1_pos_1[k], df.sc1_pos_2[k], df.sc1_pos_3[k])
        v_chaser = SVector{3, Float64}(df.sc1_vel_1[k], df.sc1_vel_2[k], df.sc1_vel_3[k])
        r_target = SVector{3, Float64}(df.sc2_pos_1[k], df.sc2_pos_2[k], df.sc2_pos_3[k])
        v_target = SVector{3, Float64}(df.sc2_vel_1[k], df.sc2_vel_2[k], df.sc2_vel_3[k])
        x_rel = SimulationModel.GuidanceHooks.inertial_to_rtn_relative_state(
            r_chaser,
            v_chaser,
            r_target,
            v_target,
        )
        actual_rtn[:, k] .= x_rel[1:3]
        ref_rtn[:, k] .= _interp_ref_columns(plan.r_ref_rtn, plan.t_ref_s, df.time[k] - plan.time_s)
        tracking_error[k] = norm(actual_rtn[:, k] - ref_rtn[:, k])
    end
    return df, actual_rtn, ref_rtn, tracking_error
end

function _save_plot(plot_obj, output_dir::AbstractString, filename::AbstractString)
    mkpath(output_dir)
    path = joinpath(output_dir, filename)
    PlotlyJS.savefig(plot_obj, path)
    return path
end

const _STATION_MESH_TRACE_CACHE = Ref{Any}(nothing)

function _build_station_mesh_trace()
    triangles = SpaceAGORA.load_rpo_station_cad_triangles(:gateway)
    ntri = size(triangles, 2) ÷ 3
    tri_idxs = collect(1:ntri)
    nfaces = length(tri_idxs)
    vertices = zeros(Float64, 3, 3 * nfaces)
    i = zeros(Int, nfaces)
    j = zeros(Int, nfaces)
    k = zeros(Int, nfaces)
    for (face_idx, tri_idx) in enumerate(tri_idxs)
        src = 3 * (tri_idx - 1)
        dst = 3 * (face_idx - 1)
        vertices[:, dst + 1] .= triangles[:, src + 1]
        vertices[:, dst + 2] .= triangles[:, src + 2]
        vertices[:, dst + 3] .= triangles[:, src + 3]
        i[face_idx] = dst
        j[face_idx] = dst + 1
        k[face_idx] = dst + 2
    end
    return mesh3d(
        x=vertices[1, :],
        y=vertices[2, :],
        z=vertices[3, :],
        i=i,
        j=j,
        k=k,
        color="rgb(150,160,170)",
        opacity=0.55,
        flatshading=true,
        name="Gateway mesh",
        hoverinfo="skip",
    )
end

function _station_mesh_trace(; refresh::Bool=false)
    if refresh || _STATION_MESH_TRACE_CACHE[] === nothing
        _STATION_MESH_TRACE_CACHE[] = _build_station_mesh_trace()
    end
    return deepcopy(_STATION_MESH_TRACE_CACHE[])
end

"""The Gateway mesh for the default station, otherwise a thinned scatter of the custom point cloud."""
function _station_geometry_trace(demo)
    demo.station.source === :gateway && return _station_mesh_trace()
    points = demo.geometry.station.points_body
    n = size(points, 2)
    idxs = n > 4000 ? unique(round.(Int, range(1, n; length=4000))) : collect(1:n)
    return scatter3d(
        x=points[1, idxs],
        y=points[2, idxs],
        z=points[3, idxs],
        mode="markers",
        marker=attr(size=1.5, color="rgb(150,160,170)"),
        name="$(demo.station.name) point cloud",
        hoverinfo="skip",
    )
end

function _cuboid_mesh_trace(centers::Matrix{Float64}, half_extents; max_cubes::Integer=12)
    n = size(centers, 2)
    n > 0 || return mesh3d(x=Float64[], y=Float64[], z=Float64[], i=Int[], j=Int[], k=Int[], name="CubeSat geometry")
    idxs = unique(round.(Int, range(1, n; length=min(max_cubes, n))))
    hx, hy, hz = Tuple(Float64.(half_extents))
    offsets = [
        (-hx, -hy, -hz), ( hx, -hy, -hz), ( hx,  hy, -hz), (-hx,  hy, -hz),
        (-hx, -hy,  hz), ( hx, -hy,  hz), ( hx,  hy,  hz), (-hx,  hy,  hz),
    ]
    faces = [
        (1, 2, 3), (1, 3, 4),
        (5, 7, 6), (5, 8, 7),
        (1, 5, 6), (1, 6, 2),
        (2, 6, 7), (2, 7, 3),
        (3, 7, 8), (3, 8, 4),
        (4, 8, 5), (4, 5, 1),
    ]
    vertices = zeros(Float64, 3, 8 * length(idxs))
    face_i = Int[]
    face_j = Int[]
    face_k = Int[]
    for (cube_idx, sample_idx) in enumerate(idxs)
        center = centers[:, sample_idx]
        vbase = 8 * (cube_idx - 1)
        for (offset_idx, offset) in enumerate(offsets)
            vertices[:, vbase + offset_idx] .= center .+ collect(offset)
        end
        for face in faces
            push!(face_i, vbase + face[1] - 1)
            push!(face_j, vbase + face[2] - 1)
            push!(face_k, vbase + face[3] - 1)
        end
    end
    return mesh3d(
        x=vertices[1, :],
        y=vertices[2, :],
        z=vertices[3, :],
        i=face_i,
        j=face_j,
        k=face_k,
        color="rgb(38,120,180)",
        opacity=0.38,
        flatshading=true,
        name="CubeSat geometry samples",
        hoverinfo="skip",
    )
end

function _save_rpo_single_case_plots(csv_path::AbstractString, demo, plan)
    df, actual_rtn, ref_rtn, tracking_error = _rpo_postprocess(csv_path, plan)
    plot_dir = joinpath(dirname(csv_path), "plots")
    mkpath(plot_dir)

    earth_plot = Plot(
        [
            _earth_surface_trace(demo.args.environment_model.planet.Rp_e),
            scatter3d(
                x=df.sc1_pos_1 ./ 1000.0,
                y=df.sc1_pos_2 ./ 1000.0,
                z=df.sc1_pos_3 ./ 1000.0,
                mode="lines+markers",
                name="CubeSat",
            ),
            scatter3d(
                x=df.sc2_pos_1 ./ 1000.0,
                y=df.sc2_pos_2 ./ 1000.0,
                z=df.sc2_pos_3 ./ 1000.0,
                mode="lines+markers",
                name="station",
            ),
        ],
        Layout(
            title="RPO Earth-Centered Trajectory",
            scene=attr(aspectmode="data", xaxis_title="x (km)", yaxis_title="y (km)", zaxis_title="z (km)"),
        ),
    )

    pso_path = SM.GuidanceHooks.rpo_sample_path(
        plan.path_rtn,
        plan.config.sample_ds_m;
        curve_type=plan.config.curve_type,
    )
    station_plot = Plot(
        [
            _station_geometry_trace(demo),
            scatter3d(
                x=pso_path[1, :],
                y=pso_path[2, :],
                z=pso_path[3, :],
                mode="lines",
                line=attr(width=6, color="rgb(210,85,45)"),
                name="PSO optimizer path",
            ),
            scatter3d(
                x=plan.r_ref_rtn[1, :],
                y=plan.r_ref_rtn[2, :],
                z=plan.r_ref_rtn[3, :],
                mode="lines",
                line=attr(dash="dot"),
                name="retimed reference",
            ),
            _cuboid_mesh_trace(actual_rtn, demo.geometry.chaser.half_extents_body),
            scatter3d(
                x=[demo.initial_relative_state_rtn[1]],
                y=[demo.initial_relative_state_rtn[2]],
                z=[demo.initial_relative_state_rtn[3]],
                mode="markers",
                marker=attr(size=7, color="rgb(45,150,90)", symbol="circle"),
                name="start",
            ),
            scatter3d(
                x=[demo.goal_rtn[1]],
                y=[demo.goal_rtn[2]],
                z=[demo.goal_rtn[3]],
                mode="markers",
                marker=attr(size=8, color="rgb(190,60,55)", symbol="diamond"),
                name="goal",
            ),
            scatter3d(
                x=actual_rtn[1, :],
                y=actual_rtn[2, :],
                z=actual_rtn[3, :],
                mode="lines",
                line=attr(width=4, color="rgb(20,90,150)"),
                name="actual CubeSat center",
            ),
        ],
        Layout(
            title="Station-Centered RPO Path and Tracking",
            scene=attr(aspectmode="data", xaxis_title="radial (m)", yaxis_title="along-track (m)", zaxis_title="cross-track (m)"),
        ),
    )

    tracking_plot = Plot(
        [
            scatter(x=df.time, y=tracking_error, mode="lines+markers", name="position tracking error"),
            scatter(x=df.time, y=actual_rtn[1, :] .- ref_rtn[1, :], mode="lines", name="radial error"),
            scatter(x=df.time, y=actual_rtn[2, :] .- ref_rtn[2, :], mode="lines", name="along-track error"),
            scatter(x=df.time, y=actual_rtn[3, :] .- ref_rtn[3, :], mode="lines", name="cross-track error"),
        ],
        Layout(title="RPO Tracking Error", xaxis_title="time (s)", yaxis_title="error (m)"),
    )

    cost_history = plan.cost_history
    pso_cost_plot = Plot(
        [scatter(x=1:length(cost_history), y=cost_history, mode="lines+markers", name="best PSO cost")],
        Layout(title="PSO Cost History", xaxis_title="iteration", yaxis_title="normalized objective"),
    )

    mass_plot = Plot(
        [scatter(x=df.time, y=df.sc1_mass, mode="lines+markers", name="CubeSat mass")],
        Layout(title="CubeSat Mass History", xaxis_title="time (s)", yaxis_title="mass (kg)"),
    )

    paths = [
        _save_plot(earth_plot, plot_dir, "earth_centered_trajectory.html"),
        _save_plot(station_plot, plot_dir, "station_centered_tracking.html"),
        _save_plot(tracking_plot, plot_dir, "tracking_error.html"),
        _save_plot(pso_cost_plot, plot_dir, "pso_cost_history.html"),
        _save_plot(mass_plot, plot_dir, "cubesat_mass_history.html"),
    ]
    return paths
end

function run_rpo_cubesat_mpc_demo(; mission_time=180.0, kwargs...)
    demo = build_rpo_cubesat_mpc_demo(; mission_time=mission_time, kwargs...)
    t = @elapsed solution = run_rpo_case(demo)
    println("COMPUTATIONAL TIME = $(t) s")
    plan = rpo_installed_plan(solution)
    settings = solution.prob.p.args.simulation_settings
    csv_path = joinpath(settings.results_directory, "simulation_results.csv")
    csv_path = settings.results && isfile(csv_path) ? csv_path : nothing
    plot_paths = csv_path === nothing ? String[] : _save_rpo_single_case_plots(csv_path, demo, plan)
    return merge(demo, (solution=solution, plan=plan, csv_path=csv_path, plot_paths=plot_paths))
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo = run_rpo_cubesat_mpc_demo()
    println("RPO single-case run complete. Results: $(abspath(demo.csv_path))")
    println("Plots:")
    for path in demo.plot_paths
        println("  ", abspath(path))
    end
end
