# One instance of every concrete built-in force/torque model for the conformance
# tests, with the fixture it runs on, and the types exempted with their reason.
# test/unit/conformance/builtin_conformance_tests.jl fails when a concrete built-in
# type is neither here nor exempted, so a new model is covered when it is added.
# The models that need SPICE ephemerides run on a SPICE fixture when the kernels of
# data/GRAMSuite.jl are present (as in CI); otherwise their checks that need no run
# are applied (`runs=false`).
module ConformanceInstances

using LinearAlgebra: I
using Random: Random
using StaticArrays: SVector
using SpaceAGORA
using SpaceAGORA.Conformance: ConformanceFixture, conformance_fixture

const SM = SpaceAGORA.SimulationModel
const C = SpaceAGORA.Conformance

export force_torque_instances, force_torque_exemptions, spice_available, sensor_instances, sensor_exemptions,
    aerobraking_mpc_controls, gnc_instances, gnc_exemptions, algorithm_instances, algorithm_exemptions, family_instances

const SPICE_PATH = joinpath(pkgdir(SpaceAGORA), "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE")
spice_available() = isdir(SPICE_PATH)

# The built-in fixture with an exponential atmosphere, for the drag models.
function atmosphere_fixture(model)
    planet = SpaceAGORA.make_no_gram_planet(:earth)
    density = SpaceAGORA.make_no_gram_density_model(planet, :exponential)
    return ConformanceFixture(C._fixture_configuration(effectors=(model,), density_model=density))
end

# The spherical-harmonics model reads its own planet, so it gets the fixture's.
function harmonics_instance()
    base = C._fixture_configuration()
    planet = base.environment_model.planet
    file = joinpath(pkgdir(SpaceAGORA), "data", "planets", "gravity", "EarthGGM05C.csv")
    model = SpaceAGORA.GravitationalHarmonicsModel(4, 4, file, planet)
    dynamics = SM.DynamicsModel(base.dynamics_model.spacecraft, (base.dynamics_model.dynamic_effectors..., model))
    return (model=model, fixture=ConformanceFixture(SM.SimConfig._with_configuration(base; dynamics_model=dynamics)),
        runs=true, reason="")
end

function mesh_surrogate_model()
    coefficients = reshape(sqrt(4pi) .* [-1.0, 0.2, 0.1, 0.05, 0.1, -0.07], 6, 1)
    surrogate = SM.MeshAeroSurrogate(0, 0, coefficients, zeros(6, 1), 1.0, 1.0, SVector(0.0, 0.0, 0.0), 1.0, 1.0, 3.0, 20.0)
    return SM.AerodynamicCoefficientMeshSurrogate(Dict(1 => surrogate), true)
end

function wheel_model()
    t = collect(0.0:1.0:60.0)
    speeds = hcat(400.0 .+ 2.0 .* t, -250.0 .+ 0.5 .* t, fill(120.0, length(t)))
    return SM.DynamicEffectors.ReactionWheelMomentum.ReactionWheelMomentumModel(t, speeds, Matrix{Float64}(I, 3, 3), 2.0e-5)
end

# The built-in fixture on SPICE ephemerides; `build` makes the model from the
# fixture's planet, which loads the kernels.
function on_spice(build)
    planet = SM.Earth("", SPICE_PATH)
    model = build(planet)
    fixture = ConformanceFixture(C._fixture_configuration(effectors=(model,), planet=planet,
        ephemerides_model=SM.SpiceEphemeridesModel()))
    return (model=model, fixture=fixture, runs=true, reason="")
end

on_fixture(model) = (model=model, fixture=conformance_fixture(effectors=(model,)), runs=true, reason="")
on_atmosphere(model) = (model=model, fixture=atmosphere_fixture(model), runs=true, reason="")
const SPICE_REASON = "needs SPICE ephemerides (the engine refuses it under SimpleEphemeridesModel before the run), " *
    "and the SPICE kernels of data/GRAMSuite.jl are absent here; checked without a run"
without_run(model) = (model=model, fixture=nothing, runs=false, reason=SPICE_REASON)

"""
A maximum-depletion aerobraking MPC controller for spacecraft 1 and a campaign
supervisor pairing it with an energy-targeting controller.
"""
function aerobraking_mpc_controls()
    config(mode) = SpaceAGORA.AerobrakingMPCConfig(
        mode=mode, bus_reference_area_m2=5.74, controllable_area_m2=5.7, mass_kg=650.0,
        drag_coefficient=2.2, qdot_max_w_cm2=0.29, heat_load_max_j_cm2=40.0, drag_max_n=7.55,
        area_slew_max_m2_s=0.2, use_constraints=true, use_slew_constraint=false,
        use_qdot_constraint=true, use_heat_load_constraint=false, use_drag_constraint=true,
        target_energy_mj_kg=-6.41, area_weight=1.0e-5, area_slew_weight=0.0, slack_weight=1.0e3,
        target_energy_weight=5.0e-10, max_depletion_energy_weight=1.0, osqp_eps_abs=1.0e-5,
        osqp_eps_rel=1.0e-5, osqp_max_iter=10_000)
    controller(mode) = SpaceAGORA.AerobrakingMPCControlModel(
        config=config(mode), spacecraft_index=1, controlled_panel_links=(2, 3),
        min_alpha_rad=1.0e-4, max_alpha_rad=pi / 2, solve_interval_s=10.0,
        build_reference_on_tick=true, qp_max_nodes=40,
        reference=SpaceAGORA.AerobrakingMPCReferenceConfig(
            h_cut_m=300.0e3, max_coast_steps=2_000_000, max_pass_steps=20_000),
        fallback_area_m2=11.44, prediction_latitude_rad=0.0, prediction_longitude_rad=0.0,
        prediction_wind=false, solve_trigger_altitude_m=300.0e3)
    single = controller(SpaceAGORA.MaxEnergyDepletionMode())
    campaign = SpaceAGORA.AerobrakingMPCCampaignControlModel(
        maximum_depletion_control=controller(SpaceAGORA.MaxEnergyDepletionMode()),
        targeting_control=controller(SpaceAGORA.TargetEnergyMode()))
    return single, campaign
end

"""
Fresh instances of the built-in force/torque models, each with its fixture.
"""
function force_torque_instances()
    aero = SM.DynamicEffectors.AerodynamicEffectors
    gravity = SM.DynamicEffectors.GravityEffectors
    perturbations = SM.DynamicEffectors.PerturbationEffectors
    return [
        on_atmosphere(aero.AerodynamicCoefficientConstant()),
        on_atmosphere(mesh_surrogate_model()),
        on_atmosphere(aero.AerodynamicCoefficientNoBallisticFlight()),
        on_atmosphere(SpaceAGORA.AerodynamicCommandedAreaDragModel()),
        on_atmosphere(SpaceAGORA.AerodynamicCoefficientfM()),
        on_atmosphere(SpaceAGORA.TelemetryVerification.ScaledAerodynamicCoefficientfM(SpaceAGORA.AerodynamicCoefficientfM(), 1.1)),
        on_fixture(gravity.ConstantGravityModel()),
        on_fixture(gravity.GravityGradientTorqueModel()),
        on_fixture(SpaceAGORA.InverseSquaredGravityModel()),
        on_fixture(SpaceAGORA.InverseSquaredJ2GravityModel()),
        harmonics_instance(),
        on_fixture(perturbations.EddyCurrentDampingModel(k_e=1.0e4)),
        on_fixture(perturbations.LVLHCascadeAttitudeControlModel(k_out=[0.1, 0.1, 0.1], w_max=0.01,
            k_rate=[1.0, 1.0, 1.0], tau_max=0.1)),
        on_fixture(perturbations.MagneticTorqueRodModel()),
        on_fixture(SpaceAGORA.PlumeSurfaceInteractionModel((; actuators=(; thrust_n=[1000.0])))),
        on_fixture(wheel_model()),
        on_fixture(SpaceAGORA.RobotArmReactionEffector()),
        (spice_available() ? spice_instances() : spice_models_without_run())...,
    ]
end

# Without the kernels: the same models, checked without a run.
spice_models_without_run() = [
    without_run(SpaceAGORA.NBodyGravityModel(body_names=("Moon", "Sun"), primary_body_name="Earth",
        planet=SM.Earth(), body_mus=(4.9028e12, 1.32712440018e20))),
    without_run(SM.DynamicEffectors.PerturbationEffectors.SolarRadiationPressureModel(1.2, 12.0)),
    without_run(SM.DynamicEffectors.PerturbationEffectors.FacetSolarRadiationPressureModel()),
]

# The N-body model declares its bodies as third bodies, so its environment check
# covers a declared third-body sample.
spice_instances() = [
    on_spice(planet -> SpaceAGORA.NBodyGravityModel(body_names=("Moon", "Sun"), primary_body_name="Earth",
        planet=planet, body_mus=(4.9028e12, 1.32712440018e20))),
    on_spice(planet -> SM.DynamicEffectors.PerturbationEffectors.SolarRadiationPressureModel(1.2, 12.0)),
    on_spice(planet -> SM.DynamicEffectors.PerturbationEffectors.FacetSolarRadiationPressureModel()),
]

"""
Concrete built-in force/torque types without an instance, each with the reason (none:
the models that cannot run here are checked without a run).
"""
force_torque_exemptions() = Dict{Any, String}()

"""
Fresh instances of the built-in sensor kinds, with noise so that the seed check
compares streams. They run on the built-in fixture.
"""
sensor_instances() = [
    SpaceAGORA.RateGyro(period_s=0.5, angle_random_walk=1.0e-4, rate_random_walk=1.0e-5, seed=3),
    SpaceAGORA.StarTracker(period_s=1.0, cross_boresight_noise_rad=1.0e-5, about_boresight_noise_rad=1.0e-4, seed=4),
    SpaceAGORA.Magnetometer(period_s=0.5, noise_t=1.0e-8, bias_t=SVector(1.0e-7, 0.0, -1.0e-7), seed=5),
    SpaceAGORA.GNSSReceiver(period_s=1.0, position_noise_m=5.0, velocity_noise_mps=0.05, seed=6),
    # The whole sky in view and no eclipse cut-off, so that every reading is finite.
    SpaceAGORA.SunSensor(period_s=0.5, noise_rad=1.0e-4, field_of_view_rad=π, min_lit_fraction=0.0, seed=7),
    SpaceAGORA.Accelerometer(period_s=0.5, noise_mps2=1.0e-5, bias_mps2=SVector(1.0e-5, 0.0, 0.0), seed=8),
]

"""
Concrete built-in sensor kinds without an instance, each with the reason.
"""
sensor_exemptions() = Dict{Any, String}()

"""
Instances of the built-in guidance and control models. `runs=false` marks a model that
cannot run in the fixture; its checks that need no run are applied, and `reason` says
why. `fixture`, when set, replaces the default one-spacecraft fixture. `BaseThrusterModel`
is not a control effector type but runs in the control model.
"""
function gnc_instances()
    planet = SpaceAGORA.make_no_gram_planet(:earth)
    braking, approach = SpaceAGORA.apollo11_descent_targets()
    apollo = SpaceAGORA.ApolloDescentConfig(reference_radius_m=planet.Rp_e - 2_000.0, site_lat_deg=30.0,
        site_lon_deg=0.0, braking=braking, approach=approach)
    edg = SpaceAGORA.AerobrakingEnergyDepletionConfig()
    edg_state() = SpaceAGORA.AerobrakingEnergyDepletionState(num_sats=1)
    full(kind, model; fixture=nothing) = (kind=kind, model=model, runs=true, reason="", fixture=fixture)
    static(kind, model, reason) = (kind=kind, model=model, runs=false, reason=reason, fixture=nothing)
    # The RPO controller ticks for its chaser and reads its target, spacecraft 1 and 2,
    # so it runs on two spacecraft; in a one-spacecraft run its scope names a spacecraft
    # the run does not have, and the engine refuses the run.
    mpc = SM.ControlHooks.RPOMPCControlModel()
    edg_reason = "energy-depletion guidance and control need the optional SpaceAGORAEDG package, an atmosphere " *
        "and each other; checked without a run"
    apollo_reason = "a lunar descent law that needs its landing-site configuration and terrain; checked without a run"
    mpc_reason = "the KS aerobraking MPC plans an atmospheric pass of a spacecraft with articulated panel links; " *
        "checked without a run"
    aerobraking_mpc, aerobraking_campaign = aerobraking_mpc_controls()
    return [
        full(:guidance_model, SM.GuidanceModels.AerobrakingCampaignPropulsiveManeuverGuidanceModel(
            maneuver_orbit_number=[2, 3], maneuver_Δv=[5.0, -3.0])),
        full(:guidance_model, SpaceAGORA.ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=planet.Rp_e + 800e3)),
        static(:guidance_model, SpaceAGORA.AerobrakingEnergyDepletionGuidanceModel(edg, edg_state()), edg_reason),
        static(:guidance_model, SpaceAGORA.ApolloDescentGuidanceModel(apollo, SpaceAGORA.ApolloDescentState(1)), apollo_reason),
        full(:control_model, SM.ControlHooks.MagneticMomentumManagerModel()),
        # Its declared control period matches the fixture's 1 s control rate.
        full(:control_model, SM.ControlHooks.RobotArmControlEffector(control_dt_s=1.0)),
        full(:control_model, mpc; fixture=conformance_fixture(control=(mpc,), spacecraft_count=2)),
        full(:control_model, SM.ThrusterModels.BaseThrusterModel(thrust=[0.5], direction=[0.0], Δv=[0.1],
            start_burn_time=[5.0], stop_burn_time=[20.0], Isp=[300.0])),
        static(:control_model, SpaceAGORA.AerobrakingEnergyDepletionControlModel(edg, edg_state()), edg_reason),
        static(:control_model, aerobraking_mpc, mpc_reason),
        static(:control_model, aerobraking_campaign, mpc_reason),
        static(:control_model, SpaceAGORA.ApolloDescentControlModel(SpaceAGORA.ApolloDescentControlConfig(), apollo,
            SpaceAGORA.ApolloDescentState(1)), apollo_reason),
    ]
end

"""
Concrete built-in guidance, navigation and control types without an instance, each
with the reason.
"""
gnc_exemptions() = Dict{Any, String}(
    SpaceAGORA.RPOPlannerLifecycle.PlannerGuidance =>
        "the binding that runs an RPO planner in a simulation; it is built by make_rpo_configuration and checked " *
        "with the algorithm contract",
    SpaceAGORA.SolarPanelAngleOfAttackControlModel =>
        "the panel actuator inside AerobrakingEnergyDepletionControlModel, with no control tick of its own",
    SpaceAGORA.AlgorithmGuidance =>
        "the generic binding that runs an algorithm of any family in a simulation; its family's preflight needs " *
        "that family's configuration (the burn-plan family: one burn scheduler sized for the run's spacecraft; " *
        "the joint-reference family: the arm controller that reads its plan buffer), and it is checked with the " *
        "algorithm contract and in test/unit/gnc/gnc_runtime_contract_tests.jl",
)

"""
The hold planner of `examples/rpo_planner_env/smoke.jl` without its HYPR imports: a
user-owned planner, valid only for an unchanged goal, that draws from its stream,
counts its proposals in its state and supports retiming.
"""
struct HoldPlanner <: SpaceAGORA.AbstractRPOPlanner end
SpaceAGORA.planner_capabilities(::HoldPlanner) =
    SpaceAGORA.RPOPlannerCapabilities(state_sources=(:truth,), frames=(:target_rtn,), retiming=true)
SpaceAGORA.initialize_planner(::HoldPlanner, context) = Ref(0)
function SpaceAGORA.plan_rpo!(state, ::HoldPlanner, r::SpaceAGORA.RPOPlanningRequest, rng::Random.AbstractRNG)
    state[] += 1
    reference = SpaceAGORA.RPOReference(t_ref_s=[0.0, r.reference_dt_s],
        r_ref_rtn_m=hcat(collect(r.x_rtn[1:3]), collect(r.goal_rtn_m)), v_ref_rtn_mps=zeros(3, 2),
        origin_time_s=r.time_s, valid_until_s=r.valid_until_s, chaser_id=r.chaser_id, target_id=r.target_id,
        geometry_revision=r.geometry_revision)
    return SpaceAGORA.RPOPlanningResult(request_id=r.request_id, status=:candidate, termination=:completed,
        reference=reference, diagnostics=(calls=state[], draw=rand(rng)))
end
SpaceAGORA.retime_rpo!(state, planner::HoldPlanner, r::SpaceAGORA.RPOPlanningRequest, active::SpaceAGORA.RPOReference,
    rng::Random.AbstractRNG) = SpaceAGORA.plan_rpo!(state, planner, r, rng)

# The small HYPR configuration of examples/rpo_planner_env/smoke.jl.
hypr_config() = SpaceAGORA.RPOPSOConfig(n_waypoints=2, n_particles=8, n_iters=3, adaptive_enable=false,
    adaptive_n_waypoints_min=2, adaptive_n_waypoints_max=2, adaptive_n_particles_min=8, adaptive_n_particles_max=8,
    adaptive_n_iters_min=3, adaptive_n_iters_max=3, sample_ds_m=0.2, retime_accel_limit_enable=true,
    retime_max_speed_mps=0.5, retime_a_max_mps2=0.5 * 0.05 / 5.2, iteration_runtime_limit_s=Inf,
    rrt_warmstart_enable=true)

# The planners that need the optional HYPR package, each as the comparison tests
# configure it (test/unit/gnc/rpo_comparison_planner_object_tests.jl), small enough
# for a short run, with the retimer limits of the RPO pilot (as there for a run
# through the RPO binding), so that they plan the pilot's requests instead of
# refusing them.
function hypr_planners()
    G = SM.GuidanceHooks
    base = G.rpo_740_mpc_final_pso_config(safe_distance_m=0.2; n_particles=8, n_iters=3, n_waypoints=2,
        adaptive_enable=false, reexplore_enable=false, cull_enable=false, sample_ds_m=0.25, retime_dt_s=0.5,
        retime_a_max_mps2=0.05, retime_max_steps=5_000, probe_candidates=4, refinement_rounds=2,
        refinement_waypoint_passes=2)
    pso = SpaceAGORA.rpo_pso_config(base; retime_accel_limit_enable=true, retime_max_speed_mps=0.5,
        retime_a_max_mps2=0.5 * 0.05 / 5.2, retime_dt_s=0.1)
    optimizer = G.RPOTrajectoryOptimizerSettings(match_hypr_runtime=false, runtime_limit_s=Inf, force_full_iters=true)
    rrt_connect = G.RPORRTConnectSettings(n_iters=120, step_size_m=0.5, shortcut_iters=6)
    return [
        SpaceAGORA.HYPRRPOPlanner(hypr_config()),
        SpaceAGORA.HYPRComparisonPlanner(pso),
        SpaceAGORA.UnrefinedPSOPlanner(pso),
        SpaceAGORA.RRTConnectPlanner(rrt_connect; config=pso, optimizer),
        SpaceAGORA.RRTConnectBezierPlanner(rrt_connect; config=pso, optimizer),
        SpaceAGORA.RRTStarPlanner(G.RPORRTStarSettings(n_iters=90, step_size_m=0.5, shortcut_iters=6); config=pso,
            optimizer),
        SpaceAGORA.CHOMPPlanner(G.RPOCHOMPSettings(n_iters=3); config=pso, optimizer, seed=rrt_connect),
        SpaceAGORA.STOMPPlanner(G.RPOSTOMPSettings(n_iters=3, n_rollouts=4); config=pso, optimizer, seed=rrt_connect),
        # A retained path from the conformance requests' start to their goal.
        SpaceAGORA.RetainedPathRPOPlanner([3.0 4.0 5.0; 0.0 0.0 0.0; 0.0 0.0 0.0], hypr_config()),
    ]
end
const HYPR_PLANNER_TYPES = (SpaceAGORA.HYPRRPOPlanner, SpaceAGORA.HYPRComparisonPlanner, SpaceAGORA.UnrefinedPSOPlanner,
    SpaceAGORA.RRTConnectPlanner, SpaceAGORA.RRTConnectBezierPlanner, SpaceAGORA.RRTStarPlanner, SpaceAGORA.CHOMPPlanner,
    SpaceAGORA.STOMPPlanner, SpaceAGORA.RetainedPathRPOPlanner)

include(joinpath(@__DIR__, "edg_guidance_conformance.jl"))

"""
The planning algorithms checked with `check_algorithm`, each with its requests, the
trusted queries its family's validation needs and a short simulation: the direct
planner (with the station clearance `make_rpo_configuration` gives it), the hold
planner above, the burn-plan family's apoapsis planner, the joint-reference
family's quintic planner, the EDG guidance decision when the optional
SpaceAGORAEDG package is loaded, and, when the optional HYPR package is
loaded, the HYPR planner, the seven comparison planners and the retained-path planner.
"""
function algorithm_instances()
    direct = SpaceAGORA.DirectRPOPlanner(segment_clearance_at=SpaceAGORA.RPOPlannerLifecycle.station_segment_clearance)
    rpo(algorithm) = (algorithm=algorithm, requests=C.rpo_conformance_requests(), trusted=C.rpo_trusted_queries(),
        simulation=C.rpo_conformance_simulation())
    entries = [
        rpo(direct),
        (algorithm=HoldPlanner(), requests=C.rpo_conformance_requests(goal_rtn_m=nothing),
            trusted=C.rpo_trusted_queries(), simulation=C.rpo_conformance_simulation(goal_rtn_m=(3.0, 0.0, 0.0))),
        (algorithm=SpaceAGORA.ApoapsisBurnPlanner(delta_v_mps=10.0), requests=C.burn_plan_conformance_requests(),
            trusted=NamedTuple(), simulation=C.burn_plan_conformance_simulation()),
        (algorithm=SpaceAGORA.QuinticJointPlanner(), requests=C.joint_reference_conformance_requests(),
            trusted=NamedTuple(), simulation=C.joint_reference_conformance_simulation()),
    ]
    SpaceAGORA.EDGServices.available() && push!(entries, (algorithm=SpaceAGORA.EDGGuidanceAlgorithm(),
        requests=edg_guidance_conformance_requests(), trusted=NamedTuple(),
        simulation=edg_guidance_conformance_simulation()))
    SpaceAGORA.hypr_available() && append!(entries, rpo.(hypr_planners()))
    return entries
end

"""
Concrete built-in planning algorithm types without an instance, each with the reason.
"""
function algorithm_exemptions()
    exemptions = SpaceAGORA.hypr_available() ? Dict{Any, String}() : Dict{Any, String}(
        T => "needs the optional HYPR package, which is not loaded here; the unit runner loads it, and " *
            "there the planner is checked" for T in HYPR_PLANNER_TYPES)
    SpaceAGORA.EDGServices.available() || (exemptions[SpaceAGORA.EDGGuidanceAlgorithm] =
        "needs the optional SpaceAGORAEDG package, which is not loaded here; the unit runner loads it, and there " *
        "the algorithm is checked")
    return exemptions
end

"""
The algorithm families checked with `check_algorithm_family`, under each rejection
policy they allow, with an algorithm of the family and its requests (which include
one the family refuses: for RPO a retiming request the direct planner's capability
check refuses, for the burn plan a request during an installed burn, for the joint
reference a request with an obstacle on the planned path).
"""
function family_instances()
    direct = SpaceAGORA.DirectRPOPlanner(segment_clearance_at=SpaceAGORA.RPOPlannerLifecycle.station_segment_clearance)
    rpo = [(family=SpaceAGORA.RPOPlannerInterfaces.RPOReferenceFamily(on_rejection=policy), algorithm=direct,
        requests=C.rpo_conformance_requests(), link=SM.RPOPlanBuffer(), trusted=C.rpo_trusted_queries())
        for policy in (:fail, :install_fallback)]
    # The burn-plan family has one policy, :keep_previous; its requests include one it refuses.
    burn = (family=SpaceAGORA.BurnPlanFamily(), algorithm=SpaceAGORA.ApoapsisBurnPlanner(delta_v_mps=10.0),
        requests=C.burn_plan_conformance_requests(),
        link=SpaceAGORA.BurnPlanSlot(SM.PropulsiveBurnPlan[SM.PropulsiveBurnPlan()], 1), trusted=NamedTuple())
    # The joint-reference family has one policy, :fail; its requests include one it refuses.
    joint = (family=SpaceAGORA.JointReferenceFamily(), algorithm=SpaceAGORA.QuinticJointPlanner(),
        requests=C.joint_reference_conformance_requests(), link=SpaceAGORA.RobotArmPlanBuffer(), trusted=NamedTuple())
    # The EDG guidance family (:fail) needs the SpaceAGORAEDG provider, which the unit runner loads.
    edg = (family=SpaceAGORA.EDGGuidanceFamily(), algorithm=SpaceAGORA.EDGGuidanceAlgorithm(),
        requests=edg_guidance_conformance_requests(), link=SpaceAGORA.AerobrakingEnergyDepletionState(num_sats=1),
        trusted=NamedTuple())
    return [rpo..., burn, joint, edg]
end

end # module ConformanceInstances
