# The joint-reference family: its representative requests, its conformance link and a
# short simulation through the generic binding (`AlgorithmGuidance`). Its validation
# needs no trusted queries.

const _JOINT = JointReferencePlanning

# The default arm, its base at the origin, and the target its joints (0.25, -0.35,
# 0.20) rad reach.
function _joint_reference_arm()
    S = SimulationModel
    arm = S.default_cloth_arm_model()
    base = S.ClothArmBasePose([0.0, 0.0, 0.0])
    return (model=arm, base_pose=base, target=S.cloth_fk(arm, base, [0.25, -0.35, 0.20]).end_effector_position)
end

"""
    joint_reference_conformance_requests() -> Vector{JointReferenceRequest}

Representative requests of the joint-reference family for [`check_algorithm`](@ref),
for the default robot arm (`default_robot_arm_model()`) with its base at the origin,
from rest at zero joint angles to a reachable target over 4 s sampled every 0.5 s:
an initial request at 0 s; a forced request at 5 s with a 2 cm sphere on the
straight joint-space path, which a planner that does not avoid it proposes through
and the family refuses (`:insufficient_clearance`) and, by its `:fail` policy, fails
the request for; and a replanning request at 10 s with the sphere moved clear of the
path. The family's validation needs no trusted queries.
"""
function joint_reference_conformance_requests()
    arm = _joint_reference_arm()
    ask(; kwargs...) = _JOINT.JointReferenceRequest(; spacecraft_index=1, model=arm.model, base_pose=arm.base_pose,
        q_start=zeros(3), target=arm.target, dt_s=0.5, duration_s=4.0, kwargs...)
    on_path = SimulationModel.RobotArmSphereObstacle([0.4536, 0.0570, 0.0669], 0.02)
    clear = SimulationModel.RobotArmSphereObstacle([-0.3, 0.0, 0.0], 0.02)
    return [ask(request_id=1, reason=:initial, time_s=0.0),
        ask(request_id=2, reason=:forced, time_s=5.0, obstacles=[on_path]),
        ask(request_id=3, reason=:replan, time_s=10.0, obstacles=[clear])]
end

conformance_link(::_JOINT.JointReferenceFamily) = SimulationModel.RobotArmPlanBuffer()

"""
    joint_reference_conformance_simulation(; mission_time_s=1.0) -> Function

A `simulation` for [`check_algorithm`](@ref): a function of a joint-reference planner
that runs it through `AlgorithmGuidance`, with one request at the run's start, on
one spacecraft in a 400 by 2000 km Earth orbit whose `RobotArmControlEffector`
tracks the default arm's plans from its plan buffer at a 0.5 s control step, the
buffer starting with a hold at zero joint angles. The
run is made twice in sequence and, with two or more threads, twice concurrently,
and the binding's run records of the first run are returned. It throws when a run
fails or the runs differ. The arm's coupled compliant state needs short solver
steps, so the default run is short: it covers the request and install at the
start and the controller's ticks at 0, 0.5 and 1 s.
"""
joint_reference_conformance_simulation(; mission_time_s::Real=1.0) =
    planner -> _joint_reference_simulation_records(planner, Float64(mission_time_s))

function _joint_reference_configuration(planner, mission_time_s::Float64)
    S = SimulationModel
    arm = _joint_reference_arm()
    planet = S.make_no_gram_planet(:earth)
    root = S.Link(root=true, m=100.0, ref_area=2.0)
    initial = S.InitialCondition(ra=planet.Rp_e + 2000e3, rp=planet.Rp_e + 400e3, i=30.0, ω=0.0, Ω=0.0, ν=10.0)
    spacecraft = S.SpacecraftModel(joints=S.Joint[], links=[root], root=root, prop_mass=50.0,
        inertia_tensor=root.inertia, initial_condition=initial, id=1)
    buffer = S.RobotArmPlanBuffer(arm.model, arm.base_pose, zeros(3))
    controller = S.RobotArmControlEffector(spacecraft_index=1, plan_buffer=buffer, control_dt_s=0.5,
        controller=S.init_robot_arm_joint_mpc(arm.model; dt_s=0.5, horizon=4))
    settings = _JOINT.JointReferenceSettings(; model=arm.model, base_pose=arm.base_pose, target=arm.target,
        dt_s=0.5, duration_s=4.0)
    binding = AlgorithmLifecycle.AlgorithmGuidance(; algorithm=planner, spacecraft_index=1, settings=settings,
        link=buffer, events=[AlgorithmLifecycle.AlgorithmEvent(0.0)])
    environment = S.EnvironmentModel(planet=planet, EI=120.0, density_model=S.NoAtmosphereModel(),
        ephemerides_model=S.SimpleEphemeridesModel(),
        thermal_model=S.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet), topography=false, wind=false)
    return S.SimulationConfiguration(
        simulation_settings=S.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=S.MissionConfiguration(mission_type=S.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_time_s, orientation_sim=false, num_steps_to_save=100),
        environment_model=environment,
        dynamics_model=S.DynamicsModel(S.SpacecraftModel[spacecraft], (S.InverseSquaredGravityModel(),)),
        guidance_model=S.GuidanceModel((binding,), [1.0e6]),
        navigation_model=S.NavigationModel((), Float64[]),
        control_model=S.ControlModel(control_effectors=(controller,), control_rates=[0.5]),
        initial_time=S.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

function _joint_reference_simulation_records(planner, mission_time_s::Float64)
    configuration = _joint_reference_configuration(planner, mission_time_s)
    base = _solve(configuration)
    _failed(base) && throw(base.error)
    again = _solve(configuration)
    gap = _model_run_difference(_with_controller(base), _with_controller(again))
    gap === nothing || error("two runs of the joint-reference configuration differ: $(gap)")
    concurrency = _concurrency_check(configuration, base, again, FIRST_RELEASE)
    concurrency.status === :fail && error(concurrency.detail)
    return only(AlgorithmLifecycle.algorithm_run_report(base.solution).bindings).records
end

# A run's outcome with the arm controller as it ended the run, for comparing its
# held command and installed plan as well as the saved states.
function _with_controller(run::_RunOutcome)
    solution = run.solution
    controller = only(solution.prob.p.args.control_model.control_effectors)
    return _RunOutcome((t=solution.t, u=solution.u, model=deepcopy(controller)), nothing)
end
