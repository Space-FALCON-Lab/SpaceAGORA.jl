# Worked example E of "Writing a GNC Algorithm" (docs/src/user/gnc_planning.md):
# a planner of the joint-reference family, run through the generic binding
# `AlgorithmGuidance` with a robot arm's controller. Each region between a
# `# snippet:` line and `# end snippet` appears on that page as it is here; the HPC
# and extensibility docs gate compares them.
using SpaceAGORA
using SpaceAGORA.Conformance: test_algorithm, joint_reference_conformance_requests,
    joint_reference_conformance_simulation
using LinearAlgebra: norm
using Test

# snippet: planner
# A joint-space path to the inverse-kinematics solution for the target, each joint
# following the smooth step 3s² - 2s³ of normalized time s. It ignores obstacles;
# the family's validation refuses a path that does not clear them.
struct SmoothStepJointPlanner <: AbstractJointReferencePlanner end
SpaceAGORA.consumes_rng(::SmoothStepJointPlanner) = false

function SpaceAGORA.propose!(state, planner::SmoothStepJointPlanner, request::JointReferenceRequest, rng)
    arm, base = request.model, request.base_pose
    q_goal = try
        robot_arm_ik(arm, base, request.target; q_seed=request.q_start)
    catch err
        err isa ErrorException || rethrow()
        return AlgorithmResult(request_id=request.request_id, status=:infeasible,
            termination=:target_unreachable)
    end
    T, Δq = request.duration_s, q_goal - request.q_start
    # The request's grid: every dt_s from 0, the last step at most dt_s, to T.
    times = [min(k * request.dt_s, T) for k in 0:ceil(Int, T / request.dt_s)]
    s = times ./ T
    q_ref = request.q_start .+ Δq * (3s .^ 2 .- 2s .^ 3)'
    dq_ref = Δq * ((6s .- 6s .^ 2) ./ T)'
    ddq_ref = Δq * ((6 .- 12s) ./ T^2)'
    ee_ref = reduce(hcat, [robot_arm_fk(arm, base, q).end_effector_position for q in eachcol(q_ref)])
    final_error_m = norm(ee_ref[:, end] - request.target)
    plan = RobotArmPlan(arm, base, times, q_ref, dq_ref, ddq_ref, ee_ref, copy(request.q_start),
                        q_goal, request.target, final_error_m, :smooth_step)
    return AlgorithmResult(request_id=request.request_id, status=:candidate,
        termination=:completed, product=plan)
end
# end snippet

# snippet: run
# The default arm on a spacecraft's bus. The plan buffer is the link between the
# binding, which installs accepted plans, and the arm's controller, which tracks
# them; it starts with the arm held at zero joint angles.
arm = default_robot_arm_model()
base = RobotArmBasePose([0.0, 0.0, 0.0])
target = robot_arm_fk(arm, base, [0.25, -0.35, 0.20]).end_effector_position
buffer = RobotArmPlanBuffer(arm, base, zeros(3))
controller = RobotArmControlEffector(spacecraft_index=1, plan_buffer=buffer, control_dt_s=0.5,
                                     controller=init_robot_arm_joint_mpc(arm; dt_s=0.5, horizon=4))
settings = JointReferenceSettings(model=arm, base_pose=base, target=target, dt_s=0.5, duration_s=4.0)
binding = AlgorithmGuidance(algorithm=SmoothStepJointPlanner(), spacecraft_index=1, settings=settings,
                            link=buffer, events=[AlgorithmEvent(0.0)])

planet = make_no_gram_planet(:earth)
root = Link(root=true, m=100.0, ref_area=2.0)
orbit = InitialCondition(ra=planet.Rp_e + 2000e3, rp=planet.Rp_e + 400e3, i=30.0, ω=0.0, Ω=0.0, ν=10.0)
spacecraft = SpacecraftModel(joints=Joint[], links=[root], root=root, prop_mass=0.0,
                             inertia_tensor=root.inertia, initial_condition=orbit, id=1)
config = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime, mission_time=1.0, number_of_orbits=1,
        keplerian=true, orientation_sim=false, num_steps_to_save=10),
    environment_model=EnvironmentModel(
        planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false),
    dynamics_model=DynamicsModel([spacecraft], (InverseSquaredGravityModel(),)),
    guidance_model=GuidanceModel(guidance_effectors=(binding,), guidance_rates=[60.0]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(controller,), control_rates=[0.5]),
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances())

solution = run_simulation(config; return_solution=true)
installed = only(filter(record -> record.event === :installed,
                        only(algorithm_run_report(solution).bindings).records))
println("installed request ", installed.request_id, ": ", installed.validation.reason,
        ", final end-effector error ", installed.validation.metrics.final_error_m, " m")
@test installed.validation.accepted
# end snippet

# snippet: conformance
@testset "SmoothStepJointPlanner conformance" begin
    test_algorithm(SmoothStepJointPlanner(), joint_reference_conformance_requests();
                   simulation=joint_reference_conformance_simulation())
end
# end snippet
