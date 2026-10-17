# The joint-reference family and its link to the robot-arm controller: trusted
# validation reasons, the quintic planner, runs through the generic binding
# (`AlgorithmGuidance`) and the shared request driver, refusals under the `:fail`
# policy, and the plan buffer's link checks.
module JointReferenceFamilyTests
using Test, Random, StaticArrays, LinearAlgebra, SpaceAGORA, SpaceAGORA.SimulationModel
const GA = SpaceAGORA.GNCAlgorithms
const SM = SpaceAGORA.SimulationModel
const C = SpaceAGORA.Conformance

const PLANET = SimulationModel.make_no_gram_planet(:earth)
const ARM = SM.default_robot_arm_model()
const BASE = SM.RobotArmBasePose([0.0, 0.0, 0.0])
const TARGET = SM.robot_arm_fk(ARM, BASE, [0.25, -0.35, 0.20]).end_effector_position
const DT, DURATION = 0.5, 4.0

settings(; kwargs...) = JointReferenceSettings(; merge((model=ARM, base_pose=BASE, target=TARGET, dt_s=DT,
    duration_s=DURATION), (; kwargs...))...)
request(; kwargs...) = JointReferenceRequest(; merge((request_id=1, reason=:forced, spacecraft_index=1, time_s=0.0,
    model=ARM, base_pose=BASE, q_start=zeros(3), target=TARGET, dt_s=DT, duration_s=DURATION), (; kwargs...))...)
planned(r) = GA.propose!(nothing, QuinticJointPlanner(), r, MersenneTwister(1))
direct_plan() = SM.plan_robot_arm_motion(ARM, BASE, zeros(3), TARGET; config=SM.RobotArmPlannerConfig(dt_s=DT, duration_s=DURATION))
reason(r, result) = validate_joint_reference(r, result).reason

# A copy of `plan` with some fields replaced.
with_plan(plan; kwargs...) = SM.RobotArmPlan((get(kwargs, f, getfield(plan, f)) for f in fieldnames(SM.RobotArmPlan))...)
candidate(r, plan) = AlgorithmResult(request_id=r.request_id, status=:candidate, termination=:planned, product=plan)

function arm_sat(; id=7)
    root = Link(root=true, m=100.0, ref_area=2.0)
    ic = InitialCondition(ra=PLANET.Rp_e + 2000e3, rp=PLANET.Rp_e + 400e3, i=30.0, ω=0.0, Ω=0.0, ν=10.0 + id)
    return SpacecraftModel(joints=Joint[], links=[root], root=root, prop_mass=50.0, inertia_tensor=root.inertia,
        initial_condition=ic, id=id)
end

mpc() = SM.init_robot_arm_joint_mpc(ARM; dt_s=DT, horizon=4)
arm_controller(; kwargs...) = SM.RobotArmControlEffector(; merge((spacecraft_index=1, control_dt_s=DT, controller=mpc()),
    (; kwargs...))...)
binding(buffer; kwargs...) = AlgorithmGuidance(; merge((algorithm=QuinticJointPlanner(), spacecraft_index=1,
    settings=settings(), link=buffer, events=[AlgorithmEvent(0.0)]), (; kwargs...))...)

# Ends at 3 s, before the 4 s reference ends, so the last tick is during the motion.
function arm_config(; controls=(arm_controller(),), guidance=(), guidance_rate=1.0e6, mission_time=3.0,
        spacecraft=SpacecraftModel[arm_sat()])
    env = EnvironmentModel(planet=PLANET, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=PLANET), topography=false, wind=false)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=false, num_steps_to_save=50),
        environment_model=env,
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(Tuple(guidance), fill(Float64(guidance_rate), length(guidance))),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel(control_effectors=Tuple(controls), control_rates=fill(DT, length(controls))),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

# A buffer holding the arm still at zero joint angles, its controller and its
# planner, in one configuration.
hold_buffer() = RobotArmPlanBuffer(ARM, BASE, zeros(3))
function linked_config(; buffer=hold_buffer(), binding_kwargs=(), guidance_rate=1.0e6)
    b = binding(buffer; binding_kwargs...)
    return arm_config(controls=(arm_controller(plan_buffer=buffer),), guidance=(b,), guidance_rate=guidance_rate), b, buffer
end

ran_controller(solution) = only(solution.prob.p.args.control_model.control_effectors)
same_states(a, b) = a.t == b.t && length(a.u) == length(b.u) && all(collect(x) == collect(y) for (x, y) in zip(a.u, b.u))
failure(f) = try f(); nothing catch e e end

@testset "Joint-reference family declarations" begin
    planner = QuinticJointPlanner()
    family = algorithm_family(planner)
    @test family === JointReferenceFamily() && planner isa AbstractJointReferencePlanner <: AbstractGNCAlgorithm
    @test on_rejection(family) === :fail
    @test family_contract_version(family) == v"1.0.0"
    @test algorithm_result_type(family) === AlgorithmResult
    @test !consumes_rng(planner)
    @test_throws ArgumentError QuinticJointPlanner(ik_tol_m=0.0)
    @test_throws ArgumentError QuinticJointPlanner(ik_max_iters=0)
    @test_throws ArgumentError settings(dt_s=0.0)
    @test_throws ArgumentError settings(clearance_m=-1.0)
    @test_throws ArgumentError settings(target=(1.0, 2.0))
    @test_throws ArgumentError request(request_id=0)
    @test_throws ArgumentError request(q_start=[0.0, NaN, 0.0])
    @test_throws ArgumentError request(goal_tolerance_m=0.0)
    # The link: one reader, at most one writer; a hold keeps the arm still.
    @test link_cardinality(RobotArmPlanBuffer()) == (producers=0:1, consumers=1:1)
    hold = hold_buffer()
    @test hold.updated_at_s == 0.0 && all(iszero, hold.plan.q_ref) && all(iszero, hold.plan.dq_ref)
    buffer = RobotArmPlanBuffer()
    @test gnc_links(arm_controller()) == ()
    @test only(gnc_links(arm_controller(plan_buffer=buffer))).link === buffer
    @test only(gnc_links(binding(buffer))) isa Produces && only(gnc_links(binding(buffer))).link === buffer
    @test gnc_links(AlgorithmGuidance(algorithm=ApoapsisBurnPlanner(delta_v_mps=1.0), spacecraft_index=1)) == ()
end

@testset "The quintic planner proposes plan_robot_arm_motion's reference" begin
    r = request()
    rng = MersenneTwister(2026); witness = copy(rng)
    result = GA.propose!(nothing, QuinticJointPlanner(), r, rng)
    @test rand(rng, UInt64) == rand(witness, UInt64)   # draws nothing
    @test result.status === :candidate && result.termination === :planned
    @test C._same_state(result.product, direct_plan())
    @test result.diagnostics.final_error_m == result.product.final_error_m
    @test validate_joint_reference(r, result).accepted
    unreachable = planned(request(target=(2.0, 0.0, 0.0)))
    @test unreachable.status === :infeasible && unreachable.termination === :target_unreachable
    @test unreachable.product === nothing
end

@testset "Trusted validation reasons" begin
    r = request()
    good = planned(r)
    plan = good.product
    accepted = validate_joint_reference(r, good)
    @test accepted.reason === :accepted && accepted.metrics.samples == length(plan.t_ref_s)
    @test keys(accepted.metrics) == (:min_clearance_m, :final_error_m, :samples)
    stale = AlgorithmResult(request_id=2, status=:candidate, termination=:planned, product=plan)
    mismatch = validate_joint_reference(r, stale)
    @test mismatch.reason === :request_mismatch && mismatch.metrics == (result_request_id=2, request_id=1)
    @test reason(r, AlgorithmResult(request_id=2, status=:infeasible, termination=:x)) === :request_mismatch
    @test reason(r, AlgorithmResult(request_id=1, status=:infeasible, termination=:x)) === :infeasible
    @test reason(r, AlgorithmResult(request_id=1, status=:unsupported, termination=:x)) === :unsupported
    @test reason(r, candidate(r, (q=zeros(3),))) === :invalid_product
    other_arm = SM.default_robot_arm_model(link_lengths_m=(0.2, 0.16, 0.12))
    @test reason(r, candidate(r, with_plan(plan; model=other_arm))) === :model_mismatch
    @test reason(r, candidate(r, with_plan(plan; base_pose=SM.RobotArmBasePose([0.0, 0.0, 0.1])))) === :base_pose_mismatch
    @test reason(r, candidate(r, with_plan(plan; q_ref=plan.q_ref[:, 1:end-1]))) === :invalid_shape
    @test reason(r, candidate(r, with_plan(plan; q_start=zeros(2)))) === :invalid_shape
    nan_rates = copy(plan.dq_ref); nan_rates[2, 3] = NaN
    @test reason(r, candidate(r, with_plan(plan; dq_ref=nan_rates))) === :nonfinite_reference
    uneven = copy(plan.t_ref_s); uneven[3] += 0.1
    @test reason(r, candidate(r, with_plan(plan; t_ref_s=uneven))) === :invalid_timing
    @test reason(request(duration_s=5.0), candidate(request(duration_s=5.0), plan)) === :invalid_timing
    @test reason(request(dt_s=0.25), candidate(request(dt_s=0.25), plan)) === :invalid_timing
    @test reason(request(q_start=[0.01, 0.0, 0.0]), candidate(request(q_start=[0.01, 0.0, 0.0]), plan)) === :start_mismatch
    beyond = copy(plan.q_ref); beyond[1, 4] = ARM.joints[1].upper_rad + 0.01
    limit = validate_joint_reference(r, candidate(r, with_plan(plan; q_ref=beyond)))
    @test limit.reason === :joint_limit && limit.metrics.joint == 1 && limit.metrics.time_s == plan.t_ref_s[4]
    @test reason(request(target=TARGET .+ 0.05), candidate(request(target=TARGET .+ 0.05), plan)) === :goal_not_reached
    blocked = request(obstacles=[SM.RobotArmSphereObstacle(plan.ee_ref[:, 5], 0.02)])
    clearance = validate_joint_reference(blocked, candidate(blocked, plan))
    @test clearance.reason === :insufficient_clearance && clearance.metrics.min_clearance_m < 0
    far = SM.RobotArmSphereObstacle([-0.3, 0.0, 0.0], 0.02)
    @test reason(request(obstacles=[far]), candidate(request(obstacles=[far]), plan)) === :accepted
    @test reason(request(obstacles=[far], clearance_m=0.5), candidate(request(obstacles=[far], clearance_m=0.5), plan)) ===
        :insufficient_clearance
    # Validation reads copies only: the request and the proposal are unchanged.
    r_copy, good_copy = deepcopy(r), deepcopy(good)
    validate_joint_reference(r, good)
    @test C._same_state(r, r_copy) && C._same_state(good, good_copy)
    @test joint_reference_is_current(request(time_s=5.0), plan, 5.0)
    @test joint_reference_is_current(request(time_s=5.0), plan, 500.0)
    @test !joint_reference_is_current(request(time_s=5.0), plan, 4.0)
end

@testset "Family through the generic binding equals the controller given the plan" begin
    # The controller tracks a plan configured on it; the binding plans the same
    # reference through the request driver and installs it in the controller's buffer.
    own = run_simulation(arm_config(controls=(arm_controller(plan=direct_plan()),)); return_solution=true)
    config, b, buffer = linked_config()
    viafamily = run_simulation(config; return_solution=true)
    @test same_states(own, viafamily)
    held(solution) = held_command(ran_controller(solution), 1)
    @test isequal(held(own), held(viafamily)) && !isempty(held(viafamily).joint_torque_nm)
    @test any(!iszero, held(viafamily).joint_torque_nm)
    report = algorithm_run_report(viafamily).bindings[1]
    @test [r.event for r in report.records] == [:request, :installed]
    @test report.records[1].request isa JointReferenceRequest && report.records[1].request.time_s == 0.0
    @test C._same_state(report.product, direct_plan())
    @test report.records[2].validation.accepted && report.records[2].validation.reason === :accepted
    @test report.stream_words == gnc_stream_words("AlgorithmGuidance/v1", 0, :JointReferenceFamily, 7)
    installed = ran_controller(viafamily).plan_buffer
    @test installed.updated_at_s == 0.0 && C._same_state(installed.plan, direct_plan())
    # The arm's coupled state follows the installed plan, as it follows a configured one.
    @test hasproperty(viafamily.u[end].sc[1], :arm_q)
    # isolate_state=true: the configured objects never run.
    @test b.runtime === nothing && buffer.plan.planner === :cloth_quintic && all(iszero, buffer.plan.q_ref)
end

@testset "A reference installed at a later tick is tracked from then" begin
    config, _, _ = linked_config(binding_kwargs=(events=[AlgorithmEvent(1.0)], ), guidance_rate=DT)
    later = run_simulation(config; return_solution=true)
    report = algorithm_run_report(later).bindings[1]
    @test [(r.event, r.time_s) for r in report.records] == [(:request, 1.0), (:installed, 1.0)]
    installed = ran_controller(later).plan_buffer
    @test installed.updated_at_s == 1.0
    # The request starts from the arm's measured joints, where its coupled state has
    # settled under the hold, and the installed plan from there.
    start = report.records[1].request.q_start
    @test maximum(abs, start) < 1e-3 && installed.plan.q_ref[:, 1] == start
    held_still = run_simulation(arm_config(controls=(arm_controller(plan_buffer=hold_buffer()),)); return_solution=true)
    @test !isequal(held_command(ran_controller(later), 1), held_command(ran_controller(held_still), 1))
end

@testset "Two runs with isolate_state=false start equal" begin
    config, b, buffer = linked_config()
    first_run = run_simulation(config; isolate_state=false, return_solution=true)
    first_end = deepcopy((buffer, held_command(only(config.control_model.control_effectors), 1)))
    @test buffer.plan !== nothing && b.runtime.request_id == 1
    second_run = run_simulation(config; isolate_state=false, return_solution=true)
    @test same_states(first_run, second_run)
    @test C._same_state((buffer, held_command(only(config.control_model.control_effectors), 1)), first_end)
    @test b.runtime.request_id == 1   # fresh run state, not a second request on top of the first
end

@testset "A buffer no model writes keeps the plan it was configured with" begin
    plan = direct_plan()
    own = run_simulation(arm_config(controls=(arm_controller(plan=plan),)); return_solution=true)
    filled = RobotArmPlanBuffer(plan=plan, updated_at_s=0.0)
    unproduced = run_simulation(arm_config(controls=(arm_controller(plan_buffer=filled),)); isolate_state=false,
        return_solution=true)
    @test filled.plan === plan
    @test isequal(held_command(ran_controller(own), 1), held_command(ran_controller(unproduced), 1))
end

# Planners with one defect each: a plan through a joint limit, a time grid of
# another step, and the first result returned for every later request.
struct OverreachPlanner <: AbstractJointReferencePlanner end
function GA.propose!(state, ::OverreachPlanner, r::JointReferenceRequest, rng)
    plan = GA.propose!(state, QuinticJointPlanner(), r, rng).product
    q = copy(plan.q_ref); q[2, end] = ARM.joints[2].lower_rad - 0.1
    return candidate(r, with_plan(plan; q_ref=q))
end
struct CoarsePlanner <: AbstractJointReferencePlanner end
GA.propose!(state, ::CoarsePlanner, r::JointReferenceRequest, rng) =
    candidate(r, SM.plan_robot_arm_motion(r.model, r.base_pose, r.q_start, r.target;
        config=SM.RobotArmPlannerConfig(dt_s=2 * r.dt_s, duration_s=r.duration_s)))
struct CachedResultPlanner <: AbstractJointReferencePlanner end
GA.initialize_algorithm(::CachedResultPlanner, run_context) = Ref{Any}(nothing)
function GA.propose!(state, ::CachedResultPlanner, r::JointReferenceRequest, rng)
    state[] === nothing && (state[] = GA.propose!(nothing, QuinticJointPlanner(), r, rng))
    return state[]
end

@testset "Negative controls: a refused reference stops the run (:fail)" begin
    plan = direct_plan()
    on_path = SM.RobotArmSphereObstacle(plan.ee_ref[:, 5], 0.02)
    cases = [
        (binding_kwargs=(settings=settings(obstacles=[on_path]),), reason=:insufficient_clearance),
        (binding_kwargs=(algorithm=OverreachPlanner(),), reason=:joint_limit),
        (binding_kwargs=(algorithm=CoarsePlanner(),), reason=:invalid_timing),
        (binding_kwargs=(settings=settings(target=(2.0, 0.0, 0.0)),), reason=:infeasible),
        (binding_kwargs=(algorithm=CachedResultPlanner(), events=[AlgorithmEvent(0.0), AlgorithmEvent(1.0; token=2)]),
            reason=:request_mismatch),
    ]
    for case in cases
        config, _, buffer = linked_config(binding_kwargs=case.binding_kwargs, guidance_rate=DT)
        hold = buffer.plan
        err = failure(() -> run_simulation(config; isolate_state=false))
        @testset let reason = case.reason
            @test err isa AlgorithmGuidanceError && err isa GNCAlgorithmError
            @test err.reason === case.reason && err.spacecraft_index == 1
            @test last(err.report.records).event === :failure
            @test !any(r -> r.event === :rejected, err.report.records)   # :fail records no rejection
            if case.reason === :request_mismatch
                # The first request's reference stays installed; the second failed.
                @test err.request_id == 2 && C._same_state(buffer.plan, plan)
                @test err.report.records[end - 1].event === :request
            else
                @test err.request_id == 1 && buffer.plan === hold   # nothing installed
            end
        end
    end
end

@testset "Configuration checks before the run" begin
    buffer = hold_buffer()
    two = SpacecraftModel[arm_sat(id=7), arm_sat(id=8)]
    # An empty buffer: the arm's state has no plan to start from.
    empty = RobotArmPlanBuffer()
    err = failure(() -> run_simulation(arm_config(controls=(arm_controller(plan_buffer=empty),), guidance=(binding(empty),))))
    @test err isa ArgumentError && occursin("reads an empty plan_buffer", sprint(showerror, err))
    # A buffer that starts with a plan for another arm than the planner's.
    other_arm = SM.default_robot_arm_model(link_lengths_m=(0.2, 0.16, 0.12))
    other = RobotArmPlanBuffer(other_arm, BASE, zeros(3))
    err = failure(() -> run_simulation(arm_config(controls=(arm_controller(plan_buffer=other),), guidance=(binding(other),))))
    @test err isa ArgumentError && occursin("another arm or base pose", sprint(showerror, err))
    # The planner and the controller act on different spacecraft.
    err = failure(() -> run_simulation(arm_config(controls=(arm_controller(plan_buffer=buffer),),
        guidance=(binding(buffer; spacecraft_index=2),), spacecraft=two)))
    @test err isa ArgumentError && occursin("plans for spacecraft 2", sprint(showerror, err))
    # A buffer nobody reads, and one read by two controllers.
    err = failure(() -> run_simulation(arm_config(guidance=(binding(RobotArmPlanBuffer()),))))
    @test err isa ArgumentError && occursin("0 models consume one RobotArmPlanBuffer", sprint(showerror, err))
    err = failure(() -> run_simulation(arm_config(controls=(arm_controller(plan_buffer=buffer),
        arm_controller(plan_buffer=buffer, spacecraft_index=2)), guidance=(binding(buffer),), spacecraft=two)))
    @test err isa ArgumentError && occursin("2 models consume one RobotArmPlanBuffer", sprint(showerror, err))
    # A controller with both a plan and a buffer.
    err = failure(() -> run_simulation(arm_config(controls=(arm_controller(plan=direct_plan(), plan_buffer=buffer),),
        guidance=(binding(buffer),))))
    @test err isa ArgumentError && occursin("both a plan and a plan_buffer", sprint(showerror, err))
    # An MPC controller built for another arm.
    two_joint = SM.default_robot_arm_model(link_lengths_m=(0.2, 0.2), link_radii_m=(0.02, 0.02), link_masses_kg=(0.1, 0.1),
        joint_axes=((0.0, 0.0, 1.0), (0.0, 1.0, 0.0)))
    err = failure(() -> run_simulation(arm_config(controls=(arm_controller(plan_buffer=buffer,
        controller=SM.init_robot_arm_joint_mpc(two_joint; dt_s=DT, horizon=4)),), guidance=(binding(buffer),))))
    @test err isa ArgumentError && occursin("built for 2 joints", sprint(showerror, err))
    # Settings of another kind.
    err = failure(() -> run_simulation(arm_config(controls=(arm_controller(plan_buffer=buffer),),
        guidance=(binding(buffer; settings=(target=TARGET,)),))))
    @test err isa ArgumentError && occursin("JointReferenceSettings", sprint(showerror, err))
    # Refused when the run starts, with the binding's error: a link that is not a
    # plan buffer. A request needs the arm's measured joints in the state.
    err = failure(() -> run_simulation(arm_config(guidance=(binding(Ref(0)),))))
    @test err isa AlgorithmGuidanceError && err.reason === :initialization_exception && err.cause isa ArgumentError
    request_error = failure(() -> GA.family_request(JointReferenceFamily(), settings(), (sc=[(pos=zeros(3),)],), nothing,
        0.0, 1, :forced, 1))
    @test request_error isa ArgumentError && occursin("carries no arm joint state", sprint(showerror, request_error))
end

@testset "The arm controller's MPC from the arm model" begin
    from_plan = SM.init_robot_arm_joint_mpc(direct_plan(); dt_s=DT, horizon=4)
    from_model = SM.init_robot_arm_joint_mpc(ARM; dt_s=DT, horizon=4)
    @test C._same_state(from_plan, from_model)
end

end # module JointReferenceFamilyTests
