# Planning Algorithms

This guide is for algorithms that propose a product for an existing family to
check and install. A controller or estimator that updates a held value at each
tick belongs in [Feedback Controllers and Estimators](gnc_feedback.md). Both
kinds follow the [shared house rules](gnc_algorithm.md#House-rules-for-both-profiles).
Defining a new trusted family is a separate
[contribution to SpaceAGORA](gnc_algorithm.md#The-family-side).

The examples are complete scripts in `examples/gnc_algorithm_env/`, a separate
project that uses only the public API. CI runs them, and a docs check fails when
a code block on this page differs from the matching region of its script. From
the SpaceAGORA checkout:

```bash
julia --startup-file=no --compiled-modules=existing --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/setup.jl
julia --startup-file=no --compiled-modules=existing --depwarn=error --threads=2 --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/rpo_planner.jl
```

Run a whole script, imports included, before you change any of its excerpts.
Each script loads `SpaceAGORA`, `Test` and the conformance functions it calls
from `SpaceAGORA.Conformance`. The RPO and joint-reference scripts also import
`norm` from `LinearAlgebra`, and the RPO script imports `AbstractRNG` from
`Random`.

## The planning profile

A planning algorithm proposes and its family decides. A family is trusted code
in SpaceAGORA that defines one problem: what a request contains, what product an
algorithm proposes, how a proposal is validated, where an accepted product is
installed, when it stops being current, and what happens on rejection. Your
algorithm never validates or installs its own result.

`run_algorithm_request!` takes every request through the same steps:

1. Trusted code builds an owned snapshot of what the algorithm needs, such as
   the state, the goal, the limits and how long the request stays valid. The
   binding records a copy of it.
2. The family checks the algorithm's declared capabilities against the request.
3. The algorithm proposes from its own copy of the request:
   `propose!(state, algorithm, request, rng)`, or
   `propose_revision!(state, algorithm, request, active_product, rng)` when the
   request asks for a revision of the installed product.
4. The family validates a copy of the proposal against the request it kept.
5. An accepted product is installed into the link the family shares with the
   product's consumer. Nothing else runs between validation and installation.
6. The binding records the outcome.

If the capability check or the validation refuses, the family's `on_rejection`
policy decides what happens. `:fail` stops the run with a `GNCAlgorithmError`
that carries the run report. `:keep_previous` keeps the installed product and
records the refusal. `:install_fallback` installs the family's
`fallback_product`, such as a hold.

An algorithm implements these functions:

| Function | Required | What it does |
|---|---|---|
| `algorithm_family(algorithm)` | yes, usually inherited from the family's abstract type | the family whose rules govern the algorithm |
| `propose!(state, algorithm, request, rng)` | yes | a result for the request: a `:candidate` with a product, or `:infeasible`, `:unsupported` or `:failed` with a termination reason |
| `propose_revision!(state, algorithm, request, active_product, rng)` | no | a revision of the installed product; by default an `:unsupported` result, never a full proposal in its place |
| `algorithm_capabilities(algorithm)` | no | what it supports, in the family's capability type; by default nothing (`default_capabilities(family)`) |
| `initialize_algorithm(algorithm, run_context)` | no | fresh run state for each run; `nothing` by default |
| `require_algorithm_support(algorithm)` | no | throw before the run if something it needs, such as an optional package, is missing |
| `consumes_rng(algorithm)` | no | `false` for an algorithm that draws nothing from `rng`; the default is `true` |

The existing families:

| Family | Request and product | Validation | Installed into | On rejection | Algorithms | Runs in a simulation through |
|---|---|---|---|---|---|---|
| RPO reference, `RPOReferenceFamily` | `RPOPlanningRequest`, `RPOReference` | `validate_rpo_result`: identity, timing, endpoints, speed and acceleration limits, sampled clearance | the `RPOPlanBuffer`, read by the `RPOMPCControlModel` | `:fail`, or `:install_fallback` (a hold) when configured | `DirectRPOPlanner`, `HYPRRPOPlanner` (with `using SpaceAGORA, HYPR`), the seven comparison planners | `make_rpo_configuration(planner=...)` |
| Burn plan, `BurnPlanFamily` | `BurnPlanRequest`, a `PropulsiveBurnPlan` | `validate_burn_plan`: the burn scheduler's own size check (`burn_plan_sizing`), thrust, specific impulse and window | the spacecraft's slot in the burn scheduler's plan buffer, `BurnPlanSlot` | `:keep_previous` | `ApoapsisBurnPlanner` | `AlgorithmGuidance` |
| Joint reference, `JointReferenceFamily` | `JointReferenceRequest`, a `RobotArmPlan` | `validate_joint_reference`: the arm and its base, the sample grid, the start joints, joint limits, the goal and obstacle clearance | the `RobotArmPlanBuffer` the arm's `RobotArmControlEffector` reads | `:fail` | `QuinticJointPlanner` | `AlgorithmGuidance` |
| EDG guidance, `EDGGuidanceFamily` | `EDGGuidanceRequest`, an `EDGGuidanceDecision` | `validate_edg_guidance_decision`: the request's spacecraft, a mode the configuration allows, a consistent targeting plan and switch time, a bracketing count that has not decreased | the run's `AerobrakingEnergyDepletionState`, read by the EDG control model | `:fail` | `EDGGuidanceAlgorithm`, with `using SpaceAGORAEDG` | `AlgorithmGuidance`, in place of the EDG guidance model |

### Worked example C: an RPO planner

The RPO names are the supported way to write an RPO planner: subtype
`AbstractRPOPlanner`, implement `planner_capabilities` and `plan_rpo!`, and
optionally `initialize_planner` and `retime_rpo!`. These forward to the general
names (`algorithm_capabilities`, `propose!`, `initialize_algorithm`,
`propose_revision!`); the
[RPO planner contract](../reference/rpo_planner_contract.md#The-general-algorithm-contract)
maps one set to the other. The script is
`examples/gnc_algorithm_env/rpo_planner.jl`.

```julia
# A reference along the straight line to the goal, in the target's RTN frame: it
# accelerates, cruises and brakes within fractions of the request's limits. It
# ignores obstacles; the family's validation refuses a line that comes too close.
Base.@kwdef struct LinePlanner <: AbstractRPOPlanner
    speed_fraction::Float64 = 0.5            # of the request's speed limit
    acceleration_fraction::Float64 = 0.5     # of the request's acceleration limit
end

SpaceAGORA.planner_capabilities(::LinePlanner) =
    RPOPlannerCapabilities(state_sources=(:truth,), frames=(:target_rtn,))
SpaceAGORA.consumes_rng(::LinePlanner) = false     # it draws nothing from `rng`

function SpaceAGORA.plan_rpo!(state, planner::LinePlanner, request::RPOPlanningRequest, rng::AbstractRNG)
    refuse(why) = RPOPlanningResult(request_id=request.request_id, status=:infeasible, termination=why)
    start, goal = collect(Float64, request.x_rtn[1:3]), collect(Float64, request.goal_rtn_m)
    distance = norm(goal - start)
    v = planner.speed_fraction * something(request.constraints.max_speed_mps, Inf)
    a = planner.acceleration_fraction * something(request.constraints.max_acceleration_mps2, Inf)
    isfinite(v) && isfinite(a) && v > 0 && a > 0 || return refuse(:needs_speed_and_acceleration_limits)
    v = min(v, sqrt(distance * a))           # too short to reach cruise speed: a triangle
    ramp_s = v / a
    duration_s = distance == 0 ? 0.0 : distance / v + ramp_s
    dt = request.reference_dt_s
    steps = max(1, ceil(Int, duration_s / dt))
    request.time_s + steps * dt <= request.valid_until_s || return refuse(:outlives_the_request)
    # Distance along the line and speed at time t from the start.
    function along(t)
        t >= duration_s && return (distance, 0.0)
        t < ramp_s && return (a * t^2 / 2, a * t)
        t <= duration_s - ramp_s && return (v * (t - ramp_s / 2), v)
        left = duration_s - t
        return (distance - a * left^2 / 2, a * left)
    end
    direction = distance == 0 ? zeros(3) : (goal - start) / distance
    times = [k * dt for k in 0:steps]
    profile = along.(times)
    reference = RPOReference(t_ref_s=times,
        r_ref_rtn_m=reduce(hcat, [start + s * direction for (s, _) in profile]),
        v_ref_rtn_mps=reduce(hcat, [speed * direction for (_, speed) in profile]),
        origin_time_s=request.time_s, valid_until_s=request.valid_until_s,
        chaser_id=request.chaser_id, target_id=request.target_id,
        geometry_revision=request.geometry_revision)
    return RPOPlanningResult(request_id=request.request_id, status=:candidate,
        termination=:completed, reference=reference, diagnostics=(duration_s=duration_s,))
end
```

Run it in the RPO pilot and read the run report:

```julia
# The RPO pilot: a chaser 3 m from the station plans to 5 m, and an MPC controller
# tracks the installed reference.
configuration = make_rpo_configuration(planner=LinePlanner(), mission_time_s=2.0)
solution = run_simulation(configuration; return_solution=true)
report = rpo_run_report(solution)
installed = only(filter(record -> record.event === :installed, only(report.planners).records))
println("installed request ", installed.request_id, ": ", installed.validation.reason,
        ", reference of ", last(installed.reference.t_ref_s), " s")
@test installed.validation.accepted
```

A planner that cannot plan returns a result with a status and a termination
reason; it does not throw. Under the RPO family's default `:fail` policy the run
then stops with an `RPOPlanningError` whose reason is that status, and the
failure record in the run report keeps the termination reason.

`LinePlanner` declares that it plans from truth states in the target's RTN frame
and does not retime. The family checks every request against these capabilities
before asking the planner, and refuses a request that needs more.

`consumes_rng = false` declares that the planner draws no random numbers, and
the conformance checks confirm that it leaves the stream untouched. A planner
that does draw uses only `rng`, the stream the binding derives from the run's
seed and the spacecraft IDs.

[Run an RPO case with interchangeable planners](rpo_planner_pilot.md) covers
replanning, holds, retiming and inspection in the pilot.

### Other families

A planner of another family runs through `AlgorithmGuidance`, the generic
binding, which you add among the guidance effectors.
`AlgorithmGuidance(; algorithm, spacecraft_index, settings, link, seed, events,
replan_interval_s)` makes a request at each `AlgorithmEvent`. When
`replan_interval_s` is finite it also makes one at the run's start and then
whenever that many seconds have passed since the last request. Requests are
made at the binding's guidance ticks. `settings` holds the family's request
settings, and `link` is the install target, for a family that does not find its
own in the run. The binding's random stream is derived from `seed`, the family
and the spacecraft's ID.

A configuration may hold only one binding per install target. `AlgorithmGuidance`
refuses checkpoint writing and resume, because its run state is not
checkpointed. On a run made with `return_solution=true`,
`algorithm_run_report(solution)` returns every binding's records: each request,
installation, rejection and failure, in order.

#### Worked example D: a burn planner

A planner of the burn-plan family subtypes `AbstractBurnPlanner` and proposes a
`PropulsiveBurnPlan`: a constant-thrust burn over a window `[start, stop)`. The
family checks it with the burn scheduler's own sizing rule, `burn_plan_sizing`,
so the planner states the propellant and impulse that rule gives and a window
as long as the burn duration it gives. The script is
`examples/gnc_algorithm_env/burn_planner.jl`.

```julia
# A prograde burn of a set Δv that starts a set time after the request.
Base.@kwdef struct BurnSoonPlanner <: AbstractBurnPlanner
    delta_v_mps::Float64 = 5.0      # the burn's Δv, m/s
    lead_s::Float64 = 10.0          # time from the request to the burn's start, s
end
SpaceAGORA.consumes_rng(::BurnSoonPlanner) = false

function SpaceAGORA.propose!(state, planner::BurnSoonPlanner, request::BurnPlanRequest, rng)
    # The scheduler's own size check, which the family's validation repeats.
    sizing = burn_plan_sizing(request.mass_kg, planner.delta_v_mps, 0.0, request.thrust_n, request.isp_s;
                              available_propellant_kg=request.available_propellant_kg)
    sizing.reason === :accepted || return AlgorithmResult(request_id=request.request_id,
        status=:infeasible, termination=sizing.reason)
    start_s = request.time_s + planner.lead_s
    plan = PropulsiveBurnPlan(valid=true, delta_v_mps=planner.delta_v_mps, direction_rad=0.0,
        source_orbit=-1, start_burn_s=start_s, stop_burn_s=start_s + sizing.burn_duration_s,
        thrust_n=request.thrust_n, isp_s=request.isp_s,
        commanded_impulse_n_s=sizing.commanded_impulse_n_s,
        propellant_required_kg=sizing.propellant_required_kg)
    return AlgorithmResult(request_id=request.request_id, status=:candidate,
        termination=:completed, product=plan)
end
```

The family installs plans into the burn scheduler, `BaseThrusterModel`. The
configuration needs exactly one scheduler, with zero `Δv` and no burn window of
its own for the binding's spacecraft, and no guidance model that commands the
scheduler's own burns. The scheduler takes up an installed plan at the end of
the step it was installed in and fires the engine inside the plan's window. The
family takes no `settings` and finds its own install target, so the binding
needs neither.

```julia
# The family installs into the burn scheduler, which plans no burn of its own here
# (zero Δv, no window) and fires a 50 N engine of 300 s specific impulse.
planet = make_no_gram_planet(:earth)
root = Link(root=true, m=100.0, ref_area=2.0)
orbit = InitialCondition(ra=planet.Rp_e + 2000e3, rp=planet.Rp_e + 400e3, i=30.0, ω=0.0, Ω=0.0, ν=10.0)
spacecraft = SpacecraftModel(joints=Joint[], links=[root], root=root, prop_mass=50.0,
                             inertia_tensor=root.inertia, initial_condition=orbit, id=1)
scheduler = BaseThrusterModel(thrust=[50.0], direction=[0.0], Δv=[0.0], start_burn_time=[-1.0],
                              stop_burn_time=[-1.0], Isp=[300.0])
binding = AlgorithmGuidance(algorithm=BurnSoonPlanner(), spacecraft_index=1,
                            events=[AlgorithmEvent(0.0)])      # one request, at the start
config = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime, mission_time=60.0, number_of_orbits=1,
        keplerian=true, orientation_sim=false, num_steps_to_save=100),
    environment_model=EnvironmentModel(
        planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false),
    dynamics_model=DynamicsModel([spacecraft], (InverseSquaredGravityModel(),)),
    guidance_model=GuidanceModel(guidance_effectors=(binding,), guidance_rates=[60.0]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(scheduler,), control_rates=[1.0]),
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances())

solution = run_simulation(config; return_solution=true)
records = only(algorithm_run_report(solution).bindings).records
plan = only(filter(record -> record.event === :installed, records)).product

# By hand: the engine takes the propellant the family validated, and that is the
# rocket equation's mass for the burn's Δv.
used_kg = (100.0 + 50.0) - solution.u[end].sc[1].mass
rocket_kg = (100.0 + 50.0) * (1 - exp(-plan.delta_v_mps / (plan.isp_s * 9.80665)))
println("burn from $(plan.start_burn_s) s to $(round(plan.stop_burn_s, digits=3)) s; propellant used ",
        "$(round(used_kg, digits=6)) kg, planned $(round(plan.propellant_required_kg, digits=6)) kg, ",
        "rocket equation $(round(rocket_kg, digits=6)) kg")
@test used_kg ≈ plan.propellant_required_kg rtol=1e-6
@test plan.propellant_required_kg ≈ rocket_kg rtol=1e-9
```

A refused proposal leaves the installed plan in place (`:keep_previous`), and the
run report records the refusal with the validation's reason, such as
`:burn_in_progress` for a request made during a burn.

#### Worked example E: a joint-reference planner

A planner of the joint-reference family subtypes `AbstractJointReferencePlanner`
and proposes a `RobotArmPlan`: joint angles, rates and accelerations, and the
end-effector position, on the request's time grid. The script is
`examples/gnc_algorithm_env/joint_planner.jl`.

```julia
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
```

The plan buffer, `RobotArmPlanBuffer`, is a link. The binding produces it and
the arm's `RobotArmControlEffector` consumes it, and before the run the engine
checks that exactly one controller reads it, for the binding's spacecraft and
arm. The binding takes the buffer as its `link` and a `JointReferenceSettings`
as its `settings`: the arm, its base pose, the target, the obstacles and the
clearance from them, the goal tolerance, and the reference's time step and
duration. Each request starts from the arm's measured joints.

```julia
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
```

The family's rejection policy is `:fail`: a refused proposal stops the run with
an `AlgorithmGuidanceError` whose reason is the validation's, such as
`:insufficient_clearance` for a path through an obstacle.

#### The EDG guidance family

The energy-depletion guidance decision runs as an algorithm of the EDG guidance
family, `EDGGuidanceAlgorithm()`. It is the `guidance_decision!` of the optional
`SpaceAGORAEDG` companion, reached through the versioned `EDGServices` facade,
so the family's contract version (`family_contract_version`) is the facade's
contract version, 3.0.0. The algorithm decides on a copy of the run's EDG state;
the family validates the decision and installs it in that state, where the EDG
control model reads it. The algorithm draws no random numbers.

`AlgorithmGuidance` runs it in place of `AerobrakingEnergyDepletionGuidanceModel`.
Its settings are `EDGGuidanceSettings(config=..., state=...)`, holding the same
configuration and state as the `AerobrakingEnergyDepletionControlModel`, and it
has no `link`. The binding handles one spacecraft; a run with several uses the
EDG guidance model. The
[EDG Service Contract](../reference/edg_adapter_boundary.md#The-EDG-guidance-family)
shows the configuration, and [Energy Depletion Guidance](edg_guidance.md) covers
installing the companion and the guidance itself.

## Test your algorithm

`SpaceAGORA.Conformance` checks a planner against its family with
`check_algorithm(algorithm, requests; trusted, simulation)`, which returns a
report; with `using Test` loaded, `test_algorithm` records the same checks as
tests. The planner proposes for the family's representative requests in order,
from a fresh state, and the checks are:

| Check | What it requires |
|---|---|
| `family_resolves`, `family_trusted` | `algorithm_family` returns a family, and the family is defined in SpaceAGORA, not in the algorithm's package |
| `result_contract` | every result has the family's result type, a known status and the request's `request_id`, and carries a product exactly when it is a `:candidate` |
| `stream_untouched_if_declared` | an algorithm that declares `consumes_rng = false` leaves the stream unchanged |
| `deterministic_given_stream` | two passes through the requests with equal streams give equal results, with Julia's global random stream seeded differently |
| `request_not_mutated` | the algorithm leaves its copy of each request unchanged |
| `initialize_independent` | two states from `initialize_algorithm` share no mutable object |
| `refuses_explicitly` | a revision request that the capability check refuses gets an `:unsupported` result from `propose_revision!` |
| `driver_round_trip` | the requests run through `run_algorithm_request!`: an accepted product is installed and current, and a refusal follows the family's `on_rejection` |
| `simulation_run` | `simulation(algorithm)`, a short simulation through the family's binding, installs a product and records no failure. Skipped without a `simulation` |

The RPO, burn-plan and joint-reference families each supply representative
requests and a short simulation for these checks. Each simulation function runs
its case twice in sequence and, with two or more threads, twice at the same
time, and fails when the runs differ. The requests include one the family can
refuse: a retiming request for the RPO family, which a planner that does not
retime cannot serve; a request made during a burn for the burn-plan family;
and, for the joint-reference family, a request with an obstacle on the straight
joint-space path, which a planner that ignores obstacles proposes through. Each
example script ends with these checks:

```julia
# The family's representative requests, its trusted validation queries and a short
# simulation through the RPO binding.
@testset "LinePlanner conformance" begin
    test_algorithm(LinePlanner(), rpo_conformance_requests(); trusted=rpo_trusted_queries(),
                   simulation=rpo_conformance_simulation())
end
```

```julia
@testset "BurnSoonPlanner conformance" begin
    test_algorithm(BurnSoonPlanner(), burn_plan_conformance_requests();
                   simulation=burn_plan_conformance_simulation())
end
```

```julia
@testset "SmoothStepJointPlanner conformance" begin
    test_algorithm(SmoothStepJointPlanner(), joint_reference_conformance_requests();
                   simulation=joint_reference_conformance_simulation())
end
```

Run the checks with two or more threads so that the simulation check compares
concurrent runs with sequential ones. A pass certifies interface behavior only.
The method, its parameters and any claims made from it remain your
responsibility.

To keep your algorithm in a repository of your own, tested against a recorded
SpaceAGORA commit, start from a contributor starter; see
[Building Your Own Research Repository](research_repository.md).
