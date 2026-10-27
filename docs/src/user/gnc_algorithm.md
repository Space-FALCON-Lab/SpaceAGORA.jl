# Writing a GNC Algorithm

Which guide you need depends on what you are adding. A controller or estimator
updates a value at each tick. A planner proposes a product, and an existing
family checks that product before anything uses it. Defining the rules of a new
family, or a service boundary for a companion package, is a change to
SpaceAGORA itself.

| What you want to add | Start here |
|---|---|
| A controller that holds a command, or an estimator that holds an estimate | [Feedback controllers and estimators](gnc_feedback.md) |
| A planner for an existing family: RPO, burn plans, joint references or EDG | [Planning algorithms](gnc_planning.md) |
| A new family, its trusted validation, or a companion service boundary | [Contributor kit integrations](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/docs/contributing/contributor_kit_integrations.md) |

To keep your work in a repository of your own, see
[Building Your Own Research Repository](research_repository.md): which starter
to choose, how to pin SpaceAGORA and how to run the shared checks. Read
[Your First Script](first_script.md) and the
[Extension Conventions](../reference/extension_conventions.md) before either
algorithm guide. The [Extension Hook Reference](../generated/extension_hook_reference.md)
lists every hook's signature and default.

## Which profile?

Use the [feedback profile](gnc_feedback.md) when each tick computes a command or
an estimate from the current state. Use the [planning profile](gnc_planning.md)
when the algorithm proposes something that can be checked before it is used: a
trajectory, a burn plan, a joint path or a guidance decision. One design can use
both. In the RPO pilot a planner proposes references and an MPC controller
tracks them, and the two share a link.

The two profiles use the same terms:

| Step | Feedback profile | Planning profile |
|---|---|---|
| Input | the state at the tick, or a measurement | an owned request snapshot built by trusted code |
| Output | a command or an estimate | a product proposed in a result |
| Check | your own limits, applied in the tick | validation by the family, which can reject the product |
| Keep | held until the next tick | installed into a link and used while it is current |

## The feedback profile

[Feedback controllers and estimators](gnc_feedback.md) explains the tick, the
held command and the output hooks, then works through a rate damper and an
estimator paired with a controller.

### Declarations

[Declarations](gnc_feedback.md#Declarations) covers the spacecraft a model acts
on, the tick period it assumes, what happens when a spacecraft becomes inactive,
reading the held command, and models that bring their own callbacks.

### Worked example A: a rate damper

[Build and check the rate damper](gnc_feedback.md#Worked-example-A:-a-rate-damper).
The page shows the same script that CI runs.

### The lifecycle hooks

[The lifecycle hooks](gnc_feedback.md#The-lifecycle-hooks) covers the checks
before a run, the reset at the start of each run, atmosphere entry and exit, and
bringing a linked producer up to date before a control tick.

### Worked example B: an estimator and a link

[Pair an estimator with a controller](gnc_feedback.md#Worked-example-B:-an-estimator-and-a-link)
through a shared estimate, with an exact gyro.

#### Pairing through a link

[Pairing through a link](gnc_feedback.md#Pairing-through-a-link) explains
producers and consumers, the checks the engine runs on a link before the run,
and how a link is reset at the start of each run.

## The planning profile

[Planning algorithms](gnc_planning.md) explains owned requests, proposals,
capability checks, validation and installation. Its table of families shows
which existing problem each planner can serve.

### Worked example C: an RPO planner

[Write an RPO planner](gnc_planning.md#Worked-example-C:-an-RPO-planner) that
proposes a reference, run it in the RPO pilot and inspect the accepted result.

### Other families

[Other families](gnc_planning.md#Other-families) shows how `AlgorithmGuidance`
runs a planner through its family's own request and installation rules.

#### Worked example D: a burn planner

[Write a burn planner](gnc_planning.md#Worked-example-D:-a-burn-planner) that
proposes a timed burn for the burn-plan family to check.

#### Worked example E: a joint-reference planner

[Write a joint-reference planner](gnc_planning.md#Worked-example-E:-a-joint-reference-planner)
that proposes a joint path for the arm's controller to track.

#### The EDG guidance family

[The EDG guidance family](gnc_planning.md#The-EDG-guidance-family) describes the
optional companion, the decision it proposes and the binding that runs it.

### The family side

The functions below belong to trusted code in SpaceAGORA. An algorithm for an
existing family does not implement them; it implements the
[planning hooks](gnc_planning.md#The-planning-profile) instead. The
[GNC Algorithm Contract API](../generated/gnc_algorithm_contract_api.md) gives
their signatures.

| Part | Functions |
|---|---|
| The family | `validate_capabilities` and `validate_proposal` check capabilities and proposals; `install!` installs an accepted product; `is_current` says whether an installed product may still be used; `on_rejection` and `fallback_product` decide what a refusal does; `family_contract_version` and `algorithm_result_type` declare the contract version and the result type; `is_revision_request` marks a request that revises the installed product; `proposed_product` reads the product from a result; `installed_record` and `rejected_record` write the run report's entries |
| The family's simulation integration, for `AlgorithmGuidance` | `family_request` builds the owned request, `family_link` resolves the install target, and `preflight_family` checks the configuration before a run |
| The binding | `algorithm_runtime`, `bound_algorithm`, `bound_family` and `bound_link` return its run state, algorithm, family and install target; `build_algorithm_request` builds an owned request; `trusted_queries` supplies queries that do not come from the algorithm; `fail_algorithm_request!` records a failed request and throws a `GNCAlgorithmError` |

A new family follows
[Adding a GNC algorithm family](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/docs/contributing/gnc_algorithm_families.md),
and a new family or companion service boundary goes through the
[Contributor kit integrations checklist](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/docs/contributing/contributor_kit_integrations.md).
Both start with a pull request to SpaceAGORA. An algorithm for an existing
family needs no such pull request and stays in your own repository.

## House rules for both profiles

- Draw random numbers only from the `rng` you are given, never from Julia's
  global stream. Planning functions take an explicit stream and have no default
  one, and the conformance checks rerun your model with the global stream
  reseeded.
- Keep run state in the model (feedback) or in the state `initialize_algorithm`
  returns (planning), never in globals. A model resets its run state in its
  `initialize_*!` hook.
- Keep a change a user makes between two runs. Some values are both user
  configuration and written by the run, such as a held command or a plan
  buffer's plan. Record a copy of each such value the model writes. In
  `initialize_*!`, reset an entry that still equals that copy to its value at
  the previous run's start. An entry that differs from it was changed by the
  user, so keep it as the new configuration. The built-in models follow this
  rule, and the burn scheduler applies it to its schedule. A robot-arm or RPO
  controller given a held command and no plan, for example, holds that command
  for the whole of every run.
- Survive `deepcopy`. `run_simulation` runs a copy of the configuration unless
  `isolate_state=false`, and a link that two models share is still shared in
  the copy.
- Give each spacecraft's copy fresh run state in `bind_spacecraft`.
- Keep the output hooks pure: no control law, no warnings, no I/O and no large
  allocations. They run at every evaluation of the equations of motion.

## Test your algorithm

Run SpaceAGORA's conformance checks on a configuration that exercises your
model. The [feedback guide](gnc_feedback.md#Test-your-algorithm) shows the
controller and estimator checks, and the
[planning guide](gnc_planning.md#Test-your-algorithm) shows each planner's
requests and simulation. Start Julia with two or more threads so that the
checks can compare concurrent runs with sequential ones. A pass shows that the
interface behaves as the contract requires; validating the science is still
your job.

## Templates

- `templates/control_hook_template.jl`: a feedback controller with its
  declarations, held command and checkpoint hooks.
- `templates/gnc_planner_template.jl`: an RPO planner written with the public
  names, with notes on writing a planner for another existing family.

## Deprecated paths

SpaceAGORA 0.5.0 removed these:

| Removed | Use instead |
|---|---|
| `RPOGuidanceModel` | `make_rpo_configuration(planner=HYPRRPOPlanner(...))`; see [Moving from RPOGuidanceModel](rpo_planner_pilot.md#Moving-from-RPOGuidanceModel) |
| selecting an RPO planner by a `Symbol`, such as `:rrt_star` | a planner object, such as `RRTStarPlanner()` |
| the old hook names (`calcControlEffect!`, `calcControlForceTorque`, `calcControlMassFlowRate`, `calcGuidanceEffect!`, `calcNavigationEffect!` and the rest) | the current hooks; see [Migrating Extension Hooks](../reference/migrating_extension_hooks.md) |
