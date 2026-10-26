# Writing a GNC Algorithm

Choose the route that matches what you are adding. A controller or estimator
updates a value at each tick. A planner proposes a product that an existing
family checks before it is used. Defining those family rules, or adding a
service boundary for a companion package, is work in SpaceAGORA itself.

| What you want to add | Start here |
|---|---|
| A controller that holds a command, or an estimator that holds an estimate | [Feedback controllers and estimators](gnc_feedback.md) |
| A planner for an existing family: RPO, burn plans, joint references or EDG | [Planning algorithms](gnc_planning.md) |
| A new family, its trusted validation, or a companion service boundary | [Contributor kit integrations](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/docs/contributing/contributor_kit_integrations.md) |

To keep your work in a separate repository, use
[Building Your Own Research Repository](research_repository.md). It explains
which starter to choose, how to pin SpaceAGORA and how to run the shared checks.
Read [Your First Script](first_script.md) and the
[Extension Conventions](../reference/extension_conventions.md) before either
algorithm guide. The [Extension Hook Reference](../generated/extension_hook_reference.md)
lists the signatures and defaults.

## Which profile?

Use the [feedback profile](gnc_feedback.md) when a tick computes a command or
estimate from the current state. Use the [planning profile](gnc_planning.md)
when an algorithm proposes a trajectory, burn plan, joint path or decision
that can be validated before use. A planner and controller can use both,
sharing a link as the RPO pilot does.

The two profiles share one vocabulary:

| Step | Feedback profile | Planning profile |
|---|---|---|
| Input | state at the tick, or a measurement | an owned request snapshot built by trusted code |
| Output | a command or an estimate | a product proposed in a result |
| Check | your own limits, in the tick | validation by the family, which can reject |
| Keep | hold until the next tick | install into a link and use while current |

## The feedback profile

[Feedback controllers and estimators](gnc_feedback.md) explains the tick,
the held command and the output hooks, then walks through a rate damper and
an estimator paired with a controller.

### Declarations

[Declarations](gnc_feedback.md#Declarations) covers spacecraft scope, assumed
period, inactive-spacecraft policy, held-command access and custom callbacks.

### Worked example A: a rate damper

[Build and check the rate damper](gnc_feedback.md#Worked-example-A:-a-rate-damper),
using the same script that the example checks run.

### The lifecycle hooks

[The lifecycle hooks](gnc_feedback.md#The-lifecycle-hooks) explains pre-run
checks, reset at the start of a run, atmosphere events and preparing a linked
producer before control.

### Worked example B: an estimator and a link

[Pair an estimator with a controller](gnc_feedback.md#Worked-example-B:-an-estimator-and-a-link)
using a shared estimate and an exact gyro.

#### Pairing through a link

[Pairing through a link](gnc_feedback.md#Pairing-through-a-link) explains
producers, consumers, their configuration checks and link initialization.

## The planning profile

[Planning algorithms](gnc_planning.md) explains owned requests, proposals,
capability checks, validation and installation. Its family table shows which
existing problem each planner can serve.

### Worked example C: an RPO planner

[Write an RPO planner](gnc_planning.md#Worked-example-C:-an-RPO-planner) that
proposes a reference, run it in the pilot and inspect the accepted result.

### Other families

[Other families](gnc_planning.md#Other-families) shows how `AlgorithmGuidance`
runs a planner through its family's existing request and installation rules.

#### Worked example D: a burn planner

[Write a burn planner](gnc_planning.md#Worked-example-D:-a-burn-planner) that
proposes a timed burn checked by the burn-plan family.

#### Worked example E: a joint-reference planner

[Write a joint-reference planner](gnc_planning.md#Worked-example-E:-a-joint-reference-planner)
that proposes a joint path for the arm's controller to track.

#### The EDG guidance family

[The EDG guidance family](gnc_planning.md#The-EDG-guidance-family) describes the
optional companion, the decision it proposes and the binding that runs it.

### The family side

These functions belong to trusted core code. An algorithm for an existing
family implements the [planning hooks](gnc_planning.md#The-planning-profile)
instead. The [GNC Algorithm Contract API](../generated/gnc_algorithm_contract_api.md)
gives the signatures of the core hooks:

- **The family** checks capabilities and proposals with `validate_capabilities`
  and `validate_proposal`, installs accepted products with `install!`, checks
  `is_current`, and selects `on_rejection` and `fallback_product`. It declares
  `family_contract_version` and `algorithm_result_type`, identifies revisions
  with `is_revision_request`, reads `proposed_product`, and produces
  `installed_record` and `rejected_record` entries for the run report.
- **The family's simulation integration** supplies the owned `family_request`,
  resolves its install target through `family_link`, and checks the configuration
  with `preflight_family` before a run.
- **The binding** supplies its `algorithm_runtime`, `bound_algorithm`,
  `bound_family` and `bound_link`, builds an owned request with
  `build_algorithm_request`, provides independent `trusted_queries`, and handles
  a stopped request with `fail_algorithm_request!`.

For a new family, follow [Adding a GNC algorithm family](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/docs/contributing/gnc_algorithm_families.md).
For a new family or companion service boundary, use the
[Contributor kit integrations checklist](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/docs/contributing/contributor_kit_integrations.md).
Both start with a core pull request; the existing-family algorithm stays in
your own repository.

## House rules for both profiles

- **Draw randomness only from the `rng` you are given**, never from Julia's
  global stream: planning functions take an explicit stream and have no default
  one, and the conformance checks rerun your model with the global stream
  reseeded.
- **Keep run state in the model** (feedback) or in the state
  `initialize_algorithm` returns (planning), never in globals, and reset a
  model's run state in its `initialize_*!` hook.
- **Keep a user's change between two runs.** A value a user may configure that
  the run also writes, such as a held command or a plan buffer's plan, goes back
  at the next run's start only while it still equals the model's last write:
  record a copy of each write in the model, and in `initialize_*!` return an
  entry still equal to it to its value at the previous run's start, and keep an
  entry changed since as the user's new configuration. Every built-in model
  follows this rule, the burn scheduler's schedule included, so a held command
  configured with no plan is the command for the whole run, in every run.
- **Survive `deepcopy`.** `run_simulation` runs a copy of the configuration unless
  `isolate_state=false`, and a link shared by two models stays shared in the copy.
- **`bind_spacecraft` gives fresh run state** to each spacecraft's copy.
- **Keep the output hooks pure:** no law, no warnings, no I/O and no heavy
  allocation. They run at every evaluation of the equations of motion.


## Test your algorithm

Run the engine's conformance checks on a representative configuration.
The [feedback guide](gnc_feedback.md#Test-your-algorithm) shows the controller
and estimator checks; the [planning guide](gnc_planning.md#Test-your-algorithm)
shows each planner's requests and simulation fixture. Run with two or more
threads to compare concurrent and sequential runs. A pass checks interface
behaviour; scientific validation remains yours.

## Templates

- `templates/control_hook_template.jl`: a feedback controller, including its
  declarations, held command and checkpoint hooks.
- `templates/gnc_planner_template.jl`: an RPO planner using the public names,
  with the steps for another existing family.

## Deprecated paths

- **`RPOGuidanceModel`** was removed in SpaceAGORA 0.5.0.
  Use `make_rpo_configuration(planner=HYPRRPOPlanner(...))`; see
  [Moving from RPOGuidanceModel](rpo_planner_pilot.md#Moving-from-RPOGuidanceModel).
- **Selecting an RPO planner by a `Symbol`**, such as `:rrt_star`, was removed in
  SpaceAGORA 0.5.0. Pass a planner object, such as `RRTStarPlanner()`.
- **The old hook names** (`calcControlEffect!`, `calcControlForceTorque`,
  `calcControlMassFlowRate`, `calcGuidanceEffect!`, `calcNavigationEffect!` and
  the rest) were removed in SpaceAGORA 0.5.0; see
  [Migrating Extension Hooks](../reference/migrating_extension_hooks.md).
