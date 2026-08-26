# Capability-Aware SpaceAGORA Integration Proposal

Status: proposal; no implementation is implied by this document

Source outline: *Capability-Aware Learning Architecture for Distributed Small Spacecraft Autonomy*

Target packages: `SpaceAGORA.jl` and `SpaceAGORA_RL.jl`

Prepared: August 26, 2026

## 1. Executive recommendation

Implement the outline as a new Earth-observation (EO) task family in
`SpaceAGORA_RL.jl`, backed by typed, public interfaces in `SpaceAGORA.jl`.
Do not place the learned policy inside the dynamics or GNC layers.

The package boundary should be:

```text
SpaceAGORA physical truth and bounded execution
    -> telemetry/resource snapshot
    -> capability estimator
    -> capability-conditioned EO MDP and fleet allocator in SpaceAGORA_RL
    -> typed mission activity command
    -> independent SpaceAGORA feasibility/safety gate
    -> conventional GNC/control execution
```

The first useful implementation should be a deterministic, event-driven EO
surrogate with an oracle capability estimate, rule and optimization baselines,
and one capability-conditioned policy. It should not begin with meta-RL or a
full battery/camera radiation model. Those later features only become useful
after the task, action, capability, and metric contracts are stable.

## 2. Assumptions and scope decisions

1. The requested `SpaceAGORA_rl` destination is the existing nested package
   `SpaceAGORA.jl/SpaceAGORA_RL.jl`.
2. The source outline is a research architecture, not a request to implement
   all physical subsystem models now.
3. EO is a new task family alongside the existing `aerobraking` and `rpo`
   task families; it does not replace them.
4. SpaceAGORA continues to own physical propagation, vehicle state, GNC,
   actuator authority, and simulation execution.
5. SpaceAGORA_RL owns learning and research-policy concerns. It may contain
   inexpensive degradation surrogates, but a physical model used by general
   SpaceAGORA simulations should ultimately be promoted to SpaceAGORA.
6. All cross-package calls should converge on documented root-level
   SpaceAGORA APIs. The EO implementation should not add new imports of private
   `SpaceAGORA.SimulationModel.*` internals.
7. A coordinated fleet is one coupled decision environment. It must not use
   `run_constellation_ensemble` as if each satellite were an independent RL
   worker; that API explicitly assumes uncoupled members. Parallelism should
   instead run complete fleet episodes across workers.
8. The learned layer proposes bounded mission activities. An independent,
   deterministic gate retains final command authority.

## 3. Current repository fit

The existing RL package already provides the main reusable seams:

- `src/core/mdp_interface.jl`: task marker plus observation normalization;
- `src/core/backend_interface.jl`: reset/step backend contract;
- `src/core/action_mask_interface.jl`: infeasible-action masking;
- `src/core/transition_schema.jl`: learner-facing transitions and summaries;
- `src/algorithms`: DDQN and A2C learners, checkpoints, and policies;
- `src/runtime`: TOML resolution, manifests, sessions, and histories;
- `src/parallel`: rollout workers and central learner orchestration;
- `src/tasks/aerobraking`: precedent for surrogate and live SpaceAGORA
  backends, baselines, metrics, and generalization campaigns; and
- `src/tasks/rpo`: precedent for a task-specific SpaceAGORA adapter.

The main package already provides multi-spacecraft propagation, vehicle and
reaction-wheel state, typed propulsive command boundaries, mission-operation
ownership under `src/mission/operations`, simulation campaigns, runtime
telemetry, and public extension rules. It does not yet provide a general EO
task scheduler, battery state-of-health model, camera quality model, data
handling model, capability contract, or crosslink-aware fleet allocator.

That gap suggests adding a task package and a narrow interface, not modifying
the existing aerobraking MDP.

## 4. Proposed packaged file structure

Legend: `[new]` is proposed, `[extend]` is an existing file to change, and
`[later]` is intentionally deferred beyond the first integration slice.

```text
SpaceAGORA.jl/
├── Project.toml                                      [existing]
├── src/
│   ├── SpaceAGORA.jl                                 [extend: public API]
│   ├── core/simulation_model.jl                      [extend: aggregator]
│   ├── mission/operations/
│   │   ├── capabilities/
│   │   │   ├── capabilities.jl                      [new: module owner]
│   │   │   ├── capability_snapshot.jl               [new: typed estimate]
│   │   │   ├── activity_commands.jl                 [new: bounded commands]
│   │   │   └── feasibility.jl                       [new: deterministic gate]
│   │   └── earth_observation/
│   │       ├── earth_observation.jl                  [new: module owner]
│   │       ├── task_types.jl                         [new: target/window/value]
│   │       └── product_requirements.jl               [new: quality contract]
│   ├── simulation/engine/
│   │   └── mission_snapshot.jl                       [new: public snapshot API]
│   └── vehicle/
│       ├── power/                                    [later]
│       │   ├── power_models.jl
│       │   └── battery_degradation.jl
│       ├── payloads/earth_observation_camera.jl      [later]
│       └── data_handling/data_models.jl              [later]
├── examples/
│   └── Earth_EO_Capability_Constellation.jl          [new after live adapter]
└── test/
    ├── unit/capability_contract_tests.jl             [new]
    ├── unit/earth_observation_contract_tests.jl      [new]
    └── integration/eo_rl_command_boundary_tests.jl   [new]

SpaceAGORA.jl/SpaceAGORA_RL.jl/
├── Project.toml                                      [extend: SpaceAGORA dep]
├── src/
│   ├── SpaceAGORA_RL.jl                              [extend: includes/exports]
│   ├── core/                                         [existing contracts reused]
│   ├── algorithms/
│   │   └── meta_rl/                                  [later]
│   │       ├── context_encoder.jl
│   │       └── adaptation.jl
│   └── tasks/earth_observation/
│       ├── EarthObservationMDP.jl                    [new: task aggregator]
│       ├── types.jl                                  [new: state/config/result]
│       ├── mdp/
│       │   ├── actions.jl                            [new]
│       │   ├── observations.jl                       [new]
│       │   ├── normalization.jl                      [new]
│       │   ├── action_masks.jl                       [new]
│       │   ├── rewards.jl                            [new]
│       │   └── termination.jl                        [new]
│       ├── capability/
│       │   ├── degradation_scenarios.jl              [new]
│       │   ├── estimator_interface.jl                [new]
│       │   └── estimators.jl                         [new]
│       ├── backend/
│       │   ├── event_surrogate.jl                    [new: first backend]
│       │   ├── spaceagora_adapter.jl                 [new: live backend]
│       │   ├── access_cache.jl                       [new]
│       │   └── fleet_environment.jl                  [new]
│       ├── coordination/
│       │   ├── allocator_interface.jl                [new]
│       │   ├── centralized_allocator.jl              [new]
│       │   ├── decentralized_messages.jl             [new]
│       │   └── crosslink_model.jl                    [new]
│       ├── baselines/
│       │   ├── fixed_schedule.jl                     [new]
│       │   ├── rule_based.jl                         [new]
│       │   └── optimization_oracle.jl                [new]
│       └── evaluation/
│           ├── metrics.jl                            [new]
│           ├── explanation_records.jl                [new]
│           ├── generalization_suites.jl              [new]
│           └── reports.jl                            [new]
├── configs/earth_observation/
│   ├── smoke.toml                                    [new]
│   ├── single_vehicle_surrogate.toml                 [new]
│   ├── fleet_surrogate.toml                          [new]
│   ├── spaceagora_physics.toml                       [new]
│   └── generalization.toml                           [new]
├── scripts/earth_observation/
│   ├── train.jl                                      [new]
│   ├── evaluate.jl                                   [new]
│   ├── run_baselines.jl                              [new]
│   └── run_generalization_suite.jl                   [new]
├── examples/earth_observation/
│   └── capability_aware_constellation.jl             [new]
└── tests/earth_observation/
    ├── contract_tests.jl                             [new]
    ├── mdp_tests.jl                                  [new]
    ├── surrogate_backend_tests.jl                    [new]
    ├── spaceagora_adapter_tests.jl                   [new]
    ├── coordination_tests.jl                         [new]
    └── campaign_tests.jl                             [new]
```

This is the target layout, not the first pull request. The smallest first
slice creates only the SpaceAGORA contracts, the EO task's `types`, `mdp`,
`event_surrogate`, three baselines, metrics, one smoke configuration, and
tests. Empty directories and placeholder implementations should not be added.

## 5. Ownership and dependency rules

| Concern | Canonical owner | Reason |
|---|---|---|
| Orbit, attitude, environment, actuator state | SpaceAGORA | Physical simulation truth |
| Battery/camera/data physical models | SpaceAGORA when implemented | Reusable beyond RL |
| Telemetry/resource snapshot | SpaceAGORA | Stable simulator-to-autonomy boundary |
| EO task and product requirement records | SpaceAGORA | Rules and optimizers need them too |
| Mission activity command and safety disposition | SpaceAGORA | Keeps authority independent of learning |
| Degradation campaign distributions | SpaceAGORA_RL | Research randomization concern |
| Oracle/noisy capability estimator variants | SpaceAGORA_RL initially | Experiment concern; promote mature flight estimator interfaces later |
| EO observation, reward, action mask, termination | SpaceAGORA_RL | MDP definition |
| DDQN/A2C/meta-learning | SpaceAGORA_RL | Learning concern |
| Fleet allocation and crosslink experiments | SpaceAGORA_RL | Decision-policy concern |
| Run manifests, checkpoints, reports | SpaceAGORA_RL | Existing experiment infrastructure |
| Final GNC/actuator execution | SpaceAGORA | Assurance boundary |

SpaceAGORA must not depend on SpaceAGORA_RL. Once the public contracts exist,
SpaceAGORA_RL should declare the parent package in its `Project.toml`:

```toml
[deps]
SpaceAGORA = "afbfb69f-5c0b-4832-b760-43725dff8540"

[sources]
SpaceAGORA = {path = ".."}

[compat]
SpaceAGORA = "0.1"
```

This replaces path/`LOAD_PATH` discovery for the new EO integration. Existing
aerobraking and RPO adapters can migrate separately; they are not prerequisites
for the EO work.

## 6. Minimum interface contracts

The names below are proposals. Field names must retain units.

### 6.1 SpaceAGORA-owned records

`CapabilitySnapshot`

- spacecraft identifier and timestamp;
- mean capability fractions for usable energy, generation, agility, payload,
  compute, storage, and downlink;
- uncertainty for each capability value;
- estimator/source identifier and freshness; and
- validity flags. Values should not silently imply certainty.

`EarthObservationTask`

- task identifier, target geometry, open/close time, mission value;
- allowed collection modes;
- minimum product quality;
- expected raw data volume and processing demand; and
- optional latency/revisit terms.

`MissionActivityCommand`

- spacecraft and task identifiers;
- activity kind (`collect`, `process`, `downlink`, `recharge`, `desaturate`,
  `defer`, or `reassign`);
- requested mode and bounded start/duration;
- resource reservation; and
- policy/explanation correlation identifier.

`CommandDisposition`

- accepted/rejected status;
- deterministic reason code;
- limiting constraint and margin; and
- optionally clipped/bounded command details.

### 6.2 SpaceAGORA root functions

```julia
mission_snapshot(simulation, spacecraft_id, epoch_s)
evaluate_activity(snapshot, task, command)
execute_activity!(simulation, accepted_command)
```

The exact implementation can use existing engine objects internally, but the
RL package should call only public root functions. `evaluate_activity` must be
deterministic and must run for learned, rule-based, optimization, and operator
commands alike.

### 6.3 SpaceAGORA_RL task records

`EODecisionState` should contain mission time, pending candidate tasks,
resource state, capability estimate, estimator belief/uncertainty, compact
fleet summaries, and communication freshness.

`EOAction` should use a stable activity kind plus a candidate slot, rather than
embedding an unbounded task identifier in the neural-network output. At each
decision, the environment sorts and selects the top `K` candidates, creates a
fixed action catalog, and masks infeasible entries through the existing
`valid_action_mask` interface.

`EOStepResult` should parallel the existing aerobraking step result: next
state, resolved action, raw and normalized observations, reward, termination
flags, command disposition, resource deltas, completed products, and an
explanation record.

### 6.4 Explanation record

Every decision path, including baselines, should emit the same record:

```text
selected action
capability estimate and uncertainty used
active constraints and margins
predicted mission value/resource outcome
material rejected alternatives
deterministic command disposition
policy/checkpoint/config/seed identifiers
```

This makes explanation fidelity testable and supports replay without requiring
neural-network weight interpretation.

## 7. MDP and fleet formulation

The first observation vector should be deliberately bounded. It should include:

- normalized mission time and next contact/window timing;
- SOC, thermal margin, free storage, compute backlog, and wheel margin;
- the seven capability estimates plus uncertainty/freshness;
- fixed features for the top `K` candidate tasks; and
- compact fleet availability and message-age summaries.

The first reward should be decomposed and logged, not returned as one opaque
number:

```text
delivered task value
- resource use
- risk/constraint-margin penalty
- avoidable degradation penalty
- uncertainty penalty
- reassignment/communication cost
```

Hard safety limits belong in action masks and the SpaceAGORA command gate, not
only in reward penalties.

The initial fleet design should use centralized training with a fleet-level
environment. Decentralized execution can be evaluated by restricting each
local policy to its own state plus timestamped neighbor summaries. Crosslink
degradation is then modeled by message delay, loss, bandwidth, and stale
capability/task summaries.

## 8. Integration phases and verification gates

### Phase 0 — Freeze boundaries and fixtures

Deliver:

- typed capability, EO task, activity command, and disposition records;
- public SpaceAGORA snapshot/feasibility stubs with deterministic fixture
  implementations;
- one hand-authored two-spacecraft, six-task scenario; and
- a versioned observation/action schema.

Verify:

- constructors reject invalid units/ranges;
- serialization round-trips;
- the same snapshot and command always produce the same disposition; and
- SpaceAGORA tests pass without loading SpaceAGORA_RL.

Exit criterion: the RL task can be implemented against fixtures without
depending on SpaceAGORA internals.

### Phase 1 — Single-vehicle event surrogate and baselines

Deliver:

- event-driven EO state transitions;
- battery, agility, payload-quality, storage, compute, and downlink resource
  ledgers;
- stochastic degradation onset/severity with an oracle estimator;
- action catalog, masks, observation normalization, rewards, and termination;
- fixed schedule, rule-based, and optimization-oracle baselines; and
- metrics/artifacts for mission value retention, oracle regret, constraint
  margins, and adaptation latency.

Verify:

- seeded episodes are bitwise repeatable on CPU;
- resource conservation and bounds hold at every step;
- infeasible actions are masked and rejected by the gate;
- zero-degradation scenarios match nominal baseline expectations; and
- the oracle never scores below the identical feasible fixed schedule due to
  an implementation artifact.

Exit criterion: one command trains/evaluates a policy and all four baselines on
the same versioned campaign.

### Phase 2 — Live SpaceAGORA adapter

Deliver:

- declared `SpaceAGORA` dependency and root API usage;
- orbit/access/attitude/resource snapshots at mission decision events;
- typed command handoff to the deterministic feasibility/safety gate;
- bounded collect/slew/desaturate/recharge execution; and
- cached exogenous access geometry for training throughput.

Verify:

- surrogate and live backends share the same task/config/result interfaces;
- unit tests compare snapshot fields and action dispositions at known states;
- an accepted action reaches only the typed command path;
- a rejected action cannot mutate physical state; and
- a short single-vehicle live episode is deterministic for a fixed seed and
  solver configuration.

Exit criterion: a frozen policy can run unchanged against both backends.

### Phase 3 — Capability estimation and generalization

Deliver:

- oracle, delayed/noisy, and telemetry-derived estimators;
- battery, wheel, and payload degradation sweeps, individually and jointly;
- configuration-conditioned DDQN/A2C models;
- held-out severity/composition campaigns; and
- uncertainty calibration and explanation-fidelity metrics.

Only after these pass should a recurrent context encoder or meta-RL adaptation
module be added.

Verify:

- true capability is absent from non-oracle policy observations;
- train/validation/test degradation families are disjoint and recorded;
- checkpoint manifests include schema, estimator, degradation, and backend
  versions; and
- generalization is reported against rule, optimization, and oracle baselines.

Exit criterion: configuration conditioning demonstrates repeatable value over
an otherwise identical policy that omits capability inputs.

### Phase 4 — Fleet allocation and degraded crosslinks

Deliver:

- fleet environment and dynamic task reassignment;
- centralized allocator baseline;
- decentralized message schema and crosslink model;
- partial payload loss and role reassignment; and
- fleet recovery and communication-efficiency metrics.

Verify:

- tasks cannot be double-counted after reassignment;
- reservations prevent two spacecraft from claiming the same exclusive task;
- local execution remains safe with stale/missing messages;
- fleet episodes are isolated across rollout workers; and
- `run_constellation_ensemble` is used only for physically uncoupled studies,
  never to split a coordinated decision episode.

Exit criterion: the fleet retains more value than independent local policies
under at least one preregistered partial-degradation regime without increasing
hard violations.

### Phase 5 — Assurance and flatsat preparation

Deliver:

- deterministic replay from config, seed, checkpoint, and event log;
- bounded-compute/decision-latency measurements;
- command-gate fault injection;
- operator-facing rationale/counterfactual/post-execution reports; and
- a hardware-interface shim that preserves the same snapshot/command records.

Verify:

- replay reproduces actions and dispositions;
- safety behavior is independent of policy implementation;
- timeout or malformed policy output resolves to a documented safe action; and
- uncertainty is calibrated rather than merely displayed.

Exit criterion: the software boundary is suitable for embedded/flatsat work;
hardware procurement and flight qualification remain separate projects.

## 9. Experiment matrix

Every campaign should specify these independent axes in TOML and record them in
the run manifest:

| Axis | Initial values |
|---|---|
| Backend | event surrogate, SpaceAGORA live |
| Policy | fixed, rule, optimization, DDQN/A2C, oracle |
| Estimator | oracle, delayed/noisy, telemetry-derived |
| Degradation | none, battery, agility, payload, joint |
| Severity | preregistered continuous bins |
| Fleet size | 1, 2, then larger only after correctness |
| Crosslink | nominal, delayed, lossy, unavailable |
| Configuration split | training, seen-test, unseen-combination test |

Required metrics from the outline are:

- mission value retention;
- regret relative to oracle;
- adaptation latency;
- hard violations and minimum energy/thermal margins;
- capability-estimation error and calibration;
- seen-to-unseen generalization ratio;
- recovered value after reassignment;
- revisit/gap and delivery latency;
- decision time, memory, and communication bytes; and
- explanation fidelity and replay completeness.

All aggregate tables should include seed count and uncertainty intervals. A
mean alone is insufficient for a stochastic degradation campaign.

## 10. Test and CI placement

Use the package's canonical Julia test entry point and include the EO suites
from it. Do not create a third competing test runner.

Fast PR tests:

- contract construction and serialization;
- observation shape/version and normalization bounds;
- action catalog and mask invariants;
- deterministic surrogate transitions;
- reward component arithmetic;
- command rejection/non-mutation;
- task reservation and reassignment invariants; and
- one single-vehicle and one two-vehicle smoke episode.

Nightly or research campaign tests:

- Monte Carlo degradation sweeps;
- live SpaceAGORA episode parity;
- held-out generalization;
- crosslink degradation; and
- performance/memory regression.

The parent-package additions must also satisfy SpaceAGORA's architecture,
canonical-path, naming, coverage, typed-config, docs, and artifact gates.

## 11. Primary risks and mitigations

| Risk | Consequence | Mitigation |
|---|---|---|
| Starting with full physics | Slow iteration and ambiguous failures | Stabilize the event surrogate and contracts first |
| Conflating health with capability | Planner consumes component-specific details | Expose mission-level capability plus uncertainty |
| Safety only in reward | Learned policy can issue invalid commands | Action masks plus independent SpaceAGORA gate |
| Variable task action space | Existing networks require fixed outputs | Top-`K` candidate slots plus stable action catalog |
| Data leakage | Policy receives latent true degradation | Separate truth, telemetry, estimate, and policy records |
| Fleet/worker confusion | Satellites become falsely independent | One complete fleet environment per rollout worker |
| Private SpaceAGORA imports | Brittle coupling | Declared dependency and root API contract |
| Reward hides behavior | High score masks resource/risk tradeoffs | Persist each reward component and constraint margin |
| Meta-RL too early | Complexity without a validated need | Require conditioning baseline and held-out evidence first |
| Schema drift | Old checkpoints silently misinterpret state | Version observation/action/capability schemas in manifests |

## 12. Recommended first implementation slice

The first development increment should stop after the following vertical
slice:

1. Add the four SpaceAGORA boundary records and deterministic fixture gate.
2. Add `EarthObservationMDP` with one spacecraft, six tasks, and event-driven
   energy/storage/downlink transitions.
3. Add only battery-capacity degradation and an oracle capability estimator.
4. Add `single image`, `downlink`, `recharge`, and `defer` actions.
5. Add fixed, rule-based, and optimization-oracle baselines.
6. Train the existing DDQN learner with the capability fraction in the
   observation.
7. Report mission value retention, regret, SOC margin, and adaptation latency.
8. Prove that the same frozen policy can step through a fixture-backed
   SpaceAGORA adapter without private imports.

This slice is deliberately narrower than the outline. It validates the
architecture and produces an end-to-end comparison before adding reaction
wheels, image-quality degradation, meta-learning, or distributed allocation.

## 13. Overall definition of done

The integrated research framework is complete when:

- SpaceAGORA and SpaceAGORA_RL have one-directional, declared dependencies;
- physical truth, telemetry, capability estimate, and policy observation are
  distinct typed records;
- learned and non-learned planners use identical task, command, gate, and
  metric paths;
- a frozen capability-conditioned policy runs on surrogate and live backends;
- single and joint degradation campaigns are repeatable and manifest-backed;
- fleet reassignment works under partial degradation and limited crosslinks;
- all hard command constraints are enforced outside the policy;
- explanations are tied to recorded inputs, alternatives, outcomes, and
  dispositions; and
- the preregistered baselines and metrics from the research outline can be
  regenerated from repository scripts.
