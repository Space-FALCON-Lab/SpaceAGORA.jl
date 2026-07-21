# SpaceAGORA_RL

Julia-native reinforcement-learning utilities for SpaceAGORA aerobraking studies.

The first task namespace implements the paper-compatible Mars Odyssey aerobraking MDP:
the 9-feature observation schema, 13-action apoapsis maneuver table, reward and
termination rules, deterministic paper-mode backend facade, baseline policies,
DDQN and A2C learner infrastructure, CPU/CUDA training-device selection,
CSV/report helpers, and example runners.

Reference material is kept under `Reference Code/` and is not on the runtime path.
Generated artifacts should go under `outputs/`; only placeholders are committed.

## Mars Odyssey aerobraking task

The current production task is the paper-compatible Mars Odyssey aerobraking
problem. Each episode is one aerobraking campaign. Each step is one passage
decision where the policy selects an apoapsis maneuver from the fixed
13-action table in `src/tasks/aerobraking/mdp/actions.jl`.

The normalized observation contains nine paper features:

- drag-passage time,
- apoapsis radius,
- periapsis altitude,
- argument of periapsis,
- right ascension of ascending node,
- inclination,
- date ordinal,
- maximum atmospheric density,
- maximum heat rate.

The action table spans lowering, zero, and raising periapsis maneuvers in
meters per second:

```text
-1.0, -0.5, -0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0
```

The reward model encourages hitting the target final apoapsis while keeping
passage heat rate inside the configured corridor. Termination captures success,
undershoot, impact, out-of-passage periapsis altitude, thermal violations, and
the configured maximum-pass cap. The default paper-surrogate backend is
deterministic when `randomization.nominal = true`; randomized runs perturb the
initial orbital state and optional process noise settings.

## Package layout

- `src/core/` defines the task-neutral MDP, backend, policy, and transition
  interfaces.
- `src/tasks/aerobraking/` contains the Mars Odyssey aerobraking task:
  observations, normalization, action mapping, reward, termination,
  randomization, baseline policies, backend adapters, evaluation metrics, and
  reports.
- `src/algorithms/ddqn/` contains the replay buffer, epsilon schedule, target
  network logic, shared dense network implementation, DDQN learner, checkpoint
  serialization, and greedy policy loader.
- `src/algorithms/a2c/` contains the rollout buffer, discounted-return
  construction, A2C actor-critic learner, and greedy policy loader.
- `src/algorithms/training_device.jl` owns CPU/CUDA device selection and array
  movement. `ext/SpaceAGORA_RLCUDAExt.jl` supplies CUDA.jl methods when CUDA is
  loaded and functional.
- `src/runtime/` loads TOML configuration, writes run manifests, and builds
  training sessions.
- `src/parallel/` owns threaded rollout collection and algorithm-specific
  training loops.

## Configuration and entry points

Commands below assume the current directory is the main `SpaceAGORA.jl` package
directory, so `SpaceAGORA_RL.jl` is a child directory.

The default configuration is `configs/aerobraking/paper_replication.toml`.
Training sessions are normally created through:

```julia
using SpaceAGORA_RL

config = resolve_config()
session = build_training_session(config)
results = train_parallel!(session)
```

The script entry point is:

```text
julia --project=SpaceAGORA_RL.jl SpaceAGORA_RL.jl/scripts/train.jl
```

`training.algorithm` selects `pr_drl`, `ddqn`, or `a2c`. `pr_drl` is the
paper-style parallel randomized DDQN path; `ddqn` remains a compatibility alias
for the same learner without the PR-DRL run label. `training.device` accepts:

- `cpu`: always keep learner arrays on CPU.
- `auto`: use CUDA when CUDA.jl is available and `CUDA.functional()` is true,
  otherwise fall back to CPU.
- `cuda` or `gpu`: require functional CUDA support and throw if it is not
  available.

The available aerobraking configs in `configs/aerobraking/` are:

| Config | Purpose |
| --- | --- |
| `paper_replication.toml` | Short PR-DRL paper-surrogate smoke run, currently 20 global steps on CPU. |
| `pr_drl_paper_surrogate.toml` | Step-limited PR-DRL paper-surrogate setting, currently 1.1M global environment steps. |
| `pr_drl_spaceagora_marsgram.toml` | Step-limited PR-DRL setting using live SpaceAGORA MarsGRAM point-density calls with the lightweight pass transition. |
| `pr_drl_spaceagora_physics.toml` | Step-limited PR-DRL setting using actual SpaceAGORA one-pass trajectory/aero/thermal propagation with the GRAM offline surrogate atmosphere. |
| `ddqn_full.toml` | Legacy step-limited DDQN-named alias for older runs. |
| `ddqn_full_campaigns.toml` | Legacy campaign-limited CUDA experiment; not the paper's 1.1M-step budget. |

The PR-DRL configs use the paper MDP, action table, DDQN hyperparameters,
parallel rollout collection, and randomized initial geometry. Three backend
modes are supported:

- `paper_surrogate`: analytic local atmosphere surrogate.
- `spaceagora_marsgram`: live SpaceAGORA MarsGRAM point calls. This requires a
  host-native `libGRAM` under `SpaceAGORA.jl/data/GRAMSuite.jl/GRAM Suite 2.0/Build/lib`.
- `spaceagora_physics`: actual SpaceAGORA propagation from apoapsis to the next
  apoapsis using the GRAM offline surrogate atmosphere, gravity harmonics, solar
  gravity/SRP, free-molecular aero, Maxwellian heat-rate evaluation, the paper's
  135 km atmosphere boundary, and ±10% aerodynamic coefficient campaign scaling.

The MarsGRAM mode replaces the analytic density sample used for each pass with
live SpaceAGORA MarsGRAM density at the current periapsis latitude/longitude.
It keeps the lightweight per-pass transition model; it is not a full ABTS
continuous-propagation replacement. The SpaceAGORA physics mode removes that
lightweight transition and propagates the pass through SpaceAGORA itself, while
using the prebuilt Mars GRAM surrogate payload instead of live GRAM/SPICE density
queries. AADS still uses bisection to target the 0.15 W/cm^2 heat-rate corridor
center; under `spaceagora_physics`, those predictions also use propagated
SpaceAGORA passes through the same surrogate-density model.

Training has two stopping modes:

- `training.global_steps > 0`: step-limited mode. One global step is one
  environment transition/action decision. The configured `episodes` value is
  not the primary stopping condition.
- `training.global_steps = 0`: campaign-limited mode. `training.episodes`
  becomes the campaign cap. One episode is one complete aerobraking campaign,
  and the total number of action steps depends on how many passes each campaign
  lasts.

Use the progress line to tell which mode is active. A step-limited run prints
`steps=current/target`; a campaign-limited run prints `ep=current/target` and
shows total `steps` separately.

## Common commands

Run the default smoke config:

```text
julia --project=SpaceAGORA_RL.jl SpaceAGORA_RL.jl/scripts/train.jl
```

Run the step-limited PR-DRL paper-surrogate config:

```text
JULIA_NUM_THREADS=16 julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/train.jl \
  SpaceAGORA_RL.jl/configs/aerobraking/pr_drl_paper_surrogate.toml
```

Run the step-limited PR-DRL SpaceAGORA-physics config:

```text
JULIA_NUM_THREADS=16 julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/train.jl \
  SpaceAGORA_RL.jl/configs/aerobraking/pr_drl_spaceagora_physics.toml
```

Run the older step-limited PR-DRL live-MarsGRAM point-call config:

```text
JULIA_NUM_THREADS=16 julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/train.jl \
  SpaceAGORA_RL.jl/configs/aerobraking/pr_drl_spaceagora_marsgram.toml
```

Run the legacy campaign-limited DDQN CUDA config:

```text
JULIA_NUM_THREADS=16 julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/train.jl \
  SpaceAGORA_RL.jl/configs/aerobraking/ddqn_full_campaigns.toml
```

Run the same paper-surrogate scenario with A2C on CUDA without editing a
config file:

```text
JULIA_NUM_THREADS=16 julia --project=SpaceAGORA_RL.jl -e '
using SpaceAGORA_RL
raw = load_config("SpaceAGORA_RL.jl/configs/aerobraking/pr_drl_paper_surrogate.toml")
raw["training"]["algorithm"] = "a2c"
raw["training"]["device"] = "cuda"
config = resolve_config(raw; source_path="SpaceAGORA_RL.jl/configs/aerobraking/a2c_pr_drl_surrogate.toml")
session = build_training_session(config; run_id="a2c_pr_drl_surrogate")
result = train_parallel!(session)
println((device=training_device_name(session.learner.device),
         campaigns=length(result.summaries),
         steps=result.global_step,
         train_steps=session.learner.train_steps,
         loss=session.learner.last_loss,
         output_dir=result.output_dir))
'
```

For shorter A2C checks, copy or override `paper_replication.toml` and set
`training.algorithm = "a2c"`. The `[a2c]` hyperparameter block is already
present in the paper configs.

## Config reference

The major TOML sections are:

- `[scenario]`: phase, backend mode, and whether the scenario is configured for
  training. The paper-style runs use `phase = "Main"` with
  `backend_mode = "paper_surrogate"`, `"spaceagora_marsgram"`, or
  `"spaceagora_physics"`.
- `[reward]`: heat-rate corridor thresholds used by `paper_reward`.
- `[termination]`: campaign termination thresholds such as impact periapsis,
  the paper's 135 km out-of-atmosphere periapsis, maximum passes, and
  thermal-violation behavior.
- `[randomization]`: nominal versus randomized initial conditions and optional
  process-noise settings. Paper-style configs also include ±10% aerodynamic
  coefficient dispersion and MarsGRAM perturbation/seed settings; the full
  physics surrogate path keeps those fields for config compatibility but uses
  the deterministic offline GRAM surrogate grid for density.
- `[ddqn]`: Q-learning hyperparameters, replay buffer size, target-network
  cadence, and hidden-network width.
- `[a2c]`: actor-critic learning rate, discount, segment length, entropy and
  value weights, gradient clipping, and hidden-network width.
- `[epsilon]`: DDQN exploration schedule. A2C samples from its categorical
  actor policy and does not use epsilon-greedy exploration.
- `[training]`: seed, algorithm, device, stopping budget, worker count,
  checkpoint cadence, progress cadence, and run output directory.
- `[reports]`: evaluation/report output directory and artifact toggles.

`n_workers` controls how many rollout campaigns are collected in parallel, but
the actual worker count is capped by `Threads.nthreads()`. Set
`JULIA_NUM_THREADS` before launching Julia when running full-scale training.

## Learning process overview

Training and evaluation share the same environment loop:

1. reset an aerobraking campaign from the configured scenario and RNG seed,
2. build the raw paper observation from the current campaign state,
3. normalize the observation with the scenario normalization bounds,
4. choose one of the 13 maneuver actions,
5. step the aerobraking backend for one passage,
6. compute reward and termination/truncation flags,
7. record the transition and update the episode summary,
8. repeat until success, failure, thermal termination, or pass cap.

The learner only sees normalized observations, integer action indices, scalar
rewards, next normalized observations, and termination/truncation flags. The
domain-specific orbital quantities stay in the aerobraking task layer and in
episode/pass logs.

Rollout collection is parallelized across Julia threads by assigning independent
campaigns to worker tasks. Each worker receives deterministic seeds derived from
the training seed, episode index, and worker id. The backend simulation remains
host-side even when the learner runs on CUDA.

The two learners consume the same task interface differently:

- DDQN is off-policy. It stores transitions in a replay buffer, samples random
  minibatches, and trains a Q-network against bootstrapped targets.
- A2C is on-policy. It keeps live worker campaigns, collects short rollout
  segments from the current actor, computes returns, and updates actor and
  critic directly from that fresh segment.

## Learning packages and dependencies

The DDQN and A2C implementations are package-local Julia code. They do not use
Flux.jl, ReinforcementLearning.jl, or another external learner framework. The
core neural-network, loss, gradient, replay-buffer, rollout-buffer, and Adam
optimizer logic lives under `src/algorithms/`.

The learning stack uses these dependencies:

- `Random`: deterministic seeds, epsilon-greedy exploration, categorical A2C
  sampling, worker RNGs, and network initialization.
- `LinearAlgebra`: orthogonal network initialization through QR factorization.
- `Statistics`: batch losses, advantage scaling, and recent training summaries.
- `CUDA`: optional GPU array backend loaded through `SpaceAGORA_RLCUDAExt`.
- `Distributed` and Julia threads: rollout-worker orchestration and parallel
  campaign collection.
- `Serialization`: learner checkpoints for DDQN and A2C.
- `TOML`, `SHA`, and `Dates`: training configuration, config hashes, run IDs,
  and run manifests.
- `Printf`: progress logging from the training loops.
- `CSV`, `DataFrames`, and `Plots`: evaluation/report artifact generation, not
  policy-gradient or Q-learning math.

## Shared network and optimizer implementation

Both learners use the lightweight `QNetwork` implementation in
`src/algorithms/ddqn/network.jl`: two ReLU hidden layers, an output layer, and
manual forward/backward passes over `Float32` arrays. We keep this network local
instead of depending on a larger neural-network stack so checkpoints,
CPU/GPU array movement, and small test fixtures remain explicit.

Weights are orthogonally initialized. Gradients are clipped by global norm and
updated with the package-local Adam implementation. Checkpoints are written with
Julia serialization and store CPU copies of networks and optimizer state, the
algorithm name, training counters, action table, device name, and run manifest.

## GPU support

GPU support is intentionally scoped to learner tensor operations. Rollout
workers and backend simulation stay on the host; learner parameters, sampled
minibatches, rollout batches, and optimizer state move through the configured
`AbstractTrainingDevice`.

This means CUDA runs still benefit from `JULIA_NUM_THREADS` for environment
rollouts. The GPU is used for the learner forward/backward passes and optimizer
state. A successful explicit CUDA run reports `device = "cuda:0"` in the final
summary printed by the examples above, while the manifest records the requested
device symbol.

The CUDA extension provides:

- `CUDA.cu` conversion for `Float32`, real-valued, and `Bool` arrays.
- `to_cpu_array(::CUDA.CuArray)` and `cpu_scalar(::CUDA.CuArray)` for losses,
  target construction, metrics, checkpoints, and CPU policy snapshots.
- `cuda_functional()` backed by `CUDA.functional()`.

When `training.device = "auto"`, `resolve_training_device` attempts to load
CUDA.jl and returns `CUDATrainingDevice(0)` only when the extension reports a
functional CUDA runtime. Explicit `cuda`/`gpu` selection is stricter and errors
if CUDA cannot be used.

## PR-DRL / DDQN implementation

`training.algorithm = "pr_drl"` uses `DDQNLearner` with the paper's PR-DRL
network and update cadence. The learner maintains an online network, target
network, replay buffer,
epsilon schedule, Adam state, global-step counter, training-step counter, and
last loss. For the SpaceAGORA physics PR-DRL backend, workers keep a full
campaign simulation alive, emit one transition at each apoapsis, and block until
the master ingests that pass-level transition, trains as needed, and sends back
the next action. This matches the paper's master/worker PR-DRL architecture.
The lightweight legacy DDQN path still collects complete worker campaigns before
ingestion. The master trains when:

- the replay buffer has at least `ddqn.batch_size` transitions,
- `global_step >= ddqn.train_start`, and
- `global_step % ddqn.train_frequency == 0`.

Each train step samples a replay minibatch, moves observations to the configured
device, evaluates online and target Q-values for next observations, computes
Double-DQN targets on the host, applies clipped gradients to the online network,
and periodically copies the online network to the target network according to
`ddqn.target_update`.

`GreedyPRDRLPolicy` and `GreedyDDQNPolicy` load the checkpoint `:online`
network and act greedily in evaluation.

PR-DRL/DDQN action selection during training is epsilon-greedy. Epsilon is read from
the configured `[epsilon]` schedule and reaches zero in test/evaluation mode.
The online network chooses the greedy next action for Double-DQN target
construction, while the target network evaluates that selected action. This
decouples action selection from target evaluation and reduces the optimistic
bias of plain DQN.

The SpaceAGORA physics PR-DRL update path is:

1. launch parallel randomized campaign workers,
2. receive one apoapsis transition from whichever worker finishes a pass,
3. push that transition into the bounded replay buffer,
4. once replay size, `train_start`, and `train_frequency` conditions are met,
   sample a minibatch,
5. evaluate online and target networks on next observations,
6. compute scalar TD targets,
7. compute network loss/gradients for the selected actions,
8. clip gradients by global norm,
9. apply the package-local Adam optimizer,
10. periodically copy online weights into the target network.

`global_step` counts environment transitions, not optimizer updates. For PR-DRL/DDQN,
`train_steps` counts replay minibatch updates. With `train_frequency = 5`, a
long run can have roughly one optimizer update per five environment steps after
`train_start`, subject to replay warmup.

PR-DRL/DDQN checkpoints contain the online network, target network, replay/training
counters, Adam state, epsilon schedule, configured device name, action table,
and the run manifest. Checkpoint filenames are keyed by global step, even for
campaign-limited runs, because training updates are scheduled by transition
count.

## A2C actor-critic implementation

`A2CLearner` uses two `QNetwork` instances:

- `actor`: maps observations to action logits over the 13 maneuver actions.
- `critic`: maps observations to a scalar value estimate.

The A2C training loop keeps one active campaign per worker thread. For each
segment it:

1. takes a CPU snapshot of the actor for rollout action selection,
2. steps active campaigns in parallel for up to `a2c.segment_length` decisions,
3. bootstraps unfinished campaigns with the critic,
4. computes discounted returns over valid worker/time entries,
5. flattens the segment into an `A2CRolloutBatch`, and
6. updates actor and critic when `global_step >= a2c.train_start`.

The actor loss uses sampled-action log probabilities weighted by normalized
advantages, plus entropy regularization controlled by `a2c.entropy_coef`. The
critic loss is a value regression term weighted by `a2c.value_coef`. Actor and
critic gradients are clipped separately and updated with independent Adam
states.

`GreedyA2CPolicy` loads the checkpoint `:actor` network and uses deterministic
argmax action selection for evaluation.

A2C action selection during training samples from the actor's softmax
distribution. Evaluation uses deterministic argmax through `GreedyA2CPolicy`.
The critic estimates state value, not action value. Its value estimate is used
both to bootstrap unfinished rollout segments and to compute advantages for the
actor update.

The A2C update path is:

1. keep one active campaign per rollout worker,
2. snapshot the actor to CPU for rollout action selection,
3. step workers for up to `a2c.segment_length` passage decisions,
4. bootstrap unfinished campaigns with the critic value estimate,
5. compute discounted returns over valid worker/time entries,
6. flatten the segment into an `A2CRolloutBatch`,
7. compute raw advantages as `returns - values`,
8. scale advantages for numerical stability,
9. update the actor with selected-action log probabilities and entropy
   regularization,
10. update the critic with value-regression gradients,
11. clip actor and critic gradients separately,
12. apply independent Adam optimizers to actor and critic.

`global_step` counts environment transitions. `train_steps` counts rollout
segment updates, so it is usually much smaller than `global_step` and depends
on `segment_length` and how many valid worker steps were collected.

A2C checkpoints contain the actor, critic, independent Adam states, training
counters, latest policy/value/entropy diagnostics, configured device name,
action table, and the run manifest. A2C does not use a replay buffer; training
updates consume fresh rollout segments.

## Evaluation and reports

Evaluation runs policies without learner updates. For each requested episode,
the evaluator resets a campaign with seed `seed + episode - 1`, repeatedly asks
the policy for an action, steps the scenario, records the transition, and
finalizes an `EpisodeSummary` when the campaign terminates or truncates.

Baseline evaluation writes CSV artifacts for `NoManeuverPolicy`,
`RandomActionPolicy`, `FixedCorridorPolicy`, and `AADSHeuristicPolicy`:

```text
julia --project=SpaceAGORA_RL.jl SpaceAGORA_RL.jl/scripts/run_baselines.jl \
  SpaceAGORA_RL.jl/configs/aerobraking/paper_replication.toml
```

The Mars Odyssey example writes CSV summaries, plots, and an interactive orbit
HTML view. By default it uses the paper-style perturbed PR-DRL/AADS comparison
distribution with full `backend_mode = :spaceagora_physics`:

```text
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/examples/aerobraking/example_SpaceAGORA_RL_MarsOdyssey.jl
```

To compare a trained PR-DRL checkpoint against the AADS heuristic from Julia:

```julia
include("SpaceAGORA_RL.jl/examples/aerobraking/example_SpaceAGORA_RL_MarsOdyssey.jl")
run_trained_pr_drl_comparison("SpaceAGORA_RL.jl/outputs/runs/<run_id>/checkpoint_final.jls")
```

The trained comparison also writes paper-shaped artifacts:

- `paper_table_v_pr_drl_vs_aads.csv`: Table V-style reward, thermal-violation,
  and goal-percentage metrics.
- `paper_fig10_heat_rate_and_delta_v.png`: Fig. 10-style heat-rate corridor and
  maneuver traces for the first ten simulations per policy.
- `paper_fig11_flight_performance_summary.png`: Fig. 11-style mean/std summary
  for maneuvers, duration, total ΔV, and thermal violations.
- `paper_fig12_heat_rate_example.png`: Fig. 12-style heat-rate trace for a
  representative simulated policy episode.

The Fig. 12-style artifact does not fabricate the 2001 Odyssey heat-rate trace;
add that flight heat-rate series before using it as a direct mission-data overlay.

This path requires the native GRAM library to be built for the host platform.
If `libGRAM.so` is missing on Linux, run
`julia --project=SpaceAGORA.jl SpaceAGORA.jl/scripts/ensure_gram_native.jl`
after installing the required build toolchain.

For the Odyssey flight-geometry comparison from the paper, use:

```julia
run_odyssey_flight_comparison("SpaceAGORA_RL.jl/outputs/runs/<run_id>/checkpoint_final.jls")
```

For A2C checkpoints, load with `load_trained_a2c_policy` and evaluate with
`evaluate_policy`, or adapt the example comparison helper to use the A2C loader.

Per-episode metrics include policy name, seed, pass count, total reward,
success/failure flags, target error, mission duration, total delta-V, maneuver
count, peak heat rate, minimum periapsis, and solver failures. `total_delta_v_mps`
and `total_mission_delta_v_mps` include ABMs, the terminal apoapsis correction,
and the final periapsis-raise maneuver. `abm_delta_v_mps`,
`apoapsis_correction_delta_v_mps`, and `periapsis_raise_delta_v_mps` expose the
components separately. Aggregate metrics summarize success rate, mean absolute
target error, total mission delta-V, its component means, mean pass count, mean
thermal violations, and total solver failures.

`pass_logs.csv` stores one row per pass with heat rate, action delta-V,
apoapsis radius, periapsis altitude, orbital angles, inclination, and reward.
Those rows drive the plotting and interactive trajectory helpers in the Mars
Odyssey example.

## Tests

Run the RL unit tests from the main `SpaceAGORA.jl` package directory:

```text
julia --project=SpaceAGORA_RL.jl SpaceAGORA_RL.jl/tests/runtests.jl
```

The test suite covers action mapping, observation normalization, reward and
termination behavior, replay-buffer wraparound, DDQN target construction, DDQN
and A2C learner updates, run-manifest naming, and backend determinism.

## Outputs

Each training run creates `outputs/runs/<run_id>/` by default, where generated
run IDs include the config name and algorithm, such as
`20260628T120000_paper_replication-pr_drl` or
`20260628T120000_paper_replication-a2c`. Runs write `manifest.toml`, periodic
`checkpoint_<step>.jls` files when
`training.checkpoint_frequency > 0`, `terminal_output.txt` from the training
script's stdout/stderr stream, and `checkpoint_final.jls` at completion.
The manifest records the generated title, source config path and hash when
available, algorithm, requested device, action table, normalization names,
reward settings, and the DDQN/A2C hyperparameters needed to interpret the run.

Important output fields:

- `manifest.toml`: immutable run context for reproducing or interpreting the
  checkpoint.
- `terminal_output.txt`: terminal output and captured exception stacktraces from
  `scripts/train.jl`.
- `checkpoint_<step>.jls`: serialized learner state at a global transition
  count. In campaign-limited runs, this number can be much larger than the
  campaign count.
- `checkpoint_final.jls`: final serialized learner state.
- `episode_metrics.csv`: per-campaign evaluation metrics when report helpers
  are used.
- `summary_metrics.csv`: aggregate evaluation metrics by policy.
- `pass_logs.csv`: per-pass traces used for plots and trajectory visualizations.

Progress lines report recent rolling metrics over completed campaigns. For
example, `recent_success=96.0%` is the recent campaign success rate, and
`recent_passes=25.4` is the recent average number of passes per completed
campaign. `steps` is always transition/action count; `ep` is completed campaign
count.
