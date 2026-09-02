# RL evaluation scripts

## HyPR-RL RPO path planning

The RPO trainer uses masked DDQN to smooth and edit raw RRT-Connect topology.
Candidate edits directly move, insert, and delete the original HyPR Bezier
control waypoints, with up to 20 internal translation waypoints, and also use
variable attitude waypoint counts.
Unsafe intermediate curves receive finite obstacle scores; feasible curves are
scored after retiming, reaction-wheel attitude propagation, and six-thruster
propellant accounting. The feasible objective contains normalized propellant
use plus a small normalized reaction-wheel energy term; it has no path-length
term. Both the RL editor and retimed-fuel PSO use the original weighted HyPR
sigmoid clearance penalty for unsafe candidates. Full LQ-MPC runs once on each
episode's best terminal candidate and
supplies a bounded terminal reward while hard clearance checks remain in force.
It uses the same PR-DRL rollout architecture as aerobraking: frozen policy
snapshots on parallel workers and a single central learner. The supplied
configuration launches 16 isolated one-thread worker processes.
The configured run lasts 100,000 episodes. Epsilon is 1.0 through episode
10,000, decays linearly to 0.01 at episode 20,000, and remains there afterward.
Each successful base scenario is replayed three additional times. Training
progress includes both elapsed time and ETA. Training and evaluation draw from
the same station-wide near-surface endpoint distribution as the HyPR comparison
cases: 0.55--1.0 m clearance, 1.5 m minimum separation, and the same surrounded
endpoint rejection rule.

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/rpo/train_hypr_rl.jl \
  SpaceAGORA_RL.jl/configs/rpo/hypr_rl.toml
```

Evaluate the newest frozen checkpoint on 100 held-out randomized RPO cases:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/rpo/evaluate_hypr_rl.jl
```

The output includes per-case and aggregate clearance, safety-violation, fuel,
tracking, actuator, attitude-wheel, and runtime metrics. `index.html` links to
an interactive 3D HTML trajectory for every case, with the existing station
CAD mesh displayed beneath the sampled Bezier plan and executed LQ-MPC path.
Control waypoints are shown as markers rather than connected as a surrogate
trajectory. The evaluation uses the
configured 16 isolated processes and applies full LQ-MPC only once after the
learned sequential editor selects its best candidate.

The RPO code and API map are documented in `src/tasks/rpo/README.md`.

Run the matched baseline HyPR evaluation on the same held-out cases with:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/rpo/evaluate_hypr_baseline.jl
```

This version retimes every PSO particle and makes PSO minimize the same
six-thruster propellant plus small reaction-wheel energy objective used by the
RL editor. It runs full LQ-MPC once on the winning path and writes the same CSV
and interactive CAD-mesh trajectory artifacts. The exact per-particle
retiming makes it considerably slower than the legacy HyPR proxy objective.

Compare that objective directly against original HyPR on paired cases with:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/rpo/compare_hypr_objectives.jl
```

The comparison uses the same scenario and PSO seed for both planners, evaluates
both selected paths with the same full terminal model, and writes paired fuel,
safety, path-length, and runtime results. Its HTML index links to a two-panel 3D
trajectory comparison for each case.

## Parallel A2C training

The A2C trainer uses synchronous, policy-versioned rollouts. The supplied
configuration matches the PR-DRL SpaceAGORA-physics setup: native MarsGRAM,
perturbed winds, and one isolated OS process per environment worker. From the
SpaceAGORA repository root, run:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/train.jl \
  SpaceAGORA_RL.jl/configs/aerobraking/a2c_spaceagora_physics_marsgram.toml
```

For live SpaceAGORA backends, `training.worker_backend` may be `threads` or
`processes`. The latter uses isolated, long-lived physics workers while the
central learner updates only after every active worker reaches the rollout
barrier.

Both A2C and A3C accept `action_space = "discrete"` or `"continuous"` in
their algorithm table. Continuous mode uses a tanh-squashed Gaussian policy
whose commanded delta-v is bounded by the outer paper-action limits of
`-1.0` and `1.0` m/s. `initial_log_std`, `log_std_min`, and `log_std_max`
control exploration. The supplied live MarsGRAM A2C and A3C configurations
currently select continuous mode; configurations without this option retain
the original discrete policy.

## Asynchronous A3C training

The A3C trainer keeps an independent actor, critic, rollout, and policy version
for each environment worker. A worker applies its n-step update as soon as its
rollout reaches `a3c.t_max` or its campaign ends, then refreshes only its own
local model. Other workers continue propagating and do not wait at a rollout
barrier. Native MarsGRAM continues to use one isolated OS process per worker,
while optimizer updates are serialized in the parent process.

From the SpaceAGORA repository root, run:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/train.jl \
  SpaceAGORA_RL.jl/configs/aerobraking/a3c_spaceagora_physics_marsgram.toml
```

## Off-policy TD3 training

TD3 uses the same continuous `-1.0` to `1.0` m/s aerobraking action range and
the same isolated live MarsGRAM worker layout. A central CUDA learner maintains
the deterministic actor, twin critics, target networks, and exact continuous
action replay buffer. Run the matched campaign with:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/train.jl \
  SpaceAGORA_RL.jl/configs/aerobraking/td3_spaceagora_physics_marsgram.toml
```

## Multi-run policy comparison

Use `--compare-run` once per additional training run to compare PR-DRL and A2C
on the same thermal-tolerant campaigns alongside AADS and the Mars Odyssey
reference:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/evaluate_rl_run.jl \
  SpaceAGORA_RL.jl/outputs/runs/PR_DRL_RUN \
  --compare-run SpaceAGORA_RL.jl/outputs/runs/A2C_RUN \
  --output SpaceAGORA_RL.jl/outputs/comparisons/pr_drl_vs_a2c
```

The first run supplies the common evaluation scenario and campaign seeds. Each
run loads its own recorded configuration and best validation checkpoint, with
the final/latest checkpoint used when no best-checkpoint record exists. The
combined output includes raw episode and pass data, a performance CSV and
figure, and a manifest identifying every source run and checkpoint.

`evaluate_rl_run.jl` turns a completed SpaceAGORA_RL run into the evaluation
metrics and figures used in Falcone and Putnam (2023), DOI
`10.1109/TAES.2022.3221697`.

From the SpaceAGORA repository root:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/evaluate_rl_run.jl \
  SpaceAGORA_RL.jl/outputs/runs/MY_RUN
```

The script reads `manifest.toml`, resolves the training config, and selects the
checkpoint recorded in
`checkpoint_validation/best_validation_checkpoint.txt` by default. For older
runs without that validation artifact, it falls back to `checkpoint_final.jls`
and then the largest numbered checkpoint. Use `--checkpoint best`,
`--checkpoint final`, or `--checkpoint PATH` to choose explicitly, and use
`--output` to put the artifacts elsewhere.

The default paper protocol runs:

- 40 IID randomized episodes for PR-DRL versus AADS (Table V and Fig. 10);
- 100 Odyssey-geometry episodes (Figs. 11-12);
- 100 episodes in the held-out IID reference and each of six SpaceAGORA-native
  generalization cases (Table VI quantities);
- 40 conservative and 40 tolerant episodes for every frozen checkpoint
  (Figs. 7-8);
- repeated-action analysis of the selected policy (Fig. 9).

The checkpoint sweep also writes
`episode_completion_and_final_target_distance.png`, a dedicated two-panel
chart of episode completion percentages and mean absolute final distance to
the target apoapsis radius with a one-standard-deviation band.

For a focused evaluation of only the selected checkpoint, run
`--final-flight-comparison`. Its 40-campaign thermal-tolerant comparison CSV,
figure, and terminal summary include the count and percentage of campaigns
finishing within ±10 km of the target apoapsis radius. The figure also compares
the mean absolute final distance from that target with one-standard-deviation
error bars.

Raw episode, summary, and pass-level CSVs are retained alongside the figures.
`evaluation_manifest.toml` records the exact run, config, checkpoint, episode
counts, case mappings, reference values, and artifact paths.

The `generalization_evaluation_suite` output directory contains the six-case
Table VI view, an all-cases table that includes the IID loss reference, raw
episode/pass data, the resolved SpaceAGORA configuration for every case, and
an automatically rendered transposed `generalization_results_table.pdf`.
The evaluated checkpoint is frozen and selected greedily; no case performs
training or fine-tuning. While the suite runs, `progress.toml` records the
current case, case episode, overall episode count and percentage, elapsed time,
coordinator process ID, worker count, and threads per worker. Each case's
campaigns run concurrently in isolated Julia worker processes; the same worker
pool is reused across cases. The default is 16 workers with one thread each.
Terminal updates follow `--progress-every N` (default: every episode).

Run that suite without the other paper figures and comparisons with:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/evaluate_rl_run.jl \
  SpaceAGORA_RL.jl/outputs/runs/MY_RUN \
  --generalization-only \
  --processes 16 --threads-per-process 1
```

For a quick smoke evaluation:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/evaluate_rl_run.jl \
  SpaceAGORA_RL.jl/outputs/runs/MY_RUN \
  --episodes 1 --flight-episodes 1 --generalization-episodes 1 \
  --skip-checkpoint-sweep
```

Run `evaluate_rl_run.jl --help` for all options.

## MarsGRAM wind modes

Training and evaluation accept `--wind-mode zero|nominal|perturbed`:

- `zero` returns a literal zero wind vector;
- `nominal` uses MarsGRAM mean winds;
- `perturbed` uses MarsGRAM perturbed winds.

The command-line option overrides `spaceagora_physics.gram_wind_mode` from the
TOML configuration. The resolved mode is written to the run or evaluation
manifest.
