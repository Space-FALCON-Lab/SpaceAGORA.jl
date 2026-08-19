# RL evaluation scripts

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
