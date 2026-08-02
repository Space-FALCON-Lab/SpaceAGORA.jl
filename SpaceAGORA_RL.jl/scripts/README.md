# RL evaluation scripts

`evaluate_rl_run.jl` turns a completed SpaceAGORA_RL run into the evaluation
metrics and figures used in Falcone and Putnam (2023), DOI
`10.1109/TAES.2022.3221697`.

From the SpaceAGORA repository root:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/evaluate_rl_run.jl \
  SpaceAGORA_RL.jl/outputs/runs/MY_RUN
```

The script reads `manifest.toml`, resolves the training config, selects
`checkpoint_final.jls` by default, and writes `MY_RUN/paper_evaluation/`.
Use `--checkpoint` to evaluate another frozen policy and `--output` to put the
artifacts elsewhere.

The default paper protocol runs:

- 40 IID randomized episodes for PR-DRL versus AADS (Table V and Fig. 10);
- 100 Odyssey-geometry episodes (Figs. 11-12);
- 100 episodes in each of six generalization cases (Table VI);
- 40 conservative and 40 tolerant episodes for every frozen checkpoint
  (Figs. 7-8);
- repeated-action analysis of the selected policy (Fig. 9).

The checkpoint sweep also writes
`episode_completion_and_final_target_distance.png`, a dedicated two-panel
chart of episode completion percentages and mean absolute final distance to
the target apoapsis radius with a one-standard-deviation band.

For a focused evaluation of only the final checkpoint, run
`--final-flight-comparison`. Its 40-campaign thermal-tolerant comparison CSV,
figure, and terminal summary include the count and percentage of campaigns
finishing within ±10 km of the target apoapsis radius. The figure also compares
the mean absolute final distance from that target with one-standard-deviation
error bars.

Raw episode, summary, and pass-level CSVs are retained alongside the figures.
`evaluation_manifest.toml` records the exact run, config, checkpoint, episode
counts, case mappings, reference values, and artifact paths.

For a quick smoke evaluation:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/evaluate_rl_run.jl \
  SpaceAGORA_RL.jl/outputs/runs/MY_RUN \
  --episodes 1 --flight-episodes 1 --generalization-episodes 1 \
  --skip-checkpoint-sweep
```

Run `evaluate_rl_run.jl --help` for all options.
