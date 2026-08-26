# HyPR-RL

HyPR-RL is a baseline RPO path planner that replaces HyPR's full PSO search
with a masked-DDQN sequential editor. RRT-Connect supplies raw path topology;
the policy then controls the internal Bezier waypoints used by the original
HyPR representation, attitude knots, and both waypoint counts. Start and goal
remain fixed. Unsafe intermediate curves receive finite obstacle scores so the
policy can correct them over multiple edits. Unsafe candidates use the original
weighted HyPR sigmoid clearance penalty; retimed fuel and reaction-wheel energy
rank candidates only after they satisfy the hard clearance requirement. The
supplied configuration permits up to 20 internal Bezier waypoints and 64
sequential edits.

During training, each editor step uses the retimed reference acceleration for
reaction-wheel attitude propagation and body-frame six-thruster pulse
allocation. Thrusters provide translation only; reaction wheels provide
attitude torque. Propellant is accumulated per thruster from impulse and Isp,
including minimum firing time. At episode termination, the best candidate is
evaluated once with the existing SpaceAGORA LQ-MPC controller and the fidelity
difference is applied to the terminal reward. Once a feasible path is found,
the best feasible candidate is retained throughout the episode.
Feasible candidates are ranked only by normalized six-thruster propellant use
and a small normalized reaction-wheel energy term; path length is not part of
the objective. Objective improvements and terminal rewards are bounded, while
clearance, tracking, and wheel-momentum requirements remain hard feasibility
checks.

The task is organized as follows:

- `types.jl`: task configuration, scenario, state, evaluation, and plan types.
- `mdp/actions.jl`: fixed action catalog and state-dependent action masks.
- `mdp/environment.jl`: observations, edits, rewards, and termination.
- `attitude.jl`: quaternion interpolation and reaction-wheel propagation.
- `backend/spaceagora_adapter.jl`: seed planning, retiming, feedforward/full-LQ-MPC evaluation, and actuator accounting.
- `planner.jl`: inference and checkpoint-backed planner entry points.
- `training/training.jl`: masked-DDQN training loop and checkpoints.

The initial implementation targets baseline path planning. Dynamic replanning
can reuse the same editor by constructing a new scenario from the current
state and remaining path; it is intentionally not mixed into this first task.

Train from the repository root with:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/rpo/train_hypr_rl.jl \
  SpaceAGORA_RL.jl/configs/rpo/hypr_rl.toml
```

Training follows the aerobraking PR-DRL execution strategy: 16 isolated
one-thread rollout processes use a frozen policy snapshot for each batch, then
a single central learner ingests the masked transitions and updates DDQN.
Successful base scenarios are replayed three additional times with the same
scenario seed, matching the aerobraking successful-case repetition semantics.
Exploration is constant at epsilon 1.0 through episode 10,000, decays linearly
to 0.01 at episode 20,000, and remains at 0.01 through episode 100,000. Progress
lines report elapsed time and an ETA based on observed episode throughput. Set
`worker_backend = "threads"` to use Julia threads instead; launch Julia with
`--threads=16` (or `JULIA_NUM_THREADS=16`) so all requested workers are active.

At inference time, construct `RPOHyPRRLMDP(config, scenario)`, load a frozen
policy with `load_hypr_rl_policy(checkpoint)`, and call
`hypr_rl_plan(mdp, policy)`. `rpo_hypr_rl_plan_path` provides the equivalent
start/goal/geometry convenience API. The returned `RPOHyPRRLPlan` contains synchronized
translation, velocity, and attitude references plus actuator diagnostics.

Evaluate the newest checkpoint on the configured 100 held-out randomized cases
with the same parallel process strategy using:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/rpo/evaluate_hypr_rl.jl
```

The evaluator freezes the policy, performs retimed-feedforward editor steps,
and runs full LQ-MPC/thruster/attitude accounting once on the selected terminal
path. It writes `cases.csv`, `summary.csv`, an evaluation manifest, an HTML
index, and one interactive 3D Plotly trajectory per case using the existing
station CAD mesh. Training and evaluation use the same area-weighted,
station-wide endpoint sampler and the HyPR comparison bounds (0.55--1.0 m
clearance, 1.5 m minimum separation, and 2.0 m surrounded-point rejection).
The positional
arguments are `CHECKPOINT CONFIG OUTPUT_DIRECTORY CASES WORKERS`; use `latest`
for the checkpoint argument to select the most recently written checkpoint.

Evaluate baseline HyPR on the same cases with PSO directly minimizing the
HyPR-RL retimed fuel and reaction-wheel objective:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/rpo/evaluate_hypr_baseline.jl
```

Every PSO particle is safety-checked, retimed, and evaluated with the same
six-thruster feedforward accounting used for RL edits. The selected path then
receives the same one-time full LQ-MPC terminal evaluation as HyPR-RL. The
baseline uses the scenario's fixed endpoint attitudes and does not add internal
attitude knots. Unsafe particles retain the original finite HyPR clearance
penalty so the swarm can move toward feasibility; they still cannot pass the
terminal safety evaluation. Its positional arguments are
`CONFIG OUTPUT_DIRECTORY CASES WORKERS`. Because retiming and actuator
accounting now occur inside the PSO loop, this matched baseline is substantially
more expensive than the original geometric/proxy-fuel HyPR evaluation.

For a paired comparison of the original and retimed-fuel objectives, run:

```sh
julia --project=SpaceAGORA_RL.jl \
  SpaceAGORA_RL.jl/scripts/rpo/compare_hypr_objectives.jl
```

Both planners receive identical scenarios, PSO seeds, swarm settings, and
one-time full LQ-MPC terminal evaluations. The output contains separate planner
CSVs, `paired_cases.csv`, fuel-savings and runtime-ratio summaries, individual
3D trajectories, and a side-by-side HTML trajectory page for every case. The
positional arguments are `CONFIG OUTPUT_DIRECTORY CASES WORKERS`.
