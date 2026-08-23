# HyPR-RL

HyPR-RL is a baseline RPO path planner that replaces HyPR's full PSO search
with a masked-DDQN sequential editor. It starts from a direct or
RRT-Connect/Bezier path, then edits translation control points, attitude knots,
and both waypoint counts. Every accepted edit is evaluated after retiming.

The coupled evaluator runs the existing SpaceAGORA LQ-MPC controller while it
propagates reaction-wheel attitude tracking and body-frame six-thruster pulse
allocation. Thrusters provide translation only; reaction wheels provide
attitude torque. Propellant is accumulated per thruster from impulse and Isp,
including minimum firing time. The best feasible path—including the original
seed—is retained throughout an episode.

The task is organized as follows:

- `types.jl`: task configuration, scenario, state, evaluation, and plan types.
- `mdp/actions.jl`: fixed action catalog and state-dependent action masks.
- `mdp/environment.jl`: observations, edits, rewards, and termination.
- `attitude.jl`: quaternion interpolation and reaction-wheel propagation.
- `backend/spaceagora_adapter.jl`: seed planning, retiming, LQ-MPC, and actuator accounting.
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

At inference time, construct `RPOHyPRRLMDP(config, scenario)`, load a frozen
policy with `load_hypr_rl_policy(checkpoint)`, and call
`hypr_rl_plan(mdp, policy)`. `rpo_hypr_rl_plan_path` provides the equivalent
start/goal/geometry convenience API. The returned `RPOHyPRRLPlan` contains synchronized
translation, velocity, and attitude references plus actuator diagnostics.
