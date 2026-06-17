# SpaceAGORA_RL

Julia-native reinforcement-learning utilities for SpaceAGORA aerobraking studies.

The first task namespace implements the paper-compatible Mars Odyssey aerobraking MDP:
the 9-feature observation schema, 13-action apoapsis maneuver table, reward and
termination rules, deterministic paper-mode backend facade, baseline policies,
DDQN learner infrastructure, CSV/report helpers, and an example runner.

Reference material is kept under `Reference Code/` and is not on the runtime path.
Generated artifacts should go under `outputs/`; only placeholders are committed.

Training uses `training.global_steps` as the primary stopping budget when that
value is positive. This matches the paper's step-based DDQN stopping condition:
one global step is one environment transition/action step, not one full
aerobraking campaign. The `training.episodes` field is still used by evaluation
scripts and by explicit episode-capped training calls.
