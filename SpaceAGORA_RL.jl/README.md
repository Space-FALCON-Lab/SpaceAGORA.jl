# SpaceAGORA_RL

Julia-native reinforcement-learning utilities for SpaceAGORA aerobraking studies.

The first task namespace implements the paper-compatible Mars Odyssey aerobraking MDP:
the 9-feature observation schema, 13-action apoapsis maneuver table, reward and
termination rules, deterministic paper-mode backend facade, baseline policies,
DDQN learner infrastructure, CSV/report helpers, and an example runner.

Reference material is kept under `Reference Code/` and is not on the runtime path.
Generated artifacts should go under `outputs/`; only placeholders are committed.
