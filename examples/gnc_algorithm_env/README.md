# GNC algorithm examples

The worked examples of
[Writing a GNC Algorithm](../../docs/src/user/gnc_algorithm.md). CI runs them
from this separate project, which uses only the public `SpaceAGORA` API.

| Script | What it shows | Page |
|---|---|---|
| `feedback_controller.jl` | A controller whose tick computes and holds a torque, an output hook that applies it, its declarations, a run checked against a hand calculation, and its conformance checks. | [Feedback Controllers and Estimators](../../docs/src/user/gnc_feedback.md) |
| `estimator.jl` | An estimator that holds a gyro's rate as its estimate, a controller that reads it through a declared link, and the conformance checks of both. | [Feedback Controllers and Estimators](../../docs/src/user/gnc_feedback.md) |
| `rpo_planner.jl` | A planner of the RPO reference family, run in the RPO pilot and checked with the family's conformance requests. | [Planning Algorithms](../../docs/src/user/gnc_planning.md) |
| `burn_planner.jl` | A planner of the burn-plan family, run through `AlgorithmGuidance` with the burn scheduler; its propellant use is checked against the rocket equation. | [Planning Algorithms](../../docs/src/user/gnc_planning.md) |
| `joint_planner.jl` | A planner of the joint-reference family, run through `AlgorithmGuidance` with a robot arm's controller. | [Planning Algorithms](../../docs/src/user/gnc_planning.md) |

From the repository root, prepare the environment once:

```sh
julia --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/setup.jl
```

Setup copies the tracked `Project.template.toml` to a local `Project.toml`,
starts `Manifest.toml` from the repository manifest so the recorded package
versions are kept, and adds this checkout by path. Git ignores both generated
files. On a fresh Julia depot setup downloads the packages it lacks, so it needs
network access the first time. After that the examples run offline, as in CI,
with no optional package, native GRAM or SPICE kernels. Run one example, or all
of them as CI does:

```sh
julia --project=examples/gnc_algorithm_env --depwarn=error examples/gnc_algorithm_env/feedback_controller.jl
JULIA_NUM_THREADS=2 julia --project=examples/gnc_algorithm_env --depwarn=error examples/gnc_algorithm_env/smoke.jl
```

Each script checks its own run, by hand or with the family's validation, and
ends with SpaceAGORA's conformance checks (`SpaceAGORA.Conformance`). With two
threads the checks also compare concurrent runs with sequential ones. Passing
them shows that the algorithm follows the interface, not that the method is
good; the examples are deliberately small.

Each region between a `# snippet:` line and its `# end snippet` line appears
verbatim on the page listed above. `test/gates/ci_hpc_extensibility_docs_gate.jl`
compares the two, so change the script and the page together.
