# GNC algorithm examples

The worked examples of
[Writing a GNC Algorithm](../../docs/src/user/gnc_algorithm.md), run by CI
from this separate project with the public `SpaceAGORA` surface only:

| Script | What it shows |
|---|---|
| `feedback_controller.jl` | a feedback controller: a tick that computes and holds a torque, an output hook that applies it, its declarations, a run checked against a closed form, and its conformance checks |
| `estimator.jl` | an estimator that holds a gyro's rate as its estimate, a controller that reads it through a declared link, and the pair's conformance checks |
| `rpo_planner.jl` | a planner of the RPO reference family, run in the RPO pilot and checked with the family's conformance requests |
| `burn_planner.jl` | a planner of the burn-plan family, run through `AlgorithmGuidance` with the burn scheduler, its propellant checked against the rocket equation |
| `joint_planner.jl` | a planner of the joint-reference family, run through `AlgorithmGuidance` with a robot arm's controller |

From the repository root, prepare the environment once:

```sh
julia --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/setup.jl
```

Setup creates a local `Project.toml` from the tracked `Project.template.toml`,
seeds `Manifest.toml` from the repository manifest and adds this checkout by
path, keeping the recorded dependency versions. Both generated files are
ignored. Setup downloads any of those versions your Julia depot does not hold
yet, so on a fresh depot it needs network access. After setup the examples run
offline, as CI runs them, and need no optional package, native GRAM or SPICE
kernels. Then run one example, or all of them as CI does:

```sh
julia --project=examples/gnc_algorithm_env --depwarn=error examples/gnc_algorithm_env/feedback_controller.jl
JULIA_NUM_THREADS=2 julia --project=examples/gnc_algorithm_env --depwarn=error examples/gnc_algorithm_env/smoke.jl
```

Each script checks its run by hand or against the family's validation, and ends
with SpaceAGORA's conformance checks (`SpaceAGORA.Conformance`). Two threads let
the checks compare concurrent runs with sequential ones. A pass certifies
interface behaviour, not the method: the examples are small on purpose.

The regions between `# snippet:` and `# end snippet` lines appear on the docs
page as they are here; `test/gates/ci_hpc_extensibility_docs_gate.jl` compares
them, so change both together.
