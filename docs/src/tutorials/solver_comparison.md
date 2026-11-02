# Compare Solvers on a Mounted-Actuator Case

This exercise measures how solver settings change one fixed simulation. You
will run a CubeSat with reaction wheels, torque rods and a laser terminal,
compare the usual Tsit5 method and multirate MRI against DP8, and refine the
reference before interpreting the differences. It is a numerical exercise on
this case, not a qualification of every spacecraft or a new propagator.

## Run the comparison

Use Julia 1.12. From the repository root, in a terminal on Windows, Linux or
macOS:

```text
julia --startup-file=no examples/setup.jl
julia --startup-file=no --project=examples examples/Mounted_Actuator_Solver_Comparison.jl
```

The first command prepares the examples environment. The second runs all nine
configurations and prints a new `output/solver_comparison/attempt-…/` directory
for this attempt. Earlier attempts remain intact. It needs no GRAM, SPICE
or other external assets. Initial compilation can take longer than the runs.
The script is included in the examples smoke suite; smoke mode uses the same
120-second comparison, not a shortened accuracy exercise.

## What stays fixed

The setup is shared with `Earth_CubeSat_Mounted_Actuators.jl` through
`examples/support/mounted_actuator_case.jl`. Each run builds fresh spacecraft,
actuator and controller objects, so no run inherits another run's wheel
momenta or held commands.

- A 10 kg rigid CubeSat starts in a roughly 500 km Earth orbit, with body rates
  `[0.05, -0.05, 0.05]` rad/s and identity attitude.
- Three wheels detumble it; three torque rods unload the wheels. An optical
  terminal transmits 1 W. The original gains, limits, epoch and force model
  are unchanged.
- Every run lasts 120 s, with control updates every 1 s. The comparison samples
  each returned solution at `0:0.25:120` s, including intermediate trajectory
  points rather than only the endpoint.
- No atmosphere, burns, contacts, interlinks or checkpoint interruptions are
  introduced. This case does not test their event behavior.

The default *method* is pinned to `:tsit5`; the exercise deliberately sets its
step cap and tolerances instead of inheriting shell overrides. For the
candidate runs, orbit, attitude and angular-rate relative tolerances are
`1e-8`, and absolute tolerances are `1e-10`. This is not a replay of every
setting of the original demo.

## Which settings change

| Run | Method | Outer/maximum step cap (s) | Fast substeps |
|---|---|---:|---:|
| Reference | DP8 | 0.03125 | not used |
| Refined reference | DP8 | 0.015625 | not used |
| Default method | Tsit5 | 1.0 | not used |
| Substep sweep | MRI | 1.0 | 1, 2, 4, 8 |
| Outer-step sweep | MRI | 0.5, 0.25 | 8 |

The reference uses relative/absolute tolerances `1e-13`/`1e-14`; the refined
reference uses `1e-14`/`1e-15`. Both use the same initial state, control schedule
and sampling grid as the candidates.

MRI's substep count is per stage, not per outer step. The outer solver still
adapts and may take steps shorter than its cap. The substep sweep holds the cap
fixed; the outer-step sweep holds the substep count fixed. A tighter slow
relative tolerance alone does not bound fast-state error. Read
[Coupled multirate integration](../user/solver_configuration.md#Coupled-multirate-integration)
for the partition, dense-output caveats and unsupported cases.

## Read the evidence

In the printed attempt directory, `summary.csv` has one row per candidate and quantity. Its `max_error` and
`rms_error` compare the candidate with the refined reference over the common
grid. For vector quantities, each sample's error is the Euclidean norm;
`rms_error` is the root mean square of those norms over time.

| Metric | Meaning and unit |
|---|---|
| `position_m` | Inertial position difference, m |
| `velocity_m_s` | Inertial velocity difference, m/s |
| `attitude_rad` | Smallest relative rotation angle, rad; equivalent quaternion signs are aligned |
| `body_rate_rad_s` | Body angular-velocity difference, rad/s |
| `wheel_momentum_n_m_s` | Integrated reaction-wheel momentum difference, N m s |

The per-run CSVs retain all sampled errors. `reference_refinement.csv` compares
the two DP8 references; its maximum for each metric is repeated in the summary
as `reference_change`. `settings.toml` records the revision, uncommitted-file
status, Julia version, manifest hashes, settings, environment and actual solver
traces. Its status starts as `running` and becomes `complete` only after every
comparison succeeds, or `failed` with the error if an exception is raised.
An interrupted process may leave `running`; do not treat partial output as a
completed exercise. Keep these files together. A dirty working tree is not fully identified
by its commit: retain your diff alongside the results.

Read `reference_resolved` before ranking the methods. It is true only when the
reference change is at most one tenth of that candidate's maximum difference.
That factor is an exercise diagnostic, **not an accepted mission tolerance or
a proof of the reference's true error**. A false value means the reference is
not fine enough to interpret that difference under this rule; refine DP8
further or report the quantity as unresolved. Differences at roundoff level
may not resolve under further refinement.

Every run must finish at 120 s with `Success` and finite samples before its
comparison is written. Those checks establish that a run completed, not that
it is accurate. Look for decreasing error in each sweep; plateaus, increases
and unresolved quantities are findings to report, not values to omit. The
recorded elapsed time includes compilation and setup effects, so it is not a
fair solver speed benchmark.

## What to submit for review

Keep the configuration and all output files, then write a short conclusion:

1. Record the code revision, any local edits, platform, Julia version and
   commands used. Confirm that all runs used the same physical case and clocks.
2. Show the reference-refinement differences before the candidate results.
   Identify quantities for which the reference is too uncertain to rank runs.
3. Compare the substep sweep and the outer-cap sweep separately, quantity by
   quantity. Include units and both maximum and RMS errors. Do not select only
   the best endpoint or silently remove a failed run.
4. Explain which settings improve accuracy on this case, which do not, and
   which questions remain. Record any instruction you could not follow without
   help.

Review the physical setup and the numerical conclusion separately. Event and
checkpoint tests need dedicated cases after this exercise. Developing or
qualifying Gauss–Jackson is separate work; completing this comparison does not
establish that a new propagator works.
