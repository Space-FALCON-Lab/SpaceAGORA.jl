# Optional KS propagation candidate

This route selects a separate KSPropagator package through
`SolverConfig(solver_mode=:ks)`. SpaceAGORA owns physical state, time, forces,
callbacks and output. KSPropagator owns its coordinates, equations and numerical
stepping. The extension uses SpaceAGORA's configured full Cartesian RHS;
standalone KS force models do not supply integrated spacecraft physics.
Mohammed's MPC is unchanged.

This is a draft development candidate. The [report](REPORT.md) and
[capability table](CAPABILITIES.md) document current-upstream trajectory, lifecycle,
installation and guard verdicts, including all eight coarse/medium misses.
Performance and final delivery review remain pending. Neither this draft nor
a package unit-test pass is a general accuracy or performance claim. The
separate checkpoint-table repair is a prerequisite.

## Installation and selection

Use Julia 1.12.7 and a checkout of this candidate, with a fresh environment.
Use normal Git/LFS checkout, then the installer below. A direct Pkg.add of the
SpaceAGORA Git URL failed with a LibGit2 submodule checkout error in the tested
environment; the working checkout/installer route is the verified route.
The installer fetches KS source `05d969c24a6595c2752f7ff947af77c504053623`
(version 0.2.0) from its maintained repository. The KS repository is private,
so normal authorized GitHub access is required. It is not a registry release.
The published dependency inventory is in KS draft PR #24; intermediate PR
#20–22 heads contain known defects and are not installation sources.

```sh
JULIA_LOAD_PATH=@:@stdlib julia --startup-file=no verification/optional_ks/install.jl /absolute/new/ks-environment
JULIA_LOAD_PATH=@:@stdlib julia --startup-file=no verification/optional_ks/install.jl /absolute/new/ordinary-environment --without-ks
```

In the first environment, use `using SpaceAGORA, KSPropagator` and set
`solver_config=SolverConfig(solver_mode=:ks)` on the normal configuration.
Physical SI state/output stays in SpaceAGORA. No KS installation is needed for
ordinary solvers. A missing package, missing stepping API or incompatible
package version must fail explicitly. OrdinaryDiffEqCore is restricted to the
4.19 minor series because the extension implements its algorithm interface.

## Initial admission boundary

One spacecraft, one rigid link, orientation propagation disabled, and exactly
one active inverse-square gravity model (optionally J2) are required. Admitted
environments are NoAtmosphere or deterministic ExponentialAtmosphere with
SimpleEphemerides, MaxwellianHeat and no topography. Constant-coefficient drag
and ideal BaseThruster use the configured force/heat/mass interfaces.

Mounted hardware including power, data stores/sources and sensors, articulated
vehicles, interlinks, external propagation and unsupported force models are
refused. Repeatable custom RHS models must explicitly implement
`SpaceAGORA.trial_rhs_safe(model)=true`; that is a purity contract, not physical
validation. Discrete accepted-step controls may retain history. Nonlinear
trials must not advance physical, random or controller history.

For state-dependent atmosphere evaluation set
`SPACEAGORA_DENSITY_FREEZE_PER_STEP=false` and
`SPACEAGORA_GRAM_TRACK_CACHE_MODE=off`; history-dependent and vacuum-cache
semantics are outside this initial route. The fixture pins these settings.

The numerical method is a canonical-reset implicit KS midpoint map, with a
positive fictitious step and a physical endpoint constraint. Step doubling
accepts the fine state, without extrapolation; auxiliary mass/heat share its
nonlinear solve and error decision. Saved dense intervals retain incoming
endpoint states/rates across command changes and impulses. Qualification of
other package solvers or sensitivities does not qualify this new map.

## Reproduction

Run with `JULIA_LOAD_PATH=@:@stdlib` and the installed environment. Use fresh
output directories; preserve failed attempts. The twelve baseline cases retain
the original physical inputs and all coarse/medium/fine settings, including
the eight historical coarse/medium accuracy misses. Targets remain maximum
1 m and 0.001 m/s on 1,441 physical epochs over 24 h. The historical fine
settings (relative 1e-12, absolute 1e-13, physical maximum step 5 s) are case
settings, not a universal prescription.

```sh
julia --startup-file=no --project=/absolute/new/ks-environment verification/optional_ks/current-checks.jl
KS_BASELINE_OUTPUT=/absolute/new/baseline julia --startup-file=no --project=/absolute/new/ks-environment verification/optional_ks/baseline-matrix.jl
julia --startup-file=no --project=/absolute/new/ks-environment verification/optional_ks/run-lifecycle.jl /absolute/new/lifecycle
julia --startup-file=no --project=/absolute/new/ordinary-environment verification/optional_ks/optional-absence.jl
```

The lifecycle workers use separate processes for uninterrupted, segmented,
interrupted and resumed runs, including coincident action/control/output/
checkpoint boundaries. Current upstream restores valid persisted output
history: resume must retain that prefix and the future nominal grid without
duplicate timestamps. Missing/incompatible output history follows core's
separate documented fallback. TrajectoryRecorder cadence and specialized
split/multirate output paths are not qualified by these table-output tests.

Integrated solar sails, multiple spacecraft, interlinks and new fixed-physical-
time sensitivities remain deferred. Expand support only with separate declared
verification. Compilation and warm execution must be reported separately in
any subsequent matched-accuracy timing study.
