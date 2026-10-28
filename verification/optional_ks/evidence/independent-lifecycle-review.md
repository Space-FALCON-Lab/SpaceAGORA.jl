# Independent current-upstream lifecycle review — attempt1

This review audits **new** `upstream-delivery-20261028/lifecycle-attempt1` results. It does not relabel the older prototype results. The adapter worktree is `c091191a1ae0dc602ae50025889d34d4b2800a96`, based on the current upstream and checkpoint repair615a447c; the normal environment installs KSPropagator from GitHub at exact revision `8cfc936fe2ba403996ed7682aa700b56d6e04c5d` (version0.2.0). The14 frozen lifecycle source/environment hashes match their current files. Ten newly produced raw artifacts—eight result tables and two interrupted-checkpoint snapshots—have separately retained SHA256 values.

**Lifecycle verdict:** supported within the frozen one-spacecraft test scope and unchanged acceptance budgets. **Overall delivery verdict remains pending** baseline completion, clean-install qualification, separate source review and matched-accuracy performance evidence. No performance result is inferred from these lifecycle runs.

## Independent method and provenance

A fresh Julia1.12.7 process decoded the new raw result files in the actual `adapter-environment`, using normal SpaceAGORA package definitions. No headless loader, model shim or substitute runtime was used. `export-current-lifecycle-attempt1.jl` exports portable numeric states and metadata. `audit-current-lifecycle-attempt1.py` independently integrates Cartesian position/velocity/current mass/heat in SciPy DOP853. It calls no SpaceAGORA/KS force or integration functions.

The reference includes central/J2 gravity, exponential co-rotating drag, independently expressed Maxwellian heat, the prescribed13s held inertial control force, trial-state/time force, current-mass prograde ideal thrust and its propellant flow, area changes and the velocity impulse. It segments exactly at all discontinuities, holding the incoming control through endpoint stages before switching. Fine/coarse reference runs independently refine numerical accuracy and descending-root epochs.

Template changes are limited to new locations, output names and the prospectively declared current-upstream resume-table contract: a valid saved prefix is restored together with future outputs. The physics, initial conditions, action/tick schedule and position/velocity/mass/heat/event acceptance targets are unchanged. The script lineage and new raw/source hashes are retained beside this review.

## New results

| Current-source result | Rows | Max position error (m) | Max velocity error (m/s) |
|---|---:|---:|---:|
| Main DP8 direct |354|3.564e−7|3.543e−10|
| Main KS direct |354|0.01093393|1.19122e−5|
| Main DP8 checkpoint segments |360|3.477e−7|3.625e−10|
| Main KS checkpoint segments |360|0.01093272|1.19111e−5|
| Main KS fresh-process resume |360|0.01093272|1.19111e−5|
| Coincident DP8 direct |27|6.265e−9|4.093e−12|
| Coincident KS direct |27|0.000167039|1.83795e−7|
| Coincident KS fresh-process resume |27|0.000167011|1.83758e−7|

All eight tables pass the frozen1m/0.001m/s physical budgets, mass and heat budgets, canonical-grid completeness/uniqueness, prescribed tick/action histories, root alias and zero clock-origin requirements. Main Cartesian-reference refinement is1.263e−7m and1.016e−10m/s; coincident refinement is6.923e−9m and7.784e−12m/s, comfortably below the frozen fraction of physical budgets.

The main reference has one descending root at2718.3509058762065s. KS root-time errors are approximately5.18e−7s; DP8 errors are a few1e−12s. Root-count, residual and descending-velocity checks pass. The short coincident case has no descending root, as expected.

The main KS heat error is at most5.09e−8J/cm², approximately0.003063of the configured heat budget. Main mass error is at most1.23e−11kg. The independently evaluated thrust propellant-flow oracle also passes. These are numerical law checks for the supplied ideal thrust/Maxwellian models, not validation of real actuator or heating fidelity.

## Resume and boundary contract

Main resumed output now contains all360rows, beginning at0and including the same six marked checkpoint endpoints as the segmented route. Coincident resumed output contains all27canonical rows, beginning at0. The old prototype's301/14resume-only row expectation is not applied to this prospectively changed normal-route contract. Every required global17s row remains present and times are strictly increasing.

The60-row restored main prefix through997.3s matches the independently executed segmented-run prefix bit for bit. Fresh resumed initialization at997.3s restores area15m², held force from988s, active burn, the committed611.25s area action and the correct root-object alias. Coincident resumed initialization at221s restores the already committed impulse and221s control tick, with no repeated committed action.

The interrupted main snapshot at997.3s agrees with Cartesian physics to0.000576805m and6.281e−7m/s; the post-impulse coincident221s snapshot agrees to0.000065961m and7.277e−8m/s. Saved ticks and actions are exact. The coincident committed trace is impulse→control→output before the checkpoint snapshot. Discarded post-checkpoint execution tails are retained separately and are not counted as duplicate committed history.

Across all matched canonical physical times, main segmented/resumed KS versus direct KS differs by at most5.719e−6m and6.193e−9m/s, far inside0.1m/0.0001m/s resume budgets. Coincident fresh resume versus direct differs by at most1.613e−7m and1.756e−10m/s. Final states also pass those budgets.

## Evidence and remaining scope

`lifecycle-current-independent-recomputation-attempt1.json` retains all per-table numerical/history checks. `current-lifecycle-extra-checks-attempt1.json` passes52additional aggregate, root-count, allowed-extra-epoch, snapshot-physics and restored-prefix checks. Portable exports, logs, source-stamp checks and raw/script SHA256 catalogs are retained in this directory. The new numerical values happen to be extremely close to the older frozen-physics experiment; their new raw paths, hashes, source stamps and changed full-prefix output layouts establish their current-source provenance.

No production code or timing work was changed by this review. These results support the tested direct/segmented/fresh-resume lifecycle of the optional route, including the specified coincident priority and observational continuous root. They do not establish arbitrary callback semantics, multi-spacecraft operation, other physical models, unsupported configuration, broad event accuracy or a speed advantage. Performance and the overall upstream-delivery decision remain pending the parent review.
