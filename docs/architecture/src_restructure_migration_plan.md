# SpaceAGORA Architecture Cleanup and Restructure Status

Date: 2026-03-05
Branch baseline: `refactor/src-restructure`

## Locked Decisions
1. `network` lives under `mission/constellation/network`.
2. No `vehicle/payloads` folder in this cycle; laser remains under `vehicle/actuators/laser_terminal`.
3. `gnc/estimation` is a sibling of `gnc/navigation`.
4. Compatibility shims stay warning-free for one release cycle.
5. Telemetry thresholds are blocking for this cleanup track.

## Current Architecture Targets
1. Canonical simulation entrypoint is `SimulationEngine.run_simulation`.
2. `SpaceAGORA.run_simulation` is a stable forwarding root export.
3. Verification consumes simulation engine APIs and does not own engine source.
4. Typed runtime configuration is centered on:
   - `ParallelConfig`
   - `SolverConfig`
   - `RuntimePolicyConfig`
   - `ArtifactConfig`
   - `SimulationEngineConfig`
5. ENV reads are adapter-boundary behavior for engine internals.

## Migration Status Table
| Area | Status | Notes |
|---|---|---|
| `src/` subsystem restructure | Implemented | Core/environment/vehicle/dynamics/gnc/mission/simulation/parallel/io/analysis layout in place with compatibility wrappers |
| Engine ownership inversion | Implemented | `SpaceAGORA.run_simulation` forwards to `SimulationEngine.run_simulation` |
| Typed runtime config boundary | Implemented | Engine config types and adapter constructors available |
| `run_simulation` split by responsibility | Implemented | Canonical split under `src/simulation/engine/` |
| Callback family split | Implemented | Canonical callback files under `src/simulation/callbacks/`; legacy path shimmed |
| Parallel routing split | Implemented | Split files under `src/parallel/routing/`; aggregator retained |
| Study/workflow migration out of `test/` | Implemented (compat window active) | Canonical study/scripts under `benchmarks/`; `test/` wrappers retained |
| `performance_runtime_analysis` monolith split | Implemented | Folderized split under `benchmarks/studies/performance_runtime_analysis/` |
| Naming contract gate | Implemented | Enforced by `test/ci_naming_contract_gate.jl` |
| Shim retirement | Pending (Wave 2) | Deferred until next release cycle per policy |

## CI and Quality Gates (Required)
1. `julia --project=.AGORA test/runtests.jl`
2. `julia --project=.AGORA --code-coverage=user test/runtests.jl`
3. `julia --project=.AGORA test/ci_coverage_quality_gate.jl`
4. `julia --project=.AGORA test/ci_p1_findings_gate.jl`
5. `julia --startup-file=no --depwarn=error --project=.AGORA test/ci_clean_depot_smoke.jl`
6. `julia --startup-file=no --depwarn=error --project=.AGORA test/ci_threaded_smoke.jl`
7. `julia --startup-file=no --project=.AGORA test/ci_examples_regression.jl`
8. `julia --startup-file=no --depwarn=error --project=.AGORA test/telemetry_orbit_accuracy_study.jl quick --enforce=true`

Additional contract gates:
1. `julia --project=.AGORA test/ci_architecture_contract_gate.jl`
2. `julia --project=.AGORA test/ci_typed_config_equivalence_gate.jl`
3. `julia --project=.AGORA test/ci_benchmark_wrapper_parity_gate.jl`
4. `julia --project=.AGORA test/ci_naming_contract_gate.jl`

## Compatibility Window Policy
1. Shim files are include-forwarders only.
2. Shims must not introduce new logic.
3. No deprecation warnings during shim window.
4. Shim retirement trigger: one tagged release after Wave 1 merge.
5. Shim inventory and removal candidates are tracked in:
   - `docs/architecture/shim_window_manifest.md`

## Wave 2 (Post-Window) Retirement Scope
1. Remove old simulation execution aliases.
2. Remove old callback path aliases.
3. Remove `test/` study wrappers once CI points to canonical benchmark paths.
4. Remove temporary naming/contract allowlists tied to shim paths.

## Assumptions
1. No intentional physics/model behavior change in this migration track.
2. Public root entrypoint signatures stay stable through compatibility window.
3. Telemetry threshold policy remains enforced and blocking.
