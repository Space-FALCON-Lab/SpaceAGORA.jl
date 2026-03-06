# SpaceAGORA Architecture Cleanup and Restructure Status

Date: 2026-03-05
Branch baseline: `refactor/src-restructure`

## Locked Decisions
1. `network` lives under `mission/constellation/network`.
2. No `vehicle/payloads` folder in this cycle; laser remains under `vehicle/actuators/laser_terminal`.
3. `gnc/estimation` is a sibling of `gnc/navigation`.
4. Telemetry thresholds are blocking for this cleanup track.

## Current Architecture Targets
1. Canonical simulation entrypoint is `SimulationEngine.run_simulation`.
2. `SpaceAGORA.run_simulation` is the stable root forwarder.
3. Verification consumes simulation engine APIs and does not own engine source.
4. Typed runtime configuration is centered on:
   - `ParallelConfig`
   - `SolverConfig`
   - `RuntimePolicyConfig`
   - `ArtifactConfig`
   - `SimulationEngineConfig`
5. ENV reads are adapter-boundary behavior for engine internals.

## Completeness Contract
1. Canonical ownership must not regress into silent stubs or legacy wrapper ownership.
2. Scaffold folders remain in place, but scaffold files must be typed interfaces with explicit `Not implemented:` failures.
3. Canonical include-chain discipline is enforced for core engine/callback/parallel and new scaffold roots.
4. See `docs/architecture/src_completeness_contract.md` for categories and allowlisted aggregators.

## Migration Status Table
| Area | Status | Notes |
|---|---|---|
| `src/` subsystem restructure | Implemented | Core/environment/vehicle/dynamics/gnc/mission/simulation/parallel/io/analysis layout in place |
| Engine ownership inversion | Implemented | `SpaceAGORA.run_simulation` forwards to `SimulationEngine.run_simulation` |
| Typed runtime config boundary | Implemented | Engine config types and adapter constructors available |
| `run_simulation` split by responsibility | Implemented | Canonical split under `src/simulation/engine/` |
| Callback family split | Implemented | Canonical callback files under `src/simulation/callbacks/` |
| Parallel routing split | Implemented | Split files under `src/parallel/routing/`; aggregator retained |
| Study/workflow migration out of `test/` | Implemented | Canonical study/scripts under `benchmarks/`; wrappers removed |
| `performance_runtime_analysis` monolith split | Implemented | Folderized split under `benchmarks/studies/performance_runtime_analysis/` |
| Naming contract gate | Implemented | Enforced by `test/ci_naming_contract_gate.jl` |
| Shim retirement | Implemented | Wave 2 removal complete |
| `src/` completeness contract | Implemented | Enforced by dedicated completeness/inert-term/scaffold gates |
| Legacy include-guard cleanup (`__legacy_*`) | Implemented | Canonical GNC/solver/report files now use shared include helpers with no `__legacy_*` markers |

## CI and Quality Gates (Required)
1. `julia --project=.AGORA test/runtests.jl`
2. `julia --project=.AGORA --code-coverage=user test/runtests.jl`
3. `julia --project=.AGORA test/ci_coverage_quality_gate.jl`
4. `julia --project=.AGORA test/ci_p1_findings_gate.jl`
5. `julia --startup-file=no --depwarn=error --project=.AGORA test/ci_clean_depot_smoke.jl`
6. `JULIA_NUM_THREADS=2 julia --startup-file=no --depwarn=error --project=.AGORA test/ci_threaded_smoke.jl`
7. `julia --startup-file=no --project=.AGORA test/ci_examples_regression.jl`
8. `julia --startup-file=no --depwarn=error --project=.AGORA benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true`

Additional contract gates:
1. `julia --project=.AGORA test/ci_architecture_contract_gate.jl`
2. `julia --project=.AGORA test/ci_typed_config_equivalence_gate.jl`
3. `julia --project=.AGORA test/ci_benchmark_wrapper_parity_gate.jl`
4. `julia --project=.AGORA test/ci_naming_contract_gate.jl`
5. `julia --project=.AGORA test/ci_no_legacy_ownership_gate.jl`
6. `julia --project=.AGORA test/ci_no_artifact_files_gate.jl`
7. `julia --project=.AGORA test/ci_canonical_path_contract_gate.jl`
8. `julia --project=.AGORA test/ci_src_completeness_contract_gate.jl`
9. `julia --project=.AGORA test/ci_no_legacy_include_chains_gate.jl`
10. `julia --project=.AGORA test/ci_no_inert_force_terms_gate.jl`
11. `julia --project=.AGORA test/ci_scaffold_interface_gate.jl`

## Wave 2 Retirement Summary
Removed:
1. legacy simulation execution wrapper file
2. benchmark/study wrappers in `test/`

Contracts tightened:
1. architecture gate fails if retired wrappers reappear
2. benchmark gate fails if retired wrappers reappear
3. docs and workflows point to canonical benchmark paths only

## Assumptions
1. No intentional physics/model behavior change in this migration track.
2. Public root entrypoint signatures remain stable.
3. Telemetry threshold policy remains enforced and blocking.
