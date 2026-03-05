# SpaceAGORA Shim Window Manifest

Date: 2026-03-05
Policy: keep wrappers for one release cycle after Wave 1 merge, then retire in Wave 2.

## Active Shim Inventory (Wave 1)
| Shim path | Canonical target | Retirement wave |
|---|---|---|
| `src/simulation/execution/run_simulation.jl` | `src/simulation/engine/*.jl` split files | Wave 2 |
| `src/simulation_model/callbacks.jl` | `src/simulation/callbacks/callbacks.jl` | Wave 2 |
| `src/parallel/ParallelProfiles.jl` | `src/parallel/routing/parallel_profiles.jl` | Wave 2 |
| `test/performance_runtime_analysis.jl` | `benchmarks/studies/performance_runtime_analysis.jl` | Wave 2 |
| `test/performance_smart_parallel_ladder.jl` | `benchmarks/studies/performance_smart_parallel_ladder.jl` | Wave 2 |
| `test/performance_smart_parallel_ladder_cross_machine.jl` | `benchmarks/studies/performance_smart_parallel_ladder_cross_machine.jl` | Wave 2 |
| `test/performance_split_imex_compare.jl` | `benchmarks/studies/performance_split_imex_compare.jl` | Wave 2 |
| `test/performance_static_vs_parallel.jl` | `benchmarks/studies/performance_static_vs_parallel.jl` | Wave 2 |
| `test/performance_effector_reduction_microbench.jl` | `benchmarks/studies/performance_effector_reduction_microbench.jl` | Wave 2 |
| `test/performance_paper_pipeline.jl` | `benchmarks/scripts/performance_paper_pipeline.jl` | Wave 2 |
| `test/gram_interpolation_vs_point_to_point_analysis.jl` | `benchmarks/studies/gram_interpolation_vs_point_to_point_analysis.jl` | Wave 2 |
| `test/gram_single_call_vs_point_to_point_analysis.jl` | `benchmarks/studies/gram_single_call_vs_point_to_point_analysis.jl` | Wave 2 |
| `test/gram_real_sim_runtime_compare.jl` | `benchmarks/studies/gram_real_sim_runtime_compare.jl` | Wave 2 |
| `test/gram_real_sim_surrogate_matrix.jl` | `benchmarks/studies/gram_real_sim_surrogate_matrix.jl` | Wave 2 |
| `test/gram_real_sim_surrogate_decision_table.jl` | `benchmarks/studies/gram_real_sim_surrogate_decision_table.jl` | Wave 2 |
| `test/gram_offline_db_accuracy_study.jl` | `benchmarks/studies/gram_offline_db_accuracy_study.jl` | Wave 2 |
| `test/gram_planet_rho_altitude_sweep.jl` | `benchmarks/studies/gram_planet_rho_altitude_sweep.jl` | Wave 2 |
| `test/telemetry_hybrid_tuner.jl` | `benchmarks/studies/telemetry_hybrid_tuner.jl` | Wave 2 |
| `test/telemetry_odyssey_tuner.jl` | `benchmarks/studies/telemetry_odyssey_tuner.jl` | Wave 2 |
| `test/telemetry_orbit_accuracy_study.jl` | `benchmarks/studies/telemetry_orbit_accuracy_study.jl` | Wave 2 |
| `test/telemetry_orbit_accuracy_plots.jl` | `benchmarks/studies/telemetry_orbit_accuracy_plots.jl` | Wave 2 |

## Retirement Preconditions
1. One tagged release completed after Wave 1 merge.
2. Canonical benchmark/script paths are used in CI workflows.
3. Contract gates are updated to canonical-only references.
4. No external downstream dependency on wrapper paths.

## Retirement Execution Checklist
1. Remove forwarding wrapper file.
2. Update docs and CI references.
3. Re-run full mandatory gates and nightly/study smoke.
