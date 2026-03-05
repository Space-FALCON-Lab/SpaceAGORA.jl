# SpaceAGORA Shim Window Manifest

Date: 2026-03-05
Policy outcome: Wave 2 shim retirement completed.

## Retired Wrapper Inventory
| Retired path | Canonical replacement |
|---|---|
| `src/simulation/execution/run_simulation.jl` | `src/simulation/engine/execution.jl` + `src/simulation/engine/public_api.jl` |
| `src/simulation_model/callbacks.jl` | `src/simulation/callbacks/callbacks.jl` |
| `test/performance_runtime_analysis.jl` | `benchmarks/studies/performance_runtime_analysis.jl` |
| `test/performance_smart_parallel_ladder.jl` | `benchmarks/studies/performance_smart_parallel_ladder.jl` |
| `test/performance_smart_parallel_ladder_cross_machine.jl` | `benchmarks/studies/performance_smart_parallel_ladder_cross_machine.jl` |
| `test/performance_split_imex_compare.jl` | `benchmarks/studies/performance_split_imex_compare.jl` |
| `test/performance_static_vs_parallel.jl` | `benchmarks/studies/performance_static_vs_parallel.jl` |
| `test/performance_effector_reduction_microbench.jl` | `benchmarks/studies/performance_effector_reduction_microbench.jl` |
| `test/performance_paper_pipeline.jl` | `benchmarks/scripts/performance_paper_pipeline.jl` |
| `test/gram_interpolation_vs_point_to_point_analysis.jl` | `benchmarks/studies/gram_interpolation_vs_point_to_point_analysis.jl` |
| `test/gram_single_call_vs_point_to_point_analysis.jl` | `benchmarks/studies/gram_single_call_vs_point_to_point_analysis.jl` |
| `test/gram_real_sim_runtime_compare.jl` | `benchmarks/studies/gram_real_sim_runtime_compare.jl` |
| `test/gram_real_sim_surrogate_matrix.jl` | `benchmarks/studies/gram_real_sim_surrogate_matrix.jl` |
| `test/gram_real_sim_surrogate_decision_table.jl` | `benchmarks/studies/gram_real_sim_surrogate_decision_table.jl` |
| `test/gram_offline_db_accuracy_study.jl` | `benchmarks/studies/gram_offline_db_accuracy_study.jl` |
| `test/gram_planet_rho_altitude_sweep.jl` | `benchmarks/studies/gram_planet_rho_altitude_sweep.jl` |
| `test/telemetry_hybrid_tuner.jl` | `benchmarks/studies/telemetry_hybrid_tuner.jl` |
| `test/telemetry_odyssey_tuner.jl` | `benchmarks/studies/telemetry_odyssey_tuner.jl` |
| `test/telemetry_orbit_accuracy_study.jl` | `benchmarks/studies/telemetry_orbit_accuracy_study.jl` |
| `test/telemetry_orbit_accuracy_plots.jl` | `benchmarks/studies/telemetry_orbit_accuracy_plots.jl` |

## Enforcement
1. `test/ci_architecture_contract_gate.jl` fails if retired wrappers reappear.
2. `test/ci_benchmark_wrapper_parity_gate.jl` fails if retired wrappers reappear.
