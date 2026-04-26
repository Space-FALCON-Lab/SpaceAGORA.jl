const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const RETIRED_WRAPPERS = (
    "test/performance_runtime_analysis.jl",
    "test/performance_smart_parallel_ladder.jl",
    "test/performance_smart_parallel_ladder_cross_machine.jl",
    "test/performance_split_imex_compare.jl",
    "test/performance_static_vs_parallel.jl",
    "test/performance_effector_reduction_microbench.jl",
    "test/performance_paper_pipeline.jl",
    "test/gram_interpolation_vs_point_to_point_analysis.jl",
    "test/gram_single_call_vs_point_to_point_analysis.jl",
    "test/gram_real_sim_runtime_compare.jl",
    "test/gram_real_sim_surrogate_matrix.jl",
    "test/gram_real_sim_surrogate_decision_table.jl",
    "test/gram_offline_db_accuracy_study.jl",
    "test/gram_planet_rho_altitude_sweep.jl",
    "test/telemetry_hybrid_tuner.jl",
    "test/telemetry_odyssey_tuner.jl",
    "test/telemetry_orbit_accuracy_study.jl",
    "test/telemetry_orbit_accuracy_plots.jl"
)

const CANONICAL_BENCHMARKS = (
    "benchmarks/studies/performance_runtime_analysis.jl",
    "benchmarks/studies/performance_smart_parallel_ladder.jl",
    "benchmarks/studies/performance_smart_parallel_ladder_cross_machine.jl",
    "benchmarks/studies/performance_split_imex_compare.jl",
    "benchmarks/studies/performance_static_vs_parallel.jl",
    "benchmarks/studies/performance_effector_reduction_microbench.jl",
    "benchmarks/scripts/performance_paper_pipeline.jl",
    "benchmarks/studies/gram_interpolation_vs_point_to_point_analysis.jl",
    "benchmarks/studies/gram_single_call_vs_point_to_point_analysis.jl",
    "benchmarks/studies/gram_real_sim_runtime_compare.jl",
    "benchmarks/studies/gram_real_sim_surrogate_matrix.jl",
    "benchmarks/studies/gram_real_sim_surrogate_decision_table.jl",
    "benchmarks/studies/gram_offline_db_accuracy_study.jl",
    "benchmarks/studies/gram_planet_rho_altitude_sweep.jl",
    "benchmarks/studies/telemetry_hybrid_tuner.jl",
    "benchmarks/studies/telemetry_odyssey_tuner.jl",
    "benchmarks/studies/telemetry_orbit_accuracy_study.jl",
    "benchmarks/studies/telemetry_orbit_accuracy_plots.jl"
)

for rel in RETIRED_WRAPPERS
    isfile(joinpath(REPO_ROOT, rel)) && error("Wave 2 wrapper still present: $rel")
end

for rel in CANONICAL_BENCHMARKS
    isfile(joinpath(REPO_ROOT, rel)) || error("Canonical benchmark missing: $rel")
end

canonical_launcher = joinpath(REPO_ROOT, "benchmarks", "studies", "performance_runtime_analysis.jl")
launcher_src = read(canonical_launcher, String)
(
    occursin(joinpath("performance_runtime_analysis", "main.jl"), launcher_src) ||
    occursin("\"performance_runtime_analysis\", \"main.jl\"", launcher_src)
) || error("Runtime-analysis canonical launcher does not include split main.jl")
launcher_lines = count(==('\n'), launcher_src) + 1
launcher_lines <= 10 || error("Runtime-analysis canonical launcher is not thin (expected <=10 lines)")

println("benchmark_wrapper_parity_gate_ok")
