const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const WRAPPER_TO_CANONICAL = Dict(
    "test/performance_runtime_analysis.jl" => "benchmarks/studies/performance_runtime_analysis.jl",
    "test/performance_smart_parallel_ladder.jl" => "benchmarks/studies/performance_smart_parallel_ladder.jl",
    "test/performance_smart_parallel_ladder_cross_machine.jl" => "benchmarks/studies/performance_smart_parallel_ladder_cross_machine.jl",
    "test/performance_split_imex_compare.jl" => "benchmarks/studies/performance_split_imex_compare.jl",
    "test/performance_static_vs_parallel.jl" => "benchmarks/studies/performance_static_vs_parallel.jl",
    "test/performance_effector_reduction_microbench.jl" => "benchmarks/studies/performance_effector_reduction_microbench.jl",
    "test/performance_paper_pipeline.jl" => "benchmarks/scripts/performance_paper_pipeline.jl",
    "test/gram_interpolation_vs_point_to_point_analysis.jl" => "benchmarks/studies/gram_interpolation_vs_point_to_point_analysis.jl",
    "test/gram_single_call_vs_point_to_point_analysis.jl" => "benchmarks/studies/gram_single_call_vs_point_to_point_analysis.jl",
    "test/gram_real_sim_runtime_compare.jl" => "benchmarks/studies/gram_real_sim_runtime_compare.jl",
    "test/gram_real_sim_surrogate_matrix.jl" => "benchmarks/studies/gram_real_sim_surrogate_matrix.jl",
    "test/gram_real_sim_surrogate_decision_table.jl" => "benchmarks/studies/gram_real_sim_surrogate_decision_table.jl",
    "test/gram_offline_db_accuracy_study.jl" => "benchmarks/studies/gram_offline_db_accuracy_study.jl",
    "test/gram_planet_rho_altitude_sweep.jl" => "benchmarks/studies/gram_planet_rho_altitude_sweep.jl",
    "test/telemetry_hybrid_tuner.jl" => "benchmarks/studies/telemetry_hybrid_tuner.jl",
    "test/telemetry_odyssey_tuner.jl" => "benchmarks/studies/telemetry_odyssey_tuner.jl",
    "test/telemetry_orbit_accuracy_study.jl" => "benchmarks/studies/telemetry_orbit_accuracy_study.jl",
    "test/telemetry_orbit_accuracy_plots.jl" => "benchmarks/studies/telemetry_orbit_accuracy_plots.jl"
)

for (wrapper_rel, canonical_rel) in WRAPPER_TO_CANONICAL
    wrapper_path = joinpath(REPO_ROOT, wrapper_rel)
    canonical_path = joinpath(REPO_ROOT, canonical_rel)

    isfile(wrapper_path) || error("Wrapper missing: $wrapper_rel")
    isfile(canonical_path) || error("Canonical benchmark missing: $canonical_rel")

    wrapper_src = read(wrapper_path, String)
    canonical_basename = basename(canonical_rel)
    occursin(canonical_basename, wrapper_src) || error("Wrapper does not point to canonical benchmark: $wrapper_rel -> $canonical_rel")

    n_lines = count(==('\n'), wrapper_src) + 1
    n_lines <= 12 || error("Wrapper is not thin (expected <=12 lines): $wrapper_rel")
end

println("benchmark_wrapper_parity_gate_ok")
