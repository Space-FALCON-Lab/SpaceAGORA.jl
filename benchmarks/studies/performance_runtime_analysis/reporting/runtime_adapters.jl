# Bind the isolated reporting module to the benchmark runtime.
function _runtime_report_context(spec::ProfileSpec)
    return (hardware=_runtime_hardware_snapshot(),
        solver_mode_default=_perf_default_solver_mode(spec.name),
        solver_mode_effective=_perf_solver_mode_env(spec.name),
        split_rollout_enabled=RuntimeReports._split_rollout_enabled(),
        split_rollout_enforce=RuntimeReports._split_rollout_enforce(),
        split_rollout_case_names=RuntimeReports._split_rollout_case_names(),
        split_rollout_solver_variants=RuntimeReports._split_rollout_solver_variants(),
        multirate_rollout_enabled=RuntimeReports._multirate_rollout_enabled(),
        multirate_rollout_enforce=RuntimeReports._multirate_rollout_enforce(),
        multirate_rollout_case_names=RuntimeReports._multirate_rollout_case_names(),
        multirate_rollout_max_slowdown_ratio=RuntimeReports._multirate_rollout_max_slowdown_ratio(),
    )
end

write_report(path, spec::ProfileSpec, args...; kwargs...) =
    RuntimeReports.write_report(path, spec, args...;
        context=_runtime_report_context(spec), baseline_scenario=PERF_BASELINE_SCENARIO, kwargs...)

evaluate_split_rollout_gate(spec::ProfileSpec, cases::Vector{BenchmarkCase}, outdir::String) =
    RuntimeReports.evaluate_split_rollout_gate(spec, cases, outdir;
        run_warmup, run_solution=_run_split_gate_solution)

evaluate_multirate_rollout_gate(spec::ProfileSpec, cases::Vector{BenchmarkCase}, outdir::String) =
    RuntimeReports.evaluate_multirate_rollout_gate(spec, cases, outdir;
        run_warmup, run_solution=_run_split_gate_solution)
