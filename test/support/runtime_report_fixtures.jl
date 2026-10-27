# Synthetic supplied attempt rows. This is test data, not a second accounting implementation.
function fixture_report_context(hardware)
    return (hardware=hardware, solver_mode_default="fixture-default", solver_mode_effective="fixture-effective",
        split_rollout_enabled=false, split_rollout_enforce=false, split_rollout_case_names=String[],
        split_rollout_solver_variants=String[], multirate_rollout_enabled=false,
        multirate_rollout_enforce=false, multirate_rollout_case_names=String[],
        multirate_rollout_max_slowdown_ratio=1.0)
end
summarize_results(raw) = RuntimeReports.summarize_results(raw; baseline_scenario=PERF_BASELINE_SCENARIO)
write_report(args...; kwargs...) = RuntimeReports.write_report(args...;
    context=fixture_report_context(FALLBACK), baseline_scenario=PERF_BASELINE_SCENARIO, kwargs...)

function raw_row(; kwargs...)
    return merge((category="baseline", scenario=PERF_BASELINE_SCENARIO,
        description="Synthetic reporting fixture", satellites=1, orientation=false,
        mission_time_s=120.0, outer_route="serial", outer_threads_safe=true,
        density_parallel_mode="serial", control_parallel_mode="serial",
        multibody_parallel_mode="serial", dt_max_orbit_s=10.0,
        dynamic_effectors="fixture", control_effectors="none", solve_success=true,
        total_time_s=2.0, solve_time_s=1.5, copy_time_s=0.5,
        copy_compile_time_s=0.0, solve_compile_time_s=0.0,
        copy_gctime_s=0.0, solve_gctime_s=0.0,
        total_bytes_mb=12.0, copy_bytes_mb=3.0, solve_bytes_mb=9.0,
        copy_alloc_calls=100, solve_alloc_calls=300, saved_points=10,
        accepted_steps=20, rejected_steps=0, solver_mode="auto",
        solver_sequence="Tsit5", solver_fallback_used=false,
        solver_fallback_count=0, solver_fallback_trigger=missing,
        policy_threads_enabled_total=0, policy_decisions_total=0,
        policy_density_threads_enabled=0, policy_control_threads_enabled=0,
        policy_multibody_threads_enabled=0, policy_other_threads_enabled=0,
        sim_seconds_per_wall_second=60.0, satellite_sim_seconds_per_wall_second=60.0,
        nbody_spkpos_runtime_calls=0, nbody_spkpos_cache_build_calls=0,
        nbody_spkpos_total_calls=0, srp_spkpos_runtime_calls=0,
        srp_spkpos_cache_build_calls=0, srp_spkpos_total_calls=0,
        planet_pxform_runtime_calls=0, planet_pxform_cache_build_calls=0,
        planet_pxform_total_calls=0, attempt=1), (; kwargs...))
end
function render_report(raw, summary=summarize_results(raw), entry=nothing; kwargs...)
    # Writing must not reorder or alter any supplied table.
    tables = Any[raw, summary, entry, values(kwargs)...]
    before = deepcopy(tables)
    text = mktempdir() do dir
        path = joinpath(dir, "report.md")
        write_report(path, SPEC, raw, summary, DataFrame(), entry; kwargs...)
        read(path, String)
    end
    @test isequal(tables, before)
    return text
end
