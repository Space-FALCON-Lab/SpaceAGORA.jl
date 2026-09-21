using Test
using DataFrames

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "benchmarks", "studies", "performance_runtime_analysis", "main.jl"))

spec = ProfileSpec(
    name="quick",
    repeats=1,
    warmup=0,
    max_attempts=1,
    mission_short_s=120.0,
    mission_long_s=600.0,
    montecarlo_samples=1,
    montecarlo_mission_s=120.0
)

raw_df = DataFrame(
    category=["baseline", "orientation"],
    scenario=[PERF_BASELINE_SCENARIO, "single_orientation_aero"],
    description=["Synthetic baseline", "Synthetic orientation"],
    satellites=[1, 1],
    orientation=[false, true],
    mission_time_s=[120.0, 120.0],
    outer_route=["serial", "serial"],
    outer_threads_safe=[true, true],
    density_parallel_mode=["serial", "serial"],
    control_parallel_mode=["serial", "serial"],
    multibody_parallel_mode=["serial", "serial"],
    dt_max_orbit_s=[10.0, 10.0],
    dynamic_effectors=["InverseSquaredGravityModel", "InverseSquaredGravityModel+Aero"],
    control_effectors=["none", "none"],
    solve_success=[true, true],
    total_time_s=[1.0, 2.0],
    solve_time_s=[0.75, 1.2],
    copy_time_s=[0.25, 0.8],
    copy_compile_time_s=[0.0, 0.0],
    solve_compile_time_s=[0.0, 0.0],
    copy_gctime_s=[0.0, 0.0],
    solve_gctime_s=[0.0, 0.0],
    total_bytes_mb=[12.0, 18.0],
    copy_bytes_mb=[3.0, 9.0],
    solve_bytes_mb=[9.0, 9.0],
    copy_alloc_calls=[100, 250],
    solve_alloc_calls=[300, 350],
    saved_points=[10, 10],
    accepted_steps=[20, 20],
    rejected_steps=[0, 0],
    solver_mode=["auto", "auto"],
    solver_sequence=["Tsit5", "Tsit5"],
    solver_fallback_used=[false, false],
    solver_fallback_count=[0, 0],
    solver_fallback_trigger=[missing, missing],
    policy_threads_enabled_total=[0, 0],
    policy_decisions_total=[0, 0],
    policy_density_threads_enabled=[0, 0],
    policy_control_threads_enabled=[0, 0],
    policy_multibody_threads_enabled=[0, 0],
    policy_other_threads_enabled=[0, 0],
    sim_seconds_per_wall_second=[120.0, 60.0],
    satellite_sim_seconds_per_wall_second=[120.0, 60.0],
    nbody_spkpos_runtime_calls=[0, 0],
    nbody_spkpos_cache_build_calls=[0, 0],
    nbody_spkpos_total_calls=[0, 0],
    srp_spkpos_runtime_calls=[0, 0],
    srp_spkpos_cache_build_calls=[0, 0],
    srp_spkpos_total_calls=[0, 0],
    planet_pxform_runtime_calls=[0, 0],
    planet_pxform_cache_build_calls=[0, 0],
    planet_pxform_total_calls=[0, 0],
    attempt=[1, 1]
)

summary_df = summarize_results(raw_df)
summary_names = Set(Symbol.(names(summary_df)))

for required in (:copy_bytes_mean_mb, :solve_bytes_mean_mb, :copy_time_share, :copy_bytes_share, :copy_alloc_share)
    @test required in summary_names
end

orbit_summary_df = DataFrame(
    scenario=[PERF_BASELINE_SCENARIO],
    category=["baseline"],
    mission_time_multiplier=[1],
    samples_success=[1],
    samples_total=[1],
    mission_time_mean_s=[120.0],
    total_time_mean_s=[1.0],
    total_time_ci95_low_s=[1.0],
    total_time_ci95_high_s=[1.0],
    total_time_p90_s=[1.0],
    time_per_baseline_period_mean_s=[1.0],
    baseline_periods_per_wall_second_mean=[1.0]
)

mktempdir() do tmp
    report_path = joinpath(tmp, "runtime_report.md")
    write_report(
        report_path,
        spec,
        raw_df,
        summary_df,
        orbit_summary_df,
        DataFrame()
    )

    report_text = read(report_path, String)
    @test occursin("## State Isolation Overhead", report_text)
    @test occursin("isolate_state=false", report_text)
end

println("runtime_analysis_copy_overhead_gate_ok")
