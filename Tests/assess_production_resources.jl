include(joinpath(@__DIR__, "..", "Simulation", "Code", "StudyConfig.jl"))
using .StudyConfig
using TOML

cases = (StudyCase(1050, 1, 300, :gve_sma), StudyCase(1050, 1, 50, :gve_inc),
    StudyCase(850, 1, 50, :gve_sma))
infos = [TOML.parsefile(joinpath(case_directory(DEFAULT_OUTPUT, c), "run_info.toml")) for c in cases]
rates = [info["output_bytes"] / ((c.helpers + 1) * mission_duration(c)) for (info, c) in zip(infos, cases)]
grid_scale = sum((c.helpers + 1) * mission_duration(c) for c in grid_cases())
low, high = extrema(rates) .* grid_scale
free = free_disk_bytes(DEFAULT_OUTPUT)
existing = sum(info["output_bytes"] for info in infos)
reserve = 5 * 1024^3
report = Dict{String,Any}("benchmark_spacecraft" => 301,
    "benchmark_orbits" => 100, "benchmark_peak_rss_bytes" => infos[1]["peak_rss_bytes"],
    "benchmark_output_bytes" => infos[1]["output_bytes"],
    "benchmark_simulation_recording_writing_s" => infos[1]["wall_elapsed_s"],
    "estimated_grid_bytes_low" => low, "estimated_grid_bytes_high" => high,
    "already_saved_full_pilot_bytes" => existing, "available_disk_bytes" => free,
    "disk_reserve_bytes" => reserve,
    "estimated_additional_free_space_needed_bytes" => max(0.0, high - existing + reserve - free),
    "estimated_grid_fits_current_disk" => high - existing + reserve <= free,
    "estimate_method" => "range of observed compressed bytes per spacecraft-second across three full-duration pilots; not a guaranteed bound")
atomic_toml(joinpath(@__DIR__, "results", "production_pipeline_compact", "resource_assessment.toml"), report)
TOML.print(stdout, report; sorted=true)
