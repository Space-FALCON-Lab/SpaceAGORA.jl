include(joinpath(@__DIR__, "..", "..", "Simulation_absolute_RN", "Code", "StudyConfig.jl"))
using .StudyConfig
using Dates, Printf, TOML

module SnapshotPreservation
include("create_snapshot.jl")
end

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RESULTS = joinpath(@__DIR__, "results")
const OLD_OUTPUT = joinpath(ROOT, "Simulation", "Output")
const PILOTS = ["sma950", "sma1050", "inc950", "inc1050", "zero_sma", "coplanar_inc",
    "large800", "inc1050_control", "baseline_sma800", "new_sma800",
    "baseline_inc800", "new_inc800", "baseline_zero"]
const COMPARISONS = ["inc1050_inc1050_control", "baseline_sma800_new_sma800",
    "baseline_inc800_new_inc800", "baseline_zero_zero_sma"]

require(condition, message) = condition || error(message)

function historical_case(c)
    path = joinpath(case_directory(OLD_OUTPUT, c), "run_info.toml")
    if isfile(path)
        info = TOML.parsefile(path)
        require(info["status"] == "verified", "Historical baseline is not verified: $path")
        trajectory_bytes = filesize(joinpath(dirname(path), "trajectory.feather"))
        return (wall=info["wall_elapsed_s"], bytes=Float64(info["output_bytes"]),
            event_bytes=Float64(info["output_bytes"] - trajectory_bytes), interpolated=false)
    end
    require(c.helpers == 700, "Missing historical baseline: $path")
    left = historical_case(StudyCase(c.altitude_km, c.inclination_deg, 600, c.policy))
    right = historical_case(StudyCase(c.altitude_km, c.inclination_deg, 800, c.policy))
    (wall=(left.wall + right.wall) / 2, bytes=(left.bytes + right.bytes) / 2,
        event_bytes=(left.event_bytes + right.event_bytes) / 2, interpolated=true)
end

function check_snapshot()
    snapshot = SnapshotPreservation.SNAPSHOT
    manifest = TOML.parsefile(joinpath(snapshot, "manifest.toml"))
    require(SnapshotPreservation.selected_hashes(joinpath(snapshot, "project")) == manifest["sha256"],
        "The frozen pre-change source snapshot was modified")
    archive_hash = first(split(read(joinpath(snapshot, "sources.tar.gz.sha256"), String)))
    require(filehash(joinpath(snapshot, "sources.tar.gz")) == archive_hash, "Snapshot archive checksum failed")
    expected = Set(joinpath(SnapshotPreservation.LIBRARY_NAME, path) for path in (
        "src/orchestration/simulation/engine/dynamics_rhs.jl",
        "src/orchestration/simulation/engine/execution.jl",
        "src/orchestration/simulation/callbacks/save_fields.jl"))
    changed = Set(path for (path, hash) in manifest["sha256"] if filehash(joinpath(ROOT, path)) != hash)
    require(changed == expected, "Unexpected changes outside the three intended library files: $changed")
    Dict("directory" => snapshot, "archive_sha256" => archive_hash,
        "file_count" => length(manifest["sha256"]), "changed_library_files" => sort!(collect(changed)),
        "original_runner_unchanged" => true)
end

function check_old_700()
    status = TOML.parsefile(joinpath(OLD_OUTPUT, "_batch_helpers700", "status.toml"))
    require(status["status"] == "stopped_by_user" && status["verified_cases"] == 4,
        "Old 700-helper campaign status changed; recheck processes before proceeding")
    require(isfile(joinpath(OLD_OUTPUT, "_batch_helpers700", "STOP_AFTER_CASE")),
        "Old campaign stop flag is no longer present")
    require(status["active_case_result"] == "deleted_unfinished" && !ispath(status["active_case"]),
        "The discarded old case unexpectedly reappeared")
    cases = String[]
    checked_bytes = 0
    for (directory, _, files) in walkdir(joinpath(OLD_OUTPUT, "paper_plot_mode", "h1000km_t850km"))
        occursin("/N700/", directory) && "run_info.toml" in files || continue
        info = TOML.parsefile(joinpath(directory, "run_info.toml"))
        require(info["status"] == "verified", "An old N700 case is no longer verified")
        for (name, checksum) in info["checksums"]
            path = joinpath(directory, name)
            require(filehash(path) == checksum, "Old N700 output checksum mismatch: $path")
            checked_bytes += filesize(path)
        end
        push!(cases, directory)
    end
    require(length(cases) == 4, "Expected four preserved old N700 cases")
    Dict("passed" => true, "verified_cases" => 4, "directories" => sort!(cases),
        "checksummed_bytes" => checked_bytes, "discarded_case_still_absent" => true)
end

function main()
    result_path = joinpath(RESULTS, "preflight.toml")
    atomic_toml(result_path, Dict("passed" => false, "status" => "checking"))
    signature = source_signature(source_hashes())
    fixtures = TOML.parsefile(joinpath(RESULTS, "fixtures.toml"))
    require(fixtures["passed"] && fixtures["source_signature"] == signature, "Tests are missing or stale")
    assertion_count = 0
    for line in eachline(joinpath(RESULTS, "fixtures_final.log"))
        match_result = match(r"\|\s+(\d+)\s+(\d+)\s+[0-9.]+s$", line)
        match_result === nothing && continue
        require(match_result[1] == match_result[2], "A fixture assertion failed")
        assertion_count += parse(Int, match_result[1])
    end
    require(assertion_count > 0, "No completed test summaries found")
    coverage = TOML.parsefile(joinpath(RESULTS, "sign_coverage.toml"))
    require(coverage["passed"], "R/N sign coverage failed")
    reports = Dict(label => TOML.parsefile(joinpath(RESULTS, label * ".toml")) for label in PILOTS)
    max_rss, max_net, max_absolute, max_quadrature = 0, 0.0, 0.0, 0.0
    for (label, report) in reports
        require(report["fresh_file_verified"], "Pilot was not freshly verified: $label")
        directory = report["output_directory"]
        require(startswith(directory, joinpath(@__DIR__, "PilotOutput") * "/"),
            "Pilot output escaped the required directory: $directory")
        info = TOML.parsefile(joinpath(directory, "run_info.toml"))
        require(info["status"] == "verified" && info["trajectory_rows"] == 1001,
            "Incomplete pilot: $label")
        require(info["julia_threads"] == 1, "Pilot was not single-threaded: $label")
        for (name, checksum) in info["checksums"]
            require(filehash(joinpath(directory, name)) == checksum, "Pilot checksum changed: $label/$name")
        end
        verification = TOML.parsefile(joinpath(directory, "verification.toml"))
        require(verification["passed"], "Pilot verification failed: $label")
        if report["kind"] == "new"
            require(report["source_signature"] == signature &&
                verification["verification_source_signature"] == signature, "Pilot code is stale: $label")
        end
        max_rss = max(max_rss, report["peak_rss_bytes"])
        if report["mode"] == "diagnostic"
            require(report["passed"] && report["fresh_reference_verified"], "Independent check failed: $label")
            max_net = max(max_net, maximum(report["net_error_mps"]))
            max_absolute = max(max_absolute, maximum(report["absolute_error_mps"]))
            max_quadrature = max(max_quadrature, maximum(report["summed_quadrature_difference_mps"]))
        end
    end
    comparisons = Dict(name => TOML.parsefile(joinpath(RESULTS, "comparison_" * name * ".toml"))
        for name in COMPARISONS)
    require(all(report["passed"] for report in values(comparisons)), "Physical comparison failed")
    snapshot = check_snapshot()
    preserved_700 = check_old_700()

    cases = grid_cases()
    require(length(cases) == 270 && all(c.orbits == 100 for c in cases), "Wrong production grid")
    baselines = historical_case.(cases)
    ratios = Dict{String,Float64}()
    projections = Dict{String,Any}()
    for (policy, suffix) in (("gve_sma", "sma800"), ("gve_inc", "inc800"))
        old, new = reports["baseline_" * suffix], reports["new_" * suffix]
        require(old["mode"] == new["mode"] == "benchmark", "Quadrature overhead must not enter benchmarks")
        period = old["settings"]["duration_s"] / old["settings"]["initial_target_periods"]
        remaining = (100 - old["settings"]["initial_target_periods"]) * period
        old_projected = old["wall_elapsed_s"] + remaining * old["timing"]["steady_wall_s_per_simulated_s"]
        new_projected = new["wall_elapsed_s"] + remaining * new["timing"]["steady_wall_s_per_simulated_s"]
        ratios[policy] = new_projected / old_projected
        projections[policy] = Dict("old_projected_100_orbit_s" => old_projected,
            "new_projected_100_orbit_s" => new_projected, "runtime_ratio" => ratios[policy],
            "method" => "10-orbit wall time plus 90 orbits at measured post-warmup propagation rate; startup counted once")
    end
    historical_hours = sum(b.wall for b in baselines) / 3600
    process_and_verifier_hours = 41.0 - historical_hours
    require(process_and_verifier_hours > 0, "Historical timings no longer match the 41-hour baseline")
    estimated_hours = sum(b.wall * ratios[string(c.policy)] for (b, c) in zip(baselines, cases)) / 3600 +
        process_and_verifier_hours
    planning_hours = [ceil(Int, 1.12estimated_hours), ceil(Int, 1.22estimated_hours)]
    expected_saved = sum(b.bytes for b in baselines) + length(cases) * 4 * TRAJECTORY_ROWS * sizeof(Float64)
    raw_trajectory = sum(TRAJECTORY_ROWS * length(trajectory_columns(c.helpers + 1)) *
        sizeof(Float64) for c in cases)
    conservative_saved = ceil(Int, raw_trajectory + 2sum(b.event_bytes for b in baselines) +
        length(cases) * 1024^2)
    max_spool = maximum((ceil(Int, mission_duration(c) / ENCOUNTER_SAMPLE_S) + 1) *
        (6(c.helpers + 1) + 1) * sizeof(Float64) for c in cases)
    available_disk = free_disk_bytes(DEFAULT_OUTPUT)
    required_disk = conservative_saved + max_spool + 5 * 1024^3
    memory_line = only(filter(line -> startswith(line, "MemAvailable:"), readlines("/proc/meminfo")))
    available_memory = parse(Int, split(memory_line)[2]) * 1024
    required_memory = max(6 * 1024^3, 2max_rss)
    require(available_disk >= required_disk, "Insufficient disk for the whole campaign plus spool/reserve")
    require(available_memory >= required_memory, "Insufficient available RAM for serial production")
    report = Dict{String,Any}("passed" => true, "checked_utc" => string(now(UTC)),
        "source_signature" => signature, "campaign" => "absolute_rn", "production_cases" => length(cases),
        "fixture_assertions" => assertion_count, "pilot_count" => length(PILOTS),
        "pilot_labels" => PILOTS, "source_snapshot" => snapshot, "comparisons" => comparisons,
        "maximum_independent_net_error_mps" => max_net,
        "maximum_independent_absolute_error_mps" => max_absolute,
        "maximum_summed_quadrature_difference_mps" => max_quadrature,
        "peak_pilot_rss_bytes" => max_rss, "available_memory_bytes" => available_memory,
        "required_memory_bytes" => required_memory, "available_disk_bytes" => available_disk,
        "required_disk_bytes" => required_disk, "expected_persistent_bytes" => round(Int, expected_saved),
        "conservative_persistent_bytes" => conservative_saved, "maximum_temporary_spool_bytes" => max_spool,
        "historical_simulation_hours" => historical_hours, "historical_baseline_hours" => 41.0,
        "historical_interpolated_cases" => count(b -> b.interpolated, baselines),
        "estimated_campaign_hours" => estimated_hours, "planning_allowance_hours" => planning_hours,
        "benchmark_projections" => projections, "old_helpers700_preservation" => preserved_700)
    atomic_toml(result_path, report)
    open(joinpath(RESULTS, "VALIDATION_REPORT.md"), "w") do io
        println(io, "# Absolute R/N campaign validation\n\nChecked UTC: $(report["checked_utc"]).")
        println(io, "\n## Scope and source preservation\n")
        println(io, "- Exactly **270 fresh cases**, 100 periods each; 1001 rows with endpoints.")
        println(io, "- Targets 900/950/1000/1050/1100 km; inclinations 0/0.5/1 degree.")
        println(io, "- Helpers initially at 1000 km; counts 1/100/200/300/400/500/600/700/800; both schedulers.")
        println(io, "- Row 801 is the 80-orbit endpoint. No separate 80-orbit runs.")
        println(io, "- Original runner files match the pre-change snapshot. Old N700 remains stopped, 4 verified cases preserved.")
        println(io, "- [Source archive](source_snapshot_pre_absolute_RN/sources.tar.gz): $(snapshot["file_count"]) checksummed files; SHA256 `$(snapshot["archive_sha256"])`.")
        println(io, "- Shared library edits are limited to diagnostic RHS/initialization, default saving and checkpoint validation.")
        println(io, "\n## Tests and independent validation\n")
        println(io, "- **$assertion_count assertions passed** in [the focused suite](fixtures_final.log), including existing interlink/scheduler tests and enabled checkpoint round trips.")
        println(io, "- **$(length(PILOTS)) pilots/baselines** passed fresh-file verification; pilot data remain in [PilotOutput](../PilotOutput).")
        println(io, "- Net and absolute integrals are independently reconstructed from dense Cartesian states and the actually held link, not saved diagnostics.")
        println(io, "- Absolute acceleration is integrated directly, with zero-crossing splitting and 16-versus-8-panel Simpson convergence checks.")
        @printf(io, "- Maximum independent net error: %.3g m/s; absolute-integral error: %.3g m/s; summed quadrature difference: %.3g m/s.\n",
            max_net, max_absolute, max_quadrature)
        println(io, "- R/N branch tolerances: 1e-8 m/s; existing T allowance: 1e-7 m/s; net-sum tolerance: 1e-8 m/s. No values were clamped and no limits were relaxed.")
        println(io, "- Independent quadrature covers the short instrumented pilots. Their measured errors are not claimed as a uniform 100-orbit cumulative-error bound; production still checks every saved row and the physical/event invariants.")
        println(io, "- Inclined lower/higher-target cases exercise both R signs for each scheduler; N has both signs. See [sign coverage](sign_coverage.toml).")
        println(io, "- Zero-power and coplanar inclination pilots have exactly zero laser diagnostics. Recorder/quadrature instrumentation leaves the same-layout dynamics and schedule bitwise unchanged.")
        println(io, "\n## Frozen-library physical comparisons\n")
        println(io, "| Comparison | Maximum position difference (m) | Maximum velocity difference (m/s) | Held-schedule disagreement (s) |")
        println(io, "|---|---:|---:|---:|")
        for name in COMPARISONS
            values = comparisons[name]
            @printf(io, "| %s | %.6g | %.6g | %.6g |\n", name,
                values["maximum_position_difference_m"], values["maximum_velocity_difference_mps"],
                values["schedule_disagreement_s"])
        end
        println(io, "\nAdded states may change adaptive steps. Let a be maximum laser acceleration, D the measured time when held links differ, T the duration and n the initial mean motion. The acceleration difference is at most 2a during those intervals. Conservative near-circular comparison tolerances are 0.02 + 10*T*a*D metres and 5e-5 + 14*a*D + 1.01*n*(10*T*a*D) m/s. The latter includes rotation of the velocity vector as phase changes. The schedule-disagreement limit is 2%. These are comparison allowances, not solver tolerances. The RHS tests separately require identical physical derivatives at identical states.")
        println(io, "\n## Serial 800-helper benchmarks\n")
        println(io, "| Implementation/policy | Orbits | Simulation wall (s) | Post-warmup wall/orbit (s) | Peak RSS (GiB) | Temporary spool (MiB) |")
        println(io, "|---|---:|---:|---:|---:|---:|")
        for label in ("baseline_sma800", "new_sma800", "baseline_inc800", "new_inc800")
            values = reports[label]
            period = values["settings"]["duration_s"] / values["settings"]["initial_target_periods"]
            @printf(io, "| %s | %.0f | %.2f | %.2f | %.3f | %.2f |\n", label,
                values["settings"]["initial_target_periods"], values["wall_elapsed_s"],
                period * values["timing"]["steady_wall_s_per_simulated_s"],
                values["peak_rss_bytes"] / 1024^3, values["temporary_dense_state_bytes"] / 1024^2)
        end
        @printf(io, "\nUpdated campaign estimate: **%.1f hours**; planning allowance: **%d-%d hours**.\n",
            estimated_hours, planning_hours[1], planning_hours[2])
        println(io, "\nThe estimate starts with 240 observed old cases plus 30 N700 interpolations, adds the existing process/verification allowance to recover the 41-hour baseline, and applies each scheduler's measured old/new 100-orbit projection ratio. A projection counts startup once, then adds 90 periods at the post-warmup rate. It does not multiply a short pilot's entire wall time by ten or fifty. Hardware load and case-to-case differences remain uncertain.")
        println(io, "\n## Capacity and launch readiness\n")
        @printf(io, "- Available disk: %.2f GiB. Expected persistent output: %.2f GiB; conservative allowance: %.2f GiB.\n",
            available_disk / 1024^3, expected_saved / 1024^3, conservative_saved / 1024^3)
        @printf(io, "- Largest one-case spool: %.2f GiB; total disk requirement including the 5 GiB reserve: %.2f GiB.\n",
            max_spool / 1024^3, required_disk / 1024^3)
        @printf(io, "- Largest measured pilot RSS: %.2f GiB; available memory: %.2f GiB; required margin: %.2f GiB.\n",
            max_rss / 1024^3, available_memory / 1024^3, required_memory / 1024^3)
        println(io, "- No existing data was deleted to obtain capacity. In-run and pre-case disk guards remain enabled.")
        println(io, "- Validated runtime source signature: `$signature`.")
        println(io, "\nPreflight passed. Actual launch and the first production verification are recorded separately in [launch confirmation](launch_confirmation.toml). Monitor/stop/resume commands are in the [runner README](../../../Simulation_absolute_RN/Code/README.md).")
    end
    TOML.print(stdout, report; sorted=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
