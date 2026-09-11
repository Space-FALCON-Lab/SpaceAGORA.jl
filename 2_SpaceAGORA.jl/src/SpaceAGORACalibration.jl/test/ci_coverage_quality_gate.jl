const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const SRC_ROOT = joinpath(REPO_ROOT, "src")

using Printf

const MIN_OVERALL = let raw = get(ENV, "SPACEAGORA_CAL_COVERAGE_MIN_OVERALL", "90.0")
    parsed = tryparse(Float64, raw)
    parsed === nothing && error("Invalid SPACEAGORA_CAL_COVERAGE_MIN_OVERALL=$raw")
    parsed
end

const MIN_FILE = let raw = get(ENV, "SPACEAGORA_CAL_COVERAGE_MIN_FILE", "70.0")
    parsed = tryparse(Float64, raw)
    parsed === nothing && error("Invalid SPACEAGORA_CAL_COVERAGE_MIN_FILE=$raw")
    parsed
end

const COVERAGE_WINDOW_SECONDS = let raw = get(ENV, "SPACEAGORA_CAL_COVERAGE_WINDOW_SECONDS", "900")
    parsed = tryparse(Float64, raw)
    parsed === nothing && error("Invalid SPACEAGORA_CAL_COVERAGE_WINDOW_SECONDS=$raw")
    Int(round(parsed))
end

struct CoverageSummary
    path::String
    covered::Int
    executable::Int
    percent::Float64
end

@inline function source_path_from_cov(cov_path::String)::String
    rel_cov = relpath(cov_path, REPO_ROOT)
    if endswith(rel_cov, ".cov")
        m = match(r"^(.*)\.[0-9]+\.cov$", rel_cov)
        if m !== nothing
            return m.captures[1]
        end
        return rel_cov[1:(end - 4)]
    end
    error("Unexpected coverage file suffix: $cov_path")
end

function list_source_files()
    src_files = String[]
    for (root, _, files) in walkdir(SRC_ROOT)
        for file in files
            endswith(file, ".jl") || continue
            push!(src_files, relpath(joinpath(root, file), REPO_ROOT))
        end
    end
    sort!(src_files)
    return src_files
end

function list_active_cov_files()
    cov_files = String[]
    mtimes = Float64[]
    for (root, _, files) in walkdir(SRC_ROOT)
        for file in files
            endswith(file, ".cov") || continue
            path = joinpath(root, file)
            push!(cov_files, path)
            push!(mtimes, stat(path).mtime)
        end
    end

    isempty(cov_files) && error("No coverage files found under $(SRC_ROOT). Run tests with --code-coverage=user first.")

    newest = maximum(mtimes)
    active = [cov_files[i] for i in eachindex(cov_files) if newest - mtimes[i] <= COVERAGE_WINDOW_SECONDS]
    isempty(active) && error("No active coverage files found in the last $(COVERAGE_WINDOW_SECONDS)s window.")

    println("calibration_coverage_files_total=$(length(cov_files)) active_window=$(length(active)) window_seconds=$(COVERAGE_WINDOW_SECONDS)")
    return sort(active)
end

function summarize_coverage(cov_files::Vector{String})
    executed_counts = Dict{Tuple{String, Int}, Int}()
    executable_lines = Set{Tuple{String, Int}}()
    files_with_cov = Set{String}()

    for cov_file in cov_files
        src_rel = source_path_from_cov(cov_file)
        push!(files_with_cov, src_rel)
        for (lineno, line) in enumerate(eachline(cov_file))
            m = match(r"^\s*([0-9]+|-)\s", line)
            m === nothing && continue
            token = m.captures[1]
            token == "-" && continue

            key = (src_rel, lineno)
            push!(executable_lines, key)
            executed_counts[key] = get(executed_counts, key, 0) + parse(Int, token)
        end
    end

    per_file_exec = Dict{String, Int}()
    per_file_cov = Dict{String, Int}()
    for (src_rel, line_no) in executable_lines
        key = (src_rel, line_no)
        per_file_exec[src_rel] = get(per_file_exec, src_rel, 0) + 1
        if get(executed_counts, key, 0) > 0
            per_file_cov[src_rel] = get(per_file_cov, src_rel, 0) + 1
        end
    end

    summaries = CoverageSummary[]
    for src_rel in sort(collect(files_with_cov))
        exec = get(per_file_exec, src_rel, 0)
        cov = get(per_file_cov, src_rel, 0)
        pct = exec == 0 ? 100.0 : 100.0 * cov / exec
        push!(summaries, CoverageSummary(src_rel, cov, exec, pct))
    end

    return summaries
end

function print_summary(summaries::Vector{CoverageSummary})
    ordered = sort(summaries; by = s -> s.percent)
    cov = sum(s.covered for s in ordered)
    exec = sum(s.executable for s in ordered)
    pct = exec == 0 ? 100.0 : 100.0 * cov / exec

    println(@sprintf("calibration_overall=%.2f%% (%d/%d) threshold=%.2f%%", pct, cov, exec, MIN_OVERALL))
    println("calibration_files:")
    for s in ordered
        println(@sprintf("  %.2f%% (%d/%d)  %s", s.percent, s.covered, s.executable, s.path))
    end
end

function enforce_gate(summaries::Vector{CoverageSummary}, src_files::Vector{String})
    failures = String[]
    missing_artifacts = String[]

    total_cov = sum(s.covered for s in summaries)
    total_exec = sum(s.executable for s in summaries)
    overall = total_exec == 0 ? 100.0 : 100.0 * total_cov / total_exec
    if overall < MIN_OVERALL
        push!(failures, @sprintf("Calibration overall coverage %.2f%% is below threshold %.2f%%", overall, MIN_OVERALL))
    end

    summary_by_path = Dict(s.path => s for s in summaries)
    for src_file in src_files
        if !haskey(summary_by_path, src_file)
            push!(missing_artifacts, src_file)
            continue
        end
        pct = summary_by_path[src_file].percent
        if pct < MIN_FILE
            push!(failures, @sprintf("Calibration file coverage %.2f%% is below %.2f%% for %s", pct, MIN_FILE, src_file))
        end
    end

    if !isempty(missing_artifacts)
        println("calibration_coverage_missing_artifacts:")
        for path in missing_artifacts
            println("  - $path")
        end
    end

    if !isempty(failures)
        println("calibration_coverage_gate_failures:")
        for msg in failures
            println("  - $msg")
        end
        error("Calibration coverage quality gate failed")
    end
end

src_files = list_source_files()
cov_files = list_active_cov_files()
summaries = summarize_coverage(cov_files)
print_summary(summaries)
enforce_gate(summaries, src_files)
println("calibration_coverage_quality_gate_ok")
