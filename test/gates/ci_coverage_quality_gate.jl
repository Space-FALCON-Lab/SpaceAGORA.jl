const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SRC_ROOT = joinpath(REPO_ROOT, "src")

using Printf

const EXCLUDED_FROM_MAIN_GATE = Set{String}()

const LEGACY_MIN_SMOKE_COVERAGE = Dict{String, Float64}()

const MIN_MAIN_OVERALL = let raw = get(ENV, "SPACEAGORA_COVERAGE_MIN_OVERALL", "90.0")
    parsed = tryparse(Float64, raw)
    parsed === nothing && error("Invalid SPACEAGORA_COVERAGE_MIN_OVERALL=$raw")
    parsed
end

const MIN_MAIN_FILE = let raw = get(ENV, "SPACEAGORA_COVERAGE_MIN_FILE", "80.0")
    parsed = tryparse(Float64, raw)
    parsed === nothing && error("Invalid SPACEAGORA_COVERAGE_MIN_FILE=$raw")
    parsed
end

const MAIN_FILE_MIN_OVERRIDES = Dict(
    # Adapter-only file is exercised indirectly via config/entrypoint tests.
    "src/simulation/engine/adapters/from_env.jl" => 70.0,
    # Dynamic RHS is branch-heavy across many mission/control combinations.
    "src/simulation/engine/dynamics_rhs.jl" => 70.0,

)

const CRITICAL_FILE_MIN_OVERRIDES = Dict(
    "src/simulation/engine/execution.jl" => 90.0,
    "src/gnc/control/propulsive_maneuvers.jl" => 90.0,
    "src/core/interfaces/reference_system.jl" => 90.0,
)

const COVERAGE_WINDOW_SECONDS = let raw = get(ENV, "SPACEAGORA_COVERAGE_WINDOW_SECONDS", "900")
    parsed = tryparse(Float64, raw)
    parsed === nothing && error("Invalid SPACEAGORA_COVERAGE_WINDOW_SECONDS=$raw")
    Int(round(parsed))
end

struct CoverageSummary
    path::String
    covered::Int
    executable::Int
    percent::Float64
    excluded::Bool
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

function list_active_cov_files()
    cov_files = String[]
    mtimes = Float64[]
    for (root, _, files) in walkdir(SRC_ROOT)
        for file in files
            if endswith(file, ".cov")
                path = joinpath(root, file)
                push!(cov_files, path)
                push!(mtimes, stat(path).mtime)
            end
        end
    end

    isempty(cov_files) && error("No coverage files found under $(SRC_ROOT). Run tests with --code-coverage=user first.")

    newest = maximum(mtimes)
    active = [cov_files[i] for i in eachindex(cov_files) if newest - mtimes[i] <= COVERAGE_WINDOW_SECONDS]
    isempty(active) && error("No active coverage files found in the last $(COVERAGE_WINDOW_SECONDS)s window.")

    println("coverage_files_total=$(length(cov_files)) active_window=$(length(active)) window_seconds=$(COVERAGE_WINDOW_SECONDS)")
    return sort(active)
end

function summarize_coverage(cov_files::Vector{String})
    executed_counts = Dict{Tuple{String, Int}, Int}()
    executable_lines = Set{Tuple{String, Int}}()

    for cov_file in cov_files
        src_rel = source_path_from_cov(cov_file)
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
    for src_rel in sort(collect(keys(per_file_exec)))
        exec = per_file_exec[src_rel]
        cov = get(per_file_cov, src_rel, 0)
        pct = exec == 0 ? 100.0 : 100.0 * cov / exec
        excluded = src_rel in EXCLUDED_FROM_MAIN_GATE
        push!(summaries, CoverageSummary(src_rel, cov, exec, pct, excluded))
    end

    return summaries
end

function print_summary(summaries::Vector{CoverageSummary})
    main = sort(filter(s -> !s.excluded, summaries); by=s -> s.percent)
    excluded = sort(filter(s -> s.excluded, summaries); by=s -> s.path)

    main_cov = sum(s.covered for s in main)
    main_exec = sum(s.executable for s in main)
    main_pct = main_exec == 0 ? 100.0 : 100.0 * main_cov / main_exec

    println(@sprintf("main_overall=%.2f%% (%d/%d) threshold=%.2f%%", main_pct, main_cov, main_exec, MIN_MAIN_OVERALL))
    println("lowest_main_files:")
    for s in Iterators.take(main, 10)
        println(@sprintf("  %.2f%% (%d/%d)  %s", s.percent, s.covered, s.executable, s.path))
    end

    if !isempty(excluded)
        println("excluded_legacy_files:")
        for s in excluded
            println(@sprintf("  %.2f%% (%d/%d)  %s", s.percent, s.covered, s.executable, s.path))
        end
    end

    summary_by_path = Dict(s.path => s for s in summaries)
    println("critical_files:")
    for path in sort(collect(keys(CRITICAL_FILE_MIN_OVERRIDES)))
        threshold = CRITICAL_FILE_MIN_OVERRIDES[path]
        summary = get(summary_by_path, path, nothing)
        if summary === nothing
            println(@sprintf("  MISSING (threshold %.2f%%)  %s", threshold, path))
        else
            println(@sprintf("  %.2f%% (threshold %.2f%%)  %s", summary.percent, threshold, path))
        end
    end
end

function enforce_gate(summaries::Vector{CoverageSummary})
    failures = String[]

    main = filter(s -> !s.excluded, summaries)
    main_cov = sum(s.covered for s in main)
    main_exec = sum(s.executable for s in main)
    main_pct = main_exec == 0 ? 100.0 : 100.0 * main_cov / main_exec

    if main_pct < MIN_MAIN_OVERALL
        push!(failures, @sprintf("Main overall coverage %.2f%% is below threshold %.2f%%", main_pct, MIN_MAIN_OVERALL))
    end

    for s in main
        threshold = get(MAIN_FILE_MIN_OVERRIDES, s.path, MIN_MAIN_FILE)
        if s.percent < threshold
            push!(failures, @sprintf("Main file coverage %.2f%% is below %.2f%% for %s", s.percent, threshold, s.path))
        end
    end

    summary_by_path = Dict(s.path => s for s in summaries)
    for (critical_path, critical_min) in CRITICAL_FILE_MIN_OVERRIDES
        summary = get(summary_by_path, critical_path, nothing)
        if summary === nothing
            push!(failures, "Critical file has no coverage artifact: $critical_path")
            continue
        end
        if summary.percent < critical_min
            push!(failures, @sprintf("Critical file coverage %.2f%% is below %.2f%% for %s", summary.percent, critical_min, critical_path))
        end
    end

    for (legacy_path, min_pct) in LEGACY_MIN_SMOKE_COVERAGE
        summary = get(summary_by_path, legacy_path, nothing)
        if summary === nothing
            push!(failures, "Legacy excluded file has no coverage artifact: $legacy_path")
            continue
        end
        if summary.percent < min_pct
            push!(failures, @sprintf("Legacy excluded file %s smoke coverage %.2f%% is below %.2f%%", legacy_path, summary.percent, min_pct))
        end
    end

    if !isempty(failures)
        println("coverage_gate_failures:")
        for msg in failures
            println("  - $msg")
        end
        error("Coverage quality gate failed")
    end
end

cov_files = list_active_cov_files()
summaries = summarize_coverage(cov_files)
print_summary(summaries)
enforce_gate(summaries)
println("coverage_quality_gate_ok")
