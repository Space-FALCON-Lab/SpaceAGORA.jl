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

@inline function source_path_from_cov(cov_path::String; repo_root::String=REPO_ROOT)::String
    rel_cov = relpath(cov_path, repo_root)
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

# Precompile-only code. PrecompileTools runs the body of `@setup_workload` and
# `@compile_workload` only while the package precompiles; a process loaded
# with `--code-coverage=user` never executes those blocks, so their lines would
# always read as uncovered. That is a measurement artifact, not untested
# behavior, so the lines inside those macro calls (the call lines included) are
# dropped from the executable-line set. Functions merely *called* from a
# workload are ordinary code and keep counting.
const PRECOMPILE_ONLY_MACROS = ("setup_workload", "compile_workload")

const _JS = Base.JuliaSyntax

# The macro's bare name ("setup_workload" for `@setup_workload` and for
# `PrecompileTools.@setup_workload`), or `nothing` when `node` is not a macro call.
function _macrocall_name(node)::Union{Nothing, String}
    _JS.kind(node) == _JS.K"macrocall" || return nothing
    kids = _JS.children(node)
    (kids === nothing || isempty(kids)) && return nothing
    head = kids[1]
    while _JS.kind(head) == _JS.K"."
        head_kids = _JS.children(head)
        (head_kids === nothing || isempty(head_kids)) && return nothing
        head = head_kids[end]
    end
    _JS.kind(head) == _JS.K"MacroName" || return nothing
    return String(_JS.sourcetext(head))
end

function _collect_precompile_spans!(spans::Vector{UnitRange{Int}}, node)
    name = _macrocall_name(node)
    if name !== nothing && name in PRECOMPILE_ONLY_MACROS
        first_line = _JS.source_line(node.source, _JS.first_byte(node))
        last_line = _JS.source_line(node.source, _JS.last_byte(node))
        push!(spans, first_line:last_line)
        return spans  # nested workload macros are inside this span already
    end
    kids = _JS.children(node)
    kids === nothing && return spans
    for kid in kids
        _collect_precompile_spans!(spans, kid)
    end
    return spans
end

"""
    precompile_only_spans(text; filename="") -> Vector{UnitRange{Int}}

Line spans of every top-most `@setup_workload` / `@compile_workload` macro call
in `text`, found by parsing with Julia's bundled parser. A file that does not
parse is an error: the gate refuses to guess which lines are precompile-only.
"""
function precompile_only_spans(text::AbstractString; filename::AbstractString="")
    tree = try
        _JS.parseall(_JS.SyntaxNode, text; filename=String(filename), ignore_errors=false, ignore_warnings=true)
    catch err
        error("Coverage gate could not parse $(isempty(filename) ? "source" : filename) to find precompile-only blocks: " *
              sprint(showerror, err))
    end
    spans = _collect_precompile_spans!(UnitRange{Int}[], tree)
    return sort!(spans; by=first)
end

"""
    exclude_precompile_only_lines!(executable_lines; repo_root, io) -> Int

Remove from `executable_lines` every `(src_rel, line)` inside a precompile-only
macro span of that source file, print one `precompile_only_excluded` line per
file that has such a span, and a total. Returns the number of lines removed.
"""
function exclude_precompile_only_lines!(executable_lines::Set{Tuple{String, Int}};
                                        repo_root::String=REPO_ROOT, io::IO=stdout)
    total_removed = 0
    files_with_spans = 0
    for src_rel in sort!(unique(first.(collect(executable_lines))))
        src_path = joinpath(repo_root, src_rel)
        if !isfile(src_path)
            println(io, "precompile_only_scan_skipped  $(src_rel)  reason=source_missing")
            continue
        end
        spans = precompile_only_spans(read(src_path, String); filename=src_rel)
        isempty(spans) && continue
        files_with_spans += 1
        removed = Int[]
        for span in spans, line_no in span
            key = (src_rel, line_no)
            if key in executable_lines
                delete!(executable_lines, key)
                push!(removed, line_no)
            end
        end
        sort!(unique!(removed))
        total_removed += length(removed)
        span_text = join((first(s) == last(s) ? string(first(s)) : "$(first(s))-$(last(s))" for s in spans), ",")
        println(io, "precompile_only_excluded  $(src_rel)  lines=$(length(removed)) ($(span_text))")
    end
    println(io, "precompile_only_excluded_total  files=$(files_with_spans) lines=$(total_removed)")
    return total_removed
end

function summarize_coverage(cov_files::Vector{String}; repo_root::String=REPO_ROOT, io::IO=stdout)
    executed_counts = Dict{Tuple{String, Int}, Int}()
    executable_lines = Set{Tuple{String, Int}}()

    for cov_file in cov_files
        src_rel = source_path_from_cov(cov_file; repo_root=repo_root)
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

    exclude_precompile_only_lines!(executable_lines; repo_root=repo_root, io=io)

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

function run_coverage_quality_gate()
    cov_files = list_active_cov_files()
    summaries = summarize_coverage(cov_files)
    print_summary(summaries)
    enforce_gate(summaries)
    println("coverage_quality_gate_ok")
    return summaries
end

# Run only when executed as a script, so tests can include this file for its
# functions. test/coverage/runtests.jl calls `run_coverage_quality_gate()`.
if abspath(PROGRAM_FILE) == @__FILE__
    run_coverage_quality_gate()
end
