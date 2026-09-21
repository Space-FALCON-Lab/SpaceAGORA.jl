const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT_DIR = joinpath(REPO_ROOT, "output", "performance", "split_imex_compare")
const DEFAULT_CASES = (
    "single_j2",
    "single_orientation_aero",
    "proximity_2sat_orientation_fullstack_gnc_highrate",
    "multi_4_gravity",
)
const SOLVER_CONFIGS = (
    (name="tsit5", solver_mode="tsit5", split_solver=nothing),
    (name="auto_stiff", solver_mode="auto_stiff", split_solver=nothing),
    (name="split_imex_kencarp4", solver_mode="split_imex", split_solver="kencarp4"),
    (name="split_imex_kencarp47", solver_mode="split_imex", split_solver="kencarp47"),
    (name="multirate", solver_mode="multirate", split_solver=nothing),
)

using CSV
using DataFrames
using Dates
using Statistics

module PerfRuntimeCompare
include(joinpath(@__DIR__, "performance_runtime_analysis.jl"))
end

const RT = PerfRuntimeCompare

@inline function _parse_positive_int(raw::AbstractString, flag::String)::Int
    txt = String(raw)
    parsed = try
        parse(Int, txt)
    catch
        throw(ArgumentError("$flag must be an integer, got '$txt'."))
    end
    parsed > 0 || throw(ArgumentError("$flag must be > 0, got $parsed."))
    return parsed
end

@inline function _parse_case_list(raw::AbstractString)::Vector{String}
    names = String[]
    for token in split(raw, ",")
        t = strip(token)
        isempty(t) && continue
        push!(names, t)
    end
    isempty(names) && throw(ArgumentError("--cases must include at least one scenario name."))
    return names
end

function parse_cli()
    profile_name = lowercase(get(ENV, "SPACEAGORA_PERF_PROFILE", "quick"))
    outdir = get(ENV, "SPACEAGORA_PERF_OUTDIR", DEFAULT_OUTPUT_DIR)
    case_names = collect(String.(DEFAULT_CASES))
    warmup_override = nothing
    repeats_override = nothing
    max_attempts_override = nothing

    for arg in ARGS
        if arg in ("quick", "full")
            profile_name = arg
        elseif startswith(arg, "--profile=")
            profile_name = lowercase(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--cases=")
            case_names = _parse_case_list(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--warmup=")
            warmup_override = _parse_positive_int(split(arg, "=", limit=2)[2], "--warmup")
        elseif startswith(arg, "--repeats=")
            repeats_override = _parse_positive_int(split(arg, "=", limit=2)[2], "--repeats")
        elseif startswith(arg, "--max-attempts=")
            max_attempts_override = _parse_positive_int(split(arg, "=", limit=2)[2], "--max-attempts")
        else
            throw(ArgumentError(
                "Unknown argument '$arg'. Supported: [quick|full], --profile=..., --outdir=..., " *
                "--cases=name1,name2,..., --warmup=..., --repeats=..., --max-attempts=..."
            ))
        end
    end

    return (
        profile_name=profile_name,
        outdir=abspath(outdir),
        case_names=case_names,
        warmup_override=warmup_override,
        repeats_override=repeats_override,
        max_attempts_override=max_attempts_override
    )
end

function _resolve_profile(profile_name::String, warmup_override, repeats_override, max_attempts_override)::RT.ProfileSpec
    base = RT._profile_from_name(profile_name)
    return RT.ProfileSpec(
        name=base.name,
        repeats=isnothing(repeats_override) ? base.repeats : repeats_override,
        warmup=isnothing(warmup_override) ? base.warmup : warmup_override,
        max_attempts=isnothing(max_attempts_override) ? base.max_attempts : max_attempts_override,
        mission_short_s=base.mission_short_s,
        mission_long_s=base.mission_long_s,
        montecarlo_samples=base.montecarlo_samples,
        montecarlo_mission_s=base.montecarlo_mission_s
    )
end

function _find_cases(all_cases::Vector{RT.BenchmarkCase}, case_names::Vector{String})::Vector{RT.BenchmarkCase}
    by_name = Dict(c.name => c for c in all_cases)
    selected = RT.BenchmarkCase[]
    for name in case_names
        haskey(by_name, name) || throw(ArgumentError("Unknown case '$name'."))
        push!(selected, by_name[name])
    end
    return selected
end

@inline function _override_solver(case::RT.BenchmarkCase, solver_mode::String)::RT.BenchmarkCase
    return RT.BenchmarkCase(
        name=case.name,
        category=case.category,
        description=case.description,
        args_template=case.args_template,
        run_in_quick=case.run_in_quick,
        solver_mode_override=solver_mode
    )
end

function _run_solver_case_grid(
    case::RT.BenchmarkCase,
    spec::RT.ProfileSpec,
    solver_cfg
)::Vector{NamedTuple}
    case_mode = _override_solver(case, solver_cfg.solver_mode)
    plan = RT.parallel_priority_plan(case_mode, :none)
    rows = NamedTuple[]
    env_pairs = Pair{String, String}[
        "SPACEAGORA_PERF_PARALLEL" => "0",
        "SPACEAGORA_PERF_PARALLEL_BACKEND" => "none"
    ]
    if !(solver_cfg.split_solver === nothing)
        push!(env_pairs, "SPACEAGORA_SPLIT_IMEX_SOLVER" => solver_cfg.split_solver)
    end

    withenv(env_pairs...) do
        RT.run_warmup(case_mode, spec.warmup, spec.name)
        for rep in 1:spec.repeats
            last_row = nothing
            for attempt in 1:spec.max_attempts
                row = RT.measure_case(case_mode, spec.name, rep; attempt=attempt, plan=plan)
                row = merge(
                    row,
                    (
                        solver_config=solver_cfg.name,
                        split_imex_solver=(solver_cfg.split_solver === nothing ? "n/a" : String(solver_cfg.split_solver))
                    )
                )
                last_row = row
                if row.solve_success
                    push!(rows, row)
                    break
                end
            end
            if !(last_row === nothing) && !last_row.solve_success
                push!(rows, last_row)
            end
        end
    end
    return rows
end

@inline _safe_mean(v) = isempty(v) ? missing : mean(v)
@inline _safe_median(v) = isempty(v) ? missing : median(v)
@inline _safe_p90(v) = isempty(v) ? missing : quantile(v, 0.9)

function _build_summary(raw_df::DataFrame)::DataFrame
    grouped = combine(
        groupby(raw_df, [:scenario, :solver_config, :split_imex_solver]),
        :solve_success => (v -> sum(Bool.(v))) => :success_count,
        :solve_success => length => :sample_count,
        :total_time_s => (v -> _safe_mean(Float64.(v))) => :total_time_mean_s,
        :total_time_s => (v -> _safe_median(Float64.(v))) => :total_time_median_s,
        :total_time_s => (v -> _safe_p90(Float64.(v))) => :total_time_p90_s,
        :solve_time_s => (v -> _safe_mean(Float64.(v))) => :solve_time_mean_s,
        :accepted_steps => (v -> _safe_mean(Float64.(collect(skipmissing(v))))) => :accepted_steps_mean,
        :rejected_steps => (v -> _safe_mean(Float64.(collect(skipmissing(v))))) => :rejected_steps_mean,
        :solver_sequence => (v -> RT._safe_unique_join(v; delimiter="|")) => :solver_sequences
    )
    grouped[!, :success_rate] = [
        n <= 0 ? missing : Float64(ok) / Float64(n)
        for (ok, n) in zip(grouped.success_count, grouped.sample_count)
    ]

    baseline = Dict{String, Float64}()
    for row in eachrow(grouped)
        if row.solver_config == "auto_stiff" &&
           !(row.total_time_mean_s isa Missing) &&
           isfinite(Float64(row.total_time_mean_s)) &&
           Float64(row.total_time_mean_s) > 0.0
            baseline[String(row.scenario)] = Float64(row.total_time_mean_s)
        end
    end
    grouped[!, :speedup_vs_auto_stiff] = [
        begin
            base = get(baseline, String(sc), NaN)
            val = tt isa Missing ? NaN : Float64(tt)
            (isfinite(base) && isfinite(val) && val > 0.0) ? (base / val) : NaN
        end
        for (sc, tt) in zip(grouped.scenario, grouped.total_time_mean_s)
    ]

    sort!(grouped, [:scenario, :solver_config])
    return grouped
end

function _write_report(
    path::String,
    spec::RT.ProfileSpec,
    case_names::Vector{String},
    raw_path::String,
    summary_path::String,
    summary_df::DataFrame
)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Split-IMEX Solver Comparison")
        println(io)
        println(io, "- Generated (UTC): $(string(now(UTC)))")
        println(io, "- Profile: `$(spec.name)`")
        println(io, "- Warmup: `$(spec.warmup)`")
        println(io, "- Repeats: `$(spec.repeats)`")
        println(io, "- Max attempts: `$(spec.max_attempts)`")
        println(io, "- Scenarios: `$(join(case_names, ", "))`")
        println(io)
        println(io, "Artifacts:")
        println(io, "- Raw CSV: `$(raw_path)`")
        println(io, "- Summary CSV: `$(summary_path)`")
        println(io)
        println(io, "| Scenario | Solver Config | Mean Total [s] | Speedup vs auto_stiff [x] | Success Rate | Solver Sequence |")
        println(io, "|---|---:|---:|---:|---:|---|")
        for row in eachrow(summary_df)
            mean_s = row.total_time_mean_s isa Missing ? "n/a" : string(round(Float64(row.total_time_mean_s); digits=4))
            speedup = isfinite(row.speedup_vs_auto_stiff) ? string(round(Float64(row.speedup_vs_auto_stiff); digits=4)) : "n/a"
            sr = row.success_rate isa Missing ? "n/a" : string(round(100.0 * Float64(row.success_rate); digits=1), "%")
            seq = row.solver_sequences isa Missing ? "n/a" : String(row.solver_sequences)
            println(io, "| $(row.scenario) | $(row.solver_config) | $(mean_s) | $(speedup) | $(sr) | $(seq) |")
        end
    end
    return nothing
end

function main()
    cli = parse_cli()
    spec = _resolve_profile(cli.profile_name, cli.warmup_override, cli.repeats_override, cli.max_attempts_override)
    mkpath(cli.outdir)

    println("Split-IMEX compare profile: $(spec.name)")
    println("Output directory: $(cli.outdir)")
    println("Cases: $(join(cli.case_names, ", "))")
    println("Solvers: $(join(String[cfg.name for cfg in SOLVER_CONFIGS], ", "))")

    planet = RT.Earth("", RT.SPICE_PATH)
    all_cases = RT.build_cases(spec, planet)
    selected = _find_cases(all_cases, cli.case_names)

    rows = NamedTuple[]
    for (idx, case) in enumerate(selected)
        println("[$idx/$(length(selected))] scenario=$(case.name)")
        for solver_cfg in SOLVER_CONFIGS
            println("  solver=$(solver_cfg.name)")
            append!(rows, _run_solver_case_grid(case, spec, solver_cfg))
        end
    end

    raw_df = DataFrame(rows)
    summary_df = _build_summary(raw_df)
    stamp = Dates.format(now(UTC), "yyyymmdd_HHMMSS")
    raw_path = joinpath(cli.outdir, "split_imex_compare_raw_$(spec.name)_$(stamp).csv")
    summary_path = joinpath(cli.outdir, "split_imex_compare_summary_$(spec.name)_$(stamp).csv")
    report_path = joinpath(cli.outdir, "split_imex_compare_report_$(spec.name)_$(stamp).md")
    CSV.write(raw_path, raw_df)
    CSV.write(summary_path, summary_df)
    _write_report(report_path, spec, cli.case_names, raw_path, summary_path, summary_df)

    println("Wrote raw CSV: $(raw_path)")
    println("Wrote summary CSV: $(summary_path)")
    println("Wrote report: $(report_path)")
    println()
    println("Run recipe:")
    println("  julia --project=. benchmarks/studies/performance_split_imex_compare.jl --profile=quick")
    println("  julia --project=. benchmarks/studies/performance_split_imex_compare.jl --profile=quick --warmup=1 --repeats=3 --cases=single_orientation_aero,proximity_2sat_orientation_fullstack_gnc_highrate,multi_4_gravity")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
