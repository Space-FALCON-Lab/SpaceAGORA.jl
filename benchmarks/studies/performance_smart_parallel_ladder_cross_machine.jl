const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_CROSS_MACHINE_OUTDIR = joinpath(REPO_ROOT, "output", "performance", "smart_parallel_ladder_cross_machine")

using CSV
using DataFrames
using Dates
using Statistics

Base.@kwdef struct CrossMachineInput
    label::String
    outdir::String
end

Base.@kwdef struct CrossMachineConfig
    profile::String
    outdir::String
    clean::Bool
    inputs::Vector{CrossMachineInput}
end

@inline function _parse_bool_token(raw::AbstractString)::Bool
    token = lowercase(strip(String(raw)))
    if token in ("1", "true", "yes", "on")
        return true
    elseif token in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid boolean token '$raw'. Use one of: 1/0, true/false, yes/no, on/off."))
end

@inline function _default_input_label(path::String)::String
    base = basename(normpath(path))
    return isempty(base) ? "machine" : base
end

@inline function _parse_input_token(raw::AbstractString)::CrossMachineInput
    token = strip(String(raw))
    isempty(token) && throw(ArgumentError("Empty cross-machine input token. Use label:/path/to/ladder_outdir"))
    if occursin(":", token)
        parts = split(token, ":", limit=2)
        label = strip(parts[1])
        outdir = strip(parts[2])
        isempty(outdir) && throw(ArgumentError("Cross-machine input '$raw' is missing output directory path after ':'."))
        if isempty(label)
            label = _default_input_label(outdir)
        end
        return CrossMachineInput(label=label, outdir=abspath(outdir))
    end
    return CrossMachineInput(label=_default_input_label(token), outdir=abspath(token))
end

@inline function _parse_input_list(raw::AbstractString)::Vector{CrossMachineInput}
    inputs = CrossMachineInput[]
    for token in split(String(raw), ",")
        stripped = strip(token)
        isempty(stripped) && continue
        push!(inputs, _parse_input_token(stripped))
    end
    return inputs
end

function parse_cross_machine_cli()::CrossMachineConfig
    profile = lowercase(strip(get(ENV, "SPACEAGORA_SMART_LADDER_PROFILE", get(ENV, "SPACEAGORA_PERF_PROFILE", "full"))))
    outdir = get(ENV, "SPACEAGORA_SMART_LADDER_CROSS_MACHINE_OUTDIR", DEFAULT_CROSS_MACHINE_OUTDIR)
    clean = _parse_bool_token(get(ENV, "SPACEAGORA_SMART_LADDER_CROSS_MACHINE_CLEAN", "1"))
    inputs = _parse_input_list(get(ENV, "SPACEAGORA_SMART_LADDER_MACHINE_INPUTS", ""))

    for arg in ARGS
        if arg in ("quick", "full")
            profile = arg
        elseif startswith(arg, "--profile=")
            profile = lowercase(strip(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--clean=")
            clean = _parse_bool_token(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--input=")
            push!(inputs, _parse_input_token(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--inputs=")
            append!(inputs, _parse_input_list(split(arg, "=", limit=2)[2]))
        else
            throw(ArgumentError(
                "Unknown argument '$arg'. Supported: [quick|full], --profile=..., --outdir=..., --clean=0|1, " *
                "--input=label:/path/to/ladder_outdir (repeatable), --inputs=label1:/path1,label2:/path2."
            ))
        end
    end

    if isempty(inputs)
        throw(ArgumentError(
            "No cross-machine inputs provided. Use --input=label:/path/to/smart_parallel_ladder_outdir " *
            "or set SPACEAGORA_SMART_LADDER_MACHINE_INPUTS."
        ))
    end

    seen = Set{String}()
    deduped = CrossMachineInput[]
    for input in inputs
        key = lowercase(strip(input.label))
        if key in seen
            throw(ArgumentError("Duplicate input label '$(input.label)'. Labels must be unique."))
        end
        push!(deduped, input)
        push!(seen, key)
    end

    return CrossMachineConfig(
        profile=profile,
        outdir=abspath(outdir),
        clean=clean,
        inputs=deduped
    )
end

function _latest_artifact_path(outdir::String, prefix::String, profile::String, suffix::String)::String
    pattern = "$(prefix)_$(profile)_"
    candidates = String[]
    for entry in readdir(outdir)
        if startswith(entry, pattern) && endswith(entry, suffix)
            push!(candidates, joinpath(outdir, entry))
        end
    end
    isempty(candidates) && error("No artifact found in '$outdir' for prefix='$prefix' profile='$profile' suffix='$suffix'")
    sort!(candidates; by=path -> stat(path).mtime)
    return last(candidates)
end

function _latest_artifact_path_optional(outdir::String, prefix::String, profile::String, suffix::String)::Union{Nothing, String}
    pattern = "$(prefix)_$(profile)_"
    candidates = String[]
    for entry in readdir(outdir)
        if startswith(entry, pattern) && endswith(entry, suffix)
            push!(candidates, joinpath(outdir, entry))
        end
    end
    isempty(candidates) && return nothing
    sort!(candidates; by=path -> stat(path).mtime)
    return last(candidates)
end

@inline function _has_column(df::DataFrame, col::Symbol)::Bool
    return col in Symbol.(names(df))
end

@inline function _first_nonmissing(df::DataFrame, col::Symbol, fallback::String)::String
    if !(_has_column(df, col)) || nrow(df) == 0
        return fallback
    end
    for v in df[!, col]
        if !(v isa Missing)
            return string(v)
        end
    end
    return fallback
end

@inline function _add_constant_column!(df::DataFrame, col::Symbol, value)::DataFrame
    df[!, col] = fill(value, nrow(df))
    return df
end

@inline function _concat_dataframes(parts::Vector{DataFrame})::DataFrame
    isempty(parts) && return DataFrame()
    if length(parts) == 1
        return parts[1]
    end
    out = parts[1]
    for part in parts[2:end]
        out = vcat(out, part; cols=:union)
    end
    return out
end

@inline function _safe_float(v)
    return v isa Missing ? missing : Float64(v)
end

@inline function _fmt_md(v)
    if v isa Missing
        return "n/a"
    elseif v isa Float64
        return isfinite(v) ? string(round(v; digits=6)) : "n/a"
    elseif v isa Bool
        return v ? "true" : "false"
    end
    return string(v)
end

@inline function _mode_token(v)::String
    return lowercase(strip(string(v)))
end

@inline function _is_adaptive_mode_token(v)::Bool
    token = _mode_token(v)
    return token in ("outer_inner_adaptive", "outer_inner_full_smart")
end

@inline function _is_fixed_candidate_mode_token(v)::Bool
    token = _mode_token(v)
    return token in ("outer_only_process", "outer_only", "inner_only", "outer_inner_static")
end

@inline function _safe_std(values)
    vec = collect(skipmissing(values))
    isempty(vec) && return missing
    return std(Float64.(vec); corrected=false)
end

@inline function _route_entropy_bits(none_pct, threads_pct, process_pct)
    probs = Float64[]
    for pct in (none_pct, threads_pct, process_pct)
        if pct isa Missing || !(pct isa Real)
            continue
        end
        p = clamp(Float64(pct) / 100.0, 0.0, 1.0)
        p > 0.0 && push!(probs, p)
    end
    isempty(probs) && return missing
    return -sum(p -> p * log2(p), probs)
end

function _write_markdown_table(io, df::DataFrame)
    names_vec = String.(names(df))
    println(io, "| ", join(names_vec, " | "), " |")
    println(io, "|", join(fill("---", length(names_vec)), "|"), "|")
    for row in eachrow(df)
        values = [_fmt_md(row[name]) for name in names(df)]
        println(io, "| ", join(values, " | "), " |")
    end
end

function _build_cross_machine_speedup_summary(df::DataFrame)::DataFrame
    nrow(df) == 0 && return DataFrame()
    valid = df[.!ismissing.(df.total_speedup_vs_r0), :]
    nrow(valid) == 0 && return DataFrame()
    out = combine(
        groupby(valid, [:rung, :mode]),
        :total_speedup_vs_r0 => mean => :mean_speedup_vs_r0,
        :total_speedup_vs_r0 => median => :median_speedup_vs_r0,
        :total_speedup_vs_r0 => minimum => :min_speedup_vs_r0,
        :total_speedup_vs_r0 => maximum => :max_speedup_vs_r0,
        nrow => :machines
    )
    sort!(out, :median_speedup_vs_r0, rev=true)
    return out
end

function _build_best_rung_table(df::DataFrame)::DataFrame
    nrow(df) == 0 && return DataFrame()
    rows = NamedTuple[]
    for grp in groupby(df, :machine_tag)
        valid = grp[
            (lowercase.(string.(grp.mode)) .!= "serial") .&
            .!ismissing.(grp.total_speedup_vs_r0), :
        ]
        nrow(valid) == 0 && continue
        sort!(valid, :total_speedup_vs_r0, rev=true)
        best = valid[1, :]
        speedups = Float64.(coalesce.(valid.total_speedup_vs_r0, NaN))
        speedups = [x for x in speedups if isfinite(x)]
        isempty(speedups) && continue
        best_speedup = first(speedups)
        worst_speedup = minimum(speedups)
        push!(rows, (
            machine_tag=String(best.machine_tag),
            machine_label=String(best.machine_label),
            hardware_class=String(best.hardware_class),
            best_rung=String(best.rung),
            best_mode=String(best.mode),
            best_speedup_vs_r0=best_speedup,
            regret_to_worst_speedup=best_speedup - worst_speedup,
            candidate_rungs=length(speedups)
        ))
    end
    out = DataFrame(rows)
    if nrow(out) > 0
        sort!(out, :best_speedup_vs_r0, rev=true)
    end
    return out
end

function _build_mission_family_summary(df::DataFrame)::DataFrame
    if nrow(df) == 0 || !(_has_column(df, :median_speedup_vs_r0))
        return DataFrame()
    end
    valid = df[.!ismissing.(df.median_speedup_vs_r0), :]
    nrow(valid) == 0 && return DataFrame()
    out = combine(
        groupby(valid, [:mission_family, :rung, :mode]),
        :median_speedup_vs_r0 => mean => :mean_median_speedup_vs_r0,
        :median_speedup_vs_r0 => median => :median_of_medians_speedup_vs_r0,
        :median_speedup_vs_r0 => minimum => :min_median_speedup_vs_r0,
        :median_speedup_vs_r0 => maximum => :max_median_speedup_vs_r0,
        nrow => :machines
    )
    sort!(out, [:mission_family, :median_of_medians_speedup_vs_r0], rev=[false, true])
    return out
end

function _build_route_mix_summary(df::DataFrame)::DataFrame
    nrow(df) == 0 && return DataFrame()
    needed = (:none_pct, :threads_pct, :process_pct)
    for col in needed
        _has_column(df, col) || return DataFrame()
    end
    out = combine(
        groupby(df, [:rung, :mode]),
        :none_pct => mean => :none_pct_mean,
        :threads_pct => mean => :threads_pct_mean,
        :process_pct => mean => :process_pct_mean,
        nrow => :machines
    )
    sort!(out, [:rung, :mode])
    return out
end

function _build_adaptive_speedup_table(df::DataFrame)::DataFrame
    nrow(df) == 0 && return DataFrame()
    mask = [_is_adaptive_mode_token(v) for v in df.mode]
    out = df[mask, :]
    nrow(out) == 0 && return DataFrame()
    keep_cols = Symbol[
        :machine_tag,
        :machine_label,
        :hardware_class,
        :rung,
        :mode,
        :total_elapsed_s,
        :total_speedup_vs_r0,
        :run_benchmarks_speedup_vs_r0,
        :run_per_orbit_speedup_vs_r0
    ]
    keep = [c for c in keep_cols if _has_column(out, c)]
    out = select(out, keep)
    sort!(out, [:machine_tag, :mode])
    return out
end

function _build_adaptive_speedup_summary(df::DataFrame)::DataFrame
    nrow(df) == 0 && return DataFrame()
    valid = df[.!ismissing.(df.total_speedup_vs_r0), :]
    nrow(valid) == 0 && return DataFrame()
    out = combine(
        groupby(valid, [:rung, :mode]),
        :total_speedup_vs_r0 => mean => :mean_speedup_vs_r0,
        :total_speedup_vs_r0 => median => :median_speedup_vs_r0,
        :total_speedup_vs_r0 => minimum => :min_speedup_vs_r0,
        :total_speedup_vs_r0 => maximum => :max_speedup_vs_r0,
        nrow => :machines
    )
    sort!(out, :median_speedup_vs_r0, rev=true)
    return out
end

function _build_adaptive_regret_table(speedup_df::DataFrame)::DataFrame
    nrow(speedup_df) == 0 && return DataFrame()
    rows = NamedTuple[]
    for grp in groupby(speedup_df, :machine_tag)
        fixed_mask = [_is_fixed_candidate_mode_token(v) for v in grp.mode]
        fixed = grp[fixed_mask, :]
        fixed = fixed[.!ismissing.(fixed.total_speedup_vs_r0) .& .!ismissing.(fixed.total_elapsed_s), :]
        nrow(fixed) == 0 && continue
        sort!(fixed, :total_speedup_vs_r0, rev=true)
        best_fixed = fixed[1, :]
        best_fixed_speedup = Float64(best_fixed.total_speedup_vs_r0)
        best_fixed_elapsed = Float64(best_fixed.total_elapsed_s)

        adaptive_mask = [_is_adaptive_mode_token(v) for v in grp.mode]
        adaptive = grp[adaptive_mask, :]
        for row in eachrow(adaptive)
            if row.total_speedup_vs_r0 isa Missing || row.total_elapsed_s isa Missing
                continue
            end
            adaptive_speedup = Float64(row.total_speedup_vs_r0)
            adaptive_elapsed = Float64(row.total_elapsed_s)
            speedup_regret = max(0.0, best_fixed_speedup - adaptive_speedup)
            time_regret_pct = if isfinite(best_fixed_elapsed) && best_fixed_elapsed > 0.0
                max(0.0, 100.0 * (adaptive_elapsed / best_fixed_elapsed - 1.0))
            else
                missing
            end
            push!(rows, (
                machine_tag=String(row.machine_tag),
                machine_label=String(row.machine_label),
                hardware_class=String(row.hardware_class),
                adaptive_rung=String(row.rung),
                adaptive_mode=String(row.mode),
                adaptive_speedup_vs_r0=adaptive_speedup,
                adaptive_elapsed_s=adaptive_elapsed,
                best_fixed_rung=String(best_fixed.rung),
                best_fixed_mode=String(best_fixed.mode),
                best_fixed_speedup_vs_r0=best_fixed_speedup,
                best_fixed_elapsed_s=best_fixed_elapsed,
                speedup_regret_vs_best_fixed=speedup_regret,
                time_regret_vs_best_fixed_pct=time_regret_pct,
                adaptive_beats_best_fixed=adaptive_speedup >= best_fixed_speedup
            ))
        end
    end
    out = DataFrame(rows)
    if nrow(out) > 0
        sort!(out, [:machine_tag, :adaptive_mode])
    end
    return out
end

function _build_adaptive_regret_summary(df::DataFrame)::DataFrame
    nrow(df) == 0 && return DataFrame()
    out = combine(
        groupby(df, [:adaptive_rung, :adaptive_mode]),
        :speedup_regret_vs_best_fixed => mean => :mean_speedup_regret,
        :speedup_regret_vs_best_fixed => median => :median_speedup_regret,
        :speedup_regret_vs_best_fixed => maximum => :max_speedup_regret,
        :time_regret_vs_best_fixed_pct => mean => :mean_time_regret_pct,
        :time_regret_vs_best_fixed_pct => maximum => :max_time_regret_pct,
        :adaptive_beats_best_fixed => (v -> 100.0 * mean(Bool.(v))) => :win_rate_pct,
        nrow => :machines
    )
    sort!(out, [:win_rate_pct, :mean_speedup_regret], rev=[true, false])
    return out
end

function _build_adaptive_route_distribution_table(route_mix_df::DataFrame)::DataFrame
    nrow(route_mix_df) == 0 && return DataFrame()
    needed = (:none_pct, :threads_pct, :process_pct, :mode, :rung, :machine_tag, :machine_label, :hardware_class)
    for col in needed
        _has_column(route_mix_df, col) || return DataFrame()
    end
    rows = NamedTuple[]
    for row in eachrow(route_mix_df)
        _is_adaptive_mode_token(row.mode) || continue
        none_pct = _safe_float(row.none_pct)
        threads_pct = _safe_float(row.threads_pct)
        process_pct = _safe_float(row.process_pct)
        route_pairs = [
            ("none", ismissing(none_pct) ? -Inf : none_pct),
            ("threads", ismissing(threads_pct) ? -Inf : threads_pct),
            ("process", ismissing(process_pct) ? -Inf : process_pct)
        ]
        sort!(route_pairs; by=x -> x[2], rev=true)
        dominant_route, dominant_pct = route_pairs[1]
        push!(rows, (
            machine_tag=String(row.machine_tag),
            machine_label=String(row.machine_label),
            hardware_class=String(row.hardware_class),
            rung=String(row.rung),
            mode=String(row.mode),
            none_pct=none_pct,
            threads_pct=threads_pct,
            process_pct=process_pct,
            dominant_route=dominant_route,
            dominant_route_pct=(isfinite(dominant_pct) ? dominant_pct : missing),
            route_entropy_bits=_route_entropy_bits(none_pct, threads_pct, process_pct)
        ))
    end
    out = DataFrame(rows)
    if nrow(out) > 0
        sort!(out, [:machine_tag, :mode])
    end
    return out
end

function _build_adaptive_route_stability_summary(df::DataFrame)::DataFrame
    nrow(df) == 0 && return DataFrame()
    out = combine(
        groupby(df, [:rung, :mode]),
        :dominant_route_pct => mean => :mean_dominant_route_pct,
        :dominant_route_pct => minimum => :min_dominant_route_pct,
        :dominant_route_pct => maximum => :max_dominant_route_pct,
        :route_entropy_bits => mean => :mean_route_entropy_bits,
        :route_entropy_bits => minimum => :min_route_entropy_bits,
        :route_entropy_bits => maximum => :max_route_entropy_bits,
        :none_pct => mean => :none_pct_mean,
        :threads_pct => mean => :threads_pct_mean,
        :process_pct => mean => :process_pct_mean,
        :none_pct => _safe_std => :none_pct_std,
        :threads_pct => _safe_std => :threads_pct_std,
        :process_pct => _safe_std => :process_pct_std,
        nrow => :machines
    )
    sort!(out, :mean_dominant_route_pct, rev=true)
    return out
end

function _build_machine_class_coverage_table(manifest_df::DataFrame)::DataFrame
    nrow(manifest_df) == 0 && return DataFrame()
    classes = sort(unique(lowercase.(string.(manifest_df.hardware_class))))
    required = ["small", "medium", "large"]
    missing = [c for c in required if !(c in classes)]
    return DataFrame([
        (
            machines=nrow(manifest_df),
            classes_present=join(classes, ","),
            required_classes=join(required, ","),
            classes_missing=isempty(missing) ? "none" : join(missing, ","),
            coverage_ok=isempty(missing)
        )
    ])
end

function _write_cross_machine_report(
    path::String,
    config::CrossMachineConfig,
    manifest_df::DataFrame,
    machine_coverage_df::DataFrame,
    speedup_df::DataFrame,
    speedup_summary_df::DataFrame,
    adaptive_speedup_df::DataFrame,
    adaptive_speedup_summary_df::DataFrame,
    adaptive_regret_df::DataFrame,
    adaptive_regret_summary_df::DataFrame,
    adaptive_route_distribution_df::DataFrame,
    adaptive_route_stability_df::DataFrame,
    best_rung_df::DataFrame,
    mission_family_summary_df::DataFrame,
    route_mix_summary_df::DataFrame;
    manifest_csv_path::String="",
    machine_coverage_csv_path::String="",
    speedup_csv_path::String="",
    speedup_summary_csv_path::String="",
    adaptive_speedup_csv_path::String="",
    adaptive_speedup_summary_csv_path::String="",
    adaptive_regret_csv_path::String="",
    adaptive_regret_summary_csv_path::String="",
    adaptive_route_distribution_csv_path::String="",
    adaptive_route_stability_csv_path::String="",
    best_rung_csv_path::String="",
    mission_family_summary_csv_path::String="",
    route_mix_summary_csv_path::String=""
)
    generated = string(now(UTC))
    open(path, "w") do io
        println(io, "# Smart Parallel Ladder Cross-Machine Replay")
        println(io)
        println(io, "- Generated (UTC): $generated")
        println(io, "- Profile: `$(config.profile)`")
        println(io, "- Machines loaded: `$(nrow(manifest_df))`")
        println(io)

        println(io, "## Replay Inputs")
        println(io)
        _write_markdown_table(io, manifest_df)
        println(io)

        println(io, "## Machine-Class Coverage")
        println(io)
        if nrow(machine_coverage_df) > 0
            _write_markdown_table(io, machine_coverage_df)
        else
            println(io, "- No machine coverage rows available.")
        end
        println(io)

        println(io, "## Adaptive Speedup vs Serial By Machine")
        println(io)
        if nrow(adaptive_speedup_df) > 0
            _write_markdown_table(io, adaptive_speedup_df)
        else
            println(io, "- No adaptive speedup rows available.")
        end
        println(io)

        println(io, "## Adaptive Speedup Summary Across Machines")
        println(io)
        if nrow(adaptive_speedup_summary_df) > 0
            _write_markdown_table(io, adaptive_speedup_summary_df)
        else
            println(io, "- No adaptive speedup summary rows available.")
        end
        println(io)

        println(io, "## Chosen Route Distribution (Adaptive Modes)")
        println(io)
        if nrow(adaptive_route_distribution_df) > 0
            _write_markdown_table(io, adaptive_route_distribution_df)
        else
            println(io, "- No adaptive route-distribution rows available.")
        end
        println(io)

        println(io, "## Adaptive-Choice Stability")
        println(io)
        if nrow(adaptive_route_stability_df) > 0
            _write_markdown_table(io, adaptive_route_stability_df)
        else
            println(io, "- No adaptive stability rows available.")
        end
        println(io)

        println(io, "## Regret vs Best Fixed Strategy")
        println(io)
        if nrow(adaptive_regret_df) > 0
            _write_markdown_table(io, adaptive_regret_df)
        else
            println(io, "- No adaptive regret rows available.")
        end
        println(io)

        println(io, "## Adaptive Regret Summary")
        println(io)
        if nrow(adaptive_regret_summary_df) > 0
            _write_markdown_table(io, adaptive_regret_summary_df)
        else
            println(io, "- No adaptive regret summary rows available.")
        end
        println(io)

        println(io, "## Best Rung By Machine")
        println(io)
        if nrow(best_rung_df) > 0
            _write_markdown_table(io, best_rung_df)
        else
            println(io, "- No best-rung rows available.")
        end
        println(io)

        println(io, "## Rung Speedup Summary Across Machines")
        println(io)
        if nrow(speedup_summary_df) > 0
            _write_markdown_table(io, speedup_summary_df)
        else
            println(io, "- No speedup summary rows available.")
        end
        println(io)

        println(io, "## Mission-Family Summary Across Machines")
        println(io)
        if nrow(mission_family_summary_df) > 0
            _write_markdown_table(io, mission_family_summary_df)
        else
            println(io, "- No mission-family summary rows available.")
        end
        println(io)

        println(io, "## Route Mix Summary Across Machines")
        println(io)
        if nrow(route_mix_summary_df) > 0
            _write_markdown_table(io, route_mix_summary_df)
        else
            println(io, "- No route-mix summary rows available.")
        end
        println(io)

        println(io, "## Output Files")
        println(io)
        println(io, "- Manifest CSV: `$(manifest_csv_path)`")
        println(io, "- Machine coverage CSV: `$(machine_coverage_csv_path)`")
        println(io, "- Combined speedup CSV: `$(speedup_csv_path)`")
        println(io, "- Speedup summary CSV: `$(speedup_summary_csv_path)`")
        println(io, "- Adaptive speedup CSV: `$(adaptive_speedup_csv_path)`")
        println(io, "- Adaptive speedup summary CSV: `$(adaptive_speedup_summary_csv_path)`")
        println(io, "- Adaptive regret CSV: `$(adaptive_regret_csv_path)`")
        println(io, "- Adaptive regret summary CSV: `$(adaptive_regret_summary_csv_path)`")
        println(io, "- Adaptive route distribution CSV: `$(adaptive_route_distribution_csv_path)`")
        println(io, "- Adaptive route stability CSV: `$(adaptive_route_stability_csv_path)`")
        println(io, "- Best rung CSV: `$(best_rung_csv_path)`")
        println(io, "- Mission-family summary CSV: `$(mission_family_summary_csv_path)`")
        println(io, "- Route-mix summary CSV: `$(route_mix_summary_csv_path)`")
        println(io)

        println(io, "## Reproducibility")
        println(io, "```bash")
        input_flags = join(
            ["--input=$(row.input_tag):$(row.input_outdir)" for row in eachrow(manifest_df)],
            " "
        )
        println(
            io,
            "julia --project=.AGORA test/performance_smart_parallel_ladder_cross_machine.jl " *
            "--profile=$(config.profile) --outdir=$(config.outdir) --clean=$(config.clean ? 1 : 0) $(input_flags)"
        )
        println(io, "```")
    end
    return nothing
end

function main_cross_machine_replay()
    config = parse_cross_machine_cli()
    if config.clean
        rm(config.outdir; recursive=true, force=true)
    end
    mkpath(config.outdir)

    println("Cross-machine ladder profile: $(config.profile)")
    println("Outdir: $(config.outdir)")
    input_tokens = [string(i.label, "=>", i.outdir) for i in config.inputs]
    println("Inputs: $(join(input_tokens, ", "))")

    manifest_rows = NamedTuple[]
    speedup_parts = DataFrame[]
    mission_family_parts = DataFrame[]
    route_mix_parts = DataFrame[]

    for input in config.inputs
        outdir = input.outdir
        speedup_path = _latest_artifact_path(outdir, "smart_parallel_ladder_speedup_vs_r0", config.profile, ".csv")
        overview_path = _latest_artifact_path(outdir, "smart_parallel_ladder_mode_overview", config.profile, ".csv")
        mission_family_path = _latest_artifact_path_optional(outdir, "smart_parallel_ladder_mission_family_speedup", config.profile, ".csv")
        route_mix_path = _latest_artifact_path_optional(outdir, "smart_parallel_ladder_route_mix", config.profile, ".csv")

        speedup_df = CSV.read(speedup_path, DataFrame)
        overview_df = CSV.read(overview_path, DataFrame)
        machine_label = _first_nonmissing(overview_df, :machine_label, input.label)
        hardware_class = _first_nonmissing(overview_df, :hardware_class, "unknown")

        _add_constant_column!(speedup_df, :machine_tag, input.label)
        _add_constant_column!(speedup_df, :machine_label, machine_label)
        _add_constant_column!(speedup_df, :hardware_class, hardware_class)
        _add_constant_column!(speedup_df, :input_outdir, outdir)
        push!(speedup_parts, speedup_df)

        if !(mission_family_path === nothing)
            mission_family_df = CSV.read(mission_family_path, DataFrame)
            _add_constant_column!(mission_family_df, :machine_tag, input.label)
            _add_constant_column!(mission_family_df, :machine_label, machine_label)
            _add_constant_column!(mission_family_df, :hardware_class, hardware_class)
            _add_constant_column!(mission_family_df, :input_outdir, outdir)
            push!(mission_family_parts, mission_family_df)
        end

        if !(route_mix_path === nothing)
            route_mix_df = CSV.read(route_mix_path, DataFrame)
            _add_constant_column!(route_mix_df, :machine_tag, input.label)
            _add_constant_column!(route_mix_df, :machine_label, machine_label)
            _add_constant_column!(route_mix_df, :hardware_class, hardware_class)
            _add_constant_column!(route_mix_df, :input_outdir, outdir)
            push!(route_mix_parts, route_mix_df)
        end

        push!(manifest_rows, (
            input_tag=input.label,
            input_outdir=outdir,
            machine_label=machine_label,
            hardware_class=hardware_class,
            speedup_csv=speedup_path,
            overview_csv=overview_path,
            mission_family_csv=coalesce(mission_family_path, ""),
            route_mix_csv=coalesce(route_mix_path, "")
        ))
    end

    manifest_df = DataFrame(manifest_rows)
    speedup_df = _concat_dataframes(speedup_parts)
    mission_family_df = _concat_dataframes(mission_family_parts)
    route_mix_df = _concat_dataframes(route_mix_parts)

    machine_coverage_df = _build_machine_class_coverage_table(manifest_df)
    speedup_summary_df = _build_cross_machine_speedup_summary(speedup_df)
    adaptive_speedup_df = _build_adaptive_speedup_table(speedup_df)
    adaptive_speedup_summary_df = _build_adaptive_speedup_summary(adaptive_speedup_df)
    adaptive_regret_df = _build_adaptive_regret_table(speedup_df)
    adaptive_regret_summary_df = _build_adaptive_regret_summary(adaptive_regret_df)
    adaptive_route_distribution_df = _build_adaptive_route_distribution_table(route_mix_df)
    adaptive_route_stability_df = _build_adaptive_route_stability_summary(adaptive_route_distribution_df)
    best_rung_df = _build_best_rung_table(speedup_df)
    mission_family_summary_df = _build_mission_family_summary(mission_family_df)
    route_mix_summary_df = _build_route_mix_summary(route_mix_df)

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    manifest_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_manifest_$(config.profile)_$(stamp).csv")
    machine_coverage_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_machine_coverage_$(config.profile)_$(stamp).csv")
    speedup_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_speedup_$(config.profile)_$(stamp).csv")
    speedup_summary_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_speedup_summary_$(config.profile)_$(stamp).csv")
    adaptive_speedup_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_adaptive_speedup_$(config.profile)_$(stamp).csv")
    adaptive_speedup_summary_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_adaptive_speedup_summary_$(config.profile)_$(stamp).csv")
    adaptive_regret_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_adaptive_regret_$(config.profile)_$(stamp).csv")
    adaptive_regret_summary_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_adaptive_regret_summary_$(config.profile)_$(stamp).csv")
    adaptive_route_distribution_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_adaptive_route_distribution_$(config.profile)_$(stamp).csv")
    adaptive_route_stability_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_adaptive_route_stability_$(config.profile)_$(stamp).csv")
    best_rung_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_best_rung_$(config.profile)_$(stamp).csv")
    mission_family_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_mission_family_$(config.profile)_$(stamp).csv")
    mission_family_summary_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_mission_family_summary_$(config.profile)_$(stamp).csv")
    route_mix_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_route_mix_$(config.profile)_$(stamp).csv")
    route_mix_summary_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_route_mix_summary_$(config.profile)_$(stamp).csv")
    report_path = joinpath(config.outdir, "smart_parallel_ladder_cross_machine_report_$(config.profile)_$(stamp).md")

    CSV.write(manifest_path, manifest_df)
    CSV.write(machine_coverage_path, machine_coverage_df)
    CSV.write(speedup_path, speedup_df)
    CSV.write(speedup_summary_path, speedup_summary_df)
    CSV.write(adaptive_speedup_path, adaptive_speedup_df)
    CSV.write(adaptive_speedup_summary_path, adaptive_speedup_summary_df)
    CSV.write(adaptive_regret_path, adaptive_regret_df)
    CSV.write(adaptive_regret_summary_path, adaptive_regret_summary_df)
    CSV.write(adaptive_route_distribution_path, adaptive_route_distribution_df)
    CSV.write(adaptive_route_stability_path, adaptive_route_stability_df)
    CSV.write(best_rung_path, best_rung_df)
    CSV.write(mission_family_path, mission_family_df)
    CSV.write(mission_family_summary_path, mission_family_summary_df)
    CSV.write(route_mix_path, route_mix_df)
    CSV.write(route_mix_summary_path, route_mix_summary_df)

    _write_cross_machine_report(
        report_path,
        config,
        manifest_df,
        machine_coverage_df,
        speedup_df,
        speedup_summary_df,
        adaptive_speedup_df,
        adaptive_speedup_summary_df,
        adaptive_regret_df,
        adaptive_regret_summary_df,
        adaptive_route_distribution_df,
        adaptive_route_stability_df,
        best_rung_df,
        mission_family_summary_df,
        route_mix_summary_df;
        manifest_csv_path=manifest_path,
        machine_coverage_csv_path=machine_coverage_path,
        speedup_csv_path=speedup_path,
        speedup_summary_csv_path=speedup_summary_path,
        adaptive_speedup_csv_path=adaptive_speedup_path,
        adaptive_speedup_summary_csv_path=adaptive_speedup_summary_path,
        adaptive_regret_csv_path=adaptive_regret_path,
        adaptive_regret_summary_csv_path=adaptive_regret_summary_path,
        adaptive_route_distribution_csv_path=adaptive_route_distribution_path,
        adaptive_route_stability_csv_path=adaptive_route_stability_path,
        best_rung_csv_path=best_rung_path,
        mission_family_summary_csv_path=mission_family_summary_path,
        route_mix_summary_csv_path=route_mix_summary_path
    )

    println()
    println("Cross-machine ladder replay summary complete.")
    println("manifest: $(manifest_path)")
    println("machine coverage: $(machine_coverage_path)")
    println("combined speedup: $(speedup_path)")
    println("speedup summary: $(speedup_summary_path)")
    println("adaptive speedup: $(adaptive_speedup_path)")
    println("adaptive speedup summary: $(adaptive_speedup_summary_path)")
    println("adaptive regret: $(adaptive_regret_path)")
    println("adaptive regret summary: $(adaptive_regret_summary_path)")
    println("adaptive route distribution: $(adaptive_route_distribution_path)")
    println("adaptive route stability: $(adaptive_route_stability_path)")
    println("best rung: $(best_rung_path)")
    println("mission-family summary: $(mission_family_summary_path)")
    println("route-mix summary: $(route_mix_summary_path)")
    println("report: $(report_path)")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_cross_machine_replay()
end
