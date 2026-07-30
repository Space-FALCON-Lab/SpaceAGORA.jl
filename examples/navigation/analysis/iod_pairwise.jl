using CSV
using DataFrames
using Statistics

if !isdefined(@__MODULE__, :NavigationPaths)
    include(joinpath(@__DIR__, "..", "paths.jl"))
end
using .NavigationPaths: navigation_output_path, resolve_repo_path

const NAVIGATION_RATE_SEC = 5.0
const NAVIGATION_DT_TOL_SEC = 0.25

@inline function _percentage(numerator::Integer, denominator::Integer)::Float64
    denominator == 0 && return NaN
    return 100.0 * numerator / denominator
end

function _group_table(pairwise::DataFrame)::DataFrame
    keys = [:time_s, :reference_observer, :reference_slot]
    return combine(
        groupby(pairwise, keys),
        :composition => first => :composition,
        :group_size => first => :group_size,
        :pair_count => first => :pair_count,
        :covariance_gate_passed => first => :covariance_gate_passed,
        :all_pairs_pass => first => :all_pairs_pass,
        :pair_passes => (values -> count(!, values)) => :violating_pairs,
        :maximum_miss_m => maximum => :largest_pairwise_miss_m
    )
end

function _mark_promoted!(groups::DataFrame, diagnostics::DataFrame)::Nothing
    promotion_times = Dict{Tuple{Int, Int}, Vector{Tuple{Float64, String}}}()
    for row in eachrow(diagnostics)
        row.outcome == "promoted" || continue
        row.stage in ("next_step_validation", "direct_promotion") || continue
        key = (Int(row.observer), Int(row.slot))
        push!(get!(promotion_times, key, Tuple{Float64, String}[]), (Float64(row.time_s), String(row.stage)))
    end

    promoted = falses(nrow(groups))
    for (index, group) in enumerate(eachrow(groups))
        key = (Int(group.reference_observer), Int(group.reference_slot))
        for (promotion_time_s, stage) in get(promotion_times, key, Tuple{Float64, String}[])
            expected_delay_s = stage == "direct_promotion" ? 0.0 : NAVIGATION_RATE_SEC
            if abs(promotion_time_s - Float64(group.time_s) - expected_delay_s) <= NAVIGATION_DT_TOL_SEC
                promoted[index] = true
                break
            end
        end
    end
    groups.promoted = promoted
    return nothing
end

function _discover_run_directories(root::String)
    discovered = Tuple{Int, String}[]
    for (directory, _, files) in walkdir(root)
        "iod_pairwise_consistency_table.csv" in files || continue
        "iod_diagnostics_table.csv" in files || continue
        basename(directory) == "proposed" || continue
        run_match = match(r"run_(\d+)", directory)
        run_match === nothing && continue
        push!(discovered, (parse(Int, run_match.captures[1]), directory))
    end
    sort!(discovered; by=first)
    return discovered
end

function _category_counts(groups::DataFrame, mask::AbstractVector{Bool})
    selected = groups[mask, :]
    total = nrow(selected)
    all_pass = total == 0 ? 0 : count(selected.all_pairs_pass)
    return (
        total=total,
        all_pass=all_pass,
        violating=total - all_pass,
        violating_pairs=total == 0 ? 0 : sum(selected.violating_pairs),
        compliance_pct=_percentage(all_pass, total),
        maximum_miss_m=total == 0 ? NaN : maximum(selected.largest_pairwise_miss_m)
    )
end

function _pooled_row(label::String, run_summary::DataFrame, prefix::String)
    total = sum(run_summary[!, Symbol(prefix * "_total")])
    all_pass = sum(run_summary[!, Symbol(prefix * "_all_pass")])
    violating = sum(run_summary[!, Symbol(prefix * "_violating")])
    violating_pairs = sum(run_summary[!, Symbol(prefix * "_violating_pairs")])
    run_compliance = [
        Float64(value) for value in run_summary[!, Symbol(prefix * "_compliance_pct")]
        if isfinite(Float64(value))
    ]
    return (
        category=label,
        groups=total,
        all_pairs_compliant=all_pass,
        groups_with_violations=violating,
        violating_pairs=violating_pairs,
        pooled_compliance_pct=_percentage(all_pass, total),
        run_compliance_median_pct=isempty(run_compliance) ? NaN : median(run_compliance),
        run_compliance_p05_pct=isempty(run_compliance) ? NaN : quantile(run_compliance, 0.05),
        run_compliance_p95_pct=isempty(run_compliance) ? NaN : quantile(run_compliance, 0.95)
    )
end

function main()
    campaign_root = isempty(ARGS) ?
        navigation_output_path("iod_pairwise") :
        resolve_repo_path(ARGS[1])
    output_root = length(ARGS) >= 2 ? resolve_repo_path(ARGS[2]) : campaign_root
    mkpath(output_root)

    run_directories = _discover_run_directories(campaign_root)
    isempty(run_directories) && error("No completed pairwise-diagnostic runs found under $(campaign_root)")

    run_rows = NamedTuple[]
    pair_type_rows = NamedTuple[]
    violating_pair_tables = DataFrame[]

    for (realization, directory) in run_directories
        pairwise = CSV.read(joinpath(directory, "iod_pairwise_consistency_table.csv"), DataFrame)
        diagnostics = CSV.read(joinpath(directory, "iod_diagnostics_table.csv"), DataFrame)
        groups = _group_table(pairwise)
        _mark_promoted!(groups, diagnostics)

        cov_pass = groups.covariance_gate_passed
        same = groups.composition .== "same_target"
        mixed = .!same
        promoted = groups.promoted
        categories = (
            all=_category_counts(groups, trues(nrow(groups))),
            same_cov=_category_counts(groups, cov_pass .& same),
            mixed_cov=_category_counts(groups, cov_pass .& mixed),
            same_promoted=_category_counts(groups, promoted .& same),
            mixed_promoted=_category_counts(groups, promoted .& mixed)
        )
        row_values = Dict{Symbol, Any}(:realization => realization)
        for (prefix, values) in pairs(categories), field in keys(values)
            row_values[Symbol(String(prefix) * "_" * String(field))] = getfield(values, field)
        end
        push!(run_rows, (; row_values...))

        for (pair_type, mask) in (
            ("reference_external", pairwise.is_reference_pair),
            ("external_external", .!pairwise.is_reference_pair)
        )
            total = count(mask)
            violations = count(.!pairwise.pair_passes[mask])
            push!(pair_type_rows, (
                realization=realization,
                pair_type=pair_type,
                pairs=total,
                compliant_pairs=total - violations,
                violating_pairs=violations,
                compliance_pct=_percentage(total - violations, total)
            ))
        end

        violating = pairwise[.!pairwise.pair_passes, :]
        if nrow(violating) > 0
            insertcols!(violating, 1, :realization => fill(realization, nrow(violating)))
            push!(violating_pair_tables, violating)
        end
    end

    run_summary = DataFrame(run_rows)
    sort!(run_summary, :realization)
    pooled = DataFrame([
        _pooled_row("All evaluated IOD groups", run_summary, "all"),
        _pooled_row("Same-target groups after covariance gate", run_summary, "same_cov"),
        _pooled_row("Mixed-target groups after covariance gate", run_summary, "mixed_cov"),
        _pooled_row("Promoted same-target initializations", run_summary, "same_promoted"),
        _pooled_row("Promoted mixed-target initializations", run_summary, "mixed_promoted")
    ])
    pair_types = DataFrame(pair_type_rows)
    pooled_pair_types = combine(
        groupby(pair_types, :pair_type),
        :pairs => sum => :pairs,
        :compliant_pairs => sum => :compliant_pairs,
        :violating_pairs => sum => :violating_pairs
    )
    pooled_pair_types.compliance_pct = 100.0 .* pooled_pair_types.compliant_pairs ./ pooled_pair_types.pairs
    violations = isempty(violating_pair_tables) ?
        DataFrame() :
        vcat(violating_pair_tables...; cols=:union)

    CSV.write(joinpath(output_root, "mc_iod_pairwise_run_summary.csv"), run_summary)
    CSV.write(joinpath(output_root, "mc_iod_pairwise_pooled_summary.csv"), pooled)
    CSV.write(joinpath(output_root, "mc_iod_pairwise_pair_type_summary.csv"), pooled_pair_types)
    CSV.write(joinpath(output_root, "mc_iod_pairwise_violations.csv"), violations)

    println("IOD pairwise campaign analysis")
    println("  completed realizations: $(nrow(run_summary))")
    show(stdout, MIME("text/plain"), pooled; allrows=true, allcols=true)
    println()
    show(stdout, MIME("text/plain"), pooled_pair_types; allrows=true, allcols=true)
    println()
end

main()
