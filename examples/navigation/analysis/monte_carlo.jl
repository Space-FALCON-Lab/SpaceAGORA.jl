using CSV
using DataFrames
using Printf
using Statistics

if !isdefined(@__MODULE__, :NavigationPaths)
    include(joinpath(@__DIR__, "..", "paths.jl"))
end
using .NavigationPaths: env_path

const OUTPUT_ROOT = env_path(
    "SPACEAGORA_MC_OUTPUT",
    "nominal"
)
const SUMMARY_PATH = joinpath(OUTPUT_ROOT, "mc_summary.csv")
const STATUS_PATH = joinpath(OUTPUT_ROOT, "mc_run_status.csv")
const WIDE_METRICS_PATH = joinpath(OUTPUT_ROOT, "mc_metrics_wide.csv")
const IOD_DIAGNOSTICS_PATH = joinpath(OUTPUT_ROOT, "mc_iod_diagnostics.csv")
const PATHOLOGY_SUMMARY_PATH = joinpath(OUTPUT_ROOT, "mc_pathology_event_summary.csv")
const PATHOLOGICAL_RUNS_PATH = joinpath(OUTPUT_ROOT, "mc_pathological_runs.csv")
const PATHOLOGY_OUTCOMES_PATH = joinpath(OUTPUT_ROOT, "mc_pathology_conditioned_outcomes.csv")
const CROWDING_SUMMARY_PATH = joinpath(OUTPUT_ROOT, "mc_crowding_pathology_summary.csv")
const IOD_EVENT_SUMMARY_PATH = joinpath(OUTPUT_ROOT, "mc_iod_event_characterization.csv")

const CASE_ORDER = [
    "proposed",
    "no_da",
    "centralized_oracle",
    "independent_local_da",
    "distributed_oracle_da",
    "baseline_da_theta_0p01",
    "baseline_da_theta_0p05",
    "baseline_da"
]

const CASE_LABELS = Dict(
    "proposed" => "Proposed",
    "no_da" => "No DA",
    "centralized_oracle" => "Centralized oracle",
    "independent_local_da" => "Independent local",
    "distributed_oracle_da" => "Distributed oracle",
    "baseline_da_theta_0p01" => "Baseline (0.01 rad)",
    "baseline_da_theta_0p05" => "Baseline (0.05 rad)",
    "baseline_da" => "Baseline"
)

@inline function _metric_row(
    summary::DataFrame,
    case::String,
    section::String,
    metric::String
)
    rows = summary[
        (summary.case .== case) .&
        (summary.section .== section) .&
        (summary.metric .== metric),
        :
    ]
    return nrow(rows) == 1 ? rows[1, :] : nothing
end

@inline function _stat(
    summary::DataFrame,
    case::String,
    section::String,
    metric::String,
    statistic::Symbol=:median
)::Float64
    row = _metric_row(summary, case, section, metric)
    row === nothing && return NaN
    value = row[statistic]
    return ismissing(value) ? NaN : Float64(value)
end

@inline _case_label(case::String) = get(CASE_LABELS, case, replace(case, '_' => ' '))

function _available_cases(summary::DataFrame)::Vector{String}
    present = Set(String.(unique(summary.case)))
    ordered = [case for case in CASE_ORDER if case in present]
    append!(ordered, sort([case for case in present if !(case in Set(CASE_ORDER))]))
    return ordered
end

function _completed_counts(status::DataFrame)::Dict{String, Int}
    counts = Dict{String, Int}()
    nrow(status) == 0 && return counts
    for row in eachrow(status)
        String(row.status) in ("ok", "resumed") || continue
        case = String(row.case)
        counts[case] = get(counts, case, 0) + 1
    end
    return counts
end

function _fmt(value::Float64; digits::Int=4)::String
    isfinite(value) || return "n/a"
    abs(value - round(value)) <= 1e-10 && return string(Int(round(value)))
    return @sprintf("%.*g", digits, value)
end

function _interval_text(summary, case, section, metric; digits=4)::String
    med = _stat(summary, case, section, metric, :median)
    p05 = _stat(summary, case, section, metric, :p05)
    p95 = _stat(summary, case, section, metric, :p95)
    isfinite(med) || return "n/a"
    return "$(_fmt(med; digits=digits)) [$(_fmt(p05; digits=digits)), $(_fmt(p95; digits=digits))]"
end

function _print_table(title::String, headers::Vector{String}, rows::Vector{Vector{String}})
    println()
    println(title)
    widths = [length(header) for header in headers]
    for row in rows, col in eachindex(headers)
        widths[col] = max(widths[col], length(row[col]))
    end
    println(join([rpad(headers[col], widths[col]) for col in eachindex(headers)], "  "))
    println(join([repeat("-", widths[col]) for col in eachindex(headers)], "  "))
    for row in rows
        println(join([rpad(row[col], widths[col]) for col in eachindex(headers)], "  "))
    end
    return nothing
end

@inline function _wide_value(row, name::String; default::Float64=0.0)::Float64
    column = Symbol(name)
    column in propertynames(row) || return default
    value = row[column]
    (ismissing(value) || !(value isa Number) || !isfinite(value)) && return default
    return Float64(value)
end

function _wide_values(data::DataFrame, name::String)::Vector{Float64}
    name in names(data) || return Float64[]
    return Float64[
        value for value in skipmissing(data[!, name])
        if value isa Number && isfinite(value)
    ]
end

@inline function _safe_ratio(numerator::Real, denominator::Real)::Float64
    denominator > 0 || return NaN
    return 100.0 * Float64(numerator) / Float64(denominator)
end

function _event_statistics(data::DataFrame, metric::String)
    values = _wide_values(data, metric)
    isempty(values) && return nothing
    return (
        total=sum(values),
        affected=count(>(0.0), values),
        incidence_pct=100.0 * count(>(0.0), values) / length(values),
        median=median(values),
        p95=quantile(values, 0.95),
        maximum=maximum(values)
    )
end

function _pathology_event_summary(wide::DataFrame, cases::Vector{String})::DataFrame
    specifications = [
        ("M2T wrong committed", "meas_assoc.committed_wrong"),
        ("T2T wrong committed", "track_assoc.tt_committed_wrong"),
        ("Mixed IOD groups", "cross_m2m.iod_groups_mixed_target"),
        ("Mixed IOD initialized", "cross_m2m.iod_initialized_mixed_target"),
        ("Wrong consensus groups", "track_assoc.consensus_group_mixed_target"),
        ("Identity switches", "track_lifecycle.id_switch_total"),
        ("Excess track fragments", "track_lifecycle.fragment_excess_total"),
        ("Pathological tracks", "track_lifecycle.bad_track_count")
    ]
    rows = NamedTuple[]
    for case in cases
        case_data = wide[wide.case .== case, :]
        for (event, metric) in specifications
            stats = _event_statistics(case_data, metric)
            stats === nothing && continue
            push!(rows, (
                case=case,
                method=_case_label(case),
                event=event,
                metric=metric,
                runs=nrow(case_data),
                total=stats.total,
                affected_runs=stats.affected,
                affected_runs_pct=stats.incidence_pct,
                median_per_run=stats.median,
                p95_per_run=stats.p95,
                maximum_per_run=stats.maximum
            ))
        end
    end
    return DataFrame(rows)
end

function _print_pathology_event_summary(events::DataFrame, cases::Vector{String})
    rows = Vector{String}[]
    for case in cases
        case_events = events[events.case .== case, :]
        for row in eachrow(case_events)
            push!(rows, [
                String(row.method),
                String(row.event),
                _fmt(Float64(row.total)),
                "$(_fmt(Float64(row.affected_runs))) / $(row.runs) ($(_fmt(Float64(row.affected_runs_pct); digits=3))%)",
                _fmt(Float64(row.median_per_run)),
                _fmt(Float64(row.p95_per_run)),
                _fmt(Float64(row.maximum_per_run))
            ])
        end
    end
    _print_table(
        "RARE AND PATHOLOGICAL EVENTS -- totals, incidence, and tails",
        ["Method", "Event", "Total", "Runs with event", "Median", "P95", "Max"],
        rows
    )
    return nothing
end

function _print_iod_validation_diagnostics(wide::DataFrame, cases::Vector{String})
    rows = Vector{String}[]
    for case in cases
        data = wide[wide.case .== case, :]
        enabled_values = _wide_values(data, "cross_m2m.iod_one_step_validation_enabled")
        validation_enabled = !isempty(enabled_values) && maximum(enabled_values) > 0.5

        mixed = sum(_wide_values(data, "cross_m2m.iod_groups_mixed_target"))
        cov_rejected = sum(_wide_values(
            data, "cross_m2m.iod_position_cov_gate_rejected_mixed_target"
        ))
        mixed_rejected = sum(_wide_values(
            data, "cross_m2m.iod_validation_rejected_mixed_target"
        ))
        mixed_confirmed = sum(_wide_values(
            data, "cross_m2m.iod_validation_confirmed_mixed_target"
        ))
        true_rejected = sum(_wide_values(
            data, "cross_m2m.iod_validation_rejected_same_target"
        ))
        true_confirmed = sum(_wide_values(
            data, "cross_m2m.iod_validation_confirmed_same_target"
        ))
        no_measure = sum(_wide_values(data, "cross_m2m.iod_validation_no_measure"))
        pending = sum(_wide_values(data, "cross_m2m.iod_validation_pending_end"))

        push!(rows, [
            _case_label(case),
            _fmt(mixed),
            _fmt(cov_rejected),
            _fmt(_safe_ratio(cov_rejected, mixed); digits=4),
            validation_enabled ? _fmt(mixed_rejected) : "n/a",
            validation_enabled ? _fmt(mixed_confirmed) : "n/a",
            validation_enabled ? _fmt(
                _safe_ratio(mixed_rejected, mixed_rejected + mixed_confirmed); digits=4
            ) : "n/a",
            validation_enabled ? _fmt(
                _safe_ratio(true_rejected, true_rejected + true_confirmed); digits=4
            ) : "n/a",
            validation_enabled ? _fmt(no_measure) : "n/a",
            validation_enabled ? _fmt(pending) : "n/a"
        ])
    end
    _print_table(
        "IOD VALIDATION DIAGNOSTICS -- pooled counts across realizations",
        ["Method", "Mixed groups", "Mixed cov rej.", "Mixed cov rej. [%]",
         "Mixed val. rej.", "Mixed promoted", "Mixed val. rej. [%]",
         "True val. rej. [%]", "No meas.", "Pending"],
        rows
    )
    return nothing
end

@inline function _finite_stat(values, statistic::Function)::Float64
    finite = Float64[
        value for value in skipmissing(values)
        if value isa Number && isfinite(value)
    ]
    return isempty(finite) ? NaN : statistic(finite)
end

function _iod_event_characterization(iod::DataFrame, case::String)::DataFrame
    data = iod[iod.case .== case, :]
    nrow(data) == 0 && return DataFrame()
    rows = NamedTuple[]
    for group in groupby(data, [:stage, :outcome, :same_target])
        d2_values = Float64[
            value for value in skipmissing(group.validation_d2)
            if value isa Number && isfinite(value)
        ]
        nonnegative_d2 = filter(>=(0.0), d2_values)
        ratios = Float64[]
        for row in eachrow(group)
            error = row.position_error_m
            sigma = row.position_rms_std_m
            if !ismissing(error) && !ismissing(sigma) && isfinite(error) &&
               isfinite(sigma) && sigma > 0
                push!(ratios, error / sigma)
            end
        end
        push!(rows, (
            case=case,
            stage=String(first(group.stage)),
            outcome=String(first(group.outcome)),
            target_class=Bool(first(group.same_target)) ? "same target" : "mixed target",
            count=nrow(group),
            position_error_median_m=_finite_stat(group.position_error_m, median),
            position_error_p95_m=_finite_stat(group.position_error_m, x -> quantile(x, 0.95)),
            formal_position_rms_median_m=_finite_stat(group.position_rms_std_m, median),
            error_to_formal_sigma_median=isempty(ratios) ? NaN : median(ratios),
            validation_d2_median=isempty(nonnegative_d2) ? NaN : median(nonnegative_d2),
            validation_d2_p95=isempty(nonnegative_d2) ? NaN : quantile(nonnegative_d2, 0.95),
            negative_validation_d2=count(<(0.0), d2_values)
        ))
    end
    return DataFrame(rows)
end

function _print_iod_event_characterization(events::DataFrame, case::String)
    nrow(events) == 0 && return nothing
    rows = Vector{String}[]
    for row in eachrow(events)
        push!(rows, [
            String(row.stage),
            String(row.outcome),
            String(row.target_class),
            string(row.count),
            _fmt(Float64(row.position_error_median_m)),
            _fmt(Float64(row.formal_position_rms_median_m)),
            _fmt(Float64(row.error_to_formal_sigma_median)),
            _fmt(Float64(row.validation_d2_median)),
            _fmt(Float64(row.validation_d2_p95)),
            string(row.negative_validation_d2)
        ])
    end
    _print_table(
        "EVENT-LEVEL IOD CHARACTERIZATION -- $(_case_label(case))",
        ["Stage", "Outcome", "Class", "N", "Err med. [m]", "Formal RMS [m]",
         "Err/RMS", "d2 med.", "d2 P95", "d2 < 0"],
        rows
    )
    negative_count = sum(events.negative_validation_d2)
    if negative_count > 0
        println()
        println(
            "WARNING: $(Int(negative_count)) saved validation statistics are negative. " *
            "A squared Mahalanobis statistic must be nonnegative; inspect the PSD " *
            "projection/pseudoinverse path before treating these events as chi-square gated."
        )
    end
    return nothing
end

function _pathology_reason(row)::String
    reasons = String[]
    _wide_value(row, "meas_assoc.committed_wrong") > 0 && push!(reasons, "M2T")
    _wide_value(row, "track_assoc.tt_committed_wrong") > 0 && push!(reasons, "T2T")
    _wide_value(row, "cross_m2m.iod_initialized_mixed_target") > 0 && push!(reasons, "IOD")
    _wide_value(row, "track_assoc.consensus_group_mixed_target") > 0 && push!(reasons, "consensus")
    _wide_value(row, "track_lifecycle.id_switch_total") > 0 && push!(reasons, "ID-switch")
    _wide_value(row, "track_lifecycle.fragment_excess_total") > 0 && push!(reasons, "fragmentation")
    _wide_value(row, "track_lifecycle.bad_track_count") > 0 && push!(reasons, "bad-track")
    return isempty(reasons) ? "none" : join(reasons, "+")
end

function _pathological_runs(wide::DataFrame)::DataFrame
    rows = NamedTuple[]
    for row in eachrow(wide)
        reasons = _pathology_reason(row)
        reasons == "none" && continue
        success = _wide_value(row, "tracking.success_rate_possible_pct"; default=NaN)
        push!(rows, (
            realization=Int(round(_wide_value(row, "realization"))),
            case=String(row.case),
            method=_case_label(String(row.case)),
            scenario_seed=Int(round(_wide_value(row, "scenario_seed"))),
            sensor_seed=Int(round(_wide_value(row, "sensor_seed"))),
            observer_od_seed=Int(round(_wide_value(row, "observer_od_seed"))),
            pathology_reasons=reasons,
            m2t_wrong=Int(round(_wide_value(row, "meas_assoc.committed_wrong"))),
            t2t_wrong=Int(round(_wide_value(row, "track_assoc.tt_committed_wrong"))),
            mixed_iod_groups=Int(round(_wide_value(row, "cross_m2m.iod_groups_mixed_target"))),
            mixed_iod_rejected=Int(round(_wide_value(
                row, "cross_m2m.iod_validation_rejected_mixed_target"
            ))),
            mixed_iod_initialized=Int(round(_wide_value(
                row, "cross_m2m.iod_initialized_mixed_target"
            ))),
            wrong_consensus_groups=Int(round(_wide_value(
                row, "track_assoc.consensus_group_mixed_target"
            ))),
            identity_switches=Int(round(_wide_value(row, "track_lifecycle.id_switch_total"))),
            excess_fragments=Int(round(_wide_value(
                row, "track_lifecycle.fragment_excess_total"
            ))),
            pathological_tracks=Int(round(_wide_value(row, "track_lifecycle.bad_track_count"))),
            initialization_success_pct=_wide_value(
                row, "track_lifecycle.initialization_success_pct"; default=NaN
            ),
            tracking_success_pct=success,
            tracking_failure_pct=isfinite(success) ? 100.0 - success : NaN,
            converged_mean_error_m=_wide_value(
                row, "tracking.converged_mean_error_m"; default=NaN
            ),
            converged_median_error_m=_wide_value(
                row, "tracking.converged_median_error_m"; default=NaN
            ),
            max_simultaneously_visible=Int(round(_wide_value(
                row, "geometry.max_simultaneously_visible_targets"
            ))),
            minimum_angular_separation_deg=_wide_value(
                row, "geometry.minimum_angular_separation_deg"; default=NaN
            )
        ))
    end
    return DataFrame(rows)
end

function _print_top_pathological_runs(pathological::DataFrame, case::String; top_n::Int=10)
    nrow(pathological) == 0 && return nothing
    data = pathological[pathological.case .== case, :]
    nrow(data) == 0 && return nothing
    sort!(
        data,
        [:wrong_consensus_groups, :identity_switches, :m2t_wrong, :t2t_wrong,
         :mixed_iod_initialized, :pathological_tracks, :excess_fragments,
         :tracking_failure_pct, :converged_mean_error_m];
        rev=fill(true, 9)
    )
    data = first(data, min(top_n, nrow(data)))
    rows = Vector{String}[]
    for row in eachrow(data)
        push!(rows, [
            string(row.realization),
            String(row.pathology_reasons),
            string(row.mixed_iod_groups),
            string(row.mixed_iod_rejected),
            string(row.mixed_iod_initialized),
            string(row.wrong_consensus_groups),
            string(row.identity_switches),
            string(row.excess_fragments),
            string(row.pathological_tracks),
            _fmt(Float64(row.tracking_success_pct)),
            _fmt(Float64(row.converged_median_error_m)),
            _fmt(Float64(row.converged_mean_error_m)),
            string(row.max_simultaneously_visible),
            _fmt(Float64(row.minimum_angular_separation_deg))
        ])
    end
    _print_table(
        "TOP PATHOLOGICAL REALIZATIONS -- $(_case_label(case)), association-first ranking",
        ["Run", "Flags", "Mixed", "Rejected", "Promoted", "Wrong cons.", "ID sw.",
         "Fragments", "Bad tracks", "Win. [%]", "Err med. [m]", "Err mean [m]",
         "Max vis.", "Min sep. [deg]"],
        rows
    )
    return nothing
end

function _pathology_conditioned_outcomes(wide::DataFrame, case::String)::DataFrame
    data = wide[wide.case .== case, :]
    nrow(data) == 0 && return DataFrame()
    fragment_values = _wide_values(data, "track_lifecycle.fragment_excess_total")
    fragment_threshold = isempty(fragment_values) ? Inf : quantile(fragment_values, 0.95)
    conditions = [
        ("Mixed IOD initialized", "cross_m2m.iod_initialized_mixed_target", 0.0),
        ("Identity switch", "track_lifecycle.id_switch_total", 0.0),
        ("Pathological track", "track_lifecycle.bad_track_count", 0.0),
        ("High fragmentation", "track_lifecycle.fragment_excess_total", fragment_threshold - eps())
    ]
    outcomes = [
        ("M2T recall [%]", "meas_assoc.recall_pct"),
        ("T2T recall [%]", "track_assoc.tt_recall_pct"),
        ("Successful windows [%]", "tracking.success_rate_possible_pct"),
        ("Converged median error [m]", "tracking.converged_median_error_m"),
        ("Converged mean error [m]", "tracking.converged_mean_error_m"),
        ("Excess fragments", "track_lifecycle.fragment_excess_total")
    ]
    rows = NamedTuple[]
    for (condition, condition_metric, threshold) in conditions
        condition_values = _wide_values(data, condition_metric)
        length(condition_values) == nrow(data) || continue
        affected = condition_values .> threshold
        any(affected) && any(.!affected) || continue
        for (outcome, outcome_metric) in outcomes
            values = _wide_values(data, outcome_metric)
            length(values) == nrow(data) || continue
            unaffected_median = median(values[.!affected])
            affected_median = median(values[affected])
            push!(rows, (
                case=case,
                condition=condition,
                condition_metric=condition_metric,
                condition_threshold=threshold,
                affected_runs=count(affected),
                unaffected_runs=count(.!affected),
                outcome=outcome,
                outcome_metric=outcome_metric,
                unaffected_median=unaffected_median,
                affected_median=affected_median,
                affected_minus_unaffected=affected_median - unaffected_median
            ))
        end
    end
    return DataFrame(rows)
end

function _print_pathology_conditioned_outcomes(conditioned::DataFrame, case::String)
    nrow(conditioned) == 0 && return nothing
    for condition in unique(conditioned.condition)
        data = conditioned[conditioned.condition .== condition, :]
        rows = Vector{String}[]
        for row in eachrow(data)
            push!(rows, [
                String(row.outcome),
                _fmt(Float64(row.unaffected_median)),
                _fmt(Float64(row.affected_median)),
                _fmt(Float64(row.affected_minus_unaffected))
            ])
        end
        affected = Int(first(data.affected_runs))
        unaffected = Int(first(data.unaffected_runs))
        _print_table(
            "CONDITIONAL PATHOLOGY EFFECT -- $(_case_label(case)), $(condition) " *
            "($affected affected, $unaffected unaffected runs)",
            ["Outcome", "Without event", "With event", "Difference"],
            rows
        )
    end
    return nothing
end

function _rankdata(values::Vector{Float64})::Vector{Float64}
    order = sortperm(values)
    ranks = similar(values)
    first_idx = 1
    while first_idx <= length(order)
        last_idx = first_idx
        while last_idx < length(order) && values[order[last_idx + 1]] == values[order[first_idx]]
            last_idx += 1
        end
        rank = (first_idx + last_idx) / 2
        for idx in first_idx:last_idx
            ranks[order[idx]] = rank
        end
        first_idx = last_idx + 1
    end
    return ranks
end

function _correlations(x::Vector{Float64}, y::Vector{Float64})
    valid = [idx for idx in eachindex(x) if isfinite(x[idx]) && isfinite(y[idx])]
    length(valid) >= 3 || return (pearson=NaN, spearman=NaN)
    xv = x[valid]
    yv = y[valid]
    std(xv) > 0 && std(yv) > 0 || return (pearson=NaN, spearman=NaN)
    return (pearson=cor(xv, yv), spearman=cor(_rankdata(xv), _rankdata(yv)))
end

function _crowding_summary(wide::DataFrame, case::String)::DataFrame
    data = wide[wide.case .== case, :]
    nrow(data) == 0 && return DataFrame()
    specifications = [
        ("Max visible", "geometry.max_simultaneously_visible_targets", false),
        ("Minimum separation", "geometry.minimum_angular_separation_deg", true)
    ]
    outcomes = [
        ("M2T recall [%]", "meas_assoc.recall_pct"),
        ("T2T recall [%]", "track_assoc.tt_recall_pct"),
        ("Identity switches", "track_lifecycle.id_switch_total"),
        ("Excess fragments", "track_lifecycle.fragment_excess_total"),
        ("Successful windows [%]", "tracking.success_rate_possible_pct"),
        ("Converged mean error [m]", "tracking.converged_mean_error_m"),
        ("Mixed IOD initialized", "cross_m2m.iod_initialized_mixed_target")
    ]
    rows = NamedTuple[]
    for (difficulty_label, difficulty_metric, reversed) in specifications
        x = _wide_values(data, difficulty_metric)
        length(x) == nrow(data) || continue
        q25, q75 = quantile(x, [0.25, 0.75])
        low_crowding = reversed ? x .>= q75 : x .<= q25
        high_crowding = reversed ? x .<= q25 : x .>= q75
        for (outcome_label, outcome_metric) in outcomes
            y = _wide_values(data, outcome_metric)
            length(y) == length(x) || continue
            correlations = _correlations(x, y)
            push!(rows, (
                case=case,
                difficulty=difficulty_label,
                difficulty_metric=difficulty_metric,
                q25=q25,
                q75=q75,
                outcome=outcome_label,
                outcome_metric=outcome_metric,
                low_crowding_median=median(y[low_crowding]),
                high_crowding_median=median(y[high_crowding]),
                high_minus_low=median(y[high_crowding]) - median(y[low_crowding]),
                pearson=correlations.pearson,
                spearman=correlations.spearman
            ))
        end
    end
    return DataFrame(rows)
end

function _print_crowding_summary(crowding::DataFrame, case::String)
    nrow(crowding) == 0 && return nothing
    for difficulty in unique(crowding.difficulty)
        data = crowding[crowding.difficulty .== difficulty, :]
        rows = Vector{String}[]
        for row in eachrow(data)
            push!(rows, [
                String(row.outcome),
                _fmt(Float64(row.low_crowding_median)),
                _fmt(Float64(row.high_crowding_median)),
                _fmt(Float64(row.high_minus_low)),
                _fmt(Float64(row.pearson)),
                _fmt(Float64(row.spearman))
            ])
        end
        q25 = Float64(first(data.q25))
        q75 = Float64(first(data.q75))
        _print_table(
            "CROWDING SENSITIVITY -- $(_case_label(case)), $(difficulty), Q25=$(_fmt(q25)), Q75=$(_fmt(q75))",
            ["Outcome", "Low crowd.", "High crowd.", "Difference", "Pearson", "Spearman"],
            rows
        )
    end
    return nothing
end

function _comparison_dataframe(summary::DataFrame, cases::Vector{String}, counts)::DataFrame
    rows = NamedTuple[]
    for case in cases
        push!(rows, (
            case=case,
            method=_case_label(case),
            runs=get(counts, case, 0),
            m2t_precision_pct=_stat(summary, case, "meas_assoc", "commit_accuracy_pct"),
            m2t_recall_pct=_stat(summary, case, "meas_assoc", "recall_pct"),
            t2t_precision_pct=_stat(summary, case, "track_assoc", "tt_accuracy_pct_known_only"),
            t2t_recall_pct=_stat(summary, case, "track_assoc", "tt_recall_pct"),
            initialization_success_pct=_stat(
                summary, case, "track_lifecycle", "initialization_success_pct"
            ),
            correctly_tracked_targets_pct=_stat(
                summary, case, "object_coverage", "successful_unique_over_jointly_detected_pct"
            ),
            successful_windows_pct=_stat(summary, case, "tracking", "success_rate_possible_pct"),
            converged_error_m=_stat(summary, case, "tracking", "converged_median_error_m"),
            initialization_error_m=_stat(
                summary, case, "track_lifecycle", "initialization_position_error_median_m"
            ),
            initialization_latency_s=_stat(
                summary, case, "track_lifecycle", "initialization_latency_median_s"
            ),
            fragmented_windows_pct=_stat(
                summary, case, "track_lifecycle", "fragmented_window_rate_pct"
            ),
            bad_tracks_per_run=_stat(summary, case, "track_lifecycle", "bad_track_count"),
            bad_track_duration_s=_stat(
                summary, case, "track_lifecycle", "bad_track_filter_duration_median_s"
            ),
            bad_track_mean_duration_s=_stat(
                summary, case, "track_lifecycle", "bad_track_filter_duration_mean_s"
            ),
            good_track_duration_s=_stat(
                summary, case, "track_lifecycle", "good_track_filter_duration_median_s"
            ),
            good_track_mean_duration_s=_stat(
                summary, case, "track_lifecycle", "good_track_filter_duration_mean_s"
            ),
            wrong_iod_promoted=_stat(
                summary, case, "cross_m2m", "iod_validation_confirmed_mixed_target"
            ),
            wrong_consensus_groups=_stat(
                summary, case, "track_assoc", "consensus_group_mixed_target"
            )
        ))
    end
    return DataFrame(rows)
end

function _print_comparisons(summary, comparison, cases, counts)
    performance_rows = Vector{String}[]
    for row in eachrow(comparison)
        push!(performance_rows, [
            String(row.method),
            string(row.runs),
            _interval_text(summary, String(row.case), "meas_assoc", "commit_accuracy_pct"),
            _interval_text(summary, String(row.case), "meas_assoc", "recall_pct"),
            _interval_text(summary, String(row.case), "track_assoc", "tt_accuracy_pct_known_only"),
            _interval_text(summary, String(row.case), "track_assoc", "tt_recall_pct"),
            _interval_text(
                summary, String(row.case), "track_lifecycle", "initialization_success_pct"
            ),
            _interval_text(
                summary, String(row.case), "object_coverage",
                "successful_unique_over_jointly_detected_pct"
            ),
            _interval_text(summary, String(row.case), "tracking", "success_rate_possible_pct"),
            _interval_text(summary, String(row.case), "tracking", "converged_median_error_m")
        ])
    end
    _print_table(
        "PRIMARY PERFORMANCE -- median [p05, p95]",
        ["Method", "N", "M2T P [%]", "M2T R [%]", "T2T P [%]", "T2T R [%]",
         "Init. [%]", "Targets [%]", "Windows [%]", "Conv. err [m]"],
        performance_rows
    )

    initialization_rows = Vector{String}[]
    for case in cases
        push!(initialization_rows, [
            _case_label(case),
            _interval_text(summary, case, "track_lifecycle", "initialization_coverage_pct"),
            _interval_text(summary, case, "track_lifecycle", "initialization_success_pct"),
            _interval_text(summary, case, "track_lifecycle", "initialized_unique_targets"),
            _interval_text(summary, case, "track_lifecycle", "correctly_initialized_unique_targets"),
            _interval_text(summary, case, "track_lifecycle", "never_initialized_unique_targets"),
            _interval_text(summary, case, "track_lifecycle", "never_correctly_initialized_unique_targets"),
            _interval_text(summary, case, "track_lifecycle", "initialization_latency_median_s")
        ])
    end
    _print_table(
        "INITIALIZATION PERFORMANCE -- median [p05, p95] per run",
        ["Method", "Coverage [%]", "Correct [%]", "Initialized", "Correct init.",
         "Never init.", "Never correct", "Latency [s]"],
        initialization_rows
    )

    iod_rows = Vector{String}[]
    for case in cases
        push!(iod_rows, [
            _case_label(case),
            _interval_text(summary, case, "cross_m2m", "iod_groups_mixed_target"),
            _interval_text(summary, case, "cross_m2m", "iod_position_cov_gate_rejected_mixed_target"),
            _interval_text(summary, case, "cross_m2m", "iod_validation_rejected_mixed_target"),
            _interval_text(summary, case, "cross_m2m", "iod_validation_confirmed_mixed_target"),
            _interval_text(summary, case, "cross_m2m", "iod_validation_rejected_same_target"),
            _interval_text(summary, case, "cross_m2m", "iod_promoted_mixed_target_error_median_m"),
            _interval_text(summary, case, "track_assoc", "consensus_group_mixed_target")
        ])
    end
    _print_table(
        "IOD ASSOCIATION OUTCOMES -- median [p05, p95] per run",
        ["Method", "Mixed groups", "Mixed cov-rej.", "Mixed val-rej.",
         "Mixed promoted", "True val-rej.", "Wrong-IOD err [m]", "Wrong cons."],
        iod_rows
    )

    lifecycle_rows = Vector{String}[]
    for case in cases
        push!(lifecycle_rows, [
            _case_label(case),
            _interval_text(summary, case, "track_lifecycle", "initialization_latency_median_s"),
            _interval_text(summary, case, "track_lifecycle", "fragmented_window_rate_pct"),
            _interval_text(summary, case, "track_lifecycle", "good_track_count"),
            _interval_text(summary, case, "track_lifecycle", "good_track_filter_duration_mean_s"),
            _interval_text(summary, case, "track_lifecycle", "good_track_filter_duration_median_s"),
            _interval_text(summary, case, "track_lifecycle", "bad_track_count"),
            _interval_text(summary, case, "track_lifecycle", "bad_track_filter_duration_mean_s"),
            _interval_text(summary, case, "track_lifecycle", "bad_track_filter_duration_median_s"),
            _interval_text(summary, case, "track_lifecycle", "id_switch_total")
        ])
    end
    _print_table(
        "TRACK LIFECYCLE -- median [p05, p95] per run",
        ["Method", "Init lat. [s]", "Frag. windows [%]", "Good tracks",
         "Good mean [s]", "Good med. [s]", "Bad tracks", "Bad mean [s]",
         "Bad med. [s]", "ID switches"],
        lifecycle_rows
    )
    return nothing
end

function main()
    isfile(SUMMARY_PATH) || error("Monte Carlo summary not found: $(SUMMARY_PATH)")
    summary = CSV.read(SUMMARY_PATH, DataFrame)
    status = isfile(STATUS_PATH) ? CSV.read(STATUS_PATH, DataFrame) : DataFrame()
    cases = _available_cases(summary)
    isempty(cases) && error("No completed Monte Carlo methods found in $(SUMMARY_PATH)")
    counts = _completed_counts(status)

    comparison = _comparison_dataframe(summary, cases, counts)
    mkpath(OUTPUT_ROOT)
    comparison_path = joinpath(OUTPUT_ROOT, "mc_method_comparison_table.csv")
    CSV.write(comparison_path, comparison)

    println("Monte Carlo method comparison")
    println("  input: $(OUTPUT_ROOT)")
    println("  methods: $([_case_label(case) for case in cases])")
    println("  values: median across realizations; brackets are p05--p95")
    _print_comparisons(summary, comparison, cases, counts)
    println()
    println("table: $(comparison_path)")

    if isfile(WIDE_METRICS_PATH)
        wide = CSV.read(WIDE_METRICS_PATH, DataFrame)
        events = _pathology_event_summary(wide, cases)
        pathological = _pathological_runs(wide)
        conditioned = _pathology_conditioned_outcomes(wide, "proposed")
        crowding = _crowding_summary(wide, "proposed")

        CSV.write(PATHOLOGY_SUMMARY_PATH, events)
        CSV.write(PATHOLOGICAL_RUNS_PATH, pathological)
        nrow(conditioned) > 0 && CSV.write(PATHOLOGY_OUTCOMES_PATH, conditioned)
        nrow(crowding) > 0 && CSV.write(CROWDING_SUMMARY_PATH, crowding)

        _print_pathology_event_summary(events, cases)
        _print_iod_validation_diagnostics(wide, cases)
        if isfile(IOD_DIAGNOSTICS_PATH)
            iod_events = _iod_event_characterization(
                CSV.read(IOD_DIAGNOSTICS_PATH, DataFrame), "proposed"
            )
            nrow(iod_events) > 0 && CSV.write(IOD_EVENT_SUMMARY_PATH, iod_events)
            _print_iod_event_characterization(iod_events, "proposed")
            nrow(iod_events) > 0 && println("table: $(IOD_EVENT_SUMMARY_PATH)")
        end
        top_n = tryparse(Int, get(ENV, "SPACEAGORA_MC_PATHOLOGY_TOP", "10"))
        _print_top_pathological_runs(
            pathological,
            "proposed";
            top_n=top_n === nothing ? 10 : max(top_n, 1)
        )
        _print_pathology_conditioned_outcomes(conditioned, "proposed")
        _print_crowding_summary(crowding, "proposed")

        println()
        println("table: $(PATHOLOGY_SUMMARY_PATH)")
        println("table: $(PATHOLOGICAL_RUNS_PATH)")
        nrow(conditioned) > 0 && println("table: $(PATHOLOGY_OUTCOMES_PATH)")
        nrow(crowding) > 0 && println("table: $(CROWDING_SUMMARY_PATH)")
    else
        println()
        println("Pathology diagnostics skipped: wide metrics not found at $(WIDE_METRICS_PATH)")
    end

    return nothing
end

main()
