using Printf

const COMPACT_NAVIGATION_METRICS = (
    ("meas_assoc", "commit_accuracy_pct", "M2T precision"),
    ("meas_assoc", "recall_pct", "M2T recall"),
    ("track_assoc", "tt_accuracy_pct_known_only", "T2T precision"),
    ("track_assoc", "tt_recall_pct", "T2T recall"),
    ("track_lifecycle", "initialization_success_pct", "Initialization success"),
    ("tracking", "success_rate_possible_pct", "Tracking success"),
    (
        "object_coverage",
        "successful_unique_over_jointly_detected_pct",
        "Correctly tracked / trackable targets"
    ),
)

function _compact_metric(
    table::DataFrame,
    section::String,
    metric::String;
    value_column::Symbol=:value
)::Float64
    value_column in propertynames(table) || return NaN
    rows = table[
        (String.(table.section) .== section) .&
        (String.(table.metric) .== metric),
        value_column
    ]
    length(rows) == 1 || return NaN
    value = rows[1]
    return ismissing(value) ? NaN : Float64(value)
end

function _compact_number(value::Float64; digits::Int=5)::String
    isfinite(value) || return "n/a"
    abs(value - round(value)) <= 1e-10 && return string(Int(round(value)))
    return @sprintf("%.*g", digits, value)
end

function _print_compact_navigation_metrics(
    table::DataFrame;
    value_column::Symbol=:value
)::Nothing
    for (section, metric, label) in COMPACT_NAVIGATION_METRICS
        value = _compact_metric(
            table,
            section,
            metric;
            value_column=value_column
        )
        value_text = isfinite(value) ? "$(_compact_number(value))%" : "n/a"
        println("    $(label): $(value_text)")
    end

    converged_error = _compact_metric(
        table,
        "tracking",
        "converged_median_error_m";
        value_column=value_column
    )
    fragments = _compact_metric(
        table,
        "track_lifecycle",
        "fragment_excess_total";
        value_column=value_column
    )
    id_switches = _compact_metric(
        table,
        "track_lifecycle",
        "id_switch_total";
        value_column=value_column
    )
    println("    Converged median error: $(_compact_number(converged_error)) m")
    println("    Track fragments / identity switches: " *
            "$(_compact_number(fragments)) / $(_compact_number(id_switches))")
    return nothing
end
