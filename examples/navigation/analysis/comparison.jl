using CSV
using DataFrames
using Dates
using Printf
using Statistics

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using Plots
using Plots.PlotMeasures: mm

if !isdefined(@__MODULE__, :NavigationPaths)
    include(joinpath(@__DIR__, "..", "paths.jl"))
end
using .NavigationPaths: env_override_path, resolve_output_path
using .NavigationPaths: navigation_output_path, transient_navigation_path

const OUTPUT_ROOT = env_override_path(
    "SPACEAGORA_COMPARISON_OUTPUT",
    transient_navigation_path("runs", "comparison")
)
const PLOT_OUTPUT_ROOT = env_override_path(
    "SPACEAGORA_PLOT_OUTPUT",
    navigation_output_path("figures")
)
const DEFAULT_BASELINE_THETA_SWEEP = "0.01,0.05"
const EARTH_MU_M3S2 = 3.986004418e14
const COMMUNICATION_RANGE_KM = 300.0
const PAPER_FONT_FAMILY = "Times"
const PAPER_GUIDE_FONT = 22
const PAPER_TICK_FONT = 16
const PAPER_LEGEND_FONT = 16
const REGENERATE_PLOTS = lowercase(strip(get(ENV, "SPACEAGORA_REGENERATE_PLOTS", "true"))) in ("1", "true", "yes", "y")
const PRINT_PROGRESS = !(lowercase(strip(get(ENV, "SPACEAGORA_PRINT_PROGRESS", "true"))) in ("0", "false", "no", "n"))
const DEFAULT_NAV_CASES = (
    :proposed,
    :centralized_oracle,
    :independent_local_da,
    :distributed_oracle_da,
    :baseline_da
)

function _progress(message::String)::Nothing
    PRINT_PROGRESS || return nothing
    mkpath(OUTPUT_ROOT)
    timestamp = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS")
    line = "[$(timestamp) UTC] $(message)"
    println(stderr, line)
    flush(stderr)
    open(joinpath(OUTPUT_ROOT, "print_navigation_comparison_metrics_progress.log"), "a") do io
        println(io, line)
    end
    return nothing
end

function _cached_plot_path(basename::String)::Union{Nothing, String}
    REGENERATE_PLOTS && return nothing
    png_path = joinpath(PLOT_OUTPUT_ROOT, "$(basename).png")
    return isfile(png_path) ? png_path : nothing
end

function _baseline_theta_values()::Vector{Float64}
    raw_values = strip(get(ENV, "SPACEAGORA_BASELINE_THETA_SWEEP", DEFAULT_BASELINE_THETA_SWEEP))
    isempty(raw_values) && return Float64[]
    return [parse(Float64, strip(value)) for value in split(raw_values, ",") if !isempty(strip(value))]
end

function _baseline_label(theta_gate_rad::Float64)::Symbol
    token = replace(replace(string(theta_gate_rad), "." => "p"), "-" => "m")
    return Symbol("baseline_da_theta_$(token)")
end

function _selected_cases()::Vector{Symbol}
    raw_cases = strip(get(ENV, "SPACEAGORA_CASES", ""))
    nav_cases = isempty(raw_cases) ? collect(DEFAULT_NAV_CASES) :
                [Symbol(strip(case)) for case in split(raw_cases, ",") if !isempty(strip(case))]

    cases = Symbol[]
    seen = Set{Symbol}()
    for nav_case in nav_cases
        labels = if nav_case === :baseline_da
            theta_values = _baseline_theta_values()
            isempty(theta_values) ? Symbol[:baseline_da] : [_baseline_label(theta) for theta in theta_values]
        else
            Symbol[nav_case]
        end

        for label in labels
            label in seen && continue
            push!(cases, label)
            push!(seen, label)
        end
    end
    return cases
end

function _case_dir(case::Symbol)::String
    return joinpath(OUTPUT_ROOT, String(case))
end

function _exit_code_hint(exit_code::Int)::String
    exit_code == 137 && return "likely SIGKILL/OOM"
    exit_code == 143 && return "likely SIGTERM/time limit"
    exit_code == 130 && return "interrupted"
    return ""
end

function _failure_reason_from_log(log_path::String, exit_code::Union{Nothing, Int}=nothing)::String
    isfile(log_path) || return ""
    lines = readlines(log_path)
    isempty(lines) && return ""
    patterns = (
        "ERROR:",
        "LoadError",
        "BoundsError",
        "DimensionMismatch",
        "OutOfMemory",
        "UndefVarError",
        "MethodError",
        "ArgumentError",
        "AssertionError",
        "Solve failed",
        "retcode",
        "Killed",
        "Terminated",
        "signal"
    )
    for line in Iterators.reverse(lines)
        if any(pattern -> occursin(pattern, line), patterns)
            return strip(line)
        end
    end
    tail_start = max(1, length(lines) - 4)
    tail = strip(join(lines[tail_start:end], " | "))
    exit_code === nothing && return "no explicit Julia error found. Last log lines: $(tail)"
    hint = _exit_code_hint(exit_code)
    prefix = isempty(hint) ?
             "process exited with code $(exit_code), but log contains no explicit Julia error" :
             "process exited with code $(exit_code) ($(hint)), but log contains no explicit Julia error"
    return "$(prefix). Last log lines: $(tail)"
end

function _read_status_by_case()
    path = joinpath(OUTPUT_ROOT, "run_status.csv")
    isfile(path) || return Dict{String, Tuple{String, String, Float64}}()
    df = CSV.read(path, DataFrame)
    (:case in propertynames(df) && :status in propertynames(df)) || return Dict{String, Tuple{String, String, Float64}}()

    status_by_case = Dict{String, Tuple{String, String, Float64}}()
    for row in eachrow(df)
        status = String(row.status)
        reason = :failure_reason in propertynames(df) && !ismissing(row.failure_reason) ? String(row.failure_reason) : ""
        theta_gate = :baseline_theta_gate_rad in propertynames(df) && !ismissing(row.baseline_theta_gate_rad) ?
                     Float64(row.baseline_theta_gate_rad) : NaN
        if isempty(reason) && status != "ok" && :log_path in propertynames(df) && !ismissing(row.log_path)
            exit_code = :exit_code in propertynames(df) && !ismissing(row.exit_code) ? Int(row.exit_code) : nothing
            log_path = resolve_output_path(String(row.log_path), OUTPUT_ROOT)
            reason = _failure_reason_from_log(log_path, exit_code)
        end
        status_by_case[String(row.case)] = (status, reason, theta_gate)
    end
    return status_by_case
end

function _read_metrics(case::Symbol)::Dict{String, Float64}
    path = joinpath(_case_dir(case), "association_quality_table.csv")
    isfile(path) || return Dict{String, Float64}()

    df = CSV.read(path, DataFrame)
    metrics = Dict{String, Float64}()
    for row in eachrow(df)
        key = string(row.section, ".", row.metric)
        metrics[key] = Float64(row.value)
    end
    return metrics
end

function _read_tracking_observer_summary(cases::Vector{Symbol})::DataFrame
    tables = DataFrame[]
    for case in cases
        path = joinpath(_case_dir(case), "tracking_observer_summary.csv")
        isfile(path) || continue
        df = CSV.read(path, DataFrame)
        if !(:case in propertynames(df))
            insertcols!(df, 1, :case => fill(String(case), nrow(df)))
        end
        for col in (:detected_unique_targets, :possible_unique_targets, :tracked_unique_targets, :successful_unique_targets)
            if !(col in propertynames(df))
                df[!, col] = fill(NaN, nrow(df))
            end
        end
        push!(tables, df)
    end
    if isempty(tables)
        return DataFrame(
            case=String[],
            nav_case=String[],
            observer=Int[],
            possible_windows=Float64[],
            tracked_windows=Float64[],
            tracking_coverage_pct=Float64[],
            successful_windows_under_1km=Float64[],
            success_rate_tracked_pct=Float64[],
            success_rate_possible_pct=Float64[],
            detected_unique_targets=Float64[],
            possible_unique_targets=Float64[],
            tracked_unique_targets=Float64[],
            successful_unique_targets=Float64[],
            mean_error_tracked_windows_m=Float64[],
            mean_error_successful_windows_m=Float64[]
        )
    end
    return vcat(tables...; cols=:union)
end

function _first_run_log(cases::Vector{Symbol})::Union{Nothing, String}
    for case in cases
        path = joinpath(_case_dir(case), "run.log")
        isfile(path) && return path
    end
    return nothing
end

function _observer_count_for_plot(cases::Vector{Symbol})::Int
    log_path = _first_run_log(cases)
    if log_path !== nothing
        for line in eachline(log_path)
            m = match(r"Observer idxs:\s*\[(.*)\]", line)
            if m !== nothing
                raw_idxs = split(m.captures[1], ",")
                return count(token -> !isempty(strip(token)), raw_idxs)
            end
        end
    end
    return 16
end

function _observer_defs_for_plot(observer_count::Int)
    e = 1e-4
    raan_deg = 10.0
    aop_deg = 14.0
    observer_a_m = 6_963.0e3

    if observer_count <= 3
        observer_i_deg = 85.0
        observer_mean_anomalies_deg = (288.0, 290.0, 292.0)
        observer_defs = Tuple((
            a_m=observer_a_m,
            e=e,
            i_deg=observer_i_deg,
            raan_deg=raan_deg,
            aop_deg=aop_deg,
            M_deg=M_deg
        ) for M_deg in observer_mean_anomalies_deg)
        return observer_defs, observer_i_deg, observer_mean_anomalies_deg
    end

    observer_i_deg = 70.0
    observer_mean_anomaly_center_deg = 290.0
    observer_raan_offsets_deg = (-1.6, -0.5, 0.5, 1.6)
    observer_mean_anomaly_offsets_deg = (-4.0, -1.3, 1.3, 4.0)
    observer_defs = Tuple([
        (
            a_m=observer_a_m,
            e=e,
            i_deg=observer_i_deg,
            raan_deg=mod(raan_deg + raan_offset_deg, 360.0),
            aop_deg=aop_deg,
            M_deg=observer_mean_anomaly_center_deg + m_offset_deg
        )
        for raan_offset_deg in observer_raan_offsets_deg
        for m_offset_deg in observer_mean_anomaly_offsets_deg
    ])
    observer_mean_anomalies_deg = Tuple(obs.M_deg for obs in observer_defs)
    return observer_defs, observer_i_deg, observer_mean_anomalies_deg
end

function _kepler_xyz_km(orbit_def, t_s::Float64=0.0)
    a = Float64(orbit_def.a_m)
    e = hasproperty(orbit_def, :e) ? Float64(orbit_def.e) : 1e-4
    inc = deg2rad(Float64(orbit_def.i_deg))
    raan = hasproperty(orbit_def, :raan_deg) ? deg2rad(Float64(orbit_def.raan_deg)) : deg2rad(10.0)
    aop = hasproperty(orbit_def, :aop_deg) ? deg2rad(Float64(orbit_def.aop_deg)) : deg2rad(14.0)
    mean0 = deg2rad(Float64(orbit_def.M_deg))
    n = sqrt(EARTH_MU_M3S2 / a^3)

    M = mod(mean0 + n * t_s, 2π)
    E = M
    for _ in 1:10
        E -= (E - e * sin(E) - M) / (1.0 - e * cos(E))
    end

    x_pf = a * (cos(E) - e)
    y_pf = a * sqrt(1.0 - e^2) * sin(E)
    cw = cos(aop)
    sw = sin(aop)
    cO = cos(raan)
    sO = sin(raan)
    ci = cos(inc)
    si = sin(inc)

    x_arg = cw * x_pf - sw * y_pf
    y_arg = sw * x_pf + cw * y_pf
    return (
        (cO * x_arg - sO * ci * y_arg) / 1000.0,
        (sO * x_arg + cO * ci * y_arg) / 1000.0,
        (si * y_arg) / 1000.0
    )
end

function _kepler_xyz_series_km(orbit_def, times_s)
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]
    sizehint!(xs, length(times_s))
    sizehint!(ys, length(times_s))
    sizehint!(zs, length(times_s))
    for t in times_s
        x, y, z = _kepler_xyz_km(orbit_def, Float64(t))
        push!(xs, x)
        push!(ys, y)
        push!(zs, z)
    end
    return xs, ys, zs
end

@inline function _metric(metrics::Dict{String, Float64}, key::String)::Float64
    return get(metrics, key, NaN)
end

@inline function _pct(num::Float64, den::Float64)::Float64
    return (isfinite(num) && isfinite(den) && den > 0.0) ? 100.0 * num / den : NaN
end

function _build_summary(cases::Vector{Symbol})::DataFrame
    status_by_case = _read_status_by_case()

    rows = NamedTuple[]
    for case in cases
        metrics = _read_metrics(case)
        status_info = get(status_by_case, String(case), (isempty(metrics) ? "missing" : "available", "", NaN))
        status = status_info[1]
        failure_reason = status_info[2]
        baseline_theta_gate_rad = _metric(metrics, "run_config.baseline_theta_gate_rad")
        if !isfinite(baseline_theta_gate_rad)
            baseline_theta_gate_rad = status_info[3]
        end
        if status == "ok" && isempty(metrics)
            status = "missing_metrics"
            isempty(failure_reason) && (failure_reason = "association_quality_table.csv missing")
        end

        m2t_total = _metric(metrics, "meas_assoc.committed_total")
        m2t_correct = _metric(metrics, "meas_assoc.committed_correct")
        m2t_wrong = _metric(metrics, "meas_assoc.committed_wrong")
        m2t_acc = _metric(metrics, "meas_assoc.commit_accuracy_pct")

        h1_created = _metric(metrics, "local_hypothesis.H1_created")
        h1_to_h2 = _metric(metrics, "local_hypothesis.H1_to_H2_created")
        h1_to_h2_same = _metric(metrics, "local_hypothesis.H1_to_H2_same_target")
        h1_to_h2_mixed = _metric(metrics, "local_hypothesis.H1_to_H2_mixed_target")
        h2_to_h3 = _metric(metrics, "local_hypothesis.H2_to_H3_attempted")
        h3_pass = _metric(metrics, "local_hypothesis.H3_los_rate_pass")
        h3_fail = _metric(metrics, "local_hypothesis.H3_los_rate_fail")
        promoted = _metric(metrics, "local_hypothesis.promoted")
        promoted_same = _metric(metrics, "local_hypothesis.promoted_same_target")
        promoted_mixed = _metric(metrics, "local_hypothesis.promoted_mixed_target")
        tracks_created_with_nonreal_measurements = _metric(metrics, "local_hypothesis.tracks_created_with_nonreal_measurements")

        xm2m_candidates = _metric(metrics, "cross_m2m.candidate_pairs")
        xm2m_gate_pass = _metric(metrics, "cross_m2m.gate_pass")
        xm2m_gate_pass_same = _metric(metrics, "cross_m2m.gate_pass_same_target")
        xm2m_gate_pass_mixed = _metric(metrics, "cross_m2m.gate_pass_mixed_target")
        xm2m_selected = _metric(metrics, "cross_m2m.selected_pairs")
        xm2m_selected_same = _metric(metrics, "cross_m2m.selected_pairs_same_target")
        xm2m_selected_mixed = _metric(metrics, "cross_m2m.selected_pairs_mixed_target")
        xm2m_iod_groups = _metric(metrics, "cross_m2m.iod_groups")
        xm2m_iod_groups_same = _metric(metrics, "cross_m2m.iod_groups_same_target")
        xm2m_iod_groups_mixed = _metric(metrics, "cross_m2m.iod_groups_mixed_target")
        iod_pos_cov_gate_threshold = _metric(metrics, "cross_m2m.iod_position_cov_gate_threshold_rms_m")
        iod_pos_cov_gate_evaluated = _metric(metrics, "cross_m2m.iod_position_cov_gate_evaluated")
        iod_pos_cov_gate_rejected = _metric(metrics, "cross_m2m.iod_position_cov_gate_rejected")
        iod_pos_cov_gate_rejected_same = _metric(metrics, "cross_m2m.iod_position_cov_gate_rejected_same_target")
        iod_pos_cov_gate_rejected_mixed = _metric(metrics, "cross_m2m.iod_position_cov_gate_rejected_mixed_target")
        iod_pos_cov_gate_rejected_pct = _metric(metrics, "cross_m2m.iod_position_cov_gate_rejected_pct")
        iod_validation_enabled = _metric(metrics, "cross_m2m.iod_one_step_validation_enabled")
        iod_validation_threshold = _metric(metrics, "cross_m2m.iod_validation_threshold_d2")
        iod_validation_attempted = _metric(metrics, "cross_m2m.iod_validation_attempted")
        iod_validation_confirmed = _metric(metrics, "cross_m2m.iod_validation_confirmed")
        iod_validation_rejected = _metric(metrics, "cross_m2m.iod_validation_rejected")
        iod_validation_confirmed_pct = _metric(metrics, "cross_m2m.iod_validation_confirmed_pct")
        iod_validation_confirmed_same = _metric(metrics, "cross_m2m.iod_validation_confirmed_same_target")
        iod_validation_confirmed_mixed = _metric(metrics, "cross_m2m.iod_validation_confirmed_mixed_target")
        iod_validation_rejected_same = _metric(metrics, "cross_m2m.iod_validation_rejected_same_target")
        iod_validation_rejected_mixed = _metric(metrics, "cross_m2m.iod_validation_rejected_mixed_target")
        iod_validation_no_measure = _metric(metrics, "cross_m2m.iod_validation_no_measure")
        iod_validation_pending_end = _metric(metrics, "cross_m2m.iod_validation_pending_end")
        xm2m_iod_initialized = _metric(metrics, "cross_m2m.iod_initialized")
        xm2m_iod_initialized_same = _metric(metrics, "cross_m2m.iod_initialized_same_target")
        xm2m_iod_initialized_mixed = _metric(metrics, "cross_m2m.iod_initialized_mixed_target")

        tt_attempt = _metric(metrics, "track_assoc.tt_attempt_total")
        tt_commit = _metric(metrics, "track_assoc.tt_committed_total")
        tt_correct = _metric(metrics, "track_assoc.tt_committed_correct")
        tt_wrong = _metric(metrics, "track_assoc.tt_committed_wrong")
        tt_acc = _metric(metrics, "track_assoc.tt_accuracy_pct_known_only")

        group_total = _metric(metrics, "track_assoc.consensus_group_total")
        group_good = _metric(metrics, "track_assoc.consensus_group_same_target")
        group_mixed = _metric(metrics, "track_assoc.consensus_group_mixed_target")
        group_good_pct = _metric(metrics, "track_assoc.consensus_group_same_target_pct_known_only")

        tracking_possible = _metric(metrics, "tracking.possible_windows")
        tracking_tracked = _metric(metrics, "tracking.tracked_windows")
        tracking_success = _metric(metrics, "tracking.successful_windows_under_1km")
        tracking_sample_weighted_mean_error = _metric(metrics, "tracking.sample_weighted_mean_error_m")
        tracking_sample_rmse_error = _metric(metrics, "tracking.sample_rmse_error_m")
        mean_estimate_duration = _metric(metrics, "tracking.mean_estimate_duration_s")
        max_estimate_duration = _metric(metrics, "tracking.max_estimate_duration_s")
        mean_successful_estimate_duration = _metric(metrics, "tracking.mean_successful_estimate_duration_s")
        window_path = joinpath(_case_dir(case), "tracking_window_table.csv")
        if isfile(window_path)
            window_df = CSV.read(window_path, DataFrame)
            if :estimate_duration_s in propertynames(window_df)
                durations = [
                    Float64(row.estimate_duration_s) for row in eachrow(window_df)
                    if row.tracked == true && isfinite(Float64(row.estimate_duration_s))
                ]
                if !isfinite(mean_estimate_duration) && !isempty(durations)
                    mean_estimate_duration = sum(durations) / length(durations)
                end
                if !isfinite(max_estimate_duration) && !isempty(durations)
                    max_estimate_duration = maximum(durations)
                end
                if !isfinite(mean_successful_estimate_duration) && :success_under_1km in propertynames(window_df)
                    good_durations = [
                        Float64(row.estimate_duration_s) for row in eachrow(window_df)
                        if row.success_under_1km == true && isfinite(Float64(row.estimate_duration_s))
                    ]
                    if !isempty(good_durations)
                        mean_successful_estimate_duration = sum(good_durations) / length(good_durations)
                    end
                end
            end
            if !isfinite(tracking_sample_weighted_mean_error) &&
                    all(name -> name in propertynames(window_df), [:tracked, :estimate_samples, :mean_error_m])
                sample_count = 0
                weighted_error = 0.0
                for row in eachrow(window_df)
                    row.tracked == true || continue
                    samples = Int(row.estimate_samples)
                    samples > 0 || continue
                    mean_err = Float64(row.mean_error_m)
                    isfinite(mean_err) || continue
                    sample_count += samples
                    weighted_error += mean_err * samples
                end
                if sample_count > 0
                    tracking_sample_weighted_mean_error = weighted_error / sample_count
                end
            end
            if !isfinite(tracking_sample_rmse_error) &&
                    all(name -> name in propertynames(window_df), [:tracked, :estimate_samples, :rmse_error_m])
                sample_count = 0
                weighted_sq_error = 0.0
                for row in eachrow(window_df)
                    row.tracked == true || continue
                    samples = Int(row.estimate_samples)
                    samples > 0 || continue
                    rmse_err = Float64(row.rmse_error_m)
                    isfinite(rmse_err) || continue
                    sample_count += samples
                    weighted_sq_error += rmse_err^2 * samples
                end
                if sample_count > 0
                    tracking_sample_rmse_error = sqrt(weighted_sq_error / sample_count)
                end
            end
        end
        detected_unique_targets = _metric(metrics, "object_coverage.detected_unique_targets")
        jointly_detected_unique_targets = _metric(metrics, "object_coverage.jointly_detected_unique_targets")
        tracked_unique_targets = _metric(metrics, "object_coverage.tracked_unique_targets")
        successful_tracked_unique_targets = _metric(metrics, "object_coverage.successful_tracked_unique_targets")
        sensor_enable_missed = _metric(metrics, "sensor_errors.enable_missed_detections")
        sensor_enable_false_alarm = _metric(metrics, "sensor_errors.enable_false_alarms")
        sensor_enable_bias = _metric(metrics, "sensor_errors.enable_measurement_bias")
        sensor_configured_missed_rate = _metric(metrics, "sensor_errors.configured_missed_detection_rate")
        sensor_configured_false_alarm_rate = _metric(metrics, "sensor_errors.configured_false_alarm_rate")
        sensor_configured_bias_rad = _metric(metrics, "sensor_errors.configured_measurement_bias_rad")
        sensor_visible_opportunities = _metric(metrics, "sensor_errors.visible_opportunities")
        sensor_true_detections = _metric(metrics, "sensor_errors.true_detections")
        sensor_missed_detections = _metric(metrics, "sensor_errors.missed_detections")
        sensor_false_alarms = _metric(metrics, "sensor_errors.false_alarms")
        sensor_biased_measurements = _metric(metrics, "sensor_errors.biased_measurements")
        sensor_detection_rate_pct = _metric(metrics, "sensor_errors.realized_detection_rate_pct")
        sensor_missed_rate_pct = _metric(metrics, "sensor_errors.realized_missed_rate_pct")
        runtime_timed_epochs = _metric(metrics, "runtime.timed_epoch_count")
        runtime_total_mean_ms = _metric(metrics, "runtime.total_epoch_mean_ms")
        runtime_total_median_ms = _metric(metrics, "runtime.total_epoch_median_ms")
        runtime_total_p95_ms = _metric(metrics, "runtime.total_epoch_p95_ms")
        runtime_total_max_ms = _metric(metrics, "runtime.total_epoch_max_ms")
        runtime_local_da_mean_ms = _metric(metrics, "runtime.local_da_mean_ms")
        runtime_cross_da_mean_ms = _metric(metrics, "runtime.cross_da_mean_ms")
        runtime_filter_mean_ms = _metric(metrics, "runtime.filter_mean_ms")
        runtime_fusion_mean_ms = _metric(metrics, "runtime.fusion_mean_ms")
        wrong_association_total = _metric(metrics, "association_health.wrong_association_total")
        id_switch_total = _metric(metrics, "track_lifecycle.id_switch_total")

        push!(
            rows,
            (
                case=String(case),
                status=status,
                failure_reason=failure_reason,
                baseline_theta_gate_rad=baseline_theta_gate_rad,

                m2t_committed=m2t_total,
                m2t_correct=m2t_correct,
                m2t_wrong=m2t_wrong,
                m2t_accuracy_pct=m2t_acc,
                m2t_ambig_committed=_metric(metrics, "meas_assoc.ambig_committed"),
                m2t_ambig_dropped=_metric(metrics, "meas_assoc.ambig_dropped"),

                h1_created=h1_created,
                h1_to_h2=h1_to_h2,
                h1_to_h2_same=h1_to_h2_same,
                h1_to_h2_mixed=h1_to_h2_mixed,
                h1_to_h2_good_pct=_pct(h1_to_h2_same, h1_to_h2),
                h2_to_h3_attempted=h2_to_h3,
                los_rate_pass=h3_pass,
                los_rate_fail=h3_fail,
                los_rate_pass_pct=_pct(h3_pass, h2_to_h3),
                promoted=promoted,
                promoted_same_target=promoted_same,
                promoted_mixed_target=promoted_mixed,
                promoted_good_pct=_pct(promoted_same, promoted),
                promoted_per_h1_pct=_pct(promoted, h1_created),
                tracks_created_with_nonreal_measurements=tracks_created_with_nonreal_measurements,

                xm2m_candidate_pairs=xm2m_candidates,
                xm2m_gate_pass=xm2m_gate_pass,
                xm2m_gate_pass_same_target=xm2m_gate_pass_same,
                xm2m_gate_pass_mixed_target=xm2m_gate_pass_mixed,
                xm2m_gate_pass_good_pct=_pct(xm2m_gate_pass_same, xm2m_gate_pass),
                xm2m_selected_pairs=xm2m_selected,
                xm2m_selected_same_target=xm2m_selected_same,
                xm2m_selected_mixed_target=xm2m_selected_mixed,
                xm2m_selected_good_pct=_pct(xm2m_selected_same, xm2m_selected),
                xm2m_iod_groups=xm2m_iod_groups,
                xm2m_iod_groups_same_target=xm2m_iod_groups_same,
                xm2m_iod_groups_mixed_target=xm2m_iod_groups_mixed,
                xm2m_iod_groups_good_pct=_pct(xm2m_iod_groups_same, xm2m_iod_groups),
                iod_position_cov_gate_threshold_rms_m=iod_pos_cov_gate_threshold,
                iod_position_cov_gate_evaluated=iod_pos_cov_gate_evaluated,
                iod_position_cov_gate_rejected=iod_pos_cov_gate_rejected,
                iod_position_cov_gate_rejected_pct=iod_pos_cov_gate_rejected_pct,
                iod_position_cov_gate_rejected_same_target=iod_pos_cov_gate_rejected_same,
                iod_position_cov_gate_rejected_mixed_target=iod_pos_cov_gate_rejected_mixed,
                iod_one_step_validation_enabled=iod_validation_enabled,
                iod_validation_threshold_d2=iod_validation_threshold,
                iod_validation_attempted=iod_validation_attempted,
                iod_validation_confirmed=iod_validation_confirmed,
                iod_validation_rejected=iod_validation_rejected,
                iod_validation_confirmed_pct=iod_validation_confirmed_pct,
                iod_validation_confirmed_same_target=iod_validation_confirmed_same,
                iod_validation_confirmed_mixed_target=iod_validation_confirmed_mixed,
                iod_validation_rejected_same_target=iod_validation_rejected_same,
                iod_validation_rejected_mixed_target=iod_validation_rejected_mixed,
                iod_validation_no_measure=iod_validation_no_measure,
                iod_validation_pending_end=iod_validation_pending_end,
                xm2m_iod_initialized=xm2m_iod_initialized,
                xm2m_iod_initialized_same_target=xm2m_iod_initialized_same,
                xm2m_iod_initialized_mixed_target=xm2m_iod_initialized_mixed,
                xm2m_iod_initialized_good_pct=_pct(xm2m_iod_initialized_same, xm2m_iod_initialized),

                t2t_attempted=tt_attempt,
                t2t_committed=tt_commit,
                t2t_correct=tt_correct,
                t2t_wrong=tt_wrong,
                t2t_commit_rate_pct=_pct(tt_commit, tt_attempt),
                t2t_accuracy_pct=tt_acc,
                t2t_ratio_fail=_metric(metrics, "track_assoc.tt_skipped_ratio_fail"),
                t2t_mutual_fail=_metric(metrics, "track_assoc.tt_skipped_mutual_fail"),
                t2t_component_conflict_rejected=_metric(metrics, "track_assoc.tt_component_conflict_rejected"),

                consensus_groups=group_total,
                consensus_same_target=group_good,
                consensus_mixed_target=group_mixed,
                consensus_good_pct=group_good_pct,

                tracking_possible_windows=tracking_possible,
                tracking_tracked_windows=tracking_tracked,
                tracking_coverage_pct=_metric(metrics, "tracking.tracking_coverage_pct"),
                tracking_successful_windows_under_1km=tracking_success,
                tracking_success_rate_tracked_pct=_metric(metrics, "tracking.success_rate_tracked_pct"),
                tracking_success_rate_possible_pct=_metric(metrics, "tracking.success_rate_possible_pct"),
                tracking_mean_error_all_m=_metric(metrics, "tracking.mean_error_tracked_windows_m"),
                tracking_mean_error_good_m=_metric(metrics, "tracking.mean_error_successful_windows_m"),
                tracking_sample_weighted_mean_error_m=tracking_sample_weighted_mean_error,
                tracking_sample_rmse_error_m=tracking_sample_rmse_error,
                tracking_success_threshold_m=_metric(metrics, "tracking.success_error_threshold_m"),
                mean_estimate_duration_s=mean_estimate_duration,
                max_estimate_duration_s=max_estimate_duration,
                mean_successful_estimate_duration_s=mean_successful_estimate_duration,

                detected_unique_targets=detected_unique_targets,
                jointly_detected_unique_targets=jointly_detected_unique_targets,
                tracked_unique_targets=tracked_unique_targets,
                successful_tracked_unique_targets=successful_tracked_unique_targets,
                tracked_unique_over_detected_pct=_metric(metrics, "object_coverage.tracked_unique_over_detected_pct"),
                successful_unique_over_detected_pct=_metric(metrics, "object_coverage.successful_unique_over_detected_pct"),
                tracked_unique_over_jointly_detected_pct=_metric(metrics, "object_coverage.tracked_unique_over_jointly_detected_pct"),
                successful_unique_over_jointly_detected_pct=_metric(metrics, "object_coverage.successful_unique_over_jointly_detected_pct"),

                sensor_enable_missed_detections=sensor_enable_missed,
                sensor_enable_false_alarms=sensor_enable_false_alarm,
                sensor_enable_measurement_bias=sensor_enable_bias,
                sensor_configured_missed_detection_rate=sensor_configured_missed_rate,
                sensor_configured_false_alarm_rate=sensor_configured_false_alarm_rate,
                sensor_configured_measurement_bias_rad=sensor_configured_bias_rad,
                sensor_visible_opportunities=sensor_visible_opportunities,
                sensor_true_detections=sensor_true_detections,
                sensor_missed_detections=sensor_missed_detections,
                sensor_false_alarms=sensor_false_alarms,
                sensor_biased_measurements=sensor_biased_measurements,
                sensor_realized_detection_rate_pct=sensor_detection_rate_pct,
                sensor_realized_missed_rate_pct=sensor_missed_rate_pct,

                runtime_timed_epochs=runtime_timed_epochs,
                runtime_total_epoch_mean_ms=runtime_total_mean_ms,
                runtime_total_epoch_median_ms=runtime_total_median_ms,
                runtime_total_epoch_p95_ms=runtime_total_p95_ms,
                runtime_total_epoch_max_ms=runtime_total_max_ms,
                runtime_local_da_mean_ms=runtime_local_da_mean_ms,
                runtime_cross_da_mean_ms=runtime_cross_da_mean_ms,
                runtime_filter_mean_ms=runtime_filter_mean_ms,
                runtime_fusion_mean_ms=runtime_fusion_mean_ms,

                wrong_association_total=wrong_association_total,
                id_switch_total=id_switch_total,
                identity_anomaly_total=_metric(metrics, "association_health.identity_anomaly_total"),
                tracks_with_id_switch=_metric(metrics, "track_lifecycle.tracks_with_id_switch"),
                track_count=_metric(metrics, "track_lifecycle.track_count"),
                track_closed_count=_metric(metrics, "track_lifecycle.closed_count"),
                mean_track_duration_s=_metric(metrics, "track_lifecycle.mean_duration_s"),
                max_track_duration_s=_metric(metrics, "track_lifecycle.max_duration_s"),
                mean_filter_duration_s=_metric(metrics, "track_lifecycle.mean_filter_duration_s"),
                max_filter_duration_s=_metric(metrics, "track_lifecycle.max_filter_duration_s")
            )
        )
    end

    return DataFrame(rows)
end

function _fmt_value(value)
    ismissing(value) && return "n/a"
    if value isa AbstractString
        return String(value)
    end
    v = Float64(value)
    if !isfinite(v)
        return "n/a"
    elseif abs(v - round(v)) < 1e-9
        return string(Int(round(v)))
    else
        return @sprintf("%.2f", v)
    end
end

function _print_table(title::String, df::DataFrame, cols::Vector{Symbol})
    println()
    println(title)
    selected = select(df, cols)
    labels = String.(names(selected))
    widths = [length(label) for label in labels]

    formatted_rows = Vector{Vector{String}}()
    for row in eachrow(selected)
        values = [_fmt_value(row[col]) for col in cols]
        push!(formatted_rows, values)
        for (idx, value) in enumerate(values)
            widths[idx] = max(widths[idx], length(value))
        end
    end

    header = join([rpad(labels[i], widths[i]) for i in eachindex(labels)], "  ")
    println(header)
    println(join([repeat("-", widths[i]) for i in eachindex(widths)], "  "))
    for values in formatted_rows
        println(join([rpad(values[i], widths[i]) for i in eachindex(values)], "  "))
    end
    flush(stdout)
    return nothing
end

function _save_paper_plot(plot_obj, basename::String)::String
    mkpath(PLOT_OUTPUT_ROOT)
    png_path = joinpath(PLOT_OUTPUT_ROOT, "$(basename).png")
    pdf_path = joinpath(PLOT_OUTPUT_ROOT, "$(basename).pdf")
    _progress("saving plot $(basename).png")
    savefig(plot_obj, png_path)
    try
        _progress("saving plot $(basename).pdf")
        savefig(plot_obj, pdf_path)
    catch err
        @warn "Could not save PDF companion for $(basename)" exception=(err, catch_backtrace())
    end
    _progress("finished plot $(basename)")
    return png_path
end

function _plot_observer_initial_formation(cases::Vector{Symbol})::Union{Nothing, String}
    cached = _cached_plot_path("observer_initial_formation_comm_graph")
    if cached !== nothing
        _progress("using cached plot observer_initial_formation_comm_graph")
        return cached
    end
    _progress("building plot observer_initial_formation_comm_graph")
    mkpath(PLOT_OUTPUT_ROOT)

    observer_count = _observer_count_for_plot(cases)
    observer_defs, _, _ = _observer_defs_for_plot(observer_count)
    positions = [_kepler_xyz_km(observer_def, 0.0) for observer_def in observer_defs]
    isempty(positions) && return nothing

    centroid = (
        mean(first.(positions)),
        mean(getindex.(positions, 2)),
        mean(getindex.(positions, 3))
    )
    rel_x = [pos[1] - centroid[1] for pos in positions]
    rel_y = [pos[2] - centroid[2] for pos in positions]

    p = plot(
        xlabel="Relative ECI x [km]",
        ylabel="Relative ECI y [km]",
        aspect_ratio=:equal,
        legend=:topright,
        size=(1150, 950),
        dpi=300,
        fontfamily=PAPER_FONT_FAMILY,
        guidefontsize=PAPER_GUIDE_FONT,
        tickfontsize=PAPER_TICK_FONT,
        legendfontsize=PAPER_LEGEND_FONT,
        left_margin=22mm,
        bottom_margin=20mm,
        right_margin=5mm,
        top_margin=5mm,
        foreground_color_legend=nothing,
        background_color_legend=:white,
        gridalpha=0.18,
        minorgrid=false
    )

    orbit_x_values = copy(rel_x)
    orbit_y_values = copy(rel_y)
    orbit_times_s = collect(range(-55.0, 55.0; length=120))
    orbit_label = true
    for observer_def in observer_defs
        xs, ys, _ = _kepler_xyz_series_km(observer_def, orbit_times_s)
        xrel = xs .- centroid[1]
        yrel = ys .- centroid[2]
        append!(orbit_x_values, xrel)
        append!(orbit_y_values, yrel)
        plot!(
            p,
            xrel,
            yrel;
            color=:gray45,
            alpha=0.42,
            linewidth=2.4,
            label=orbit_label ? "Observer local orbit arcs" : false
        )
        orbit_label = false
    end

    link_label = true
    link_count = 0
    for i in eachindex(positions)
        for j in (i + 1):length(positions)
            dist_ij = sqrt(
                (positions[i][1] - positions[j][1])^2 +
                (positions[i][2] - positions[j][2])^2 +
                (positions[i][3] - positions[j][3])^2
            )
            dist_ij <= COMMUNICATION_RANGE_KM || continue
            link_count += 1
            plot!(
                p,
                [rel_x[i], rel_x[j]],
                [rel_y[i], rel_y[j]];
                color=:steelblue,
                alpha=0.68,
                linewidth=2.5,
                label=link_label ? "Communication links" : false
            )
            link_label = false
        end
    end

    scatter!(
        p,
        rel_x,
        rel_y;
        color=:navy,
        markersize=8.5,
        markerstrokecolor=:white,
        markerstrokewidth=1.4,
        label="Observers"
    )

    label_x_span = max(maximum(rel_x) - minimum(rel_x), 1.0)
    label_y_span = max(maximum(rel_y) - minimum(rel_y), 1.0)
    label_dx = 0.018 * label_x_span
    label_dy = 0.018 * label_y_span
    for idx in eachindex(rel_x)
        annotate!(p, rel_x[idx] + label_dx, rel_y[idx] + label_dy, text(string(idx), 11, :left))
    end

    x_span = max(maximum(orbit_x_values) - minimum(orbit_x_values), 1.0)
    y_span = max(maximum(orbit_y_values) - minimum(orbit_y_values), 1.0)
    margin = max(35.0, 0.18 * max(x_span, y_span))
    x_min = minimum(orbit_x_values) - margin
    x_max = maximum(orbit_x_values) + margin
    y_min = minimum(orbit_y_values) - margin
    y_max = maximum(orbit_y_values) + margin

    scale_y = y_min + 0.10 * (y_max - y_min)
    scale_x0 = x_min + 0.08 * (x_max - x_min)
    scale_x1 = scale_x0 + COMMUNICATION_RANGE_KM
    if scale_x1 < x_max - 0.05 * (x_max - x_min)
        plot!(p, [scale_x0, scale_x1], [scale_y, scale_y]; color=:black, linewidth=3.0, label=false)
        annotate!(
            p,
            (scale_x0 + scale_x1) / 2.0,
            scale_y + 0.045 * (y_max - y_min),
            text("300 km", 13, :center)
        )
    end

    annotate!(
        p,
        x_min + 0.02 * (x_max - x_min),
        y_max - 0.06 * (y_max - y_min),
        text("3D links: $(link_count), R_com = 300 km", 14, :left)
    )
    plot!(p; xlims=(x_min, x_max), ylims=(y_min, y_max))

    return _save_paper_plot(p, "observer_initial_formation_comm_graph")
end

function main()
    _progress("starting navigation comparison metrics printer")
    cases = _selected_cases()
    _progress("selected cases: $(join(String.(cases), ", "))")
    summary = _build_summary(cases)
    _progress("loaded summary metrics")
    tracking_by_observer = _read_tracking_observer_summary(cases)
    _progress("loaded tracking-by-observer metrics")
    mkpath(OUTPUT_ROOT)

    selected_path = joinpath(OUTPUT_ROOT, "comparison_selected_metrics.csv")
    tracking_observer_path = joinpath(OUTPUT_ROOT, "comparison_tracking_by_observer.csv")
    CSV.write(selected_path, summary)
    CSV.write(tracking_observer_path, tracking_by_observer)
    _progress("wrote comparison CSV files")
    observer_formation_plot_path = _plot_observer_initial_formation(cases)
    _progress("finished plot generation")

    println("Navigation comparison metrics")
    println("  output_root: $(OUTPUT_ROOT)")
    println("  selected csv: $(selected_path)")
    println("  tracking by observer csv: $(tracking_observer_path)")
    observer_formation_plot_path !== nothing && println("  observer formation plot: $(observer_formation_plot_path)")

    _print_table(
        "Run Status",
        summary,
        [:case, :status, :baseline_theta_gate_rad, :failure_reason]
    )
    _progress("printed run-status table")

    _print_table(
        "LOS Sensor Error Injection",
        summary,
        [
            :case,
            :sensor_enable_missed_detections,
            :sensor_configured_missed_detection_rate,
            :sensor_missed_detections,
            :sensor_realized_missed_rate_pct,
            :sensor_enable_false_alarms,
            :sensor_configured_false_alarm_rate,
            :sensor_false_alarms,
            :sensor_enable_measurement_bias,
            :sensor_configured_measurement_bias_rad,
            :sensor_biased_measurements
        ]
    )
    _progress("printed LOS sensor error table")

    _print_table(
        "Measurement-to-Track Association",
        summary,
        [
            :case,
            :m2t_committed,
            :m2t_correct,
            :m2t_wrong,
            :m2t_accuracy_pct,
            :m2t_ambig_committed,
            :m2t_ambig_dropped
        ]
    )
    _progress("printed M2T table")

    _print_table(
        "Local M2M / Hypothesis Initialization",
        summary,
        [
            :case,
            :h1_created,
            :h1_to_h2,
            :h1_to_h2_good_pct,
            :h2_to_h3_attempted,
            :los_rate_pass,
            :los_rate_pass_pct,
            :promoted,
            :promoted_good_pct,
            :tracks_created_with_nonreal_measurements
        ]
    )
    _progress("printed local M2M table")

    _print_table(
        "Cross-Observer M2M / IOD",
        summary,
        [
            :case,
            :xm2m_candidate_pairs,
            :xm2m_gate_pass,
            :xm2m_gate_pass_good_pct,
            :xm2m_gate_pass_mixed_target,
            :xm2m_selected_pairs,
            :xm2m_selected_good_pct,
            :xm2m_selected_mixed_target,
            :xm2m_iod_groups,
            :xm2m_iod_groups_good_pct,
            :iod_position_cov_gate_threshold_rms_m,
            :iod_position_cov_gate_evaluated,
            :iod_position_cov_gate_rejected,
            :iod_position_cov_gate_rejected_pct,
            :iod_position_cov_gate_rejected_same_target,
            :iod_position_cov_gate_rejected_mixed_target,
            :iod_one_step_validation_enabled,
            :iod_validation_threshold_d2,
            :iod_validation_attempted,
            :iod_validation_confirmed,
            :iod_validation_rejected,
            :iod_validation_confirmed_pct,
            :iod_validation_confirmed_same_target,
            :iod_validation_confirmed_mixed_target,
            :iod_validation_rejected_same_target,
            :iod_validation_rejected_mixed_target,
            :iod_validation_no_measure,
            :iod_validation_pending_end,
            :xm2m_iod_initialized,
            :xm2m_iod_initialized_good_pct,
            :xm2m_iod_initialized_mixed_target
        ]
    )
    _progress("printed cross-observer M2M/IOD table")

    _print_table(
        "Track-to-Track / Consensus",
        summary,
        [
            :case,
            :t2t_attempted,
            :t2t_committed,
            :t2t_commit_rate_pct,
            :t2t_accuracy_pct,
            :t2t_ratio_fail,
            :t2t_mutual_fail,
            :t2t_component_conflict_rejected,
            :consensus_groups,
            :consensus_good_pct,
            :consensus_mixed_target
        ]
    )
    _progress("printed T2T/consensus table")

    _print_table(
        "Tracking Windows Overall",
        summary,
        [
            :case,
            :tracking_possible_windows,
            :tracking_tracked_windows,
            :tracking_coverage_pct,
            :tracking_successful_windows_under_1km,
            :tracking_success_rate_possible_pct,
            :tracking_mean_error_all_m,
            :tracking_mean_error_good_m,
            :tracking_sample_weighted_mean_error_m,
            :tracking_sample_rmse_error_m,
            :tracking_success_threshold_m
        ]
    )
    _progress("printed tracking-windows overall table")

    _print_table(
        "Unique Object Coverage",
        summary,
        [
            :case,
            :detected_unique_targets,
            :jointly_detected_unique_targets,
            :tracked_unique_targets,
            :successful_tracked_unique_targets,
            :tracked_unique_over_detected_pct,
            :successful_unique_over_detected_pct,
            :tracked_unique_over_jointly_detected_pct,
            :successful_unique_over_jointly_detected_pct
        ]
    )
    _progress("printed unique-object coverage table")

    _print_table(
        "Runtime Per Navigation Epoch",
        summary,
        [
            :case,
            :runtime_timed_epochs,
            :runtime_total_epoch_mean_ms,
            :runtime_total_epoch_median_ms,
            :runtime_total_epoch_p95_ms,
            :runtime_local_da_mean_ms,
            :runtime_cross_da_mean_ms,
            :runtime_filter_mean_ms
        ]
    )
    _progress("printed runtime table")

    _print_table(
        "Association Health / Estimate Duration",
        summary,
        [
            :case,
            :wrong_association_total,
            :id_switch_total,
            :identity_anomaly_total,
            :tracks_created_with_nonreal_measurements,
            :tracks_with_id_switch,
            :track_count,
            :track_closed_count,
            :mean_estimate_duration_s,
            :max_estimate_duration_s,
            :mean_successful_estimate_duration_s
        ]
    )
    _progress("printed association-health table")

    if nrow(tracking_by_observer) > 0
        _progress("printing tracking-windows-by-observer table")
        _print_table(
            "Tracking Windows By Observer",
            tracking_by_observer,
            [
                :case,
                :observer,
                :possible_windows,
                :tracked_windows,
                :tracking_coverage_pct,
                :successful_windows_under_1km,
                :success_rate_possible_pct,
                :detected_unique_targets,
                :tracked_unique_targets,
                :successful_unique_targets,
                :mean_error_tracked_windows_m,
                :mean_error_successful_windows_m
            ]
        )
        _progress("printed tracking-windows-by-observer table")
    end

    _progress("finished navigation comparison metrics printer")

    return nothing
end

function _run_main_with_diagnostics()::Nothing
    try
        main()
    catch err
        bt = catch_backtrace()
        log_path = joinpath(OUTPUT_ROOT, "print_navigation_comparison_metrics_error.log")
        mkpath(dirname(log_path))

        open(log_path, "w") do io
            println(io, "Error while executing $(@__FILE__).")
            showerror(io, err, bt)
            println(io)
        end

        println(stderr, "Error while executing $(@__FILE__).")
        showerror(stderr, err, bt)
        println(stderr)
        println(stderr, "Full traceback saved to: $(log_path)")
        rethrow()
    end

    return nothing
end

_run_main_with_diagnostics()
