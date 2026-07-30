module MonteCarloIODGeometryAnalysis

using CSV
using DataFrames
using Dates
using Statistics

if !isdefined(@__MODULE__, :NavigationPaths)
    include(joinpath(@__DIR__, "..", "paths.jl"))
end
using .NavigationPaths: REPO_ROOT, env_path, resolve_repo_path

const DETECTION_RANGE_M = 500.0e3
const OBSERVER_COUNT = 16
const NAVIGATION_RATE_SEC = 5.0

@inline _truthy(value::AbstractString)::Bool =
    lowercase(strip(value)) in ("1", "true", "yes", "y", "on")

@inline function _finite_values(values)
    return Float64[Float64(value) for value in values if !ismissing(value) && isfinite(Float64(value))]
end

@inline function _percentile(values::Vector{Float64}, probability::Float64)
    isempty(values) && return NaN
    return quantile(values, probability)
end

function _validation_lookup(diagnostics::DataFrame)
    lookup = Dict{Tuple{Int, Int, Int, Int}, NamedTuple}()
    for row in eachrow(diagnostics)
        String(row.stage) == "next_step_validation" || continue
        key = (
            Int(row.realization), Int(row.observer), Int(row.slot),
            round(Int, 1_000.0 * Float64(row.time_s))
        )
        lookup[key] = (
            outcome=String(row.outcome),
            validation_d2=Float64(row.validation_d2),
            validation_position_error_m=Float64(row.position_error_m)
        )
    end
    return lookup
end

function _attach_validation!(events::DataFrame, diagnostics::DataFrame)
    lookup = _validation_lookup(diagnostics)
    last_iod_time = Dict(
        Int(group.realization[1]) => maximum(Float64.(group.time_s))
        for group in groupby(events, :realization)
    )
    final_outcomes = String[]
    validation_d2 = Float64[]
    validation_errors = Float64[]
    for row in eachrow(events)
        if String(row.outcome) == "rejected"
            push!(final_outcomes, "covariance_rejected")
            push!(validation_d2, NaN)
            push!(validation_errors, NaN)
            continue
        end
        validation_time_ms = round(Int, 1_000.0 * (Float64(row.time_s) + NAVIGATION_RATE_SEC))
        key = (Int(row.realization), Int(row.observer), Int(row.slot), validation_time_ms)
        validation = get(lookup, key, nothing)
        if validation === nothing
            pending_at_end = Float64(row.time_s) >=
                get(last_iod_time, Int(row.realization), Inf) - NAVIGATION_RATE_SEC / 2.0
            push!(final_outcomes, pending_at_end ? "pending_at_mission_end" : "validation_record_missing")
            push!(validation_d2, NaN)
            push!(validation_errors, NaN)
        else
            push!(final_outcomes, validation.outcome)
            push!(validation_d2, validation.validation_d2)
            push!(validation_errors, validation.validation_position_error_m)
        end
    end
    events.final_outcome = final_outcomes
    events.next_step_validation_d2 = validation_d2
    events.next_step_position_error_m = validation_errors
    return events
end

function _window_proxy_lookup(windows::DataFrame, case::String)
    lookup = Dict{Tuple{Int, Int}, Vector{Tuple{Float64, Float64}}}()
    nrow(windows) == 0 && return lookup
    selected = :case in propertynames(windows) ? windows[windows.case .== case, :] : windows
    for row in eachrow(selected)
        key = (Int(row.realization), Int(row.observer))
        push!(
            get!(lookup, key, Tuple{Float64, Float64}[]),
            (Float64(row.window_start_t_s), Float64(row.window_end_t_s))
        )
    end
    return lookup
end

function _attach_window_proxy!(events::DataFrame, windows::DataFrame, case::String)
    lookup = _window_proxy_lookup(windows, case)
    events.joint_visibility_window_proxy_count = [
        count(interval -> interval[1] <= Float64(row.time_s) <= interval[2],
            get(lookup, (Int(row.realization), Int(row.observer)), Tuple{Float64, Float64}[]))
        for row in eachrow(events)
    ]
    return events
end

function _ensure_geometry_columns!(events::DataFrame)
    definitions = (
        :visible_target_count => fill(-1, nrow(events)),
        :minimum_angular_separation_deg => fill(NaN, nrow(events)),
        :closest_target_1 => fill(0, nrow(events)),
        :closest_target_2 => fill(0, nrow(events)),
        :reference_minimum_angular_separation_deg => fill(NaN, nrow(events)),
        :reference_closest_target => fill(0, nrow(events)),
        :visible_target_ids => fill("", nrow(events)),
        :minimum_3nn_radius_deg => fill(NaN, nrow(events)),
        :maximum_3nn_density_targets_per_sr => fill(NaN, nrow(events)),
        :reference_3nn_radius_deg => fill(NaN, nrow(events)),
        :reference_3nn_density_targets_per_sr => fill(NaN, nrow(events)),
        :geometry_source => fill("unavailable", nrow(events)),
        :composition => ifelse.(events.same_target, "same_target", "mixed_real_target"),
        :group_observers => fill("", nrow(events)),
        :group_target_ids => fill("", nrow(events))
    )
    for (name, values) in definitions
        name in propertynames(events) || (events[!, name] = values)
    end
    return events
end

function _truth_replay_output_is_valid(path::String, event_ids::Vector{Int})::Bool
    isfile(path) || return false
    try
        table = CSV.read(path, DataFrame)
        return nrow(table) == length(event_ids) && sort(Int.(table.event_id)) == sort(event_ids)
    catch
        return false
    end
end

function _replay_missing_geometry_truth!(
    events::DataFrame,
    target_populations::DataFrame,
    output_root::String;
    max_realizations::Int=0
)
    _ensure_geometry_columns!(events)
    missing_rows = findall(row -> Int(row.visible_target_count) < 0, eachrow(events))
    isempty(missing_rows) && return 0

    realizations = sort!(unique(Int(events.realization[row]) for row in missing_rows))
    if max_realizations > 0
        realizations = realizations[1:min(max_realizations, length(realizations))]
    end
    population = target_populations
    if nrow(population) > 0 && :case in propertynames(population) && any(population.case .== "proposed")
        population = population[population.case .== "proposed", :]
    end

    cache_root = joinpath(output_root, "iod_geometry_analysis", "truth_replay_cache")
    mkpath(cache_root)
    simulator_path = joinpath(@__DIR__, "..", "..", "navigation.jl")
    julia_executable = joinpath(Sys.BINDIR, Base.julia_exename())
    heap_hint = strip(get(ENV, "SPACEAGORA_JULIA_HEAP_SIZE_HINT", "4G"))
    resume = _truthy(get(ENV, "SPACEAGORA_MC_IOD_GEOMETRY_TRUTH_RESUME", "true"))
    status_path = joinpath(cache_root, "truth_replay_status.csv")
    status = DataFrame()
    replayed_event_count = 0

    # Use the same robust, sequential child-process pattern as the nominal
    # Monte Carlo runner. The child runs truth dynamics and the geometry sampler
    # only; no measurement generation, association, IOD, or navigation filter.
    for (realization_index, realization) in enumerate(realizations)
        row_indices = [row for row in missing_rows if Int(events.realization[row]) == realization]
        event_ids = Int.(row_indices)
        input_path = joinpath(cache_root, "run_$(lpad(realization, 4, '0'))_events.csv")
        output_path = joinpath(cache_root, "run_$(lpad(realization, 4, '0'))_truth_geometry.csv")
        log_path = joinpath(cache_root, "run_$(lpad(realization, 4, '0')).log")
        replay_input = DataFrame(
            event_id=event_ids,
            time_s=Float64.(events.time_s[row_indices]),
            observer=Int.(events.observer[row_indices]),
            reference_target=Int.(events.reference_target[row_indices])
        )
        CSV.write(input_path, replay_input)

        resumed = resume && _truth_replay_output_is_valid(output_path, event_ids)
        exit_code = 0
        elapsed_sec = 0.0
        run_status = "resumed"

        if !resumed
            scenario_seed = Int(events.scenario_seed[first(row_indices)])
            target_count = if nrow(population) > 0 && :realization in propertynames(population)
                count(==(realization), Int.(population.realization))
            else
                300
            end
            target_count > 0 || error("No saved target population for realization $(realization)")
            mission_time = maximum(Float64.(events.time_s[row_indices])) + NAVIGATION_RATE_SEC
            child_environment = copy(ENV)
            child_environment["SPACEAGORA_TRUTH_GEOMETRY_REPLAY_ONLY"] = "true"
            child_environment["SPACEAGORA_TRUTH_GEOMETRY_REPLAY_INPUT"] = input_path
            child_environment["SPACEAGORA_TRUTH_GEOMETRY_REPLAY_OUTPUT"] = output_path
            child_environment["SPACEAGORA_SCENARIO_SEED"] = string(scenario_seed)
            child_environment["SPACEAGORA_SENSOR_SEED"] = string(scenario_seed + 1)
            child_environment["SPACEAGORA_N_CLUSTER_TARGETS"] = string(target_count)
            child_environment["SPACEAGORA_MISSION_TIME_SEC"] = string(mission_time)
            child_environment["SPACEAGORA_NAV_CASE"] = "proposed"
            child_environment["SPACEAGORA_NAV_OUTPUT_LABEL"] = "truth_geometry_replay"
            child_environment["SPACEAGORA_COMPARISON_OUTPUT"] = cache_root
            child_environment["SPACEAGORA_SAVE_SIMULATION_RESULTS"] = "false"
            child_environment["SPACEAGORA_SAVE_TARGET_ESTIMATE_FIELDS"] = "false"
            child_environment["SPACEAGORA_SAVE_COMPARISON_DETAILED_TABLES"] = "false"
            child_environment["SPACEAGORA_ENABLE_NAV_TIMING"] = "false"
            command = isempty(heap_hint) ?
                setenv(`$julia_executable --project=$REPO_ROOT $simulator_path`, child_environment) :
                setenv(`$julia_executable --heap-size-hint=$heap_hint --project=$REPO_ROOT $simulator_path`, child_environment)
            println(
                "  truth geometry replay: realization $(realization) started " *
                "($(realization_index)/$(length(realizations)))"
            )
            flush(stdout)
            start_time = time()
            open(log_path, "w") do log_io
                process = run(
                    pipeline(ignorestatus(command); stdout=log_io, stderr=log_io);
                    wait=false
                )
                last_heartbeat = time()
                while process_running(process)
                    sleep(1.0)
                    if time() - last_heartbeat >= 30.0
                        println(
                            "    realization $(realization): still propagating " *
                            "($(round((time() - start_time) / 60; digits=1)) min)"
                        )
                        flush(stdout)
                        last_heartbeat = time()
                    end
                end
                wait(process)
                exit_code = try
                    Int(process.exitcode)
                catch
                    success(process) ? 0 : 1
                end
            end
            elapsed_sec = time() - start_time
            valid_output = _truth_replay_output_is_valid(output_path, event_ids)
            if exit_code == 0 && valid_output
                run_status = "ok"
                println(
                    "    realization $(realization): completed in " *
                    "$(round(elapsed_sec; digits=1)) s"
                )
            else
                run_status = "failed"
                log_tail = if isfile(log_path)
                    lines = readlines(log_path)
                    join(lines[max(1, length(lines) - 3):end], " | ")
                else
                    "log not created"
                end
                println(
                    "    realization $(realization): FAILED (exit_code=$(exit_code)); " *
                    "$(log_tail)"
                )
            end
            flush(stdout)
        else
            println(
                "  truth geometry replay: realization $(realization) resumed " *
                "($(realization_index)/$(length(realizations)))"
            )
            flush(stdout)
        end

        if run_status in ("ok", "resumed")
            replay = CSV.read(output_path, DataFrame)
            replay_by_id = Dict(Int(row.event_id) => row for row in eachrow(replay))
            for row_index in row_indices
                geometry = replay_by_id[row_index]
                events.visible_target_count[row_index] = Int(geometry.visible_target_count)
                events.minimum_angular_separation_deg[row_index] = Float64(geometry.minimum_angular_separation_deg)
                events.closest_target_1[row_index] = Int(geometry.closest_target_1)
                events.closest_target_2[row_index] = Int(geometry.closest_target_2)
                events.reference_minimum_angular_separation_deg[row_index] =
                    Float64(geometry.reference_minimum_angular_separation_deg)
                events.reference_closest_target[row_index] = Int(geometry.reference_closest_target)
                events.visible_target_ids[row_index] = String(geometry.visible_target_ids)
                events.minimum_3nn_radius_deg[row_index] = Float64(geometry.minimum_3nn_radius_deg)
                events.maximum_3nn_density_targets_per_sr[row_index] =
                    Float64(geometry.maximum_3nn_density_targets_per_sr)
                events.reference_3nn_radius_deg[row_index] = Float64(geometry.reference_3nn_radius_deg)
                events.reference_3nn_density_targets_per_sr[row_index] =
                    Float64(geometry.reference_3nn_density_targets_per_sr)
                events.geometry_source[row_index] = "full_truth_dynamics_replay"
            end
            replayed_event_count += length(row_indices)
        end

        status_row = (
            realization=realization,
            status=run_status,
            exit_code=exit_code,
            elapsed_sec=elapsed_sec,
            event_count=length(row_indices),
            output_path=output_path,
            log_path=log_path
        )
        if ncol(status) == 0
            status = DataFrame([status_row])
        else
            push!(status, status_row)
        end
        CSV.write(status_path, status)
    end
    return replayed_event_count
end

function _mark_direct_geometry!(events::DataFrame)
    _ensure_geometry_columns!(events)
    for row_index in 1:nrow(events)
        if Int(events.visible_target_count[row_index]) >= 0
            events.geometry_source[row_index] = "direct_truth_geometry"
        end
    end
    return events
end

function _class_summary(events::DataFrame)::DataFrame
    rejected = in.(events.final_outcome, Ref(Set((
        "covariance_rejected", "rejected_score", "rejected_invalid_state"
    ))))
    subsets = Pair{String, Vector{Int}}[
        "all" => collect(1:nrow(events)),
        "same_target" => findall(events.same_target),
        "mixed_target" => findall(.!events.same_target),
        "same_target_promoted" => findall(events.same_target .& (events.final_outcome .== "promoted")),
        "same_target_rejected" => findall(events.same_target .& rejected),
        "same_target_no_measurement" => findall(events.same_target .& (events.final_outcome .== "no_measurement")),
        "same_target_pending_at_end" => findall(events.same_target .& (events.final_outcome .== "pending_at_mission_end")),
        "mixed_target_promoted" => findall((.!events.same_target) .& (events.final_outcome .== "promoted")),
        "mixed_target_rejected" => findall((.!events.same_target) .& rejected),
        "mixed_target_no_measurement" => findall((.!events.same_target) .& (events.final_outcome .== "no_measurement")),
        "mixed_target_pending_at_end" => findall((.!events.same_target) .& (events.final_outcome .== "pending_at_mission_end"))
    ]
    rows = NamedTuple[]
    for (label, indices) in subsets
        isempty(indices) && continue
        visible = _finite_values(events.visible_target_count[indices])
        filter!(>=(0.0), visible)
        separation = _finite_values(events.minimum_angular_separation_deg[indices])
        reference_separation = _finite_values(events.reference_minimum_angular_separation_deg[indices])
        cluster_radius = _finite_values(events.minimum_3nn_radius_deg[indices])
        cluster_density = _finite_values(events.maximum_3nn_density_targets_per_sr[indices])
        reference_density = _finite_values(events.reference_3nn_density_targets_per_sr[indices])
        position_error = _finite_values(events.position_error_m[indices])
        formal_sigma = _finite_values(events.position_rms_std_m[indices])
        push!(rows, (
            event_class=label,
            count=length(indices),
            geometry_event_count=length(visible),
            visible_count_mean=isempty(visible) ? NaN : mean(visible),
            visible_count_median=isempty(visible) ? NaN : median(visible),
            visible_count_p05=_percentile(visible, 0.05),
            visible_count_p95=_percentile(visible, 0.95),
            minimum_separation_median_deg=isempty(separation) ? NaN : median(separation),
            minimum_separation_p05_deg=_percentile(separation, 0.05),
            reference_separation_median_deg=isempty(reference_separation) ? NaN : median(reference_separation),
            reference_separation_p05_deg=_percentile(reference_separation, 0.05),
            minimum_3nn_radius_median_deg=isempty(cluster_radius) ? NaN : median(cluster_radius),
            maximum_3nn_density_median_targets_per_sr=isempty(cluster_density) ? NaN : median(cluster_density),
            reference_3nn_density_median_targets_per_sr=isempty(reference_density) ? NaN : median(reference_density),
            position_error_median_m=isempty(position_error) ? NaN : median(position_error),
            position_error_p95_m=_percentile(position_error, 0.95),
            formal_position_rms_median_m=isempty(formal_sigma) ? NaN : median(formal_sigma)
        ))
    end
    return DataFrame(rows)
end

function _cluster_density_summary(events::DataFrame)::DataFrame
    finite_rows = findall(isfinite, Float64.(events.maximum_3nn_density_targets_per_sr))
    isempty(finite_rows) && return DataFrame()
    values = Float64.(events.maximum_3nn_density_targets_per_sr[finite_rows])
    q25, q50, q75 = quantile(values, (0.25, 0.50, 0.75))
    labels = [
        value <= q25 ? "Q1: lowest density" :
        value <= q50 ? "Q2" : value <= q75 ? "Q3" : "Q4: highest density"
        for value in values
    ]
    usable = events[finite_rows, :]
    usable.cluster_density_bin = labels
    order = Dict(
        "Q1: lowest density" => 1, "Q2" => 2,
        "Q3" => 3, "Q4: highest density" => 4
    )
    rows = NamedTuple[]
    for group in groupby(usable, :cluster_density_bin)
        mixed = .!group.same_target
        mixed_promoted = mixed .& (group.final_outcome .== "promoted")
        group_values = Float64.(group.maximum_3nn_density_targets_per_sr)
        push!(rows, (
            cluster_density_bin=String(group.cluster_density_bin[1]),
            lower_bound_targets_per_sr=minimum(group_values),
            upper_bound_targets_per_sr=maximum(group_values),
            total_iod_events=nrow(group),
            mixed_iod_events=count(mixed),
            mixed_iod_rate_pct=100.0 * count(mixed) / nrow(group),
            mixed_iod_promoted=count(mixed_promoted),
            promoted_fraction_of_mixed_pct=count(mixed) > 0 ?
                100.0 * count(mixed_promoted) / count(mixed) : NaN
        ))
    end
    result = DataFrame(rows)
    sort!(result, :cluster_density_bin; by=value -> order[String(value)])
    return result
end

@inline function _visible_bin(value::Int)::String
    value <= 5 && return "0--5"
    value <= 8 && return "6--8"
    value <= 11 && return "9--11"
    value <= 14 && return "12--14"
    return ">=15"
end

function _visible_count_summary(events::DataFrame)::DataFrame
    usable = events[events.visible_target_count .>= 0, :]
    nrow(usable) == 0 && return DataFrame()
    usable.visible_count_bin = _visible_bin.(Int.(usable.visible_target_count))
    rows = NamedTuple[]
    order = Dict("0--5" => 1, "6--8" => 2, "9--11" => 3, "12--14" => 4, ">=15" => 5)
    for group in groupby(usable, :visible_count_bin)
        mixed = .!group.same_target
        mixed_promoted = mixed .& (group.final_outcome .== "promoted")
        push!(rows, (
            visible_count_bin=String(group.visible_count_bin[1]),
            total_iod_events=nrow(group),
            mixed_iod_events=count(mixed),
            mixed_iod_rate_pct=100.0 * count(mixed) / nrow(group),
            mixed_iod_promoted=count(mixed_promoted),
            promoted_fraction_of_mixed_pct=count(mixed) > 0 ?
                100.0 * count(mixed_promoted) / count(mixed) : NaN
        ))
    end
    result = DataFrame(rows)
    sort!(result, :visible_count_bin; by=value -> order[String(value)])
    return result
end

function _separation_summary(events::DataFrame)::DataFrame
    finite_rows = findall(isfinite, Float64.(events.minimum_angular_separation_deg))
    isempty(finite_rows) && return DataFrame()
    values = Float64.(events.minimum_angular_separation_deg[finite_rows])
    q25, q50, q75 = quantile(values, (0.25, 0.50, 0.75))
    labels = [
        value <= q25 ? "Q1: smallest separation" :
        value <= q50 ? "Q2" : value <= q75 ? "Q3" : "Q4: largest separation"
        for value in values
    ]
    usable = events[finite_rows, :]
    usable.separation_bin = labels
    order = Dict(
        "Q1: smallest separation" => 1, "Q2" => 2,
        "Q3" => 3, "Q4: largest separation" => 4
    )
    rows = NamedTuple[]
    for group in groupby(usable, :separation_bin)
        mixed = .!group.same_target
        mixed_promoted = mixed .& (group.final_outcome .== "promoted")
        group_values = Float64.(group.minimum_angular_separation_deg)
        push!(rows, (
            separation_bin=String(group.separation_bin[1]),
            lower_bound_deg=minimum(group_values),
            upper_bound_deg=maximum(group_values),
            total_iod_events=nrow(group),
            mixed_iod_events=count(mixed),
            mixed_iod_rate_pct=100.0 * count(mixed) / nrow(group),
            mixed_iod_promoted=count(mixed_promoted),
            promoted_fraction_of_mixed_pct=count(mixed) > 0 ?
                100.0 * count(mixed_promoted) / count(mixed) : NaN
        ))
    end
    result = DataFrame(rows)
    sort!(result, :separation_bin; by=value -> order[String(value)])
    return result
end

function _observer_summary(events::DataFrame)::DataFrame
    rows = NamedTuple[]
    for group in groupby(events, :observer)
        mixed = .!group.same_target
        mixed_promoted = mixed .& (group.final_outcome .== "promoted")
        visible_mixed = _finite_values(group.visible_target_count[mixed])
        filter!(>=(0.0), visible_mixed)
        separation_mixed = _finite_values(group.minimum_angular_separation_deg[mixed])
        push!(rows, (
            observer=Int(group.observer[1]),
            total_iod_events=nrow(group),
            mixed_iod_events=count(mixed),
            mixed_iod_rate_pct=100.0 * count(mixed) / nrow(group),
            mixed_iod_promoted=count(mixed_promoted),
            mixed_visible_count_median=isempty(visible_mixed) ? NaN : median(visible_mixed),
            mixed_minimum_separation_median_deg=isempty(separation_mixed) ? NaN : median(separation_mixed)
        ))
    end
    result = DataFrame(rows)
    sort!(result, :observer)
    return result
end

function _outcome_summary(events::DataFrame)::DataFrame
    categorized = DataFrame(events)
    categorized.target_class = ifelse.(categorized.same_target, "same_target", "mixed_target")
    rows = NamedTuple[]
    for group in groupby(categorized, [:target_class, :outcome, :final_outcome])
        visible = _finite_values(group.visible_target_count)
        filter!(>=(0.0), visible)
        separation = _finite_values(group.minimum_angular_separation_deg)
        reference_separation = _finite_values(group.reference_minimum_angular_separation_deg)
        error_values = _finite_values(group.position_error_m)
        push!(rows, (
            target_class=String(group.target_class[1]),
            covariance_gate_outcome=String(group.outcome[1]),
            final_outcome=String(group.final_outcome[1]),
            count=nrow(group),
            visible_count_median=isempty(visible) ? NaN : median(visible),
            minimum_separation_median_deg=isempty(separation) ? NaN : median(separation),
            reference_separation_median_deg=isempty(reference_separation) ? NaN : median(reference_separation),
            position_error_median_m=isempty(error_values) ? NaN : median(error_values)
        ))
    end
    result = DataFrame(rows)
    sort!(result, [:target_class, :covariance_gate_outcome, :final_outcome])
    return result
end

function _realization_summary(events::DataFrame)::DataFrame
    rows = NamedTuple[]
    for group in groupby(events, :realization)
        mixed = .!group.same_target
        same = group.same_target
        mixed_promoted = mixed .& (group.final_outcome .== "promoted")
        mixed_visible = _finite_values(group.visible_target_count[mixed])
        same_visible = _finite_values(group.visible_target_count[same])
        filter!(>=(0.0), mixed_visible)
        filter!(>=(0.0), same_visible)
        mixed_separation = _finite_values(group.minimum_angular_separation_deg[mixed])
        same_separation = _finite_values(group.minimum_angular_separation_deg[same])
        push!(rows, (
            realization=Int(group.realization[1]),
            scenario_seed=Int(group.scenario_seed[1]),
            total_iod_events=nrow(group),
            mixed_iod_events=count(mixed),
            mixed_iod_rate_pct=100.0 * count(mixed) / nrow(group),
            mixed_iod_promoted=count(mixed_promoted),
            same_visible_count_median=isempty(same_visible) ? NaN : median(same_visible),
            mixed_visible_count_median=isempty(mixed_visible) ? NaN : median(mixed_visible),
            same_minimum_separation_median_deg=isempty(same_separation) ? NaN : median(same_separation),
            mixed_minimum_separation_median_deg=isempty(mixed_separation) ? NaN : median(mixed_separation)
        ))
    end
    result = DataFrame(rows)
    sort!(result, :realization)
    return result
end

function _write_metadata(
    output_directory::String,
    event_count::Int,
    reconstructed_count::Int,
    direct_count::Int
)
    metadata = DataFrame(
        field=[
            "created_at", "event_count", "direct_truth_event_count",
            "replayed_event_count", "replay_model", "replay_purpose",
            "replay_limitation", "legacy_group_membership",
            "detection_range_m", "observer_count"
        ],
        value=[
            string(now()), string(event_count), string(direct_count), string(reconstructed_count),
            "same full truth dynamics as the nominal simulation: central gravity, J2--J4, Sun/Moon and SRP",
            "exact truth-geometry recovery for campaigns that predate direct IOD geometry logging",
            "measurement-space quantities and legacy IOD group membership are not reconstructed",
            "not recoverable; generator observer, reference target, and reconstructed visible targets remain available",
            string(DETECTION_RANGE_M), string(OBSERVER_COUNT)
        ]
    )
    CSV.write(joinpath(output_directory, "iod_geometry_metadata.csv"), metadata)
    return nothing
end

function run_iod_geometry_analysis(
    output_root::AbstractString;
    case::String="proposed",
    replay_missing::Bool=true,
    max_realizations::Int=0
)
    root = resolve_repo_path(output_root)
    diagnostics_path = joinpath(root, "mc_iod_diagnostics.csv")
    isfile(diagnostics_path) || error("Missing Monte Carlo IOD diagnostics: $(diagnostics_path)")
    diagnostics = CSV.read(diagnostics_path, DataFrame)
    selected = diagnostics[
        (diagnostics.case .== case) .& (diagnostics.stage .== "covariance_gate"),
        :
    ]
    output_directory = joinpath(root, "iod_geometry_analysis")
    mkpath(output_directory)
    if nrow(selected) == 0
        CSV.write(joinpath(output_directory, "iod_event_geometry.csv"), selected)
        _write_metadata(output_directory, 0, 0, 0)
        println("No covariance-gate IOD events found for case=$(case)")
        println("IOD geometry analysis written to $(output_directory)")
        return output_directory
    end
    events = DataFrame(selected)
    if max_realizations > 0
        selected_realizations = sort!(unique(Int.(events.realization)))
        selected_realizations = selected_realizations[1:min(max_realizations, length(selected_realizations))]
        events = events[in.(Int.(events.realization), Ref(Set(selected_realizations))), :]
    end
    _ensure_geometry_columns!(events)
    _mark_direct_geometry!(events)
    _attach_validation!(events, diagnostics[diagnostics.case .== case, :])

    windows_path = joinpath(root, "mc_tracking_windows.csv")
    windows = isfile(windows_path) ? CSV.read(windows_path, DataFrame) : DataFrame()
    _attach_window_proxy!(events, windows, case)

    reconstructed_count = 0
    if replay_missing
        populations_path = joinpath(root, "mc_target_populations.csv")
        populations = isfile(populations_path) ? CSV.read(populations_path, DataFrame) : DataFrame()
        reconstructed_count = _replay_missing_geometry_truth!(
            events,
            populations,
            root;
            max_realizations=0
        )
    end
    direct_count = count(==("direct_truth_geometry"), events.geometry_source)

    sort!(events, [:realization, :time_s, :observer, :slot])
    CSV.write(joinpath(output_directory, "iod_event_geometry.csv"), events)
    CSV.write(joinpath(output_directory, "iod_geometry_by_class.csv"), _class_summary(events))
    CSV.write(joinpath(output_directory, "iod_geometry_by_outcome.csv"), _outcome_summary(events))
    CSV.write(joinpath(output_directory, "iod_mixed_rate_by_visible_count.csv"), _visible_count_summary(events))
    CSV.write(joinpath(output_directory, "iod_mixed_rate_by_angular_separation.csv"), _separation_summary(events))
    CSV.write(joinpath(output_directory, "iod_mixed_rate_by_cluster_density.csv"), _cluster_density_summary(events))
    CSV.write(joinpath(output_directory, "iod_geometry_by_observer.csv"), _observer_summary(events))
    CSV.write(joinpath(output_directory, "iod_geometry_by_realization.csv"), _realization_summary(events))
    _write_metadata(output_directory, nrow(events), reconstructed_count, direct_count)

    println("IOD geometry analysis written to $(output_directory)")
    println("  covariance-gate events: $(nrow(events))")
    println("  direct truth geometry: $(direct_count)")
    println("  full-dynamics truth replay: $(reconstructed_count)")
    unresolved = count(==("unavailable"), events.geometry_source)
    unresolved > 0 && println("  unavailable geometry: $(unresolved)")
    return output_directory
end

function main()
    output_root = env_path(
        "SPACEAGORA_MC_OUTPUT",
        "nominal"
    )
    replay_missing = _truthy(get(ENV, "SPACEAGORA_MC_IOD_GEOMETRY_REPLAY", "true"))
    max_realizations = parse(Int, get(ENV, "SPACEAGORA_MC_IOD_GEOMETRY_MAX_RUNS", "0"))
    run_iod_geometry_analysis(
        output_root;
        case=get(ENV, "SPACEAGORA_MC_IOD_GEOMETRY_CASE", "proposed"),
        replay_missing=replay_missing,
        max_realizations=max_realizations
    )
end

end # module MonteCarloIODGeometryAnalysis

if abspath(PROGRAM_FILE) == @__FILE__
    MonteCarloIODGeometryAnalysis.main()
end
