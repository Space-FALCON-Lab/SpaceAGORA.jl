"""
The post-processing tool. It reads saved Feather trajectories and encounter/laser records, then exports:

Simulation data as CSV.
Encounter start/end times, durations, and laser-on durations.
Spacecraft states during encounters.
Minimum/maximum relative speed, range rate, and durations.
"""

using Arrow, CSV, DataFrames
using LinearAlgebra

function encounter_csv_directory(feather_path)
    directory = isdir(feather_path) ? abspath(feather_path) : dirname(abspath(feather_path))
    root = basename(dirname(directory)) == "feather" ? dirname(dirname(directory)) :
        joinpath(@__DIR__, "output")
    return joinpath(root, "CSV", basename(directory))
end

function contact_windows(times, ranges, maximum_range)
    windows = Tuple{Float64,Float64}[]
    for sample in 1:length(times)-1
        left_range, right_range = ranges[sample], ranges[sample+1]
        min(left_range, right_range) > maximum_range && continue
        start_time, end_time = times[sample], times[sample+1]
        if left_range > maximum_range
            start_time += (end_time - start_time) *
                          (left_range - maximum_range) / (left_range - right_range)
        elseif right_range > maximum_range
            end_time = start_time + (end_time - start_time) *
                       (maximum_range - left_range) / (right_range - left_range)
        end
        if !isempty(windows) && last(windows)[2] == start_time
            windows[end] = (last(windows)[1], end_time)
        else
            push!(windows, (start_time, end_time))
        end
    end
    return windows
end

function encounter_state_at(times, states, time)
    sample = min(searchsortedlast(times, time), length(times)-1)
    fraction = (time - times[sample]) / (times[sample+1] - times[sample])
    return (1 - fraction) .* states[:, sample] .+ fraction .* states[:, sample+1]
end

function pair_motion(state_a, state_b)
    relative_position = state_b[1:3] - state_a[1:3]
    relative_velocity = state_b[4:6] - state_a[4:6]
    range = norm(relative_position)
    range_rate = iszero(range) ? missing : dot(relative_position, relative_velocity) / range
    return range, norm(relative_velocity), range_rate
end

"""
    extract_feather_encounters(feather_path; maximum_range_m, output_dir)

Read a three-file scenario bundle and use recorded encounter/laser intervals.
Write direct simulation tables under output_dir/sim_output and derived reports
under output_dir/analysis. The output_dir keyword specifies the scenario root.
The range threshold is read from metadata unless supplied explicitly. Legacy
single-file input remains supported with interpolated contacts and estimated
laser durations; trajectory.feather requires both associated interval files.
"""
function extract_feather_encounters(feather_path;
                                    maximum_range_m=nothing,
                                    output_dir=encounter_csv_directory(feather_path))
    isdir(feather_path) && (feather_path = joinpath(feather_path, "trajectory.feather"))
    table = Arrow.Table(feather_path)
    metadata = Arrow.getmetadata(table)
    stored_range = metadata === nothing ? nothing : get(metadata, "maximum_range_m", nothing)
    if maximum_range_m === nothing
        stored_range === nothing && throw(ArgumentError("Supply maximum_range_m for legacy Feather files"))
        maximum_range_m = parse(Float64, stored_range)
    end
    isfinite(maximum_range_m) && maximum_range_m > 0 ||
        throw(ArgumentError("maximum_range_m must be finite and positive"))
    directory = dirname(feather_path)
    geometry_path = joinpath(directory, "geometry_encounters.feather")
    laser_path = joinpath(directory, "laser_on.feather")
    event_bundle = isfile(geometry_path) && isfile(laser_path)
    if basename(feather_path) == "trajectory.feather" && !event_bundle
        throw(ArgumentError("trajectory.feather requires geometry_encounters.feather and laser_on.feather"))
    end
    geometry = event_bundle ? DataFrame(Arrow.Table(geometry_path)) : nothing
    laser_intervals = event_bundle ? DataFrame(Arrow.Table(laser_path)) : nothing
    if event_bundle
        for path in (geometry_path, laser_path)
            interval_metadata = Arrow.getmetadata(Arrow.Table(path))
            interval_metadata["scenario"] == metadata["scenario"] ||
                throw(ArgumentError("Mismatched scenario in $path"))
            parse(Float64, interval_metadata["maximum_range_m"]) == maximum_range_m ||
                throw(ArgumentError("Range threshold does not match the recorded events"))
        end
        for interval in eachrow(laser_intervals)
            matches = filter(row -> row.encounter_id == interval.encounter_id, geometry)
            nrow(matches) == 1 || throw(ArgumentError("Laser interval has an invalid encounter ID"))
            encounter = matches[1, :]
            (interval.sc_a, interval.sc_b) == (encounter.sc_a, encounter.sc_b) ||
                throw(ArgumentError("Laser interval pair does not match its encounter"))
            encounter.start_time_s <= interval.start_time_s <= interval.end_time_s <= encounter.end_time_s ||
                throw(ArgumentError("Laser interval extends beyond its geometry encounter"))
        end
    end
    times = Float64.(table.time)
    length(times) >= 2 && all(isfinite, times) && all(diff(times) .> 0) ||
        throw(ArgumentError("At least two finite, strictly increasing times are required"))
    spacecraft_ids = sort!([parse(Int, matched.captures[1])
        for name in propertynames(table)
        for matched in (match(r"^sc(\d+)_pos_1$", string(name)),) if matched !== nothing])
    length(spacecraft_ids) >= 2 && first(spacecraft_ids) == 1 ||
        throw(ArgumentError("At least two spacecraft starting at sc1 are required"))
    target_id = metadata === nothing ? 1 : parse(Int, get(metadata, "target_id", "1"))
    target_id in spacecraft_ids || throw(ArgumentError("Invalid target ID in metadata"))
    states = Dict{Int,Matrix{Float64}}()
    for spacecraft in spacecraft_ids
        state = Matrix{Float64}(undef, 6, length(times))
        for (field, offset) in (("pos", 0), ("vel", 3)), component in 1:3
            state[offset+component, :] = getproperty(table, Symbol("sc$(spacecraft)_$(field)_$(component)"))
        end
        all(isfinite, state) || throw(ArgumentError("Nonfinite spacecraft state for sc$spacecraft"))
        states[spacecraft] = state
    end
    active_helpers = hasproperty(table, :laser_active_helper) ? table.laser_active_helper : nothing
    if active_helpers !== nothing
        all(helper -> ismissing(helper) || helper == 0 || (helper in spacecraft_ids && helper != target_id), active_helpers) ||
            throw(ArgumentError("Invalid laser_active_helper spacecraft ID"))
    end

    encounters = DataFrame(
        encounter_id=Int[], sc_a=Int[], sc_b=Int[], start_time_s=Float64[], end_time_s=Float64[],
        geometry_duration_s=Float64[], laser_on_estimate_s=Union{Missing,Float64}[],
        start_clipped=Bool[], end_clipped=Bool[], maximum_range_m=Float64[])
    event_bundle && rename!(encounters, :laser_on_estimate_s => :laser_on_s)
    samples = DataFrame(
        encounter_id=Int[], sc_a=Int[], sc_b=Int[], time_s=Float64[], interpolated=Bool[],
        geometry_range_m=Float64[], relative_speed_mps=Float64[], range_rate_mps=Union{Missing,Float64}[],
        laser_on_held=Union{Missing,Bool}[])
    event_bundle && rename!(samples, :laser_on_held => :laser_on)
    for label in ("a", "b"), field in ("r", "v"), axis in ("x", "y", "z")
        samples[!, "sc_$(label)_$(field)_$(axis)"] = Float64[]
    end
    extrema = DataFrame(
        metric=String[], extremum=String[], value=Union{Missing,Float64}[],
        sc_a=Union{Missing,Int}[], sc_b=Union{Missing,Int}[],
        time_s=Union{Missing,Float64}[], encounter_id=Union{Missing,Int}[])
    for metric in ("relative_speed_mps", "range_rate_mps"), extremum in ("min", "max")
        push!(extrema, (metric, extremum, missing, missing, missing, missing, missing))
    end

    for first_index in 1:length(spacecraft_ids)-1, second_index in first_index+1:length(spacecraft_ids)
        sc_a, sc_b = spacecraft_ids[first_index], spacecraft_ids[second_index]
        pair_geometry = event_bundle ? filter(row -> row.sc_a == sc_a && row.sc_b == sc_b, geometry) : nothing
        event_bundle && isempty(pair_geometry) && continue
        state_a, state_b = states[sc_a], states[sc_b]
        motions = [pair_motion(state_a[:, sample], state_b[:, sample]) for sample in eachindex(times)]
        ranges = first.(motions)
        windows = event_bundle ? [(row.start_time_s, row.end_time_s) for row in eachrow(pair_geometry)] :
            contact_windows(times, ranges, maximum_range_m)
        isempty(windows) && continue

        function laser_on_at(sample)
            active_helpers === nothing && return missing
            helper = active_helpers[sample]
            ismissing(helper) && return missing
            return (sc_a == target_id && helper == sc_b) || (sc_b == target_id && helper == sc_a)
        end

        for (window_index, (start_time, end_time)) in enumerate(windows)
            event = event_bundle ? pair_geometry[window_index, :] : nothing
            encounter_id = event_bundle ? event.encounter_id : nrow(encounters) + 1
            laser_duration = active_helpers === nothing ? missing : 0.0
            linked_laser = event_bundle ? filter(row -> row.encounter_id == encounter_id, laser_intervals) : nothing
            if event_bundle
                laser_duration = sum(linked_laser.duration_s; init=0.0)
            else
              for sample in 1:length(times)-1
                overlap = max(0.0, min(end_time, times[sample+1]) - max(start_time, times[sample]))
                overlap == 0 && continue
                on = laser_on_at(sample)
                if ismissing(on)
                    laser_duration = missing
                    break
                elseif on
                    laser_duration += overlap
                end
              end
            end
            push!(encounters, (encounter_id, sc_a, sc_b, start_time, end_time,
                end_time-start_time, laser_duration,
                event_bundle ? event.start_clipped : start_time == first(times),
                event_bundle ? event.end_clipped : end_time == last(times), Float64(maximum_range_m)))
            sample_times = sort!(unique(vcat(start_time, times[start_time .< times .< end_time], end_time)))
            for time in sample_times
                interpolated = !event_bundle && !(time in times)
                current_a = encounter_state_at(times, state_a, time)
                current_b = encounter_state_at(times, state_b, time)
                if event_bundle && (time == start_time || time == end_time)
                    boundary = time == start_time ? "start" : "end"
                    current_a = [event["$(boundary)_sc_a_$(field)_$(axis)"] for field in ("r", "v") for axis in ("x", "y", "z")]
                    current_b = [event["$(boundary)_sc_b_$(field)_$(axis)"] for field in ("r", "v") for axis in ("x", "y", "z")]
                end
                _, speed, range_rate = pair_motion(current_a, current_b)
                for (metric_index, value) in enumerate((speed, range_rate)), extremum_index in 1:2
                    ismissing(value) && continue
                    row = 2*(metric_index-1) + extremum_index
                    previous = extrema.value[row]
                    if ismissing(previous) || (extremum_index == 1 ? value < previous : value > previous)
                        extrema[row, [:value, :sc_a, :sc_b, :time_s, :encounter_id]] =
                            [value, sc_a, sc_b, time, encounter_id]
                    end
                end
                sample = min(searchsortedlast(times, time), length(times)-1)
                fraction = (time-times[sample]) / (times[sample+1]-times[sample])
                geometry_range = (1-fraction)*ranges[sample] + fraction*ranges[sample+1]
                event_bundle && (geometry_range = norm(current_b[1:3]-current_a[1:3]))
                held_sample = searchsortedlast(times, time)
                on = event_bundle ? any(row -> row.start_time_s <= time < row.end_time_s ||
                    (row.end_clipped && time == row.end_time_s), eachrow(linked_laser)) : laser_on_at(held_sample)
                push!(samples, (encounter_id, sc_a, sc_b, time, interpolated,
                    geometry_range, speed, range_rate, on, current_a..., current_b...))
            end
        end
    end
    laser_metric = event_bundle ? "laser_on_s" : "laser_on_estimate_s"
    for metric in ("geometry_duration_s", laser_metric), extremum in ("min", "max")
        metric_label = metric == "geometry_duration_s" ? "encounter_duration_s" : metric
        valid_rows = findall(!ismissing, encounters[!, metric])
        if isempty(valid_rows)
            push!(extrema, (metric_label, extremum, missing, missing, missing, missing, missing))
        else
            values = encounters[valid_rows, metric]
            selected = extremum == "min" ? argmin(values) : argmax(values)
            encounter = encounters[valid_rows[selected], :]
            push!(extrema, (metric_label, extremum, encounter[metric], encounter.sc_a, encounter.sc_b,
                encounter.start_time_s, encounter.encounter_id))
        end
    end
    analysis_dir = joinpath(output_dir, "analysis")
    mkpath(analysis_dir)
    extrema_path = joinpath(analysis_dir, "extrema.csv")
    encounters_path = joinpath(analysis_dir, "encounters.csv")
    samples_path = joinpath(analysis_dir, "encounter_states.csv")
    CSV.write(extrema_path, extrema)
    sort!(encounters, :encounter_id)
    CSV.write(encounters_path, encounters)
    CSV.write(samples_path, samples)
    encounter_states_dir = joinpath(analysis_dir, "encounter_states")
    mkpath(encounter_states_dir)
    state_columns = [Symbol("sc_$(label)_$(field)_$(axis)") =>
        Symbol("spacecraft_$(label)_$(field)_$(axis)_", field == "r" ? "m" : "m_s")
        for label in ("a", "b") for field in ("r", "v") for axis in ("x", "y", "z")]
    encounter_filenames = Set("encounter_$(encounter_id).csv" for encounter_id in encounters.encounter_id)
    for filename in readdir(encounter_states_dir)
        if occursin(r"^encounter_\d+\.csv$", filename) && filename ∉ encounter_filenames
            rm(joinpath(encounter_states_dir, filename))
        end
    end
    for encounter_samples in groupby(samples, :encounter_id; sort=true)
        encounter_id = first(encounter_samples.encounter_id)
        states_table = select(encounter_samples, :sc_a => :spacecraft_a,
            :sc_b => :spacecraft_b, :time_s, state_columns...)
        states_table.laser_on = map(states_table.time_s) do time
            if event_bundle
                return any(interval -> interval.start_time_s <= time < interval.end_time_s ||
                    (interval.end_clipped && time == interval.end_time_s), eachrow(laser_intervals))
            end
            active_helpers === nothing && return missing
            helper = active_helpers[searchsortedlast(times, time)]
            return ismissing(helper) ? missing : helper != 0
        end
        CSV.write(joinpath(encounter_states_dir, "encounter_$(encounter_id).csv"), states_table)
    end
    longest_encounter = only(filter(row -> row.metric == "encounter_duration_s" &&
        row.extremum == "max", eachrow(extrema)))
    if !ismissing(longest_encounter.encounter_id)
        filename = "encounter_$(longest_encounter.encounter_id).csv"
        cp(joinpath(encounter_states_dir, filename), joinpath(analysis_dir, filename); force=true)
    end
    if event_bundle
        sim_output_dir = joinpath(output_dir, "sim_output")
        mkpath(sim_output_dir)
        CSV.write(joinpath(sim_output_dir, "trajectory.csv"), table)
        CSV.write(joinpath(sim_output_dir, "encounter_record.csv"), geometry)
        CSV.write(joinpath(sim_output_dir, "laser_link_record.csv"), laser_intervals)
    end
    println("Saved $(nrow(encounters)) encounters and $(nrow(samples)) states to $output_dir")
    return (; extrema, encounters, samples, extrema_path, encounters_path, samples_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (1, 2, 3) || error("Usage: extract_feather_encounters.jl SCENARIO_DIRECTORY_OR_FEATHER [MAXIMUM_RANGE_KM [OUTPUT_DIR]]")
    feather_path = abspath(ARGS[1])
    maximum_range_m = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) * 1000 : nothing
    output_dir = length(ARGS) == 3 ? abspath(ARGS[3]) : encounter_csv_directory(feather_path)
    extract_feather_encounters(feather_path; maximum_range_m, output_dir)
end