module FeatherRecording

using Arrow, DataFrames, DiffEqCallbacks, LinearAlgebra, Printf
using DiffEqBase: CallbackSet, DiscreteCallback, VectorContinuousCallback

export output_times, scenario_paths, IntervalRecorder, recording_callbacks, finish_recording!, write_intervals

function output_times(duration, interval=10.0)
    duration > 0 || throw(ArgumentError("duration must be positive"))
    isfinite(interval) && interval > 0 || throw(ArgumentError("output interval must be finite and positive"))
    times = collect(0.0:Float64(interval):Float64(duration))
    last(times) < duration && push!(times, Float64(duration))
    return times
end

function scenario_paths(root, opts, duration; source, smoke=false, schedule=opts.schedule, use_J2=true)
    interval = hasproperty(opts, :output_interval_s) ? opts.output_interval_s : 10.0
    name = @sprintf("%s_N%d_h%.6gkm_t%.6gkm_ih%.6g_it%.6g_e%.6g_nu%.6g_T%.9gs_R%.6gkm_P%.6gW_B%.6g_m%.6gkg_%s_J2%s_dt%.6gs",
        source, opts.helpers, opts.helper_altitude_km, opts.target_altitude_km,
        opts.helper_inclination_deg, opts.target_inclination_deg, opts.target_ecc,
        opts.target_nu_deg, duration, opts.laser_range_km, opts.laser_power_w,
        opts.magnification, opts.mass_kg, string(schedule), string(use_J2), interval)
    if hasproperty(opts, :beta) && (opts.beta != 1.0 || opts.eta != 1.0)
        name *= @sprintf("_beta%.6g_eta%.6g", opts.beta, opts.eta)
    end
    base = smoke ? joinpath(root, "smoke") : root
    return (feather=joinpath(base, "feather", name), csv=joinpath(base, "CSV", name), scenario=name)
end

function interval_table(; laser=false)
    table = DataFrame(encounter_id=Int[], sc_a=Int[], sc_b=Int[],
        start_time_s=Float64[], end_time_s=Float64[], duration_s=Float64[],
        start_clipped=Bool[], end_clipped=Bool[])
    laser && insertcols!(table, 1, :laser_id => Int[])
    for boundary in ("start", "end"), label in ("a", "b"), field in ("r", "v"), axis in ("x", "y", "z")
        table[!, "$(boundary)_sc_$(label)_$(field)_$(axis)"] = Float64[]
    end
    return table
end

mutable struct IntervalRecorder{State, Laser, Gates}
    state::State
    laser::Laser
    gates::Gates
    maximum_range::Float64
    pairs::Vector{Tuple{Int,Int}}
    geometry::DataFrame
    laser_on::DataFrame
    open_geometry::Dict{Tuple{Int,Int},NamedTuple}
    open_laser::Dict{Tuple{Int,Int},NamedTuple}
    next_geometry::Int
    next_laser::Int
    initialized::Bool
    gate_count::Int
end

function IntervalRecorder(state, laser, initial_state, maximum_range;
                          gates=(states, time)->Float64[])
    isfinite(maximum_range) && maximum_range > 0 || throw(ArgumentError("A positive finite laser range is required"))
    count = size(state(initial_state), 2)
    pairs = [(first, second) for first in 1:count-1 for second in first+1:count]
    return IntervalRecorder(state, laser, gates, Float64(maximum_range), pairs,
        interval_table(), interval_table(; laser=true),
        Dict{Tuple{Int,Int},NamedTuple}(), Dict{Tuple{Int,Int},NamedTuple}(),
        0, 0, false, length(gates(state(initial_state), 0.0)))
end

function close_interval!(table, open, pair, time, states; clipped=false)
    entry = pop!(open, pair)
    identifiers = hasproperty(entry, :laser_id) ? (entry.laser_id, entry.encounter) : (entry.encounter,)
    push!(table, (identifiers..., pair..., entry.time, Float64(time),
        Float64(time)-entry.time, entry.clipped, clipped,
        entry.states[:, pair[1]]..., entry.states[:, pair[2]]...,
        states[:, pair[1]]..., states[:, pair[2]]...))
end

function observe!(recorder, states, time; initial=false, geometry_pair=nothing, entering=false)
    if initial
        recorder.initialized = true
        for pair in recorder.pairs
            if norm(states[1:3, pair[2]] - states[1:3, pair[1]]) <= recorder.maximum_range
                recorder.next_geometry += 1
                recorder.open_geometry[pair] = (encounter=recorder.next_geometry,
                    time=Float64(time), states=copy(states), clipped=true)
            end
        end
    end
    if geometry_pair !== nothing && entering && !haskey(recorder.open_geometry, geometry_pair)
        recorder.next_geometry += 1
        recorder.open_geometry[geometry_pair] = (encounter=recorder.next_geometry,
            time=Float64(time), states=copy(states), clipped=false)
    end
    probe = copy(states)
    probe[1:3, :] .+= 1e-6 .* states[4:6, :]
    active = Set{Tuple{Int,Int}}(recorder.laser(probe, time+1e-6))
    if geometry_pair !== nothing && !entering
        delete!(active, geometry_pair)
    end
    for pair in collect(keys(recorder.open_laser))
        if !(pair in active)
            close_interval!(recorder.laser_on, recorder.open_laser, pair, time, states)
        end
    end
    for pair in sort!(collect(active))
        haskey(recorder.open_geometry, pair) || continue
        if !haskey(recorder.open_laser, pair)
            recorder.next_laser += 1
            recorder.open_laser[pair] = (laser_id=recorder.next_laser,
                encounter=recorder.open_geometry[pair].encounter, time=Float64(time),
                states=copy(states), clipped=initial)
        end
    end
    if geometry_pair !== nothing && !entering && haskey(recorder.open_geometry, geometry_pair)
        close_interval!(recorder.geometry, recorder.open_geometry, geometry_pair, time, states)
    end
    return nothing
end

function recording_callbacks(recorder)
    function condition!(values, state, time, integrator)
        states = recorder.state(state)
        for (index, pair) in enumerate(recorder.pairs)
            values[index] = norm(states[1:3, pair[2]] - states[1:3, pair[1]]) - recorder.maximum_range
        end
        values[length(recorder.pairs)+1:end] .= recorder.gates(states, time)
    end
    function crossing!(integrator, index, entering)
        states = recorder.state(integrator.u)
        transitions = Tuple{Tuple{Int,Int},Bool}[]
        for (pair_index, pair) in enumerate(recorder.pairs)
            if pair_index == index
                push!(transitions, (pair, entering))
                continue
            end
            relative_position = states[1:3, pair[2]] - states[1:3, pair[1]]
            if abs(norm(relative_position) - recorder.maximum_range) <= 1e-6
                relative_velocity = states[4:6, pair[2]] - states[4:6, pair[1]]
                direction = dot(relative_position, relative_velocity)
                direction != 0 && push!(transitions, (pair, direction < 0))
            end
        end
        sort!(transitions; by=transition -> !transition[2])
        for (pair, pair_entering) in transitions
            observe!(recorder, states, integrator.t; geometry_pair=pair, entering=pair_entering)
        end
        isempty(transitions) && observe!(recorder, states, integrator.t)
    end
    geometry_cb = VectorContinuousCallback(condition!,
        (integrator, index)->crossing!(integrator, index, false),
        (integrator, index)->crossing!(integrator, index, true),
        length(recorder.pairs)+recorder.gate_count;
        abstol=1e-8, reltol=0.0, interp_points=20, save_positions=(false, false))
    initialize = (callback, state, time, integrator) -> begin
        recorder.initialized || observe!(recorder, recorder.state(state), time; initial=true)
    end
    observer_cb = DiscreteCallback((state, time, integrator)->true,
        integrator->observe!(recorder, recorder.state(integrator.u), integrator.t);
        initialize=initialize, save_positions=(false, false))
    return geometry_cb, observer_cb
end

function finish_recording!(recorder, state, time)
    states = recorder.state(state)
    for pair in sort!(collect(keys(recorder.open_laser)))
        close_interval!(recorder.laser_on, recorder.open_laser, pair, time, states; clipped=true)
    end
    for pair in sort!(collect(keys(recorder.open_geometry)))
        close_interval!(recorder.geometry, recorder.open_geometry, pair, time, states; clipped=true)
    end
    sort!(recorder.geometry, :encounter_id)
    sort!(recorder.laser_on, :laser_id)
    return recorder
end

function write_intervals(recorder, directory; source, scenario, laser_timing)
    mkpath(directory)
    metadata = Dict("source"=>source, "scenario"=>scenario,
        "maximum_range_m"=>string(recorder.maximum_range),
        "geometry_timing"=>"continuous range roots; abstol=1e-8; interp_points=20",
        "simultaneous_boundary_tolerance_m"=>"1e-6",
        "laser_timing"=>laser_timing,
        "boundary_side_probe_s"=>"1e-6", "time_units"=>"s", "state_frame"=>"ECI")
    Arrow.write(joinpath(directory, "geometry_encounters.feather"), recorder.geometry; metadata)
    Arrow.write(joinpath(directory, "laser_on.feather"), recorder.laser_on; metadata)
end

end