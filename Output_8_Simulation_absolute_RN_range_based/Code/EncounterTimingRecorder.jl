module EncounterTimingRecorder

using LinearAlgebra
using StaticArrays
using DataFrames

export Recorder, snapshot, observe_geometry!, observe_firing!, finish!,
    hermite_state, crossing_brackets, bisect_crossing, analyze_intervals, range_segments

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

function close_interval!(table, open, pair, time, states; clipped=false)
    entry = pop!(open, pair)
    time >= entry.time || error("Interval ends before it starts")
    time == entry.time && return nothing
    identifiers = hasproperty(entry, :laser_id) ? (entry.laser_id, entry.encounter) : (entry.encounter,)
    push!(table, (identifiers..., pair..., entry.time, Float64(time),
        Float64(time) - entry.time, entry.clipped, clipped,
        entry.states[:, pair[1]]..., entry.states[:, pair[2]]...,
        states[:, pair[1]]..., states[:, pair[2]]...))
end

snapshot(u) = hcat([vcat(SVector{3, Float64}(s.pos), SVector{3, Float64}(s.vel))
    for s in u.sc]...)

function hermite_state(left, right, duration, fraction)
    s = fraction
    position = (2s^3 - 3s^2 + 1) * left[1:3, :] +
        (s^3 - 2s^2 + s) * duration * left[4:6, :] +
        (-2s^3 + 3s^2) * right[1:3, :] + (s^3 - s^2) * duration * right[4:6, :]
    velocity = ((6s^2 - 6s) * left[1:3, :] +
        (-6s^2 + 6s) * right[1:3, :]) / duration +
        (3s^2 - 4s + 1) * left[4:6, :] + (3s^2 - 2s) * right[4:6, :]
    return vcat(position, velocity)
end

function stationary_fractions(coefficients)
    n = length(coefficients)
    scale = maximum(abs, coefficients)
    scale == 0 && return Float64[]
    while n > 1 && abs(coefficients[n]) <= 16eps(scale)
        n -= 1
    end
    n == 1 && return Float64[]
    degree = n - 1
    companion = zeros(degree, degree)
    for i in 2:degree
        companion[i, i - 1] = 1.0
    end
    companion[:, end] .= -coefficients[1:degree] / coefficients[n]
    return sort!(unique([real(z) for z in eigvals(companion)
        if abs(imag(z)) < 1e-9 && 0 < real(z) < 1]))
end

function crossing_brackets(r0, v0, r1, v1, duration, maximum_range; inside=norm(r0) <= maximum_range)
    controls = (r0, r0 + duration * v0 / 3, r1 - duration * v1 / 3, r1)
    lower_bound_sq = sum(begin
        low, high = extrema(p[k] for p in controls)
        low > 0 ? low^2 : high < 0 ? high^2 : 0.0
    end for k in 1:3)
    # Bezier bounds avoid solving polynomials for distant or wholly interior pairs.
    lower_bound_sq > (maximum_range + 0.01)^2 && return Tuple{Float64,Float64,Bool}[]
    all(norm(p) < maximum_range - 0.01 for p in controls) &&
        return Tuple{Float64,Float64,Bool}[]
    c = (r0, duration * v0, 3(r1 - r0) - duration * (2v0 + v1),
        2(r0 - r1) + duration * (v0 + v1))
    squared = zeros(7)
    for i in 1:4, j in 1:4
        squared[i + j - 1] += dot(c[i], c[j])
    end
    critical = stationary_fractions([i * squared[i + 1] for i in 1:6])
    knots = vcat(0.0, critical, 1.0)
    relative(s) = c[1] + s * (c[2] + s * (c[3] + s * c[4]))
    brackets = Tuple{Float64,Float64,Bool}[]
    for i in 2:length(knots)
        left, right = knots[i - 1], knots[i]
        gap = norm(relative(right)) - maximum_range
        after = if abs(gap) < 1e-7
            if i < length(knots)
                norm(relative((right + knots[i + 1]) / 2)) < maximum_range
            else
                dot(r1, v1) < 0
            end
        else
            gap < 0
        end
        if after != inside
            push!(brackets, (left, right, after))
            inside = after
        end
    end
    return brackets
end

function bisect_crossing(gap, left, right; tolerance_s=1e-9)
    fl, fr = gap(left), gap(right)
    abs(fl) <= 1e-7 && return left
    abs(fr) <= 1e-7 && return right
    signbit(fl) != signbit(fr) ||
        error("Range crossing is not bracketed: [$left, $right], gaps=($fl, $fr)")
    for _ in 1:80
        right - left <= tolerance_s && return (left + right) / 2
        middle = (left + right) / 2
        fm = gap(middle)
        if signbit(fm) == signbit(fl)
            left, fl = middle, fm
        else
            right = middle
        end
    end
    error("Range-root bisection did not converge")
end

mutable struct Recorder
    maximum_range::Float64
    pairs::Vector{Tuple{Int,Int}}
    geometry::DataFrame
    laser::DataFrame
    open_geometry::Dict{Tuple{Int,Int},NamedTuple}
    open_laser::Dict{Tuple{Int,Int},NamedTuple}
    geometry_id::Int
    laser_id::Int
    max_boundary_residual_m::Float64
end

function Recorder(initial, maximum_range)
    isfinite(maximum_range) && maximum_range > 0 ||
        throw(ArgumentError("Encounter range must be finite and positive"))
    count = size(initial, 2)
    laser = interval_table(; laser=true)
    rename!(laser, :encounter_id => :start_encounter_id)
    recorder = Recorder(maximum_range,
        [(a, b) for a in 1:count-1 for b in a+1:count],
        interval_table(), laser, Dict{Tuple{Int,Int},NamedTuple}(),
        Dict{Tuple{Int,Int},NamedTuple}(), 0, 0, 0.0)
    for pair in recorder.pairs
        if norm(initial[1:3, pair[2]] - initial[1:3, pair[1]]) <= maximum_range
            recorder.geometry_id += 1
            recorder.open_geometry[pair] = (encounter=recorder.geometry_id,
                time=0.0, states=copy(initial), clipped=true)
        end
    end
    return recorder
end

@inline function may_cross(r0, v0, r1, v1, dt, range)
    c1 = r0 + (dt / 3) * v0
    c2 = r1 - (dt / 3) * v1
    lower = 0.0
    @inbounds for k in 1:3
        low = min(r0[k], c1[k], c2[k], r1[k])
        high = max(r0[k], c1[k], c2[k], r1[k])
        lower += low > 0 ? low^2 : high < 0 ? high^2 : 0.0
    end
    lower > (range + 0.01)^2 && return false
    inner = max(range - 0.01, 0.0)^2
    return !(dot(r0, r0) < inner && dot(c1, c1) < inner &&
        dot(c2, c2) < inner && dot(r1, r1) < inner)
end

function observe_geometry!(recorder, previous, current, start, stop, dense_state)
    dt = stop - start
    dt > 0 || error("Recorder needs strictly increasing step times")
    for pair in recorder.pairs
        a, b = pair
        r0 = SVector{3, Float64}(ntuple(k -> previous[k, b] - previous[k, a], 3))
        r1 = SVector{3, Float64}(ntuple(k -> current[k, b] - current[k, a], 3))
        v0 = SVector{3, Float64}(ntuple(k -> previous[k + 3, b] - previous[k + 3, a], 3))
        v1 = SVector{3, Float64}(ntuple(k -> current[k + 3, b] - current[k + 3, a], 3))
        may_cross(r0, v0, r1, v1, dt, recorder.maximum_range) || continue
        brackets = crossing_brackets(r0, v0, r1, v1, dt, recorder.maximum_range;
            inside=haskey(recorder.open_geometry, pair))
        isempty(brackets) && continue
        gap(t) = begin
            state = dense_state(t)
            norm(state[1:3, b] - state[1:3, a]) - recorder.maximum_range
        end
        for (left, right, entering) in brackets
            time = bisect_crossing(gap, start + left * dt, start + right * dt)
            state = dense_state(time)
            recorder.max_boundary_residual_m = max(recorder.max_boundary_residual_m, abs(gap(time)))
            if entering
                haskey(recorder.open_geometry, pair) && error("Duplicate encounter entry for $pair")
                recorder.geometry_id += 1
                recorder.open_geometry[pair] = (encounter=recorder.geometry_id,
                    time=time, states=state, clipped=false)
            else
                haskey(recorder.open_geometry, pair) || error("Encounter exit without entry for $pair")
                close_interval!(recorder.geometry, recorder.open_geometry, pair, time, state)
            end
        end
    end
end

function observe_firing!(recorder, active, state, time; initial=false)
    selected = Set(active)
    for pair in collect(keys(recorder.open_laser))
        pair in selected ||
            close_interval!(recorder.laser, recorder.open_laser, pair, time, state)
    end
    for pair in sort!(collect(selected))
        haskey(recorder.open_laser, pair) && continue
        recorder.laser_id += 1
        encounter = get(recorder.open_geometry, pair, nothing)
        recorder.open_laser[pair] = (laser_id=recorder.laser_id,
            encounter=encounter === nothing ? 0 : encounter.encounter,
            time=Float64(time), states=copy(state), clipped=initial)
    end
end

function finish!(recorder, state, time)
    for pair in sort!(collect(keys(recorder.open_laser)))
        close_interval!(recorder.laser, recorder.open_laser, pair, time, state; clipped=true)
    end
    for pair in sort!(collect(keys(recorder.open_geometry)))
        close_interval!(recorder.geometry, recorder.open_geometry, pair, time, state; clipped=true)
    end
    sort!(recorder.geometry, [:sc_a, :sc_b, :start_time_s])
    sort!(recorder.laser, :start_time_s)
    return recorder
end

overlap(a, b, c, d) = max(0.0, min(b, d) - max(a, c))

function analyze_intervals(geometry, laser)
    encounters = select(geometry, :encounter_id, :sc_a, :sc_b, :start_time_s,
        :end_time_s, :duration_s, :start_clipped, :end_clipped)
    encounters.laser_on_s = Float64[sum((overlap(e.start_time_s, e.end_time_s,
        l.start_time_s, l.end_time_s) for l in eachrow(laser)
        if (l.sc_a, l.sc_b) == (e.sc_a, e.sc_b)); init=0.0) for e in eachrow(encounters)]
    audit = DataFrame(laser_id=Int[], sc_a=Int[], sc_b=Int[], actual_firing_s=Float64[],
        within_range_s=Float64[], outside_range_s=Float64[])
    for l in eachrow(laser)
        within = sum((overlap(e.start_time_s, e.end_time_s, l.start_time_s, l.end_time_s)
            for e in eachrow(geometry) if (e.sc_a, e.sc_b) == (l.sc_a, l.sc_b)); init=0.0)
        push!(audit, (l.laser_id, l.sc_a, l.sc_b, l.duration_s, within, l.duration_s - within))
    end
    return encounters, audit
end

function range_segments(geometry, laser; run_id::String)
    isempty(run_id) && throw(ArgumentError("A nonempty run ID is required"))
    for (table, id) in ((geometry, :encounter_id), (laser, :laser_id))
        length(unique(table[!, id])) == nrow(table) ||
            throw(ArgumentError("Duplicate interval IDs"))
        for row in eachrow(table)
            isfinite(row.start_time_s) && isfinite(row.end_time_s) &&
                0 <= row.start_time_s < row.end_time_s ||
                throw(ArgumentError("Invalid interval times"))
            row.sc_a < row.sc_b || throw(ArgumentError("Pairs must be ordered"))
        end
        for group in groupby(table, [:sc_a, :sc_b])
            ordered = sort(group, :start_time_s)
            all(ordered.end_time_s[1:end-1] .<= ordered.start_time_s[2:end]) ||
                throw(ArgumentError("Overlapping intervals for the same pair"))
        end
    end
    segments = DataFrame(run_id=String[], laser_id=Int[], sc_a=Int[], sc_b=Int[],
        encounter_id=Union{Missing,Int}[], start_time_s=Float64[], end_time_s=Float64[],
        duration_s=Float64[], within_range=Bool[])
    for firing in eachrow(laser)
        contacts = filter(e -> (e.sc_a, e.sc_b) == (firing.sc_a, firing.sc_b) &&
            e.start_time_s < firing.end_time_s && e.end_time_s > firing.start_time_s, geometry)
        boundaries = [firing.start_time_s, firing.end_time_s]
        for e in eachrow(contacts)
            push!(boundaries, max(firing.start_time_s, e.start_time_s),
                min(firing.end_time_s, e.end_time_s))
        end
        sort!(unique!(boundaries))
        for i in 1:length(boundaries)-1
            start, stop = boundaries[i], boundaries[i + 1]
            found = findfirst(e -> e.start_time_s <= start && stop <= e.end_time_s,
                eachrow(contacts))
            encounter = found === nothing ? missing : contacts.encounter_id[found]
            push!(segments, (run_id, firing.laser_id, firing.sc_a, firing.sc_b,
                encounter, start, stop, stop - start, found !== nothing))
        end
    end
    return segments
end

end
