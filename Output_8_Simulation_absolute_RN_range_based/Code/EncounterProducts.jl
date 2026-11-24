module EncounterProducts

using DataFrames, LinearAlgebra, StaticArrays

export write_dense_sample, build_encounter_products!, longest_ids, boundary_state,
    encounter_quantities

function write_dense_sample(io, time, state)
    write(io, Float64(time))
    write(io, state)
    nothing
end

function boundary_state(event, boundary, endpoint)
    SVector{6,Float64}(ntuple(k -> getproperty(event,
        Symbol("$(boundary)_sc_$(endpoint)_$(k <= 3 ? "r" : "v")_$("xyz"[mod1(k, 3)])")), 6))
end

function encounter_quantities(a, b)
    r = SVector{3,Float64}(b[1:3] - a[1:3])
    v = SVector{3,Float64}(b[4:6] - a[4:6])
    distance = norm(r)
    distance > 0 || error("Coincident spacecraft have undefined range rate")
    (norm(v), dot(r, v) / distance, distance)
end

function longest_ids(geometry, target)
    isempty(geometry) && return (0, 0)
    ordered = sort(geometry, [:duration_s, :start_time_s, :encounter_id]; rev=[true, false, false])
    overall = first(ordered.encounter_id)
    subset = filter(row -> row.sc_a == target || row.sc_b == target, ordered)
    (overall, isempty(subset) ? 0 : first(subset.encounter_id))
end

function update_extrema!(geometry, index, time, a, b)
    values = encounter_quantities(a, b)
    for (quantity, value) in zip(("relative_speed_mps", "range_rate_mps", "separation_m"), values)
        for (kind, compare) in (("min", <), ("max", >))
            field = Symbol("$(kind)_$quantity")
            if compare(value, geometry[index, field])
                geometry[index, field] = value
                geometry[index, Symbol("$(kind)_$(quantity)_time_s")] = time
            end
        end
    end
    geometry.sample_count[index] += 1
end

function track_table()
    table = DataFrame(run_id=String[], encounter_id=Int[], sc_a=Int[], sc_b=Int[],
        time=Float64[], is_all_pair_longest=Bool[], is_target_helper_longest=Bool[])
    for endpoint in ("a", "b"), field in ("pos", "vel"), axis in 1:3
        table[!, Symbol("sc_$(endpoint)_$(field)_$axis")] = Float64[]
    end
    table
end

function append_track!(tracks, event, time, a, b, selected, run_id)
    push!(tracks, (run_id, event.encounter_id, event.sc_a, event.sc_b, time,
        event.encounter_id == selected[1], event.encounter_id == selected[2], a..., b...))
end

function build_encounter_products!(geometry, spool_path, spacecraft_count, run_id, duration; cadence=10.0)
    for quantity in ("relative_speed_mps", "range_rate_mps", "separation_m")
        for kind in ("min", "max")
            geometry[!, Symbol("$(kind)_$quantity")] = fill(kind == "min" ? Inf : -Inf, nrow(geometry))
            geometry[!, Symbol("$(kind)_$(quantity)_time_s")] = zeros(nrow(geometry))
        end
    end
    geometry.sample_count = zeros(Int, nrow(geometry))
    selected = longest_ids(geometry, spacecraft_count)
    tracks = track_table()
    for (index, event) in enumerate(eachrow(geometry)), boundary in ("start", "end")
        a, b = boundary_state(event, boundary, "a"), boundary_state(event, boundary, "b")
        time = getproperty(event, Symbol("$(boundary)_time_s"))
        update_extrema!(geometry, index, time, a, b)
        event.encounter_id in selected && append_track!(tracks, event, time, a, b, selected, run_id)
    end
    if isempty(geometry)
        return tracks
    end
    order = sortperm(geometry.start_time_s)
    active = Set{Int}()
    cursor, rows = 1, 0
    state = Matrix{Float64}(undef, 6, spacecraft_count)
    last_time = -Inf
    open(spool_path, "r") do io
        while !eof(io)
            time = read(io, Float64)
            read!(io, state)
            time > last_time || error("Nonmonotonic dense encounter spool")
            isapprox(time, min(rows * cadence, duration); atol=2eps(duration), rtol=0) ||
                error("Dense encounter spool cadence mismatch")
            while cursor <= length(order) && geometry.start_time_s[order[cursor]] < time
                push!(active, order[cursor])
                cursor += 1
            end
            filter!(index -> geometry.end_time_s[index] > time, active)
            for index in active
                event = geometry[index, :]
                a, b = SVector{6,Float64}(state[:, event.sc_a]), SVector{6,Float64}(state[:, event.sc_b])
                update_extrema!(geometry, index, time, a, b)
                event.encounter_id in selected && append_track!(tracks, event, time, a, b, selected, run_id)
            end
            rows += 1
            last_time = time
        end
    end
    last_time == duration || error("Dense encounter spool lacks its final sample")
    sort!(tracks, [:encounter_id, :time])
    tracks
end

end
