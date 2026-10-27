# Generated from SpaceAGORA.jl scripts/companion_support/mission_history.jl.
# Source SHA-256: 67c6110e0471e0588cd613f582c43e326812c2b7fda482361f2c6c21a5f1c687; regenerate with scripts/sync_companion_support.jl.
# Mission history is optional plain data shared with the independent viewer.
# Validate structure without changing times, names or the underlying results.
function _validated_mission_history(mission)
    mission === nothing && return nothing
    mission isa AbstractDict || throw(ArgumentError("mission must be a JSON object."))
    schema = get(mission, "schema", nothing)
    schema isa Integer && !(schema isa Bool) && schema == 1 ||
        throw(ArgumentError("Unsupported mission history schema; expected 1."))
    name_ok(x) = x isa AbstractString && !isempty(strip(x))
    time_ok(x) = x isa Real && !(x isa Bool) && isfinite(x) && x >= 0
    name_ok(get(mission, "name", nothing)) || throw(ArgumentError("mission.name must be a nonempty string."))
    segments = get(mission, "segments", nothing)
    transitions = get(mission, "transitions", nothing)
    segments isa AbstractVector || throw(ArgumentError("mission.segments must be an array."))
    transitions isa AbstractVector || throw(ArgumentError("mission.transitions must be an array."))
    names = Set{String}()
    for (i, segment) in enumerate(segments)
        segment isa AbstractDict || throw(ArgumentError("mission segment must be an object."))
        name = get(segment, "name", nothing)
        name_ok(name) && !(name in names) || throw(ArgumentError("Mission segment names must be nonempty and unique."))
        push!(names, name)
        start_s, end_s = get(segment, "start_s", nothing), get(segment, "end_s", nothing)
        time_ok(start_s) && time_ok(end_s) && start_s <= end_s ||
            throw(ArgumentError("Mission segment times must be finite and satisfy 0 <= start_s <= end_s."))
        i == 1 || segments[i - 1]["end_s"] == start_s ||
            throw(ArgumentError("Mission segment intervals must be contiguous and ordered."))
    end
    n = length(segments)
    (n == 0 ? isempty(transitions) : length(transitions) in (n - 1, n)) ||
        throw(ArgumentError("Mission history needs one transition per adjacent segment and an optional final exit."))
    for (i, transition) in enumerate(transitions)
        transition isa AbstractDict || throw(ArgumentError("mission transition must be an object."))
        t = get(transition, "time_s", nothing)
        time_ok(t) && t == segments[i]["end_s"] ||
            throw(ArgumentError("Mission transition time must match its segment boundary."))
        get(transition, "from", nothing) == segments[i]["name"] ||
            throw(ArgumentError("Mission transition.from must name its outgoing segment."))
        haskey(transition, "to") && transition["to"] == (i < n ? segments[i + 1]["name"] : nothing) ||
            throw(ArgumentError("Mission transition.to must name the next segment, or be null at final exit."))
        name_ok(get(transition, "reason", nothing)) ||
            throw(ArgumentError("Mission transition.reason must be a nonempty string."))
    end
    return mission
end
