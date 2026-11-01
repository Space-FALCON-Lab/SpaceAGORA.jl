# The feasible-edge graph of a constellation: which spacecraft pairs could carry
# a laser link at one sample time, under line of sight past the planet, range
# against transmitted power, each end's field of regard, terminal count and duty
# cycle. It is a planning view, evaluated on propagated states for an edge
# schedule over a horizon; a run decides its own links with
# `update_availability!`. Included by `InterLinkModels` (interlink_graph.jl),
# after laser_beam_links.jl, whose occultation and beam geometry it reuses.

"""
    OpticalTerminal(; boresight, field_of_regard_rad, frame=:rtn)

One steerable optical terminal of a spacecraft in an `EdgeFeasibilityModel`. It
can point anywhere within `field_of_regard_rad` (in (0, π]) of `boresight`, a
direction in `frame`:

- `:rtn`: the spacecraft's own radial, along-track and cross-track axes, from
  its position and velocity, for a spacecraft that keeps a local-orbital
  attitude. The convention is that of the RPO frame utilities: radial along the
  position, cross-track along the orbit normal, along-track completing the set.
- `:body`: its body frame, through the attitude of its state sample, which must
  have one. `terminal_pose(spacecraft, terminal)` gives a mounted terminal's
  boresight in that frame.
- `:inertial`: the inertial frame.

A direction exactly `field_of_regard_rad` from the boresight is inside, as for
`TerminalPointing`. A terminal carries at most one active link at a time.
"""
struct OpticalTerminal
    boresight::SVector{3, Float64}
    field_of_regard_rad::Float64
    frame::Symbol

    function OpticalTerminal(boresight, field_of_regard_rad::Real, frame::Symbol)
        axis = SVector{3, Float64}(boresight)
        n = norm(axis)
        isfinite(n) && n > 0 || throw(ArgumentError("An OpticalTerminal's boresight must be finite and nonzero; got $(axis)."))
        isfinite(field_of_regard_rad) && 0 < field_of_regard_rad <= pi ||
            throw(ArgumentError("field_of_regard_rad must be in (0, π]; got $(field_of_regard_rad)."))
        frame in (:rtn, :body, :inertial) ||
            throw(ArgumentError("An OpticalTerminal's frame is :rtn, :body or :inertial; got $(repr(frame))."))
        new(axis / n, Float64(field_of_regard_rad), frame)
    end
end

OpticalTerminal(; boresight, field_of_regard_rad::Real, frame::Symbol=:rtn) =
    OpticalTerminal(boresight, field_of_regard_rad, frame)

"""
    MinimumForceRange(; min_force_n)

A range against power for an `EdgeFeasibilityModel`: the largest separation
`ρ_max(P)` at which its laser link, transmitting `P` watts, still exerts at
least `min_force_n` newtons on each end (`power_limited_range_m`).

- **With a `GaussianBeam`.** The force is `magnification × P / c` times the
  beam's capture fraction at the far aperture (`beam_capture_fraction`), which
  falls with range. With `η = F_min c / (magnification × P)`, the waist radius
  `w0`, the aperture radius `a` and the Rayleigh range `z_R = π w0² / λ`,
  `ρ_max = z_R √(w² / w0² − 1)` for the spot radius `w² = 2a² / (−ln(1 − η))`
  at which the capture fraction is `η`.
- **Limits.** `ρ_max` is zero when even the waist captures less than `η` (or
  `η ≥ 1`), and unbounded when `η ≤ 0`.
- **Without a beam** the force does not fall with range: `ρ_max` is unbounded
  when `magnification × P / c ≥ F_min` and zero otherwise.

The link's own `range_m` caps `ρ_max` in the feasible-edge graph.
`min_force_n` must be finite and nonnegative.
"""
struct MinimumForceRange
    min_force_n::Float64

    function MinimumForceRange(min_force_n::Real)
        isfinite(min_force_n) && min_force_n >= 0 ||
            throw(ArgumentError("min_force_n must be finite and nonnegative; got $(min_force_n)."))
        new(Float64(min_force_n))
    end
end

MinimumForceRange(; min_force_n::Real) = MinimumForceRange(min_force_n)

"""
    power_limited_range_m(limit, link, power_w) -> Float64

The range `ρ_max(P)` (m) up to which `link` can close while transmitting
`power_w` watts, by the range-against-power rule `limit` of an
`EdgeFeasibilityModel`: a `MinimumForceRange`, any function of the power in
watts returning metres, or `nothing` for no power dependence (`Inf`). The
feasible-edge graph uses `min(link_range_m(link), ρ_max(P))`.
"""
power_limited_range_m(::Nothing, link, power_w::Real) = Inf
power_limited_range_m(limit, link, power_w::Real) = Float64(limit(Float64(power_w)))

function power_limited_range_m(limit::MinimumForceRange, link::LaserBeamLink, power_w::Real)
    power = Float64(power_w)
    isfinite(power) && power >= 0 || throw(ArgumentError("power_w must be finite and nonnegative; got $(power_w)."))
    scale = link.magnification * power / SPEED_OF_LIGHT_MPS
    link.beam === nothing && return scale >= limit.min_force_n ? Inf : 0.0
    beam = link.beam
    beam isa GaussianBeam || throw(ArgumentError(
        "MinimumForceRange needs a GaussianBeam or no beam; for a $(nameof(typeof(beam))), give the range " *
        "against power as a function of the power in watts."))
    scale == 0 && return limit.min_force_n == 0 ? Inf : 0.0
    fraction = limit.min_force_n / scale
    fraction <= 0 && return Inf
    fraction >= 1 && return 0.0
    spot_sq = 2.0 * beam.aperture_radius_m^2 / -log1p(-fraction)
    ratio = spot_sq / beam.waist_radius_m^2 - 1.0
    ratio <= 0 && return 0.0
    return (pi * beam.waist_radius_m^2 / beam.wavelength_m) * sqrt(ratio)
end

"""
    DutyCycleLimit(; on_steps, window_steps, scope=:spacecraft)

A duty-cycle limit on laser links over a grid of planning steps: in any
`window_steps` consecutive steps, at most `on_steps` active link-steps

- at each spacecraft (`scope=:spacecraft`): each active link of the spacecraft
  counts once in each step it is active, so two links in one step count two;
- on each link (`scope=:edge`): the steps in which that link is active.

`on_steps` must be a nonnegative integer and `window_steps` a positive one.
"""
struct DutyCycleLimit
    on_steps::Int
    window_steps::Int
    scope::Symbol

    function DutyCycleLimit(on_steps::Integer, window_steps::Integer, scope::Symbol)
        on_steps >= 0 || throw(ArgumentError("on_steps must be nonnegative; got $(on_steps)."))
        window_steps >= 1 || throw(ArgumentError("window_steps must be positive; got $(window_steps)."))
        scope in (:spacecraft, :edge) ||
            throw(ArgumentError("A DutyCycleLimit's scope is :spacecraft or :edge; got $(repr(scope))."))
        new(Int(on_steps), Int(window_steps), scope)
    end
end

DutyCycleLimit(; on_steps::Integer, window_steps::Integer, scope::Symbol=:spacecraft) =
    DutyCycleLimit(on_steps, window_steps, scope)

"""
    EdgeFeasibilityModel(link, terminals; pairs=all_pairs(terminals), power_w=link.power_w,
                         max_range=nothing, duty_cycle=nothing)

The rule that decides which spacecraft pairs of a constellation could carry a
laser link at a sample time (`feasible_edge_graph`). An edge `(i, j)` is
feasible when

1. **line of sight**: with the link's `occultation`, the segment between the
   two spacecraft clears the planet's sphere raised by its
   `min_ray_altitude_m`, as for the link in a run;
2. **range against power**: `0 < ρ_ij ≤ min(link_range_m(link), ρ_max(P_ij))`,
   with `ρ_max` from `max_range` (`power_limited_range_m`) and `P_ij` the pair's
   transmitted power;
3. **field of regard at both ends**: some terminal of `i` sees `ℓ̂_ij`, the unit
   line of sight from `i` to `j`, within its field of regard, and some terminal
   of `j` sees `−ℓ̂_ij`;
4. **terminal and duty-cycle limits**: an edge needs a terminal at each end
   (3); a set of edges active together needs distinct terminals
   (`terminal_assignment`); and with `duty_cycle`, the edge must stay within the
   `DutyCycleLimit` given the steps before (`previous`).

- `link` is a `LaserBeamLink` (a `LaserMomentumLink` is taken as the beam link
  with its defaults) without terminal pointing or mounted recoil: the field of
  regard comes from `terminals`, and the geometry is the spacecraft
  body-reference positions.
- `terminals` holds, for each spacecraft by vector index, its `OpticalTerminal`s.
- `pairs` lists the candidate pairs `(i, j)`; each is kept as `(min, max)`.
- `power_w` is the transmitted power (W) of every pair, a vector with one entry
  per pair, or a function `(i, j) -> P`.
"""
struct EdgeFeasibilityModel{L <: LaserBeamLink, R}
    link::L
    terminals::Vector{Vector{OpticalTerminal}}
    pairs::Vector{Tuple{Int, Int}}
    power_w::Vector{Float64}
    max_range::R
    duty_cycle::Union{Nothing, DutyCycleLimit}
end

_feasibility_link(link::LaserBeamLink) = link
_feasibility_link(link::LaserMomentumLink) = LaserBeamLink(; power_w=link.power_w, magnification=link.magnification,
    range_m=link.range_m)
_feasibility_link(link) = throw(ArgumentError(
    "An EdgeFeasibilityModel takes a LaserBeamLink or a LaserMomentumLink; got $(typeof(link))."))

function EdgeFeasibilityModel(link, terminals; pairs=all_pairs(terminals), power_w=nothing, max_range=nothing,
                              duty_cycle::Union{Nothing, DutyCycleLimit}=nothing)
    beam_link = _feasibility_link(link)
    (beam_link.pointing === nothing && beam_link.recoil === :central) || throw(ArgumentError(
        "An EdgeFeasibilityModel takes its field of regard from its OpticalTerminals; give it a link without " *
        "terminal pointing or mounted recoil."))
    ends = [collect(OpticalTerminal, set) for set in terminals]
    n = length(ends)
    keys = Tuple{Int, Int}[]
    for (i, j) in pairs
        (1 <= i <= n && 1 <= j <= n) || throw(ArgumentError("Pair ($i, $j) names a spacecraft outside 1:$n."))
        i == j && throw(ArgumentError("Pair ($i, $j) joins a spacecraft to itself."))
        key = minmax(Int(i), Int(j))
        key in keys && throw(ArgumentError("Pair $key is listed twice."))
        push!(keys, key)
    end
    powers = if power_w === nothing
        fill(beam_link.power_w, length(keys))
    elseif power_w isa Real
        fill(Float64(power_w), length(keys))
    elseif power_w isa AbstractVector
        length(power_w) == length(keys) || throw(ArgumentError(
            "power_w lists $(length(power_w)) powers for $(length(keys)) pairs."))
        collect(Float64, power_w)
    else
        [Float64(power_w(i, j)) for (i, j) in keys]
    end
    all(p -> isfinite(p) && p >= 0, powers) || throw(ArgumentError("Every pair's power_w must be finite and nonnegative."))
    return EdgeFeasibilityModel(beam_link, ends, keys, powers, max_range, duty_cycle)
end

"""
    FeasibleEdgeGraph

The feasible-edge graph `G_k` at one sample time, from `feasible_edge_graph`.
Each field is a vector over `pairs`, in the model's order, unless noted.

- `time_s`: the sample time (s), from the `LinkContext`.
- `pairs`: the candidate pairs `(i, j)`, `i < j`.
- `feasible`: whether each edge is feasible: all of the flags below.
- `line_of_sight`, `in_range`, `duty_cycle`: each constraint's flag.
- `field_of_regard`: `(at i, at j)`: whether some terminal at that end sees the
  other end within its field of regard.
- `range_m`, `max_range_m`: the separation and its limit `min(link_range_m,
  ρ_max(P))`.
- `direction`: `ℓ̂_ij`, the inertial unit line of sight from `i` to `j`.
- `terminal_options`: `(at i, at j)`: the numbers of the terminals at each end
  that see the other end.
- `terminal_count`: each spacecraft's number of terminals (by vector index).
- `duty_limit`: the model's `DutyCycleLimit`, or `nothing`.
"""
struct FeasibleEdgeGraph
    time_s::Float64
    pairs::Vector{Tuple{Int, Int}}
    feasible::Vector{Bool}
    line_of_sight::Vector{Bool}
    in_range::Vector{Bool}
    field_of_regard::Vector{NTuple{2, Bool}}
    duty_cycle::Vector{Bool}
    range_m::Vector{Float64}
    max_range_m::Vector{Float64}
    direction::Vector{SVector{3, Float64}}
    terminal_options::Vector{NTuple{2, Vector{Int}}}
    terminal_count::Vector{Int}
    duty_limit::Union{Nothing, DutyCycleLimit}
end

"""
    feasible_edge_graph(model::EdgeFeasibilityModel, samples; context=LinkContext(), previous=()) -> FeasibleEdgeGraph
    feasible_edge_graph(model::EdgeFeasibilityModel, state; spacecraft=nothing, context=LinkContext(), previous=())

The feasible-edge graph `G_k` of `model` at one sample time.

- `samples` holds one `StateSample` per spacecraft, by vector index; `state` is
  a run state, whose spacecraft `i` is `state.sc[i]` (as a solution `sol(t)` of
  `run_simulation`), with the run's spacecraft models as `spacecraft` when
  known. A `:body` terminal needs the sample's attitude.
- `context` gives the time and the planet (any value with an equatorial radius
  `Rp_e`); without a planet, nothing blocks the line of sight.
- `previous` lists the active edge sets of the steps just before this one,
  oldest first, each a collection of pairs. With a `DutyCycleLimit`, an edge is
  feasible only if activating it at this step, alone, keeps every window that
  ends here within the limit. Edges active together at this step share the
  budget: `admissible_edge_set` checks a whole set.
"""
function feasible_edge_graph(model::EdgeFeasibilityModel, samples::AbstractVector{<:StateSample};
                             context::LinkContext=LinkContext(), previous=())
    n = length(model.terminals)
    length(samples) == n || throw(ArgumentError(
        "The feasible-edge graph needs one state sample per spacecraft: $n terminal sets, $(length(samples)) samples."))
    link = model.link
    history = [_pair_set(set) for set in previous]
    count = length(model.pairs)
    graph = FeasibleEdgeGraph(Float64(context.time_s), copy(model.pairs), fill(false, count), fill(false, count),
        fill(false, count), fill((false, false), count), fill(true, count), zeros(count), zeros(count),
        fill(SVector(0.0, 0.0, 0.0), count), [(Int[], Int[]) for _ in 1:count], length.(model.terminals),
        model.duty_cycle)
    for (e, (i, j)) in enumerate(model.pairs)
        a, b = samples[i], samples[j]
        separation = b.pos_ii - a.pos_ii
        distance = norm(separation)
        sight = distance > 0 ? separation / distance : separation
        limit = min(link_range_m(link), power_limited_range_m(model.max_range, link, model.power_w[e]))
        clear = !(link.occultation && _segment_occulted(a.pos_ii, b.pos_ii, context.planet, link.min_ray_altitude_m))
        options = (_seeing_terminals(model.terminals[i], a, sight), _seeing_terminals(model.terminals[j], b, -sight))
        duty = model.duty_cycle === nothing ||
            _duty_window_ok(model.duty_cycle, vcat(history, [Set([(i, j)])]), length(history) + 1)
        graph.line_of_sight[e] = clear
        graph.in_range[e] = 0 < distance <= limit
        graph.field_of_regard[e] = (!isempty(options[1]), !isempty(options[2]))
        graph.duty_cycle[e] = duty
        graph.range_m[e] = distance
        graph.max_range_m[e] = limit
        graph.direction[e] = sight
        graph.terminal_options[e] = options
        graph.feasible[e] = clear && graph.in_range[e] && all(graph.field_of_regard[e]) && duty
    end
    return graph
end

function feasible_edge_graph(model::EdgeFeasibilityModel, state; spacecraft=nothing, context::LinkContext=LinkContext(),
                             previous=())
    samples = [link_endpoint(state.sc[i], _spacecraft_at(spacecraft, i), i) for i in eachindex(model.terminals)]
    return feasible_edge_graph(model, samples; context, previous)
end

"""
    feasible_edge_graphs(model::EdgeFeasibilityModel, states, times; planet=nothing, spacecraft=nothing) -> Vector{FeasibleEdgeGraph}

The feasible-edge graphs `G_k` of a propagated constellation at the sample
times `times` (s): `states(t)` returns the constellation's run state at time
`t`, as the solution of `run_simulation(...; return_solution=true)` does.
`planet` is the run's planet (for the line of sight), and `spacecraft` its
spacecraft models when known. No activity precedes the samples, so each graph's
duty-cycle flags are all true; the scheduler applies the limit to its schedule.
"""
function feasible_edge_graphs(model::EdgeFeasibilityModel, states, times; planet=nothing, spacecraft=nothing)
    return [feasible_edge_graph(model, states(t); spacecraft, context=LinkContext(time_s=t, planet=planet))
            for t in times]
end

"""
    feasible_pairs(graph::FeasibleEdgeGraph) -> Vector{Tuple{Int, Int}}

The pairs of `graph` that are feasible, in the graph's order.
"""
feasible_pairs(graph::FeasibleEdgeGraph) = graph.pairs[graph.feasible]

"""
    terminal_assignment(graph::FeasibleEdgeGraph, active) -> Union{Nothing, Dict{Tuple{Int, Int}, NTuple{2, Int}}}

A terminal for each end of every pair in `active` such that no terminal
carries two links, each terminal seeing the other end within its field of
regard (`graph.terminal_options`): a dictionary from each pair to its terminal
numbers `(at i, at j)`. `nothing` when no such assignment exists, as when a
spacecraft has fewer terminals able to point than active links. Each
spacecraft's terminals are matched to its active links by augmenting paths, so
an assignment is found whenever one exists.
"""
function terminal_assignment(graph::FeasibleEdgeGraph, active)
    pairs = collect(_pair_set(active))
    index = Dict(pair => e for (e, pair) in enumerate(graph.pairs))
    for pair in pairs
        haskey(index, pair) || throw(ArgumentError("Pair $pair is not a candidate of the feasible-edge graph."))
    end
    assignment = Dict(pair => (0, 0) for pair in pairs)
    for craft in sort!(unique!([s for pair in pairs for s in pair]))
        incident = [(pair, pair[1] == craft ? 1 : 2) for pair in pairs if craft in pair]
        options = [graph.terminal_options[index[pair]][side] for (pair, side) in incident]
        matched = _match_terminals(options, graph.terminal_count[craft])
        matched === nothing && return nothing
        for ((pair, side), terminal) in zip(incident, matched)
            assignment[pair] = side == 1 ? (terminal, assignment[pair][2]) : (assignment[pair][1], terminal)
        end
    end
    return assignment
end

"""
    admissible_edge_set(graph::FeasibleEdgeGraph, active; max_active=nothing, previous=()) -> Bool

Whether the edges `active` (a collection of pairs) can be active together at
the graph's step: each is a candidate whose line of sight, range and field of
regard hold, there are at most `max_active` of them, a `terminal_assignment`
exists, and with the graph's `DutyCycleLimit`, the step's edges, after the
steps `previous` (oldest first), keep every window that ends here within it.
"""
function admissible_edge_set(graph::FeasibleEdgeGraph, active; max_active::Union{Nothing, Integer}=nothing,
                             previous=())
    pairs = _pair_set(active)
    max_active === nothing || length(pairs) <= max_active || return false
    index = Dict(pair => e for (e, pair) in enumerate(graph.pairs))
    for pair in pairs
        e = get(index, pair, 0)
        e == 0 && throw(ArgumentError("Pair $pair is not a candidate of the feasible-edge graph."))
        graph.line_of_sight[e] && graph.in_range[e] && all(graph.field_of_regard[e]) || return false
    end
    terminal_assignment(graph, pairs) === nothing && return false
    graph.duty_limit === nothing && return true
    history = vcat([_pair_set(set) for set in previous], [pairs])
    return _duty_window_ok(graph.duty_limit, history, length(history))
end

# The active edges `set` as a set of `(min, max)` pairs.
_pair_set(set) = Set{Tuple{Int, Int}}(minmax(Int(i), Int(j)) for (i, j) in set)

# The terminal numbers in `terminals` that see the unit inertial direction
# `sight` from `sample` within their field of regard.
function _seeing_terminals(terminals::Vector{OpticalTerminal}, sample, sight::SVector{3, Float64})
    seeing = Int[]
    iszero(sight) && return seeing
    for (number, terminal) in enumerate(terminals)
        local_sight = _terminal_frame_sight(terminal.frame, sample, sight)
        angle = atan(norm(cross(local_sight, terminal.boresight)), dot(local_sight, terminal.boresight))
        angle <= terminal.field_of_regard_rad && push!(seeing, number)
    end
    return seeing
end

_terminal_frame_sight(frame::Symbol, sample, sight) =
    frame === :inertial ? sight : frame === :body ? _body_sight(sample, sight) : _rtn_sight(sample, sight)

function _body_sight(sample, sight)
    sample.q_ib === nothing && throw(ArgumentError(
        "A :body OpticalTerminal needs the attitude of spacecraft $(sample.spacecraft_index)'s state sample."))
    return rot(sample.q_ib) * sight
end

# The inertial direction `sight` in the radial, along-track and cross-track
# axes of `sample`'s own orbit.
function _rtn_sight(sample, sight)
    r, v = sample.pos_ii, sample.vel_ii
    h = cross(r, v)
    (norm(r) > 0 && isfinite(norm(h)) && norm(h) > 0) || throw(ArgumentError(
        "An :rtn OpticalTerminal needs a nonzero position and orbital angular momentum for spacecraft " *
        "$(sample.spacecraft_index)."))
    radial = r / norm(r)
    normal = h / norm(h)
    along = normalize(cross(normal, radial))
    return SVector(dot(radial, sight), dot(along, sight), dot(normal, sight))
end

# A distinct terminal for each link from its `options`, by augmenting paths, or
# `nothing` when there is none.
function _match_terminals(options::Vector{Vector{Int}}, terminal_count::Int)
    owner = zeros(Int, terminal_count)
    for link in eachindex(options)
        _augment!(owner, options, link, falses(terminal_count)) || return nothing
    end
    matched = zeros(Int, length(options))
    for (terminal, link) in enumerate(owner)
        link == 0 || (matched[link] = terminal)
    end
    return matched
end

function _augment!(owner::Vector{Int}, options, link::Int, seen::BitVector)
    for terminal in options[link]
        seen[terminal] && continue
        seen[terminal] = true
        if owner[terminal] == 0 || _augment!(owner, options, owner[terminal], seen)
            owner[terminal] = link
            return true
        end
    end
    return false
end

# Whether the active edge sets `steps` (one per step, oldest first) keep every
# `limit.window_steps` window that contains step `k` within `limit`. Windows that
# run past the listed steps count only the listed ones.
function _duty_window_ok(limit::DutyCycleLimit, steps::AbstractVector, k::Int)
    for start in max(1, k - limit.window_steps + 1):k
        counts = Dict{Tuple{Int, Int}, Int}()
        for step in start:min(start + limit.window_steps - 1, length(steps)), pair in steps[step]
            for key in _duty_keys(limit, pair)
                counts[key] = get(counts, key, 0) + 1
                counts[key] > limit.on_steps && return false
            end
        end
    end
    return true
end

# What one active link-step counts against: each end spacecraft `(i, 0)`, or the link itself.
_duty_keys(limit::DutyCycleLimit, (i, j)) = limit.scope === :edge ? ((i, j),) : ((i, 0), (j, 0))
