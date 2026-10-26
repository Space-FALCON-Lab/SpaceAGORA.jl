module InterlinkTerminals

# Interlink terminal hardware on a spacecraft: which mounted actuators are
# terminals and what they are numbered, a run's hardware link endpoints, the
# levels and power their links set, and the shedding that closes their links.
# The link contract and the graph are in InterLinkModels (interlink_graph.jl),
# which comes before the mounted-actuator and power models this module needs;
# its terminal methods for a SpacecraftModel are added here.

using StaticArrays: SVector
using ..AbstractTypes: AbstractActuator
using ..SpacecraftModels: SpacecraftModel
using ..InterLinkModels: InterlinkModel, LinkContext, link_terminal, link_terminal_levels, terminal_supports,
    link_directed, link_received_power_w, link_endpoint, PowerBeamLink, terminal_shed_tier, link_data_rate_bps,
    _endpoint_sample, _node_active, _station_antenna
import ..InterLinkModels: terminal_count, _terminal_device, _terminal_devices, _terminal_number, _battery_capacity_wh,
    _battery_within_full, terminal_pose, _data_store_capacity_bits, _data_store_within_empty,
    _data_store_within_full, _data_store_production_bps, _data_store_outflows_bps
using ..ActuatorMounting: mounted_actuators, actuator_context
using ..Actuators: actuator_consumption, actuator_limits
using ..PowerSystems: battery, within_full_tolerance
using ..DataStorage: has_data_store, data_store, _production_bps, _within_store_empty, _within_store_full

export mounted_link_terminals, terminal_number, data_link_rates_bps

"""
    mounted_link_terminals(spacecraft::SpacecraftModel) -> Vector

The terminal hardware mounted on `spacecraft` (actuators with `link_terminal`),
in mounting order. Each entry gives:

- `number`: its terminal number, after the `n_terminal` virtual terminals;
- `position`: its position among the mounted actuators, the index of its level;
- `link`: the link it is mounted on;
- `terminal`: the actuator.
"""
function mounted_link_terminals(spacecraft::SpacecraftModel)
    out = NamedTuple{(:number, :position, :link, :terminal), Tuple{Int, Int, Any, Any}}[]
    for (position, entry) in enumerate(mounted_actuators(spacecraft))
        link_terminal(entry.actuator) || continue
        push!(out, (number=spacecraft.n_terminal + length(out) + 1, position=position,
                    link=entry.link, terminal=entry.actuator))
    end
    return out
end

terminal_count(spacecraft::SpacecraftModel) =
    spacecraft.n_terminal + count(entry -> link_terminal(entry.actuator), mounted_actuators(spacecraft))

_terminal_device(spacecraft::SpacecraftModel, number::Integer) =
    number <= spacecraft.n_terminal ? nothing : mounted_link_terminals(spacecraft)[number - spacecraft.n_terminal].terminal
_terminal_devices(spacecraft::SpacecraftModel) =
    Any[fill(nothing, spacecraft.n_terminal); [entry.terminal for entry in mounted_link_terminals(spacecraft)]]

"""
    terminal_number(spacecraft::SpacecraftModel, terminal) -> Int

The terminal number of `terminal`, mounted on `spacecraft`, for a link endpoint
`(spacecraft index, terminal number)`. `register_link!` also takes the terminal
itself. Two mounted terminals with equal fields cannot be told apart, so for
those give the number.
"""
function terminal_number(spacecraft::SpacecraftModel, terminal)::Int
    matches = [entry.number for entry in mounted_link_terminals(spacecraft) if entry.terminal === terminal]
    isempty(matches) && throw(ArgumentError(
        "The $(nameof(typeof(terminal))) given as a link endpoint is not terminal hardware mounted on spacecraft " *
        "$(spacecraft.id) (see link_terminal and add_actuator!)."))
    length(matches) == 1 || throw(ArgumentError(
        "Spacecraft $(spacecraft.id) has $(length(matches)) mounted terminals equal to the one given (numbers " *
        "$(join(matches, ", "))): give the terminal number."))
    return only(matches)
end

_terminal_number(spacecraft::SpacecraftModel, terminal::AbstractActuator) = terminal_number(spacecraft, terminal)

function terminal_pose(spacecraft::SpacecraftModel, terminal)
    pose = nothing
    for link in spacecraft.links, mounted in Iterators.flatten((link.laser_terminals, link.custom_actuators))
        mounted === terminal || continue
        here = _mounted_terminal_pose(link, terminal)
        pose === nothing || pose == here || throw(ArgumentError(
            "Spacecraft $(spacecraft.id) has equal $(nameof(typeof(terminal)))s mounted at different poses, so " *
            "terminal_pose cannot tell them apart; give them distinct fields."))
        pose = here
    end
    pose === nothing && throw(ArgumentError(
        "The $(nameof(typeof(terminal))) given is not terminal hardware mounted on spacecraft $(spacecraft.id)."))
    return pose
end

function _mounted_terminal_pose(link, terminal)
    ctx = actuator_context(link)
    axis = hasproperty(terminal, :boresight) ? terminal.boresight :
        hasproperty(terminal, :normal) ? terminal.normal : throw(ArgumentError(
            "terminal_pose: a $(nameof(typeof(terminal))) has no boresight or normal to point with."))
    location = hasproperty(terminal, :location) ? SVector{3, Float64}(terminal.location) : SVector{3, Float64}(0.0, 0.0, 0.0)
    return (ctx.link_to_body * SVector{3, Float64}(axis), ctx.link_to_body * location + ctx.link_r_body)
end

function _battery_capacity_wh(spacecraft::SpacecraftModel)
    b = battery(spacecraft)
    return b === nothing ? NaN : Float64(b.capacity_wh)
end

function _data_store_capacity_bits(spacecraft::SpacecraftModel)
    store = data_store(spacecraft)
    return store === nothing ? NaN : store.capacity_bits
end

_data_store_production_bps(spacecraft::SpacecraftModel, sample, context, t::Float64) =
    _production_bps(spacecraft, sample, context, t)

function _data_store_within_empty(spacecraft::SpacecraftModel, bits::Real)
    store = data_store(spacecraft)
    return store !== nothing && _within_store_empty(store, bits)
end

function _data_store_within_full(spacecraft::SpacecraftModel, bits::Real)
    store = data_store(spacecraft)
    return store !== nothing && _within_store_full(store, bits)
end

function _battery_within_full(spacecraft::SpacecraftModel, energy_wh::Real)
    b = battery(spacecraft)
    return b !== nothing && within_full_tolerance(b, energy_wh)
end

# ---- A run's hardware link endpoints -----------------------------------------

"""
    prepare_terminal_endpoints!(model::InterlinkModel, spacecraft) -> Bool

At run setup, record on `model` the hardware endpoints of every candidate link
of the run's spacecraft. For each spacecraft, each endpoint is
`(edge, side, mounted actuator position)`; for each edge, the terminal hardware
at its two ends (`nothing` for a virtual terminal), which the link hooks
receive. Each hardware terminal must still be
able to carry its link (`terminal_supports`), at a level within its
`actuator_limits`; so must a ground station's antenna, which is not recorded (a
station has no power budget). Returns whether there are any.
"""
function prepare_terminal_endpoints!(model::InterlinkModel, spacecraft)::Bool
    graph = model.graph
    endpoints = [NTuple{3, Int}[] for _ in eachindex(spacecraft)]
    terminals = [Any[nothing, nothing] for _ in eachindex(graph.edges)]
    hardware = Dict{Int, Any}()
    for (edge, vertices) in enumerate(graph.edges), (side, vertex) in enumerate(vertices)
        satellite, number = graph.vertices[vertex]
        if satellite < 0
            _check_station_antenna(model, edge, side, (satellite, number))
            terminals[edge][side] = _station_antenna(model, (satellite, number))
            continue
        end
        model_s = spacecraft[satellite]
        model_s isa SpacecraftModel && number > model_s.n_terminal || continue
        mounted = get!(() -> mounted_link_terminals(model_s), hardware, satellite)
        number - model_s.n_terminal <= length(mounted) || throw(ArgumentError(
            "Interlinks: spacecraft $satellite has no terminal $number (terminal_count $(terminal_count(model_s)))."))
        entry = mounted[number - model_s.n_terminal]
        link = graph.edge_properties[edge]
        role = link_directed(link) ? (side == 1 ? :transmitter : :receiver) : :either
        terminal_supports(entry.terminal, link, role) || throw(ArgumentError(
            "Interlinks: terminal $number of spacecraft $satellite, a $(nameof(typeof(entry.terminal))), cannot " *
            "carry a $(nameof(typeof(link)))$(role === :either ? "" : " as its $(role)") (terminal_supports)."))
        tier = terminal_shed_tier(entry.terminal)
        tier isa Integer && 0 <= tier <= 3 || throw(ArgumentError(
            "Interlinks: terminal $number of spacecraft $satellite, a $(nameof(typeof(entry.terminal))), has shedding " *
            "tier $(repr(tier)); terminal_shed_tier must be 0, 1, 2 or 3."))
        level = Float64(link_terminal_levels(link)[side])
        low, high = actuator_limits(entry.terminal)
        low <= level <= high || throw(ArgumentError(
            "Interlinks: the $(nameof(typeof(link))) on terminal $number of spacecraft $satellite sets it to level " *
            "$level, outside its limits [$low, $high] (for a LaserTerminal, max_power_w); build the link from the " *
            "terminal, LaserMomentumLink(terminal), or lower its power."))
        push!(endpoints[satellite], (edge, side, entry.position))
        terminals[edge][side] = entry.terminal
    end
    model.terminal_endpoints = endpoints
    model.edge_terminals = Tuple{Any, Any}[(pair[1], pair[2]) for pair in terminals]
    return any(!isempty, endpoints)
end

# A ground station's antenna at end `side` of `edge`: able to carry its link, at
# a level within its limits.
function _check_station_antenna(model::InterlinkModel, edge::Int, side::Int, endpoint)
    s, number = -endpoint[1], endpoint[2]
    s <= length(model.ground_stations) && number <= length(model.ground_stations[s].antennas) || throw(ArgumentError(
        "Interlinks: ground station $s has no antenna $number."))
    antenna = _station_antenna(model, endpoint)
    link = model.graph.edge_properties[edge]
    role = link_directed(link) ? (side == 1 ? :transmitter : :receiver) : :either
    terminal_supports(antenna, link, role) || throw(ArgumentError(
        "Interlinks: antenna $number of ground station $s, a $(nameof(typeof(antenna))), cannot carry a " *
        "$(nameof(typeof(link)))$(role === :either ? "" : " as its $(role)") (terminal_supports)."))
    level = Float64(link_terminal_levels(link)[side])
    low, high = actuator_limits(antenna)
    low <= level <= high || throw(ArgumentError(
        "Interlinks: the $(nameof(typeof(link))) on antenna $number of ground station $s sets it to level $level, " *
        "outside its limits [$low, $high]."))
    return nothing
end

# Whether links drive terminal hardware on spacecraft `satellite`.
@inline link_driven(model, satellite::Int) =
    model !== nothing && satellite <= length(model.terminal_endpoints) && !isempty(model.terminal_endpoints[satellite])

"""
    without_link_driven(levels, model, satellite) -> AbstractVector

`levels` (one per mounted actuator of spacecraft `satellite`) with its
link-driven terminals at zero. Their link drives them for the whole run, so a
controller's levels for them apply nothing: no recoil and no draw. This is a
copy, or `levels` itself when no terminal of that spacecraft is link-driven.
"""
function without_link_driven(levels, model, satellite::Int)
    (levels === nothing || !link_driven(model, satellite)) && return levels
    out = collect(Float64, levels)
    for (_, _, position) in model.terminal_endpoints[satellite]
        out[position] = 0.0
    end
    return out
end

"""
    link_driven_levels!(levels, model, satellite) -> levels

Write into `levels` (one per mounted actuator) the levels that the links set on
the link-driven terminals of spacecraft `satellite`. That is the
`link_terminal_levels` of a terminal's active link, and 0 when none of its
links is active. Other entries are left as they are.
"""
function link_driven_levels!(levels, model, satellite::Int)
    link_driven(model, satellite) || return levels
    graph = model.graph
    for (_, _, position) in model.terminal_endpoints[satellite]
        levels[position] = 0.0
    end
    for (edge, side, position) in model.terminal_endpoints[satellite]
        graph.active[edge] || continue
        levels[position] = Float64(link_terminal_levels(graph.edge_properties[edge])[side])
    end
    return levels
end

"""
    link_terminal_power(model, satellite, spacecraft::SpacecraftModel) -> (by_tier, load_w)

The power (W) the link-driven terminals of spacecraft `satellite` (model
`spacecraft`) ask for and draw.

- `by_tier` is the demand, split as `actuator_power_by_tier` does: laser
  terminals, torque rods, reaction wheels, every other kind. Each terminal's
  demand counts in its `terminal_shed_tier` (a tier-0 terminal with every other
  kind). A terminal with an active link asks for its link's level. A terminal
  the battery's shedding has closed keeps asking for the highest level of its
  links, so the tier clears only when the battery can supply it again, not
  because its link closed.
- `load_w` is what is drawn: the active links' levels, less what the shedding
  sheds (tier 1 sheds laser and radio terminals).
"""
function link_terminal_power(model, satellite::Int, spacecraft::SpacecraftModel)
    link_driven(model, satellite) || return ((0.0, 0.0, 0.0, 0.0), 0.0)
    graph = model.graph
    endpoints = model.terminal_endpoints[satellite]
    b = battery(spacecraft)
    tier = b === nothing ? 0 : b.shed_tier
    entries = mounted_actuators(spacecraft)
    first_tier = second_tier = third_tier = never = load = 0.0
    for position in unique(last(endpoint) for endpoint in endpoints)
        active = false
        level = wanted = 0.0
        for (edge, side, at) in endpoints
            at == position || continue
            this = Float64(link_terminal_levels(graph.edge_properties[edge])[side])
            wanted = max(wanted, this)
            graph.active[edge] && (active = true; level = this)
        end
        entry = entries[position]
        own = terminal_shed_tier(entry.actuator)
        shed = own > 0 && tier >= own
        demand = active ? level : (shed ? wanted : 0.0)
        demand == 0.0 && continue
        w = Float64(actuator_consumption(entry.actuator, demand, actuator_context(entry.link)).power_w)
        own == 1 ? (first_tier += w) : own == 2 ? (second_tier += w) : own == 3 ? (third_tier += w) : (never += w)
        active && !shed && (load += w)
    end
    return ((first_tier, second_tier, third_tier, never), load)
end

# ---- Received and beamed power ------------------------------------------------

# Whether spacecraft `satellite` is the receiving end of a directed link on its
# terminal hardware.
@inline link_receives(model, satellite::Int) = link_driven(model, satellite) &&
    any(endpoint -> endpoint[2] == 2 && link_directed(model.graph.edge_properties[endpoint[1]]),
        model.terminal_endpoints[satellite])

"""
    received_link_power_w(model, satellite, spacecraft, state; context=LinkContext()) -> Float64

The electrical power (W) that the active directed links ending on spacecraft
`satellite`'s terminal hardware deliver to its bus, each from
`link_received_power_w` at `state`. `spacecraft` is the run's spacecraft and
`state` the run state, with every spacecraft. A ground station transmitter
reads `planet_frame`, the planet-fixed rotation at the state's time (`nothing`
for the planet's current `L_PI`). Returns 0 when it receives none.
"""
function received_link_power_w(model, satellite::Int, spacecraft, state; context=LinkContext(),
    planet_frame=nothing)::Float64
    link_receives(model, satellite) || return 0.0
    graph = model.graph
    entries = mounted_actuators(spacecraft[satellite])
    total = 0.0
    for (edge, side, position) in model.terminal_endpoints[satellite]
        side == 2 && graph.active[edge] || continue
        link = graph.edge_properties[edge]
        link_directed(link) || continue
        first, second = map(vertex -> graph.vertices[vertex], graph.edges[edge])
        a = first[1] > 0 ? link_endpoint(state.sc[first[1]], spacecraft[first[1]], first[1]) :
            _endpoint_sample(model, first[1], state, spacecraft, context, planet_frame)
        b = link_endpoint(state.sc[second[1]], spacecraft[second[1]], second[1])
        total += Float64(link_received_power_w(link, entries[position].actuator, a, b, context))
    end
    return total
end

"""
    beamed_link_power_w(model, satellite) -> Float64

The optical power (W) that spacecraft `satellite` transmits on its active power
beams (`PowerBeamLink`). Returns 0 when it transmits none.
"""
function beamed_link_power_w(model, satellite::Int)::Float64
    link_driven(model, satellite) || return 0.0
    graph = model.graph
    total = 0.0
    for (edge, side, _) in model.terminal_endpoints[satellite]
        side == 1 && graph.active[edge] || continue
        link = graph.edge_properties[edge]
        link isa PowerBeamLink && (total += link.power_w)
    end
    return total
end

# ---- Data rates ------------------------------------------------------------------

# The terminal hardware at end `side` of `edge` on spacecraft `satellite` (model
# `spacecraft`), from the run's hardware endpoints, or `nothing` for a virtual
# terminal.
function _endpoint_terminal(model, satellite::Int, edge::Int, side::Int, spacecraft)
    link_driven(model, satellite) || return nothing
    for (at_edge, at_side, position) in model.terminal_endpoints[satellite]
        at_edge == edge && at_side == side && return mounted_actuators(spacecraft)[position].actuator
    end
    return nothing
end

"""
    data_link_rates_bps(model, spacecraft, state; is_active=nothing, context=LinkContext()) -> (sent, received)

Per spacecraft, the data rates (bits per second) that the run's active links
carry out of it (`sent`) and into it (`received`), each link's
`link_data_rate_bps` at `state`, the run state with every spacecraft.
`spacecraft` is the run's spacecraft. Optional `is_active` flags drop links with
an inactive end. Undirected links carry none. A link with a ground station end
counts at its spacecraft end, and the station end reads `planet_frame`, the
planet-fixed rotation at the state's time (`nothing` for the planet's current
`L_PI`).
"""
function data_link_rates_bps(model, spacecraft, state; is_active=nothing, context=LinkContext(), planet_frame=nothing)
    rates = _data_link_rates(model, spacecraft, state; is_active, context, planet_frame)
    return (rates.sent, rates.received)
end

# The rates of `data_link_rates_bps`, with each ground station's too:
# `station_sent` and `station_received`, by station number.
function _data_link_rates(model, spacecraft, state; is_active=nothing, context=LinkContext(), planet_frame=nothing)
    sent = zeros(Float64, length(spacecraft))
    received = zeros(Float64, length(spacecraft))
    station_sent = zeros(Float64, model === nothing ? 0 : length(model.ground_stations))
    station_received = zeros(Float64, length(station_sent))
    model === nothing && return (; sent, received, station_sent, station_received)
    graph = model.graph
    for edge in eachindex(graph.edges)
        graph.active[edge] || continue
        link = graph.edge_properties[edge]
        link_directed(link) || continue
        first, second = map(vertex -> graph.vertices[vertex], graph.edges[edge])
        is_active === nothing || (_node_active(is_active, first[1]) && _node_active(is_active, second[1])) || continue
        if first[1] > 0 && second[1] > 0
            transmitter = _endpoint_terminal(model, first[1], edge, 1, spacecraft[first[1]])
            receiver = _endpoint_terminal(model, second[1], edge, 2, spacecraft[second[1]])
            a = link_endpoint(state.sc[first[1]], spacecraft[first[1]], first[1])
            b = link_endpoint(state.sc[second[1]], spacecraft[second[1]], second[1])
            rate = Float64(link_data_rate_bps(link, transmitter, receiver, a, b, context))
            sent[first[1]] += rate
            received[second[1]] += rate
        else
            rate = _station_link_rate_bps(model, edge, first, second, spacecraft, state, context, planet_frame)
            first[1] > 0 ? (sent[first[1]] += rate) : (station_sent[-first[1]] += rate)
            second[1] > 0 ? (received[second[1]] += rate) : (station_received[-second[1]] += rate)
        end
    end
    return (; sent, received, station_sent, station_received)
end

# The rate of the active data link of `edge` between a spacecraft and a ground station.
function _station_link_rate_bps(model, edge::Int, first, second, spacecraft, state, context, planet_frame)
    link = model.graph.edge_properties[edge]
    terminal(endpoint, side) = endpoint[1] < 0 ? _station_antenna(model, endpoint) :
        _endpoint_terminal(model, endpoint[1], edge, side, spacecraft[endpoint[1]])
    a = _endpoint_sample(model, first[1], state, spacecraft, context, planet_frame)
    b = _endpoint_sample(model, second[1], state, spacecraft, context, planet_frame)
    return Float64(link_data_rate_bps(link, terminal(first, 1), terminal(second, 2), a, b, context))
end

# Per spacecraft, the data flows at `state` in a run with data stores (bit/s):
# `sent` and `received`, what its active links carry out and in; `produced`, its
# sources' production; `stored`, the rate of the data its store holds; and
# `lost`, the production a full store cannot keep. A link carries at most its
# `link_data_rate_bps` (`data_link_rates_bps`), and less at a store: an empty
# store sends no more than comes in (its production and its incoming links),
# shared among its outgoing links in proportion to their rates, and a full store
# takes in no more than it sends, so links lose nothing and only its production
# beyond that is lost. A spacecraft without a store neither limits nor is
# limited, and neither does a ground station, an unbounded store: its links'
# rates are `station_sent` and `station_received`, by station number, and its
# end reads `planet_frame` (see `data_link_rates_bps`). An inactive spacecraft's
# links carry nothing and it produces nothing.
function _data_flows(model, spacecraft, state; is_active=nothing, context=LinkContext(), planet_frame=nothing,
                     available_only::Bool=false, active_edges=nothing)
    n = length(spacecraft)
    sent, received, produced, stored, lost = (zeros(Float64, n) for _ in 1:5)
    station_sent = zeros(Float64, model === nothing ? 0 : length(model.ground_stations))
    station_received = zeros(Float64, length(station_sent))
    # The flow nodes: spacecraft `i` is node `i`, and station `s` node `n + s`.
    transmitters, receivers, capacities = Int[], Int[], Float64[]
    if model !== nothing
        graph = model.graph
        for edge in eachindex(graph.edges)
            (active_edges === nothing ? graph.active[edge] : active_edges[edge]) || continue
            available_only && !graph.available[edge] && continue
            link = graph.edge_properties[edge]
            link_directed(link) || continue
            first, second = map(vertex -> graph.vertices[vertex], graph.edges[edge])
            is_active === nothing || (_node_active(is_active, first[1]) && _node_active(is_active, second[1])) || continue
            if first[1] > 0 && second[1] > 0
                transmitter = _endpoint_terminal(model, first[1], edge, 1, spacecraft[first[1]])
                receiver = _endpoint_terminal(model, second[1], edge, 2, spacecraft[second[1]])
                a = link_endpoint(state.sc[first[1]], spacecraft[first[1]], first[1])
                b = link_endpoint(state.sc[second[1]], spacecraft[second[1]], second[1])
                push!(transmitters, first[1])
                push!(receivers, second[1])
                push!(capacities, Float64(link_data_rate_bps(link, transmitter, receiver, a, b, context)))
            else
                push!(transmitters, first[1] > 0 ? first[1] : n - first[1])
                push!(receivers, second[1] > 0 ? second[1] : n - second[1])
                push!(capacities, _station_link_rate_bps(model, edge, first, second, spacecraft, state, context, planet_frame))
            end
        end
    end
    for i in 1:n
        has_data_store(spacecraft[i]) || continue
        is_active === nothing || is_active[i] || continue
        x = link_endpoint(state.sc[i], spacecraft[i], i)
        produced[i] = _production_bps(spacecraft[i], x, context, context.time_s)
    end
    rates = _carried_data_rates(transmitters, receivers, capacities, produced, spacecraft, length(station_sent))
    for e in eachindex(rates)
        i, j = transmitters[e], receivers[e]
        i <= n ? (sent[i] += rates[e]) : (station_sent[i - n] += rates[e])
        j <= n ? (received[j] += rates[e]) : (station_received[j - n] += rates[e])
    end
    for i in 1:n
        store = data_store(spacecraft[i])
        store === nothing && continue
        is_active === nothing || is_active[i] || continue
        kept = store.full ? min(produced[i], max(0.0, sent[i] - received[i])) : produced[i]
        lost[i] = produced[i] - kept
        net = kept + received[i] - sent[i]
        # An empty store does not go below nothing, nor a full one above its
        # capacity, by the rounding of the shares.
        stored[i] = store.empty ? max(net, 0.0) : store.full ? min(net, 0.0) : net
    end
    return (sent=sent, received=received, produced=produced, stored=stored, lost=lost,
        station_sent=station_sent, station_received=station_received)
end

# The rates (bit/s) the active data links carry, from link `e`'s transmitter
# `transmitters[e]`, receiver `receivers[e]` and rate `capacities[e]`, and each
# spacecraft's production: each link at most its rate; out of an empty store, a
# share of what comes in; into a full store, a share of what goes out. Starting
# from no flow, each pass sets every link's rate from the previous pass's
# inflows and outflows, so the rates only grow and every pass keeps within the
# stores' limits. The passes stop when nothing changes, or after one more pass
# than there are stores, which a chain of relays needs. Nodes past the
# spacecraft are the `stations` ground stations, which never limit a link.
function _carried_data_rates(transmitters, receivers, capacities, produced, spacecraft, stations::Int=0)
    limits_out = [(store = data_store(sc); store !== nothing && store.empty) for sc in spacecraft]
    limits_in = [(store = data_store(sc); store !== nothing && store.full) for sc in spacecraft]
    (any(limits_out) || any(limits_in)) || return copy(capacities)
    if stations > 0
        append!(limits_out, fill(false, stations))
        append!(limits_in, fill(false, stations))
        produced = [produced; zeros(stations)]
    end
    n = length(spacecraft) + stations
    out_capacity, in_capacity, inflow, outflow = (zeros(Float64, n) for _ in 1:4)
    for e in eachindex(capacities)
        out_capacity[transmitters[e]] += capacities[e]
        in_capacity[receivers[e]] += capacities[e]
    end
    rates = zeros(Float64, length(capacities))
    next = similar(rates)
    for _ in 1:(count(has_data_store, spacecraft) + 1)
        copyto!(inflow, produced)
        fill!(outflow, 0.0)
        for e in eachindex(rates)
            inflow[receivers[e]] += rates[e]
            outflow[transmitters[e]] += rates[e]
        end
        for e in eachindex(rates)
            i, j = transmitters[e], receivers[e]
            rate = capacities[e]
            limits_out[i] && out_capacity[i] > inflow[i] && (rate = min(rate, capacities[e] * (inflow[i] / out_capacity[i])))
            limits_in[j] && in_capacity[j] > outflow[j] && (rate = min(rate, capacities[e] * (outflow[j] / in_capacity[j])))
            next[e] = rate
        end
        next == rates && break
        rates, next = next, rates
    end
    return rates
end

"""
    shed_link_driven_levels!(levels, model, satellite, spacecraft::SpacecraftModel) -> levels

Set to zero, in `levels` (one per mounted actuator), the link-driven terminals
of spacecraft `satellite` that its battery's shedding sheds (see
`terminal_shed_tier`). Other entries are left as they are.
"""
function shed_link_driven_levels!(levels, model, satellite::Int, spacecraft::SpacecraftModel)
    link_driven(model, satellite) || return levels
    entries = mounted_actuators(spacecraft)
    for (_, _, position) in model.terminal_endpoints[satellite]
        _terminal_shed(spacecraft, entries[position].actuator) && (levels[position] = 0.0)
    end
    return levels
end

"""
    terminal_availability(model, spacecraft) -> Union{Nothing, Vector{Bool}}

One flag per graph vertex, false for terminal hardware that its spacecraft's
power budget sheds: a terminal while the battery sheds its
`terminal_shed_tier` or more (tier 1 for laser and radio terminals).
Links through such a terminal are unavailable until the tier clears. Virtual
terminals and ground station antennas are always available. Returns `nothing` when the run has no hardware
endpoints.
"""
function terminal_availability(model, spacecraft)
    any(!isempty, model.terminal_endpoints) || return nothing
    graph = model.graph
    flags = fill(true, length(graph.vertices))
    for (vertex, (satellite, number)) in enumerate(graph.vertices)
        satellite > 0 || continue # a ground station's antennas are never shed
        model_s = spacecraft[satellite]
        model_s isa SpacecraftModel && number > model_s.n_terminal || continue
        flags[vertex] = !_terminal_shed(model_s, _terminal_device(model_s, number))
    end
    return flags
end

function _terminal_shed(spacecraft::SpacecraftModel, terminal)
    own = terminal_shed_tier(terminal)
    own > 0 || return false
    b = battery(spacecraft)
    return b !== nothing && b.shed_tier >= own
end



# Snapshot current carried outflow before a new matching replaces graph.active.
# Availability has already been refreshed at this scheduling state.
_data_store_outflows_bps(model::InterlinkModel, spacecraft, state, context; is_active=nothing, active_edges=nothing) =
    _data_flows(model, spacecraft, state; context, is_active, active_edges, available_only=true).sent

end # module InterlinkTerminals
