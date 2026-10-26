module SchedulingPolicies

using LinearAlgebra
using ..InterLinkModels
using ..InterLinkModels: link_endpoint, _terminal_devices, _spacecraft_at, _battery_capacity_wh, _battery_within_full,
    _endpoint_sample, _station_antenna, _data_store_capacity_bits, _data_store_within_empty,
    _data_store_within_full, _data_store_production_bps, _data_store_outflows_bps
using ..EffectorSampling: StateSample

export AbstractSchedulingPolicy, LinkCandidate, SchedulingContext, candidate_score
export SemimajorAxisRatePolicy, EccentricityRatePolicy, PowerDeliveryPolicy, DataThroughputPolicy
export EnergyBalancingPolicy, WeightedPolicy, PolicyTerm, state_of_charge
export DataBacklogPolicy, stored_bits, store_room_bits
export SchedulingPolicyModel, score_candidates!
export select_interlinks!, schedule_interlinks!

# 1. Matching workspace for dynamic programming
# 1.1. Scheduler-owned scores and matching scratch (reused between scheduling calls)
struct InterLinkMatchingWorkspace
	graph::InterLinkGraph
	score::Vector{Float64}
	candidates::Vector{Int}
	exact_scores::Vector{Rational{BigInt}}
	weights::Vector{BigInt}
	parent::Vector{Int}
	order::Vector{Int}
	stack::Vector{Int}
	free_score::Vector{BigInt}
	blocked_score::Vector{BigInt}
	choice::Vector{Int}
	blocked::Vector{Bool}
	cyclic::Vector{Bool}
	# A weighted policy's term scores, one row per top-level term and one column per
	# edge, when the interlink model records them (`record_terms`); otherwise `nothing`.
	terms::Union{Nothing, Matrix{Float64}}
end

# 2. Scheduling policies: what the scheduler maximizes. A policy is a value of a
# policy type; it scores each available candidate link (`candidate_score`), and
# the scheduler activates the terminal matching with the largest total score.
"""
    AbstractSchedulingPolicy

Supertype of the interlink scheduling policies. A policy is an immutable value
whose fields say what it scores, and it holds no run state. A policy type
implements `candidate_score`. The built-in policies are
`SemimajorAxisRatePolicy`, `EccentricityRatePolicy`, `PowerDeliveryPolicy`,
`DataThroughputPolicy`, `DataBacklogPolicy`, `EnergyBalancingPolicy` and
`WeightedPolicy`, a weighted mix of policies. Give a run its policy with
`SimulationConfiguration(scheduling_policy_model=policy)`.
"""
abstract type AbstractSchedulingPolicy end

"""
    LinkCandidate

One available candidate link, as a scheduling policy scores it.

- `key`: its endpoints, `(spacecraft vector index, terminal number)` each, in key
  order: the transmitter first for a directed link, the lower endpoint first
  otherwise. A ground station end is `(-s, antenna number)`.
- `link`: its link kind value.
- `a`, `b`: the endpoints' state samples (`StateSample`), in key order, with
  their spacecraft models when the run's are known. A station end's sample has
  the station as its model, no mass and no attitude.
- `terminal_a`, `terminal_b`: the terminal hardware at the two ends, `nothing`
  for a virtual terminal or when the spacecraft models are not known; a
  station's antenna at a station end.
"""
struct LinkCandidate
    key::InterLinkKey
    link::AbstractInterlink
    a::StateSample
    b::StateSample
    terminal_a::Any
    terminal_b::Any
end

"""
    SchedulingContext

What a scheduling policy may read besides the candidate:

- `mu`: the planet's gravitational parameter (m³/s²);
- `link_context`: the time and the planet (`LinkContext`), to pass to the link
  hooks;
- `spacecraft`: the run's spacecraft models, or `nothing` when a caller scores
  without them.

Its other fields are the engine's.
"""
struct SchedulingContext{L <: LinkContext, S, X, F, T, O}
    mu::Float64
    link_context::L
    spacecraft::S
    # The run state, which the built-in orbital policies read as the engine
    # holds it; a caller's `candidate_force(target, partner, link)`; and each
    # spacecraft's terminal hardware by terminal number, looked up once per
    # scoring call (`nothing` without the spacecraft models).
    _state::X
    _candidate_force::F
    _terminals::T
    _outflows::O # carried rates from the previous matching, sampled once per pass
end
SchedulingContext(mu::Real, link_context::LinkContext, spacecraft, state, candidate_force) =
    SchedulingContext(Float64(mu), link_context, spacecraft, state, candidate_force,
        spacecraft === nothing ? nothing : map(_terminal_devices, spacecraft), nothing)
SchedulingContext(mu::Real, link_context::LinkContext, spacecraft, state, candidate_force, terminals) =
    SchedulingContext(Float64(mu), link_context, spacecraft, state, candidate_force, terminals, nothing)

"""
    candidate_score(policy, candidate::LinkCandidate, ctx::SchedulingContext, t::Float64) -> Float64

The score `policy` gives `candidate`, an available candidate link, at `t` s into
the run (`ctx.link_context.time_s`; NaN outside a run). The scheduler activates
the terminal matching with the largest total score, and only links with a
positive score can be selected. The engine calls it at the start of the run and
after every accepted step, once for each available candidate; an unavailable
candidate scores 0 without a call. A nonfinite score is refused. Every policy
type implements it: a run refuses at setup a policy without a method.
"""
function candidate_score end

# The candidate of graph edge `edge`, at the context's state. A ground station
# end reads the planet's current frame, as the availability does.
function _link_candidate(model, edge::Int, ctx::SchedulingContext)
    graph = model.graph
    first, second = map(vertex -> graph.vertices[vertex], graph.edges[edge])
    spacecraft, state = ctx.spacecraft, ctx._state
    first[1] > 0 && second[1] > 0 || return _station_link_candidate(model, edge, ctx)
    a = link_endpoint(state.sc[first[1]], _spacecraft_at(spacecraft, first[1]), first[1])
    b = link_endpoint(state.sc[second[1]], _spacecraft_at(spacecraft, second[1]), second[1])
    return LinkCandidate((first, second), graph.edge_properties[edge], a, b,
        _candidate_terminal(ctx._terminals, first), _candidate_terminal(ctx._terminals, second))
end
function _station_link_candidate(model, edge::Int, ctx::SchedulingContext)
    graph = model.graph
    first, second = map(vertex -> graph.vertices[vertex], graph.edges[edge])
    a = _endpoint_sample(model, first[1], ctx._state, ctx.spacecraft, ctx.link_context)
    b = _endpoint_sample(model, second[1], ctx._state, ctx.spacecraft, ctx.link_context)
    terminal(endpoint) = endpoint[1] < 0 ? _station_antenna(model, endpoint) : _candidate_terminal(ctx._terminals, endpoint)
    return LinkCandidate((first, second), graph.edge_properties[edge], a, b, terminal(first), terminal(second))
end
_candidate_terminal(::Nothing, endpoint) = nothing
_candidate_terminal(terminals, (satellite, number)) = _candidate_terminal(terminals[satellite], number)
_candidate_terminal(::Nothing, number::Int) = nothing
_candidate_terminal(devices::AbstractVector, number::Int) = devices[number]

# Whether a built-in policy can give graph edge `edge` a score other than 0. The
# loop scores an edge it cannot as 0, without building its candidate: an edge
# with no end among an orbital policy's targets, or one that the power and data
# policies do not score. Any other policy is called for every available edge.
_may_score(policy, model, edge::Int, ctx::SchedulingContext) = true

# Checks a policy makes of the run before scoring (its targets, the gravitational
# parameter), and checks across the candidates' scores after it. Both default to
# none.
_prepare_scoring(policy, model, ctx::SchedulingContext) = nothing
_check_scores(policy, model, scores, ctx::SchedulingContext) = nothing

# The run's policy has a `candidate_score` method: checked at setup.
function _check_candidate_score(policy)
    hasmethod(candidate_score, Tuple{typeof(policy), LinkCandidate, SchedulingContext, Float64}) || throw(ArgumentError(
        "The scheduling policy $(nameof(typeof(policy))) has no candidate_score method; implement " *
        "SpaceAGORA.candidate_score(policy::$(nameof(typeof(policy))), candidate, ctx, t), the score of one candidate link."))
    return nothing
end

# ---- The orbital policies ---------------------------------------------------------

# The link's forces on the candidate's two ends (`interlink_wrench`'s, by default
# `link_force`'s), evaluated once for both, or `nothing` when a caller's
# `candidate_force` gives the forces.
function _link_forces(ctx::SchedulingContext, candidate::LinkCandidate)
    ctx._candidate_force === nothing || return nothing
    force_a, _, force_b, _ = interlink_wrench(candidate.link, candidate.terminal_a, candidate.terminal_b,
        candidate.a, candidate.b, ctx.link_context)
    return (force_a, force_b)
end

# The force the candidate puts on spacecraft `target`: from the link's forces,
# or from a caller's `candidate_force(target, partner, link)`. The caller's
# function is called dynamically, so it is compiled on its own as it was before
# the context had field types, and its arithmetic rounds the same.
function _target_force(ctx::SchedulingContext, candidate::LinkCandidate, target::Int, partner::Int, forces)
    forces === nothing || return target == candidate.key[1][1] ? forces[1] : forces[2]
    return Base.inferencebarrier(ctx._candidate_force)(ctx._state.sc[target], _partner_state(ctx, candidate, partner),
        candidate.link)
end
# The partner a caller's `candidate_force` receives: its run state, or for a
# ground station, its state sample.
_partner_state(ctx::SchedulingContext, candidate::LinkCandidate, partner::Int) = partner > 0 ? ctx._state.sc[partner] :
    candidate.key[1][1] == partner ? candidate.a : candidate.b

# The rate (m/s) at which `force` changes the osculating semimajor axis of the
# spacecraft state `target`, as the engine holds it. The rate functions are
# compiled on their own (`@noinline`), not into each caller: StaticArrays' dot
# product rounds by its compile context, so every caller gets the same bits.
@noinline function _semimajor_axis_rate(target, force, mu::Real)
    semimajor_axis = inv(2 / norm(target.pos) - dot(target.vel, target.vel) / mu)
    return 2 * semimajor_axis^2 / mu * dot(target.vel, force) / target.mass
end

_check_targets(name, targets) = throw(ArgumentError(
    "$(name) targets must be a positive spacecraft vector index$(name === :SemimajorAxisRatePolicy ?
        ", a nonempty vector of distinct positive indices, :both or :lower_index" : " or :both"); got $(repr(targets))."))

"""
    SemimajorAxisRatePolicy(; targets=1)

Score each available link by the rate (m/s) at which its force raises the
osculating semimajor axis of its target spacecraft:
`da/dt = 2a² (v · F) / (μ m)` at each target end, summed. The force is the
link's, from `interlink_wrench` (by default its `link_force`). `targets` picks the ends:

- a spacecraft vector index: that spacecraft's end only;
- a vector or range of indices: every end in the set, stored sorted;
- `:both`: both ends of every link;
- `:lower_index`: the lower spacecraft index of each pair.

A ground station end is never a target: a station has no orbit, and a link to
one scores at its spacecraft end alone (`:lower_index` takes that end).

A link that lowers the semimajor axis scores negative, and links that apply no
force (power beams, data links) score 0, so this policy never selects them.
"""
struct SemimajorAxisRatePolicy <: AbstractSchedulingPolicy
    targets::Union{Int, Vector{Int}, Symbol}

    function SemimajorAxisRatePolicy(targets)
        name = :SemimajorAxisRatePolicy
        if targets isa Symbol
            targets in (:both, :lower_index) || _check_targets(name, targets)
            return new(targets)
        elseif targets isa Integer && !(targets isa Bool)
            targets > 0 || _check_targets(name, targets)
            return new(Int(targets))
        elseif targets isa AbstractVector{<:Integer} && !(eltype(targets) <: Bool)
            selected = sort!(Int[targets...])
            !isempty(selected) && all(>(0), selected) && allunique(selected) || _check_targets(name, targets)
            return new(selected)
        end
        _check_targets(name, targets)
    end
end
SemimajorAxisRatePolicy(; targets=1) = SemimajorAxisRatePolicy(targets)

_is_target(targets::Int, satellite::Int) = targets == satellite
# The `:lower_index` target of a link and its partner: the lower spacecraft
# index, or the spacecraft end of a link to a ground station.
_lower_spacecraft_end(first::Int, second::Int) = first < 0 ? (second, first) : second < 0 ? (first, second) :
    minmax(first, second)
_is_target(targets::Vector{Int}, satellite::Int) = satellite in targets
_is_target(targets::Symbol, satellite::Int) = targets === :both && satellite > 0

function _prepare_scoring(policy::SemimajorAxisRatePolicy, model, ctx::SchedulingContext)
    isfinite(ctx.mu) && ctx.mu > 0 || throw(ArgumentError("The gravitational parameter must be finite and positive."))
    count = length(ctx._state.sc)
    targets = policy.targets
    targets isa Vector{Int} && !all(index -> index <= count, targets) &&
        throw(ArgumentError("SemimajorAxisRatePolicy targets contain an index outside the spacecraft state vector."))
    targets isa Int && targets > count &&
        throw(ArgumentError("The SemimajorAxisRatePolicy target is outside the spacecraft state vector."))
    return nothing
end

function _may_score(policy::SemimajorAxisRatePolicy, model, edge::Int, ctx::SchedulingContext)
    policy.targets isa Symbol && return true
    graph = model.graph
    first, second = graph.edges[edge]
    return _is_target(policy.targets, graph.vertices[first][1]) || _is_target(policy.targets, graph.vertices[second][1])
end

function candidate_score(policy::SemimajorAxisRatePolicy, candidate::LinkCandidate, ctx::SchedulingContext, t::Float64)
    first_satellite, second_satellite = candidate.key[1][1], candidate.key[2][1]
    if policy.targets === :lower_index
        target, partner = _lower_spacecraft_end(first_satellite, second_satellite)
        force = _target_force(ctx, candidate, target, partner, _link_forces(ctx, candidate))
        return _semimajor_axis_rate(ctx._state.sc[target], force, ctx.mu)
    end
    score, forces = 0.0, nothing
    for (target, partner) in ((first_satellite, second_satellite), (second_satellite, first_satellite))
        _is_target(policy.targets, target) || continue
        forces === nothing && (forces = _link_forces(ctx, candidate))
        score += _semimajor_axis_rate(ctx._state.sc[target], _target_force(ctx, candidate, target, partner, forces), ctx.mu)
    end
    return score
end

"""
    EccentricityRatePolicy(; targets=1)

Score each available link by the rate (1/s) at which its force raises the
eccentricity of its target spacecraft, `de/dt = ê · ė`, at each target end,
summed. `targets` is a spacecraft vector index, or `:both` for both ends of
every link. The force is the link's, from `interlink_wrench` (by default its
`link_force`). A ground station end is never a target.

At an eccentricity no greater than `sqrt(eps(Float64))` the score is the norm of
the eccentricity-vector rate, `‖ė‖`. Such norms are not additive, so a candidate
set that could put two positive-score links on one such target, through
different terminals, is refused. A link that lowers the eccentricity scores
negative, and links that apply no force score 0.
"""
struct EccentricityRatePolicy <: AbstractSchedulingPolicy
    targets::Union{Int, Symbol}

    function EccentricityRatePolicy(targets)
        name = :EccentricityRatePolicy
        targets === :both && return new(targets)
        targets isa Integer && !(targets isa Bool) && targets > 0 && return new(Int(targets))
        _check_targets(name, targets)
    end
end
EccentricityRatePolicy(; targets=1) = EccentricityRatePolicy(targets)

function _prepare_scoring(policy::EccentricityRatePolicy, model, ctx::SchedulingContext)
    isfinite(ctx.mu) && ctx.mu > 0 || throw(ArgumentError("The gravitational parameter must be finite and positive."))
    policy.targets isa Int && policy.targets > length(ctx._state.sc) &&
        throw(ArgumentError("The EccentricityRatePolicy target is outside the spacecraft state vector."))
    return nothing
end

# The eccentricity rate `force` gives the spacecraft state `target`, and whether
# its orbit is near circular (the rate is then the norm of the vector rate).
@noinline function _eccentricity_rate(target, force, mu::Real)
    position = target.pos
    velocity = target.vel
    acceleration = force / target.mass
    angular_momentum = cross(position, velocity)
    eccentricity_vector = cross(velocity, angular_momentum) / mu - position / norm(position)
    eccentricity_vector_rate = (
        cross(acceleration, angular_momentum) + cross(velocity, cross(position, acceleration))
    ) / mu
    eccentricity = norm(eccentricity_vector)
    near_circular = eccentricity <= sqrt(eps(Float64))
    rate = near_circular ? norm(eccentricity_vector_rate) :
        dot(eccentricity_vector, eccentricity_vector_rate) / eccentricity
    return rate, near_circular
end

function _may_score(policy::EccentricityRatePolicy, model, edge::Int, ctx::SchedulingContext)
    policy.targets === :both && return true
    graph = model.graph
    first, second = graph.edges[edge]
    return policy.targets == graph.vertices[first][1] || policy.targets == graph.vertices[second][1]
end

candidate_score(policy::EccentricityRatePolicy, candidate::LinkCandidate, ctx::SchedulingContext, t::Float64) =
    _eccentricity_score(policy, candidate, ctx, nothing, 0)

# The score of the candidate of graph edge `edge`. With `circular`, it also
# records each target spacecraft on a near-circular orbit that the candidate
# gives a positive rate.
function _eccentricity_score(policy::EccentricityRatePolicy, candidate::LinkCandidate, ctx::SchedulingContext, circular, edge::Int)
    first_satellite, second_satellite = candidate.key[1][1], candidate.key[2][1]
    score, forces = 0.0, nothing
    for (target, partner) in ((first_satellite, second_satellite), (second_satellite, first_satellite))
        target > 0 && (policy.targets === :both || policy.targets == target) || continue
        forces === nothing && (forces = _link_forces(ctx, candidate))
        rate, near_circular = _eccentricity_rate(ctx._state.sc[target], _target_force(ctx, candidate, target, partner, forces), ctx.mu)
        score += rate
        circular !== nothing && near_circular && rate > 0 && push!(get!(Vector{Int}, circular, target), edge)
    end
    return score
end

# The scheduler's own loop records the near-circular candidates while it scores
# them; the check of a policy inside another (a mix) computes them again. Only a
# mix has terms to record.
function _score_and_check!(scores, policy::EccentricityRatePolicy, model, ctx::SchedulingContext; terms=nothing)
    circular = Dict{Int, Vector{Int}}()
    _score_candidates!(scores, policy, model, ctx, circular)
    return _check_circular_candidates(model, scores, circular)
end

function _check_scores(policy::EccentricityRatePolicy, model, scores, ctx::SchedulingContext)
    graph = model.graph
    circular = Dict{Int, Vector{Int}}()
    for edge in eachindex(graph.edges)
        graph.available[edge] && _may_score(policy, model, edge, ctx) || continue
        _eccentricity_score(policy, _link_candidate(model, edge, ctx), ctx, circular, edge)
    end
    return _check_circular_candidates(model, scores, circular)
end

# Near circularity the score is the norm of the eccentricity-vector rate. Norms
# from two simultaneous forces are not additive: opposite forces can receive
# positive individual scores even though their combined rate is zero. Admit the
# objective only when a scored near-circular spacecraft cannot receive two
# positive-score candidate links in the same terminal matching.
function _check_circular_candidates(model, scores, circular_candidates)
    graph = model.graph
    for (satellite, edges) in circular_candidates
        candidates = Tuple{Int, Int}[]
        for edge in edges
            graph.available[edge] && scores[edge] > 0 || continue
            endpoints = graph.edges[edge]
            for previous in candidates
                if all(vertex -> vertex != endpoints[1] && vertex != endpoints[2], previous)
                    throw(ArgumentError("The EccentricityRatePolicy objective does not support simultaneous links " *
                        "on near-circular target spacecraft $satellite: per-link eccentricity-rate norms " *
                        "are not additive. Use candidate links that share one terminal on that target."))
                end
            end
            push!(candidates, endpoints)
        end
    end
    return nothing
end

# ---- Power and data ----------------------------------------------------------------

"""
    PowerDeliveryPolicy()

Score each available directed link whose receiving end is terminal hardware by
the electrical power (W) it would deliver to the receiver's bus,
`link_received_power_w` at the current state. It selects power beams, which
apply no force. Other links score 0, and so does every link when the run's
spacecraft models are not given.
"""
struct PowerDeliveryPolicy <: AbstractSchedulingPolicy end

_may_score(::PowerDeliveryPolicy, model, edge::Int, ctx::SchedulingContext) =
    ctx.spacecraft !== nothing && link_directed(model.graph.edge_properties[edge])

function candidate_score(::PowerDeliveryPolicy, candidate::LinkCandidate, ctx::SchedulingContext, t::Float64)
    ctx.spacecraft === nothing && return 0.0
    link_directed(candidate.link) || return 0.0
    candidate.terminal_b === nothing && return 0.0
    return Float64(link_received_power_w(candidate.link, candidate.terminal_b, candidate.a, candidate.b, ctx.link_context))
end

"""
    DataThroughputPolicy()

Score each available directed link by the data rate (bit/s) it would carry,
`link_data_rate_bps` at the current state, with the terminal hardware at its two
ends (`nothing` for a virtual terminal). It selects data links. Undirected links
score 0, and so does every link when the run's spacecraft models are not given.
"""
struct DataThroughputPolicy <: AbstractSchedulingPolicy end

_may_score(::DataThroughputPolicy, model, edge::Int, ctx::SchedulingContext) =
    ctx.spacecraft !== nothing && link_directed(model.graph.edge_properties[edge])

function candidate_score(::DataThroughputPolicy, candidate::LinkCandidate, ctx::SchedulingContext, t::Float64)
    ctx.spacecraft === nothing && return 0.0
    link_directed(candidate.link) || return 0.0
    return Float64(link_data_rate_bps(candidate.link, candidate.terminal_a, candidate.terminal_b,
        candidate.a, candidate.b, ctx.link_context))
end

"""
    stored_bits(ctx::SchedulingContext, spacecraft_index::Integer) -> Float64

The data (bit) spacecraft `spacecraft_index`'s data store holds at the
scheduling state, the `data_stored_bits` column's value. NaN without a store,
without the run's spacecraft models, or for a state without `data_stored_bits`.
Inf for a ground station (a negative index), an unbounded store that always
has data to send.
"""
function stored_bits(ctx::SchedulingContext, spacecraft_index::Integer)
    spacecraft_index < 0 && return Inf
    ctx.spacecraft === nothing && return NaN
    isnan(_data_store_capacity_bits(ctx.spacecraft[spacecraft_index])) && return NaN
    view = ctx._state.sc[spacecraft_index]
    hasproperty(view, :data_stored_bits) || return NaN
    return Float64(view.data_stored_bits)
end

"""
    store_room_bits(ctx::SchedulingContext, spacecraft_index::Integer) -> Float64

The room (bit) left in spacecraft `spacecraft_index`'s data store at the
scheduling state: its capacity less `stored_bits`. Inf without a store or
without the run's spacecraft models, since nothing then limits what the
spacecraft takes in; NaN for a state without `data_stored_bits`. Inf for a
ground station (a negative index), an unbounded store that always has room.
"""
function store_room_bits(ctx::SchedulingContext, spacecraft_index::Integer)
    spacecraft_index < 0 && return Inf
    ctx.spacecraft === nothing && return Inf
    capacity = _data_store_capacity_bits(ctx.spacecraft[spacecraft_index])
    isnan(capacity) && return Inf
    return capacity - stored_bits(ctx, spacecraft_index)
end

"""
    DataBacklogPolicy(; horizon_s)

Score each available directed link by the data (as a rate, bit/s) it could move
in the next `horizon_s` seconds: the least of the rate it would carry
(`link_data_rate_bps`, as `DataThroughputPolicy` scores it), the data its
transmitter's store holds over `horizon_s` (`stored_bits`) plus its current
source production rate, and the room in its receiver's store over `horizon_s`
(`store_room_bits`) plus its current carried outflow. Outflow is sampled once
before selection from links that were already active and remain available,
using the same flow caps as the runtime. It does not predict the new matching.
Production is sampled at the candidate's current state and
time; this score estimates that rate continuing over the horizon, not future
source changes or incoming relay traffic. It prefers fuller sending stores and
emptier receiving stores, while a full store that is draining can keep an incoming link. An empty sender with ongoing production can stay
selected; one with neither held data nor production scores zero. A store the
engine counts as empty (within `1e-9` of its capacity of nothing) holds nothing
here; a full receiver has only its current outflow as scoring room. A transmitter without a store counts as always holding data, and a
receiver without one as always having room; so does a ground station, which
never limits a link: the policy drains fuller stores to the ground first. Undirected
links score 0, and so does every link when the run's spacecraft models are not
given. The unit is that of `DataThroughputPolicy`, so the two mix without
rescaling.
"""
struct DataBacklogPolicy <: AbstractSchedulingPolicy
    horizon_s::Float64

    function DataBacklogPolicy(horizon_s::Real)
        isfinite(horizon_s) && horizon_s > 0 || throw(ArgumentError(
            "A DataBacklogPolicy's horizon_s must be finite and positive, in seconds; got $(horizon_s)."))
        return new(Float64(horizon_s))
    end
end
DataBacklogPolicy(; horizon_s::Real) = DataBacklogPolicy(horizon_s)

_may_score(::DataBacklogPolicy, model, edge::Int, ctx::SchedulingContext) =
    ctx.spacecraft !== nothing && link_directed(model.graph.edge_properties[edge])

function candidate_score(policy::DataBacklogPolicy, candidate::LinkCandidate, ctx::SchedulingContext, t::Float64)
    ctx.spacecraft === nothing && return 0.0
    link_directed(candidate.link) || return 0.0
    rate = Float64(link_data_rate_bps(candidate.link, candidate.terminal_a, candidate.terminal_b,
        candidate.a, candidate.b, ctx.link_context))
    transmitter, receiver = candidate.key[1][1], candidate.key[2][1]
    held = _backlog_held_bits(ctx, transmitter)
    production = isfinite(held) ?
        _data_store_production_bps(ctx.spacecraft[transmitter], candidate.a, ctx.link_context, t) : 0.0
    outflow = receiver > 0 && ctx._outflows !== nothing ? ctx._outflows[receiver] : 0.0
    return min(rate, held / policy.horizon_s + production,
        _backlog_room_bits(ctx, receiver) / policy.horizon_s + outflow)
end

# The data a transmitter holds and the room a receiver has, for the backlog
# score: as the store's own events count it, so the hair a full or empty event
# leaves does not count; Inf without a store and at a ground station.
function _backlog_held_bits(ctx::SchedulingContext, i::Int)
    i < 0 && return Inf
    held = stored_bits(ctx, i)
    isnan(held) && return Inf
    return _data_store_within_empty(ctx.spacecraft[i], held) ? 0.0 : held
end
function _backlog_room_bits(ctx::SchedulingContext, i::Int)
    i < 0 && return Inf
    room = store_room_bits(ctx, i)
    isnan(room) && return Inf
    return _data_store_within_full(ctx.spacecraft[i], stored_bits(ctx, i)) ? 0.0 : room
end

"""
    state_of_charge(ctx::SchedulingContext, spacecraft_index::Integer) -> Float64

The state of charge (a fraction of its capacity) of spacecraft
`spacecraft_index`'s battery at the scheduling state: its `battery_wh` over its
`capacity_wh`, the `battery_soc` column's value. NaN without a battery, without
the run's spacecraft models, for a state without `battery_wh`, or for a ground
station (a negative index), which has no battery.
"""
function state_of_charge(ctx::SchedulingContext, spacecraft_index::Integer)
    ctx.spacecraft === nothing && return NaN
    spacecraft_index < 0 && return NaN
    view = ctx._state.sc[spacecraft_index]
    hasproperty(view, :battery_wh) || return NaN
    return Float64(view.battery_wh) / _battery_capacity_wh(ctx.spacecraft[spacecraft_index])
end

# ---- Energy balancing ----------------------------------------------------------

"""
    EnergyBalancingPolicy()

Charge the lowest state of charge. Score each available directed link whose
receiving end is terminal hardware, between two spacecraft with a battery, by
`P × (s_tx − s_rx)` (W): `P` is the electrical power it would deliver to the
receiver's bus (`link_received_power_w`, as `PowerDeliveryPolicy` scores it),
and `s_tx` and `s_rx` are the two batteries' states of charge
(`state_of_charge`). Power flows only from a fuller battery to an emptier one,
and a beam into a full battery never scores: a battery the power system counts
as full (its stored energy within `1e-9` of its capacity) counts as 1, so two
full batteries tie at 0. Of receivers competing for one sender's terminal, the
emptier wins at equal delivered power. The score is signed: a beam from an
emptier battery to a fuller one scores negative. Links with no battery at
either end, a link that would deliver no power (or, for a link kind of your
own, a negative power), and every link without the run's spacecraft models
score 0.
"""
struct EnergyBalancingPolicy <: AbstractSchedulingPolicy end

# The state of charge energy balancing compares: 1 for a battery the power
# system counts as full (`within_full_tolerance`, in energy, as its flags
# decide it). The full event leaves the stored energy just below the capacity,
# so a full battery's raw state of charge can be a hair under 1.
function _balancing_charge(ctx::SchedulingContext, spacecraft_index::Int)
    charge = state_of_charge(ctx, spacecraft_index)
    isnan(charge) && return charge
    energy_wh = Float64(ctx._state.sc[spacecraft_index].battery_wh)
    return _battery_within_full(ctx.spacecraft[spacecraft_index], energy_wh) ? 1.0 : charge
end

function candidate_score(::EnergyBalancingPolicy, candidate::LinkCandidate, ctx::SchedulingContext, t::Float64)
    ctx.spacecraft === nothing && return 0.0
    link_directed(candidate.link) || return 0.0
    candidate.terminal_b === nothing && return 0.0
    gap = _balancing_charge(ctx, candidate.key[1][1]) - _balancing_charge(ctx, candidate.key[2][1])
    isnan(gap) && return 0.0
    delivered = Float64(link_received_power_w(candidate.link, candidate.terminal_b, candidate.a, candidate.b, ctx.link_context))
    delivered <= 0 && return 0.0
    return delivered * gap
end

# ---- Weighted mixes ------------------------------------------------------------

"""
    PolicyTerm(policy::AbstractSchedulingPolicy; scale::Real, weight::Real=1.0)

One term of a `WeightedPolicy`: `policy`'s score divided by `scale` and times
`weight`. `scale` is in the policy's own unit (m/s, W, bit/s, ...), the amount
of its score that counts as one, so it makes the term dimensionless; it must be
finite and positive. `weight` must be finite; a negative weight penalizes the
term.
"""
struct PolicyTerm{P <: AbstractSchedulingPolicy}
    policy::P
    scale::Float64
    weight::Float64

    function PolicyTerm(policy::P, scale::Real, weight::Real) where {P <: AbstractSchedulingPolicy}
        isfinite(scale) && scale > 0 || throw(ArgumentError("A PolicyTerm's scale must be finite and positive, in its policy's unit; got $(scale)."))
        isfinite(weight) || throw(ArgumentError("A PolicyTerm's weight must be finite; got $(weight)."))
        new{P}(policy, Float64(scale), Float64(weight))
    end
end
PolicyTerm(policy::AbstractSchedulingPolicy; scale::Real, weight::Real=1.0) = PolicyTerm(policy, scale, weight)

"""
    WeightedPolicy(terms::PolicyTerm...)

A weighted mix of policies: each available link scores
`Σ weight × (score / scale)` over the terms, summed in term order from 0.0. Each
term's policy applies its own rules (an orbital term scores power beams 0, a
power term scores momentum links 0), so one run can weigh link kinds against
each other. Each term's checks run, and an eccentricity term keeps its
near-circular refusal over the mix's own candidates. A nonfinite term is
refused. Mixes nest. A term with `weight=1` and `scale=1` alone schedules as its
policy does.

```julia
WeightedPolicy(PolicyTerm(SemimajorAxisRatePolicy(targets=1); scale=1.0),   # 1 m/s counts as one
               PolicyTerm(EnergyBalancingPolicy(); scale=100.0))            # and so do 100 W
```
"""
struct WeightedPolicy{T <: Tuple{Vararg{PolicyTerm}}} <: AbstractSchedulingPolicy
    terms::T

    function WeightedPolicy(terms::T) where {T <: Tuple{Vararg{PolicyTerm}}}
        isempty(terms) && throw(ArgumentError("A WeightedPolicy needs at least one PolicyTerm."))
        new{T}(terms)
    end
end
WeightedPolicy(terms::PolicyTerm...) = WeightedPolicy(terms)

candidate_score(policy::WeightedPolicy, candidate::LinkCandidate, ctx::SchedulingContext, t::Float64) =
    _weighted_score(policy, candidate, ctx, t, nothing, 0)

# The mix's score, and with `terms` each term's own score in column `edge`: one
# loop for both, so recording the terms leaves the mix's arithmetic unchanged.
function _weighted_score(policy::WeightedPolicy, candidate::LinkCandidate, ctx::SchedulingContext, t::Float64,
    terms, edge::Int)
    score = 0.0
    for (k, term) in enumerate(policy.terms)
        value = Float64(candidate_score(term.policy, candidate, ctx, t))
        isfinite(value) || throw(ArgumentError("Nonfinite $(nameof(typeof(term.policy))) score in a WeightedPolicy " *
            "term for candidate $(candidate.key)."))
        terms === nothing || (terms[k, edge] = value)
        score += term.weight * (value / term.scale)
    end
    return score
end

function _prepare_scoring(policy::WeightedPolicy, model, ctx::SchedulingContext)
    foreach(term -> _prepare_scoring(term.policy, model, ctx), policy.terms)
    return nothing
end
# The terms' checks across candidates use the mix's scores: an eccentricity term
# refuses a near-circular target that two positive-score links of the mix could reach.
function _check_scores(policy::WeightedPolicy, model, scores, ctx::SchedulingContext)
    foreach(term -> _check_scores(term.policy, model, scores, ctx), policy.terms)
    return nothing
end
function _check_candidate_score(policy::WeightedPolicy)
    foreach(term -> _check_candidate_score(term.policy), policy.terms)
    return nothing
end

# ---- The scheduling policies chosen by a symbol, removed in 0.6.0 ---------------

# The typed policy a removed symbol form stands for: its type and its targets
# (`nothing` for a type without targets), from the old `target_index` and
# `target_indices`. One table for the refusal message and the checkpoint reader.
function _symbol_policy_spec(policy::Symbol, index, indices)
    targets = indices !== nothing ? indices : index === nothing ? :both : index
    policy === :gve_sma && return (SemimajorAxisRatePolicy, targets)
    policy === :gve_sma_first_index && return (SemimajorAxisRatePolicy, :lower_index)
    policy === :gve_eccentricity && return (EccentricityRatePolicy, targets)
    policy === :power_delivery && return (PowerDeliveryPolicy, nothing)
    policy === :data_throughput && return (DataThroughputPolicy, nothing)
    return nothing
end

# The typed policy a checkpoint's recorded symbol and targets stand for, or
# `nothing` when they stand for none.
function _symbol_scheduling_policy(policy::Symbol, index, indices)
    spec = _symbol_policy_spec(policy, index, indices)
    spec === nothing && return nothing
    T, targets = spec
    try
        return targets === nothing ? T() : T(targets)
    catch err
        err isa ArgumentError || rethrow()
        return nothing
    end
end

# The symbol, target index and target set that a checkpoint records for a policy
# that had a symbol form, as 0.5 recorded them; `nothing` for any other policy.
_symbol_policy_fields(policy) = nothing
function _symbol_policy_fields(policy::SemimajorAxisRatePolicy)
    targets = policy.targets
    targets === :lower_index && return (:gve_sma_first_index, 1, nothing)
    targets === :both && return (:gve_sma, nothing, nothing)
    targets isa Int && return (:gve_sma, targets, nothing)
    return (:gve_sma, 1, copy(targets))
end
_symbol_policy_fields(policy::EccentricityRatePolicy) =
    (:gve_eccentricity, policy.targets === :both ? nothing : policy.targets, nothing)
_symbol_policy_fields(::PowerDeliveryPolicy) = (:power_delivery, 1, nothing)
_symbol_policy_fields(::DataThroughputPolicy) = (:data_throughput, 1, nothing)

const _TARGET_KEYWORDS = (:target_index, :target_indices, :target_idx, :target_idxs)

# The message for a removed symbol form: the call as given, and its typed replacement.
function _removed_symbol_message(policy::Symbol, kwargs)
    given(name) = haskey(kwargs, name)
    index = given(:target_index) ? kwargs[:target_index] : given(:target_idx) ? kwargs[:target_idx] : 1
    indices = given(:target_indices) ? kwargs[:target_indices] : given(:target_idxs) ? kwargs[:target_idxs] : nothing
    call = "SchedulingPolicyModel($(repr(policy))" *
        (isempty(kwargs) ? "" : "; " * join(("$(name)=$(repr(value))" for (name, value) in pairs(kwargs)), ", ")) * ")"
    spec = _symbol_policy_spec(policy, index, indices)
    if spec === nothing
        return "$(call): $(repr(policy)) is not a scheduling policy, and policies chosen by a symbol were removed in " *
            "SpaceAGORA 0.6.0. Use a policy type: SemimajorAxisRatePolicy, EccentricityRatePolicy, PowerDeliveryPolicy, " *
            "DataThroughputPolicy, DataBacklogPolicy, EnergyBalancingPolicy or WeightedPolicy, or one of your own (a subtype of " *
            "AbstractSchedulingPolicy with a candidate_score method)."
    end
    T, targets = spec
    shown = targets !== nothing && (policy === :gve_sma_first_index || any(given, _TARGET_KEYWORDS))
    replacement = "$(nameof(T))(" * (shown ? "targets=$(repr(targets))" : "") * ")"
    others = [name for name in keys(kwargs) if !(name in _TARGET_KEYWORDS)]
    return "$(call) was removed in SpaceAGORA 0.6.0: scheduling policies are types. Use $(replacement), for example " *
        "scheduling_policy_model=$(replacement)." *
        (isempty(others) ? "" : " The keyword$(length(others) == 1 ? "" : "s") $(join(others, ", ")) " *
            "$(length(others) == 1 ? "has" : "have") no replacement.")
end

# 2.1. The scheduler
"""
    SchedulingPolicyModel(policy::AbstractSchedulingPolicy)

The run's interlink scheduler: the scheduling policy whose total score it
maximizes (`objective`) and the matching workspace it reuses between
scheduling calls (`matching`). Pass it, or the policy alone, as
`SimulationConfiguration(scheduling_policy_model=...)`; the default is
`SchedulingPolicyModel(SemimajorAxisRatePolicy())`.

The policies chosen by a symbol, `SchedulingPolicyModel(:gve_sma)` and the
others, were removed in SpaceAGORA 0.6.0, with their `target_index` and
`target_indices` keywords and properties and the older `target_idx` and
`target_idxs`. A call with a symbol throws an `ArgumentError` that names its
typed replacement, such as `SemimajorAxisRatePolicy(targets=[31, 32])` for
`SchedulingPolicyModel(:gve_sma; target_indices=[31, 32])`, and so does a read
or write of a removed property, naming `.objective` or `.objective.targets`.
"""
mutable struct SchedulingPolicyModel
	const objective::AbstractSchedulingPolicy
	matching::Union{Nothing, InterLinkMatchingWorkspace}

	SchedulingPolicyModel(objective::AbstractSchedulingPolicy) = new(objective, nothing)
end
SchedulingPolicyModel(policy::Symbol; kwargs...) = throw(ArgumentError(_removed_symbol_message(policy, kwargs)))

# The properties removed in 0.6.0 refuse with their replacement; every other name
# is a field. With a constant name, as in `scheduler.objective`, the check folds
# away at compile time.
const _REMOVED_PROPERTIES = (:policy, :target_index, :target_indices, :target_idx, :target_idxs)
@inline function Base.getproperty(model::SchedulingPolicyModel, name::Symbol)
	name in _REMOVED_PROPERTIES && _removed_property(name)
	return getfield(model, name)
end
@inline function Base.setproperty!(model::SchedulingPolicyModel, name::Symbol, value)
	name in _REMOVED_PROPERTIES && _removed_property(name)
	return setfield!(model, name, convert(fieldtype(SchedulingPolicyModel, name), value))
end
@noinline _removed_property(name::Symbol) = throw(ArgumentError(name === :policy ?
	"SchedulingPolicyModel.policy was removed in SpaceAGORA 0.6.0: the scheduler's policy is a typed policy, " *
	"`.objective` (for example SemimajorAxisRatePolicy(targets=1))." :
	"SchedulingPolicyModel.$(name) was removed in SpaceAGORA 0.6.0: a policy's targets are `.objective.targets` " *
	"(SemimajorAxisRatePolicy, EccentricityRatePolicy)."))

# 1.2. Create and maintain the dynamic-programming workspace.
function matching_workspace!(model::InterlinkModel, policy::SchedulingPolicyModel)
	graph = model.graph
	vertices, edges = length(graph.vertices), length(graph.edges)
	work = policy.matching
	term_count = model.record_terms && policy.objective isa WeightedPolicy ? length(policy.objective.terms) : 0
	if work !== nothing && work.graph === graph &&
		length(work.parent) == vertices && length(work.weights) == edges &&
		(work.terms === nothing ? term_count == 0 : size(work.terms) == (term_count, edges))
		return work
	end
	scores = zeros(Float64, edges)
	if work !== nothing && work.graph === graph
		copyto!(scores, 1, work.score, 1, min(edges, length(work.score)))
	end
	policy.matching = InterLinkMatchingWorkspace(
		graph, scores, sizehint!(Int[], edges), sizehint!(Rational{BigInt}[], edges), zeros(BigInt, edges),
		zeros(Int, vertices), sizehint!(Int[], vertices), sizehint!(Int[], vertices),
		zeros(BigInt, vertices), zeros(BigInt, vertices), zeros(Int, vertices),
		fill(false, vertices), fill(false, vertices), term_count == 0 ? nothing : zeros(term_count, edges))
	return policy.matching
end

# 2.2. Score every candidate link with the scheduler's policy
"""
	score_candidates!(model::InterlinkModel, policy::SchedulingPolicyModel, state, mu::Real;
	                  candidate_force=nothing, spacecraft=nothing, context=LinkContext())

Score every candidate link of `model` with the scheduler's policy
(`policy.objective`) into its matching workspace: `candidate_score` for each
available candidate, at the run state `state`, and 0 for the others. `mu` is the
planet's gravitational parameter (m³/s²), `spacecraft` the run's spacecraft
models (the power and data policies score nothing without them) and `context`
the link context. The orbital policies take each link's force from its
`interlink_wrench` (by default its `link_force`), or from a caller's
`candidate_force(target, partner, link)`. The target is a spacecraft run-state
view. A spacecraft partner is also a run-state view; a ground-station partner
is a `StateSample` with its station model, inertial position and velocity and
negative station index. A callback that handles station links must accept
that sample rather than assume both arguments are spacecraft state views.
"""
function score_candidates!(model::InterlinkModel, policy::SchedulingPolicyModel, state, mu::Real;
	candidate_force=nothing, spacecraft=nothing, context=LinkContext(), is_active=nothing, previous_active=nothing)
	ctx = SchedulingContext(Float64(mu), context, spacecraft, state, candidate_force)
	objective = policy.objective
	if spacecraft !== nothing && _uses_backlog_outflow(objective) &&
		any((previous_active === nothing ? model.graph.active : previous_active) .& model.graph.available)
		outflows = _data_store_outflows_bps(model, spacecraft, state, context; is_active, active_edges=previous_active)
		ctx = SchedulingContext(ctx.mu, context, spacecraft, state, candidate_force, ctx._terminals, outflows)
	end
	_prepare_scoring(objective, model, ctx)
	work = matching_workspace!(model, policy)
	_score_and_check!(work.score, objective, model, ctx; terms=work.terms)
	return nothing
end

_uses_backlog_outflow(objective) = false
_uses_backlog_outflow(::DataBacklogPolicy) = true
_uses_backlog_outflow(policy::WeightedPolicy) = any(term -> _uses_backlog_outflow(term.policy), policy.terms)

# Score every edge, then make the policy's checks across the scores.
function _score_and_check!(scores, objective, model, ctx::SchedulingContext; terms=nothing)
	_score_candidates!(scores, objective, model, ctx; terms)
	return _check_scores(objective, model, scores, ctx)
end

# The scores of every edge: one loop compiled for each policy type. With
# `circular`, the eccentricity policy also records its near-circular candidates.
function _score_candidates!(scores, objective::P, model::InterlinkModel, ctx::SchedulingContext, circular=nothing;
	terms=nothing) where {P}
	graph = model.graph
	t = ctx.link_context.time_s
	for edge in eachindex(graph.edges)
		score = 0.0
		terms === nothing || (terms[:, edge] .= 0.0)
		if graph.available[edge] && _may_score(objective, model, edge, ctx)
			candidate = _link_candidate(model, edge, ctx)
			score = Float64(circular !== nothing ? _eccentricity_score(objective, candidate, ctx, circular, edge) :
				terms === nothing ? candidate_score(objective, candidate, ctx, t) :
				_weighted_score(objective, candidate, ctx, t, terms, edge))
			isfinite(score) || throw(ArgumentError("Nonfinite $(nameof(P)) score for candidate $(candidate.key)."))
		end
		scores[edge] = score
	end
	return scores
end

# 3. Dynamic programming:  Select the best set of interlinks based on the current scores
# Helper Function: given one vertex and a related edge, this helper function of find the other vertex of the edge
other_endpoint(model, edge, vertex) = model.graph.edges[edge][1] == vertex ?
	model.graph.edges[edge][2] : model.graph.edges[edge][1]

# Helper function: select the best set of interlinks in a tree structure
function select_tree!(model, work, first_index, last_index, selected)
	for index in last_index:-1:first_index
		vertex = work.order[index]
		base = BigInt(0)
		for edge in model.graph.vertex_neighbors[vertex]
			work.weights[edge] > 0 || continue
			child = other_endpoint(model, edge, vertex)
			work.parent[child] == vertex || continue
			base += work.free_score[child]
		end
		work.blocked_score[vertex] = base
		gain, chosen = BigInt(0), 0
		for edge in model.graph.vertex_neighbors[vertex]
			work.weights[edge] > 0 || continue
			child = other_endpoint(model, edge, vertex)
			work.parent[child] == vertex || continue
			candidate = work.weights[edge] + work.blocked_score[child] - work.free_score[child]
			if candidate > gain
				gain, chosen = candidate, edge
			end
		end
		work.free_score[vertex] = base + gain
		work.choice[vertex] = chosen
	end
	for index in first_index:last_index
		vertex = work.order[index]
		edge = work.choice[vertex]
		(work.blocked[vertex] || edge == 0) && continue
		work.blocked[other_endpoint(model, edge, vertex)] = true
		push!(selected, map(v -> model.graph.vertices[v], model.graph.edges[edge]))
	end
	return nothing
end


function _select_cyclic!(model, work, selected)
	extension = Base.get_extension(parentmodule(parentmodule(@__MODULE__)), :SpaceAGORAPythonCallExt)
	extension === nothing && throw(ArgumentError(
		"Cyclic interlink matching requires the optional PythonCall backend. " *
		"Install PythonCall in the active Julia environment and load it with `using PythonCall`; " *
		"install NetworkX in the Python environment used by PythonCall. Tree matching works without Python."))
	return extension.select_cyclic!(model, work, selected)
end


"""
	select_interlinks!(model::InterlinkModel, policy::SchedulingPolicyModel)

Activate a maximum-total-score matching of available, positive-score links,
with at most one link per terminal. Of parallel links (directed links both ways
between the same two terminals), only the highest-scoring is a candidate, the
first in endpoint order on a tie. Topology belongs to `model.graph`; `policy.matching`
contains raw per-edge `score` values and reusable solver scratch arrays. Pass the same
policy used by `score_candidates!`; selection uses its stored scores, not its scoring
objective. A fresh workspace starts with zero scores. Appending links preserves existing
scores and initializes new scores to zero when the workspace is refreshed. Switching
graphs or discarding the workspace resets scores; rescore before selecting again.
Use separate policies for independent scheduling runs. Deep-copy the model and policy
together to preserve their workspace binding. Edge order is sorted locally on each
selection, not stored in either structure. Registered topology and iterative DFS
identify tree components, solved with linear-time dynamic programming. Cyclic components
use NetworkX's weighted Blossom algorithm through PythonCall, with O(N^3)
arithmetic operations for N cyclic endpoints. Both solvers use exact BigInt
weights scaled from Float64 scores; no floating-point score sums are compared.
Ties are deterministic for a given backend version, but tree ties may choose a
different optimum than Blossom alone. Register topology through `register_allowed_link!`.
For manual scores, refresh with `matching_workspace!(model, policy)` before editing
its `score` vector. Graph availability and workspace scores are reread on every call.
Tree components work without Python. Cyclic components require the optional
PythonCall package loaded with `using PythonCall`, and NetworkX installed in the
Python environment used by PythonCall. Core loading does not install either backend.
"""
function select_interlinks!(model::InterlinkModel, policy::SchedulingPolicyModel)
	# Step 1. Collect available, positive-score edges in canonical endpoint order.
	work = matching_workspace!(model, policy)
	empty!(work.candidates)
	empty!(work.exact_scores)
	edge_order = sortperm([map(v -> model.graph.vertices[v], endpoints) for endpoints in model.graph.edges])
	for edge in edge_order
		if model.graph.available[edge] && work.score[edge] > 0.0
			isfinite(work.score[edge]) || throw(ArgumentError("Interlink matching scores must be finite."))
			push!(work.candidates, edge)
			push!(work.exact_scores, Rational{BigInt}(work.score[edge]))
		end
	end
	_drop_parallel_candidates!(model, work)
	selected = InterLinkKey[]
	if !isempty(work.candidates)
		# Step 2. Scale exact Float64 rationals to integer weights without rounding score sums.
		# Their denominators are powers of two, so the largest is a common denominator.
		common_denominator = maximum(denominator, work.exact_scores)
		fill!(work.weights, 0)
		for (edge, score) in zip(work.candidates, work.exact_scores)
			work.weights[edge] = numerator(score) * div(common_denominator, denominator(score))
		end
		# Step 3. Traverse positive-weight components with depth-first search; solve trees and mark cycles if detected.
		fill!(work.parent, 0)
		fill!(work.blocked, false)
		fill!(work.cyclic, false)
		empty!(work.order)
		empty!(work.stack)
		for root in sortperm(model.graph.vertices)
			work.parent[root] == 0 || continue
			work.parent[root] = -1
			first_index = length(work.order) + 1
			degree_sum = 0
			push!(work.stack, root)
			while !isempty(work.stack)
				vertex = pop!(work.stack)
				push!(work.order, vertex)
				for edge in model.graph.vertex_neighbors[vertex]
					work.weights[edge] > 0 || continue
					degree_sum += 1
					partner = other_endpoint(model, edge, vertex)
					work.parent[partner] == 0 || continue
					work.parent[partner] = vertex
					push!(work.stack, partner)
				end
			end
			last_index = length(work.order)
			degree_sum == 0 && continue
			# Each undirected edge is counted twice; a connected tree has E = V - 1.
			if degree_sum ÷ 2 == last_index - first_index
				select_tree!(model, work, first_index, last_index, selected)
			else
				for index in first_index:last_index
					work.cyclic[work.order[index]] = true
				end
			end
		end
	end
	# Step 4. Solve cyclic components with NetworkX's maximum-weight Blossom matching.
	if !isempty(work.candidates) && any(work.cyclic)
		_select_cyclic!(model, work, selected)
	end
	# Step 5. Replace the active-link flags with the combined tree and cyclic selections.
	fill!(model.graph.active, false)
	for key in selected
		model.graph.active[model.graph.edge_index_lookup[key]] = true
	end
	return sort!(selected)
end

# Parallel candidates (directed links both ways between the same two terminals)
# share both terminals, so a matching holds at most one: keep the highest-scoring,
# the first in endpoint order on a tie. The components then stay simple graphs,
# and a two-way pair alone stays a tree.
function _drop_parallel_candidates!(model, work)
	length(work.candidates) < 2 && return nothing
	edges = model.graph.edges
	best = Dict{Tuple{Int, Int}, Int}()   # vertex pair => position in the candidates
	keep = fill(true, length(work.candidates))
	for (position, edge) in enumerate(work.candidates)
		pair = minmax(edges[edge]...)
		other = get(best, pair, 0)
		if other == 0
			best[pair] = position
		elseif work.exact_scores[position] > work.exact_scores[other]
			keep[other] = false
			best[pair] = position
		else
			keep[position] = false
		end
	end
	all(keep) && return nothing
	keepat!(work.candidates, keep)
	keepat!(work.exact_scores, keep)
	return nothing
end

# 4. Overall scheduling function to schedule interlinks based on the current state and scoring policy
function schedule_interlinks!(model::InterlinkModel, policy::SchedulingPolicyModel,
	spacecraft, state, mu::Real; candidate_force=nothing, is_active=nothing, context=LinkContext(),
	terminal_available=nothing, previous_active=nothing)
	# step 1. update the availability of each interlink based on the current state
	update_availability!(model, spacecraft, state; is_active, context, terminal_available)

	# step 2. compute the score of each interlink based on the scheduling policy / scoring policy
	score_candidates!(model, policy, state, mu; candidate_force, spacecraft, context, is_active, previous_active)

	# step 3. select the best set of interlinks based on the current scores
	return select_interlinks!(model, policy)
end

end # module SchedulingPolicies
