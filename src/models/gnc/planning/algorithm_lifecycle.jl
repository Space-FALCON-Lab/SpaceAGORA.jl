# The generic simulation binding of the GNC algorithm contract (gnc_algorithms.jl):
# a guidance model that runs one algorithm of any family through the shared
# request driver. The family supplies the request, the install target and its
# configuration checks; nothing here is specific to a family.
module AlgorithmLifecycle

using ..SimulationModel
import ..GNCAlgorithms as GA
import ..SimulationLifecycle as L
import ..GNCRuntime as GR
const S = SimulationModel

export AlgorithmEvent, AlgorithmGuidance, AlgorithmGuidanceError, algorithm_run_report

# The versioned namespace of the binding's random streams (`gnc_stream_words`).
const STREAM_NAMESPACE = "AlgorithmGuidance/v1"

"""
    AlgorithmEvent(time_s; reason=:forced, token=1)

A request that an `AlgorithmGuidance` makes at its first tick at or after
`time_s` (elapsed seconds, finite and nonnegative), with `reason` handed to the
family's request. A positive `token` identifies the event, which is delivered
once; events with distinct tokens due at one tick are separate requests, in the
order they are listed.
"""
struct AlgorithmEvent
    time_s::Float64
    reason::Symbol
    token::Int
    function AlgorithmEvent(time_s::Real; reason::Symbol=:forced, token::Integer=1)
        (isfinite(time_s) && time_s >= 0) || throw(ArgumentError("AlgorithmEvent time must be finite and nonnegative."))
        token > 0 || throw(ArgumentError("AlgorithmEvent token must be positive."))
        new(Float64(time_s), reason, Int(token))
    end
end

"""
    AlgorithmGuidanceError <: GNCAlgorithmError

A request of an `AlgorithmGuidance` failed: its family's rejection policy is
`:fail`, or algorithm or family code threw. Fields: `spacecraft_index`,
`request_id`, `reason`, `report` (a copy of the binding's run report, as in
`algorithm_run_report`), and the original `cause` and `backtrace`.
"""
struct AlgorithmGuidanceError <: GA.GNCAlgorithmError
    spacecraft_index::Int
    request_id::Int
    reason::Symbol
    report::Any
    cause::Any
    backtrace::Any
end
function Base.showerror(io::IO, e::AlgorithmGuidanceError)
    print(io, "Algorithm request failed for spacecraft ", e.spacecraft_index,
        ", request ", e.request_id, ": ", e.reason)
    e.cause === nothing || (print(io, ". Cause: "); showerror(io, e.cause))
end

"""
    AlgorithmGuidance(; algorithm, spacecraft_index, settings=nothing, link=nothing,
                      seed=0, events=AlgorithmEvent[], replan_interval_s=Inf)

The guidance model that runs `algorithm`, of any family, for the spacecraft at
`spacecraft_index` (its position in the run, from 1). Each request goes through
the shared request driver, `run_algorithm_request!`: the family builds the
owned request (`family_request`), the algorithm proposes (`propose!`), the
family validates the proposal (`validate_proposal`) and installs an accepted
product into the binding's install target (`family_link`, `install!`), and a
refusal follows the family's `on_rejection`.

- `settings`: the family's request settings, or `nothing` when it needs none.
- `link`: the install target, shared with the product's consumer, for a family
  that does not find its own in the run. The binding declares it as a GNC link it
  produces (`gnc_links`), so the engine checks and initializes it with its
  consumers (joint-reference family: the arm controller's plan buffer).
- `seed`: the seed of the binding's random stream,
  `MersenneTwister(gnc_stream_words("AlgorithmGuidance/v1", seed, F, id))`, with
  `F` the family's type name and `id` the spacecraft's ID; fresh in every run.
- `events`: requests at given times (`AlgorithmEvent`). An event due when the
  run starts is requested before the first tick.
- `replan_interval_s`: when finite, a request with reason `:interval` at the
  run's start and then whenever this many seconds have passed since the last
  request. The default `Inf` makes none.

Requests happen at the binding's guidance ticks, at the rate its
`GuidanceModel` gives it, and before the same tick's control, so the consumer
sees a product at the tick it is installed. The engine also runs every
guidance model outside those ticks when the burn scheduler plans (as a
spacecraft leaves the atmosphere); a request that is due by then is made
there. The binding keeps a record of every
request, installation, rejection and failure (`algorithm_run_report`). A
configuration may hold one binding per install target: one per family and
spacecraft for a family that finds its own target, one per configured `link`
otherwise. Checkpoint writing and resume are refused, because the binding's run
state is not checkpointed. On a spacecraft, `bind_spacecraft` gives a copy for
that spacecraft that shares the configuration and has fresh run state.
"""
Base.@kwdef mutable struct AlgorithmGuidance{A<:GA.AbstractGNCAlgorithm, T} <: S.AbstractGuidanceModel
    algorithm::A
    spacecraft_index::Int
    settings::T = nothing
    link::Any = nothing
    seed::Int = 0
    events::Vector{AlgorithmEvent} = AlgorithmEvent[]
    replan_interval_s::Float64 = Inf
    # Run state, set when each run starts.
    target::Any = nothing
    last_request_s::Float64 = -Inf
    runtime::Union{Nothing, GA.AlgorithmRuntime} = nothing
end

# A binding acts on its one spacecraft (GNCRuntime scope), and produces its
# configured link, which the product's consumer reads.
GR.spacecraft_scope(g::AlgorithmGuidance) = GR.BoundSpacecraft((g.spacecraft_index,))
GR.gnc_links(g::AlgorithmGuidance) = g.link === nothing ? () : (GR.Produces(g.link),)

# Guidance ticks accumulate their period, so an event or interval due at a tick
# may fall a few rounding units after it.
_time_atol(t) = 8 * eps(max(abs(Float64(t)), 1.0))

function _binding_report(g::AlgorithmGuidance)
    r = g.runtime
    return deepcopy((spacecraft_index=g.spacecraft_index, algorithm=string(typeof(g.algorithm)),
        settings=g.settings, seed=g.seed, generator="MersenneTwister",
        stream_mapping="SHA256/SpaceAGORA/$(STREAM_NAMESPACE)/little-endian",
        stream_words=r === nothing ? nothing : r.stream_words,
        request_id=r === nothing ? 0 : r.request_id,
        installed_request=r === nothing ? nothing : r.installed_request,
        product=r === nothing ? nothing : r.product,
        records=r === nothing ? NamedTuple[] : r.records))
end

function _fail!(g::AlgorithmGuidance, reason; cause=nothing, backtrace=nothing, details=NamedTuple())
    r = g.runtime
    r === nothing || push!(r.records, (event=:failure, request_id=r.request_id,
        reason=reason, details=deepcopy(details)))
    throw(AlgorithmGuidanceError(g.spacecraft_index, r === nothing ? 0 : r.request_id,
        reason, _binding_report(g), cause, backtrace))
end

# The binding functions of the request driver.
GA.algorithm_runtime(g::AlgorithmGuidance) = g.runtime
GA.bound_algorithm(g::AlgorithmGuidance) = g.algorithm
GA.bound_link(g::AlgorithmGuidance) = g.target
GA.build_algorithm_request(g::AlgorithmGuidance, u, p, t, reason, request_id) =
    GA.family_request(GA.algorithm_family(g.algorithm), g.settings, u, p, t, g.spacecraft_index, reason, request_id)
GA.fail_algorithm_request!(g::AlgorithmGuidance, reason; kwargs...) = _fail!(g, reason; kwargs...)

# Two bindings would install into one object: the same configured link, or, for
# a family that finds its own target, the same family and spacecraft.
_same_target(a::AlgorithmGuidance, b::AlgorithmGuidance) = a.link === nothing ?
    b.link === nothing && a.spacecraft_index == b.spacecraft_index &&
        typeof(GA.algorithm_family(a.algorithm)) === typeof(GA.algorithm_family(b.algorithm)) :
    a.link === b.link

function L.preflight_guidance(g::AlgorithmGuidance, args; isolate_state)
    n = length(args.dynamics_model.spacecraft)
    1 <= g.spacecraft_index <= n || throw(ArgumentError(
        "AlgorithmGuidance spacecraft_index $(g.spacecraft_index) is not a spacecraft of the run (1:$(n))."))
    g.replan_interval_s > 0 || throw(ArgumentError("AlgorithmGuidance replan_interval_s must be positive."))
    allunique(e.token for e in g.events) || throw(ArgumentError("AlgorithmGuidance event tokens must be distinct."))
    ss = args.simulation_settings
    (ss.checkpoint_enabled || ss.resume_from_checkpoint) && throw(ArgumentError(
        "AlgorithmGuidance does not support checkpoint writing or resume: its run state is not checkpointed."))
    count(x -> x isa AlgorithmGuidance && _same_target(x, g), args.guidance_model.guidance_effectors) == 1 ||
        throw(ArgumentError("Two AlgorithmGuidance bindings install into one target (spacecraft $(g.spacecraft_index)); keep one."))
    GA.require_algorithm_support(g.algorithm)
    GA.preflight_family(GA.algorithm_family(g.algorithm), g.algorithm, g.settings, args, g.spacecraft_index;
        isolate_state=isolate_state)
    return nothing
end

function L.initialize_guidance!(g::AlgorithmGuidance, u, p, t)
    t = Float64(t)
    g.runtime, g.target, g.last_request_s = nothing, nothing, -Inf
    try
        family = GA.algorithm_family(g.algorithm)
        id = p.args.dynamics_model.spacecraft[g.spacecraft_index].id
        words = GA.gnc_stream_words(STREAM_NAMESPACE, g.seed, nameof(typeof(family)), id)
        g.runtime = GA.AlgorithmRuntime(words)
        g.target = GA.family_link(family, g.link, g.settings, p, g.spacecraft_index)
        g.runtime.state = GA.initialize_algorithm(g.algorithm,
            (seed=g.seed, spacecraft_index=g.spacecraft_index, spacecraft_id=id,
             epoch=deepcopy(p.args.initial_time), stream_words=copy(words)))
        # A tick at this same time (the burn scheduler's pass at the run's start)
        # then makes no second request.
        g.runtime.last_guidance_time = t
        _request_due!(g, u, p, t)
    catch e
        e isa AlgorithmGuidanceError && rethrow()
        _fail!(g, :initialization_exception; cause=e, backtrace=catch_backtrace())
    end
    return nothing
end

_due(e::AlgorithmEvent, r, t) = e.time_s <= t + _time_atol(t) && !(e.token in r.processed_tokens)

function _request!(g::AlgorithmGuidance, u, p, t, reason)
    g.last_request_s = t
    GA.run_algorithm_request!(g, u, p, t, reason)
    return nothing
end

function _request_due!(g::AlgorithmGuidance, u, p, t)
    r = g.runtime
    requested = false
    for e in g.events
        _due(e, r, t) || continue
        push!(r.processed_tokens, e.token)
        _request!(g, u, p, t, e.reason)
        requested = true
    end
    if !requested && isfinite(g.replan_interval_s) &&
            t - g.last_request_s >= g.replan_interval_s - _time_atol(t)
        _request!(g, u, p, t, :interval)
    end
    return nothing
end

function S.GuidanceHooks.update_guidance!(g::AlgorithmGuidance, u, p, t::Float64, spacecraft_index::Int)
    spacecraft_index == g.spacecraft_index || return nothing
    r = g.runtime
    r === nothing && _fail!(g, :uninitialized)
    # One tick can arrive through the periodic callback and an event-driven
    # controller's scheduling pass; only a due event makes a second request.
    t == r.last_guidance_time && !any(e -> _due(e, r, t), g.events) && return nothing
    t >= r.last_guidance_time || _fail!(g, :nonmonotonic_guidance_time)
    r.last_guidance_time = t
    _request_due!(g, u, p, t)
    return nothing
end

# A binding installs its product before the same tick's control, which consumes it.
S.SimulationCallbacks._guidance_precedes_control(::AlgorithmGuidance) = true

L.bind_spacecraft(g::AlgorithmGuidance, sat_idx::Int) = AlgorithmGuidance(;
    algorithm=g.algorithm, spacecraft_index=sat_idx, settings=g.settings, link=g.link,
    seed=g.seed, events=g.events, replan_interval_s=g.replan_interval_s)

"""
    algorithm_run_report(solution) -> NamedTuple

Copied run records of every `AlgorithmGuidance` in an isolated `run_simulation`
result (`return_solution=true`), as `bindings`, one entry per binding: its
`spacecraft_index`, `algorithm` type, `settings`, `seed`, random-stream identity
(`generator`, `stream_mapping`, `stream_words`), `request_id` (the number of
requests made), the `installed_request` and `product`, and `records`: one per
request (`:request`, with a copy of the request), installation (`:installed`),
rejection (`:rejected`, with the family's validation result and reason) and
failure (`:failure`), in order.
"""
function algorithm_run_report(solution)
    args = solution.prob.p.args
    return (bindings=[_binding_report(g) for g in args.guidance_model.guidance_effectors if g isa AlgorithmGuidance],)
end

end # module AlgorithmLifecycle
