# Algorithms of the general GNC algorithm contract (`GNCAlgorithms`) and the families
# they declare. This file needs only `GNCAlgorithms`, `reports.jl`, `walker.jl` and
# standard libraries, so the standalone contract tests include it as they include
# the contract itself.

const _GA = GNCAlgorithms
const _ALGORITHM_STATUSES = (:candidate, :infeasible, :unsupported, :failed)
const _REQUEST_DRIVER_FAILURES = (:invalid_result_type, :algorithm_exception)

"""
    check_algorithm(algorithm, requests; trusted=NamedTuple(), run_context=(seed=1,), simulation=nothing)
        -> ConformanceReport

Check a planning algorithm of the general GNC algorithm contract against its
declared family, on `requests`, the family's representative requests (for RPO,
[`rpo_conformance_requests`](@ref)), proposed in order from a fresh state.
`trusted` holds the queries the family's validation needs (for RPO,
[`rpo_trusted_queries`](@ref)); `run_context` is passed to `initialize_algorithm`
with the stream's seed words added. A request the family's `is_revision_request`
marks is sent to `propose_revision!` with the product of the latest accepted
proposal. The checks, in report order:

- `family_resolves`: `algorithm_family` returns an `AbstractAlgorithmFamily`.
- `family_trusted`: the family is defined in core (a family module or a service
  facade), not in the algorithm's own package.
- `result_contract`: every proposal has the family's `algorithm_result_type`, a
  known status, the request's `request_id`, and a product (`proposed_product`)
  exactly when its status is `:candidate`.
- `stream_untouched_if_declared`: an algorithm that declares `consumes_rng = false`
  leaves the stream unchanged. Skipped for one that declares `true`.
- `deterministic_given_stream`: two walks through the requests from fresh states and
  equal streams, with the global random number generator seeded 1 and 2, give equal
  results.
- `request_not_mutated`: the copy of each request handed to the algorithm is
  unchanged after its proposal.
- `initialize_independent`: two states from `initialize_algorithm` share no mutable
  object.
- `refuses_explicitly`: a revision request the family's capability check refuses
  gets an `:unsupported` result from `propose_revision!`, never a product. Skipped
  when no request is such a revision.
- `driver_round_trip`: the requests run in order through `run_algorithm_request!`
  with the kit's test binding end as the family decides on proposals made from
  copies of the same state and stream: an accepted proposal is installed (the link
  changes, the installed product equals the proposal and is current at the
  request's time), a refusal follows the family's `on_rejection`. No request may end
  as `:invalid_result_type` or `:algorithm_exception`. Skipped for a family without
  a conformance link.
- `simulation_run`: `simulation(algorithm)`, a short simulation that returns the
  binding's run records (for RPO, [`rpo_conformance_simulation`](@ref)), completes
  with a product installed and no failure record. Skipped when `simulation` is
  `nothing`.

A pass certifies interface behaviour, not the algorithm's method.
"""
check_algorithm(algorithm, requests; trusted::NamedTuple=NamedTuple(), run_context::NamedTuple=(seed=1,),
    simulation=nothing) = _preserving_global_rng(() -> _check_algorithm(algorithm, requests, trusted, run_context, simulation))

function _check_algorithm(algorithm, requests, trusted, run_context, simulation)
    since = FIRST_RELEASE
    requests = collect(Any, requests)
    target = string(nameof(typeof(algorithm)))
    family = try
        _GA.algorithm_family(algorithm)
    catch err
        err
    end
    if !(family isa _GA.AbstractAlgorithmFamily)
        why = family isa Exception ? "algorithm_family threw: " * _summary(family) :
            "algorithm_family returned $(repr(family)), not an AbstractAlgorithmFamily"
        checks = [i == 1 ? _check(id, since, false, why) :
            _check(id, since, false, "the algorithm's family did not resolve (family_resolves)")
            for (i, id) in enumerate(_ALGORITHM_CHECK_IDS)]
        return ConformanceReport(target, :algorithm, checks)
    end
    walks = map(seed -> _with_global_seed(() -> _walk(family, algorithm, requests, trusted, run_context), seed), (1, 2))
    checks = [
        _check(:family_resolves, since, true),
        _family_trusted_check(family, since),
        _result_contract_check(family, walks[1], since),
        _stream_check(algorithm, walks[1], since),
        _deterministic_check(walks, since),
        _request_check(walks[1], since),
        _initialize_check(algorithm, run_context, since),
        _refusal_check(family, algorithm, walks[1], run_context, since),
        _driver_check(family, algorithm, requests, trusted, run_context, conformance_link(family), since),
        _simulation_check(algorithm, simulation, since),
    ]
    return ConformanceReport(target, :algorithm, checks)
end

const _ALGORITHM_CHECK_IDS = (:family_resolves, :family_trusted, :result_contract, :stream_untouched_if_declared,
    :deterministic_given_stream, :request_not_mutated, :initialize_independent, :refuses_explicitly,
    :driver_round_trip, :simulation_run)
const _ALGORITHM_FAMILY_CHECK_IDS = (:family_resolves, :family_obligations, :driver_round_trip)

"""
    check_algorithm_family(family, algorithm, requests; link, trusted=NamedTuple(), run_context=(seed=1,))
        -> ConformanceReport

Check a family of the general GNC algorithm contract (a new family's core half, or
a family a service facade serves), with an algorithm of the family and its
representative `requests`, which should include one the family refuses. `link` is
a fresh object of the kind products are installed into; each run uses a copy. The
checks, in report order:

- `family_resolves`: `algorithm_family(algorithm)` has the type of `family`, as the
  request driver requires.
- `family_obligations`: `validate_proposal`, `install!` and `is_current` have
  methods for the family's requests and products; `fallback_product` has one when
  `on_rejection(family)` is `:install_fallback`; `on_rejection` is a known policy;
  validation leaves the request, the result, the trusted queries and the global
  random number generator unchanged; `install!` changes none of its arguments but
  the link; `is_current` returns a `Bool`; `algorithm_result_type` is a subtype of
  `AbstractAlgorithmResult`; and a family a service facade serves (its module
  defines `CONTRACT_VERSION`) states `family_contract_version`.
- `driver_round_trip`: as for [`check_algorithm`](@ref), with `family` as the bound
  family, so its rejection policy is the one exercised.
"""
check_algorithm_family(family::_GA.AbstractAlgorithmFamily, algorithm, requests; link,
    trusted::NamedTuple=NamedTuple(), run_context::NamedTuple=(seed=1,)) =
    _preserving_global_rng(() -> _check_algorithm_family(family, algorithm, requests, link, trusted, run_context))

function _check_algorithm_family(family, algorithm, requests, link, trusted, run_context)
    since = FIRST_RELEASE
    requests = collect(Any, requests)
    declared = try
        _GA.algorithm_family(algorithm)
    catch err
        err
    end
    resolves = if declared isa Exception
        _check(:family_resolves, since, false, "algorithm_family threw: " * _summary(declared))
    else
        _check(:family_resolves, since, Base.typename(typeof(declared)) === Base.typename(typeof(family)),
            "the algorithm declares a $(typeof(declared)), not a $(typeof(family)); the request driver refuses " *
            "the pair (:family_mismatch)")
    end
    checks = [
        resolves,
        _obligations_check(family, algorithm, requests, trusted, link, run_context, since),
        _driver_check(family, algorithm, requests, trusted, run_context, link, since; copy_link=true),
    ]
    return ConformanceReport(string(nameof(typeof(family))), :algorithm_family, checks)
end

"""
    test_algorithm(algorithm, requests; kwargs...)

Run [`check_algorithm`](@ref) and record its report with [`test_report`](@ref).
Needs the `Test` standard library: `using Test` loads the method.
"""
function test_algorithm end

"""
    test_algorithm_family(family, algorithm, requests; kwargs...)

Run [`check_algorithm_family`](@ref) and record its report with [`test_report`](@ref).
Needs the `Test` standard library: `using Test` loads the method.
"""
function test_algorithm_family end

# A fresh object of the kind `family`'s products are installed into, for the kit's
# test binding (for RPO, an empty `RPOPlanBuffer`). A core family defines a method
# next to its conformance requests; the fallback `nothing` makes `check_algorithm`
# skip `driver_round_trip`.
conformance_link(family) = nothing

# ---- Walking the requests ----------------------------------------------------

_request_time(request) = hasproperty(request, :time_s) ? Float64(getproperty(request, :time_s)) : 0.0
_request_reason(request) = hasproperty(request, :reason) ? getproperty(request, :reason) : :initial
_stream_words(run_context) = _GA.gnc_stream_words("Conformance/v1", get(run_context, :seed, 1))
_initialize(algorithm, run_context) =
    _GA.initialize_algorithm(algorithm, deepcopy(merge(run_context, (stream_words=_stream_words(run_context),))))

# One proposal: what the driver would ask for, on a copy of the request.
function _propose(family, state, algorithm, request, active, rng)
    given = deepcopy(request)
    if _GA.is_revision_request(family, request)
        active === nothing && return (given=given, result=nothing, error=nothing, revision=true, skipped=true)
        result = _GA.propose_revision!(state, algorithm, given, deepcopy(active), rng)
    else
        result = _GA.propose!(state, algorithm, given, rng)
    end
    return (given=given, result=result, error=nothing, revision=_GA.is_revision_request(family, request), skipped=false)
end

# The requests in order, from a fresh state and stream, as the driver would see
# them: a request the capability check refuses is not proposed; an accepted
# proposal becomes the active product of later revision requests.
function _walk(family, algorithm, requests, trusted, run_context)
    entries = NamedTuple[]
    result_type = _result_type(family)
    result_type isa Type || return (entries=entries, failure=result_type, rng=nothing, state=nothing)
    state = try
        _initialize(algorithm, run_context)
    catch err
        return (entries=entries, failure="initialize_algorithm threw: " * _summary(err), rng=nothing, state=nothing)
    end
    rng = Random.MersenneTwister(_stream_words(run_context))
    active = nothing
    for (k, request) in enumerate(requests)
        capabilities = try
            _GA.validate_capabilities(family, algorithm, request)
        catch err
            push!(entries, (index=k, request=request, accepted=false, capabilities=nothing, given=nothing,
                result=nothing, error=err, untouched=true, revision=false, validation=nothing))
            continue
        end
        if !capabilities.accepted
            push!(entries, (index=k, request=request, accepted=false, capabilities=capabilities, given=nothing,
                result=nothing, error=nothing, untouched=true, revision=_GA.is_revision_request(family, request),
                validation=nothing))
            continue
        end
        witness = copy(rng)
        proposal = try
            _propose(family, state, algorithm, request, active, rng)
        catch err
            (given=nothing, result=nothing, error=err, revision=false, skipped=false)
        end
        untouched = rand(copy(rng), UInt64) == rand(copy(witness), UInt64)
        proposal.skipped && continue
        validation = nothing
        if proposal.error === nothing && proposal.result isa result_type
            validation = try
                _GA.validate_proposal(family, deepcopy(request), proposal.result, trusted)
            catch err
                err
            end
            validation isa Exception || !validation.accepted || (active = _GA.proposed_product(proposal.result))
        end
        push!(entries, (index=k, request=request, accepted=true, capabilities=capabilities, given=proposal.given,
            result=proposal.result, error=proposal.error, untouched=untouched, revision=proposal.revision,
            validation=validation))
    end
    return (entries=entries, failure=nothing, rng=rng, state=state)
end

_walk_failure(walk) = walk.failure

# The family's result type, or why there is none: `algorithm_result_type` must return
# a subtype of `AbstractAlgorithmResult`, and every check that compares proposals
# with it stops when it does not.
function _result_type(family)
    T = try
        _GA.algorithm_result_type(family)
    catch err
        return "algorithm_result_type threw: " * _summary(err)
    end
    T isa Type && T <: _GA.AbstractAlgorithmResult && return T
    return "algorithm_result_type returns $(repr(T)), not a subtype of AbstractAlgorithmResult"
end

# ---- The algorithm checks ----------------------------------------------------

# The module tree that holds the contract (SpaceAGORA, or the standalone tests).
_contract_root() = parentmodule(_GA)

function _defined_under(m::Module, root::Module)
    while true
        m === root && return true
        parent = parentmodule(m)
        parent === m && return false
        m = parent
    end
end

function _family_trusted_check(family, since)
    home = parentmodule(typeof(family))
    return _check(:family_trusted, since, _defined_under(home, _contract_root()),
        "the family $(typeof(family)) is defined in $(home), outside $(_contract_root()); families are trusted and " *
        "core-owned: declare a family core defines (for RPO, any AbstractRPOPlanner), or add the family to core " *
        "with a review")
end

function _result_contract_check(family, walk, since)
    problem = _walk_failure(walk)
    problem === nothing || return _check(:result_contract, since, false, problem)
    expected = _result_type(family)
    problems = String[]
    for entry in walk.entries
        entry.accepted || continue
        k = entry.index
        if entry.error !== nothing
            push!(problems, "request $k: the proposal threw: " * _summary(entry.error))
            continue
        end
        result = entry.result
        if !(result isa expected)
            push!(problems, "request $k: the proposal is a $(typeof(result)), not the family's $(expected)")
            continue
        end
        status = hasproperty(result, :status) ? result.status : nothing
        status in _ALGORITHM_STATUSES ||
            push!(problems, "request $k: status $(repr(status)) is not one of $(_ALGORITHM_STATUSES)")
        if hasproperty(entry.request, :request_id)
            id = hasproperty(result, :request_id) ? result.request_id : nothing
            isequal(id, entry.request.request_id) ||
                push!(problems, "request $k: the result's request_id is $(repr(id)), the request's $(entry.request.request_id)")
        end
        product = try
            _GA.proposed_product(result)
        catch err
            push!(problems, "request $k: proposed_product threw: " * _summary(err))
            continue
        end
        (status === :candidate) == (product !== nothing) || push!(problems, status === :candidate ?
            "request $k: a :candidate result without a product" :
            "request $k: a $(repr(status)) result with a product; only a :candidate carries one")
    end
    return _check(:result_contract, since, isempty(problems), join(problems, "\n"))
end

function _stream_check(algorithm, walk, since)
    _GA.consumes_rng(algorithm) &&
        return _skip(:stream_untouched_if_declared, since, "the algorithm declares consumes_rng = true")
    problem = _walk_failure(walk)
    problem === nothing || return _check(:stream_untouched_if_declared, since, false, problem)
    moved = [entry.index for entry in walk.entries if !entry.untouched]
    return _check(:stream_untouched_if_declared, since, isempty(moved),
        "the algorithm declares consumes_rng = false, yet its proposal for request(s) $(join(moved, ", ")) " *
        "advanced the stream")
end

function _deterministic_check(walks, since)
    a, b = walks
    for walk in walks
        problem = _walk_failure(walk)
        problem === nothing || return _check(:deterministic_given_stream, since, false, problem)
    end
    length(a.entries) == length(b.entries) || return _check(:deterministic_given_stream, since, false,
        "the walks proposed for $(length(a.entries)) and $(length(b.entries)) requests")
    for (x, y) in zip(a.entries, b.entries)
        same = x.accepted == y.accepted && (x.error === nothing) == (y.error === nothing) &&
            (x.error === nothing ? _same_state(x.result, y.result) : typeof(x.error) === typeof(y.error))
        same || return _check(:deterministic_given_stream, since, false,
            "request $(x.index): equal streams and fresh states gave different results with the global random " *
            "number generator seeded 1 and 2; draw only from the stream passed to propose!, and keep run state in " *
            "the state initialize_algorithm returns")
    end
    _same_state(a.state, b.state) || return _check(:deterministic_given_stream, since, false,
        "the algorithm's state after the walks differs with the global random number generator seeded 1 and 2")
    return _check(:deterministic_given_stream, since, true)
end

function _request_check(walk, since)
    problem = _walk_failure(walk)
    problem === nothing || return _check(:request_not_mutated, since, false, problem)
    changed = [entry.index for entry in walk.entries
        if entry.given !== nothing && !_same_state(entry.given, entry.request)]
    return _check(:request_not_mutated, since, isempty(changed),
        "the algorithm changed its copy of request(s) $(join(changed, ", ")); treat a request as read-only")
end

function _initialize_check(algorithm, run_context, since)
    states = try
        (_initialize(algorithm, run_context), _initialize(algorithm, run_context))
    catch err
        return _check(:initialize_independent, since, false, "initialize_algorithm threw: " * _summary(err))
    end
    shared = _shared_mutables(states...)
    return _check(:initialize_independent, since, isempty(shared),
        "two states from initialize_algorithm share $(length(shared)) mutable object(s) (the first is a " *
        "$(isempty(shared) ? "" : typeof(first(shared)))); create every mutable part of the state in initialize_algorithm")
end

function _refusal_check(family, algorithm, walk, run_context, since)
    problem = _walk_failure(walk)
    problem === nothing || return _check(:refuses_explicitly, since, false, problem)
    refused = [entry for entry in walk.entries if !entry.accepted && entry.revision && entry.capabilities !== nothing]
    isempty(refused) && return _skip(:refuses_explicitly, since,
        "no revision request among the requests is refused by the family's capability check")
    active = nothing
    for entry in walk.entries
        entry.accepted && entry.validation isa _GA.AbstractValidationResult && entry.validation.accepted &&
            (active = _GA.proposed_product(entry.result))
    end
    active === nothing && return _skip(:refuses_explicitly, since,
        "no proposal was accepted, so there is no product to ask a revision of")
    problems = String[]
    for entry in refused
        result = try
            _GA.propose_revision!(_initialize(algorithm, run_context), algorithm, deepcopy(entry.request),
                deepcopy(active), Random.MersenneTwister(_stream_words(run_context)))
        catch err
            push!(problems, "request $(entry.index): propose_revision! threw: " * _summary(err))
            continue
        end
        status = hasproperty(result, :status) ? result.status : nothing
        product = try
            _GA.proposed_product(result)
        catch err
            push!(problems, "request $(entry.index): proposed_product threw on the refused revision's result: " *
                _summary(err))
            continue
        end
        status === :unsupported && product === nothing && continue
        returned = repr(status) * (product === nothing ? "" : " with a product")
        push!(problems, "request $(entry.index): the capability check refuses this revision " *
            "($(entry.capabilities.reason)), yet propose_revision! returned $(returned); return :unsupported with " *
            "no product, never a full proposal in its place")
    end
    return _check(:refuses_explicitly, since, isempty(problems), join(problems, "\n"))
end

# ---- The request driver with the kit's test binding ---------------------------

# The kit's binding: it hands the driver copies of the caller's requests, in order,
# and fails with `_TestBindingFailure`.
mutable struct _TestBinding
    algorithm::Any
    family::Any
    runtime::_GA.AlgorithmRuntime
    link::Any
    trusted::NamedTuple
    requests::Vector{Any}
end
struct _TestBindingFailure <: _GA.GNCAlgorithmError
    reason::Symbol
    cause::Any
    details::Any
end
_GA.algorithm_runtime(binding::_TestBinding) = binding.runtime
_GA.bound_algorithm(binding::_TestBinding) = binding.algorithm
_GA.bound_family(binding::_TestBinding) = binding.family
_GA.bound_link(binding::_TestBinding) = binding.link
_GA.trusted_queries(binding::_TestBinding) = binding.trusted
_GA.build_algorithm_request(binding::_TestBinding, u, p, t, reason, request_id) = deepcopy(binding.requests[u])
function _GA.fail_algorithm_request!(binding::_TestBinding, reason; cause=nothing, backtrace=nothing,
        details=NamedTuple())
    push!(binding.runtime.records, (event=:failure, request_id=binding.runtime.request_id, reason=reason))
    throw(_TestBindingFailure(reason, cause, details))
end

# What the family decides on one request, from copies of the binding's state and
# stream: (:refused, validation), (:accepted, product), or (:driver_failure, reason).
function _verdict(family, algorithm, request, state, rng, active, trusted)
    try
        capabilities = _GA.validate_capabilities(family, algorithm, request)
        capabilities.accepted || return (:refused, capabilities)
        result = _GA.is_revision_request(family, request) ?
            _GA.propose_revision!(state, algorithm, deepcopy(request), deepcopy(active), rng) :
            _GA.propose!(state, algorithm, deepcopy(request), rng)
        T = _result_type(family)
        T isa Type && result isa T || return (:driver_failure, :invalid_result_type)
        product = deepcopy(_GA.proposed_product(result))
        validation = _GA.validate_proposal(family, deepcopy(request), result, trusted)
        validation.accepted || return (:refused, validation)
        return (:accepted, product)
    catch
        return (:driver_failure, :algorithm_exception)
    end
end

function _driver_check(family, algorithm, requests, trusted, run_context, link, since; copy_link::Bool=false)
    link === nothing && return _skip(:driver_round_trip, since,
        "the family has no conformance link (conformance_link) for the kit's test binding")
    runtime = _GA.AlgorithmRuntime(_stream_words(run_context))
    binding = _TestBinding(algorithm, family, runtime, copy_link ? deepcopy(link) : link, trusted, requests)
    try
        runtime.state = _initialize(algorithm, run_context)
    catch err
        return _check(:driver_round_trip, since, false, "initialize_algorithm threw: " * _summary(err))
    end
    policy = _GA.on_rejection(family)
    problems = String[]
    for (k, request) in enumerate(requests)
        t = _request_time(request)
        if _GA.is_revision_request(family, request) && runtime.product === nothing
            continue   # the driver has no installed product to revise
        end
        kind, value = _verdict(family, algorithm, request, deepcopy(runtime.state), copy(runtime.rng), runtime.product,
            trusted)
        link_before = deepcopy(binding.link)
        outcome = try
            _GA.run_algorithm_request!(binding, k, nothing, t, _request_reason(request))
        catch err
            err
        end
        last_event = isempty(runtime.records) ? nothing : last(runtime.records).event
        if outcome isa _TestBindingFailure && outcome.reason in _REQUEST_DRIVER_FAILURES
            push!(problems, "request $k: the request driver failed with $(repr(outcome.reason))" *
                (outcome.cause === nothing ? "" : ": " * _summary(outcome.cause)))
        elseif outcome isa Exception && !(outcome isa _TestBindingFailure)
            push!(problems, "request $k: the request driver threw outside the binding's failure: " * _summary(outcome))
        elseif kind === :driver_failure
            push!(problems, "request $k: on copies of the binding's state and stream, the request ends as " *
                "$(repr(value)), yet the driver's run gave $(repr(outcome))")
        elseif kind === :accepted
            if outcome !== true
                push!(problems, "request $k: on copies of the binding's state and stream the family accepts the " *
                    "proposal, yet the driver's run gave $(_outcome_text(outcome)); the algorithm answers differently " *
                    "for a copy of its state or stream")
            elseif !_same_state(runtime.product, value)
                push!(problems, "request $k: the installed product differs from the proposal made on copies of the " *
                    "binding's state and stream")
            elseif _same_state(binding.link, link_before)
                push!(problems, "request $k: the proposal was accepted, yet install! left the link unchanged")
            elseif !_current(family, request, runtime.product, t)
                push!(problems, "request $k: the installed product is not current at the request's time " *
                    "(is_current returned false)")
            end
        else   # refused
            reason = value.reason
            if policy === :fail
                outcome isa _TestBindingFailure && outcome.reason === reason ||
                    push!(problems, "request $k: the family refuses it ($(repr(reason))) under :fail, yet the " *
                        "driver's run gave $(_outcome_text(outcome))")
                _same_state(binding.link, link_before) ||
                    push!(problems, "request $k: a refused request under :fail changed the link")
            elseif policy === :keep_previous
                outcome === false && last_event === :rejected && _same_state(binding.link, link_before) ||
                    push!(problems, "request $k: the family refuses it ($(repr(reason))) under :keep_previous, yet " *
                        "the driver's run gave $(_outcome_text(outcome)) or changed the link")
            elseif policy === :install_fallback
                fallback = try
                    deepcopy(_GA.fallback_product(family, deepcopy(request)))
                catch err
                    err
                end
                outcome === true && last_event === :fallback_installed && _same_state(runtime.product, fallback) ||
                    push!(problems, "request $k: the family refuses it ($(repr(reason))) under :install_fallback, " *
                        "yet the driver's run gave $(_outcome_text(outcome)) without installing fallback_product")
            end
        end
    end
    return _check(:driver_round_trip, since, isempty(problems), join(problems, "\n"))
end

_outcome_text(outcome) = outcome === true ? "an installed product" : outcome === false ? "a kept product" :
    outcome isa _TestBindingFailure ? "a failure with reason $(repr(outcome.reason))" : _summary(outcome)

function _current(family, request, product, t)
    try
        return _GA.is_current(family, request, product, t) === true
    catch
        return false
    end
end

function _simulation_check(algorithm, simulation, since)
    simulation === nothing && return _skip(:simulation_run, since, "no simulation was given")
    records = try
        simulation(deepcopy(algorithm))
    catch err
        return _check(:simulation_run, since, false, "the simulation failed: " * _summary(err))
    end
    events = [record.event for record in records if hasproperty(record, :event)]
    failures = [record for record in records if hasproperty(record, :event) && record.event === :failure]
    isempty(failures) || return _check(:simulation_run, since, false,
        "the run recorded a failure: $(first(failures))")
    return _check(:simulation_run, since, any(in((:installed, :fallback_installed)), events),
        "the run installed no product (events: $(join(unique(events), ", ")))")
end

# ---- The family obligations --------------------------------------------------

_defines_contract_version(m::Module) = isdefined(m, :CONTRACT_VERSION)

_global_draw_state() = copy(Random.default_rng())
_global_stream_unchanged(before) = rand(copy(before), UInt64) == rand(copy(Random.default_rng()), UInt64)

function _obligations_check(family, algorithm, requests, trusted, link, run_context, since)
    problems = String[]
    policy = _GA.on_rejection(family)
    policy in (:fail, :keep_previous, :install_fallback) ||
        push!(problems, "on_rejection returns $(repr(policy)), not :fail, :keep_previous or :install_fallback")
    result_type = _result_type(family)
    result_type isa Type || push!(problems, result_type)
    version = _GA.family_contract_version(family)
    version === nothing || version isa VersionNumber ||
        push!(problems, "family_contract_version returns $(repr(version)), not a VersionNumber or nothing")
    home = parentmodule(typeof(family))
    _defines_contract_version(home) && version === nothing &&
        push!(problems, "the family is served by a facade ($(home) defines CONTRACT_VERSION), so " *
            "family_contract_version must state the family's contract version")
    if policy === :install_fallback
        for request in requests
            hasmethod(_GA.fallback_product, Tuple{typeof(family), typeof(request)}) && continue
            push!(problems, "on_rejection is :install_fallback, yet fallback_product has no method for a " *
                "$(typeof(request))")
            break
        end
    end
    # Without a result type no proposal can be checked against the family.
    result_type isa Type || return _check(:family_obligations, since, false, join(unique(problems), "\n"))
    walk = _walk(family, algorithm, requests, trusted, run_context)
    walk.failure === nothing || push!(problems, walk.failure)
    candidates = [entry for entry in walk.entries if entry.accepted && entry.error === nothing &&
        entry.result isa result_type]
    isempty(candidates) && push!(problems, "no request gave a proposal of the family's result type to validate")
    accepted = false
    for entry in candidates
        k, request, result = entry.index, entry.request, entry.result
        # The hooks a validated product meets, whether or not this one is accepted.
        offered = try
            _GA.proposed_product(result)
        catch err
            push!(problems, "request $k: proposed_product threw: " * _summary(err))
            continue
        end
        if offered !== nothing
            P, R = typeof(offered), typeof(request)
            hasmethod(_GA.install!, Tuple{typeof(link), typeof(family), P, R, Float64}) || push!(problems,
                "install! has no method for a $(typeof(link)) link, this family, a $(P) product and a $(R)")
            hasmethod(_GA.is_current, Tuple{typeof(family), R, P, Float64}) ||
                push!(problems, "is_current has no method for this family, a $(R) and a $(P) product")
        end
        request_copy, result_copy, trusted_copy = deepcopy(request), deepcopy(result), deepcopy(trusted)
        # A validation that draws from the global generator is reported, and the
        # generator is restored, whether or not the validation returns.
        before = _global_draw_state()
        drew = Ref(false)
        validation = try
            _preserving_global_rng() do
                checked = _GA.validate_proposal(family, request_copy, result_copy, trusted_copy)
                drew[] = !_global_stream_unchanged(before)
                checked
            end
        catch err
            push!(problems, "request $k: validate_proposal threw: " * _summary(err))
            continue
        end
        drew[] && push!(problems, "request $k: validate_proposal drew from the global random number generator")
        _same_state(request_copy, request) && _same_state(result_copy, result) && _same_state(trusted_copy, trusted) ||
            push!(problems, "request $k: validate_proposal changed the request, the result or the trusted queries")
        validation isa _GA.AbstractValidationResult ||
            (push!(problems, "request $k: validate_proposal returned a $(typeof(validation)), not an " *
                "AbstractValidationResult"); continue)
        validation.accepted || continue
        accepted = true
        product = offered
        target = deepcopy(link)
        diagnostics = hasproperty(result, :diagnostics) ? result.diagnostics : NamedTuple()
        arguments = (deepcopy(product), deepcopy(request), deepcopy(validation), deepcopy(diagnostics))
        snapshot = deepcopy(arguments)
        try
            _GA.install!(target, family, arguments[1], arguments[2], _request_time(request);
                validation=arguments[3], diagnostics=arguments[4])
        catch err
            push!(problems, "request $k: install! threw: " * _summary(err))
            continue
        end
        _same_state(arguments, snapshot) ||
            push!(problems, "request $k: install! changed the product, request, validation or diagnostics it was given")
        current = try
            _GA.is_current(family, request, product, _request_time(request))
        catch err
            push!(problems, "request $k: is_current threw: " * _summary(err))
            continue
        end
        current isa Bool || push!(problems, "request $k: is_current returned a $(typeof(current)), not a Bool")
    end
    isempty(candidates) || accepted || push!(problems, "validate_proposal accepted none of the proposals, so " *
        "install! and is_current never ran; include a request the family accepts")
    return _check(:family_obligations, since, isempty(problems), join(unique(problems), "\n"))
end
