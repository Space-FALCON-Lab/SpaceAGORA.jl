# General GNC algorithm contract and its shared lifecycle machinery. This module
# loads with Julia standard libraries alone. The RPO planner contract
# (rpo_planner.jl) is its first family.
module GNCAlgorithms

using Random: MersenneTwister
using SHA: sha256

"""
    AbstractGNCAlgorithm

An algorithm whose output is consumed later or can be rejected: a planner, a
scheduler or a decision service. It runs in the planning profile of the GNC
contract: its family builds an owned request, the algorithm proposes a product,
trusted code validates the proposal, and only a validated copy is installed.
Feedback control laws that act at once use the feedback profile instead (a tick
that computes and holds a command).

An algorithm declares `algorithm_family` and implements
`propose!`; the other algorithm-side functions have defaults.
`AbstractRPOPlanner` is a subtype.
"""
abstract type AbstractGNCAlgorithm end

"""
    AbstractAlgorithmFamily

The trusted definition of one problem: its request and product types, how a
proposal is validated (`validate_proposal`), where an accepted product
is installed (`install!`), when an installed product stops being usable
(`is_current`) and what happens on rejection (`on_rejection`).
Algorithms declare a family; they never validate or install their own products.
`RPOReferenceFamily` is a subtype.
"""
abstract type AbstractAlgorithmFamily end

"""
    AbstractAlgorithmRequest

An owned snapshot handed to an algorithm, in SI units, carrying a positive
`request_id`. The lifecycle keeps its own copy for validation and passes the
algorithm a separate deep copy. `RPOPlanningRequest` is a subtype.
"""
abstract type AbstractAlgorithmRequest end

"""
    AbstractAlgorithmResult

An algorithm's outcome for one request: `request_id`, `status` (`:candidate`,
`:infeasible`, `:unsupported` or `:failed`), `termination`, the proposed product
(present exactly for `:candidate`) and diagnostics. A candidate still needs
trusted validation. `AlgorithmResult` and `RPOPlanningResult` are subtypes.
"""
abstract type AbstractAlgorithmResult end

"""
    AbstractValidationResult

The outcome of trusted validation, with fields `accepted::Bool`, `reason::Symbol`
and `metrics`. `ValidationResult` and `RPOValidationResult` are subtypes.
"""
abstract type AbstractValidationResult end

"""
    GNCAlgorithmError <: Exception

Supertype of the errors an algorithm lifecycle raises when a request fails under
the `:fail` policy. Subtypes carry a copied run report, the reason, and the
original cause and backtrace. `RPOPlanningError` is a subtype.
"""
abstract type GNCAlgorithmError <: Exception end

const ALGORITHM_STATUSES = (:candidate, :infeasible, :unsupported, :failed)

"""
    AlgorithmResult(; request_id, status, termination, product=nothing, diagnostics=NamedTuple())

A generic algorithm outcome for families that need no result type of their own.
`status` is one of `:candidate`, `:infeasible`, `:unsupported` or `:failed`; only a
`:candidate` carries a `product`, and it must carry one. Construction deep-copies
the product and the diagnostics, so later edits through an algorithm-held handle
cannot change them. Diagnostics are never read by validation or installation.
"""
struct AlgorithmResult{P, D} <: AbstractAlgorithmResult
    request_id::Int
    status::Symbol
    termination::Symbol
    product::P
    diagnostics::D
    function AlgorithmResult(; request_id::Integer, status::Symbol, termination::Symbol,
            product=nothing, diagnostics=NamedTuple())
        status in ALGORITHM_STATUSES || throw(ArgumentError(
            "Unknown algorithm result status $(repr(status)); use one of $(ALGORITHM_STATUSES)."))
        (status === :candidate) == (product !== nothing) || throw(ArgumentError(
            "A :candidate result carries a product and no other status does."))
        p, d = deepcopy(product), deepcopy(diagnostics)
        new{typeof(p), typeof(d)}(Int(request_id), status, termination, p, d)
    end
end

"""
    ValidationResult(accepted, reason, metrics=NamedTuple())

A generic validation outcome: whether the proposal is accepted, a machine-readable
reason and the measured values behind it.
"""
struct ValidationResult{M} <: AbstractValidationResult
    accepted::Bool
    reason::Symbol
    metrics::M
end
ValidationResult(accepted::Bool, reason::Symbol) = ValidationResult(accepted, reason, NamedTuple())

_request_id(request) = hasproperty(request, :request_id) ? Int(getproperty(request, :request_id)) : 0

# ---- Algorithm side: what an algorithm author implements ----------------------

"""
    algorithm_family(algorithm) -> AbstractAlgorithmFamily

Required. The family whose request, validation, installation and failure rules
govern `algorithm`. The fallback throws an `ArgumentError` that names the type.
"""
function algorithm_family(algorithm)
    throw(ArgumentError("$(nameof(typeof(algorithm))) declares no algorithm family: implement SpaceAGORA.algorithm_family."))
end

"""
    default_capabilities(family)

The capabilities an algorithm of `family` has when it declares none: support for
nothing. The fallback is an empty `NamedTuple`; a family with a capability type
returns its empty instance.
"""
default_capabilities(family) = NamedTuple()

"""
    algorithm_capabilities(algorithm)

What `algorithm` supports, in its family's capability type. Defaults to
`default_capabilities(algorithm_family(algorithm))`, which supports nothing.
Declaring capabilities never loads algorithm configuration.
"""
algorithm_capabilities(algorithm) = default_capabilities(algorithm_family(algorithm))

"""
    require_algorithm_support(algorithm) -> nothing

Throw before a run starts if something `algorithm` needs is missing, such as an
optional package. The default checks nothing.
"""
require_algorithm_support(algorithm) = nothing

"""
    initialize_algorithm(algorithm, run_context) -> state

Create fresh run-owned state for `algorithm`; stateless algorithms return
`nothing` (the default). `run_context` is a copied `NamedTuple` (for example
`seed`, `stream_words`, `epoch` and spacecraft IDs).
"""
initialize_algorithm(algorithm, run_context) = nothing

"""
    propose!(state, algorithm, request, rng) -> AbstractAlgorithmResult

Required. Propose a product for `request` using run-owned `state`. `rng` is the
only random stream the algorithm may draw from; there is no default global stream.
The algorithm receives its own copy of `request` and must not install anything.
The fallback returns an `:unsupported` result and draws nothing.
"""
function propose!(state, algorithm, request, rng)
    return AlgorithmResult(request_id=_request_id(request), status=:unsupported,
        termination=:proposal_not_implemented)
end

"""
    propose_revision!(state, algorithm, request, active_product, rng) -> AbstractAlgorithmResult

Optional. Revise the installed `active_product` instead of proposing from scratch
(RPO calls this retiming). The fallback returns an `:unsupported` result and
draws nothing; an unsupported revision never falls back to a full proposal.
"""
function propose_revision!(state, algorithm, request, active_product, rng)
    return AlgorithmResult(request_id=_request_id(request), status=:unsupported,
        termination=:revision_not_implemented)
end

"""
    consumes_rng(algorithm) -> Bool

Whether `algorithm` draws from the stream passed to `propose!`. The default
is `true`. A deterministic algorithm declares `false`; its contract tests then
check that a proposal leaves the stream unchanged.
"""
consumes_rng(algorithm) = true

# ---- Family side: trusted, core-owned -----------------------------------------

"""
    validate_capabilities(family, algorithm, request; checkpoint=false) -> AbstractValidationResult

Check `algorithm`'s declared capabilities against `request` before asking for a
proposal; `checkpoint=true` also requires restart support. The default accepts,
for families that have no capability rules.
"""
validate_capabilities(family, algorithm, request; checkpoint::Bool=false) =
    ValidationResult(true, :supported)

"""
    validate_proposal(family, request, result, trusted) -> AbstractValidationResult

Required for each family. Validate `result` against the lifecycle's own copy of
`request`. `trusted` holds queries the lifecycle supplies independently of the
algorithm (for RPO, `clearance_at`). Validation never installs a product, never
advances a random stream and never reads diagnostics as geometry.
"""
function validate_proposal end

"""
    install!(link, family, product, request, t; validation, diagnostics) -> nothing

Required for each family. Install a validated copy of `product` into `link`, the
object the family shares with the product's consumer, at time `t`. Installation
is atomic: no yield or user callback runs between validation and installation.
The lifecycle passes the accepted `validation` result and the algorithm's copied
`diagnostics` as keywords, so a family may keep them with the product.
"""
function install! end

"""
    is_current(family, request, product, t) -> Bool

Required for each family. Whether an installed `product` (accepted for `request`)
may still be used at time `t`, including any preview the consumer needs.
"""
function is_current end

"""
    on_rejection(family) -> Symbol

What the lifecycle does when capability checks or validation reject a proposal:
`:fail` raises a `GNCAlgorithmError` (the default), `:keep_previous`
keeps the installed product, and `:install_fallback` installs
`fallback_product`.
"""
on_rejection(family) = :fail

"""
    fallback_product(family, request)

The product installed on rejection by a family whose `on_rejection` is
`:install_fallback`, such as a hold. Families with another policy need no method.
"""
function fallback_product end

"""
    family_contract_version(family) -> Union{Nothing, VersionNumber}

The version of `family`'s request, product and validation contract, so that a
companion package can require the one it was built for. The default `nothing`
means unversioned.
"""
family_contract_version(family) = nothing

# ---- Shared lifecycle machinery ----------------------------------------------
#
# A binding is the guidance model that runs one algorithm in a simulation: the
# generic `AlgorithmGuidance` for any family, or for RPO `PlannerGuidance`. It
# implements the binding functions below; the request driver then runs the
# contract's steps in a fixed order for every family.

"""
    gnc_stream_words(namespace, seed, ids...) -> Vector{UInt32}

The seed words of a binding's random stream, for `MersenneTwister(words)`: SHA-256
over the ASCII text `"SpaceAGORA/<namespace>/<seed>/<id>/..."`, decoded
little-endian into eight 32-bit words (held as unsigned integers; the RPO mapping
has always produced them this way). The namespace carries a version, so a
changed mapping is a new namespace. RPO uses `"RPO/v1"` with the chaser and target
IDs. A stream depends only on these values, never on run order or threads.
"""
function gnc_stream_words(namespace::AbstractString, seed, ids...)
    digest = sha256(join(("SpaceAGORA", namespace, string(seed), map(string, ids)...), "/"))
    return [sum(UInt32(digest[i+j]) << (8j) for j in 0:3) for i in 1:4:32]
end

"""
    AlgorithmRuntime(stream_words; state=nothing)

Run-owned record of one binding, created fresh when the run starts and never
shared between bindings:
- the algorithm's `state` and its random stream `rng` (seeded from `stream_words`);
- the request counter `request_id`;
- the `installed_request` and `product`, which stay `nothing` until a product is
  installed;
- the delivered event tokens (`processed_tokens`);
- the last guidance and request times;
- the run `records`.
"""
mutable struct AlgorithmRuntime
    state::Any
    rng::MersenneTwister
    stream_words::Vector{UInt32}
    request_id::Int
    installed_request::Any
    product::Any
    processed_tokens::Set{Int}
    last_guidance_time::Float64
    last_request_time::Float64
    records::Vector{NamedTuple}
end
function AlgorithmRuntime(stream_words::AbstractVector{<:Unsigned}; state=nothing)
    # Own one copy, so a caller that later reuses its vector cannot change what the
    # runtime reports as the stream's seed.
    words = collect(UInt32, stream_words)
    return AlgorithmRuntime(state, MersenneTwister(words), words, 0, nothing, nothing,
        Set{Int}(), -Inf, -Inf, NamedTuple[])
end

"""
    algorithm_runtime(binding) -> AlgorithmRuntime

Required for each binding: its run-owned `AlgorithmRuntime`.
"""
function algorithm_runtime end

"""
    bound_algorithm(binding)

Required for each binding: the algorithm it runs.
"""
function bound_algorithm end

"""
    bound_family(binding) -> AbstractAlgorithmFamily

The family instance a binding runs its algorithm under. Defaults to
`algorithm_family(bound_algorithm(binding))`. A binding returns its own instance to
configure the family's policy, such as the RPO rejection policy. The instance must
have the same type as the family the algorithm declares; otherwise the request
fails with reason `:family_mismatch`.
"""
bound_family(binding) = algorithm_family(bound_algorithm(binding))

"""
    bound_link(binding)

Required for each binding: the object accepted products are installed into, shared
with their consumer (for RPO, the plan buffer).
"""
function bound_link end

"""
    build_algorithm_request(binding, u, p, t, reason, request_id)

Required for each binding: the owned request snapshot for this tick, built by
trusted code from the simulation state `u` and parameters `p`.
"""
function build_algorithm_request end

"""
    trusted_queries(binding) -> NamedTuple

The queries the lifecycle hands to `validate_proposal`, supplied independently of
the algorithm (for RPO, `clearance_at`). Defaults to none.
"""
trusted_queries(binding) = NamedTuple()

"""
    fail_algorithm_request!(binding, reason; cause=nothing, backtrace=nothing, details=NamedTuple())

Required for each binding: record the failure and throw a `GNCAlgorithmError`
carrying the binding's copied run report. The request driver calls it when a
request fails under the `:fail` policy, and with the driver's `exception_reason`
when algorithm code throws.
"""
function fail_algorithm_request! end

"""
    algorithm_result_type(family) -> Type

The type every proposal of `family` must have. A proposal of any other type fails
the request with reason `:invalid_result_type`. Defaults to `AbstractAlgorithmResult`.
"""
algorithm_result_type(family) = AbstractAlgorithmResult

"""
    is_revision_request(family, request) -> Bool

Whether `request` asks for a revision of the installed product, which the driver
sends to `propose_revision!` instead of `propose!` (RPO: `reason === :retime`).
Defaults to `false`.
"""
is_revision_request(family, request) = false

"""
    proposed_product(result)

The product a candidate `result` proposes: `result.product` for an
`AlgorithmResult` (RPO: the reference). The request driver calls it on every
proposal of the family's result type before validation, so it must also return,
not throw, for a result that validation will refuse.
"""
proposed_product(result::AlgorithmResult) = result.product

"""
    installed_record(family, request, result, validation, product, reason, t) -> NamedTuple

The run record kept when `product` is installed. The driver calls it after
installation, with copies of its arguments, and stores a deep copy of the result.
The default names the product `product`.
"""
installed_record(family, request, result, validation, product, reason, t) =
    (event=:installed, time_s=t, request_id=_request_id(request), reason=reason,
     termination=result.termination, validation=validation, product=product,
     diagnostics=result.diagnostics)

"""
    rejected_record(family, request, result, validation, reason, t) -> NamedTuple

The run record kept when a request is rejected under the `:keep_previous` or
`:install_fallback` policy. `result` is `nothing` when the capability check refused
the request before any proposal.
"""
rejected_record(family, request, result, validation, reason, t) =
    (event=:rejected, time_s=t, request_id=_request_id(request), reason=reason,
     termination=result === nothing ? :not_proposed : result.termination,
     validation=validation)

# Marks the failure this request's binding raised, so the driver passes it on
# unwrapped; any other error, including a GNCAlgorithmError thrown by algorithm
# code, is wrapped with the driver's `exception_reason`. Never escapes the driver.
struct _BindingFailure <: Exception
    error::GNCAlgorithmError
end

function _fail_request!(binding, reason; kwargs...)
    try
        fail_algorithm_request!(binding, reason; kwargs...)
    catch e
        e isa GNCAlgorithmError && throw(_BindingFailure(e))
        rethrow()
    end
    throw(ArgumentError("fail_algorithm_request! for $(nameof(typeof(binding))) must throw a GNCAlgorithmError."))
end

# Record hooks receive isolated copies, so a family's record code cannot change the
# request, result, validation or product the lifecycle keeps or installs.
_rejection(family, request, result, validation, reason, t) =
    deepcopy(rejected_record(family, deepcopy(request), deepcopy(result), deepcopy(validation), reason, t))

function _reject_request!(binding, runtime, family, request, result, validation, reason, t)
    policy = on_rejection(family)
    if policy === :fail
        result === nothing && _fail_request!(binding, validation.reason)
        _fail_request!(binding, validation.reason; details=(termination=result.termination,
            diagnostics=result.diagnostics, validation=validation))
    elseif policy === :keep_previous
        push!(runtime.records, _rejection(family, request, result, validation, reason, t))
        return false
    elseif policy === :install_fallback
        push!(runtime.records, _rejection(family, request, result, validation, reason, t))
        link = bound_link(binding)
        product = deepcopy(fallback_product(family, deepcopy(request)))
        install!(link, family, product, request, t; validation=validation, diagnostics=NamedTuple())
        runtime.product, runtime.installed_request, runtime.last_request_time = product, request, t
        push!(runtime.records, deepcopy((event=:fallback_installed, time_s=t,
            request_id=_request_id(request), reason=reason, product=product)))
        return true
    end
    throw(ArgumentError("Unknown rejection policy $(repr(policy)) for $(nameof(typeof(family)))."))
end

"""
    run_algorithm_request!(binding, u, p, t, reason; exception_reason=:algorithm_exception) -> Bool

Run one request of `binding` at time `t` through the contract's steps:
1. advance the request counter, resolve the algorithm and its family
   (`bound_algorithm`, `bound_family`, checked against `algorithm_family`) and
   build the request (`build_algorithm_request`);
2. record a copy of it;
3. check capabilities (`validate_capabilities`);
4. ask for a proposal, passing a separate copy of the request (`propose!`, or
   `propose_revision!` with a copy of the installed product when
   `is_revision_request`);
5. check the proposal's type and copy it, then resolve the link (`bound_link`) and
   copy the product (`proposed_product`);
6. validate the copy against the retained request (`validate_proposal` with the
   binding's `trusted_queries`);
7. install the product (`install!`) and update the runtime. Nothing else runs
   between validation and installation: no yield, no hook and no user callback;
8. record the install (`installed_record`, called after installation on copies).

A refusal follows the family's `on_rejection` policy. Any error thrown once the
counter has advanced, by algorithm or family code (resolving the family included)
and including a `GNCAlgorithmError`, becomes a failure of this binding with reason
`exception_reason`; only the binding's own failure passes through unwrapped. Every
failure counts as a request: the counter stays advanced. Returns whether a product
was installed.
"""
function run_algorithm_request!(binding, u, p, t, reason; exception_reason::Symbol=:algorithm_exception)
    runtime = algorithm_runtime(binding)::AlgorithmRuntime
    runtime.request_id += 1
    try
        algorithm = bound_algorithm(binding)
        family = bound_family(binding)
        # A binding configures the declared family; it cannot substitute another one.
        Base.typename(typeof(family)) === Base.typename(typeof(algorithm_family(algorithm))) ||
            _fail_request!(binding, :family_mismatch)
        request = build_algorithm_request(binding, u, p, t, reason, runtime.request_id)
        # Retain a trusted snapshot. Algorithm code only receives another copy.
        push!(runtime.records, (event=:request, time_s=t, request=deepcopy(request)))
        capabilities = validate_capabilities(family, algorithm, request)
        capabilities.accepted ||
            return _reject_request!(binding, runtime, family, request, nothing, capabilities, reason, t)
        result = is_revision_request(family, request) ?
            propose_revision!(runtime.state, algorithm, deepcopy(request), deepcopy(runtime.product), runtime.rng) :
            propose!(runtime.state, algorithm, deepcopy(request), runtime.rng)
        result isa algorithm_result_type(family) || _fail_request!(binding, :invalid_result_type)
        result = deepcopy(result)
        # Resolve the link and the product before validating, so that only install!
        # and the runtime update run between validation and installation.
        link = bound_link(binding)
        product = deepcopy(proposed_product(result))
        checked = validate_proposal(family, request, result, trusted_queries(binding))
        checked.accepted ||
            return _reject_request!(binding, runtime, family, request, result, checked, reason, t)
        # No yield or user callback between validation and installation.
        install!(link, family, product, request, t; validation=checked, diagnostics=result.diagnostics)
        runtime.product, runtime.installed_request, runtime.last_request_time = product, request, t
        # The record hook runs after installation, on isolated copies.
        push!(runtime.records, deepcopy(installed_record(family, deepcopy(request), deepcopy(result),
            deepcopy(checked), deepcopy(product), reason, t)))
        return true
    catch e
        e isa _BindingFailure && throw(e.error)
        backtrace = catch_backtrace()
        try
            _fail_request!(binding, exception_reason; cause=e, backtrace=backtrace)
        catch failure
            failure isa _BindingFailure && throw(failure.error)
            rethrow()
        end
    end
end

# ---- Family side for the generic binding -------------------------------------
#
# `AlgorithmGuidance` runs an algorithm of any family through the request driver,
# with no family-specific orchestration. It asks the family for the three pieces
# below; the burn-plan family (`BurnPlanFamily`) implements them.

"""
    family_request(family, settings, u, p, t, spacecraft_index, reason, request_id) -> AbstractAlgorithmRequest

Required for a family that `AlgorithmGuidance` runs: the owned request snapshot
for the spacecraft at `spacecraft_index` (its position in the run, from 1) at
elapsed time `t` in seconds. Trusted family code builds it from the simulation
state `u`, the run's parameters `p` and the binding's `settings`, which the
family defines (`nothing` when it needs none). `reason` is the request's reason
and `request_id` its positive number. The fallback throws an `ArgumentError`
that names the family.
"""
function family_request(family, settings, u, p, t, spacecraft_index, reason, request_id)
    throw(ArgumentError("$(nameof(typeof(family))) builds no requests for AlgorithmGuidance: implement SpaceAGORA.family_request."))
end

"""
    family_link(family, link, settings, p, spacecraft_index)

The object `AlgorithmGuidance` installs accepted products into for the
spacecraft at `spacecraft_index`, resolved once per run before the first
request. `link` is the binding's configured link and `p` the run's parameters,
so a family whose install target belongs to the run (the burn-plan family: the
spacecraft's slot in the burn scheduler's plan buffer) finds it there. The
fallback returns `link`, and throws an `ArgumentError` when it is `nothing`.
"""
function family_link(family, link, settings, p, spacecraft_index)
    link === nothing && throw(ArgumentError(
        "$(nameof(typeof(family))) has no install target of its own: configure the binding's `link`."))
    return link
end

"""
    preflight_family(family, algorithm, settings, args, spacecraft_index; isolate_state) -> nothing

Check, before a run starts, that the configuration `args` can run `algorithm`
of `family` for the spacecraft at `spacecraft_index` with the binding's
`settings`, and throw an `ArgumentError` naming the problem otherwise (for
example a missing consumer of the family's products). `isolate_state` is the
run's `isolate_state` setting. The default checks nothing.
"""
preflight_family(family, algorithm, settings, args, spacecraft_index; isolate_state::Bool=true) = nothing

end # module GNCAlgorithms
