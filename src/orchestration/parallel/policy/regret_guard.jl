# Regret guard (OFF by default): bound the policy's worst case against the
# static baseline using only data the policy already holds. The baseline cannot
# be run alongside, so the caller supplies `baseline_ns`, a prediction of what
# the baseline would cost in the same unit as the `chosen_ns` samples (per-RHS
# or per-solve; the unit is the caller's, the rule only compares like with like).
#
#   SPACEAGORA_POLICY_REGRET_GUARD       0/1, default 0
#   SPACEAGORA_POLICY_REGRET_MARGIN      relative margin (0.10 = 10% slower than
#                                        baseline). No default: required when
#                                        the guard is on.
#   SPACEAGORA_POLICY_REGRET_MIN_SOLVES  samples required before the guard may
#                                        fire (default 3 = the repo's
#                                        SPACEAGORA_RHS_CALIBRATE_HEURISTIC_VOTES
#                                        default; assumed, not measured)

"""
    regret_guard_decision(baseline_ns, chosen_ns, margin, min_solves) -> (fire::Bool, reason::Symbol)

Pure rule. Fires only when at least `min_solves` finite positive samples exist
and the MINIMUM sample exceeds `baseline_ns * (1 + margin)`. The minimum, not
the mean or latest: one warm sample under the line vetoes the guard, so JIT and
first-repeat slowdowns (a cold first sample is routinely 1.5x the warm ones)
cannot fire it alone. Reasons: `:fire`, `:disabled_input`, `:too_few_solves`,
`:within_margin`.
"""
function regret_guard_decision(baseline_ns::Real, chosen_ns::AbstractVector{<:Real},
                               margin::Real, min_solves::Integer)
    (isfinite(baseline_ns) && baseline_ns > 0 && isfinite(margin) && margin >= 0) ||
        return (false, :disabled_input)
    good = [Float64(x) for x in chosen_ns if isfinite(x) && x > 0]
    length(good) >= max(1, min_solves) || return (false, :too_few_solves)
    return minimum(good) > baseline_ns * (1 + margin) ? (true, :fire) : (false, :within_margin)
end

regret_guard_enabled()::Bool = parse_bool_env("SPACEAGORA_POLICY_REGRET_GUARD", false)

function regret_guard_margin()::Float64
    raw = strip(get(ENV, "SPACEAGORA_POLICY_REGRET_MARGIN", ""))
    isempty(raw) && throw(ArgumentError(
        "SPACEAGORA_POLICY_REGRET_GUARD=1 requires SPACEAGORA_POLICY_REGRET_MARGIN " *
        "(relative margin over the baseline, e.g. 0.10); there is no default."))
    m = tryparse(Float64, raw)
    (m !== nothing && isfinite(m) && m >= 0) || throw(ArgumentError(
        "SPACEAGORA_POLICY_REGRET_MARGIN must be a finite float >= 0, got '$raw'"))
    return m
end

function regret_guard_min_solves()::Int
    raw = strip(get(ENV, "SPACEAGORA_POLICY_REGRET_MIN_SOLVES", "3"))
    n = tryparse(Int, raw)
    (n !== nothing && n >= 1) || throw(ArgumentError(
        "SPACEAGORA_POLICY_REGRET_MIN_SOLVES must be an integer >= 1, got '$raw'"))
    return n
end

# Process-lifetime state: the freeze outlives reset_policy_telemetry!.
const _regret_guard_lock = ReentrantLock()
const _regret_guard_samples = Dict{String, Vector{Float64}}()
const _regret_guard_frozen = Ref(false)

regret_guard_frozen()::Bool = _regret_guard_frozen[]

function reset_regret_guard!()::Nothing
    lock(_regret_guard_lock) do
        empty!(_regret_guard_samples)
        _regret_guard_frozen[] = false
    end
    return nothing
end

"""
    regret_guard_observe!(signature, chosen_ns, baseline_ns) -> Bool

Record one solve of the policy-chosen plan for `signature` and return true once
the guard has frozen the process to the baseline plan. No-op (false) when the
guard is off. Fires at most once; the firing is mirrored into the policy
telemetry (`regret_guard_*` fields).
"""
function regret_guard_observe!(signature::String, chosen_ns::Real, baseline_ns::Real)::Bool
    regret_guard_enabled() || return false
    _regret_guard_frozen[] && return true
    margin = regret_guard_margin()
    min_solves = regret_guard_min_solves()
    fire, reason, n = lock(_regret_guard_lock) do
        s = get!(Vector{Float64}, _regret_guard_samples, signature)
        push!(s, Float64(chosen_ns))
        f, r = regret_guard_decision(baseline_ns, s, margin, min_solves)
        f && (_regret_guard_frozen[] = true)
        (f, r, length(s))
    end
    if fire
        ctx = _active_policy_context()
        lock(ctx.lock) do
            t = ctx.telemetry
            t.regret_guard_fired = true
            t.regret_guard_reason = reason
            t.regret_guard_margin = margin
            t.regret_guard_count = Int64(n)
        end
    end
    return fire
end

function regret_guard_snapshot()
    ctx = _active_policy_context()
    lock(ctx.lock) do
        t = ctx.telemetry
        (fired=t.regret_guard_fired, reason=String(t.regret_guard_reason),
         margin=t.regret_guard_margin, count=t.regret_guard_count)
    end
end
