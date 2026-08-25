# Run-start measurement of effectors the model cannot count analytically.
#
# The reason this exists rather than a required trait: a mandatory
# `effector_cost_terms` method would have to do something on the `::Any`
# fallback, and every option is bad. Returning zero makes the model silently
# underestimate a user's expensive effector. Throwing breaks every custom
# effector already in the wild. Guessing is confidently wrong. Measuring is the
# only honest fallback, and it is affordable: one call per effector at setup is
# strictly cheaper than a single point of the timed sweep this replaces, which
# evaluated the whole RHS fifteen times per candidate.
#
# It also collapses the cache-invalidation problem. The rule is to cache only
# what is expensive to obtain and stable, and re-derive per run whatever is
# cheap and volatile. Machine constants are expensive and stable, so they are
# cached and fingerprinted. Per-effector costs are cheap and volatile -- a user
# edits their effector, changes a degree from 20 to 50, swaps a coefficient
# table -- so they are never cached and always re-probed. There is no need to
# hash user source, track method world-age, or reason about whether an
# effector's parameters changed, because the axis that is hard to invalidate
# correctly is the one that is not cached at all.

"""
    EffectorProbe

Per-effector measurement taken at run start.

- `ns_per_sat`      measured nanoseconds for one satellite's evaluation.
- `reference_ratio` `ns_per_sat` divided by the reference kernel's per-lane cost
                    measured in the same window. Dimensionless, so a probe taken
                    on a throttled or contended machine stays comparable to
                    constants calibrated on a quiet one.
- `declared`        whether the effector also declared analytic terms, in which
                    case the probe was used to validate rather than to supply
                    the cost.
"""
struct EffectorProbe
    ns_per_sat::Float64
    reference_ratio::Float64
    declared::Bool
end

"""
    validate_declaration(measured_ns, predicted_ns; factor = 2.0) -> Bool

Whether a declared cost agrees with a probe of the same effector.

Declaring `effector_cost_terms` is optional, so a declaration cannot be trusted
by default -- but requiring it to be *right* would defeat the purpose of making
it optional. The compromise is to validate it once, cheaply, against a single
probe, and warn on gross disagreement. A declaration within a factor of two is
accepted: the model consumes rankings, not absolute times, and a factor-of-two
error in one effector among several rarely flips a candidate ordering. An
order-of-magnitude error usually does.

Returns true when the declaration is usable.
"""
@inline function validate_declaration(
    measured_ns::Float64,
    predicted_ns::Float64;
    factor::Float64 = 2.0,
)::Bool
    (isfinite(measured_ns) && measured_ns > 0.0) || return true
    (isfinite(predicted_ns) && predicted_ns > 0.0) || return false
    ratio = measured_ns / predicted_ns
    return (1.0 / factor) <= ratio <= factor
end

"""
    probe_ns_map(probes) -> Dict{Int, Float64}

Reduce a probe table to the `probe => ns` mapping `constellation_work_counts`
consumes. Declared effectors are omitted: their cost comes from their counts,
and including them would double-count.
"""
function probe_ns_map(probes::Dict{Int, EffectorProbe})::Dict{Int, Float64}
    out = Dict{Int, Float64}()
    for (idx, probe) in probes
        probe.declared && continue
        out[idx] = probe.ns_per_sat
    end
    return out
end
