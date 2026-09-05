"""
    record_outer_route_feedback!(state, features; route, successes, failures, kwargs...) -> Nothing

Record observed execution feedback for a selected outer route so future adaptive
route selection can use empirical runtime and success statistics.
"""
function record_outer_route_feedback!(
    state::OuterRouteState,
    f::OuterRouteFeatures;
    route::Symbol,
    successes::Int,
    failures::Int,
    elapsed_success_s::Float64=0.0,
    elapsed_success_sq_sum_s::Float64=NaN,
    tuning::OuterRouteTuning=OuterRouteTuning(),
    signature_prefix::String="",
    discard_cold_observation::Bool=true,
    # How many campaigns this observation counts for. The in-campaign split
    # race measures every width after a warm-up sample, in one campaign, and
    # records each as `adaptive_min_samples` observations so the selector can
    # exploit the winner on the next campaign rather than re-trying each width
    # once more. Weight > 1 skips cold eviction: the reading was already warm.
    weight::Int=1,
)::Nothing
    # Arms are not restricted to the three route symbols any more. The split
    # selector (`select_outer_split!`) records under `split_<route>_w<N>` arms in
    # the same per-signature bucket, so it inherits this function's statistics,
    # its confidence handling and its persistence rather than duplicating them.
    # Route arms and split arms never collide because the split arms are
    # prefixed, and each selector only scores the arms it enumerated.
    success_count = max(0, successes)
    failure_count = max(0, failures)
    samples = success_count + failure_count
    samples <= 0 && return nothing

    success_elapsed_sum_s = max(0.0, elapsed_success_s)
    success_elapsed_sq_sum_s = if success_count > 0
        approx_sq = (success_elapsed_sum_s^2) / max(1, success_count)
        provided_sq = max(0.0, elapsed_success_sq_sum_s)
        if isfinite(provided_sq)
            max(provided_sq, approx_sq)
        else
            approx_sq
        end
    else
        0.0
    end

    failure_penalty_s = max(0.0, tuning.failure_penalty_s)
    elapsed_sum_s = success_elapsed_sum_s + failure_count * failure_penalty_s
    elapsed_sq_sum_s = success_elapsed_sq_sum_s + failure_count * failure_penalty_s^2
    elapsed_sq_sum_s = max(elapsed_sq_sum_s, (elapsed_sum_s^2) / samples)

    signatures = _outer_route_signature_hierarchy(f)
    if !isempty(signature_prefix)
        signatures = String[signature_prefix * sig for sig in signatures]
    end
    lock(state.lock) do
        for signature in signatures
            bucket = get!(state.history, signature) do
                Dict{Symbol, OuterRouteStats}()
            end
            stats = get!(bucket, route) do
                OuterRouteStats()
            end
            # Evict this arm's first, cold timing as soon as a warm one exists.
            # See OuterRouteStats.observations for why it is evicted rather than
            # skipped, and for the measurement that motivates it at all.
            stats.observations += max(1, weight)
            stats.campaigns += max(1, weight)
            if discard_cold_observation && weight == 1 && stats.observations == 2
                stats.samples = samples
                stats.successes = success_count
                stats.failures = failure_count
                stats.elapsed_sum_s = elapsed_sum_s
                stats.elapsed_sq_sum_s = elapsed_sq_sum_s
                continue
            end
            stats.samples += samples
            stats.successes += success_count
            stats.failures += failure_count
            stats.elapsed_sum_s += elapsed_sum_s
            stats.elapsed_sq_sum_s += elapsed_sq_sum_s
        end
    end
    return nothing
end
