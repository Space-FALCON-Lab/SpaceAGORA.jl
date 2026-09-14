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
    tuning::OuterRouteTuning=OuterRouteTuning()
)::Nothing
    route in (:none, :threads, :process) || return nothing
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
    lock(state.lock) do
        for signature in signatures
            bucket = get!(state.history, signature) do
                Dict{Symbol, OuterRouteStats}()
            end
            stats = get!(bucket, route) do
                OuterRouteStats()
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
