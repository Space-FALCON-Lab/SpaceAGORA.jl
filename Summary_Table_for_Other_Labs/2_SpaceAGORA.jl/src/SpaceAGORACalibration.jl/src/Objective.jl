module Objective

using Statistics
using ..Backend: BackendEvaluation

export primary_score, robust_statistics, robust_score, sort_indices_by_score

@inline function primary_score(ev::BackendEvaluation)::Float64
    if !ev.success || !isfinite(ev.score)
        return Inf
    end
    return ev.score
end

function robust_score(
    samples::Vector{BackendEvaluation};
    method::String="weighted",
    cvar_alpha::Float64=0.9,
    p95_weight::Float64=0.5,
    fail_weight::Float64=5.0
)::Float64
    stats = robust_statistics(samples; cvar_alpha=cvar_alpha)
    if !isfinite(stats.mean)
        return Inf
    end

    if method == "weighted"
        return stats.mean + p95_weight * stats.p95 + fail_weight * stats.fail_rate
    elseif method == "cvar"
        return stats.cvar + fail_weight * stats.fail_rate
    elseif method == "p95"
        return stats.p95 + fail_weight * stats.fail_rate
    end
    throw(ArgumentError("Unsupported robust scoring method '$method'."))
end

function robust_statistics(
    samples::Vector{BackendEvaluation};
    cvar_alpha::Float64=0.9
)::NamedTuple
    isempty(samples) && return (
        mean=Inf,
        p95=Inf,
        cvar=Inf,
        fail_rate=1.0,
        n_success=0,
        n_total=0
    )
    scores = Float64[]
    failures = 0
    for ev in samples
        if ev.success && isfinite(ev.score)
            push!(scores, ev.score)
        else
            failures += 1
        end
    end
    isempty(scores) && return (
        mean=Inf,
        p95=Inf,
        cvar=Inf,
        fail_rate=1.0,
        n_success=0,
        n_total=length(samples)
    )

    sorted = sort(scores)
    p95_idx = max(1, ceil(Int, 0.95 * length(sorted)))
    p95 = sorted[p95_idx]
    cvar_idx = clamp(ceil(Int, cvar_alpha * length(sorted)), 1, length(sorted))
    cvar = mean(sorted[cvar_idx:end])
    fail_rate = failures / length(samples)
    return (
        mean=mean(scores),
        p95=p95,
        cvar=cvar,
        fail_rate=fail_rate,
        n_success=length(scores),
        n_total=length(samples)
    )
end

function sort_indices_by_score(evals::Vector{BackendEvaluation})::Vector{Int}
    idx = collect(eachindex(evals))
    sort!(idx, by=i -> primary_score(evals[i]))
    return idx
end

end # module Objective
