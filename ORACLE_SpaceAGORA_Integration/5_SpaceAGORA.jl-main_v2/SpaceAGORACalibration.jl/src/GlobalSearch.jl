module GlobalSearch

using LinearAlgebra
using Random
using Statistics

using ..Spec: CalibrationSpec, ParameterSpec, continuous, integer, categorical
using ..ParamSpace: Candidate, candidate_signature, sample_initial_population, sample_random_candidate
using ..Backend: AbstractBackend, BackendEvaluation, evaluate_candidate

export GlobalSearchResult
export global_target_count, plan_initial_design, propose_bo_candidate, run_global_search

Base.@kwdef struct GlobalSearchResult
    candidates::Vector{Candidate}
    evaluations::Vector{BackendEvaluation}
    planned_count::Int
end

Base.@kwdef struct _SurrogateModel
    X::Matrix{Float64}
    y::Vector{Float64}
    Ainv::Matrix{Float64}
    alpha::Vector{Float64}
    length_scale::Float64
end

@inline function global_target_count(spec::CalibrationSpec)::Int
    return spec.budgets.initial_samples + spec.budgets.global_iters * spec.budgets.batch_size
end

function plan_initial_design(
    spec::CalibrationSpec;
    start_id::Int=1,
    stage::String="global_search_quick"
)::Vector{Candidate}
    n = spec.budgets.initial_samples
    n > 0 || return Candidate[]
    rng = MersenneTwister(hash((spec.seed, stage, "initial_design")))
    return sample_initial_population(
        rng,
        spec,
        n;
        start_id=start_id,
        stage=stage,
        design=spec.budgets.initial_design
    )
end

@inline function _norm_feature(p::ParameterSpec, raw)::Float64
    if p.kind == continuous || p.kind == integer
        span = p.upper - p.lower
        if span <= 0.0
            return 0.5
        end
        return clamp((Float64(raw) - p.lower) / span, 0.0, 1.0)
    elseif p.kind == categorical
        n = length(p.choices)
        n <= 1 && return 0.0
        val = String(raw)
        idx = findfirst(==(val), p.choices)
        idx === nothing && (idx = 1)
        return (idx - 1) / (n - 1)
    end
    return 0.0
end

function _candidate_vector(spec::CalibrationSpec, cand::Candidate)::Vector{Float64}
    vec = Vector{Float64}(undef, length(spec.parameters))
    for i in eachindex(spec.parameters)
        p = spec.parameters[i]
        haskey(cand.values, p.name) || throw(ArgumentError("Candidate $(cand.id) missing parameter $(p.name)."))
        vec[i] = _norm_feature(p, cand.values[p.name])
    end
    return vec
end

@inline function _rbf_kernel(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, length_scale::Float64)::Float64
    d2 = 0.0
    @inbounds for i in eachindex(x)
        d = x[i] - y[i]
        d2 += d * d
    end
    return exp(-d2 / (2.0 * length_scale * length_scale))
end

@inline function _phi(z::Float64)::Float64
    return exp(-0.5 * z * z) / sqrt(2.0 * π)
end

@inline function _Phi(z::Float64)::Float64
    # Smooth normal-CDF approximation (GELU-style) to avoid extra deps.
    return 0.5 * (1.0 + tanh(sqrt(2.0 / π) * (z + 0.044715 * z^3)))
end

function _failure_penalty(evals::Vector{BackendEvaluation})::Float64
    finite_scores = [ev.score for ev in evals if ev.success && isfinite(ev.score)]
    if isempty(finite_scores)
        return 1_000.0
    end
    worst = maximum(finite_scores)
    return worst + max(1.0, 0.25 * abs(worst))
end

function _fit_surrogate(
    spec::CalibrationSpec,
    candidates::Vector{Candidate},
    evals::Vector{BackendEvaluation}
)::Union{Nothing, _SurrogateModel}
    n = length(candidates)
    n == 0 && return nothing

    d = length(spec.parameters)
    X = Matrix{Float64}(undef, n, d)
    for i in 1:n
        X[i, :] = _candidate_vector(spec, candidates[i])
    end

    penalty = _failure_penalty(evals)
    y = Vector{Float64}(undef, n)
    for i in 1:n
        ev = evals[i]
        y[i] = (ev.success && isfinite(ev.score)) ? ev.score : penalty
    end

    ls = spec.budgets.bo_length_scale
    noise = max(spec.budgets.bo_noise, 1.0e-9)

    K = Matrix{Float64}(undef, n, n)
    for i in 1:n
        xi = @view X[i, :]
        for j in i:n
            xj = @view X[j, :]
            kval = _rbf_kernel(xi, xj, ls)
            K[i, j] = kval
            K[j, i] = kval
        end
    end

    A = K + noise * I
    Ainv = Matrix{Float64}(inv(Matrix(A)))
    alpha = Ainv * y

    return _SurrogateModel(X=X, y=y, Ainv=Ainv, alpha=alpha, length_scale=ls)
end

function _predict(model::_SurrogateModel, x::Vector{Float64})::Tuple{Float64, Float64}
    n = size(model.X, 1)
    k = Vector{Float64}(undef, n)
    for i in 1:n
        xi = @view model.X[i, :]
        k[i] = _rbf_kernel(x, xi, model.length_scale)
    end

    mu = dot(k, model.alpha)
    var = 1.0 - dot(k, model.Ainv * k)
    sigma = sqrt(max(var, 1.0e-12))
    return mu, sigma
end

function _ei(mu::Float64, sigma::Float64, best::Float64, xi::Float64)::Float64
    sigma <= 1.0e-12 && return 0.0
    improve = best - mu - xi
    z = improve / sigma
    return improve * _Phi(z) + sigma * _phi(z)
end

@inline function _lcb(mu::Float64, sigma::Float64, kappa::Float64)::Float64
    return mu - kappa * sigma
end

function _candidate_pool(
    rng::AbstractRNG,
    spec::CalibrationSpec,
    n::Int,
    stage::String,
    excluded::Set{String}
)::Vector{Candidate}
    out = Candidate[]
    seen = Set{String}()
    attempts = 0
    max_attempts = max(64, n * 40)

    while length(out) < n && attempts < max_attempts
        c = sample_random_candidate(rng, spec, 0; stage=stage)
        sig = candidate_signature(c)
        if sig ∉ excluded && sig ∉ seen
            push!(out, c)
            push!(seen, sig)
        end
        attempts += 1
    end

    if isempty(out)
        c = sample_random_candidate(rng, spec, 0; stage=stage)
        push!(out, c)
    end

    return out
end

function propose_bo_candidate(
    spec::CalibrationSpec,
    observed_candidates::Vector{Candidate},
    observed_evals::Vector{BackendEvaluation},
    next_id::Int;
    stage::String="global_search_quick",
    excluded_signatures::Set{String}=Set{String}()
)::Candidate
    n_obs = length(observed_candidates)

    pool_rng = MersenneTwister(hash((spec.seed, stage, "bo_pool", next_id)))
    pool = _candidate_pool(
        pool_rng,
        spec,
        spec.budgets.bo_pool_size,
        stage,
        excluded_signatures
    )

    if n_obs == 0
        best = pool[1]
        return Candidate(id=next_id, values=Dict(best.values), stage=stage)
    end

    model = _fit_surrogate(spec, observed_candidates, observed_evals)
    model === nothing && begin
        best = pool[1]
        return Candidate(id=next_id, values=Dict(best.values), stage=stage)
    end

    acq = spec.budgets.global_acquisition
    best_obs = minimum([ev.score for ev in observed_evals if ev.success && isfinite(ev.score)]; init=Inf)
    isfinite(best_obs) || (best_obs = _failure_penalty(observed_evals))

    best_idx = 1
    best_val = acq == "ei" ? -Inf : Inf
    best_sig = candidate_signature(pool[1])

    for i in eachindex(pool)
        cand = pool[i]
        x = _candidate_vector(spec, cand)
        mu, sigma = _predict(model, x)

        value = if acq == "ei"
            _ei(mu, sigma, best_obs, spec.budgets.bo_xi)
        else
            _lcb(mu, sigma, spec.budgets.bo_kappa)
        end

        sig = candidate_signature(cand)
        if acq == "ei"
            if value > best_val || (value == best_val && sig < best_sig)
                best_idx = i
                best_val = value
                best_sig = sig
            end
        else
            if value < best_val || (value == best_val && sig < best_sig)
                best_idx = i
                best_val = value
                best_sig = sig
            end
        end
    end

    best = pool[best_idx]
    return Candidate(id=next_id, values=Dict(best.values), stage=stage)
end

function run_global_search(
    spec::CalibrationSpec,
    backend::AbstractBackend;
    stage::String="global_search_quick",
    prior_candidates::Vector{Candidate}=Candidate[],
    prior_evaluations::Vector{BackendEvaluation}=BackendEvaluation[]
)::GlobalSearchResult
    length(prior_candidates) == length(prior_evaluations) || throw(ArgumentError(
        "prior_candidates and prior_evaluations must have equal lengths."
    ))

    total = global_target_count(spec)
    candidates = copy(prior_candidates)
    evals = copy(prior_evaluations)

    initial_plan = plan_initial_design(spec; start_id=1, stage=stage)
    seen = Set{String}(candidate_signature(c) for c in candidates)

    next_id = isempty(candidates) ? 1 : maximum(c.id for c in candidates) + 1

    while length(candidates) < total
        idx = length(candidates) + 1

        cand = if idx <= length(initial_plan)
            base = initial_plan[idx]
            sig = candidate_signature(base)
            if sig ∈ seen
                propose_bo_candidate(spec, candidates, evals, next_id; stage=stage, excluded_signatures=seen)
            else
                Candidate(id=next_id, values=Dict(base.values), stage=stage)
            end
        else
            propose_bo_candidate(spec, candidates, evals, next_id; stage=stage, excluded_signatures=seen)
        end

        ev = evaluate_candidate(backend, spec, cand; stage=stage)
        push!(candidates, cand)
        push!(evals, ev)
        push!(seen, candidate_signature(cand))
        next_id += 1
    end

    new_count = total - length(prior_candidates)
    start = length(candidates) - new_count + 1
    out_c = new_count > 0 ? candidates[start:end] : Candidate[]
    out_e = new_count > 0 ? evals[start:end] : BackendEvaluation[]

    return GlobalSearchResult(candidates=out_c, evaluations=out_e, planned_count=total)
end

end # module GlobalSearch
