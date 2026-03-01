module LocalRefine

using Random

using ..Spec: CalibrationSpec
using ..ParamSpace: Candidate, candidate_signature, perturb_candidate
using ..Backend: AbstractBackend, BackendEvaluation, evaluate_candidate
using ..Objective: primary_score, sort_indices_by_score

export LocalRefineResult, plan_local_candidates, run_local_refine

Base.@kwdef struct LocalRefineResult
    candidates::Vector{Candidate}
    evaluations::Vector{BackendEvaluation}
    planned_count::Int
end

function plan_local_candidates(
    rng::AbstractRNG,
    spec::CalibrationSpec,
    seed_candidates::Vector{Candidate},
    seed_evals::Vector{BackendEvaluation};
    start_id::Int,
    stage::String="local_refine_full"
)::Vector{Candidate}
    isempty(seed_candidates) && return Candidate[]

    ranked = sort_indices_by_score(seed_evals)
    topk = min(spec.budgets.local_refine_topk, length(ranked))

    out = Candidate[]
    next_id = start_id

    for i in 1:topk
        base_idx = ranked[i]
        base = seed_candidates[base_idx]
        for _ in 1:spec.budgets.local_refine_steps
            for _ in 1:max(1, spec.budgets.local_refine_neighbors)
                cand = perturb_candidate(
                    rng,
                    spec,
                    base,
                    next_id;
                    stage=stage,
                    scale=spec.budgets.local_refine_init_scale,
                    perturb_discrete=false
                )
                push!(out, cand)
                next_id += 1
            end
        end
    end

    return out
end

@inline function _planned_count(spec::CalibrationSpec, n_seed::Int)::Int
    topk = min(spec.budgets.local_refine_topk, n_seed)
    return topk * spec.budgets.local_refine_steps * max(1, spec.budgets.local_refine_neighbors)
end

function run_local_refine(
    rng::AbstractRNG,
    spec::CalibrationSpec,
    backend::AbstractBackend,
    seed_candidates::Vector{Candidate},
    seed_evals::Vector{BackendEvaluation};
    start_id::Int,
    stage::String="local_refine_full",
    prior_candidates::Vector{Candidate}=Candidate[],
    prior_evals::Vector{BackendEvaluation}=BackendEvaluation[],
    run_dir::Union{Nothing, String}=nothing
)::LocalRefineResult
    _ = rng
    isempty(seed_candidates) && return LocalRefineResult(
        candidates=Candidate[],
        evaluations=BackendEvaluation[],
        planned_count=0
    )

    length(prior_candidates) == length(prior_evals) || throw(ArgumentError(
        "prior_candidates and prior_evals must have equal lengths."
    ))

    ranked = sort_indices_by_score(seed_evals)
    topk = min(spec.budgets.local_refine_topk, length(ranked))
    n_neighbors = max(1, spec.budgets.local_refine_neighbors)

    all_candidates = copy(prior_candidates)
    all_evals = copy(prior_evals)
    prior_idx = 1

    next_id = isempty(all_candidates) ? start_id : max(start_id, maximum(c.id for c in all_candidates) + 1)

    for i in 1:topk
        base_idx = ranked[i]
        base = seed_candidates[base_idx]
        center = base
        center_score = primary_score(seed_evals[base_idx])
        scale = spec.budgets.local_refine_init_scale

        for step in 1:spec.budgets.local_refine_steps
            step_candidates = Candidate[]
            step_evals = BackendEvaluation[]

            for nb in 1:n_neighbors
                if prior_idx <= length(prior_candidates)
                    cand = prior_candidates[prior_idx]
                    ev = prior_evals[prior_idx]
                    prior_idx += 1
                else
                    hash_seed = hash((
                        spec.seed,
                        stage,
                        base.id,
                        step,
                        nb,
                        round(scale, digits=8),
                        candidate_signature(center),
                        next_id
                    ))
                    local_rng = MersenneTwister(hash_seed)
                    parent = spec.budgets.local_refine_strategy == "trust_region" ? center : base
                    cand = perturb_candidate(
                        local_rng,
                        spec,
                        parent,
                        next_id;
                        stage=stage,
                        scale=scale,
                        perturb_discrete=false
                    )
                    ev = evaluate_candidate(backend, spec, cand; stage=stage, run_dir=run_dir)
                end

                next_id = max(next_id, cand.id + 1)
                push!(step_candidates, cand)
                push!(step_evals, ev)
                push!(all_candidates, cand)
                push!(all_evals, ev)
            end

            if spec.budgets.local_refine_strategy == "trust_region"
                step_scores = [primary_score(ev) for ev in step_evals]
                best_idx = argmin(step_scores)
                best_score = step_scores[best_idx]
                if isfinite(best_score) && (best_score + spec.budgets.local_refine_min_improvement < center_score)
                    center = step_candidates[best_idx]
                    center_score = best_score
                    scale = min(1.0, scale * spec.budgets.local_refine_expand)
                else
                    scale = max(1.0e-6, scale * spec.budgets.local_refine_shrink)
                end
            end
        end
    end

    return LocalRefineResult(
        candidates=all_candidates,
        evaluations=all_evals,
        planned_count=_planned_count(spec, length(seed_candidates))
    )
end

end # module LocalRefine
