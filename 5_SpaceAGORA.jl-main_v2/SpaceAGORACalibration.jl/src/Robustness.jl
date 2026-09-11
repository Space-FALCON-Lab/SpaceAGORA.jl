module Robustness

using Random

using ..Spec: CalibrationSpec
using ..ParamSpace: Candidate
using ..Backend: AbstractBackend, BackendEvaluation, evaluate_candidate
using ..Objective: robust_statistics, robust_score

export RobustnessRecord, RobustnessResult
export plan_robustness_candidates, summarize_robustness, run_robustness

Base.@kwdef struct RobustnessRecord
    candidate::Candidate
    samples::Vector{BackendEvaluation}
    robust_score::Float64
    mean_score::Float64
    p95_score::Float64
    cvar_score::Float64
    fail_rate::Float64
end

Base.@kwdef struct RobustnessResult
    records::Vector{RobustnessRecord}
end

@inline function _sample_with_replica(c::Candidate, replica::Int, stage::String)::Candidate
    vals = Dict{String, Any}(c.values)
    vals["robustness_replica"] = replica
    return Candidate(id=c.id, values=vals, stage=stage)
end

function plan_robustness_candidates(
    spec::CalibrationSpec,
    finalists::Vector{Candidate};
    stage::String="robustness_validation"
)::Vector{Candidate}
    n_draw = max(spec.budgets.robustness_samples, 1)
    out = Candidate[]
    for cand in finalists
        for replica in 1:n_draw
            push!(out, _sample_with_replica(cand, replica, stage))
        end
    end
    return out
end

function summarize_robustness(
    spec::CalibrationSpec,
    finalists::Vector{Candidate},
    sample_evals::Dict{Int, Vector{BackendEvaluation}}
)::RobustnessResult
    records = RobustnessRecord[]
    for cand in finalists
        samples = get(sample_evals, cand.id, BackendEvaluation[])
        stats = robust_statistics(samples; cvar_alpha=spec.budgets.robust_cvar_alpha)
        score = robust_score(
            samples;
            method=spec.budgets.robust_ranking,
            cvar_alpha=spec.budgets.robust_cvar_alpha,
            p95_weight=spec.budgets.robust_p95_weight,
            fail_weight=spec.budgets.robust_fail_weight
        )
        push!(records, RobustnessRecord(
            candidate=cand,
            samples=samples,
            robust_score=score,
            mean_score=stats.mean,
            p95_score=stats.p95,
            cvar_score=stats.cvar,
            fail_rate=stats.fail_rate
        ))
    end

    return RobustnessResult(records=records)
end

function run_robustness(
    rng::AbstractRNG,
    spec::CalibrationSpec,
    backend::AbstractBackend,
    finalists::Vector{Candidate};
    stage::String="robustness_validation"
)::RobustnessResult
    _ = rng
    sample_evals = Dict{Int, Vector{BackendEvaluation}}()
    for replica in plan_robustness_candidates(spec, finalists; stage=stage)
        ev = evaluate_candidate(backend, spec, replica; stage=stage)
        bucket = get!(sample_evals, replica.id, BackendEvaluation[])
        push!(bucket, ev)
    end
    return summarize_robustness(spec, finalists, sample_evals)
end

end # module Robustness
