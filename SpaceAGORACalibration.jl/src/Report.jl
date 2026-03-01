module Report

using Dates

using ..Spec: CalibrationSpec
using ..ParamSpace: Candidate
using ..Backend: BackendEvaluation
using ..Robustness: RobustnessResult
using ..Store: RunStore, evaluations_path, state_path, best_manifest_path, report_path

export write_report!

@inline function _candidate_values_text(c::Candidate)::String
    keys_sorted = sort(collect(keys(c.values)))
    return join([string(k, "=", c.values[k]) for k in keys_sorted], ", ")
end

function write_report!(
    store::RunStore,
    spec::CalibrationSpec,
    global_evals::Vector{BackendEvaluation},
    local_evals::Vector{BackendEvaluation},
    robust::RobustnessResult,
    best::Candidate,
    best_score::Float64;
    promote_eval::Union{Nothing, BackendEvaluation}=nothing
)::String
    path = report_path(store)
    open(path, "w") do io
        println(io, "# SpaceAGORACalibration Report")
        println(io)
        println(io, "- Generated (UTC): $(now(UTC))")
        println(io, "- Run ID: `$(store.run_id)`")
        println(io, "- Spec ID: `$(spec.id)`")
        println(io, "- Schema Version: `$(spec.schema_version)`")
        println(io, "- Verification Script: `$(spec.verification_script)`")
        println(io, "- Manifest Paths: `$(join(spec.manifest_paths, ", "))`")
        println(io, "- Robust Ranking: `$(spec.budgets.robust_ranking)`")
        println(io, "- Robustness Uncertainty: `$(spec.budgets.robustness_uncertainty)`")
        println(io, "- Objective: `J(θ) = Σ_s w_s Σ_e [Huber(rmse_km/limit_rmse) + 0.5*Huber(max_abs_km/limit_abs)] + λfail*I(run_failed) + λtime*max(0, runtime/t_budget - 1)`")
        println(io, "- Objective Weights: `lambda_fail=$(spec.budgets.objective_lambda_fail), lambda_time=$(spec.budgets.objective_lambda_time), huber_delta=$(spec.budgets.objective_huber_delta)`")
        println(io)

        println(io, "## Lifecycle")
        println(io)
        println(io, "- prepare")
        println(io, "- global_search_quick")
        println(io, "- local_refine_full")
        println(io, "- robustness_validation")
        println(io, "- promote")
        println(io)

        println(io, "## Search Summary")
        println(io)
        println(io, "- Global evaluations: `$(length(global_evals))`")
        println(io, "- Local refinement evaluations: `$(length(local_evals))`")
        println(io, "- Robustness candidates: `$(length(robust.records))`")
        println(io)

        println(io, "## Best Candidate")
        println(io)
        best_record = nothing
        for rec in robust.records
            if rec.candidate.id == best.id
                best_record = rec
                break
            end
        end
        println(io, "- Candidate ID: `$(best.id)`")
        println(io, "- Robust score: `$(best_score)`")
        println(io, "- Values: `$( _candidate_values_text(best) )`")
        if best_record !== nothing
            println(io, "- Mean score: `$(best_record.mean_score)`")
            println(io, "- P95 score: `$(best_record.p95_score)`")
            println(io, "- CVaR score: `$(best_record.cvar_score)`")
            println(io, "- Fail rate: `$(best_record.fail_rate)`")
        end
        if promote_eval !== nothing
            println(io, "- Full confirmation score: `$(promote_eval.score)`")
            println(io, "- Full confirmation success: `$(promote_eval.success)`")
            println(io, "- Full confirmation runtime [s]: `$(promote_eval.runtime_s)`")
        end
        println(io)

        println(io, "## Artifacts")
        println(io)
        println(io, "- Run directory: `$(store.run_dir)`")
        println(io, "- Spec: `$(store.spec_path)`")
        println(io, "- Evaluations: `$(evaluations_path(store))`")
        println(io, "- State: `$(state_path(store))`")
        println(io, "- Best manifest: `$(best_manifest_path(store))`")
    end
    return path
end

end # module Report
