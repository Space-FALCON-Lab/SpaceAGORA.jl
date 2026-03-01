module Engine

using Random
using TOML

using ..Spec: CalibrationSpec, primary_manifest_path, validate_spec
using ..ParamSpace: Candidate, candidate_signature
using ..Backend: AbstractBackend, BackendEvaluation, apply_candidate_to_manifest!, evaluate_candidate
using ..GlobalSearch: global_target_count, plan_initial_design, propose_bo_candidate
using ..LocalRefine: run_local_refine
using ..Robustness: RobustnessResult, plan_robustness_candidates, summarize_robustness
using ..Objective: sort_indices_by_score
using ..Store: RunStore, RunState, STAGE_SEQUENCE
using ..Store: init_store, append_evaluation!, load_ledger_entries, load_stage_entries
using ..Store: load_state, save_state!, stage_is_completed!, report_path, best_manifest_path
using ..Report: write_report!

export CalibrationResult, run_calibration

Base.@kwdef struct CalibrationResult
    run_id::String
    run_dir::String
    report_path::String
    best_candidate::Candidate
    best_score::Float64
end

@inline _stage_done(state::RunState, stage::String) = get(state.stage_status, stage, "pending") == "completed"

@inline function _entries_to_vectors(entries)
    candidates = Candidate[]
    evals = BackendEvaluation[]
    for e in entries
        push!(candidates, e.candidate)
        push!(evals, e.evaluation)
    end
    return candidates, evals
end

@inline function _next_candidate_id(entries)::Int
    next_id = 1
    for e in entries
        next_id = max(next_id, e.candidate.id + 1)
    end
    return next_id
end

@inline function _top_candidates(
    candidates::Vector{Candidate},
    evals::Vector{BackendEvaluation},
    topk::Int
)::Vector{Candidate}
    isempty(candidates) && return Candidate[]
    ranked = sort_indices_by_score(evals)
    k = min(topk, length(ranked))
    return [candidates[ranked[i]] for i in 1:k]
end

function _execute_stage_plan!(
    store::RunStore,
    state::RunState,
    spec::CalibrationSpec,
    backend::AbstractBackend,
    stage::String,
    plan::Vector{Candidate};
    advance_candidate_id::Bool=true
)
    if _stage_done(state, stage)
        return _entries_to_vectors(load_stage_entries(store, stage))
    end

    state.stage_status[stage] = "in_progress"
    state.current_stage = stage
    save_state!(store, state)

    existing_entries = load_stage_entries(store, stage)
    skip = length(existing_entries)

    if skip < length(plan)
        for i in (skip + 1):length(plan)
            cand = plan[i]
            ev = evaluate_candidate(backend, spec, cand; stage=stage, run_dir=store.run_dir)
            append_evaluation!(store, cand, ev)

            if advance_candidate_id
                state.next_candidate_id = max(state.next_candidate_id, cand.id + 1)
            end
            save_state!(store, state)
        end
    end

    if advance_candidate_id && !isempty(plan)
        state.next_candidate_id = max(state.next_candidate_id, maximum(c.id for c in plan) + 1)
    end

    stage_is_completed!(store, state, stage; next_candidate_id=state.next_candidate_id)
    return _entries_to_vectors(load_stage_entries(store, stage))
end

function _group_evals_by_candidate_id(entries)::Dict{Int, Vector{BackendEvaluation}}
    grouped = Dict{Int, Vector{BackendEvaluation}}()
    for e in entries
        bucket = get!(grouped, e.candidate.id, BackendEvaluation[])
        push!(bucket, e.evaluation)
    end
    return grouped
end

function _write_best_manifest!(store::RunStore, spec::CalibrationSpec, best::Candidate)::String
    base_manifest = primary_manifest_path(spec)
    if isfile(base_manifest)
        doc = TOML.parsefile(base_manifest)
        apply_candidate_to_manifest!(doc, spec, best)
        open(best_manifest_path(store), "w") do io
            TOML.print(io, doc)
        end
    else
        payload = Dict(
            "source_manifest_missing" => base_manifest,
            "candidate_id" => best.id,
            "candidate_values" => Dict(best.values)
        )
        open(best_manifest_path(store), "w") do io
            TOML.print(io, payload)
        end
    end
    return best_manifest_path(store)
end

function _prepare_state!(store::RunStore, state::RunState)
    if _stage_done(state, "prepare")
        return state
    end

    state.stage_status["prepare"] = "in_progress"
    state.current_stage = "prepare"
    all_entries = load_ledger_entries(store)
    state.next_candidate_id = _next_candidate_id(all_entries)
    save_state!(store, state)

    stage_is_completed!(store, state, "prepare"; next_candidate_id=state.next_candidate_id)
    return state
end

function _run_global_stage!(
    store::RunStore,
    state::RunState,
    spec::CalibrationSpec,
    backend::AbstractBackend
)::Tuple{Vector{Candidate}, Vector{BackendEvaluation}}
    stage = "global_search_quick"
    if _stage_done(state, stage)
        return _entries_to_vectors(load_stage_entries(store, stage))
    end

    state.stage_status[stage] = "in_progress"
    state.current_stage = stage
    save_state!(store, state)

    entries = load_stage_entries(store, stage)
    candidates, evals = _entries_to_vectors(entries)
    seen = Set{String}(candidate_signature(c) for c in candidates)

    target = global_target_count(spec)
    initial_plan = plan_initial_design(spec; start_id=1, stage=stage)
    next_id = isempty(candidates) ? 1 : maximum(c.id for c in candidates) + 1

    while length(candidates) < target
        remaining = target - length(candidates)
        batch = Candidate[]

        if length(candidates) < length(initial_plan)
            init_take = min(remaining, length(initial_plan) - length(candidates))
            for _ in 1:init_take
                base = initial_plan[length(candidates) + length(batch) + 1]
                sig = candidate_signature(base)
                cand = if sig ∈ seen
                    propose_bo_candidate(
                        spec,
                        candidates,
                        evals,
                        next_id;
                        stage=stage,
                        excluded_signatures=union(seen, Set(candidate_signature(c) for c in batch))
                    )
                else
                    Candidate(id=next_id, values=Dict(base.values), stage=stage)
                end
                push!(batch, cand)
                next_id += 1
            end
        else
            q = min(remaining, max(1, spec.budgets.batch_size))
            batch_seen = Set{String}()
            for _ in 1:q
                cand = propose_bo_candidate(
                    spec,
                    candidates,
                    evals,
                    next_id;
                    stage=stage,
                    excluded_signatures=union(seen, batch_seen)
                )
                push!(batch, cand)
                push!(batch_seen, candidate_signature(cand))
                next_id += 1
            end
        end

        for cand in batch
            ev = evaluate_candidate(backend, spec, cand; stage=stage, run_dir=store.run_dir)
            append_evaluation!(store, cand, ev)

            push!(candidates, cand)
            push!(evals, ev)
            push!(seen, candidate_signature(cand))

            state.next_candidate_id = max(state.next_candidate_id, cand.id + 1)
            save_state!(store, state)
        end
    end

    stage_is_completed!(store, state, stage; next_candidate_id=next_id)
    return candidates, evals
end

function _run_local_stage!(
    store::RunStore,
    state::RunState,
    spec::CalibrationSpec,
    backend::AbstractBackend,
    seed_candidates::Vector{Candidate},
    seed_evals::Vector{BackendEvaluation},
    start_id::Int
)::Tuple{Vector{Candidate}, Vector{BackendEvaluation}}
    stage = "local_refine_full"
    if _stage_done(state, stage)
        return _entries_to_vectors(load_stage_entries(store, stage))
    end

    state.stage_status[stage] = "in_progress"
    state.current_stage = stage
    save_state!(store, state)

    existing_entries = load_stage_entries(store, stage)
    prior_candidates, prior_evals = _entries_to_vectors(existing_entries)

    rng = MersenneTwister(hash((spec.seed, stage)))
    local_result = run_local_refine(
        rng,
        spec,
        backend,
        seed_candidates,
        seed_evals;
        start_id=start_id,
        stage=stage,
        prior_candidates=prior_candidates,
        prior_evals=prior_evals,
        run_dir=store.run_dir
    )

    for i in (length(prior_candidates) + 1):length(local_result.candidates)
        cand = local_result.candidates[i]
        ev = local_result.evaluations[i]
        append_evaluation!(store, cand, ev)
        state.next_candidate_id = max(state.next_candidate_id, cand.id + 1)
        save_state!(store, state)
    end

    stage_is_completed!(store, state, stage; next_candidate_id=state.next_candidate_id)
    return _entries_to_vectors(load_stage_entries(store, stage))
end

function _restore_completed_result(store::RunStore, state::RunState)::Union{Nothing, CalibrationResult}
    if !_stage_done(state, "promote")
        return nothing
    end
    if state.best_candidate_id <= 0 || !isfinite(state.best_score)
        return nothing
    end

    best = Candidate(
        id=state.best_candidate_id,
        values=Dict{String, Any}(state.best_candidate_values),
        stage="promote"
    )

    return CalibrationResult(
        run_id=store.run_id,
        run_dir=store.run_dir,
        report_path=report_path(store),
        best_candidate=best,
        best_score=state.best_score
    )
end

function run_calibration(spec::CalibrationSpec, backend::AbstractBackend)::CalibrationResult
    validate_spec(spec)

    store = init_store(spec)
    state = load_state(store)

    restored = _restore_completed_result(store, state)
    restored !== nothing && return restored

    _prepare_state!(store, state)

    global_candidates, global_evals = _run_global_stage!(store, state, spec, backend)

    isempty(global_candidates) && throw(ArgumentError("No global candidates available for calibration."))

    local_start_id = isempty(global_candidates) ? 1 : maximum(c.id for c in global_candidates) + 1
    local_seed_candidates = _top_candidates(global_candidates, global_evals, spec.budgets.local_refine_topk)
    local_seed_evals = begin
        ranked = sort_indices_by_score(global_evals)
        k = min(spec.budgets.local_refine_topk, length(ranked))
        [global_evals[ranked[i]] for i in 1:k]
    end

    local_candidates, local_evals = _run_local_stage!(
        store,
        state,
        spec,
        backend,
        local_seed_candidates,
        local_seed_evals,
        local_start_id
    )

    finalists = if isempty(local_candidates)
        _top_candidates(global_candidates, global_evals, spec.budgets.local_refine_topk)
    else
        _top_candidates(local_candidates, local_evals, spec.budgets.local_refine_topk)
    end

    robust_stage = "robustness_validation"
    robust_plan = plan_robustness_candidates(spec, finalists; stage=robust_stage)
    _execute_stage_plan!(
        store,
        state,
        spec,
        backend,
        robust_stage,
        robust_plan;
        advance_candidate_id=false
    )

    robust_entries = load_stage_entries(store, robust_stage)
    grouped = _group_evals_by_candidate_id(robust_entries)
    robust_result = summarize_robustness(spec, finalists, grouped)

    isempty(robust_result.records) && throw(ArgumentError("No robustness records were produced."))
    best_idx = argmin([r.robust_score for r in robust_result.records])
    best_record = robust_result.records[best_idx]

    promote_stage = "promote"
    if !_stage_done(state, promote_stage)
        state.stage_status[promote_stage] = "in_progress"
        state.current_stage = promote_stage
        save_state!(store, state)

        promote_candidate = Candidate(
            id=best_record.candidate.id,
            values=Dict{String, Any}(best_record.candidate.values),
            stage=promote_stage
        )
        promote_entries = load_stage_entries(store, promote_stage)
        promote_eval = if isempty(promote_entries)
            ev = evaluate_candidate(backend, spec, promote_candidate; stage=promote_stage, run_dir=store.run_dir)
            append_evaluation!(store, promote_candidate, ev)
            ev
        else
            promote_entries[end].evaluation
        end

        _write_best_manifest!(store, spec, best_record.candidate)

        report = write_report!(
            store,
            spec,
            global_evals,
            local_evals,
            robust_result,
            best_record.candidate,
            best_record.robust_score;
            promote_eval=promote_eval
        )

        state.best_candidate_id = best_record.candidate.id
        state.best_candidate_values = Dict{String, Any}(best_record.candidate.values)
        state.best_score = best_record.robust_score
        save_state!(store, state)

        stage_is_completed!(store, state, promote_stage; next_candidate_id=state.next_candidate_id)

        return CalibrationResult(
            run_id=store.run_id,
            run_dir=store.run_dir,
            report_path=report,
            best_candidate=best_record.candidate,
            best_score=best_record.robust_score
        )
    end

    return CalibrationResult(
        run_id=store.run_id,
        run_dir=store.run_dir,
        report_path=report_path(store),
        best_candidate=Candidate(
            id=state.best_candidate_id,
            values=Dict{String, Any}(state.best_candidate_values),
            stage="promote"
        ),
        best_score=state.best_score
    )
end

end # module Engine
