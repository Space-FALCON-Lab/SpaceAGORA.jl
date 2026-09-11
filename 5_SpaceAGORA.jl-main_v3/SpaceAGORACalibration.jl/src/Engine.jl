module Engine

using Random
using TOML

using ..Spec: CalibrationSpec, primary_manifest_path, validate_spec
using ..ParamSpace: Candidate, candidate_signature
using ..Backend: AbstractBackend, BackendEvaluation, CandidateRuntimePolicy
using ..Backend: apply_candidate_to_manifest!, backend_full_auto_requested, backend_parallel_profile, evaluate_candidate
using ..GlobalSearch: global_target_count, plan_initial_design, propose_bo_candidate
using ..LocalRefine: run_local_refine
using ..Robustness: RobustnessResult, plan_robustness_candidates, summarize_robustness
using ..Objective: sort_indices_by_score
using ..AdaptiveRuntime: RuntimeController, RuntimeDecision
using ..AdaptiveRuntime: init_runtime_controller, choose_runtime_decision!, record_runtime_feedback!, save_runtime_controller!
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

@inline function _progress_enabled()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_CALIBRATION_PROGRESS", "1")))
    return !(raw in ("0", "false", "off", "no"))
end

@inline function _log_progress(msg::AbstractString)::Nothing
    _progress_enabled() || return nothing
    println(msg)
    flush(stdout)
    return nothing
end

@inline function _fmt_eval_log(stage::String, cand::Candidate, ev::BackendEvaluation)::String
    return "[calibration] stage=$(stage) candidate=$(cand.id) success=$(ev.success) score=$(round(ev.score; digits=6)) runtime_s=$(round(ev.runtime_s; digits=2))"
end

@inline function _parallel_evaluations(spec::CalibrationSpec)::Int
    default = max(1, spec.budgets.parallel_evaluations)
    raw = strip(get(ENV, "SPACEAGORA_CALIBRATION_PARALLEL_EVALUATIONS", ""))
    isempty(raw) && return default
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_CALIBRATION_PARALLEL_EVALUATIONS must be a positive integer, got '$raw'"))
    end
    parsed > 0 || throw(ArgumentError(
        "SPACEAGORA_CALIBRATION_PARALLEL_EVALUATIONS must be a positive integer, got $parsed"
    ))
    return parsed
end

function _evaluate_candidates_ordered(
    backend::AbstractBackend,
    spec::CalibrationSpec,
    candidates::Vector{Candidate};
    stage::String,
    run_dir::Union{Nothing, String},
    workers::Union{Nothing, Int}=nothing,
    runtime_policy::Union{Nothing, CandidateRuntimePolicy}=nothing
)::Vector{BackendEvaluation}
    isempty(candidates) && return BackendEvaluation[]
    base_workers = workers === nothing ? _parallel_evaluations(spec) : max(1, workers)
    worker_budget = min(base_workers, length(candidates))
    if worker_budget <= 1
        out = Vector{BackendEvaluation}(undef, length(candidates))
        for i in eachindex(candidates)
            out[i] = evaluate_candidate(
                backend,
                spec,
                candidates[i];
                stage=stage,
                run_dir=run_dir,
                runtime_policy=runtime_policy
            )
        end
        return out
    end

    gate = Base.Semaphore(worker_budget)
    tasks = Vector{Task}(undef, length(candidates))
    for i in eachindex(candidates)
        cand = candidates[i]
        tasks[i] = @async begin
            Base.acquire(gate)
            try
                return evaluate_candidate(
                    backend,
                    spec,
                    cand;
                    stage=stage,
                    run_dir=run_dir,
                    runtime_policy=runtime_policy
                )
            finally
                Base.release(gate)
            end
        end
    end

    out = Vector{BackendEvaluation}(undef, length(candidates))
    for i in eachindex(tasks)
        out[i] = fetch(tasks[i])
    end
    return out
end

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

@inline function _runtime_policy_from_decision(decision::RuntimeDecision)::CandidateRuntimePolicy
    return CandidateRuntimePolicy(
        outer_parallel_active=decision.workers > 1,
        outer_backend=decision.outer_backend,
        inner_thread_budget=max(1, decision.inner_thread_budget)
    )
end

@inline function _fmt_runtime_decision(stage::String, decision::RuntimeDecision)::String
    return "[calibration] stage=$(stage) runtime_policy workers=$(decision.workers) batch_size=$(decision.batch_size) outer_backend=$(decision.outer_backend) inner_thread_budget=$(decision.inner_thread_budget)"
end

function _execute_stage_plan!(
    store::RunStore,
    state::RunState,
    spec::CalibrationSpec,
    backend::AbstractBackend,
    stage::String,
    plan::Vector{Candidate};
    advance_candidate_id::Bool=true,
    controller::Union{Nothing, RuntimeController}=nothing
)
    if _stage_done(state, stage)
        return _entries_to_vectors(load_stage_entries(store, stage))
    end

    state.stage_status[stage] = "in_progress"
    state.current_stage = stage
    save_state!(store, state)
    _log_progress("[calibration] stage=$(stage) status=in_progress")

    existing_entries = load_stage_entries(store, stage)
    skip = length(existing_entries)

    if skip < length(plan)
        cursor = skip + 1
        while cursor <= length(plan)
            remaining = length(plan) - cursor + 1
            decision = controller === nothing ? nothing : choose_runtime_decision!(controller, stage, remaining)
            worker_budget = decision === nothing ? max(1, _parallel_evaluations(spec)) : decision.workers
            runtime_policy = decision === nothing ? nothing : _runtime_policy_from_decision(decision)
            decision === nothing || _log_progress(_fmt_runtime_decision(stage, decision))

            chunk_size = max(1, min(worker_budget, remaining))
            hi = min(length(plan), cursor + chunk_size - 1)
            chunk = plan[cursor:hi]
            ids = join((string(c.id) for c in chunk), ",")
            _log_progress("[calibration] stage=$(stage) evaluating candidates=[$(ids)]")
            chunk_evals = _evaluate_candidates_ordered(
                backend,
                spec,
                chunk;
                stage=stage,
                run_dir=store.run_dir,
                workers=worker_budget,
                runtime_policy=runtime_policy
            )
            decision === nothing || record_runtime_feedback!(controller, stage, decision, chunk_evals)
            for i in eachindex(chunk)
                cand = chunk[i]
                ev = chunk_evals[i]
                append_evaluation!(store, cand, ev)

                if advance_candidate_id
                    state.next_candidate_id = max(state.next_candidate_id, cand.id + 1)
                end
                save_state!(store, state)
                _log_progress(_fmt_eval_log(stage, cand, ev))
            end
            cursor = hi + 1
        end
    end

    if advance_candidate_id && !isempty(plan)
        state.next_candidate_id = max(state.next_candidate_id, maximum(c.id for c in plan) + 1)
    end

    stage_is_completed!(store, state, stage; next_candidate_id=state.next_candidate_id)
    _log_progress("[calibration] stage=$(stage) status=completed")
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
    _log_progress("[calibration] stage=prepare status=in_progress")

    stage_is_completed!(store, state, "prepare"; next_candidate_id=state.next_candidate_id)
    _log_progress("[calibration] stage=prepare status=completed")
    return state
end

function _run_global_stage!(
    store::RunStore,
    state::RunState,
    spec::CalibrationSpec,
    backend::AbstractBackend;
    controller::Union{Nothing, RuntimeController}=nothing
)::Tuple{Vector{Candidate}, Vector{BackendEvaluation}}
    stage = "global_search_quick"
    if _stage_done(state, stage)
        return _entries_to_vectors(load_stage_entries(store, stage))
    end

    state.stage_status[stage] = "in_progress"
    state.current_stage = stage
    save_state!(store, state)
    _log_progress("[calibration] stage=$(stage) status=in_progress")

    entries = load_stage_entries(store, stage)
    candidates, evals = _entries_to_vectors(entries)
    seen = Set{String}(candidate_signature(c) for c in candidates)

    target = global_target_count(spec)
    initial_plan = plan_initial_design(spec; start_id=1, stage=stage)
    next_id = isempty(candidates) ? 1 : maximum(c.id for c in candidates) + 1

    while length(candidates) < target
        remaining = target - length(candidates)
        batch = Candidate[]
        decision = controller === nothing ? nothing : choose_runtime_decision!(controller, stage, remaining)
        max_parallel = decision === nothing ? max(1, _parallel_evaluations(spec)) : decision.workers
        bo_batch_size = decision === nothing ? max(1, spec.budgets.batch_size) : decision.batch_size
        runtime_policy = decision === nothing ? nothing : _runtime_policy_from_decision(decision)
        decision === nothing || _log_progress(_fmt_runtime_decision(stage, decision))

        if length(candidates) < length(initial_plan)
            init_take = min(remaining, length(initial_plan) - length(candidates), max_parallel)
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
            q = min(remaining, bo_batch_size)
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

        ids = join((string(c.id) for c in batch), ",")
        _log_progress("[calibration] stage=$(stage) evaluating candidates=[$(ids)]")
        batch_evals = _evaluate_candidates_ordered(
            backend,
            spec,
            batch;
            stage=stage,
            run_dir=store.run_dir,
            workers=max_parallel,
            runtime_policy=runtime_policy
        )
        decision === nothing || record_runtime_feedback!(controller, stage, decision, batch_evals)

        for i in eachindex(batch)
            cand = batch[i]
            ev = batch_evals[i]
            append_evaluation!(store, cand, ev)

            push!(candidates, cand)
            push!(evals, ev)
            push!(seen, candidate_signature(cand))

            state.next_candidate_id = max(state.next_candidate_id, cand.id + 1)
            save_state!(store, state)
            _log_progress(_fmt_eval_log(stage, cand, ev))
        end
    end

    stage_is_completed!(store, state, stage; next_candidate_id=next_id)
    _log_progress("[calibration] stage=$(stage) status=completed")
    return candidates, evals
end

function _run_local_stage!(
    store::RunStore,
    state::RunState,
    spec::CalibrationSpec,
    backend::AbstractBackend,
    seed_candidates::Vector{Candidate},
    seed_evals::Vector{BackendEvaluation},
    start_id::Int;
    controller::Union{Nothing, RuntimeController}=nothing
)::Tuple{Vector{Candidate}, Vector{BackendEvaluation}}
    stage = "local_refine_full"
    if _stage_done(state, stage)
        return _entries_to_vectors(load_stage_entries(store, stage))
    end

    state.stage_status[stage] = "in_progress"
    state.current_stage = stage
    save_state!(store, state)
    _log_progress("[calibration] stage=$(stage) status=in_progress")

    existing_entries = load_stage_entries(store, stage)
    prior_candidates, prior_evals = _entries_to_vectors(existing_entries)
    decision = controller === nothing ? nothing : choose_runtime_decision!(controller, stage, max(1, spec.budgets.local_refine_neighbors))
    workers_override = decision === nothing ? nothing : decision.workers
    runtime_policy = decision === nothing ? nothing : _runtime_policy_from_decision(decision)
    decision === nothing || _log_progress(_fmt_runtime_decision(stage, decision))

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
        run_dir=store.run_dir,
        workers_override=workers_override,
        runtime_policy=runtime_policy
    )

    new_evals = BackendEvaluation[]
    for i in (length(prior_candidates) + 1):length(local_result.candidates)
        cand = local_result.candidates[i]
        ev = local_result.evaluations[i]
        push!(new_evals, ev)
        append_evaluation!(store, cand, ev)
        state.next_candidate_id = max(state.next_candidate_id, cand.id + 1)
        save_state!(store, state)
        _log_progress(_fmt_eval_log(stage, cand, ev))
    end
    decision === nothing || record_runtime_feedback!(controller, stage, decision, new_evals)

    stage_is_completed!(store, state, stage; next_candidate_id=state.next_candidate_id)
    _log_progress("[calibration] stage=$(stage) status=completed")
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
    _log_progress("[calibration] run_id=$(store.run_id) run_dir=$(store.run_dir)")
    profile = backend_parallel_profile(backend)
    controller = init_runtime_controller(
        spec,
        profile;
        base_parallel=_parallel_evaluations(spec),
        base_batch_size=spec.budgets.batch_size,
        full_auto_requested=backend_full_auto_requested(backend)
    )
    if controller !== nothing
        _log_progress(
            "[calibration] runtime_policy_controller enabled profile=$(controller.profile_name) machine=$(controller.machine_key) cache=$(controller.cache_path)"
        )
    end

    restored = _restore_completed_result(store, state)
    restored !== nothing && return restored

    _prepare_state!(store, state)

    global_candidates, global_evals = _run_global_stage!(store, state, spec, backend; controller=controller)

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
        local_start_id;
        controller=controller
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
        advance_candidate_id=false,
        controller=controller
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
        controller === nothing || save_runtime_controller!(controller)

        return CalibrationResult(
            run_id=store.run_id,
            run_dir=store.run_dir,
            report_path=report,
            best_candidate=best_record.candidate,
            best_score=best_record.robust_score
        )
    end

    controller === nothing || save_runtime_controller!(controller)
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
