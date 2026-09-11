function _evaluate_tasks_batch(
    cfg::TunerConfig,
    stage::String,
    profile::Symbol,
    tasks::Vector{EvalTask},
    base_manifest::Dict{String, Any}
)::DataFrame
    n = length(tasks)
    n == 0 && return DataFrame()

    rows = Vector{NamedTuple}(undef, n)
    progress_lock = ReentrantLock()

    worker_count = min(cfg.parallel_candidates, n)
    if worker_count <= 1
        for idx in 1:n
            t = tasks[idx]
            _progress_print(progress_lock, @sprintf("[%s %d/%d] candidate=%d start", stage, idx, n, t.candidate.id))
            row = evaluate_candidate(
                cfg,
                stage,
                profile,
                t.candidate,
                base_manifest;
                uncertainty_draw=t.draw,
                progress_lock=progress_lock
            )
            rows[idx] = row
            _progress_print(progress_lock, @sprintf("[%s %d/%d] candidate=%d done success=%s objective=%.6f", stage, idx, n, t.candidate.id, string(row.success), Float64(row.objective)))
        end
        return DataFrame(rows)
    end

    jobs = Channel{Int}(n)
    for i in 1:n
        put!(jobs, i)
    end
    close(jobs)

    @sync for _ in 1:worker_count
        Threads.@spawn begin
            for idx in jobs
                t = tasks[idx]
                _progress_print(progress_lock, @sprintf("[%s %d/%d] candidate=%d start", stage, idx, n, t.candidate.id))
                row = evaluate_candidate(
                    cfg,
                    stage,
                    profile,
                    t.candidate,
                    base_manifest;
                    uncertainty_draw=t.draw,
                    progress_lock=progress_lock
                )
                rows[idx] = row
                _progress_print(progress_lock, @sprintf("[%s %d/%d] candidate=%d done success=%s objective=%.6f", stage, idx, n, t.candidate.id, string(row.success), Float64(row.objective)))
            end
        end
    end

    return DataFrame(rows)
end

function _candidate_from_row(row)::TuneCandidate
    values = Dict{String, Any}(
        "epoch_shift_s" => Int(row.epoch_shift_s),
        "ra_scale" => Float64(row.ra_scale),
        "rp_altitude_offset_m" => Float64(row.rp_altitude_offset_m),
        "i_offset_deg" => Float64(row.i_offset_deg),
        "aop_offset_deg" => Float64(row.aop_offset_deg),
        "raan_offset_deg" => Float64(row.raan_offset_deg),
        "ta_offset_deg" => Float64(row.ta_offset_deg),
        "bus_mass_scale" => Float64(row.bus_mass_scale),
        "prop_mass_scale" => Float64(row.prop_mass_scale),
        "panel_mass_scale" => Float64(row.panel_mass_scale),
        "bus_dims_scale" => Float64(row.bus_dims_scale),
        "panel_dims_scale" => Float64(row.panel_dims_scale),
        "panel_offset_scale" => Float64(row.panel_offset_scale),
        "srp_cr_scale" => Float64(row.srp_cr_scale),
        "srp_area_scale" => Float64(row.srp_area_scale),
        "gravity_variant" => String(row.gravity_variant),
        "dv_global_scale" => Float64(row.dv_global_scale),
        "dv_early_scale" => Float64(row.dv_early_scale),
        "dv_orbit7_bias_mps" => Float64(row.dv_orbit7_bias_mps),
        "solver_mode" => String(row.solver_mode_requested),
        "dt_max_orbit_s" => Float64(row.dt_max_orbit_requested_s),
        "dt_max_atm_s" => Float64(row.dt_max_atm_requested_s)
    )
    return TuneCandidate(id=Int(row.candidate_id), values=values, stage=String(row.stage))
end

function _sort_by_objective(df::DataFrame)::DataFrame
    good = df[df.success .== true, :]
    nrow(good) == 0 && return DataFrame()
    sort!(good, [:pass, :objective], rev=[true, false])
    return good
end

function _p95(vals::Vector{Float64})::Float64
    isempty(vals) && return Inf
    sorted = sort(vals)
    idx = max(1, ceil(Int, 0.95 * length(sorted)))
    return sorted[idx]
end

function _run_global_bo(
    cfg::TunerConfig,
    params::Vector{ParameterSpec},
    base_manifest::Dict{String, Any}
)
    candidates = TuneCandidate[]
    scores = Float64[]
    all_rows = DataFrame()

    seen = Set{String}()
    next_id = 1

    init_candidates = _initial_design(cfg, params, cfg.n_init; start_id=next_id, stage="quick_global_init")
    next_id += length(init_candidates)
    foreach(c -> push!(seen, _candidate_signature(c)), init_candidates)

    init_tasks = [EvalTask(candidate=c) for c in init_candidates]
    init_df = _evaluate_tasks_batch(cfg, "quick_global_init", :quick, init_tasks, base_manifest)
    all_rows = vcat(all_rows, init_df; cols=:union)

    for row in eachrow(init_df)
        push!(candidates, _candidate_from_row(row))
        push!(scores, Float64(row.objective))
    end

    for k in 1:cfg.n_bo_iters
        model = _fit_surrogate(cfg, params, candidates, scores)
        model === nothing && break

        batch = _propose_batch(
            cfg,
            params,
            model,
            scores,
            seen,
            next_id,
            k,
            cfg.batch_size
        )
        isempty(batch) && break
        next_id += length(batch)

        stage = @sprintf("quick_global_bo_%03d", k)
        for i in eachindex(batch)
            batch[i] = TuneCandidate(id=batch[i].id, values=batch[i].values, stage=stage)
        end
        df = _evaluate_tasks_batch(cfg, stage, :quick, [EvalTask(candidate=c) for c in batch], base_manifest)
        all_rows = vcat(all_rows, df; cols=:union)

        for row in eachrow(df)
            push!(candidates, _candidate_from_row(row))
            push!(scores, Float64(row.objective))
        end
    end

    return all_rows, candidates, scores, next_id
end

function _continuous_params(params::Vector{ParameterSpec})::Vector{ParameterSpec}
    return [p for p in params if p.kind == :continuous]
end

function _run_local_refinement(
    cfg::TunerConfig,
    params::Vector{ParameterSpec},
    base_manifest::Dict{String, Any},
    quick_df::DataFrame,
    next_id::Int
)
    sorted = _sort_by_objective(quick_df)
    nrow(sorted) == 0 && return DataFrame(), next_id

    seed_count = min(cfg.local_seed_topk, nrow(sorted))
    seeds = [_candidate_from_row(sorted[i, :]) for i in 1:seed_count]
    cont = _continuous_params(params)
    isempty(cont) && return DataFrame(), next_id

    local_rows = DataFrame()
    stage_prefix = "quick_local"

    for (seed_idx, seed) in enumerate(seeds)
        center = seed
        center_score = Float64(sorted[seed_idx, :objective])
        step_scale = cfg.local_step_init

        for step in 1:cfg.local_steps
            neighbors = TuneCandidate[]
            for p in cont
                base_val = Float64(center.values[p.name])
                span = p.upper - p.lower
                delta = span * step_scale
                for sgn in (-1.0, 1.0)
                    vals = deepcopy(center.values)
                    vals[p.name] = clamp(base_val + sgn * delta, p.lower, p.upper)
                    cand = TuneCandidate(
                        id=next_id,
                        values=vals,
                        stage=@sprintf("%s_seed%02d_step%02d", stage_prefix, seed_idx, step)
                    )
                    push!(neighbors, cand)
                    next_id += 1
                end
            end

            isempty(neighbors) && break
            df = _evaluate_tasks_batch(cfg, neighbors[1].stage, :quick, [EvalTask(candidate=c) for c in neighbors], base_manifest)
            local_rows = vcat(local_rows, df; cols=:union)

            valid = df[df.success .== true, :]
            if nrow(valid) == 0
                step_scale *= cfg.local_step_shrink
                continue
            end
            sort!(valid, :objective)
            best_row = valid[1, :]
            best_score = Float64(best_row.objective)
            if best_score + cfg.local_min_improve < center_score
                center = _candidate_from_row(best_row)
                center_score = best_score
                step_scale *= cfg.local_step_expand
            else
                step_scale *= cfg.local_step_shrink
            end
        end
    end

    return local_rows, next_id
end

function _run_robustness_validation(
    cfg::TunerConfig,
    base_manifest::Dict{String, Any},
    finalists_df::DataFrame
)
    robust_samples = DataFrame()
    robust_rank_rows = NamedTuple[]

    for i in 1:nrow(finalists_df)
        cand = _candidate_from_row(finalists_df[i, :])
        tasks = [EvalTask(candidate=TuneCandidate(id=cand.id, values=deepcopy(cand.values), stage="robust_mc"), draw=j) for j in 1:cfg.robust_draws]
        sample_df = _evaluate_tasks_batch(cfg, "robust_mc", :quick, tasks, base_manifest)
        robust_samples = vcat(robust_samples, sample_df; cols=:union)

        success_mask = sample_df.success .== true
        sample_scores = Float64[sample_df.objective[j] for j in 1:nrow(sample_df) if success_mask[j]]
        fail_rate = nrow(sample_df) == 0 ? 1.0 : count(!, Bool.(sample_df.success)) / nrow(sample_df)

        if isempty(sample_scores)
            mean_j = Inf
            p95_j = Inf
            robust_j = Inf
        else
            mean_j = mean(sample_scores)
            p95_j = _p95(sample_scores)
            robust_j = mean_j + cfg.robust_alpha * p95_j + cfg.robust_beta * fail_rate
        end

        push!(robust_rank_rows, merge(
            _candidate_row_nt(cand),
            (
                stage="robust_rank",
                draws=cfg.robust_draws,
                mean_j=mean_j,
                p95_j=p95_j,
                fail_rate=fail_rate,
                j_robust=robust_j
            )
        ))
    end

    robust_rank_df = DataFrame(robust_rank_rows)
    if nrow(robust_rank_df) > 0
        sort!(robust_rank_df, :j_robust)
    end
    return robust_samples, robust_rank_df
end

function _write_report(
    path::String,
    cfg::TunerConfig,
    quick_df::DataFrame,
    robust_rank_df::DataFrame,
    full_df::DataFrame,
    best_manifest_path::String
)
    open(path, "w") do io
        println(io, "# Telemetry Hybrid Tuner Report")
        println(io)
        println(io, "- Generated (UTC): $(now(UTC))")
        println(io, "- Base manifest: `$(cfg.base_manifest)`")
        println(io, "- Outdir: `$(cfg.outdir)`")
        println(io, "- Global design: `n_init=$(cfg.n_init), n_bo_iters=$(cfg.n_bo_iters), batch=$(cfg.batch_size), acq=$(String(cfg.acquisition))`")
        println(io, "- Local refine: `topk=$(cfg.local_seed_topk), steps=$(cfg.local_steps)`")
        println(io, "- Robustness: `draws=$(cfg.robust_draws), alpha=$(cfg.robust_alpha), beta=$(cfg.robust_beta)`")
        println(io, "- Objective weights: `lambda_fail=$(cfg.lambda_fail), lambda_time=$(cfg.lambda_time), huber_delta=$(cfg.huber_delta)`")
        println(io)

        println(io, "## Objective")
        println(io)
        println(io, "`J(θ) = Σ_s w_s Σ_e [Huber(rmse_km/limit_rmse) + 0.5*Huber(max_abs_km/limit_abs)] + λfail*I(run_failed) + λtime*max(0, runtime/t_budget - 1)`")
        println(io)

        println(io, "## Stage Counts")
        println(io)
        println(io, "- Quick evaluations: `$(nrow(quick_df))`")
        println(io, "- Robustness samples: `$(nrow(robust_rank_df) == 0 ? 0 : nrow(robust_rank_df) * cfg.robust_draws)`")
        println(io, "- Full confirmations: `$(nrow(full_df))`")
        println(io)

        if nrow(robust_rank_df) > 0
            best = robust_rank_df[1, :]
            println(io, "## Best Robust Candidate")
            println(io)
            println(io, "- Candidate ID: `$(best.candidate_id)`")
            println(io, "- J_robust: `$(best.j_robust)`")
            println(io, "- Mean(J): `$(best.mean_j)`")
            println(io, "- P95(J): `$(best.p95_j)`")
            println(io, "- Fail rate: `$(best.fail_rate)`")
            println(io)
        end

        if nrow(full_df) > 0
            best_full = full_df[1, :]
            println(io, "## Full Profile Confirmation")
            println(io, "- Success: `$(best_full.success)`")
            println(io, "- Objective: `$(best_full.objective)`")
            println(io, "- Pass: `$(best_full.pass)`")
            println(io)
        end

        println(io, "## Artifacts")
        println(io)
        println(io, "- Quick results: `$(joinpath(cfg.outdir, "telemetry_hybrid_quick_results.csv"))`")
        println(io, "- Robust sample results: `$(joinpath(cfg.outdir, "telemetry_hybrid_robust_samples.csv"))`")
        println(io, "- Robust ranking: `$(joinpath(cfg.outdir, "telemetry_hybrid_robust_rank.csv"))`")
        println(io, "- Full results: `$(joinpath(cfg.outdir, "telemetry_hybrid_full_results.csv"))`")
        println(io, "- Combined results: `$(joinpath(cfg.outdir, "telemetry_hybrid_all_results.csv"))`")
        println(io, "- Best manifest: `$(best_manifest_path)`")
        println(io)

        println(io, "## Reproduce Best (Full)")
        println(io)
        println(io, "```bash")
        println(io, "JULIA_NUM_THREADS=$(Threads.nthreads()) julia --project=.AGORA --startup-file=no benchmarks/studies/telemetry_orbit_accuracy_study.jl --profile=full --manifest=$(best_manifest_path) --enforce=$(cfg.enforce ? "1" : "0")")
        println(io, "```")
    end
end

function main_hybrid_tuner()
    cfg = parse_cli(copy(ARGS))
    mkpath(cfg.outdir)

    base_manifest = TOML.parsefile(cfg.base_manifest)
    params = _parameter_space(cfg)

    println("Telemetry hybrid tuner")
    println("outdir=$(cfg.outdir)")
    println("manifest=$(cfg.base_manifest)")
    println("global: n_init=$(cfg.n_init) n_bo_iters=$(cfg.n_bo_iters) batch=$(cfg.batch_size) acq=$(cfg.acquisition)")

    quick_global_df, _, _, next_id = _run_global_bo(cfg, params, base_manifest)

    println("running local refinement...")
    local_df, next_id = _run_local_refinement(cfg, params, base_manifest, quick_global_df, next_id)

    quick_df = vcat(quick_global_df, local_df; cols=:union)
    quick_csv = joinpath(cfg.outdir, "telemetry_hybrid_quick_results.csv")
    CSV.write(quick_csv, quick_df)

    valid_quick = _sort_by_objective(quick_df)
    nrow(valid_quick) > 0 || error("No successful quick evaluations produced by hybrid tuner.")

    finalist_count = min(cfg.finalists, nrow(valid_quick))
    finalists_df = valid_quick[1:finalist_count, :]

    println("running robustness Monte Carlo on $(finalist_count) finalists...")
    robust_samples_df, robust_rank_df = _run_robustness_validation(cfg, base_manifest, finalists_df)
    robust_samples_csv = joinpath(cfg.outdir, "telemetry_hybrid_robust_samples.csv")
    robust_rank_csv = joinpath(cfg.outdir, "telemetry_hybrid_robust_rank.csv")
    CSV.write(robust_samples_csv, robust_samples_df)
    CSV.write(robust_rank_csv, robust_rank_df)

    best_row = if nrow(robust_rank_df) > 0
        robust_rank_df[1, :]
    else
        finalists_df[1, :]
    end
    best_candidate = _candidate_from_row(best_row)

    full_df = DataFrame()
    if cfg.full_confirm
        println("confirming best robust candidate on full profile...")
        full_task = [EvalTask(candidate=TuneCandidate(id=best_candidate.id, values=deepcopy(best_candidate.values), stage="full_confirm"))]
        full_df = _evaluate_tasks_batch(cfg, "full_confirm", :full, full_task, base_manifest)
        sort!(full_df, :objective)
    end
    full_csv = joinpath(cfg.outdir, "telemetry_hybrid_full_results.csv")
    CSV.write(full_csv, full_df)

    best_manifest = _apply_candidate_to_manifest(base_manifest, best_candidate, cfg)
    best_manifest_path = joinpath(cfg.outdir, "telemetry_hybrid_best_manifest.toml")
    open(best_manifest_path, "w") do io
        TOML.print(io, best_manifest)
    end

    all_df = vcat(quick_df, robust_samples_df, full_df; cols=:union)
    all_csv = joinpath(cfg.outdir, "telemetry_hybrid_all_results.csv")
    CSV.write(all_csv, all_df)

    report_path = joinpath(cfg.outdir, "telemetry_hybrid_report.md")
    _write_report(report_path, cfg, quick_df, robust_rank_df, full_df, best_manifest_path)

    println()
    println("Hybrid tuning complete.")
    println("Quick results      : $quick_csv")
    println("Robust samples     : $robust_samples_csv")
    println("Robust rank        : $robust_rank_csv")
    println("Full results       : $full_csv")
    println("All results        : $all_csv")
    println("Best manifest      : $best_manifest_path")
    println("Report             : $report_path")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_hybrid_tuner()
end
