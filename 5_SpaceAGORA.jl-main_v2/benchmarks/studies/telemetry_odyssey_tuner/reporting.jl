function _evaluate_candidates_batch(
    cfg::TunerConfig,
    stage::Symbol,
    candidates::Vector{TuneCandidate},
    base_scenario::Dict{String, Any}
)::DataFrame
    n = length(candidates)
    n == 0 && return DataFrame()
    stage_name = String(stage)
    rows = Vector{NamedTuple}(undef, n)
    progress_lock = ReentrantLock()

    worker_count = min(cfg.parallel_candidates, n)
    if worker_count <= 1
        for idx in 1:n
            cand = candidates[idx]
            _progress_print(progress_lock, @sprintf("[%s %d/%d] candidate=%d start", stage_name, idx, n, cand.id))
            tuned = _apply_candidate(base_scenario, cand; calibration_enabled=cfg.calibration_enabled)
            row = evaluate_candidate(cfg, stage, cand, tuned; progress_lock=progress_lock)
            rows[idx] = row
            _progress_print(
                progress_lock,
                @sprintf(
                    "[%s %d/%d] candidate=%d done success=%s pass=%s score=%s wall_s=%.1f",
                    stage_name,
                    idx,
                    n,
                    cand.id,
                    string(row.success),
                    string(row.pass),
                    string(row.score),
                    Float64(row.wall_runtime_s)
                )
            )
        end
        return DataFrame(rows)
    end

    jobs = Channel{Int}(n)
    for idx in 1:n
        put!(jobs, idx)
    end
    close(jobs)

    @sync for _ in 1:worker_count
        Threads.@spawn begin
            for idx in jobs
                cand = candidates[idx]
                _progress_print(progress_lock, @sprintf("[%s %d/%d] candidate=%d start", stage_name, idx, n, cand.id))
                tuned = _apply_candidate(base_scenario, cand; calibration_enabled=cfg.calibration_enabled)
                row = evaluate_candidate(cfg, stage, cand, tuned; progress_lock=progress_lock)
                rows[idx] = row
                _progress_print(
                    progress_lock,
                    @sprintf(
                        "[%s %d/%d] candidate=%d done success=%s pass=%s score=%s wall_s=%.1f",
                        stage_name,
                        idx,
                        n,
                        cand.id,
                        string(row.success),
                        string(row.pass),
                        string(row.score),
                        Float64(row.wall_runtime_s)
                    )
                )
            end
        end
    end
    return DataFrame(rows)
end

function _candidate_from_row(row)::TuneCandidate
    return TuneCandidate(
        id=Int(row.candidate_id),
        epoch_shift_s=Int(row.epoch_shift_s),
        ra_scale=Float64(row.ra_scale),
        rp_altitude_offset_m=Float64(row.rp_altitude_offset_m),
        i_offset_deg=Float64(row.i_offset_deg),
        aop_offset_deg=Float64(row.aop_offset_deg),
        raan_offset_deg=Float64(row.raan_offset_deg),
        ta_offset_deg=Float64(row.ta_offset_deg),
        bus_mass_scale=Float64(row.bus_mass_scale),
        prop_mass_scale=Float64(row.prop_mass_scale),
        panel_mass_scale=Float64(row.panel_mass_scale),
        bus_dims_scale=Float64(row.bus_dims_scale),
        panel_dims_scale=Float64(row.panel_dims_scale),
        panel_offset_scale=Float64(row.panel_offset_scale),
        srp_cr_scale=Float64(row.srp_cr_scale),
        srp_area_scale=Float64(row.srp_area_scale),
        gravity_variant=String(row.gravity_variant),
        dv_global_scale=Float64(row.dv_global_scale),
        dv_early_scale=Float64(row.dv_early_scale),
        dv_orbit7_bias_mps=Float64(row.dv_orbit7_bias_mps),
        solver_mode=String(row.solver_mode_requested),
        dt_max_orbit_s=Float64(row.dt_max_orbit_requested_s),
        dt_max_atm_s=Float64(row.dt_max_atm_requested_s)
    )
end

function _write_report(
    path::String,
    cfg::TunerConfig,
    quick_df::DataFrame,
    full_df::DataFrame,
    best_row,
    best_manifest::String;
    plot_paths::Vector{String}=String[]
)
    generated = string(now(UTC))
    open(path, "w") do io
        println(io, "# Odyssey Telemetry Tuning Report")
        println(io)
        println(io, "- Generated (UTC): $generated")
        println(io, "- Base manifest: `$(cfg.base_manifest)`")
        println(io, "- Scenario: `$(cfg.scenario_name)`")
        println(io, "- Quick candidates: `$(cfg.quick_candidates)`")
        println(io, "- Full top-k: `$(cfg.full_topk)`")
        println(io, "- Parallel candidates: `$(cfg.parallel_candidates)`")
        println(io, "- Calibration enabled: `$(cfg.calibration_enabled)`")
        println(io, "- Best manifest: `$(best_manifest)`")
        println(io)

        println(io, "## Best Candidate")
        println(io)
        println(io, "- Stage selected: `$(best_row.stage)`")
        println(io, "- Score: `$(best_row.score)`")
        println(io, "- Passes paper limits: `$(best_row.pass)`")
        println(io, "- Peri RMSE [km]: `$(best_row.peri_rmse_km)`")
        println(io, "- Apo RMSE [km]: `$(best_row.apo_rmse_km)`")
        println(io, "- Peri max abs [km]: `$(best_row.peri_max_abs_km)`")
        println(io, "- Apo max abs [km]: `$(best_row.apo_max_abs_km)`")
        println(io)

        println(io, "## Reproduce Best")
        println(io)
        println(io, "```bash")
        println(io, "JULIA_NUM_THREADS=$(Threads.nthreads()) julia --project=.AGORA --startup-file=no benchmarks/studies/telemetry_orbit_accuracy_study.jl --profile=full --manifest=$(best_manifest) --enforce=$(cfg.enforce ? "1" : "0")")
        println(io, "```")
        println(io)
        println(io, "## Artifacts")
        println(io)
        println(io, "- Quick results: `$(joinpath(cfg.outdir, "odyssey_tuning_quick_results.csv"))`")
        println(io, "- Full results: `$(joinpath(cfg.outdir, "odyssey_tuning_full_results.csv"))`")
        println(io, "- Combined results: `$(joinpath(cfg.outdir, "odyssey_tuning_all_results.csv"))`")
        println(io, "- Best manifest: `$(best_manifest)`")
        if !isempty(plot_paths)
            println(io)
            println(io, "## Plots")
            println(io)
            for pth in plot_paths
                println(io, "- `$(pth)`")
            end
        end
    end
end

function main_tuner()
    cfg = parse_tuner_cli(ARGS)
    mkpath(cfg.outdir)
    base_manifest = TOML.parsefile(cfg.base_manifest)
    base_scenario = _find_base_scenario(base_manifest, cfg.scenario_name)
    candidates = sample_candidates(cfg)

    println("Odyssey tuner")
    println("outdir=$(cfg.outdir)")
    println("scenario=$(cfg.scenario_name)")
    println("quick_candidates=$(length(candidates)) full_topk=$(cfg.full_topk) parallel_candidates=$(cfg.parallel_candidates) quick_only=$(cfg.quick_only)")

    quick_df = _evaluate_candidates_batch(cfg, :quick, candidates, base_scenario)
    quick_csv = joinpath(cfg.outdir, "odyssey_tuning_quick_results.csv")
    CSV.write(quick_csv, quick_df)

    valid_quick = quick_df[quick_df.success .== true, :]
    if nrow(valid_quick) == 0
        error("All quick candidates failed. Check logs under $(joinpath(cfg.outdir, "logs")).")
    end
    sort!(valid_quick, [:pass, :score, :peri_rmse_km, :apo_rmse_km], rev=[true, false, false, false])

    full_candidates = TuneCandidate[]
    if !cfg.quick_only
        k = min(cfg.full_topk, nrow(valid_quick))
        for i in 1:k
            row = valid_quick[i, :]
            push!(full_candidates, _candidate_from_row(row))
        end
    end
    full_df = _evaluate_candidates_batch(cfg, :full, full_candidates, base_scenario)
    full_csv = joinpath(cfg.outdir, "odyssey_tuning_full_results.csv")
    CSV.write(full_csv, full_df)

    combined_df = vcat(quick_df, full_df; cols=:union)
    combined_csv = joinpath(cfg.outdir, "odyssey_tuning_all_results.csv")
    CSV.write(combined_csv, combined_df)

    best_row = nothing
    if nrow(full_df) > 0
        valid_full = full_df[full_df.success .== true, :]
        if nrow(valid_full) > 0
            sort!(valid_full, [:pass, :score, :peri_rmse_km, :apo_rmse_km], rev=[true, false, false, false])
            best_row = valid_full[1, :]
        else
            sort!(full_df, :score)
            best_row = full_df[1, :]
        end
    else
        best_row = valid_quick[1, :]
    end

    best_cand = _candidate_from_row(best_row)
    best_scenario = _apply_candidate(base_scenario, best_cand; calibration_enabled=cfg.calibration_enabled)
    best_manifest_doc = Dict{String, Any}("version" => 1, "scenarios" => Any[best_scenario])
    best_manifest_path = joinpath(cfg.outdir, "odyssey_tuning_best_manifest.toml")
    open(best_manifest_path, "w") do io
        TOML.print(io, best_manifest_doc)
    end

    plot_paths = _generate_tuning_plots(cfg.outdir, quick_df, full_df)
    report_path = joinpath(cfg.outdir, "odyssey_tuning_report.md")
    _write_report(report_path, cfg, quick_df, full_df, best_row, best_manifest_path; plot_paths=plot_paths)

    println()
    println("Tuning complete.")
    println("Quick results: $quick_csv")
    println("Full results : $full_csv")
    println("Combined CSV : $combined_csv")
    println("Best manifest: $best_manifest_path")
    println("Plots        : $(length(plot_paths))")
    println("Report       : $report_path")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_tuner()
end
