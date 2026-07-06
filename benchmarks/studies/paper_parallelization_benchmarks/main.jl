function _ppb_active_phases(ppb::PPBConfig)::Vector{PPBPhase}
    ids = isempty(ppb.phases) ? Set(p.id for p in PAPER_BENCHMARK_PHASES) : Set(ppb.phases)
    ppb.preview && setdiff!(ids, PPB_PREVIEW_SKIP_PHASES)
    phases = [p for p in PAPER_BENCHMARK_PHASES if p.id in ids]
    ppb.preview && (phases = _ppb_preview_phase.(phases))
    return phases
end

function _ppb_thread_ladder(phase::PPBPhase, ppb::PPBConfig)::Vector{Int}
    base = isempty(ppb.threads) ? _ppc_full_thread_ladder() : ppb.threads
    phase.thread_mode == :max_only && return [maximum(base)]
    phase.thread_mode == :single   && return [1]
    return base
end

function _ppb_build_ppc_config(
    phase::PPBPhase,
    ppb::PPBConfig,
    outdir::String;
    process_workers::Int = ppb.process_workers,
)::PPCConfig
    effective_workers = ppb.preview ? min(process_workers, PPB_PREVIEW_MAX_WORKERS) : process_workers
    return PPCConfig(
        profile         = "full",
        outdir          = outdir,
        modes           = phase.modes,
        cases           = phase.cases,
        parity_cases    = phase.parity_cases,
        threads         = _ppb_thread_ladder(phase, ppb),
        repeats         = phase.repeats,
        warmup          = phase.warmup,
        seed            = ppb.seed,
        solver_mode     = ppb.solver_mode,
        process_workers = effective_workers,
        mc_samples      = phase.mc_samples,
        parity_samples  = 512,
        cpu_pinning     = ppb.cpu_pinning,
    )
end

function _ppb_fmt_elapsed(s::Float64)::String
    s < 60.0  && return "$(round(Int, s))s"
    m  = div(round(Int, s), 60)
    rs = round(Int, s) - 60 * m
    m < 60    && return "$(m)m$(rs)s"
    h  = div(m, 60)
    rm = m - 60 * h
    return "$(h)h$(rm)m$(rs)s"
end

function _ppb_dry_print(phase::PPBPhase, ppb::PPBConfig, outdir::String; process_workers::Int=ppb.process_workers)
    println("[dry-run]   outdir          = $(outdir)")
    println("[dry-run]   cases           = $(join(phase.cases, ", "))")
    println("[dry-run]   parity_cases    = $(join(phase.parity_cases, ", "))")
    println("[dry-run]   modes           = $(join(phase.modes, ", "))")
    println("[dry-run]   threads         = $(_ppb_thread_ladder(phase, ppb))")
    println("[dry-run]   mc_samples      = $(phase.mc_samples)")
    println("[dry-run]   process_workers = $(process_workers)")
    println("[dry-run]   repeats         = $(phase.repeats), warmup = $(phase.warmup)")
end

function _ppb_run_phase(phase::PPBPhase, ppb::PPBConfig, phase_dir::String)::NamedTuple
    started  = time()
    errors   = String[]

    if isempty(phase.worker_ladder)
        # Single run — no process-worker sweep.
        if ppb.dry_run
            println("[dry-run] phase=$(phase.id) — single run")
            _ppb_dry_print(phase, ppb, phase_dir)
        else
            try
                mkpath(phase_dir)
                ppc_run_controller(_ppb_build_ppc_config(phase, ppb, phase_dir))
            catch err
                push!(errors, sprint(showerror, err))
            end
        end
        runs = 1
    else
        # Process-worker sweep — one controller call per worker count.
        runs = length(phase.worker_ladder)
        for w in phase.worker_ladder
            sub_dir = joinpath(phase_dir, "workers_$(lpad(w, 2, '0'))")
            if ppb.dry_run
                println("[dry-run] phase=$(phase.id) — workers=$(w)")
                _ppb_dry_print(phase, ppb, sub_dir; process_workers=w)
            else
                try
                    mkpath(sub_dir)
                    cfg = _ppb_build_ppc_config(phase, ppb, sub_dir; process_workers=w)
                    ppc_run_controller(cfg)
                catch err
                    push!(errors, "workers=$(w): $(sprint(showerror, err))")
                end
            end
        end
    end

    elapsed = time() - started
    status  = isempty(errors) ? "ok" : "error"
    return (; id=phase.id, label=phase.label, status, runs, elapsed, errors)
end

function _ppb_print_summary(results::Vector{NamedTuple})
    sep = "=" ^ 82
    println()
    println(sep)
    println("  Paper Benchmark Run Summary")
    println(sep)
    println(rpad("Phase", 6), rpad("Status", 9), rpad("Runs", 6), rpad("Elapsed", 12), "Label")
    println("-" ^ 82)
    for r in results
        status_str = r.status == "ok" ? "ok" : "FAILED"
        println(
            rpad(r.id, 6),
            rpad(status_str, 9),
            rpad(string(r.runs), 6),
            rpad(_ppb_fmt_elapsed(r.elapsed), 12),
            r.label,
        )
        for msg in r.errors
            println("    Error: ", first(msg, 110))
        end
    end
    println(sep)
    n_ok = count(r -> r.status == "ok", results)
    println("  $(n_ok)/$(length(results)) phases completed successfully.")
    println(sep)
end

function main_paper_benchmarks()
    ppb    = ppb_parse_cli()
    active = _ppb_active_phases(ppb)
    stamp  = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    root   = ppb.dry_run ? joinpath(ppb.outdir, "dry_run_$(stamp)") : joinpath(ppb.outdir, stamp)
    ppb.dry_run || mkpath(root)

    println("[paper-benchmarks] outdir          = $(root)")
    println("[paper-benchmarks] phases           = $(join([p.id for p in active], ", "))")
    println("[paper-benchmarks] thread_ladder    = $(isempty(ppb.threads) ? "auto ($(Sys.CPU_THREADS) CPU threads)" : join(ppb.threads, ","))")
    println("[paper-benchmarks] process_workers  = $(ppb.preview ? "≤$(PPB_PREVIEW_MAX_WORKERS) (preview)" : string(ppb.process_workers))")
    println("[paper-benchmarks] solver_mode      = $(ppb.solver_mode)")
    println("[paper-benchmarks] cpu_pinning      = $(isempty(ppb.cpu_pinning) ? "off" : join(ppb.cpu_pinning, ","))")
    println("[paper-benchmarks] seed             = $(ppb.seed)")
    println("[paper-benchmarks] preview          = $(ppb.preview)")
    println("[paper-benchmarks] dry_run          = $(ppb.dry_run)")
    println()

    results = NamedTuple[]
    for phase in active
        phase_dir = joinpath(root, phase.id)
        println("[phase $(phase.id)] $(phase.label)")
        result = _ppb_run_phase(phase, ppb, phase_dir)
        push!(results, result)
        println("[phase $(phase.id)] done — status=$(result.status)  elapsed=$(_ppb_fmt_elapsed(result.elapsed))")
        println()
    end

    _ppb_print_summary(results)

    if !ppb.dry_run
        println("[paper-benchmarks] Collecting results from $(root) ...")
        raw = _ppb_collect_raw_csvs(root)
        if nrow(raw) > 0
            agg      = _ppb_aggregate(raw)
            raw_path = joinpath(root, "paper_benchmarks_raw_$(stamp).csv")
            agg_path = joinpath(root, "paper_benchmarks_aggregated_$(stamp).csv")
            CSV.write(raw_path, raw)
            CSV.write(agg_path, agg)
            println("[paper-benchmarks] raw CSV        = $(raw_path)")
            println("[paper-benchmarks] aggregated CSV = $(agg_path)")
            plot_paths = _ppb_write_plots(root, agg)
            println("[paper-benchmarks] plots          = $(length(plot_paths))")
            for p in plot_paths
                println("  $(p)")
            end
            try
                report_path = _ppb_write_report(root, agg, active, results, stamp)
                println("[paper-benchmarks] report         = $(report_path)")
            catch err
                @warn "Report generation failed; raw/aggregated CSVs are still on disk" exception=(err, catch_backtrace())
            end
        else
            println("[paper-benchmarks] No raw CSV data found; skipping aggregation and plots.")
        end
    end

    return results
end
