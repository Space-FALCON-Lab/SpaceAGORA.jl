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
    # Two points (1 and the machine max) rather than the full ladder: B6
    # covers many workload axes at once (case count x mode count already
    # multiplies out large), so it takes a coarse cut of the thread-budget
    # axis specifically to keep total cost bounded -- the *full* ladder is
    # already covered in depth by B1/B2's dedicated thread-scaling phases.
    phase.thread_mode == :low_high && return sort(unique([1, maximum(base)]))
    return base
end

# Caps a phase's MC sample ladder at ppb.mc_samples_max (e.g. B4's default
# [1, 4, 16, 64, 256, 1024]), keeping at least the smallest value if the cap
# would otherwise remove everything. No-op when mc_samples_max is unset.
function _ppb_capped_mc_samples(phase::PPBPhase, ppb::PPBConfig)::Vector{Int}
    ppb.mc_samples_max === nothing && return phase.mc_samples
    samples = filter(s -> s <= ppb.mc_samples_max, phase.mc_samples)
    isempty(samples) && (samples = [minimum(phase.mc_samples)])
    return samples
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
        mc_samples      = _ppb_capped_mc_samples(phase, ppb),
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
    println("[dry-run]   mc_samples      = $(_ppb_capped_mc_samples(phase, ppb))")
    println("[dry-run]   process_workers = $(process_workers)")
    println("[dry-run]   repeats         = $(phase.repeats), warmup = $(phase.warmup)")
end

function _ppb_run_phase(
    phase::PPBPhase, ppb::PPBConfig, phase_dir::String;
    on_run_complete::Union{Nothing, Function}=nothing,
)::NamedTuple
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
                ppc_run_controller(_ppb_build_ppc_config(phase, ppb, phase_dir); on_run_complete)
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
                    ppc_run_controller(cfg; on_run_complete)
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

# Rebuilds the combined raw/aggregated CSVs from whatever is on disk right
# now (finished per-phase CSVs plus in-progress worker_rows for anything
# still running) and overwrites the canonical CSV files with that snapshot.
# Called after every individual worker run completes, so a hang or crash
# partway through a multi-hour/multi-day phase still leaves up-to-date CSVs
# for everything that finished so far — not just whatever was captured the
# last time the whole suite completed. Deliberately CSV-only: plots/report
# regeneration is comparatively expensive (Plots.jl rendering) and only
# happens once, at the end of the full run (see main_paper_benchmarks).
#
# Wrapped in a top-level try/catch because this runs unattended, potentially
# thousands of times over a multi-day run; a transient hiccup here (e.g.
# reading a file mid-write) must never take down the benchmark run itself.
function _ppb_snapshot_csvs!(root::String, stamp::String)
    try
        collected = _ppb_collect_partial(root)
        nrow(collected.raw) == 0 && return nothing
        agg = _ppb_aggregate(collected.raw)
        nrow(agg) == 0 && return nothing

        CSV.write(joinpath(root, "paper_benchmarks_raw_$(stamp).csv"), collected.raw)
        CSV.write(joinpath(root, "paper_benchmarks_aggregated_$(stamp).csv"), agg)
    catch err
        @warn "Incremental CSV snapshot failed; will retry after the next run" exception=(err, catch_backtrace())
    end
    return nothing
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

    # --resume points at an existing run directory (e.g. one a prior invocation
    # aborted mid-phase) rather than starting a fresh timestamped one. Reuse
    # that directory's own stamp for the final CSV/report names when it
    # matches the expected pattern, so a resumed run's outputs read as a
    # continuation of the original rather than a separate run; fall back to a
    # fresh stamp otherwise (e.g. a hand-picked directory name).
    resuming = !isempty(ppb.resume)
    stamp = if resuming
        m = match(r"^\d{8}_\d{6}", basename(ppb.resume))
        m === nothing ? Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS") : String(m.match)
    else
        Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    end
    root = if resuming
        ppb.resume
    elseif ppb.dry_run
        joinpath(ppb.outdir, "dry_run_$(stamp)")
    else
        joinpath(ppb.outdir, stamp)
    end
    ppb.dry_run || mkpath(root)

    resuming && println("[paper-benchmarks] resuming        = $(root)")
    println("[paper-benchmarks] outdir          = $(root)")
    println("[paper-benchmarks] phases           = $(join([p.id for p in active], ", "))")
    println("[paper-benchmarks] thread_ladder    = $(isempty(ppb.threads) ? "auto ($(Sys.CPU_THREADS) CPU threads)" : join(ppb.threads, ","))")
    println("[paper-benchmarks] process_workers  = $(ppb.preview ? "≤$(PPB_PREVIEW_MAX_WORKERS) (preview)" : string(ppb.process_workers))")
    println("[paper-benchmarks] mc_samples_max   = $(ppb.mc_samples_max === nothing ? "unset" : (ppb.preview ? "$(min(ppb.mc_samples_max, PPB_PREVIEW_MAX_SAMPLES)) (preview-capped)" : string(ppb.mc_samples_max)))")
    println("[paper-benchmarks] solver_mode      = $(ppb.solver_mode)")
    println("[paper-benchmarks] cpu_pinning      = $(isempty(ppb.cpu_pinning) ? "off" : join(ppb.cpu_pinning, ","))")
    println("[paper-benchmarks] seed             = $(ppb.seed)")
    println("[paper-benchmarks] preview          = $(ppb.preview)")
    println("[paper-benchmarks] dry_run          = $(ppb.dry_run)")
    println()

    results = NamedTuple[]
    on_run_complete = ppb.dry_run ? nothing : () -> _ppb_snapshot_csvs!(root, stamp)

    for phase in active
        phase_dir = joinpath(root, phase.id)
        println("[phase $(phase.id)] $(phase.label)")
        result = _ppb_run_phase(phase, ppb, phase_dir; on_run_complete)
        push!(results, result)
        println("[phase $(phase.id)] done — status=$(result.status)  elapsed=$(_ppb_fmt_elapsed(result.elapsed))")
        println()
    end

    _ppb_print_summary(results)

    if !ppb.dry_run
        println("[paper-benchmarks] Collecting final results from $(root) ...")
        collected = _ppb_collect_partial(root)
        if nrow(collected.raw) > 0
            agg      = _ppb_aggregate(collected.raw)
            raw_path = joinpath(root, "paper_benchmarks_raw_$(stamp).csv")
            agg_path = joinpath(root, "paper_benchmarks_aggregated_$(stamp).csv")
            CSV.write(raw_path, collected.raw)
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
