# Validate the paper harness's precompile workload on representative points:
# for each point, launch the harness worker twice, back to back -- once as the
# controller launches it without the workload ("stock") and once with it
# ("workload") -- and record
#
#   import_s        process start -> harness files and packages loaded
#   first_solve_s   process start -> the coordinator's first solve returns
#   first_timed_s   process start -> the first timed repeat starts (after the
#                   warm-ups and, for a campaign, the pool's provisioning)
#   total_s         the worker process's wall time
#   timed_median_s  median of the point's timed repeats (the harness row's wall_time_s)
#
# and whether every solve's step times and final state are byte-identical
# between the two runs (raw Float64 dumps; see SPACEAGORA_PPC_DUMP_STATE_DIR in
# parallelization_performance/execution.jl).
#
# Usage (from the repository root, after build_workload.sh):
#   julia --project=. --startup-file=no \
#     benchmarks/studies/paper_parallelization_benchmarks/workload/validate_workload.jl <out_dir> [--points=a,b] [--pairs=N]
#
# --pairs=N launches each point N times per variant, alternating which variant
# goes first (stock/workload, then workload/stock, ...), because a timed median
# moves by several percent between two launches of the same process, and a
# single back-to-back pair cannot separate that from an effect of the image.
#
# One Julia worker runs at a time. Numbers are this machine's; nothing here is
# a threshold.

include(joinpath(@__DIR__, "..", "..", "paper_parallelization_benchmarks.jl"))

const PPBV_OUT = abspath(get(ARGS, 1, joinpath(PPC_REPO_ROOT, "output", "paper_workload", "validate")))
const PPBV_ONLY = let a = findfirst(s -> startswith(s, "--points="), ARGS)
    a === nothing ? String[] : String.(split(split(ARGS[a], "=", limit=2)[2], ","))
end
const PPBV_PAIRS = let a = findfirst(s -> startswith(s, "--pairs="), ARGS)
    a === nothing ? 1 : parse(Int, split(ARGS[a], "=", limit=2)[2])
end

# Representative points: one per launch shape the harness uses. Thread and
# worker counts are capped at 8 so a point fits this workstation's measurement
# rules; the case, mode, sample count, repeats and warm-up are the phase's own.
const PPBV_POINTS = [
    (id="p1_1sat_serial",   phase="P1", case=_ppb_paper_size_case(1),    mode="serial",
     threads=8, workers=1, samples=1,   repeats=3, warmup=1),
    (id="p2_4096_inner",    phase="P2", case=_ppb_paper_size_case(4096), mode="inner_only",
     threads=8, workers=1, samples=1,   repeats=3, warmup=1),
    (id="p3_process",       phase="P3", case="independent_1sat_1hr",    mode="outer_process",
     threads=4, workers=4, samples=256, repeats=5, warmup=1),
    (id="p3_policy_v2",     phase="P3", case="independent_1sat_1hr",    mode="policy_v2",
     threads=4, workers=4, samples=256, repeats=5, warmup=1),
    (id="p5_split_static",  phase="P5", case="mcgrid_16sat_8mc",        mode="outer_inner_static",
     threads=4, workers=2, samples=8,   repeats=5, warmup=1),
]

function ppbv_trace(log::String)
    events = Dict{Tuple{String, Int}, Float64}()   # (label, myid) -> first epoch
    for line in eachline(log)
        m = match(r"PPC_STARTUP_TRACE (\S+) (\d+) (\d+) ([0-9.eE+-]+)", line)
        m === nothing && continue
        key = (String(m.captures[1]), parse(Int, m.captures[3]))
        haskey(events, key) || (events[key] = parse(Float64, m.captures[4]))
    end
    return events
end

function ppbv_dumps_identical(a::String, b::String)
    fa = isdir(a) ? sort(readdir(a)) : String[]
    fb = isdir(b) ? sort(readdir(b)) : String[]
    (isempty(fa) || fa != fb) && return (false, length(fa), length(fb))
    same = all(f -> read(joinpath(a, f)) == read(joinpath(b, f)), fa)
    return (same, length(fa), length(fb))
end

function ppbv_run(p, variant::String, workload_env, pair::Int)
    dir = joinpath(PPBV_OUT, p.id, "pair$(pair)", variant)
    rm(dir; recursive=true, force=true)
    mkpath(dir)
    cfg = PPCConfig(profile="full", outdir=dir, process_workers=p.workers, warmup=p.warmup,
                    repeats=p.repeats, solver_mode=PPBConfig().solver_mode, parity_samples=512)
    outfile = joinpath(dir, "rows.csv")
    cmd = ppc_worker_cmd(cfg; case=p.case, mode=p.mode, threads=p.threads, repeat=1,
                         repeats=p.repeats, seed=cfg.seed + 1, mc_samples=p.samples,
                         outfile=outfile, parity=false,
                         workload_env=variant == "workload" ? workload_env : nothing)
    cmd = addenv(cmd, "SPACEAGORA_PPC_STARTUP_TRACE" => "1",
                 "SPACEAGORA_PPC_DUMP_STATE_DIR" => joinpath(dir, "states"))
    log = joinpath(dir, "worker.log")
    println("[validate] $(p.id) pair $(pair) $(variant): $(p.case) mode=$(p.mode) threads=$(p.threads) workers=$(p.workers)")
    t0 = time()
    ok = success(pipeline(cmd; stdout=log, stderr=log))
    total = time() - t0
    ev = ppbv_trace(log)
    rel(label) = haskey(ev, (label, 1)) ? round(ev[(label, 1)] - t0; digits=2) : missing
    pool_first = [v for ((l, id), v) in ev if l == "first_solve" && id != 1]
    rows = isfile(outfile) ? CSV.read(outfile, DataFrame) : DataFrame()
    walls = nrow(rows) > 0 ? Float64.(rows.wall_time_s) : Float64[]
    return (
        point=p.id, pair=pair, phase=p.phase, case=p.case, mode=p.mode, threads=p.threads,
        process_workers=p.workers, mc_samples=p.samples, variant=variant, success=ok,
        import_s=rel("loaded"), first_solve_s=rel("first_solve"), first_timed_s=rel("timed_start"),
        pool_first_solve_s=isempty(pool_first) ? missing : round(maximum(pool_first) - t0; digits=2),
        total_s=round(total; digits=2),
        timed_median_s=isempty(walls) ? missing : median(walls),
        timed_min_s=isempty(walls) ? missing : minimum(walls),
        repeats=length(walls),
        timed_walls_s=join(string.(walls), ";"),
    )
end

function main_validate()
    workload_env = withenv("SPACEAGORA_PPB_WORKLOAD" => "1") do
        ppc_workload_env()
    end
    mkpath(PPBV_OUT)
    points = isempty(PPBV_ONLY) ? PPBV_POINTS : filter(p -> p.id in PPBV_ONLY, PPBV_POINTS)
    # Bring both images into the page cache once, so neither variant's first
    # point pays a cold disk read the other does not.
    julia_bin = Base.julia_cmd().exec[1]
    run(addenv(`$(julia_bin) --startup-file=no --project=$(PPC_REPO_ROOT) -e "using SpaceAGORA; import SpaceAGORAPaperWorkload"`,
               "JULIA_LOAD_PATH" => ppc_workload_load_path(workload_env)))
    rows = NamedTuple[]
    for p in points, pair in 1:PPBV_PAIRS
        if isodd(pair)
            stock = ppbv_run(p, "stock", workload_env, pair)
            wl = ppbv_run(p, "workload", workload_env, pair)
        else
            wl = ppbv_run(p, "workload", workload_env, pair)
            stock = ppbv_run(p, "stock", workload_env, pair)
        end
        same, na, nb = ppbv_dumps_identical(joinpath(PPBV_OUT, p.id, "pair$(pair)", "stock", "states"),
                                            joinpath(PPBV_OUT, p.id, "pair$(pair)", "workload", "states"))
        for r in (stock, wl)
            push!(rows, merge(r, (states_identical=same, state_files=r.variant == "stock" ? na : nb)))
        end
        println("[validate] $(p.id) pair $(pair): import $(stock.import_s) -> $(wl.import_s) s, first solve $(stock.first_solve_s) -> $(wl.first_solve_s) s, " *
                "first timed $(stock.first_timed_s) -> $(wl.first_timed_s) s, total $(stock.total_s) -> $(wl.total_s) s, " *
                "timed median $(stock.timed_median_s) -> $(wl.timed_median_s) s, states identical=$(same) ($(na) files)")
        CSV.write(joinpath(PPBV_OUT, "workload_validation.csv"), DataFrame(rows))
    end
    return rows
end

main_validate()
