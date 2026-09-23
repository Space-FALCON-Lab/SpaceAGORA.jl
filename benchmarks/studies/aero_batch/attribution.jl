# WS11e deliverable 1: where the time goes in an exponential-atmosphere aero
# constellation, and how much it allocates per spacecraft per derivative
# evaluation.
#
# Serial, one case, one process. Every number this prints is a share of one
# profile or a ratio between two solves timed back to back in the same process.
# It is not a machine benchmark and must not be quoted as one.
#
# The case is built through the parallelization_performance study's own
# builders, so its constellation geometry, tolerances, step cap and mission are
# the catalog's. The vacuum counterpart is solved alongside it as the baseline
# the P6 calibration used, so the atmosphere's marginal cost is a ratio between
# two rungs rather than an absolute.
#
# Usage (one Julia process at a time; check `uptime` first):
#   julia --project=. --threads=1 benchmarks/studies/aero_batch/attribution.jl \
#       --n=1024 --mission=100 [--repeats=1] [--csv=<path>] [--no-allocs]
#
# ROUTE, and why this script takes `--mode` and `--rhs`. At a thread budget of
# one, `_rhs_execution_plan` routes a (harmonics, aero) stack to
# `:satellite_batch`, not to the flat constellation queue: the budget <= 1
# admission to the flat route requires every effector to be served by a
# pre-pass, and aerodynamics is not. So a serial run profiles the per-satellite
# route, which is the right attribution for the serial P6 number, and a
# multi-thread run (`--mode=inner_only --threads=8`) profiles the flat route,
# which is where the pre-passes live. Run both; they are different code.

using Printf
using Profile

const AB_STUDY_DIR = @__DIR__
const AB_REPO_ROOT = normpath(joinpath(AB_STUDY_DIR, "..", "..", ".."))
const AB_PPC_DIR = joinpath(AB_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(AB_PPC_DIR, "cli.jl"))
include(joinpath(AB_PPC_DIR, "modes.jl"))
include(joinpath(AB_PPC_DIR, "cases.jl"))

ab_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end
ab_flag(args, key)::Bool = any(a -> a == "--$key", args)

const AB_N        = parse(Int, ab_arg(ARGS, "n", "1024"))
const AB_MISSION  = parse(Int, ab_arg(ARGS, "mission", "100"))
const AB_REPEATS  = parse(Int, ab_arg(ARGS, "repeats", "1"))
const AB_CSV      = ab_arg(ARGS, "csv", "")
const AB_ALLOCS   = !ab_flag(ARGS, "no-allocs")
const AB_MODE     = ab_arg(ARGS, "mode", "serial")
const AB_RHS      = ab_arg(ARGS, "rhs", "auto")

# Trace 4 of the P6 density ladder at this spacecraft count, and the vacuum
# rung the ladder is read against. The `aero_<N>sat_l50_expatm_<S>s` builder
# branch accepts any N by regex; the catalog registers only 16 and 4096, so the
# builder is called directly rather than through ppc_resolve_cases.
#
# The two rungs are NOT iso-step: the density ladder caps the step at 5 s and
# the vacuum ladder at 20 s, which is the catalog's choice and not this
# script's. So the wall ratio between them is "what the atmosphere costs,
# including the steps it forces", and the per-RHS-evaluation ratio printed
# beside it is the step-count-normalized form. Read both.
const AB_AERO_CASE   = "aero_$(AB_N)sat_l50_expatm_$(AB_MISSION)s"
const AB_VACUUM_CASE = "gravity_$(AB_N)sat_l50_vacuum_$(AB_MISSION)s"

# Nested phase markers. A sample is attributed to the INNERMOST marker on its
# backtrace, so an outer phase's share is the work it does outside every inner
# phase it contains -- `flat_effector_slots` is the queue/pre-pass bookkeeping
# that is neither harmonics nor aero, not the sum of everything under it.
#
# Order in this vector is documentation only; attribution is by stack depth.
const AB_MARKERS = [
    ("aero_hart_coefficients",   ("aerodynamic_coefficient_fM",)),
    ("aero_wrench",              ("_aero_pure_wrench", "wrench_caching!")),
    ("aero_prepass",             ("_accumulate_aero_flat_batch!",)),
    ("flat_queue_dispatch",      ("_evaluate_dynamic_effector",)),
    ("harmonics_prepass",        ("_accumulate_harmonics_flat_batch!", "_harmonics_flat_batch_kernel!")),
    ("flat_reduce",              ("_reduce_flat_effector_slots_range!",)),
    ("flat_effector_slots",      ("_accumulate_dynamic_effectors_flat_slots!", "_accumulate_dynamic_effectors_flat_batch!")),
    ("rhs_atmosphere_prefill",   ("_prefill_environment_samples!",)),
    ("density_callback",         ("update_density_sat!", "getDensityBatch!", "_gram_isolated_pool_batch_eval!")),
    ("rhs_other",                ("_spacecraft_dynamics_flat_constellation_effector_queue!",
                                  "_spacecraft_dynamics_satellite_batch!")),
]

function ab_marker_for(name::String)::Union{Nothing, String}
    for (bucket, funcs) in AB_MARKERS
        for f in funcs
            occursin(f, name) && return bucket
        end
    end
    return nothing
end

# Split the raw profile buffer into one backtrace per sample. Metadata is
# excluded at fetch time so the buffer is plain instruction pointers separated
# by a zero word, which is the layout this splitter assumes.
function ab_backtraces(data::Vector{UInt})
    traces = Vector{Vector{UInt}}()
    current = UInt[]
    for ip in data
        if ip == 0
            isempty(current) || push!(traces, reverse(current))  # reverse: root first
            current = UInt[]
        else
            push!(current, ip)
        end
    end
    isempty(current) || push!(traces, reverse(current))
    return traces
end

function ab_profile_shares(data::Vector{UInt}, lidict)
    counts = Dict{String, Int}()
    total = 0
    for trace in ab_backtraces(data)
        total += 1
        bucket = "solver_and_other"
        for ip in trace
            frames = get(lidict, ip, nothing)
            frames === nothing && continue
            # `frames` holds the inlined frames for one instruction pointer,
            # leaf first; reversed here so the whole scan runs root to leaf and
            # the last marker it sees is the innermost one.
            for frame in Iterators.reverse(frames)
                m = ab_marker_for(String(frame.func))
                m === nothing || (bucket = m)   # innermost marker wins
            end
        end
        counts[bucket] = get(counts, bucket, 0) + 1
    end
    return counts, total
end

function ab_solve(case_name::String, profile::String)
    cfg = PPCConfig(profile=profile, solver_mode="auto_stiff")
    args = ppc_single_config(case_name, cfg)
    timed = @timed SimulationEngine.run_simulation(
        args; isolate_state=false, return_solution=true, return_solver_metadata=true
    )
    r = timed.value
    sol = r.solution
    geti(name) = try Int(getproperty(sol.stats, name)) catch; -1 end
    return (
        wall_s=Float64(timed.time), gc_s=Float64(timed.gctime), bytes=Int(timed.bytes),
        nf=geti(:nf), naccept=geti(:naccept), nreject=geti(:nreject),
        steps=length(sol.t), nsats=length(args.dynamics_model.spacecraft),
        retcode=string(sol.retcode),
    )
end

function main()
    mode = ppc_mode_specs()[AB_MODE]
    cfg = PPCConfig(profile="full", solver_mode="auto_stiff")
    envpairs = copy(ppc_mode_env_pairs(mode, cfg; outer_tasks=1))
    if AB_RHS != "auto"
        push!(envpairs, "SPACEAGORA_RHS_EXECUTION_MODE" => AB_RHS)
    end

    @printf("host=%s threads=%d mode=%s rhs_mode=%s N=%d mission=%ds aero_case=%s vacuum_case=%s\n",
            gethostname(), Threads.nthreads(), AB_MODE, AB_RHS,
            AB_N, AB_MISSION, AB_AERO_CASE, AB_VACUUM_CASE)
    println("load-at-start: ", strip(read(`uptime`, String)))
    flush(stdout)

    for case in (AB_AERO_CASE, AB_VACUUM_CASE)
        withenv(envpairs...) do
            ab_solve(case, "test")
        end
    end
    println("warmed both cases"); flush(stdout)

    rows = Dict{String, Vector{NamedTuple}}()
    for rep in 1:AB_REPEATS, case in (AB_AERO_CASE, AB_VACUUM_CASE)
        r = withenv(envpairs...) do
            ab_solve(case, "full")
        end
        push!(get!(rows, case, NamedTuple[]), r)
        @printf("rep%d %-44s N=%5d wall=%9.3fs gc=%6.3fs alloc=%8.3fGiB nf=%8d acc=%6d rej=%5d %s\n",
                rep, case, r.nsats, r.wall_s, r.gc_s, r.bytes / 2^30,
                r.nf, r.naccept, r.nreject, r.retcode)
        flush(stdout)
    end

    aero = rows[AB_AERO_CASE][1]
    vac = rows[AB_VACUUM_CASE][1]
    aero_w = minimum(r.wall_s for r in rows[AB_AERO_CASE])
    vac_w = minimum(r.wall_s for r in rows[AB_VACUUM_CASE])
    println("\n-- exponential aero against the vacuum rung (same process, back to back) --")
    @printf("wall x%.3f   nf x%.3f   per-RHS-eval x%.3f   marginal wall %.3f s (%.1f%% of the aero case)\n",
            aero_w / vac_w, aero.nf / max(1, vac.nf),
            (aero_w / max(1, aero.nf)) / (vac_w / max(1, vac.nf)),
            aero_w - vac_w, 100 * (aero_w - vac_w) / aero_w)

    println("\n-- allocation per spacecraft per derivative evaluation --")
    for (label, r) in (("exponential aero", aero), ("vacuum", vac))
        denom = max(1, r.nf * r.nsats)
        @printf("%-18s %10.1f B/sat/eval   (total %.3f GiB, nf=%d, N=%d)\n",
                label, r.bytes / denom, r.bytes / 2^30, r.nf, r.nsats)
    end

    # Time profile of the aero case only. delay is deliberately short and the
    # buffer large: one solve of a 1024-spacecraft case is tens of thousands of
    # derivative evaluations and the buckets below are small.
    Profile.clear()
    Profile.init(n = 40_000_000, delay = 0.0005)
    withenv(envpairs...) do
        @profile ab_solve(AB_AERO_CASE, "full")
    end
    data = Profile.fetch(include_meta = false)
    lidict = Profile.getdict(data)
    counts, total = ab_profile_shares(data, lidict)
    println("\n-- time profile of $(AB_AERO_CASE), by innermost phase marker ($(total) samples) --")
    for (bucket, n) in sort(collect(counts); by = kv -> -last(kv))
        @printf("%-26s %7d  %6.2f%%\n", bucket, n, 100 * n / max(1, total))
    end
    profile_path = joinpath(AB_STUDY_DIR, "profile_$(AB_AERO_CASE)_$(AB_MODE)_$(Threads.nthreads())t.txt")
    open(profile_path, "w") do io
        Profile.print(IOContext(io, :displaysize => (24, 2000)); format=:flat, C=true, sortedby=:count, mincount=20)
    end
    println("flat profile written: ", profile_path)

    if AB_ALLOCS
        # Allocation attribution. sample_rate below 1 keeps the recorder's own
        # overhead from dominating a solve that allocates tens of GiB; the
        # per-site SHARES are what this is read for, not the absolute bytes.
        Profile.Allocs.clear()
        withenv(envpairs...) do
            Profile.Allocs.@profile sample_rate = 0.0005 ab_solve(AB_AERO_CASE, "full")
        end
        allocs = Profile.Allocs.fetch()
        sites = Dict{String, Tuple{Int, Int}}()
        for a in allocs.allocs
            frames = a.stacktrace
            isempty(frames) && continue
            key = string(frames[1].func, " @ ", basename(String(frames[1].file)), ":", frames[1].line)
            n, b = get(sites, key, (0, 0))
            sites[key] = (n + 1, b + a.size)
        end
        total_b = sum(last(v) for v in values(sites); init = 0)
        println("\n-- allocation sites of $(AB_AERO_CASE) (sampled at 5e-4; shares, not totals) --")
        for (key, v) in first(sort(collect(sites); by = kv -> -last(last(kv))), 25)
            @printf("%9d  %10.2f MiB  %6.2f%%  %s\n",
                    v[1], v[2] / 2^20, 100 * v[2] / max(1, total_b), key)
        end
    end

    if !isempty(AB_CSV)
        mkpath(dirname(AB_CSV))
        open(AB_CSV, "w") do io
            println(io, "host,threads,mode,rhs_mode,case,rep,nsats,wall_s,gc_s,bytes,nf,naccept,nreject,steps,retcode")
            for case in (AB_AERO_CASE, AB_VACUUM_CASE)
                for (rep, r) in enumerate(rows[case])
                    @printf(io, "%s,%d,%s,%s,%s,%d,%d,%.6f,%.6f,%d,%d,%d,%d,%d,%s\n",
                            gethostname(), Threads.nthreads(), AB_MODE, AB_RHS, case, rep, r.nsats,
                            r.wall_s, r.gc_s, r.bytes, r.nf, r.naccept, r.nreject,
                            r.steps, r.retcode)
                end
            end
        end
        println("\ncsv written: ", AB_CSV)
    end

    println("load-at-end: ", strip(read(`uptime`, String)))
end

main()
