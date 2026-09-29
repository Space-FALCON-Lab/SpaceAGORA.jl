# Per-sample memory of one process-pool worker running P6 trace 6
# (aero_<N>sat_l50_gram_process_100s): what a Distributed pool worker executes
# for every sample it is handed, run in this one process, with its resident
# memory, Julia heap and malloc arena recorded after every sample.
#
# Question it answers: when a pool worker's memory grows across samples, is the
# growth on the Julia GC heap (bounded, if at all, by --heap-size-hint) or
# outside it (native allocations the collector cannot see)? And how many native
# GRAM atmospheres are alive at each point?
#
# Every figure is read, not estimated: VmRSS/VmHWM from /proc/self/status,
# Base.gc_live_bytes() and Base.gc_num() from the runtime, glibc mallinfo2()
# for the malloc arena. The live-atmosphere count is atmospheres created minus
# atmospheres finalized, both counted exactly (a counting finalizer is attached
# to each native atmosphere this sample's run constructs; it does not replace
# the wrapper's own finalizer, which still frees the native object).
#
# Each sample goes through the same steps as ppc_run_sample_once: build the case
# configuration, apply the mode environment, run_simulation. The density-model
# epoch alignment run_simulation starts with is performed here first (the same
# function), so the atmosphere it builds can be counted; run_simulation's own
# alignment then finds the epoch already matching and returns the same object.
#
# Usage (one process, always capped):
#   systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
#     julia --project=. --threads=1 [--heap-size-hint=7.4G] \
#     benchmarks/studies/gram_thread_scaling/worker_growth.jl \
#     --case=aero_4096sat_l50_gram_process_100s --samples=300 --out=results/x.csv \
#     [--definition=new|old] [--model=harness|run_epoch] [--gc-every=K] \
#     [--dump-dir=DIR --dump-samples=8 (0 = every sample)] [--order=forward|reverse]
#
#   --definition=old  the earlier trace-6 sample: the default ppc_spacecraft at
#                     500-550 km with the default 120 km entry interface
#                     (never inside the atmosphere), as worker_rss.jl measured.
#   --model=harness   whatever ppc_gram_atmosphere_model builds (since the fix:
#                     a model at the harness epoch, so no sample rebuilds it).
#   --model=run_epoch pre-seed the harness's GRAM model cache with a model built
#                     at the run's own epoch, so no sample rebuilds it.
#   --model=default_epoch  pre-seed it with a model at GRAMSuite's default
#                     construction epoch, which is what the harness built before
#                     the fix: every sample's run rebuilds it at the run epoch.
#   --gc-every=K      a full collection (GC.gc(true)) after every K-th sample.

using Printf

const GW_DIR = @__DIR__
const GW_REPO_ROOT = normpath(joinpath(GW_DIR, "..", "..", ".."))
const GW_PPC_DIR = joinpath(GW_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

gw_arg(key, default) = begin
    for a in ARGS
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end

const GW_CASE = gw_arg("case", "aero_4096sat_l50_gram_process_100s")
# cases.jl loads GRAMSuite and sets the trace-6 density env from --case on ARGS.
any(a -> startswith(a, "--case="), ARGS) || push!(ARGS, "--case=$(GW_CASE)")

include(joinpath(GW_PPC_DIR, "cli.jl"))
include(joinpath(GW_PPC_DIR, "modes.jl"))
include(joinpath(GW_PPC_DIR, "cases.jl"))
include(joinpath(GW_PPC_DIR, "execution.jl"))

const GW_SAMPLES = parse(Int, gw_arg("samples", "300"))
const GW_OUT = gw_arg("out", joinpath(GW_DIR, "results", "worker_growth.csv"))
const GW_DEFINITION = gw_arg("definition", "new")
const GW_MODEL = gw_arg("model", "harness")
const GW_GC_EVERY = parse(Int, gw_arg("gc-every", "0"))
const GW_DUMP_DIR = gw_arg("dump-dir", "")
const GW_DUMP_SAMPLES = parse(Int, gw_arg("dump-samples", "8"))
const GW_PROFILE = gw_arg("profile", "full")
const GW_MODE = gw_arg("mode", "outer_process")
# forward (1..N) or reverse (N..1): the order samples are run in, so a model
# reused across samples can be checked for any dependence on its predecessor.
const GW_ORDER = gw_arg("order", "forward")

function proc_status_kb(key::String)::Int
    for line in eachline("/proc/self/status")
        startswith(line, key * ":") || continue
        return parse(Int, split(line)[2])
    end
    return -1
end

# glibc mallinfo2: arena, ordblks, smblks, hblks, hblkhd, usmblks, fsmblks,
# uordblks, fordblks, keepcost (all size_t).
function malloc_in_use_bytes()::Tuple{Int, Int}
    mi = ccall(:mallinfo2, NTuple{10, Csize_t}, ())
    arena, hblkhd, uordblks = Int(mi[1]), Int(mi[5]), Int(mi[8])
    return (uordblks + hblkhd, arena + hblkhd)
end

const GW_CREATED = Ref(0)
const GW_FINALIZED = Threads.Atomic{Int}(0)

function count_atmosphere!(model)
    model isa SpaceAGORA.SimulationModel.EnvironmentModels.GRAMAtmosphereModel || return nothing
    atm = model.core.gram_atmosphere
    GW_CREATED[] += 1
    finalizer(_ -> Threads.atomic_add!(GW_FINALIZED, 1), atm)
    return nothing
end

function old_definition_config(cfg::PPCConfig)
    planet = Earth("", PPC_SPICE_PATH)
    return ppc_build_config(
        planet=planet,
        spacecraft=[ppc_spacecraft(planet; id=1)],
        mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=90.0, full=100.0),
        orientation_sim=false,
        dynamic_effectors=(ppc_harmonics_model(planet, 50), AerodynamicCoefficientfM()),
        density_model=ppc_gram_atmosphere_model("earth"),
        dt_max_orbit=5.0
    )
end

sample_config(cfg, i, seed) = GW_DEFINITION == "old" ? old_definition_config(cfg) :
    ppc_single_config(GW_CASE, cfg; seed=seed, mc_index=i)

function dump_history(dir, i, sol)
    mkpath(dir)
    open(joinpath(dir, "history_i$(i).bin"), "w") do io
        write(io, Float64.(sol.t))
        for u in sol.u
            write(io, Float64[x for x in u])
        end
    end
end

function main()
    cfg = PPCConfig(profile=GW_PROFILE, process_workers=1, worker_seed=20260924)
    mode = ppc_mode_specs()[GW_MODE]
    hint = Base.JLOptions().heap_size_hint
    @printf("pid=%d case=%s definition=%s model=%s gc_every=%d samples=%d heap_size_hint=%d\n",
            getpid(), GW_CASE, GW_DEFINITION, GW_MODEL, GW_GC_EVERY, GW_SAMPLES, hint)
    println("density env: FREEZE_PER_STEP=", get(ENV, "SPACEAGORA_DENSITY_FREEZE_PER_STEP", "<unset>"),
            " VACUUM_GRAM_CACHE=", get(ENV, "SPACEAGORA_VACUUM_GRAM_CACHE", "<unset>"))
    flush(stdout)

    if GW_MODEL == "run_epoch"
        epoch = sample_config(cfg, 1, 1).initial_time
        lock(_PPC_GRAM_MODEL_CACHE_LOCK) do
            _PPC_GRAM_MODEL_CACHE["earth"] = GRAMAtmosphereModel(planet_name="earth", initial_time=epoch)
        end
    elseif GW_MODEL == "default_epoch"
        lock(_PPC_GRAM_MODEL_CACHE_LOCK) do
            _PPC_GRAM_MODEL_CACHE["earth"] = GRAMAtmosphereModel(planet_name="earth")
        end
    end
    println("SPACEAGORA_GRAM_NATIVE_COLLECT_LIMIT=", get(ENV, "SPACEAGORA_GRAM_NATIVE_COLLECT_LIMIT", "<unset>"))

    mkpath(dirname(GW_OUT))
    open(GW_OUT, "w") do io
        println(io, "sample,member,wall_s,vm_rss_kb,vm_hwm_kb,gc_live_bytes,gc_total_bytes,gc_pauses,gc_full_sweeps,",
                "malloc_in_use_bytes,malloc_arena_bytes,atmospheres_created,atmospheres_finalized,",
                "heap_size_hint,definition,model,gc_every,gc_time_ns,elapsed_s")
        loop_t0 = time()
        for i in (GW_ORDER == "reverse" ? (GW_SAMPLES:-1:1) : (1:GW_SAMPLES))
            seed = cfg.worker_seed + i
            t0 = time()
            sol = withenv(ppc_mode_env_pairs(mode, cfg)...) do
                args = sample_config(cfg, i, seed)
                aligned = SimulationEngine._with_density_model_epoch(args)
                aligned.environment_model.density_model === args.environment_model.density_model ||
                    count_atmosphere!(aligned.environment_model.density_model)
                ppc_solve_once(aligned, cfg).solution
            end
            wall = time() - t0
            string(sol.retcode) == "Success" || error("sample $i: retcode $(sol.retcode)")
            (GW_DUMP_DIR != "" && (GW_DUMP_SAMPLES <= 0 || i <= GW_DUMP_SAMPLES)) && dump_history(GW_DUMP_DIR, i, sol)
            sol = nothing
            GW_GC_EVERY > 0 && i % GW_GC_EVERY == 0 && GC.gc(true)
            num = Base.gc_num()
            in_use, arena = malloc_in_use_bytes()
            @printf(io, "%d,%d,%.4f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,%s,%d,%d,%.4f\n",
                    i, GW_DEFINITION == "old" ? 0 : mod1(i, 4096), wall,
                    proc_status_kb("VmRSS"), proc_status_kb("VmHWM"),
                    Base.gc_live_bytes(), Base.gc_total_bytes(num), num.pause, num.full_sweep,
                    in_use, arena, GW_CREATED[], GW_FINALIZED[], hint,
                    GW_DEFINITION, GW_MODEL, GW_GC_EVERY, num.total_time, time() - loop_t0)
            flush(io)
            if i == 1 || i % 25 == 0
                @printf("sample %4d  wall %.2f s  RSS %7.1f MB  gc_live %7.1f MB  malloc %7.1f MB  atm live %d\n",
                        i, wall, proc_status_kb("VmRSS") / 1024, Base.gc_live_bytes() / 2^20,
                        in_use / 2^20, GW_CREATED[] - GW_FINALIZED[])
                flush(stdout)
            end
        end
    end
end

main()
