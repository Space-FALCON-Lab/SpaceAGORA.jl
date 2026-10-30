# Probe GC driver: ONE run_monte_carlo campaign per call, ONE worker (timing), process route asserted. Memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. probe_gc_run.jl --phase=plan|pilot|calls|main --results-root=DIR [--horizon-days=D] [--host=NAME] [--workers=1]
#   plan   print the flight list (no Julia flights, loads nothing heavy)
#   pilot  1-day warm-up + 3-day direct 20x20 perturbed flight (wall seconds per day -> horizon rule of the plan)
#   calls  untimed 3-day direct 20x20 flights (mean, perturbed) with the GRAM profile counter on (native mean-density call counts)
#   main   warm-up + 18 flights (order A, reverse order B, two alternative-seed flights), horizon from --horizon-days
# Reads the Stage I calibration from <root>/stagei/mars/{cells_base.csv,cal/<tag>}. Writes <root>/probe_gc/{runs/<run_id>/,summary_<phase>_<host>.csv,mc_<phase>_<host>.toml}.
# TIMING JOBS MUST RUN ALONE on the host. Plan: tasks/2026-10-15_periapsis-region-border-targeting/probe_gc_plan.md.
const PHASE_ARG = let m = [match(r"^--phase=(.*)$", a) for a in ARGS]; i = findfirst(!isnothing, m); i === nothing ? "" : m[i].captures[1]; end
PHASE_ARG == "plan" || include(joinpath(@__DIR__, "probe_gc_common.jl"))
using TOML
const CFGI = joinpath(@__DIR__, "stagei_config.toml")
const TAG = "I0252_R000_W090_E1"

function parse_cli()
    kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2]); end
    return kv
end
const CLI = parse_cli(); const ROOT = abspath(get(CLI, "results-root", ".")); const PG = joinpath(ROOT, "probe_gc"); const HOST = get(CLI, "host", "local")
const CELL = (inc=25.2, raan=0.0, omega=90.0, epoch_idx=1)   # SOURCED: <root>/stagei/mars/cells_base.csv row I0252_R000_W090_E1

mkseed(id, cfg, deg, arm, days; rep=1, stag=TAG, profile=false, warm=false) = (; config_path=CFGI, run_id=id, cfg, degree=deg, arm, seed_tag=stag, tag=TAG, inc=CELL.inc, raan=CELL.raan, omega=CELL.omega,
    epoch_idx=CELL.epoch_idx, horizon_days=days, out=joinpath(PG, "runs"), cal_root=joinpath(ROOT, "stagei", "mars", "cal"), profile, warm, rep)
rid(cfg, deg, arm, rep, sfx="") = "$(cfg)$(deg)_$(arm)_r$(rep)$(sfx)"
const CONFIGS = [("direct", 20), ("cache", 20), ("direct", 10), ("cache", 10)]

function flights(phase, days)
    if phase == "pilot"; return [mkseed("warm", "direct", 20, "pert", 1.0; warm=true), mkseed("pilot_direct20_pert", "direct", 20, "pert", 3.0)]; end
    if phase == "calls"; return [mkseed("calls_direct20_$(a)", "direct", 20, a, 3.0; profile=true) for a in ("mean", "pert")]; end
    phase == "main" || error("unknown phase $phase")
    A = [(c, d, a) for a in ("mean", "pert") for (c, d) in CONFIGS]
    out = [mkseed("warm", "direct", 20, "pert", 1.0; warm=true)]
    for (c, d, a) in A; push!(out, mkseed(rid(c, d, a, 1), c, d, a, days; rep=1)); end
    for (c, d, a) in reverse(A); push!(out, mkseed(rid(c, d, a, 2), c, d, a, days; rep=2)); end
    for c in ("direct", "cache"); push!(out, mkseed(rid(c, 20, "pert", 1, "_alt"), c, 20, "pert", days; rep=1, stag=TAG * "|alt")); end
    return out
end

phase = PHASE_ARG; days = parse(Float64, get(CLI, "horizon-days", "30"))
seeds = flights(phase == "plan" ? "main" : phase, days)
if phase == "plan"
    foreach(s -> println(s.run_id, "  cfg=", s.cfg, " deg=", s.degree, " arm=", s.arm, " seed_tag=", s.seed_tag, " days=", s.horizon_days, " warm=", s.warm), seeds)
    println(length(seeds), " flights incl. warm-up; calls phase: ", length(flights("calls", 3.0)), "; pilot: ", length(flights("pilot", 3.0))); exit(0)
end
# One worker, one flight at a time. run_monte_carlo with workers=1 plans route :none on this host class (observed error: "GRAM campaign ran on route none, not :process"), and workers=2 would
# overlap two flights, which a timing run must not. So the flights go one by one to the single process worker of SpaceAGORA's own pool (--threads=1, its own native GRAM), the same pool
# run_gram_campaign provisions. Recorded in the plan as the single-worker fallback.
function sequential_single_worker(task, seeds; worker_files)
    pool = campaign_process_pool(); ids = ensure_process_workers!(pool, 1); length(ids) == 1 || error("expected one process worker, got $(length(ids))")
    Distributed.remotecall_eval(Main, ids, :(ENV["PC_SWEEP_WORKER"] = "1"))
    for f in worker_files; Distributed.remotecall_eval(Main, ids, :(include($f))); end
    t0 = time(); samples = NamedTuple[]; failed = NamedTuple[]
    for s in seeds
        try; push!(samples, (; seed=s, value=Distributed.remotecall_fetch(task, ids[1], s)))
        catch err; push!(failed, (; seed=s, error=err)); end
    end
    return (; route=:process, threads=1, samples, successful=samples, failed, elapsed_s=time() - t0)
end
res = sequential_single_worker(probe_gc_task, seeds; worker_files=[joinpath(@__DIR__, "probe_gc_common.jl")])
println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1))")
foreach(f -> println("FAILED $(f.seed.run_id): $(first(sprint(showerror, f.error), 600))"), res.failed)
rows = [merge(s.value, (; host=HOST)) for s in res.successful]
write_rows_m(joinpath(PG, "summary_$(phase)_$(HOST).csv"), rows)
open(joinpath(PG, "mc_$(phase)_$(HOST).toml"), "w") do io
    println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nhost = \"$HOST\"\nuptime = \"$(strip(read(`uptime`, String)))\"\ncpu_model = \"$(try strip(split(first(filter(l -> startswith(l, "model name"), readlines("/proc/cpuinfo"))), ":"; limit=2)[2]) catch; "unknown" end)\"\nlogical_cpus = $(Sys.CPU_THREADS)")
end
res.route === :process || error("route was $(res.route), not :process")
for r in rows; println(rpad(r.run_id, 28), " ", r.status, " passes=", r.passes, " wall_s=", round(r.run_wall_s; digits=1), " ", r.error); end
