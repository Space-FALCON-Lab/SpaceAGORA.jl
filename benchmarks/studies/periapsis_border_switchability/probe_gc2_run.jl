# Probe GC-2 driver (plan: tasks/2026-10-15_periapsis-region-border-targeting/probe_gc2_plan.md). One Julia job per phase, memory capped:
#   systemd-run --user --scope -p MemoryMax=11G julia --project=. probe_gc2_run.jl --phase=plan|smoke|states|r1|r4|ens|timed --results-root=DIR [--workers=3] [--host=NAME]
#   smoke   3-day states flight, 2-day R1 sampler flight, one restart task (plumbing check, own output folder)
#   states  direct 20x20 mean-arm trajectory, 12 days -> <root>/probe_gc2/states.csv (apoapsis states)
#   r1      cache 20x20 and 10x10 mean-arm 30-day flights with the live-spline sampler
#   r4      20 restart tasks (4 configurations x (mean + 3 seeds), one orbit ahead)
#   ens     16 seeds x 4 configurations, perturbed, 30 days, plus the 4 mean-arm flights (up to 3 workers; NOT timing)
#   timed   single worker, alone: direct, cache, cache, direct at seed k = 0 (the only timing measurements)
# Reads the Stage I calibration from <root>/stagei/mars/cal/<tag>. Writes <root>/probe_gc2/.
const PHASE_ARG = let m = [match(r"^--phase=(.*)$", a) for a in ARGS]; i = findfirst(!isnothing, m); i === nothing ? "" : m[i].captures[1]; end
PHASE_ARG == "plan" || include(joinpath(@__DIR__, "probe_gc2_common.jl"))
using TOML
const CFGI = joinpath(@__DIR__, "stagei_config.toml"); const TAG = "I0252_R000_W090_E1"
function parse_cli()
    kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2]); end
    return kv
end
const CLI = parse_cli(); const ROOT = abspath(get(CLI, "results-root", ".")); const HOST = get(CLI, "host", "local")
const PG = joinpath(ROOT, PHASE_ARG == "smoke" ? "probe_gc2_smoke" : "probe_gc2")
const CELL = (inc=25.2, raan=0.0, omega=90.0, epoch_idx=1)
const N_ENS = 16; const N_RESTARTS = 20; const CFGS = (("direct", 20), ("cache", 20), ("direct", 10), ("cache", 10))

mk(id, cfg, deg, arm, days; stag=TAG, extra=(;)) = merge((; config_path=CFGI, run_id=id, cfg, degree=deg, arm, seed_tag=stag, tag=TAG, inc=CELL.inc, raan=CELL.raan, omega=CELL.omega, epoch_idx=CELL.epoch_idx,
    horizon_days=days, out=joinpath(PG, "runs"), cal_root=joinpath(ROOT, "stagei", "mars", "cal"), keep_log=false), extra)
ek(k) = "|e" * lpad(k, 2, '0')

states_seed(days) = merge(mk("states_run", "direct", 20, "mean", days), (; out=PG))
r1_seeds(days) = [mk("dens_$(c)$(d)_mean", c, d, "mean", days; extra=(; sampler=true)) for (c, d) in (("cache", 20), ("cache", 10))]
restart_seeds(n) = [(; restart_i=i, states_path=joinpath(PG, "states.csv"), out=joinpath(PG, "runs"), cal_root=joinpath(ROOT, "stagei", "mars", "cal"), config_path=CFGI, tag=TAG, inc=CELL.inc, raan=CELL.raan,
    omega=CELL.omega, epoch_idx=CELL.epoch_idx, n_seeds=3) for i in 1:n]
function ens_seeds()
    out = [mk("m_$(c)$(d)_mean", c, d, "mean", 30.0) for (c, d) in CFGS]
    for k in 0:(N_ENS - 1), (c, d) in CFGS; push!(out, mk("e$(lpad(k, 2, '0'))_$(c)$(d)_pert", c, d, "pert", 30.0; stag=ek(k), extra=(; keep_log=true))); end
    return out
end
timed_seeds() = [mk("warm", "direct", 20, "pert", 1.0; stag=ek(0)); [mk("timed_$(i)_$(c)20_pert", c, 20, "pert", 30.0; stag=ek(0)) for (i, c) in enumerate(("direct", "cache", "cache", "direct"))]]

if PHASE_ARG == "plan"
    println("states 12 d; r1 ", length(r1_seeds(30.0)), " flights; r4 ", N_RESTARTS, " restart tasks x 16 one-orbit flights = ", N_RESTARTS * 16, "; ens ", length(ens_seeds()), " flights; timed ", length(timed_seeds()), " (incl. warm-up)")
    exit(0)
end

function sequential_single_worker(task, seeds; worker_files)
    pool = campaign_process_pool(); ids = ensure_process_workers!(pool, 1)
    Distributed.remotecall_eval(Main, ids, :(ENV["PC_SWEEP_WORKER"] = "1"))
    for f in worker_files; Distributed.remotecall_eval(Main, ids, :(include($f))); end
    t0 = time(); samples = NamedTuple[]; failed = NamedTuple[]
    for s in seeds
        try; push!(samples, (; seed=s, value=Distributed.remotecall_fetch(task, ids[1], s)))
        catch err; push!(failed, (; seed=s, error=err)); end
    end
    return (; route=:process, threads=1, samples, successful=samples, failed, elapsed_s=time() - t0)
end
VARIANTS = [("cache", 20, 40), ("cache", 20, 80)]
r4b_seeds() = [merge(x, (; n_orbits=2, configs=[("direct", 20, 20), ("cache", 20, 20), ("direct", 10, 20), ("cache", 10, 20), VARIANTS...])) for x in restart_seeds(N_RESTARTS)]
r4v_seeds() = [merge(x, (; n_orbits=1, configs=VARIANTS)) for x in restart_seeds(N_RESTARTS)]
r1v_seeds() = [mk("dens_cache20n$(np)_mean", c, d, "mean", 30.0; extra=(; sampler=true, npoints=np)) for (c, d, np) in VARIANTS]
timedv_seeds() = [mk("warm", "direct", 20, "pert", 1.0; stag=ek(0)); [mk("timedv_$(np)_cache20n$(np)_pert", "cache", 20, "pert", 30.0; stag=ek(0), extra=(; npoints=np)) for (_, _, np) in VARIANTS]]
const WF = [joinpath(@__DIR__, "probe_gc2_common.jl")]
campaign(task, seeds, w; days=30.0) = run_gram_campaign(task, seeds, w; mission_time_s=days * 86400.0, worker_files=WF)
workers = parse(Int, get(CLI, "workers", "3"))

function report(res, name)
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1))")
    foreach(f -> println("FAILED $(f.seed): $(first(sprint(showerror, f.error), 600))"), res.failed)
    res.route === :process || error("route was $(res.route), not :process")
    open(joinpath(PG, "mc_$(name)_$(HOST).toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nhost = \"$HOST\"\nuptime = \"$(strip(read(`uptime`, String)))\"")
        println(io, "cpu_model = \"$(try strip(split(first(filter(l -> startswith(l, "model name"), readlines("/proc/cpuinfo"))), ":"; limit=2)[2]) catch; "unknown" end)\"\nlogical_cpus = $(Sys.CPU_THREADS)")
    end
end
function flat(res, name)
    rows = [merge(s.value, (; host=HOST)) for s in res.successful]; mkpath(PG); write_rows_m(joinpath(PG, "summary_$(name)_$(HOST).csv"), rows)
    for r in rows; println(rpad(r.run_id, 30), " ", r.status, " passes=", r.passes, " wall_s=", round(r.run_wall_s; digits=1), " ", r.error); end
    return rows
end
mkpath(PG); phase = PHASE_ARG
if phase == "smoke"
    r = sequential_single_worker(gc2_states_task, [states_seed(3.0)]; worker_files=WF); report(r, "smoke_states"); flat(r, "smoke_states")
    r = sequential_single_worker(gc2_task, r1_seeds(2.0)[1:1]; worker_files=WF); report(r, "smoke_r1"); rr = flat(r, "smoke_r1")
    r = sequential_single_worker(gc2_restart_task, [merge(restart_seeds(1)[1], (; n_seeds=1))]; worker_files=WF); report(r, "smoke_r4")
    rows = reduce(vcat, [collect(s.value.rows) for s in r.successful]); write_rows_m(joinpath(PG, "restarts_smoke_$(HOST).csv"), rows)
    foreach(x -> println(x), rows); all(x -> x.status == "ok", rows) || error("smoke restart rows not ok")
    all(x -> x.status == "ok" && x.n_samples > 0, rr) || error("smoke R1 sampler produced no samples: $(rr)")
elseif phase == "states"
    r = sequential_single_worker(gc2_states_task, [states_seed(parse(Float64, get(CLI, "horizon-days", "12")))]; worker_files=WF); report(r, "states"); flat(r, "states")
elseif phase == "r1"
    r = campaign(gc2_task, r1_seeds(30.0), 2); report(r, "r1"); flat(r, "r1")
elseif phase == "r4"
    r = campaign(gc2_restart_task, restart_seeds(N_RESTARTS), workers; days=2.0); report(r, "r4")
    rows = reduce(vcat, [collect(s.value.rows) for s in r.successful]); write_rows_m(joinpath(PG, "restarts_$(HOST).csv"), rows); println("restart rows ", length(rows), " not ok ", count(x -> x.status != "ok", rows))
elseif phase == "r1v"
    r = campaign(gc2_task, r1v_seeds(), 2); report(r, "r1v"); flat(r, "r1v")
elseif phase == "r4b" || phase == "r4v"
    r = campaign(gc2_restart_task, phase == "r4b" ? r4b_seeds() : r4v_seeds(), workers; days=3.0); report(r, phase)
    rows = reduce(vcat, [collect(s.value.rows) for s in r.successful]); write_rows_m(joinpath(PG, "restarts_$(phase)_$(HOST).csv"), rows); println("restart rows ", length(rows), " not ok ", count(x -> x.status != "ok", rows))
elseif phase == "timedv"
    r = sequential_single_worker(gc2_task, timedv_seeds(); worker_files=WF); report(r, "timedv"); flat(r, "timedv")
elseif phase == "ens"
    r = campaign(gc2_task, ens_seeds(), workers); report(r, "ens"); flat(r, "ens")
elseif phase == "timed"
    r = sequential_single_worker(gc2_task, timed_seeds(); worker_files=WF); report(r, "timed"); flat(r, "timed")
else
    error("unknown phase $phase")
end
