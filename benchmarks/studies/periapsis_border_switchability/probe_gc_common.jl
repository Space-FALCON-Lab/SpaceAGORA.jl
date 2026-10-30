# Probe GC: one open-loop Mars flight with a chosen gravity degree and a direct or look-ahead-cache GRAM sampling path (plan: tasks/2026-10-15_periapsis-region-border-targeting/probe_gc_plan.md).
# New file only. Reuses the Stage I/F case helpers unchanged. The corridor guard `_corridor_gram_require_direct_sampling` is not edited: it runs at model construction inside `stagef_gram_model`, so this task
# clears the sampling variables, builds the model through the unchanged helper, and only afterwards turns the vacuum cache on in its own process (restored in `finally`). run_simulation snapshots the env at run build.
include(joinpath(@__DIR__, "stagei_common.jl"))

const PROBE_GC_ENV = ("SPACEAGORA_GRAM_TRACK_CACHE", "SPACEAGORA_GRAM_SEGMENT_CACHE", "SPACEAGORA_VACUUM_GRAM_CACHE", "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS",
    "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S", "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M", "SPACEAGORA_DENSITY_FREEZE_PER_STEP")
const PROBE_GC_CACHE = ("SPACEAGORA_VACUUM_GRAM_CACHE" => "1", "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "20", "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => "600", "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "5000")

_probe_gc_load() = try parse(Float64, split(read("/proc/loadavg", String))[1]) catch; NaN end   # 1-minute load average (Linux)

"""seed fields: run_id, cfg ("direct"|"cache"), degree, arm ("mean"|"pert"), seed_tag, tag, inc, raan, omega, epoch_idx, horizon_days, out, cal_root, config_path, profile::Bool, warm::Bool, rep.
Plain-data row; never throws."""
function probe_gc_task(seed)
    t0 = time(); saved = Dict(k => get(ENV, k, nothing) for k in (PROBE_GC_ENV..., "SPACEAGORA_GRAM_PROFILE"))
    try
        for k in PROBE_GC_ENV; delete!(ENV, k); end
        meta = stagef_meta(seed.config_path); sc = load_stagef("mars"; path=seed.config_path)
        sc["gravity"]["degree"] = seed.degree; sc["gravity"]["order"] = seed.degree; ms0 = stagef_setup(sc)
        epoch = parse_epoch(sc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms0.planet, epoch)
        cdir = joinpath(seed.cal_root, seed.tag); G = CSV.read(joinpath(cdir, "guess.csv"), DataFrame)[1, :]; N = CSV.read(joinpath(cdir, "nodes.csv"), DataFrame)
        pwf = pw_fit(N.h_m, N.rho_kg_m3, Float64(sc["atmosphere"]["temperature_k"]), "piecewise"); ms = with_fit(ms0, pwf)
        pert = seed.arm == "pert"; stagef_env!(pert, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed("mars" * seed.seed_tag))   # guard runs here with every sampling variable cleared
        core_t = model.core.initial_time
        seed.cfg == "cache" && for (k, v) in PROBE_GC_CACHE; ENV[k] = v; end   # bypass: cache on strictly after the guarded construction
        seed.profile && (ENV["SPACEAGORA_GRAM_PROFILE"] = "1")
        dir = joinpath(seed.out, seed.run_id); mkpath(dir); horizon = Float64(seed.horizon_days) * SECONDS_PER_DAY
        state0 = mars_initial_state(ms, rot, Float64(G.seed_alt_m), seed.inc, seed.raan, seed.omega)
        args = gram_args(ms, sc, model, epoch, horizon, dir; state=state0)
        c0 = seed.profile ? _gram_density_evaluations() : 0
        load0 = _probe_gc_load(); tstart = time(); tr = time(); sol = run_simulation(args; return_solution=true, isolate_state=false); run_wall = time() - tr; load1 = _probe_gc_load()
        calls = seed.profile ? _gram_density_evaluations() - c0 : -1
        logpath = joinpath(dir, "gram_density_perturbation_log.csv")
        rows = stagef_table(sol, args, ms, rot, sc, pwf, logpath)   # the perturbation log and its .summary.toml are kept (stagef_run_task deletes them)
        write_rows_m(joinpath(dir, "passes.csv"), rows)
        nf = try Int(sol.stats.nf) catch; -1 end; nacc = try Int(sol.stats.naccept) catch; -1 end; nrej = try Int(sol.stats.nreject) catch; -1 end
        ret = string(sol.retcode)
        return (; run_id=seed.run_id, cfg=seed.cfg, degree=seed.degree, arm=seed.arm, seed_tag=seed.seed_tag, rep=seed.rep, warm=seed.warm, profile=seed.profile, status="ok", retcode=ret,
            passes=length(rows), t_end_days=sol.t[end] / SECONDS_PER_DAY, horizon_days=seed.horizon_days, n_sol_t=length(sol.t), nf=nf, naccept=nacc, nreject=nrej, density_calls=calls,
            gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second),
            t_start_unix=tstart, load1_before=load0, load1_after=load1, run_wall_s=run_wall, wall_s=time() - t0, error="")
    catch err
        return (; run_id=seed.run_id, cfg=seed.cfg, degree=seed.degree, arm=seed.arm, seed_tag=seed.seed_tag, rep=seed.rep, warm=seed.warm, profile=seed.profile, status="run_error", retcode="",
            passes=0, t_end_days=NaN, horizon_days=seed.horizon_days, n_sol_t=0, nf=-1, naccept=-1, nreject=-1, density_calls=-1, gram_core_epoch="", t_start_unix=NaN, load1_before=NaN, load1_after=NaN, run_wall_s=NaN,
            wall_s=time() - t0, error=replace(first(sprint(showerror, err), 500), "\n" => " "))
    finally
        for (k, v) in saved; v === nothing ? delete!(ENV, k) : (ENV[k] = v); end
    end
end
