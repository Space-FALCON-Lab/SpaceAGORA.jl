# Probe GC-2 flight machinery (plan: tasks/2026-10-15_periapsis-region-border-targeting/probe_gc2_plan.md). New file only; reuses probe_gc_common.jl (env handling, guard bypass) and the Stage F/I helpers unchanged.
# Guard bypass: as probe_gc_common.jl: the sampling variables are cleared, the GRAM model is built through the unchanged `stagef_gram_model`, then the cache variables are set in this process and restored in `finally`.
using OrdinaryDiffEq
include(joinpath(@__DIR__, "probe_gc_common.jl"))

"""Keep the perturbation knots (kind 3) as pass,t,alt,r (every knot at or below 300 km, every 5th above) in knots.csv (the chain gzips it) and delete the full log (kept: the summary toml)."""
function gc2_compact_log!(dir::AbstractString)
    log = joinpath(dir, "gram_density_perturbation_log.csv"); isfile(log) || return nothing
    open(joinpath(dir, "knots.csv"), "w") do io
        println(io, "pass,t_s,alt_m,r"); n = 0
        for line in eachline(log)
            startswith(line, "3,") || continue
            f = split(line, ','); alt = parse(Float64, f[5]); n += 1
            (alt <= 300e3 || n % 5 == 0) && println(io, f[3], ",", f[4], ",", f[5], ",", f[6])
        end
    end
    rm(log; force=true); return nothing
end

"""Discrete callback: every `dt` s in the atmosphere, the engine's live look-ahead spline value vs a direct GRAM mean sample at the identical state (read-only; never rebuilds the cache)."""
function gc2_sampler(rows::Vector{NamedTuple}, args, model, dt::Float64, dev_m::Float64)
    S = SimulationModel.SimulationCallbacks; planet = args.environment_model.planet; pargs = (args=args,)
    nxt = Ref(0.0); lastt = Ref(-1.0e9); pid = Ref(0); miss = Ref(0)
    cond(u, t, i) = t >= nxt[]
    function aff(i)
        t = Float64(i.t); nxt[] = t + dt; p = i.p; sb = p.shared_buffers
        (length(sb.in_atmosphere) >= 1 && sb.in_atmosphere[1]) || return nothing
        c = sb.vacuum_gram_caches[1]; c === nothing && return nothing
        sc1 = i.u.sc[1]; pos = SVector{3, Float64}(sc1.pos); vel = SVector{3, Float64}(sc1.vel)
        t - lastt[] > 600.0 && (pid[] += 1); lastt[] = t
        if !(c.valid && t >= c.t0 && t <= c.t1 && norm(pos - S._interp_vacuum_position(c, t)) <= dev_m)
            miss[] += 1; return nothing
        end
        rho_c = exp(S._eval_natural_cubic_spline(t, c.t0, c.h, c.log_rhos, c.Ms_rho))
        l_pi = S._planet_lpi_at(p, t); pos_pp, vel_pp = S._planet_relative_state(pos, vel, planet, l_pi)
        alt, lat, lon = S.rtolatlong(pos_pp, planet)
        rho_d = gram_rho(model, pargs, alt, lat, lon, t)
        push!(rows, (; pass=pid[], t_s=t, alt_m=Float64(alt), rho_cache=rho_c, rho_direct=rho_d, vrel_mps=Float64(norm(vel_pp))))
        return nothing
    end
    return DiscreteCallback(cond, aff; save_positions=(false, false)), miss
end

"""seed: run_id, cfg, degree, arm, seed_tag, tag, inc, raan, omega, epoch_idx, out, cal_root, config_path, plus horizon_days (or horizon_s), optional restart=(r, v, t_s), sampler::Bool, keep_log::Bool.
Returns (row, table, sol, args, wall); the caller turns it into plain data."""
function gc2_fly(seed)
    saved = Dict(k => get(ENV, k, nothing) for k in (PROBE_GC_ENV..., "SPACEAGORA_GRAM_PROFILE"))
    try
        for k in PROBE_GC_ENV; delete!(ENV, k); end
        meta = stagef_meta(seed.config_path); sc = load_stagef("mars"; path=seed.config_path)
        sc["gravity"]["degree"] = seed.degree; sc["gravity"]["order"] = seed.degree; ms0 = stagef_setup(sc)
        epoch0 = parse_epoch(sc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms0.planet, epoch0)
        cdir = joinpath(seed.cal_root, seed.tag); G = CSV.read(joinpath(cdir, "guess.csv"), DataFrame)[1, :]; N = CSV.read(joinpath(cdir, "nodes.csv"), DataFrame)
        pwf = pw_fit(N.h_m, N.rho_kg_m3, Float64(sc["atmosphere"]["temperature_k"]), "piecewise"); ms = with_fit(ms0, pwf)
        pert = seed.arm == "pert"; stagef_env!(pert, meta)
        restart = haskey(seed, :restart) ? seed.restart : nothing
        epoch = restart === nothing ? epoch0 : _shift_initial_time(epoch0, Float64(restart.t_s))
        model = stagef_gram_model(sc, epoch, stagef_seed("mars" * seed.seed_tag))   # guard runs here with every sampling variable cleared
        seed.cfg == "cache" && for (k, v) in PROBE_GC_CACHE; ENV[k] = v; end
        seed.cfg == "cache" && haskey(seed, :npoints) && (ENV["SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS"] = string(seed.npoints))   # exploratory variants only
        dir = joinpath(seed.out, seed.run_id); mkpath(dir)
        horizon = haskey(seed, :horizon_s) ? Float64(seed.horizon_s) : Float64(seed.horizon_days) * SECONDS_PER_DAY
        state0 = restart === nothing ? mars_initial_state(ms, rot, Float64(G.seed_alt_m), seed.inc, seed.raan, seed.omega) : (SVector{3, Float64}(restart.r), SVector{3, Float64}(restart.v))
        args = gram_args(ms, sc, model, epoch, horizon, dir; state=state0)
        srows = NamedTuple[]; extra = (); miss = Ref(0)
        if get(seed, :sampler, false)
            cb, miss = gc2_sampler(srows, args, model, 2.0, 5000.0); extra = (cb,)
        end
        load0 = _probe_gc_load(); tstart = time(); sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=extra); wall = time() - tstart; load1 = _probe_gc_load()
        logpath = joinpath(dir, "gram_density_perturbation_log.csv")
        rows = stagef_table(sol, args, ms, rot, sc, pwf, logpath)
        if get(seed, :keep_log, false)
            gc2_compact_log!(dir)
        else
            rm(logpath; force=true); rm(logpath * ".summary.toml"; force=true)
        end
        isempty(srows) || write_rows_m(joinpath(dir, "r1_samples.csv"), srows)
        get(seed, :write_passes, true) && write_rows_m(joinpath(dir, "passes.csv"), rows)
        nf = try Int(sol.stats.nf) catch; -1 end
        row = (; run_id=seed.run_id, cfg=seed.cfg, degree=seed.degree, arm=seed.arm, seed_tag=seed.seed_tag, status="ok", retcode=string(sol.retcode), passes=length(rows),
            t_end_days=sol.t[end] / SECONDS_PER_DAY, n_sol_t=length(sol.t), nf=nf, n_samples=length(srows), n_miss=miss[], load1_before=load0, load1_after=load1, run_wall_s=wall, error="")
        return (; row, table=rows, sol, args)
    finally
        for (k, v) in saved; v === nothing ? delete!(ENV, k) : (ENV[k] = v); end
    end
end

gc2_err(seed, err) = (; run_id=seed.run_id, cfg=seed.cfg, degree=seed.degree, arm=seed.arm, seed_tag=seed.seed_tag, status="run_error", retcode="", passes=0, t_end_days=NaN, n_sol_t=0, nf=-1,
    n_samples=0, n_miss=0, load1_before=NaN, load1_after=NaN, run_wall_s=NaN, error=replace(first(sprint(showerror, err), 500), "\n" => " "))

"""Ensemble / R1 / timed flight: plain-data row."""
function gc2_task(seed)
    try; return gc2_fly(seed).row; catch err; return gc2_err(seed, err); end
end

"""Direct 20x20 mean-arm trajectory: apsis events -> states.csv (apoapsis state per row with the time of the next apoapsis and the next periapsis)."""
function gc2_states_task(seed)
    try
        f = gc2_fly(seed); ev = physical_apsis_events(f.sol, _geodetic_of(f.args)); apo = [e for e in ev if e.kind == "apoapsis"]; per = [e for e in ev if e.kind == "periapsis"]
        rows = NamedTuple[]
        for i in 1:(length(apo) - 1)
            a = apo[i]; nx = apo[i + 1]; pp = [e for e in per if e.time_s > a.time_s && e.time_s < nx.time_s]
            isempty(pp) && continue
            push!(rows, (; i, t_s=a.time_s, rx=a.position_ii_m[1], ry=a.position_ii_m[2], rz=a.position_ii_m[3], vx=a.velocity_ii_mps[1], vy=a.velocity_ii_mps[2], vz=a.velocity_ii_mps[3],
                mass_kg=a.mass_kg, apo_alt_m=a.physical_altitude_m, t_next_apo_s=nx.time_s, next_apo_alt_m=nx.physical_altitude_m, next_peri_alt_m=first(pp).physical_altitude_m))
        end
        write_rows_m(joinpath(seed.out, "states.csv"), rows)
        return merge(f.row, (; n_states=length(rows)))
    catch err; return merge(gc2_err(seed, err), (; n_states=0)); end
end

const GC2_CONFIGS = (("direct", 20), ("cache", 20), ("direct", 10), ("cache", 10))

"""One restart: configurations x (mean arm + n_seeds seeds), n_orbits orbits ahead. seed: restart_i, states_path, out, cal_root, config_path, tag, inc, raan, omega, epoch_idx, n_seeds,
optional n_orbits (default 1) and configs (vector of (cfg, degree, npoints), default the four plan configurations at 20 points). Returns (; rows)."""
function gc2_restart_task(seed)
    rows = NamedTuple[]; no = get(seed, :n_orbits, 1); cfgs = get(seed, :configs, [(c, d, 20) for (c, d) in GC2_CONFIGS])
    row(; kw...) = merge((; restart=seed.restart_i, n_orbits=no), (; kw...))
    try
        S = CSV.read(seed.states_path, DataFrame); s = S[S.i .== seed.restart_i, :][1, :]; sl = S[S.i .== seed.restart_i + no - 1, :][1, :]
        rs = (r=SVector(s.rx, s.ry, s.rz), v=SVector(s.vx, s.vy, s.vz), t_s=s.t_s); horizon = (sl.t_next_apo_s - s.t_s) + 600.0
        for (arm, nseed) in (("mean", 1), ("pert", seed.n_seeds)), j in 1:nseed, (cfg, deg, np) in cfgs
            id = "r$(seed.restart_i)o$(no)_$(cfg)$(deg)n$(np)_$(arm)_s$(j)"
            sd = (; config_path=seed.config_path, run_id=id, cfg, degree=deg, arm, seed_tag=seed.tag * "|r$(seed.restart_i)|s$(j)", tag=seed.tag, inc=seed.inc, raan=seed.raan, omega=seed.omega,
                epoch_idx=seed.epoch_idx, horizon_s=horizon, out=seed.out, cal_root=seed.cal_root, restart=rs, write_passes=false, npoints=np)
            try
                f = gc2_fly(sd); t = f.table; rm(joinpath(seed.out, id); recursive=true, force=true)
                if length(t) < no; push!(rows, row(; cfg, degree=deg, np, arm, seed=j, status="no_periapsis", q_eff=NaN, alt_m=NaN, apo_change_m=NaN, wall_s=f.row.run_wall_s)); continue; end
                pn = t[no]; push!(rows, row(; cfg, degree=deg, np, arm, seed=j, status="ok", q_eff=pn.q_eff_w_cm2, alt_m=pn.alt_m, apo_change_m=pn.apo_alt_m - s.apo_alt_m, wall_s=f.row.run_wall_s))
            catch err
                push!(rows, row(; cfg, degree=deg, np, arm, seed=j, status="error:" * replace(first(sprint(showerror, err), 200), "," => ";"), q_eff=NaN, alt_m=NaN, apo_change_m=NaN, wall_s=NaN))
            end
        end
    catch err
        push!(rows, row(; cfg="", degree=0, np=0, arm="", seed=0, status="error:" * replace(first(sprint(showerror, err), 200), "," => ";"), q_eff=NaN, alt_m=NaN, apo_change_m=NaN, wall_s=NaN))
    end
    return (; rows)
end
