# Part C: switchability by predictor finite differences (impulsive kicks, full truth model).
#   --phase=states   extract the unkicked states at the three kick locations (campaign start and mid-campaign)
#   --phase=kicks    reference + kicked runs per (geometry, start, location); writes switch_fd.csv
#   --phase=steps    direct periapsis step (tangential kick at apoapsis, secant on the step) and hold behavior
# Geometries come from <out>/geometries.csv (written by select_geometries.py from the Part B grid).
include(joinpath(@__DIR__, "common.jl"))

const SIGNS = (1.0, -1.0)

function write_rows(path, rows)
    mkpath(dirname(path))
    CSV.write(path, reduce((x, y) -> vcat(x, y; cols=:union), (DataFrame([r]) for r in rows)))  # rows may differ in fields (reference vs kick)
end

function geometry_config(study, g, out)
    cfg = make_study_config(study; raan_deg=g.raan, omega_deg=g.omega)
    cal = CSV.read(joinpath(out, "grid", g.tag, "calibration.csv"), DataFrame)
    cfg.values["orbit"]["initial_periapsis_altitude_m"] = cal.seed_alt_m[end]  # the accepted calibration seed
    return cfg, parse_epoch(study["grid"]["epochs_utc"][g.epoch_idx])
end

state_at(sol, t) = (r=SVector{3, Float64}(sol(t).sc[1].pos), v=SVector{3, Float64}(sol(t).sc[1].vel), m=Float64(sol(t).sc[1].mass))

"""States of the unkicked 3B run at apoapsis, approach and nu=135 deg of the leg ending at the first periapsis at or after `target_day`."""
function extract_states(g, start::String, study, out)
    cfg, epoch = geometry_config(study, g, out)
    sw = study["switch"]
    target_s = start == "campaign" ? 0.0 : Float64(sw["mid_campaign_day"]) * SECONDS_PER_DAY
    horizon = target_s + 2.2 * 2π * sqrt(7.0e7^3 / 3.986004418e14)  # ~2 initial periods past the target
    res = run_truth(cfg, study, epoch, horizon, joinpath(out, "states", g.tag); third_bodies=true)
    geo = _geodetic_of(res.args); mu = res.args.environment_model.planet.μ
    events = physical_apsis_events(res.sol, geo)
    peri = [e for e in events if e.kind == "periapsis"]; apo = [e for e in events if e.kind == "apoapsis"]
    kp = findfirst(e -> e.time_s >= target_s - 1.0, peri)
    kp === nothing && throw(ErrorException("no periapsis after day $(target_s / SECONDS_PER_DAY)"))
    tp = peri[kp].time_s
    prev_apo = [e for e in apo if e.time_s < tp - 1.0]
    t_apo = isempty(prev_apo) ? 0.0 : prev_apo[end].time_s
    next_apo = [e for e in apo if e.time_s > tp + 1.0]
    isempty(next_apo) && throw(ErrorException("run too short: no apoapsis after periapsis $kp"))
    t_apo_next = next_apo[1].time_s
    t_appr = find_zero(t -> geo(t, state_at(res.sol, t).r)[1] - Float64(sw["approach_altitude_m"]), (t_apo + 1.0, tp - 1.0), Roots.Bisection(); xatol=1e-3)
    nu_target = deg2rad(Float64(sw["true_anomaly_deg"]))
    t_nu = find_zero(t -> (s = state_at(res.sol, t); _corridor_element_summary(s.r, s.v, mu).nu - nu_target),
        (tp + 1.0, t_apo_next - 1.0), Roots.Bisection(); xatol=1e-3)
    rows = NamedTuple[]
    for (loc, t) in (("apoapsis", t_apo), ("approach", t_appr), ("nu135", t_nu))
        s = state_at(res.sol, t)
        push!(rows, (; tag=g.tag, start, location=loc, periapsis_index=kp, t_s=t, mass_kg=s.m,
            rx=s.r[1], ry=s.r[2], rz=s.r[3], vx=s.v[1], vy=s.v[2], vz=s.v[3]))
    end
    return rows, res.wall_s
end

function states_task(seed)
    study = load_study(seed.study_path); nbody_planet()
    rows, wall = extract_states(seed.g, seed.start, study, seed.out)
    return (; rows, wall_s=wall)
end

direction_vector(kind, r, v) = kind == "tangential" ? normalize(v) : kind == "radial" ? normalize(r) : normalize(cross(r, v))

"""Least-squares line through (x, y): returns (intercept, slope)."""
function linfit(x, y)
    n = length(x); mx, my = sum(x) / n, sum(y) / n
    sxx = sum((x .- mx) .^ 2)
    slope = sxx > 0 ? sum((x .- mx) .* (y .- my)) / sxx : NaN
    return (my - slope * mx, slope)
end

angdiff(a, b) = mod(a - b + π, 2π) - π

"""Run the unkicked restart (reference) and every kick for one (geometry, start, location). Returns plain rows."""
function kicks_task(seed)
    study = load_study(seed.study_path); nbody_planet()
    sw = study["switch"]; g = seed.g
    cfg, epoch0 = geometry_config(study, g, seed.out)
    st = seed.state
    r0 = SVector(st.rx, st.ry, st.rz); v0 = SVector(st.vx, st.vy, st.vz)
    epoch = _shift_initial_time(epoch0, Float64(st.t_s))
    horizon = Float64(sw["propagate_days"]) * SECONDS_PER_DAY
    dir = joinpath(seed.out, "kicks", g.tag, seed.start, seed.location)
    ref = run_truth(cfg, study, epoch, horizon, dir; third_bodies=true, state=(r0, v0), mass_kg=st.mass_kg)
    mu = ref.args.environment_model.planet.μ
    fit = fitted_atmosphere_parameters(cfg)
    with_band(rows) = [merge(r, (h_qhigh_m=band_altitudes(fit, r.vrel_mps, study)[1], h_qlow_m=band_altitudes(fit, r.vrel_mps, study)[2])) for r in rows]
    write_rows(joinpath(dir, "passes_ref.csv"), with_band(ref.rows))
    tref = [r.time_s for r in ref.rows] ./ SECONDS_PER_DAY
    aref = [r.alt_m for r in ref.rows] ./ 1000
    ref_icpt, ref_slope = linfit(tref, aref)
    el0 = _corridor_element_summary(r0, v0, mu)
    out = NamedTuple[(; tag=g.tag, role=g.role, start=seed.start, location=seed.location, direction="reference", dv_mps=0.0,
        ref_slope_km_day=ref_slope, ref_alt0_km=aref[1], ref_q0=ref.rows[1].q_w_cm2, ref_passes=length(ref.rows),
        wall_s=ref.wall_s, retcode=ref.retcode)]
    dirs = haskey(seed, :dirs) ? seed.dirs : String.(sw["directions"])
    dvs = haskey(seed, :dvs) ? seed.dvs : Float64.(sw["probe_dv_mps"])
    for kind in dirs, dv in dvs, sgn in SIGNS
        v1 = v0 + sgn * dv * direction_vector(kind, r0, v0)
        el1 = _corridor_element_summary(r0, v1, mu)
        run = run_truth(cfg, study, epoch, horizon, dir; third_bodies=true, state=(r0, v1), mass_kg=st.mass_kg)
        n = min(length(run.rows), length(ref.rows))
        d = [(run.rows[j].alt_m - ref.rows[j].alt_m) / 1000 for j in 1:n]
        icpt, slope = linfit(tref[1:n], d)
        push!(out, (; tag=g.tag, role=g.role, start=seed.start, location=seed.location, direction=kind, dv_mps=sgn * dv,
            ref_slope_km_day=ref_slope, ref_alt0_km=aref[1], ref_q0=ref.rows[1].q_w_cm2, ref_passes=length(ref.rows),
            passes=n, d_first_km=d[1], d_last_km=d[n], d_intercept_km=icpt, delta_rate_km_day=slope,
            d_apo_first_km=(run.rows[1].apo_alt_m - ref.rows[1].apo_alt_m) / 1000,
            dq_first=run.rows[1].q_w_cm2 - ref.rows[1].q_w_cm2,
            kick_da_m=el1.a - el0.a, kick_de=el1.e - el0.e, kick_domega_rad=angdiff(el1.omega, el0.omega),
            kick_di_rad=el1.i - el0.i, kick_draan_rad=angdiff(el1.raan, el0.raan),
            end_da_m=run.rows[n].a_m - ref.rows[n].a_m, end_de=run.rows[n].e - ref.rows[n].e,
            end_domega_rad=angdiff(run.rows[n].omega_rad, ref.rows[n].omega_rad), end_di_rad=run.rows[n].i_rad - ref.rows[n].i_rad,
            wall_s=run.wall_s, retcode=run.retcode))
    end
    return out
end

"""Direct-step cost: tangential kick at apoapsis sized (secant) to move the next-pass periapsis by +/- step_km."""
function steps_task(seed)
    study = load_study(seed.study_path); nbody_planet()
    sw = study["switch"]; g = seed.g
    cfg, epoch0 = geometry_config(study, g, seed.out)
    st = seed.state
    r0 = SVector(st.rx, st.ry, st.rz); v0 = SVector(st.vx, st.vy, st.vz)
    epoch = _shift_initial_time(epoch0, Float64(st.t_s))
    horizon = Float64(seed.hold_days) * SECONDS_PER_DAY
    dir = joinpath(seed.out, "steps", g.tag, seed.start)
    fit = fitted_atmosphere_parameters(cfg)
    ref = run_truth(cfg, study, epoch, horizon, dir; third_bodies=true, state=(r0, v0), mass_kg=st.mass_kg)
    tref = [r.time_s for r in ref.rows] ./ SECONDS_PER_DAY; aref = [r.alt_m for r in ref.rows] ./ 1000
    out = NamedTuple[]
    for sgn in SIGNS
        target = sgn * Float64(sw["direct_step_km"])
        dv = target / seed.km_per_mps       # linear first guess from the finite-difference table
        local run, n, d, dv_applied
        for it in 1:4  # secant on the step: rescale dv by target / achieved
            dv_applied = dv
            run = run_truth(cfg, study, epoch, horizon, dir; third_bodies=true, state=(r0, v0 + dv * normalize(v0)), mass_kg=st.mass_kg)
            n = min(length(run.rows), length(ref.rows))
            d = [(run.rows[j].alt_m - ref.rows[j].alt_m) / 1000 for j in 1:n]
            abs(d[1] - target) < 0.05 && break
            dv *= target / d[1]
        end
        stats = band_stats(run.rows, study)
        push!(out, (; tag=g.tag, role=g.role, start=seed.start, target_step_km=target, dv_applied_mps=dv_applied, achieved_first_km=d[1],
            cost_mps_per_km=abs(dv_applied) / abs(d[1]), retained_day5_km=interp_series(tref[1:n], d, 5.0),
            retained_day10_km=interp_series(tref[1:n], d, 10.0), retained_last_km=d[n], last_day=tref[n],
            ref_drift_km_day=linfit(tref, aref)[2], stepped_drift_km_day=linfit([r.time_s for r in run.rows] ./ SECONDS_PER_DAY, [r.alt_m for r in run.rows] ./ 1000)[2],
            stepped_q_first=run.rows[1].q_w_cm2, ref_q_first=ref.rows[1].q_w_cm2,
            stepped_first_exit_days=stats.first_exit_time_days, stepped_first_exit_edge=stats.first_exit_edge,
            ref_first_exit_days=band_stats(ref.rows, study).first_exit_time_days, wall_s=run.wall_s))
        write_rows(joinpath(dir, "passes_step_$(sgn > 0 ? "up" : "down").csv"), run.rows)
    end
    write_rows(joinpath(dir, "passes_ref.csv"), ref.rows)
    return out
end

function interp_series(t, y, t0)
    (t0 < t[1] || t0 > t[end]) && return NaN
    j = searchsortedlast(t, t0); j >= length(t) && return y[end]
    w = (t0 - t[j]) / (t[j + 1] - t[j]); return (1 - w) * y[j] + w * y[j + 1]
end

function parse_args(args)
    o = Dict{String, Any}("out" => nothing, "phase" => "kicks", "workers" => 2, "hold_days" => 30.0, "only" => nothing, "dv" => nothing)
    for a in args
        startswith(a, "--out=") && (o["out"] = abspath(split(a, "=", limit=2)[2]))
        startswith(a, "--phase=") && (o["phase"] = split(a, "=", limit=2)[2])
        startswith(a, "--workers=") && (o["workers"] = parse(Int, split(a, "=", limit=2)[2]))
        startswith(a, "--hold-days=") && (o["hold_days"] = parse(Float64, split(a, "=", limit=2)[2]))
        startswith(a, "--only=") && (o["only"] = split(split(a, "=", limit=2)[2], ","))  # geometry tags
    end
    return o
end

function geometries(out, only)
    df = CSV.read(joinpath(out, "geometries.csv"), DataFrame)
    only === nothing || (df = df[in.(df.tag, Ref(only)), :])
    return [(; tag=r.tag, raan=r.raan_deg, omega=r.omega_deg, epoch_idx=r.epoch_idx, role=r.role) for r in eachrow(df)]
end

function switch_main(args=ARGS)
    o = parse_args(args); study = load_study(); out = something(o["out"], results_root(study))
    gs = geometries(out, o["only"]); starts = ("campaign", "mid")
    files = [joinpath(@__DIR__, "run_switch.jl")]
    mt = Float64(study["switch"]["propagate_days"]) * SECONDS_PER_DAY
    if o["phase"] == "states"
        seeds = [(; study_path=STUDY_CONFIG_PATH, out, g, start=s) for g in gs for s in starts]
        res = run_process_campaign(states_task, seeds, o["workers"]; mission_time_s=mt, worker_files=files)
        report(res)
        write_rows(joinpath(out, "kick_states.csv"), [r for s in res.successful for r in s.value.rows])
    elseif o["phase"] == "kicks"
        states = CSV.read(joinpath(out, "kick_states.csv"), DataFrame)
        seeds = [(; study_path=STUDY_CONFIG_PATH, out, g, start=r.start, location=r.location, state=NamedTuple(r))
                 for g in gs for r in eachrow(states) if r.tag == g.tag]
        res = run_process_campaign(kicks_task, seeds, o["workers"]; mission_time_s=mt, worker_files=files)
        report(res)
        write_rows(joinpath(out, "switch_fd.csv"), [r for s in res.successful for r in s.value])
    elseif o["phase"] == "verify"  # the best lever's kick at the extrapolated switch_dv values (checks the linear extrapolation)
        states = CSV.read(joinpath(out, "kick_states.csv"), DataFrame)
        sw_df = CSV.read(joinpath(out, "switchability.csv"), DataFrame)
        seeds = NamedTuple[]
        for g in gs, s in starts
            sub = sw_df[(sw_df.tag .== g.tag) .& (sw_df.start .== s), :]
            loc, dirn = sub.best_location[1], sub.best_direction[1]   # the step-compensated best lever
            row = states[(states.tag .== g.tag) .& (states.start .== s) .& (states.location .== loc), :]
            push!(seeds, (; study_path=STUDY_CONFIG_PATH, out, g, start=s, location=loc, state=NamedTuple(row[1, :]),
                dirs=[dirn], dvs=unique(round.(sub.kick_dv_mps; digits=3))))
        end
        res = run_process_campaign(kicks_task, seeds, o["workers"]; mission_time_s=mt, worker_files=files)
        report(res)
        write_rows(joinpath(out, "switch_verify.csv"), [r for s in res.successful for r in s.value])
    elseif o["phase"] == "steps"
        states = CSV.read(joinpath(out, "kick_states.csv"), DataFrame)
        fd = CSV.read(joinpath(out, "switch_fd.csv"), DataFrame)
        seeds = NamedTuple[]
        for g in gs, s in starts
            row = states[(states.tag .== g.tag) .& (states.start .== s) .& (states.location .== "apoapsis"), :]
            sel = fd[(fd.tag .== g.tag) .& (fd.start .== s) .& (fd.location .== "apoapsis") .& (fd.direction .== "tangential") .& (fd.dv_mps .== Float64(study["switch"]["probe_dv_mps"][end])), :]
            push!(seeds, (; study_path=STUDY_CONFIG_PATH, out, g, start=s, state=NamedTuple(row[1, :]), hold_days=o["hold_days"],
                km_per_mps=sel.d_first_km[1] / sel.dv_mps[1]))
        end
        res = run_process_campaign(steps_task, seeds, o["workers"]; mission_time_s=o["hold_days"] * SECONDS_PER_DAY, worker_files=files)
        report(res)
        write_rows(joinpath(out, "direct_step.csv"), [r for s in res.successful for r in s.value])
    else
        throw(ArgumentError("--phase must be states, kicks, verify or steps"))
    end
end

function report(res)
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1))")
    for f in res.failed
        println("FAILED $(f.seed.g.tag) $(get(f.seed, :start, "")) $(get(f.seed, :location, "")): $(first(sprint(showerror, f.error), 400))")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    switch_main()
end
