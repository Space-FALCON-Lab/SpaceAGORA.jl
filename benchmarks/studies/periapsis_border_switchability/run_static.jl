# Static-region sweep. Phases (one Julia job at a time, process route via run_monte_carlo):
#   --pilot                     one cell, serial, timed
#   --phase=grid                inclination x RAAN x omega x epoch, 60 d, Sun+Moon      -> grid/<tag>/passes_3b.csv
#   --phase=control             no-third-body control for <out>/control_cells.csv       -> grid/<tag>/passes_ctl.csv
#   --phase=race                long passive runs for <out>/race_cells.csv + FD lever   -> race/<tag>/passes_race.csv
# Selection files are written by analyze_static.py between phases. common.jl is unchanged.
include(joinpath(@__DIR__, "common.jl"))

const STATIC_CONFIG_PATH = joinpath(@__DIR__, "static_config.toml")
cell_tag(inc, raan, omega, e) = "I$(lpad(round(Int, 10inc), 4, '0'))_R$(lpad(round(Int, raan), 3, '0'))_W$(lpad(round(Int, omega), 3, '0'))_E$(e)"

write_rows(path, rows) = (mkpath(dirname(path)); CSV.write(path, DataFrame(rows)))
state_at(sol, t) = (r=SVector{3, Float64}(sol(t).sc[1].pos), v=SVector{3, Float64}(sol(t).sc[1].vel), m=Float64(sol(t).sc[1].mass))

function cell_config(study, inc, raan, omega; seed_alt=nothing)
    cfg = make_study_config(study; raan_deg=raan, omega_deg=omega)
    cfg.values["orbit"]["inclination_deg"] = Float64(inc)
    seed_alt === nothing || (cfg.values["orbit"]["initial_periapsis_altitude_m"] = Float64(seed_alt))
    return cfg
end

seed_alt_of(dir) = CSV.read(joinpath(dir, "calibration.csv"), DataFrame).seed_alt_m[end]  # the accepted calibration seed

"""Add the heat-limited altitudes (altitude where peak Q = q_high / q_low at the pass speed) to each row."""
function with_band(rows, cfg, study)
    fit = fitted_atmosphere_parameters(cfg)
    return [merge(r, (h_qhigh_m=b[1], h_qlow_m=b[2])) for r in rows for b in (band_altitudes(fit, r.vrel_mps, study),)]
end

function grid_cell_task(seed)
    study = load_study(seed.study_path); nbody_planet()
    cfg0 = cell_config(study, seed.inc, seed.raan, seed.omega)
    epoch = parse_epoch(study["grid"]["epochs_utc"][seed.epoch_idx])
    tag = cell_tag(seed.inc, seed.raan, seed.omega, seed.epoch_idx)
    dir = joinpath(seed.out, "grid", tag)
    t0 = time()
    cal, recs = calibrate_peak_q(cfg0, study, epoch, dir; third_bodies=true)
    cal_wall = time() - t0
    r3 = run_truth(cal, study, epoch, Float64(seed.horizon_days) * SECONDS_PER_DAY, dir; third_bodies=true, planes=true)
    write_rows(joinpath(dir, "passes_3b.csv"), with_band(r3.rows, cal, study))
    write_rows(joinpath(dir, "calibration.csv"), recs)
    return (; tag, inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, epoch_idx=seed.epoch_idx, cal_iterations=length(recs),
        seed_alt_m=recs[end].seed_alt_m, cal_wall_s=cal_wall, wall_3b_s=r3.wall_s, retcode_3b=r3.retcode, passes_3b=length(r3.rows))
end

function control_cell_task(seed)
    study = load_study(seed.study_path); nbody_planet()
    dir = joinpath(seed.out, "grid", seed.tag)
    cfg = cell_config(study, seed.inc, seed.raan, seed.omega; seed_alt=seed_alt_of(dir))
    epoch = parse_epoch(study["grid"]["epochs_utc"][seed.epoch_idx])
    rc = run_truth(cfg, study, epoch, Float64(seed.horizon_days) * SECONDS_PER_DAY, dir; third_bodies=false, planes=true)
    write_rows(joinpath(dir, "passes_ctl.csv"), with_band(rc.rows, cfg, study))
    return (; tag=seed.tag, wall_ctl_s=rc.wall_s, retcode_ctl=rc.retcode, passes_ctl=length(rc.rows))
end

"""Tangential kick at the last apoapsis that has a following periapsis; first-pass periapsis change per kick (finite difference)."""
function end_lever(res, cfg, study, epoch, dir, dv)
    geo = _geodetic_of(res.args); ev = physical_apsis_events(res.sol, geo)
    ka = findlast(e -> e.kind == "apoapsis" && any(p -> p.kind == "periapsis" && p.time_s > e.time_s, ev), ev)
    ka === nothing && return (; lever_ok=false)
    ta = ev[ka].time_s
    tp = first(p.time_s for p in ev if p.kind == "periapsis" && p.time_s > ta)
    s = state_at(res.sol, ta)
    ep = _shift_initial_time(epoch, ta)
    horizon = tp - ta + 600.0
    ref = run_truth(cfg, study, ep, horizon, dir; third_bodies=true, state=(s.r, s.v), mass_kg=s.m)
    up = run_truth(cfg, study, ep, horizon, dir; third_bodies=true, state=(s.r, s.v + dv * normalize(s.v)), mass_kg=s.m)
    dn = run_truth(cfg, study, ep, horizon, dir; third_bodies=true, state=(s.r, s.v - dv * normalize(s.v)), mass_kg=s.m)
    (isempty(ref.rows) || isempty(up.rows) || isempty(dn.rows)) && return (; lever_ok=false)
    d = (up.rows[1].alt_m - dn.rows[1].alt_m) / 2   # central difference, m of periapsis altitude per dv
    el = _corridor_element_summary(s.r, s.v, ref.args.environment_model.planet.μ)
    return (; lever_ok=true, lever_t_apo_days=ta / SECONDS_PER_DAY, lever_a_m=el.a, lever_e=el.e, lever_dv_mps=dv,
        lever_dalt_m=d, lever_mps_per_km=dv / (abs(d) / 1000), lever_ref_alt_m=ref.rows[1].alt_m)
end

function race_cell_task(seed)
    study = load_study(seed.study_path); nbody_planet()
    gdir = joinpath(seed.out, "grid", seed.tag); dir = joinpath(seed.out, "race", seed.tag)
    cfg = cell_config(study, seed.inc, seed.raan, seed.omega; seed_alt=seed_alt_of(gdir))
    epoch = parse_epoch(study["grid"]["epochs_utc"][seed.epoch_idx])
    r = run_truth(cfg, study, epoch, Float64(seed.horizon_days) * SECONDS_PER_DAY, dir; third_bodies=true, planes=true)
    write_rows(joinpath(dir, "passes_race.csv"), with_band(r.rows, cfg, study))
    lev = try end_lever(r, cfg, study, epoch, dir, Float64(seed.fd_dv_mps)) catch err; (; lever_ok=false, lever_error=first(sprint(showerror, err), 200)) end
    return merge((; tag=seed.tag, wall_race_s=r.wall_s, retcode_race=r.retcode, passes_race=length(r.rows)), lev)
end

function parse_args(args)
    o = Dict{String, Any}("out" => "output/periapsis_static_region", "pilot" => false, "phase" => "grid", "workers" => nothing)
    for a in args
        a == "--pilot" && (o["pilot"] = true)
        startswith(a, "--out=") && (o["out"] = abspath(split(a, "=", limit=2)[2]))
        startswith(a, "--phase=") && (o["phase"] = split(a, "=", limit=2)[2])
        startswith(a, "--workers=") && (o["workers"] = parse(Int, split(a, "=", limit=2)[2]))
    end
    return o
end

report(res) = (println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1))");
    foreach(f -> println("FAILED $(f.seed): $(first(sprint(showerror, f.error), 400))"), res.failed))

function save_mc(res, out, name)
    open(joinpath(out, "$(name)_mc.toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_elapsed_s = $(res.elapsed_s)")
    end
end

function static_main(args=ARGS)
    o = parse_args(args); study = load_study(); sc = TOML.parsefile(STATIC_CONFIG_PATH); out = o["out"]; mkpath(out)
    workers = something(o["workers"], Int(sc["grid"]["workers"]))
    files = [joinpath(@__DIR__, "run_static.jl")]
    hd = Float64(sc["grid"]["horizon_days"])
    mk(c; kw...) = (; study_path=STUDY_CONFIG_PATH, out, kw...)
    if o["pilot"]
        t0 = time(); row = grid_cell_task(mk(nothing; inc=63.4, raan=0.0, omega=90.0, epoch_idx=1, horizon_days=hd))
        println("PILOT wall_s=$(round(time() - t0; digits=1)) cal=$(round(row.cal_wall_s; digits=1)) run3b=$(round(row.wall_3b_s; digits=1)) retcode=$(row.retcode_3b) passes=$(row.passes_3b)")
        return
    end
    if o["phase"] == "grid"
        seeds = [mk(nothing; inc, raan, omega, epoch_idx=e, horizon_days=hd) for inc in sc["grid"]["inclination_deg"] for e in eachindex(study["grid"]["epochs_utc"])
                 for raan in sc["grid"]["raan_deg"] for omega in sc["grid"]["omega_deg"]]
        res = run_process_campaign(grid_cell_task, seeds, workers; mission_time_s=hd * SECONDS_PER_DAY, worker_files=files)
        report(res); write_rows(joinpath(out, "static_run_summary.csv"), [s.value for s in res.successful]); save_mc(res, out, "grid")
    elseif o["phase"] == "control"
        cells = CSV.read(joinpath(out, "control_cells.csv"), DataFrame)
        seeds = [mk(nothing; tag=r.tag, inc=r.inc_deg, raan=r.raan_deg, omega=r.omega_deg, epoch_idx=Int(r.epoch_idx), horizon_days=hd) for r in eachrow(cells)]
        res = run_process_campaign(control_cell_task, seeds, workers; mission_time_s=hd * SECONDS_PER_DAY, worker_files=files)
        report(res); write_rows(joinpath(out, "control_run_summary.csv"), [s.value for s in res.successful]); save_mc(res, out, "control")
    elseif o["phase"] == "race"
        cells = CSV.read(joinpath(out, "race_cells.csv"), DataFrame)
        rd = Float64(sc["race"]["horizon_days"])
        seeds = [mk(nothing; tag=r.tag, inc=r.inc_deg, raan=r.raan_deg, omega=r.omega_deg, epoch_idx=Int(r.epoch_idx), horizon_days=rd, fd_dv_mps=Float64(sc["race"]["fd_dv_mps"])) for r in eachrow(cells)]
        res = run_process_campaign(race_cell_task, seeds, workers; mission_time_s=rd * SECONDS_PER_DAY, worker_files=files)
        report(res); write_rows(joinpath(out, "race_run_summary.csv"), [s.value for s in res.successful]); save_mc(res, out, "race")
    else
        throw(ArgumentError("--phase must be grid, control or race"))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    static_main()
end
