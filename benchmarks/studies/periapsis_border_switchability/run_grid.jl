# Part B: passive 60-day grid, with and without Sun+Moon.
#   julia --project=. benchmarks/studies/periapsis_border_switchability/run_grid.jl --pilot [--out=DIR]
#   julia --project=. benchmarks/studies/periapsis_border_switchability/run_grid.jl [--out=DIR] [--cells=raan:omega:epochIndex,...] [--workers=N]
include(joinpath(@__DIR__, "common.jl"))

cell_tag(raan, omega, epoch_idx) = "R$(lpad(round(Int, raan), 3, '0'))_W$(lpad(round(Int, omega), 3, '0'))_E$(epoch_idx)"

function write_rows(path, rows)
    mkpath(dirname(path)); CSV.write(path, DataFrame(rows))
end

"""One grid cell: calibrate, then run with and without Sun+Moon from the same state. Returns a plain summary row."""
function grid_cell_task(seed)
    study = load_study(seed.study_path)
    cfg0 = make_study_config(study; raan_deg=seed.raan, omega_deg=seed.omega)
    epoch = parse_epoch(study["grid"]["epochs_utc"][seed.epoch_idx])
    dir = joinpath(seed.out, "grid", cell_tag(seed.raan, seed.omega, seed.epoch_idx))
    nbody_planet()  # furnish SPICE kernels
    t0 = time()
    cal, cal_records = calibrate_peak_q(cfg0, study, epoch, dir; third_bodies=true)
    cal_wall = time() - t0
    horizon = Float64(study["grid"]["horizon_days"]) * SECONDS_PER_DAY
    r3 = run_truth(cal, study, epoch, horizon, dir; third_bodies=true, planes=true)
    rc = run_truth(cal, study, epoch, horizon, dir; third_bodies=false, planes=true)
    write_rows(joinpath(dir, "passes_3b.csv"), r3.rows)
    write_rows(joinpath(dir, "passes_ctl.csv"), rc.rows)
    write_rows(joinpath(dir, "calibration.csv"), cal_records)
    s3, sc = band_stats(r3.rows, study), band_stats(rc.rows, study)
    t_end = min(r3.rows[end].time_s, rc.rows[end].time_s)  # last pass common to both runs (just under the horizon)
    d3 = interp_col(r3.rows, :alt_m, t_end) - r3.rows[1].alt_m
    dc = interp_col(rc.rows, :alt_m, t_end) - rc.rows[1].alt_m
    r3rad = interp_col(r3.rows, :radius_m, t_end) - r3.rows[1].radius_m
    rcrad = interp_col(rc.rows, :radius_m, t_end) - rc.rows[1].radius_m
    kozai = sum((r.moon_kozai_drp_mps + r.sun_kozai_drp_mps) * ((j < length(r3.rows) ? r3.rows[j + 1].time_s : t_end) - r.time_s)
                for (j, r) in enumerate(r3.rows) if r.time_s <= t_end; init=0.0)
    avg(col) = sum(getfield(r, col) for r in r3.rows) / length(r3.rows)
    return (; tag=cell_tag(seed.raan, seed.omega, seed.epoch_idx), raan_deg=seed.raan, omega_deg=seed.omega, epoch_idx=seed.epoch_idx,
        epoch=study["grid"]["epochs_utc"][seed.epoch_idx], cal_iterations=length(cal_records), cal_wall_s=cal_wall,
        initial_periapsis_alt_m=r3.rows[1].alt_m, initial_vrel_mps=r3.rows[1].vrel_mps, initial_q=r3.rows[1].q_w_cm2,
        wall_3b_s=r3.wall_s, wall_ctl_s=rc.wall_s, retcode_3b=r3.retcode, retcode_ctl=rc.retcode,
        t_end_days=t_end / SECONDS_PER_DAY, passes_3b=s3.passes, passes_ctl=sc.passes,
        drp_geodetic_3b_m=d3, drp_geodetic_ctl_m=dc, third_body_effect_geodetic_m=d3 - dc,
        drp_radius_3b_m=r3rad, drp_radius_ctl_m=rcrad, third_body_effect_radius_m=r3rad - rcrad,
        kozai_predicted_radius_change_m=kozai,
        mean_sin2wp_moon=sum(sin(2r.moon_omega_p_rad) for r in r3.rows) / length(r3.rows),
        mean_sin2wp_sun=sum(sin(2r.sun_omega_p_rad) for r in r3.rows) / length(r3.rows),
        omega_p_moon_start_deg=rad2deg(r3.rows[1].moon_omega_p_rad), omega_p_moon_end_deg=rad2deg(r3.rows[end].moon_omega_p_rad),
        omega_p_sun_start_deg=rad2deg(r3.rows[1].sun_omega_p_rad), omega_p_sun_end_deg=rad2deg(r3.rows[end].sun_omega_p_rad),
        i_p_moon_start_deg=rad2deg(r3.rows[1].moon_i_p_rad), i_p_sun_start_deg=rad2deg(r3.rows[1].sun_i_p_rad),
        min_alt_3b_m=minimum(r.alt_m for r in r3.rows), max_alt_3b_m=maximum(r.alt_m for r in r3.rows),
        min_alt_ctl_m=minimum(r.alt_m for r in rc.rows), max_alt_ctl_m=maximum(r.alt_m for r in rc.rows),
        apo_start_m=r3.rows[1].apo_alt_m, apo_end_3b_m=r3.rows[end].apo_alt_m, apo_end_ctl_m=rc.rows[end].apo_alt_m,
        in_band_passes_3b=s3.passes_in_band, in_band_days_3b=s3.time_in_band_days, first_exit_days_3b=s3.first_exit_time_days, first_exit_edge_3b=s3.first_exit_edge,
        q_max_3b=s3.q_max, q_min_3b=s3.q_min, exceeds_hard_3b=s3.exceeds_hard,
        in_band_passes_ctl=sc.passes_in_band, in_band_days_ctl=sc.time_in_band_days, first_exit_days_ctl=sc.first_exit_time_days, first_exit_edge_ctl=sc.first_exit_edge,
        q_max_ctl=sc.q_max, q_min_ctl=sc.q_min)
end

function parse_args(args)
    o = Dict{String, Any}("out" => nothing, "pilot" => false, "cells" => nothing, "workers" => nothing)
    for a in args
        a == "--pilot" && (o["pilot"] = true)
        startswith(a, "--out=") && (o["out"] = abspath(split(a, "=", limit=2)[2]))
        startswith(a, "--workers=") && (o["workers"] = parse(Int, split(a, "=", limit=2)[2]))
        startswith(a, "--cells=") && (o["cells"] = [Tuple(parse.(Float64, split(c, ":"))) for c in split(split(a, "=", limit=2)[2], ",")])
    end
    return o
end

function grid_main(args=ARGS)
    o = parse_args(args)
    study = load_study()
    o["out"] = something(o["out"], results_root(study)); mkpath(o["out"])
    cells = o["pilot"] ? [(0.0, 90.0, 1.0)] : o["cells"] !== nothing ? o["cells"] :
        [(r, w, Float64(e)) for e in eachindex(study["grid"]["epochs_utc"]) for r in study["grid"]["raan_deg"] for w in study["grid"]["omega_deg"]]
    seeds = [(; study_path=STUDY_CONFIG_PATH, out=o["out"], raan=c[1], omega=c[2], epoch_idx=Int(c[3])) for c in cells]
    if o["pilot"]  # one cell, serial in-process, to time it and surface any SPICE problem
        t0 = time(); row = grid_cell_task(seeds[1])
        println("PILOT wall_s=$(round(time() - t0; digits=1)) (3b run $(round(row.wall_3b_s; digits=1)) s, control $(round(row.wall_ctl_s; digits=1)) s, calibration $(round(row.cal_wall_s; digits=1)) s)")
        display(pairs(row)); println()
        write_rows(joinpath(o["out"], "pilot_row.csv"), [row]); return nothing
    end
    workers = something(o["workers"], Int(study["grid"]["workers"]))
    res = run_process_campaign(grid_cell_task, seeds, workers; mission_time_s=Float64(study["grid"]["horizon_days"]) * SECONDS_PER_DAY,
        worker_files=[joinpath(@__DIR__, "run_grid.jl")])
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1))")
    for f in res.failed
        println("FAILED $(f.seed): $(first(sprint(showerror, f.error), 400))")
    end
    write_rows(joinpath(o["out"], "grid_summary.csv"), [s.value for s in res.successful])
    open(joinpath(o["out"], "grid_mc.toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_elapsed_s = $(res.elapsed_s)")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    grid_main()
end
