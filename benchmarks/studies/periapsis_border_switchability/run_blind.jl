# Part 2: blind truth runs of the estimator-v1 picks (45 cells from blind_cells.csv, written by select_blind.py).
# Same task as the static sweep (run_static.jl grid_cell_task: calibrate to the band-center peak Q, then a 60-day run with
# Sun+Moon via SPICE), on the process route through run_monte_carlo.
#   julia --project=. run_blind.jl [--pilot] [--results-root=DIR] [--workers=N]
# Reads <root>/blind/blind_cells.csv; writes <root>/blind/grid/<tag>/{calibration,passes_3b}.csv, blind_run_summary.csv, blind_mc.toml.
include(joinpath(@__DIR__, "run_static.jl"))

function blind_cli(args)
    o = Dict{String, Any}("root" => nothing, "pilot" => false, "workers" => nothing, "only" => nothing)
    for a in args
        a == "--pilot" && (o["pilot"] = true)
        startswith(a, "--results-root=") && (o["root"] = abspath(split(a, "=", limit=2)[2]))
        startswith(a, "--workers=") && (o["workers"] = parse(Int, split(a, "=", limit=2)[2]))
        startswith(a, "--only=") && (o["only"] = split(split(a, "=", limit=2)[2], ","))
    end
    return o
end

function blind_main(args=ARGS)
    o = blind_cli(args); study = load_study(); bc = TOML.parsefile(joinpath(@__DIR__, "blind_config.toml"))
    root = something(o["root"], results_root(study)); out = joinpath(root, "blind"); mkpath(out)
    cells = CSV.read(joinpath(out, "blind_cells.csv"), DataFrame)
    o["only"] === nothing || (cells = cells[in.(cells.tag, Ref(o["only"])), :])
    hd = Float64(bc["run"]["horizon_days"]); workers = something(o["workers"], Int(bc["run"]["workers"]))
    mk(r) = (; study_path=STUDY_CONFIG_PATH, out, inc=r.inc_deg, raan=r.raan_deg, omega=r.omega_deg, epoch_idx=Int(r.epoch), horizon_days=hd)
    seeds = [mk(r) for r in eachrow(cells)]
    if o["pilot"]
        t0 = time(); row = grid_cell_task(first(seeds))
        println("PILOT tag=$(row.tag) wall_s=$(round(time() - t0; digits=1)) cal=$(round(row.cal_wall_s; digits=1)) run3b=$(round(row.wall_3b_s; digits=1)) retcode=$(row.retcode_3b) passes=$(row.passes_3b)")
        return
    end
    res = run_process_campaign(grid_cell_task, seeds, workers; mission_time_s=hd * SECONDS_PER_DAY, worker_files=[joinpath(@__DIR__, "run_static.jl")])
    report(res); write_rows(joinpath(out, "blind_run_summary.csv"), [s.value for s in res.successful]); save_mc(res, out, "blind")
end

if abspath(PROGRAM_FILE) == @__FILE__
    blind_main()
end
