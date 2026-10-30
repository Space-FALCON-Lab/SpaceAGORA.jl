# Earth truth at higher gravity degree/order (user directive 2026-10-17: every truth study uses at least 20x20 harmonics).
# Same cells, same calibration rule, same everything else as run_static.jl / run_grid.jl; only the truth gravity model differs
# (common.jl `build_args` honors an optional [gravity] table of the study config). One Julia job at a time, memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. run_earth_g20.jl --phase=cells|pilot|grid --results-root=DIR [--cells=CSV] [--epochs=ISO,ISO,ISO] [--subdir=NAME] [--only=tag,tag] [--workers=N] [--degree=20] [--order=20]
#   cells   writes <root>/earth_cells_336.csv (the 336 cells of the earlier Earth truth: 50 deg border grid + 288 static-sweep cells)
#   pilot   one cell, serial, timed
#   grid    cells of --cells (default earth_cells_336.csv) on the process route (run_monte_carlo, route asserted)  -> <root>/earth_g20/grid/<tag>/{calibration,passes_3b}.csv
# The study config used by the workers is <root>/study_config_g20.toml = study_config.toml + a [gravity] table (written here, not committed).
include(joinpath(@__DIR__, "run_static.jl"))
using Statistics

function g20_cli()
    kv = Dict{String, String}()
    for a in ARGS
        m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2])
    end
    return kv
end
const G20 = g20_cli()

function earth_cells_336()
    study = load_study(); sc = TOML.parsefile(STATIC_CONFIG_PATH); rows = NamedTuple[]
    for inc in vcat(Float64(study["grid"]["inclination_deg"]), Float64.(sc["grid"]["inclination_deg"])), e in eachindex(study["grid"]["epochs_utc"])
        raans = inc == 50.0 ? study["grid"]["raan_deg"] : sc["grid"]["raan_deg"]; omegas = inc == 50.0 ? study["grid"]["omega_deg"] : sc["grid"]["omega_deg"]
        for r in raans, w in omegas
            push!(rows, (tag=cell_tag(inc, r, w, e), inc_deg=inc, raan_deg=Float64(r), omega_deg=Float64(w), epoch_idx=e))
        end
    end
    return rows
end

function write_g20_study(root, degree, order)
    path = joinpath(root, "study_config_g20.toml")
    txt = read(STUDY_CONFIG_PATH, String)
    if haskey(G20, "epochs")   # blind cells: replace the three start epochs (comma-separated ISO strings)
        txt = replace(txt, r"epochs_utc = \[[^\]]*\]" => "epochs_utc = [" * join(["\"$x\"" for x in split(G20["epochs"], ",")], ", ") * "]")
        occursin(split(G20["epochs"], ",")[1], txt) || error("epochs_utc line not replaced")
    end
    txt *= "\n[gravity]\n# Truth gravity override read by common.jl build_args (absent in study_config.toml: degree 3 order 0 of the example).\n" *
        "coefficient_file = \"data/planets/gravity/EarthGGM05C.csv\"\ndegree = $degree\norder = $order\n"
    mkpath(root); write(path, txt); return path
end

function g20_main()
    root = abspath(get(G20, "results-root", results_root(load_study()))); phase = get(G20, "phase", "cells")
    degree = parse(Int, get(G20, "degree", "20")); order = parse(Int, get(G20, "order", "20")); mkpath(root)
    cells_csv = abspath(get(G20, "cells", joinpath(root, "earth_cells_336.csv")))
    if phase == "cells"
        write_rows(joinpath(root, "earth_cells_336.csv"), earth_cells_336()); println("wrote $(length(earth_cells_336())) cells"); return
    end
    spath = write_g20_study(root, degree, order); cells = CSV.read(cells_csv, DataFrame)
    haskey(G20, "only") && (cells = cells[in.(cells.tag, Ref(split(G20["only"], ","))), :])
    out = joinpath(root, get(G20, "subdir", "earth_g20")); mkpath(out); hd = 60.0
    seeds = [(; study_path=spath, out, inc=r.inc_deg, raan=r.raan_deg, omega=r.omega_deg, epoch_idx=Int(r.epoch_idx), horizon_days=hd) for r in eachrow(cells)]
    if phase == "pilot"
        for s in seeds
            t0 = time(); row = grid_cell_task(s); println("PILOT $(row.tag) wall_s=$(round(time() - t0; digits=1)) cal=$(round(row.cal_wall_s; digits=1)) run3b=$(round(row.wall_3b_s; digits=1)) retcode=$(row.retcode_3b) passes=$(row.passes_3b)")
        end
        return
    end
    workers = parse(Int, get(G20, "workers", "4"))
    res = run_process_campaign(grid_cell_task, seeds, workers; mission_time_s=hd * SECONDS_PER_DAY, worker_files=[joinpath(@__DIR__, "run_static.jl")])
    report(res); write_rows(joinpath(out, "grid_run_summary.csv"), [s.value for s in res.successful]); save_mc(res, out, "grid")
end
g20_main()
