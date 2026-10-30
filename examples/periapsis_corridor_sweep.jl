# Periapsis corridor sweep: inclination x argument of periapsis (--atmosphere=us76_piecewise|gram, --cells=i:omega,...).
#
#   julia --project=. examples/periapsis_corridor_sweep.jl --phase=A [--output=DIR]
#   julia --project=. examples/periapsis_corridor_sweep.jl --phase=B [--margin-threshold-m=500] [--policy=frequent|near_passive]
#
# Phase A: passive coast for every cell (runner case "favorable", no controller).
# Phase B: the `--policy` controller (default `frequent`; `near_passive` writes *_near_passive files) on cells that breached or whose passive minimum
# margin to the lower bound is below the threshold. Every cell re-calibrates the
# initial periapsis to the reference altitude (geodetic) with the runner's own
# calibration path. Physics, policy and outputs are the runner's; nothing is
# duplicated here. Grid defaults come from the [sweep] table of the config.
include(joinpath(@__DIR__, "AGORA_Earth_Periapsis_Corridor.jl"))
using Distributed
# Pool workers set PC_SWEEP_WORKER=1: they only run cells, so skip the plotting stack.
get(ENV, "PC_SWEEP_WORKER", "") == "" && @eval using Plots

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

function _sweep_cli(args)
    config_path = joinpath(@__DIR__, "periapsis_corridor_config.toml")
    for arg in args
        startswith(arg, "--config=") && (config_path = abspath(split(arg, "=", limit=2)[2]))
    end
    config = load_corridor_config(config_path)
    sweep = get(config.values, "sweep", Dict{String, Any}())
    opt = Dict{String, Any}("phase" => "A", "output" => joinpath(REPO_ROOT, "output", "periapsis_corridor_sweep"),
        "inclinations" => Float64.(get(sweep, "inclinations_deg", [50.0])),
        "omegas" => Float64.(get(sweep, "argument_of_periapsis_deg", [90.0])),
        "margin" => Float64(get(sweep, "margin_threshold_m", 500.0)),
        "horizon" => nothing, "backend" => "us76_piecewise", "workers" => 4, "cells" => nothing, "output_given" => false,
        "calib" => Int(config.values["calibration"]["max_iterations"]), "policy" => :frequent, "mc_us76" => false)
    parse_list(s) = [parse(Float64, x) for x in split(s, ",")]
    for arg in args
        startswith(arg, "--") && occursin("=", arg) || continue
        key, value = split(arg[3:end], "=", limit=2)
        key == "phase" && (value in ("A", "B", "both", "plot") || throw(ArgumentError("--phase must be A, B, both or plot")); opt["phase"] = value)
        key == "output" && ((opt["output"] = abspath(value); opt["output_given"] = true))
        key == "atmosphere" && (value in ("us76_piecewise", "gram") || throw(ArgumentError("--atmosphere must be us76_piecewise or gram")); opt["backend"] = value)
        key == "policy" && (value in ("frequent", "near_passive") || throw(ArgumentError("--policy must be frequent or near_passive")); opt["policy"] = Symbol(value))
        key == "us76-process-route" && (opt["mc_us76"] = value == "true")  # --us76-process-route=true: US76 cells through run_monte_carlo, process route
        key == "workers" && (opt["workers"] = parse(Int, value))
        key == "cells" && (opt["cells"] = [Tuple(parse.(Float64, split(c, ":"))) for c in split(value, ",")])
        key == "inclinations-deg" && (opt["inclinations"] = parse_list(value))
        key == "omegas-deg" && (opt["omegas"] = parse_list(value))
        key == "margin-threshold-m" && (opt["margin"] = parse(Float64, value))
        key == "horizon-s" && (opt["horizon"] = parse(Float64, value))
        key == "calibration-max-iterations" && (opt["calib"] = parse(Int, value))
    end
    opt["backend"] == "gram" && !opt["output_given"] && (opt["output"] = joinpath(REPO_ROOT, "output", "periapsis_corridor_sweep_gram"))
    config.values["atmosphere"]["backend"] = opt["backend"]  # the runner reads and enforces it per case dir
    if opt["cells"] !== nothing  # explicit cells replace the grid; heat maps shrink to the listed cells
        opt["inclinations"] = unique(first.(opt["cells"])); opt["omegas"] = unique(last.(opt["cells"]))
    end
    return config, opt
end

_cells(opt) = opt["cells"] === nothing ? [(i, w) for i in opt["inclinations"] for w in opt["omegas"]] : opt["cells"]

function _write_manifest(opt)  # backend-stamped record of what this output directory holds
    mkpath(opt["output"])
    open(_policy_file(opt, "sweep_manifest_$(opt["phase"])", "toml"), "w") do io
        println(io, "atmosphere_backend = \"", opt["backend"], "\"\nphase = \"", opt["phase"], "\"")
        opt["phase"] == "B" && println(io, "policy = \"", opt["policy"], "\"")
        println(io, "cells = [", join(["\"$(i):$(w)\"" for (i, w) in _cells(opt)], ", "), "]")
        opt["backend"] == "gram" && println(io, "gram_density = \"mean only; deterministic; single epoch\"")
    end
end

# Phase B outputs of a non-default policy carry its name so they never overwrite or mix with `frequent`.
_policy_file(opt, stem, ext) = joinpath(opt["output"],
    opt["phase"] == "B" && opt["policy"] !== :frequent ? "$(stem)_$(opt["policy"]).$(ext)" : "$stem.$ext")

# A CSV from another backend (older US76 CSVs have no column) must never be read back as this one's.
function _read_grid(path, opt)
    df = CSV.read(path, DataFrame)
    b = "backend" in names(df) ? unique(df.backend) : ["us76_piecewise"]
    b == [opt["backend"]] || throw(ErrorException("$path holds backend $b, run is $(opt["backend"]); refusing to mix"))
    return df
end

_tag(i, w) = "i" * replace(string(i), "." => "p") * "_w" * lpad(string(round(Int, w)), 3, '0')

# Uniform row layout so every cell (success or failure) stacks into one table.
_blank(inc, omega, backend, policy="") = (; backend, policy=String(policy), i_deg=Float64(inc), omega_deg=Float64(omega), status="", reason="",
    initial_periapsis_lat_deg=NaN, final_periapsis_lat_deg=NaN, omega_first_deg=NaN, omega_last_deg=NaN,
    retcode="", elapsed_s=NaN, event_count=0, coverage_complete=false, breach_count=0,
    first_breach_time_s=NaN, min_peri_alt_m=NaN, max_peri_alt_m=NaN, min_margin_to_floor_m=NaN,
    final_apoapsis_alt_m=NaN, delivered_dv_mps=NaN, burn_count=0, propellant_used_kg=NaN,
    full_mission_verified=false, controller_success=false, controller_failures="", mission_result="", wall_s=NaN)

_df(rows) = reduce((x, y) -> vcat(x, y; cols=:union), (DataFrame([r]) for r in rows))

"""Run one cell; never throws. `policy` is :passive, :frequent or :near_passive (the runner case of the same name;
deadband and margin come from the config: [comparison] near_passive_deadband_m / frequent_deadband_m, [mission] safety_margin_m)."""
function run_cell(config, inc, omega, root, policy, opt)
    cfg = CorridorConfig(config.path, deepcopy(config.values))
    cfg.values["orbit"]["inclination_deg"] = inc
    cfg.values["orbit"]["argument_of_periapsis_deg"] = omega
    cfg.values["calibration"]["max_iterations"] = opt["calib"]
    floor_m = Float64(cfg.values["mission"]["corridor_min_altitude_m"])
    case = policy === :passive ? "favorable" : String(policy)
    dir = joinpath(root, _tag(inc, omega))
    t0 = time()
    try
        s = run_corridor_case(cfg, case, dir; horizon_s=opt["horizon"], phase_override=:favorable)
        passes = CSV.read(joinpath(dir, case, "passes.csv"), DataFrame)
        horizon_done = s.elapsed_horizon_s + 1e-6 >= something(opt["horizon"], Float64(cfg.values["mission"]["horizon_s"]))
        status = !s.retcode_success || !horizon_done ? "run_incomplete" :
            !s.event_coverage_complete ? "coverage_incomplete" :
            s.physical_breach_count > 0 ? "breach" : "pass"
        return merge(_blank(inc, omega, String(config.values["atmosphere"]["backend"]), policy), (; status,
            initial_periapsis_lat_deg=rad2deg(passes.latitude_rad[1]),
            final_periapsis_lat_deg=rad2deg(passes.latitude_rad[end]),
            omega_first_deg=rad2deg(passes.omega_rad[1]), omega_last_deg=rad2deg(passes.omega_rad[end]),
            retcode=s.solver_retcode, elapsed_s=s.elapsed_horizon_s,
            event_count=s.physical_periapsis_events, coverage_complete=s.event_coverage_complete,
            breach_count=s.physical_breach_count, first_breach_time_s=s.first_physical_breach_time_s,
            min_peri_alt_m=s.minimum_physical_periapsis_altitude_m,
            max_peri_alt_m=s.maximum_physical_periapsis_altitude_m,
            min_margin_to_floor_m=s.minimum_physical_periapsis_altitude_m - floor_m,
            final_apoapsis_alt_m=passes.physical_apoapsis_altitude_m[end],
            delivered_dv_mps=s.delivered_delta_v_mps, burn_count=s.maintenance_burn_count,
            propellant_used_kg=Float64(cfg.values["vehicle"]["propellant_mass_kg"]) - s.propellant_remaining_kg,
            full_mission_verified=s.full_mission_verified, controller_success=s.controller_success, controller_failures=s.controller_failures,
            mission_result=s.mission_result, wall_s=time() - t0))
    catch err
        msg = replace(first(sprint(showerror, err), 300), "\n" => " ")
        return merge(_blank(inc, omega, String(config.values["atmosphere"]["backend"]), policy), (;
            status=occursin("calibration", msg) ? "calibration_failed" : "run_error", reason=msg,
            wall_s=time() - t0))
    end
end

# One sweep cell as a Monte Carlo seed; a plain global function so pool workers resolve it by name.
function sweep_cell_task(seed)
    config = load_corridor_config(seed.config_path)
    config.values["atmosphere"]["backend"] = seed.backend
    return run_cell(config, seed.i, seed.omega, seed.root, seed.policy,
        Dict{String, Any}("horizon" => seed.horizon, "calib" => seed.calib))
end

"""Run `cells` and return rows in order. US76 stays serial in-process; GRAM goes through
`run_monte_carlo` forced onto the process route (native GRAM is not thread-safe, so each
worker is its own `--threads=1` OS process with its own GRAM instance)."""
function _run_cells(config, cells, root, policy, opt, on_row)
    rows = NamedTuple[]
    if opt["backend"] != "gram" && !opt["mc_us76"]  # default: serial in-process, as the original sweep ran
        for (i, w) in cells
            push!(rows, run_cell(config, i, w, root, policy, opt)); on_row(rows)
        end
        return rows
    end
    seeds = [(; i, omega=w, policy, root, horizon=opt["horizon"], calib=opt["calib"],
              config_path=config.path, backend=opt["backend"]) for (i, w) in cells]
    pool = campaign_process_pool()
    # Create the workers up front (no warmup_fn, which would run a whole cell per worker), then load
    # the example code into them; the campaign finds them already present.
    ids = ensure_process_workers!(pool, opt["workers"])
    Distributed.remotecall_eval(Main, ids, :(ENV["PC_SWEEP_WORKER"] = "1"))
    Distributed.remotecall_eval(Main, ids, :(include($(@__FILE__))))
    # density_family="gram_point" makes the router's candidate set [:none, :process] (native GRAM);
    # adaptive off returns the default (:process) with no exploration trial on :none.
    features = campaign_route_features(samples=length(seeds), density_family="gram_point",
        mission_time_s=Float64(config.values["mission"]["horizon_s"]))
    tuning = OuterRouteTuning(adaptive_enabled=false, process_workers_resident=length(ids),
        process_max_workers=opt["workers"])
    res = run_monte_carlo(sweep_cell_task, seeds; threads=:auto, route_features=features,
        route_state=OuterRouteState(), route_tuning=tuning)
    println("MC route=$(res.route) workers=$(res.threads) local_slots=$(res.local_slots) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time)")
    open(_policy_file(opt, "sweep_manifest_$(opt["phase"])", "toml"), "a") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_elapsed_s = $(res.elapsed_s)")
    end
    for smp in res.samples  # failed sample: keep a row so the cell is not silently dropped
        push!(rows, smp.success ? smp.value : merge(_blank(smp.seed.i, smp.seed.omega, smp.seed.backend, smp.seed.policy),
            (; status="run_error", reason=replace(first(sprint(showerror, smp.error), 300), "\n" => " "))))
    end
    on_row(rows)
    res.route === :process || throw(ErrorException("GRAM cells ran on route $(res.route), not :process"))
    return rows
end

function _heat(df, col, grid_i, grid_w, path, title; mark_breach=true, band=(95_000.0, 100_000.0))
    grid_i, grid_w = sort(unique(grid_i)), sort(unique(grid_w))
    z = fill(NaN, length(grid_i), length(grid_w))
    for r in eachrow(df)
        a, b = findfirst(==(r.i_deg), grid_i), findfirst(==(r.omega_deg), grid_w)
        (a === nothing || b === nothing) && continue
        z[a, b] = r[col]
    end
    plt = heatmap(1:length(grid_w), 1:length(grid_i), z; xticks=(1:length(grid_w), string.(Int.(round.(grid_w)))),
        yticks=(1:length(grid_i), string.(grid_i)), xlabel="argument of periapsis (deg)",
        ylabel="inclination (deg)", title=title, color=:viridis, size=(900, 560),
        titlefontsize=10, left_margin=6Plots.mm, bottom_margin=4Plots.mm)
    for r in eachrow(df)
        a, b = findfirst(==(r.i_deg), grid_i), findfirst(==(r.omega_deg), grid_w)
        (a === nothing || b === nothing) && continue
        txt = r.status == "calibration_failed" ? "CAL FAIL" :
            r.status == "run_error" ? "ERROR" :
            mark_breach && r.status == "breach" ?
                "BREACH " * (r.min_peri_alt_m < band[1] ? "LOW" : "HIGH") * "\n" * string(round(Int, r[col])) :
            isfinite(r[col]) ? string(round(r[col]; digits=abs(r[col]) < 1 ? 3 : 0)) : "n/a"
        color = r.status in ("breach", "calibration_failed", "run_error") ? :red : :white
        annotate!(plt, b, a, Plots.text(txt, 7, color))
    end
    savefig(plt, path)
end

function phase_a(config, opt)
    root = joinpath(opt["output"], "cells_passive")
    rows = NamedTuple[]
    csv_path = joinpath(opt["output"], "passive_grid.csv")
    mkpath(opt["output"])
    _write_manifest(opt)
    rows = _run_cells(config, _cells(opt), root, :passive, opt, rs -> CSV.write(csv_path, _df(rs)))
    for row in rows
        println("A i=$(row.i_deg) omega=$(row.omega_deg) status=$(row.status) min_margin_m=$(row.min_margin_to_floor_m) wall_s=$(round(row.wall_s; digits=1))")
    end
    df = _df(rows)
    label = opt["horizon"] === nothing ? "" : "\n[reduced_horizon_diagnostic, $(opt["horizon"]) s]"
    _heat(df, :min_margin_to_floor_m, opt["inclinations"], opt["omegas"],
        joinpath(opt["output"], "passive_margin_heatmap.png"),
        "Passive coast, $(opt["backend"]): minimum periapsis margin to the 95 km geodetic floor (m)\n" *
        "BREACH LOW = below 95 km, HIGH = above 100 km (number shown is still the margin to the floor)" * label;
        band=(Float64(config.values["mission"]["corridor_min_altitude_m"]),
              Float64(config.values["mission"]["corridor_max_altitude_m"])))
    return df
end

function phase_b(config, opt)
    csv_path = joinpath(opt["output"], "passive_grid.csv")
    isfile(csv_path) || throw(ErrorException("Phase B needs $csv_path from Phase A"))
    passive = _read_grid(csv_path, opt)
    _write_manifest(opt)
    # explicit --cells run exactly those cells (as passive row context); otherwise the margin rule selects
    select_mask = [opt["cells"] !== nothing ? (r.i_deg, r.omega_deg) in opt["cells"] :
                   r.status in ("breach", "pass") && (r.breach_count > 0 || r.min_margin_to_floor_m < opt["margin"])
                   for r in eachrow(passive)]
    chosen = passive[select_mask, :]
    println("Phase B: $(nrow(chosen)) of $(nrow(passive)) cells below margin threshold $(opt["margin"]) m or breaching")
    root = joinpath(opt["output"], opt["policy"] === :frequent ? "cells_controlled" : "cells_controlled_$(opt["policy"])")
    rows = NamedTuple[]
    controlled_path = _policy_file(opt, "controlled_grid", "csv")
    raw = _run_cells(config, [(p.i_deg, p.omega_deg) for p in eachrow(chosen)], root, opt["policy"], opt, _ -> nothing)
    for (p, row) in zip(eachrow(chosen), raw)
        # full_mission_verified is false by construction at a reduced horizon, so judge pass + controller success
        feasible = row.status == "pass" && row.controller_success
        cause = feasible ? "" : row.status != "pass" ? "$(row.status): $(row.reason)$(row.mission_result) $(row.controller_failures)" :
            "not verified: $(row.mission_result) $(row.controller_failures)"
        push!(rows, merge(row, (; passive_status=p.status, passive_min_margin_m=p.min_margin_to_floor_m,
            passive_breach_count=p.breach_count, classification=feasible ? "feasible" : "infeasible",
            infeasible_cause=cause)))
        println("B i=$(p.i_deg) omega=$(p.omega_deg) $(rows[end].classification) dv=$(row.delivered_dv_mps) burns=$(row.burn_count) wall_s=$(round(row.wall_s; digits=1))")
    end
    CSV.write(controlled_path, _df(rows))
    isempty(rows) && return DataFrame()
    df = _df(rows)
    label = opt["horizon"] === nothing ? "" : "\n[reduced_horizon_diagnostic]"
    _heat(df, :delivered_dv_mps, opt["inclinations"], opt["omegas"],
        _policy_file(opt, "controlled_dv_heatmap", "png"),
        "$(opt["policy"]) policy: delivered delta-v (m/s), cells run in Phase B only" * label; mark_breach=false)
    return df
end

function sweep_main(args=ARGS)
    config, opt = _sweep_cli(args)
    if opt["phase"] == "plot"  # redraw the Phase A heat map from an existing passive_grid.csv
        df = _read_grid(joinpath(opt["output"], "passive_grid.csv"), opt)
        _heat(df, :min_margin_to_floor_m, opt["inclinations"], opt["omegas"],
            joinpath(opt["output"], "passive_margin_heatmap.png"),
            "Passive coast, $(opt["backend"]): minimum periapsis margin to the 95 km geodetic floor (m)\n" *
            "BREACH LOW = below 95 km, HIGH = above 100 km (number shown is the margin to the floor)")
    end
    opt["phase"] in ("A", "both") && phase_a(config, opt)
    opt["phase"] in ("B", "both") && phase_b(config, opt)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    sweep_main()
end
