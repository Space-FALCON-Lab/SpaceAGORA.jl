# Stage F at 20x20 driver (copy of stagef_run.jl reading stagef20_config.toml; output below <root>/stagef20/). Stage F driver: ONE run_monte_carlo campaign per phase over all planets together (process route asserted). One Julia job, memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. stagef_run.jl --phase=feasibility|cal|f1|f2|selftest --results-root=DIR [--workers=N] [--groups=venus:1,mars:2,...]
#        [--planets=venus,mars,earth] [--host=NAME] [--horizon-days=D] [--only=tag,tag] [--ctrls=open,fixed] [--cells-dir=DIR]
# Rows are stamped with --host and written to <root>/stagef/summary_<phase>_<host>.csv. Plan: tasks/2026-10-15_periapsis-region-border-targeting/stagef20_plan.md.
#   cal  Pi-rule altitude, start fit, mean-density calibration per (planet, cell)           -> <root>/stagef/<planet>/cal/<tag>/
#   f1   per cell: {open, fixed} x {mean, pert} runs (paired mean/perturbed share the calibrated start)  -> <root>/stagef/<planet>/runs/<ctrl>_<arm>_<tag>/passes.csv
#   f2   adaptive-margin arm (pert) on held-out cells with the dev-derived parameters in <root>/stagef/<planet>/adapt_params.csv -> runs/adaptive_pert_<tag>/
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "stagef_common.jl"))
const CFG20 = joinpath(@__DIR__, "stagef20_config.toml")   # Stage F at 20x20: the only difference to stagef_run.jl is that every config read goes to this file

function parse_cli()
    kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2]); end
    return kv
end
const CLI = parse_cli(); const META = stagef_meta(CFG20)
const ROOT = abspath(CLI["results-root"]); const SF = joinpath(ROOT, String(META["stagef"]["results_subdir"]))
const HOST = get(CLI, "host", "local")
const PLANETS = haskey(CLI, "planets") ? String.(split(CLI["planets"], ",")) : String.(META["stagef"]["planets"])
const GROUPS = haskey(CLI, "groups") ? Set(String.(split(CLI["groups"], ","))) : nothing
ingroup(pl, e) = GROUPS === nothing || "$(pl):$(e)" in GROUPS
workers() = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(META["stagef"]["workers"])
config_path() = CFG20
const CTRLS = haskey(CLI, "ctrls") ? String.(split(CLI["ctrls"], ",")) : ["open", "fixed"]   # --ctrls=fixed restricts f1 (cross-host check)

function cells_of(pl)
    D = CSV.read(joinpath(get(CLI, "cells-dir", SF), pl, "cells_base.csv"), DataFrame)   # --cells-dir: where stagef_cells.py wrote the cell lists (remote hosts have no pandas)
    haskey(CLI, "only") && (D = D[in.(D.tag, Ref(split(CLI["only"], ","))), :])
    return D[[ingroup(pl, e) for e in D.epoch_idx], :]
end

"""Per-planet feasibility: construct the planet, the GRAM model and sample a density. A planet that fails is reported with the exact error and dropped."""
function feasible_planets()
    ok = String[]
    for pl in PLANETS
        try
            sc = load_stagef(pl; path=CFG20); ms = stagef_setup(sc); epoch = parse_epoch(sc["grid"]["epochs_utc"][1]); model = stagef_gram_model(sc, epoch, 1001)
            pargs = gram_pargs(ms.planet, 2000.0); rho = gram_rho(model, pargs, 130e3, 0.3, 1.0, 0.0); isfinite(rho) && rho > 0 || error("non-positive density $rho")
            println("FEASIBLE $pl: GRAM mean density at 130 km = $rho kg/m3"); push!(ok, pl)
        catch err
            println("INFEASIBLE $pl: ", sprint(showerror, err)); 
        end
    end
    return ok
end

function finish(res, phase)
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time, not a performance number)")
    foreach(f -> println("FAILED $(f.seed.run_id): $(first(sprint(showerror, f.error), 600))"), res.failed)
    rows = [merge(s.value, (; host=HOST)) for s in res.successful]
    write_rows_m(joinpath(SF, "summary_$(phase)_$(HOST).csv"), rows)
    open(joinpath(SF, "mc_$(phase)_$(HOST).toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_density_family = \"gram\"\nhost = \"$HOST\"")
    end
    bad = count(r -> r.status != "ok", rows); println("rows=$(length(rows)) not ok=$bad")
    for r in rows; r.status != "ok" && println("  $(r.planet) $(get(r, :run_id, r.tag)): $(r.error)"); end
end

campaign(task, seeds, mission_s) = run_gram_campaign(task, seeds, workers(); mission_time_s=mission_s, worker_files=[joinpath(@__DIR__, "stagef_common.jl")])

function phase_cal(pls)
    seeds = NamedTuple[]
    for pl in pls, r in eachrow(cells_of(pl))
        push!(seeds, (; config_path=config_path(), planet=pl, tag=r.tag, run_id="cal_$(pl)_$(r.tag)", inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch_idx), out=joinpath(SF, pl, "cal")))
    end
    println("cal: $(length(seeds)) cases over $(pls)")
    finish(campaign(stagef_cal_task, seeds, 0.65 * 86400.0), "cal")
end

hd(pl) = haskey(CLI, "horizon-days") ? parse(Float64, CLI["horizon-days"]) : Float64(load_stagef(pl; path=CFG20)["grid"]["horizon_days"])
function run_seed(pl, r, arm, ctrl, extra=(;); label=ctrl)
    merge((; config_path=config_path(), planet=pl, tag=r.tag, inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch_idx), arm, ctrl, run_id="$(label)_$(arm)_$(r.tag)",
        out=joinpath(SF, pl, "runs"), cal_root=joinpath(SF, pl, "cal"), horizon_days=hd(pl), adapt_pred="", adapt_e_hi=0.0, adapt_e_lo=0.0), extra)
end
calibrated(pl, tag) = isfile(joinpath(SF, pl, "cal", tag, "guess.csv"))

function phase_f1(pls)
    seeds = NamedTuple[]
    for pl in pls, r in eachrow(cells_of(pl))
        calibrated(pl, r.tag) || (println("skip uncalibrated $pl $(r.tag)"); continue)
        for ctrl in CTRLS, arm in ("mean", "pert"); push!(seeds, run_seed(pl, r, arm, ctrl)); end
    end
    println("f1: $(length(seeds)) cases over $(pls)")
    finish(campaign(stagef_run_task, seeds, maximum(hd.(pls)) * 86400.0), "f1")
end

function phase_f2(pls)
    seeds = NamedTuple[]
    for pl in pls
        P = CSV.read(joinpath(SF, pl, "adapt_params.csv"), DataFrame)[1, :]
        for r in eachrow(cells_of(pl))
            r.split == "heldout" && calibrated(pl, r.tag) || continue
            # --variant=track (UNFROZEN diagnostic, plan Addendum F): same predictor, quantile offsets zero (density tracking only, no buffer); run ids prefixed "track"
            track = get(CLI, "variant", "") == "track"
            push!(seeds, run_seed(pl, r, "pert", "adaptive", (; adapt_pred=String(P.predictor), adapt_e_hi=track ? 0.0 : Float64(P.e_hi), adapt_e_lo=track ? 0.0 : Float64(P.e_lo)); label=track ? "track" : "adaptive"))
        end
    end
    println("f2: $(length(seeds)) cases over $(pls)")
    finish(campaign(stagef_run_task, seeds, maximum(hd.(pls)) * 86400.0), get(CLI, "variant", "") == "track" ? "f2track" : "f2")
end

function selftest()
    ys = [0.0, 0.1, 0.25, 0.2, 0.4, 0.5, 0.45, 0.7, 0.9, 1.0, 1.1, 1.3]; lons = cumsum(fill(0.1, 12)) .+ 0.01 .* sin.(1:12)
    for n in (1, 3, 8, 12)
        println("n=$n ", join(["$(p)=$(round(stagef_predict_next(p, ys[1:n], lons[1:n]); digits=12))" for p in STAGEF_PREDICTORS], " "))
    end
end

function main()
    mkpath(SF); phase = get(CLI, "phase", "")
    phase == "selftest" && return selftest()
    pls = feasible_planets(); phase == "feasibility" && return
    isempty(pls) && error("no feasible planet")
    phase == "cal" && phase_cal(pls); phase == "f1" && phase_f1(pls); phase == "f2" && phase_f2(pls)
end
main()
