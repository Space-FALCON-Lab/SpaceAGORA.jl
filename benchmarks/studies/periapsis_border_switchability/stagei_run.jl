# Stage I driver: the paper corridor (stagei_config.toml). One run_monte_carlo campaign per call (process route asserted); memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. stagei_run.jl --phase=cal|dev|held|selftest --results-root=DIR [--workers=N] [--planets=mars,earth,venus] [--host=NAME] [--only=tag,tag] [--horizon-days=D]
#   cal   Pi-rule altitude, start fit and first-pass calibration to Q = 0.15 for every cell (dev and held out)        -> <root>/stagei/<planet>/cal/<tag>/
#   dev   fixed 1000 m keeper, perturbed arm, on the development cells (feeds the e_hi / e_lo refit, stagei_analysis.py --phase=params) -> runs/fixed_pert_<tag>/
#   held  held-out cells: fixed, Stage F adaptive, Stage G asymmetric, Stage H (e_lo_frac) arms, perturbed, from <root>/stagei/<planet>/adapt_params.csv -> runs/<ctrl>_pert_<tag>/
# Rows: <root>/stagei/summary_<phase>_<host>.csv. Plan: tasks/2026-10-15_periapsis-region-border-targeting/stagei_plan.md.
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "stagei_common.jl"))
const CFGI = joinpath(@__DIR__, "stagei_config.toml")

function parse_cli()
    kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2]); end
    return kv
end
const CLI = parse_cli(); const META = stagef_meta(CFGI)
const ROOT = abspath(get(CLI, "results-root", ".")); const SF = joinpath(ROOT, String(META["stagef"]["results_subdir"]))
const HOST = get(CLI, "host", "local")
const PLANETS = haskey(CLI, "planets") ? String.(split(CLI["planets"], ",")) : String.(META["stagef"]["planets"])
workers() = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(META["stagef"]["workers"])
hd(pl) = haskey(CLI, "horizon-days") ? parse(Float64, CLI["horizon-days"]) : Float64(load_stagef(pl; path=CFGI)["grid"]["horizon_days"])
const FRAC = Float64(META["stagei"]["e_lo_frac"])
const HCTRL = "adaptive_h" * lpad(string(round(Int, 100 * FRAC)), 3, '0')

function cells_of(pl)
    D = CSV.read(joinpath(SF, pl, "cells_base.csv"), DataFrame)
    haskey(CLI, "only") && (D = D[in.(D.tag, Ref(split(CLI["only"], ","))), :])
    return D
end
calibrated(pl, tag) = isfile(joinpath(SF, pl, "cal", tag, "guess.csv"))

function finish(res, phase)
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time, not a performance number)")
    foreach(f -> println("FAILED $(f.seed.run_id): $(first(sprint(showerror, f.error), 600))"), res.failed)
    rows = [merge(s.value, (; host=HOST)) for s in res.successful]
    write_rows_m(joinpath(SF, "summary_$(phase)_$(HOST).csv"), rows)
    open(joinpath(SF, "mc_$(phase)_$(HOST).toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_density_family = \"gram\"\nhost = \"$HOST\"")
    end
    res.route === :process || error("route was $(res.route), not :process")
    bad = count(r -> r.status != "ok", rows); println("rows=$(length(rows)) not ok=$bad")
    for r in rows; r.status != "ok" && println("  $(r.planet) $(get(r, :run_id, r.tag)): $(r.error)"); end
end
campaign(task, seeds, mission_s) = run_gram_campaign(task, seeds, workers(); mission_time_s=mission_s, worker_files=[joinpath(@__DIR__, "stagei_common.jl")])

function phase_cal(pls)
    seeds = NamedTuple[]
    for pl in pls, r in eachrow(cells_of(pl))
        push!(seeds, (; config_path=CFGI, planet=pl, tag=r.tag, run_id="cal_$(pl)_$(r.tag)", inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch_idx), out=joinpath(SF, pl, "cal")))
    end
    println("cal: $(length(seeds)) cases over $(pls)")
    finish(campaign(stagei_cal_task, seeds, 0.65 * 86400.0), "cal")
end

seedf(pl, r, ctrl, extra=(;)) = merge((; config_path=CFGI, planet=pl, tag=r.tag, inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch_idx), arm="pert", ctrl, run_id="$(ctrl)_pert_$(r.tag)",
    out=joinpath(SF, pl, "runs"), cal_root=joinpath(SF, pl, "cal"), horizon_days=hd(pl), adapt_pred="", adapt_e_hi=0.0, adapt_e_lo=0.0), extra)

function phase_dev(pls)
    seeds = NamedTuple[]
    for pl in pls, r in eachrow(cells_of(pl))
        r.split == "dev" && calibrated(pl, r.tag) && push!(seeds, seedf(pl, r, "fixed"))
    end
    println("dev: $(length(seeds)) cases over $(pls)")
    finish(campaign(stagei_task, seeds, maximum(hd.(pls)) * 86400.0), "dev")
end

function phase_held(pls)
    seeds = NamedTuple[]
    for pl in pls
        P = CSV.read(joinpath(SF, pl, "adapt_params.csv"), DataFrame)[1, :]   # Stage I development refit (stagei_analysis.py --phase=params); never the Stage F20 values
        pred = String(P.predictor); ehi = Float64(P.e_hi); g = FRAC * abs(Float64(P.e_lo))
        for r in eachrow(cells_of(pl))
            r.split == "heldout" && calibrated(pl, r.tag) || continue
            push!(seeds, seedf(pl, r, "fixed"))
            push!(seeds, seedf(pl, r, "adaptive", (; adapt_pred=pred, adapt_e_hi=ehi, adapt_e_lo=Float64(P.e_lo))))
            push!(seeds, seedf(pl, r, "adaptive_asym", (; adapt_pred=pred, adapt_e_hi=ehi, adapt_e_lo=0.0)))
            push!(seeds, seedf(pl, r, HCTRL, (; adapt_pred=pred, adapt_e_hi=ehi, adapt_e_lo=g)))
        end
    end
    println("held: $(length(seeds)) cases over $(pls)")
    finish(campaign(stagei_task, seeds, maximum(hd.(pls)) * 86400.0), "held")
end

"""Exact-case check of the calibration view: the geometric mean of (q_low, q_high) of the view equals the paper target, and the real band is untouched."""
function selftest()
    sc = load_stagef("mars"; path=CFGI); b = sc["band"]; v = stagei_cal_view(sc)
    @assert b["q_low_w_cm2"] == 0.05 && b["q_high_w_cm2"] == 0.25 && b["q_target_w_cm2"] == 0.15 && b["q_medium_w_cm2"] == 0.30 && b["q_hard_w_cm2"] == 0.45
    gm = sqrt(Float64(v["band"]["q_low_w_cm2"]) * Float64(v["band"]["q_high_w_cm2"])); @assert abs(gm - 0.15) < 1e-15
    @assert sc["band"]["q_low_w_cm2"] == 0.05   # the original dict is not mutated
    println("stagei selftest ok: calibration view geometric mean $(gm) W/cm2 = q_target; hctrl=$HCTRL")
end

function main()
    mkpath(SF); phase = get(CLI, "phase", "")
    phase == "selftest" && return selftest()
    pls = String[]
    for pl in PLANETS
        try; sc = load_stagef(pl; path=CFGI); stagef_setup(sc); push!(pls, pl); catch err; println("INFEASIBLE $pl: ", sprint(showerror, err)); end
    end
    isempty(pls) && error("no feasible planet")
    phase == "cal" && phase_cal(pls); phase == "dev" && phase_dev(pls); phase == "held" && phase_held(pls)
end
main()
