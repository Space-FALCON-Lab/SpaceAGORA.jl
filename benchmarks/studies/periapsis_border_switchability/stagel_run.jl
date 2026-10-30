# Stage L driver: the dV / campaign-time / violation frontier (stagel_config.toml). One run_monte_carlo campaign per call (process route asserted); memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. stagel_run.jl --phase=repro|dev|held|heldall|selftest --results-root=DIR [--workers=N] [--planets=mars] --host=NAME [--only=tag,tag] [--dry-run]
# Trees: calibrations, cells and the Stage I development refit are READ from <root>/stagei/<planet>/{cells_base.csv,cal}; family parameters from <root>/stagel/<planet>/family_params.csv
# (stagel_analysis.py --phase=params); selected held-out members from <root>/stagel/<planet>/selected.csv (stagel_analysis.py --phase=select). Output: <root>/stagel/<planet>/runs/<code>_pert_<tag>/.
#   repro    one held-out cell per planet, arms fm01000, fam_p095_sp000, fam_p095_sm100 -> <planet>/repro/ (stagel_analysis.py --phase=repro compares them with the Stage I fixed / G / F runs)
#   dev      all arms of the grid on the development half-cells (first dev cell of each epoch x inclination group)
#   held     the dev-selected members (selected.csv) on the 24 held-out cells
#   heldall  OPTIONAL, descriptive: the whole grid on the held-out cells (arms already on disk are skipped)
# Rows: <root>/stagel/summary_<phase>_<planets>_<host>.csv (use a DISTINCT --host stamp per call). Plan: tasks/2026-10-15_periapsis-region-border-targeting/stagel_plan.md.
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "stagel_common.jl"))
const CFGL = joinpath(@__DIR__, "stagel_config.toml")

function parse_cli()
    kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)(?:=(.*))?$", a); m === nothing || (kv[m.captures[1]] = something(m.captures[2], "1")); end
    return kv
end
const CLI = parse_cli(); const META = stagef_meta(CFGL); const SP = META["stagel"]
const ROOT = abspath(get(CLI, "results-root", ".")); const SFL = joinpath(ROOT, String(META["stagef"]["results_subdir"])); const SFI = joinpath(ROOT, "stagei")
const HOST = get(CLI, "host", "local")
const PLANETS = haskey(CLI, "planets") ? String.(split(CLI["planets"], ",")) : String.(META["stagef"]["planets"])
workers() = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(META["stagef"]["workers"])
hd(pl) = Float64(load_stagef(pl; path=CFGL)["grid"]["horizon_days"])

"""Cells of one planet from the Stage I tree. which: "dev" (the dev half-cells: first dev cell of every (epoch, inclination) group), "heldout" (all), "repro" (the first held-out cell)."""
function cells_of(pl, which)
    D = CSV.read(joinpath(SFI, pl, "cells_base.csv"), DataFrame)
    if which == "dev"
        D = D[D.split .== "dev", :]; seen = Set{Tuple{Int, Float64}}(); keep = Bool[]
        for r in eachrow(D); k = (Int(r.epoch_idx), Float64(r.inc)); push!(keep, !(k in seen)); push!(seen, k); end
        D = D[keep, :]
    elseif which == "heldout"; D = D[D.split .== "heldout", :]
    elseif which == "repro"; D = first(D[D.split .== "heldout", :], 1)
    end
    haskey(CLI, "only") && (D = D[in.(D.tag, Ref(split(CLI["only"], ","))), :])
    return D
end
calibrated(pl, tag) = isfile(joinpath(SFI, pl, "cal", tag, "guess.csv"))

"""Full grid of arms of one planet: fixed margins and the family p x s (e_hi from family_params.csv, g = s * |e_lo|). NamedTuple per arm."""
function grid_arms(pl)
    P = CSV.read(joinpath(SFL, pl, "family_params.csv"), DataFrame)[1, :]; arms = NamedTuple[]
    for m in Float64.(SP["fixed_margins_m"]); push!(arms, (; kind="fixedm", margin_m=m, p_hi=NaN, s_frac=NaN, pred="", e_hi=0.0, e_lo=0.0, code=stagel_arm_code("fixedm", m, NaN, NaN))); end
    for p in Float64.(SP["p_hi"]), s in Float64.(SP["s_frac"])
        e_hi = Float64(P[Symbol("e_hi_" * lpad(string(round(Int, 100p)), 3, '0'))])
        push!(arms, (; kind="adaptive_l", margin_m=NaN, p_hi=p, s_frac=s, pred=String(P.predictor), e_hi=e_hi, e_lo=s * abs(Float64(P.e_lo)), code=stagel_arm_code("adaptive_l", NaN, p, s)))
    end
    return arms
end

"""Arms listed in <root>/stagel/<planet>/selected.csv (columns code, kind, margin_m, p_hi, s_frac), resolved against the grid."""
function selected_arms(pl)
    S = CSV.read(joinpath(SFL, pl, "selected.csv"), DataFrame); G = grid_arms(pl); out = NamedTuple[]
    for c in S.code; i = findfirst(a -> a.code == c, G); i === nothing && error("selected arm $c is not in the grid"); push!(out, G[i]); end
    return out
end

seedl(pl, r, a, outdir) = (; config_path=CFGL, planet=pl, tag=r.tag, inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch_idx), arm="pert", kind=a.kind, code=a.code, ctrl=a.kind,
    run_id="$(a.code)_pert_$(r.tag)", out=outdir, cal_root=joinpath(SFI, pl, "cal"), horizon_days=hd(pl), margin_m=a.margin_m, p_hi=a.p_hi, s_frac=a.s_frac,
    adapt_pred=a.pred, adapt_e_hi=a.e_hi, adapt_e_lo=a.e_lo)

function finish(res, phase, pls)
    tag = phase * "_" * join(pls, "-") * "_" * HOST
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time, not a performance number)")
    foreach(f -> println("FAILED $(f.seed.run_id): $(first(sprint(showerror, f.error), 600))"), res.failed)
    rows = [merge(s.value, (; host=HOST)) for s in res.successful]
    write_rows_m(joinpath(SFL, "summary_$(tag).csv"), rows)
    open(joinpath(SFL, "mc_$(tag).toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_density_family = \"gram\"\nhost = \"$HOST\"")
    end
    res.route === :process || error("route was $(res.route), not :process")
    bad = count(r -> r.status != "ok", rows); println("rows=$(length(rows)) not ok=$bad")
    for r in rows; r.status != "ok" && println("  $(r.planet) $(r.run_id): $(r.error)"); end
end
campaign(seeds, mission_s) = run_gram_campaign(stagel_run_task, seeds, workers(); mission_time_s=mission_s, worker_files=[joinpath(@__DIR__, "stagel_common.jl")])

function build_seeds(phase, pls)
    seeds = NamedTuple[]
    for pl in pls
        runs = joinpath(SFL, pl, "runs")
        if phase == "repro"
            arms = filter(a -> a.code in ("fm01000", "fam_p095_sp000", "fam_p095_sm100"), grid_arms(pl))
            for r in eachrow(cells_of(pl, "repro")), a in arms; calibrated(pl, r.tag) && push!(seeds, seedl(pl, r, a, joinpath(SFL, pl, "repro"))); end
        elseif phase == "dev"
            for r in eachrow(cells_of(pl, "dev")), a in grid_arms(pl); calibrated(pl, r.tag) && push!(seeds, seedl(pl, r, a, runs)); end
        elseif phase == "held"
            for r in eachrow(cells_of(pl, "heldout")), a in selected_arms(pl); calibrated(pl, r.tag) && push!(seeds, seedl(pl, r, a, runs)); end
        elseif phase == "heldall"
            for r in eachrow(cells_of(pl, "heldout")), a in grid_arms(pl)
                calibrated(pl, r.tag) && !isfile(joinpath(runs, "$(a.code)_pert_$(r.tag)", "passes.csv")) && push!(seeds, seedl(pl, r, a, runs))
            end
        end
    end
    return seeds
end

"""Exact checks without any flight: arm codes, signed edge offsets of the family, dev half-cell rule on a synthetic table."""
function selftest()
    @assert stagel_arm_code("fixedm", 500.0, NaN, NaN) == "fm00500" && stagel_arm_code("fixedm", 4000.0, NaN, NaN) == "fm04000"
    @assert stagel_arm_code("adaptive_l", NaN, 0.8, -1.0) == "fam_p080_sm100" && stagel_arm_code("adaptive_l", NaN, 0.99, 0.5) == "fam_p099_sp050" && stagel_arm_code("adaptive_l", NaN, 0.95, 0.0) == "fam_p095_sp000"
    # signed g: s = -1 reproduces Stage F's e_lo (narrow band), s = +0.5 widens, s = 0 is Stage G, all through the unchanged stageh_edges
    H = 6000.0; rho0 = 1e-8; h0 = 120000.0; fit = exp_fit(rho0, h0, H, 400000.0, 170.0, "selftest")
    mk(e_hi, g) = (; fit, vrel_first=4000.0, q_low_w_m2=0.05e4, q_high_w_m2=0.25e4, meas=[(; y=0.1, lon_rad=0.0)], pred_name="P1", e_hi, e_lo=g, min_band_m=200.0)
    elo = -0.58; e0 = stageh_edges(mk(0.5, 0.0)); en = stageh_edges(mk(0.5, 1.0 * elo)); ew = stageh_edges(mk(0.5, 0.5 * abs(elo)))
    @assert en.hi < e0.hi < ew.hi && abs((ew.hi - e0.hi) - H * 0.5 * abs(elo)) < 1e-6 * H && abs((e0.hi - en.hi) - H * abs(elo)) < 1e-6 * H && en.lo == e0.lo == ew.lo
    println("stagel selftest ok: codes, signed edges (band widths narrow/G/wide = $(round(en.hi - en.lo)) / $(round(e0.hi - e0.lo)) / $(round(ew.hi - ew.lo)) m)")
end

function main()
    mkpath(SFL); phase = get(CLI, "phase", "")
    phase == "selftest" && return selftest()
    phase in ("repro", "dev", "held", "heldall") || error("--phase must be repro|dev|held|heldall|selftest")
    pls = String[]
    for pl in PLANETS
        try; sc = load_stagef(pl; path=CFGL); stagef_setup(sc); push!(pls, pl); catch err; println("INFEASIBLE $pl: ", sprint(showerror, err)); end
    end
    isempty(pls) && error("no feasible planet")
    seeds = build_seeds(phase, pls); println("$phase: $(length(seeds)) cases over $(pls)")
    haskey(CLI, "dry-run") && (foreach(c -> println("  ", c, " x ", count(x -> x.code == c, seeds)), unique(getfield.(seeds, :code))); return println("dry run: nothing flown"))
    isempty(seeds) && error("no cases (missing family_params.csv / selected.csv / calibrations?)")
    finish(campaign(seeds, maximum(hd.(pls)) * 86400.0), phase, pls)
end
main()
