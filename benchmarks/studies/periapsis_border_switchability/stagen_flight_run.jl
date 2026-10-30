# Stage N-flight driver: closed-loop flight of the chance-constrained edge rule (stagen_flight_config.toml). One run_monte_carlo campaign per call (process route asserted); memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. stagen_flight_run.jl --phase=smoke|dev|held|repro|selftest --results-root=DIR [--workers=N] --host=NAME [--only=tag,tag] [--dry-run]
# Trees: calibrations and cells are READ from <root>/stagei/<planet>/{cells_base.csv,cal}; filter constants from <root>/stagen_flight/<planet>/kalman_params.csv (stagen_flight_analysis.py --phase params);
# the dev-selected band width from <root>/stagen_flight/<planet>/selected.csv (--phase select). Output: <root>/stagen_flight/<planet>/runs/<code>_pert_<tag>/ (runs already on disk are skipped).
#   smoke  the first dev half-cell, arms w = none and w = 2000 m (2 flights, 2 workers)
#   dev    the 12 dev half-cells x the band-width grid (the only tuned knob)
#   held   the selected arm on the 24 held-out cells
#   repro  the first two held-out cells of the selected arm into <planet>/repro (bitwise determinism check against the held runs)
# Rows: <root>/stagen_flight/summary_<phase>_<planet>_<host>.csv (use a DISTINCT --host stamp per call).
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "stagen_flight_common.jl"))
const CFGN = joinpath(@__DIR__, "stagen_flight_config.toml")

function parse_cli()
    kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)(?:=(.*))?$", a); m === nothing || (kv[m.captures[1]] = something(m.captures[2], "1")); end
    return kv
end
const CLI = parse_cli(); const META = stagef_meta(CFGN); const SP = META["stagen"]
const ROOT = abspath(get(CLI, "results-root", ".")); const SFN = joinpath(ROOT, String(META["stagef"]["results_subdir"])); const SFI = joinpath(ROOT, "stagei")
const HOST = get(CLI, "host", "local")
const PLANETS = haskey(CLI, "planets") ? String.(split(CLI["planets"], ",")) : ["mars"]
workers() = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(META["stagef"]["workers"])
hd(pl) = Float64(load_stagef(pl; path=CFGN)["grid"]["horizon_days"])

"""Cells of one planet from the Stage I tree (same rule as stagel_run.jl): "dev" = first dev cell of every (epoch, inclination) group, "heldout" = all, "first2" = the first two held-out cells, "smoke" = the first dev half-cell."""
function cells_of(pl, which)
    D = CSV.read(joinpath(SFI, pl, "cells_base.csv"), DataFrame)
    if which in ("dev", "smoke")
        D = D[D.split .== "dev", :]; seen = Set{Tuple{Int, Float64}}(); keep = Bool[]
        for r in eachrow(D); k = (Int(r.epoch_idx), Float64(r.inc)); push!(keep, !(k in seen)); push!(seen, k); end
        D = D[keep, :]; which == "smoke" && (D = first(D, 1))
    elseif which == "heldout"; D = D[D.split .== "heldout", :]
    elseif which == "first2"; D = first(D[D.split .== "heldout", :], 2)
    end
    haskey(CLI, "only") && (D = D[in.(D.tag, Ref(split(CLI["only"], ","))), :])
    return D
end
calibrated(pl, tag) = isfile(joinpath(SFI, pl, "cal", tag, "guess.csv"))

"""Filter constants and quantile multiplier of one planet (kalman_params.csv: r, q, m0, v0, t_eps_<100 eps>)."""
function kparams(pl)
    P = CSV.read(joinpath(SFN, pl, "kalman_params.csv"), DataFrame)[1, :]; eps = Float64(SP["epsilon"])
    return (; kr=Float64(P.r), kq=Float64(P.q), km0=Float64(P.m0), kv0=Float64(P.v0), t_mult=Float64(P[Symbol("t_emp_" * lpad(string(round(Int, 100eps)), 2, '0'))]))
end

arm_of(pl, w) = merge((; eps=Float64(SP["epsilon"]), w_m=Float64(w), code=stagen_arm_code(SP["epsilon"], w)), kparams(pl))
function selected_w(pl)
    S = CSV.read(joinpath(SFN, pl, "selected.csv"), DataFrame); return Float64(S.w_m[1])
end

seedn(pl, r, a, outdir) = (; config_path=CFGN, planet=pl, tag=r.tag, inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch_idx), arm="pert", ctrl="cc", code=a.code, eps=a.eps, w_m=a.w_m,
    run_id="$(a.code)_pert_$(r.tag)", out=outdir, cal_root=joinpath(SFI, pl, "cal"), horizon_days=hd(pl), kr=a.kr, kq=a.kq, km0=a.km0, kv0=a.kv0, t_mult=a.t_mult)

function finish(res, phase, pls)
    tag = phase * "_" * join(pls, "-") * "_" * HOST
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time, not a performance number)")
    foreach(f -> println("FAILED $(f.seed.run_id): $(first(sprint(showerror, f.error), 600))"), res.failed)
    rows = [merge(s.value, (; host=HOST)) for s in res.successful]
    write_rows_m(joinpath(SFN, "summary_$(tag).csv"), rows)
    open(joinpath(SFN, "mc_$(tag).toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_density_family = \"gram\"\nhost = \"$HOST\"")
    end
    res.route === :process || error("route was $(res.route), not :process")
    bad = count(r -> r.status != "ok", rows); println("rows=$(length(rows)) not ok=$bad")
    for r in rows; r.status != "ok" && println("  $(r.planet) $(r.run_id): $(r.error)"); end
    (isempty(res.failed) && bad == 0) || error("failed or not-ok cases present")
end
campaign(seeds, mission_s) = run_gram_campaign(stagen_run_task, seeds, workers(); mission_time_s=mission_s, worker_files=[joinpath(@__DIR__, "stagen_flight_common.jl")])

function build_seeds(phase, pls)
    seeds = NamedTuple[]
    for pl in pls
        runs = joinpath(SFN, pl, "runs"); todo(a, r, dir) = calibrated(pl, r.tag) && !isfile(joinpath(dir, "$(a.code)_pert_$(r.tag)", "passes.csv")) && push!(seeds, seedn(pl, r, a, dir))
        if phase == "smoke"
            for r in eachrow(cells_of(pl, "smoke")), w in (0.0, 2000.0); todo(arm_of(pl, w), r, runs); end
        elseif phase == "dev"
            for r in eachrow(cells_of(pl, "dev")), w in Float64.(SP["band_width_grid_m"]); todo(arm_of(pl, w), r, runs); end
        elseif phase == "held"
            a = arm_of(pl, selected_w(pl)); for r in eachrow(cells_of(pl, "heldout")); todo(a, r, runs); end
        elseif phase == "repro"
            a = arm_of(pl, selected_w(pl)); for r in eachrow(cells_of(pl, "first2")); todo(a, r, joinpath(SFN, pl, "repro")); end
        end
    end
    return seeds
end

"""Exact checks without any flight: arm codes, the online filter against hand-computed values of the same recursion, edge geometry (cap, narrow rule, monotone in the quantile multiplier)."""
function selftest()
    @assert stagen_arm_code(0.01, 0.0) == "cc01_wnone" && stagen_arm_code(0.01, 2000.0) == "cc01_w02000" && stagen_arm_code(0.001, 4000.0) == "cc00_w04000"
    # filter: r = 0.1, q = 0.01, prior (0, 0.2). No data: m = 0, s2 = 0.3. One y = 1.0: K = 0.2/0.3, m = 2/3, p = 0.2/3*... p1 = (1-K)*0.2 = 1/15; predictive p = 1/15 + 0.01, s2 = p + 0.1
    pr0 = stagen_predict(Float64[], 0.1, 0.01, 0.0, 0.2); @assert pr0.m == 0.0 && abs(pr0.s2 - 0.3) < 1e-15
    pr1 = stagen_predict([1.0], 0.1, 0.01, 0.0, 0.2); @assert abs(pr1.m - 2 / 3) < 1e-15 && abs(pr1.s2 - (1 / 15 + 0.01 + 0.1)) < 1e-15
    # second y = 0: p before update = 1/15 + 0.01 = 0.076666..., s2 = p + 0.1, K = p / s2, m = 2/3 (1 - K)
    p2 = 1 / 15 + 0.01; K2 = p2 / (p2 + 0.1); pr2 = stagen_predict([1.0, 0.0], 0.1, 0.01, 0.0, 0.2); @assert abs(pr2.m - (2 / 3) * (1 - K2)) < 1e-15 && abs(pr2.s2 - ((1 - K2) * p2 + 0.01 + 0.1)) < 1e-15
    H = 6000.0; fit = exp_fit(1e-8, 120000.0, H, 400000.0, 170.0, "selftest"); ys = [0.1, 0.2, 0.0]
    mk(t, w, meas=ys) = (; f=(; fit, vrel_first=4000.0, q_low_w_m2=0.05e4, q_high_w_m2=0.25e4, meas=[(; y=y) for y in meas], min_band_m=200.0), kr=0.1, kq=0.01, km0=0.0, kv0=0.2, t_mult=t, w_m=w)
    e0 = stagen_edges(mk(2.0, 0.0)); e1 = stagen_edges(mk(3.0, 0.0)); ec = stagen_edges(mk(2.0, 2000.0))
    @assert e1.lo > e0.lo "a larger t is more conservative: the thermal edge moves to higher altitude"
    println("stagen selftest ok: codes, filter recursion; lo(t=2) = $(round(e0.lo)), lo(t=3) = $(round(e1.lo)), band(none) = $(round(e0.hi - e0.lo)), band(cap 2000) = $(round(ec.hi - ec.lo)) m")
    @assert abs((ec.hi - ec.lo) - 2000.0) < 1e-6 && e0.hi > ec.hi
end

function main()
    mkpath(SFN); phase = get(CLI, "phase", "")
    phase == "selftest" && return selftest()
    phase in ("smoke", "dev", "held", "repro") || error("--phase must be smoke|dev|held|repro|selftest")
    pls = String[]
    for pl in PLANETS
        try; sc = load_stagef(pl; path=CFGN); stagef_setup(sc); push!(pls, pl); catch err; println("INFEASIBLE $pl: ", sprint(showerror, err)); end
    end
    isempty(pls) && error("no feasible planet")
    seeds = build_seeds(phase, pls); println("$phase: $(length(seeds)) cases over $(pls)")
    haskey(CLI, "dry-run") && (foreach(c -> println("  ", c, " x ", count(x -> x.code == c, seeds)), unique(getfield.(seeds, :code))); return println("dry run: nothing flown"))
    isempty(seeds) && error("no cases (missing kalman_params.csv / selected.csv / calibrations, or all runs already on disk)")
    finish(campaign(seeds, maximum(hd.(pls)) * 86400.0), phase, pls)
end
main()
