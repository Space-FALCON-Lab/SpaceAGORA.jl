# Stage N2 driver: edge-aimed chance-constrained controller with a calibrated aim offset (stagen2_config.toml). One run_monte_carlo campaign per call (process route asserted); memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. stagen2_run.jl --phase=smoke|dev|held|repro|selftest --results-root=DIR [--workers=N] --host=NAME [--only=tag,tag] [--dry-run]
# Trees: calibrations and cells are READ from <root>/stagei; filter constants from <root>/stagen_flight/<planet>/kalman_params.csv (unchanged Stage N-flight values); the calibrated arms from <root>/stagen2/<planet>/aim_grid.csv
# (stagen2_analysis.py --phase calibrate) and the selected arm from selected.csv (--phase select). Output: <root>/stagen2/<planet>/runs/<code>_pert_<tag>/ (runs already on disk are skipped).
#   smoke  the first dev half-cell, first two grid arms (2 flights: one case plans route none, so two)    dev  the 12 dev half-cells x the grid arms    held  the selected arm on the 24 held-out cells
#   repro  the first two held-out cells of the selected arm into <planet>/repro (bitwise determinism check against the held runs)
# Rows: <root>/stagen2/summary_<phase>_<planet>_<host>.csv (use a DISTINCT --host stamp per call).
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "stagen2_common.jl"))
const CFGN = joinpath(@__DIR__, "stagen2_config.toml")

function parse_cli()
    kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)(?:=(.*))?$", a); m === nothing || (kv[m.captures[1]] = something(m.captures[2], "1")); end
    return kv
end
const CLI = parse_cli(); const META = stagef_meta(CFGN); const SP = META["stagen2"]
const ROOT = abspath(get(CLI, "results-root", ".")); const SFN = joinpath(ROOT, String(META["stagef"]["results_subdir"])); const SFI = joinpath(ROOT, "stagei"); const SFK = joinpath(ROOT, "stagen_flight")
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

"""Filter constants (unchanged Stage N-flight values) from <root>/stagen_flight; the calibrated arms (g_a, t_mult) from <root>/stagen2/<planet>/aim_grid.csv (stagen2_analysis.py --phase calibrate)."""
function kconst(pl)
    P = CSV.read(joinpath(SFK, pl, "kalman_params.csv"), DataFrame)[1, :]; return (; kr=Float64(P.r), kq=Float64(P.q), km0=Float64(P.m0), kv0=Float64(P.v0))
end
arm_of(pl, g, t) = merge((; eps=Float64(SP["epsilon"]), w_m=Float64(SP["w_cap_m"]), aim=Float64(g), t_mult=Float64(t), code=stagen2_arm_code(SP["epsilon"], t, g, SP["w_cap_m"])), kconst(pl))
grid_arms(pl) = (G = CSV.read(joinpath(SFN, pl, "aim_grid.csv"), DataFrame); [arm_of(pl, r.g_a, r.t_mult) for r in eachrow(G)])
selected_arm(pl) = (S = CSV.read(joinpath(SFN, pl, "selected.csv"), DataFrame)[1, :]; arm_of(pl, S.g_a, S.t_mult))

seedn(pl, r, a, outdir) = (; config_path=CFGN, planet=pl, tag=r.tag, inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch_idx), arm="pert", ctrl="cc", code=a.code, eps=a.eps, w_m=a.w_m, aim=a.aim,
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
campaign(seeds, mission_s) = run_gram_campaign(stagen2_run_task, seeds, workers(); mission_time_s=mission_s, worker_files=[joinpath(@__DIR__, "stagen2_common.jl")])

function build_seeds(phase, pls)
    seeds = NamedTuple[]
    for pl in pls
        runs = joinpath(SFN, pl, "runs"); todo(a, r, dir) = calibrated(pl, r.tag) && !isfile(joinpath(dir, "$(a.code)_pert_$(r.tag)", "passes.csv")) && push!(seeds, seedn(pl, r, a, dir))
        if phase == "smoke"
            for r in eachrow(cells_of(pl, "smoke")), a in first(grid_arms(pl), 2); todo(a, r, runs); end
        elseif phase == "dev"
            for r in eachrow(cells_of(pl, "dev")), a in grid_arms(pl); todo(a, r, runs); end
        elseif phase == "held"
            a = selected_arm(pl); for r in eachrow(cells_of(pl, "heldout")); todo(a, r, runs); end
        elseif phase == "repro"
            a = selected_arm(pl); for r in eachrow(cells_of(pl, "first2")); todo(a, r, joinpath(SFN, pl, "repro")); end
        end
    end
    return seeds
end

"""Exact checks without any flight: arm code, aim target geometry (target inside the band, monotone in g_a, equal to the edge + 10 m at g_a = 0, below the edge in Q for g_a > 0)."""
function selftest()
    @assert stagen2_arm_code(0.01, 1.4, 0.1, 4000.0) == "cc01_a010_t140_w04000" && stagen2_arm_code(0.01, 1.9024, 0.0, 4000.0) == "cc01_a000_t190_w04000"
    H = 6000.0; fit = exp_fit(1e-8, 120000.0, H, 400000.0, 170.0, "selftest"); ys = [0.1, 0.2, 0.0]
    mk(t, w, g) = (; f=(; fit, vrel_first=4000.0, q_low_w_m2=0.05e4, q_high_w_m2=0.25e4, meas=[(; y=y) for y in ys], min_band_m=200.0), kr=0.1, kq=0.01, km0=0.0, kv0=0.2, t_mult=t, w_m=w, aim_ln=g)
    e0 = stagen2_edges(mk(1.9, 4000.0, 0.0)); e1 = stagen2_edges(mk(1.9, 4000.0, 0.1)); e2 = stagen2_edges(mk(1.9, 4000.0, 0.3)); e3 = stagen2_edges(mk(1.9, 4000.0, 5.0))
    @assert abs(e0.target - (e0.lo + 10.0)) < 1e-6 && e0.lo < e1.target < e2.target < e0.hi && abs(e1.target - e0.lo - H * 0.1) < 20.0 && abs(e3.target - (e3.hi - 10.0)) < 1e-6
    @assert stagen2_edges(mk(1.4, 4000.0, 0.1)).lo < e1.lo    # a smaller multiplier moves the edge to lower altitude (higher Q)
    println("stagen2 selftest ok: target - lo = $(round(e0.target - e0.lo)) / $(round(e1.target - e0.lo)) / $(round(e2.target - e0.lo)) m for g_a = 0 / 0.1 / 0.3 (scale height $H m)")
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
