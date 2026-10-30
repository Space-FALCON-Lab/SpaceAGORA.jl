# Stage H driver: asymmetric band with a partial low-heat relaxation e_lo_frac * |e_lo| (see stageh_common.jl). Reads stagef20_config.toml and <root>/stagef20/<planet>/{cells_base.csv,cal,adapt_params.csv}.
# One run_monte_carlo campaign per call (process route asserted); memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. stageh_run.jl --phase=dev|held|selftest --frac=0.25 --results-root=DIR [--workers=N] [--planets=mars] [--host=NAME] [--horizon-days=D] [--only=tag,tag]
# dev = development cells (split "dev", used for choosing e_lo_frac), held = held-out cells (evaluated once at the chosen frac).
# Output: <root>/stagef20/<planet>/runs/adaptive_h<FFF>_pert_<tag>/ (FFF = 100*frac, three digits), rows in <root>/stagef20/summary_h<phase>FFF_<host>.csv. Plan: tasks/2026-10-15_periapsis-region-border-targeting/stageh_plan.md.
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "stageh_common.jl"))
const CFG20 = joinpath(@__DIR__, "stagef20_config.toml")
function parse_cli()
    kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2]); end
    return kv
end
const CLI = parse_cli(); const META = stagef_meta(CFG20)
const ROOT = abspath(get(CLI, "results-root", ".")); const SF = joinpath(ROOT, String(META["stagef"]["results_subdir"]))
const HOST = get(CLI, "host", "local")
const PLANETS = haskey(CLI, "planets") ? String.(split(CLI["planets"], ",")) : String.(META["stagef"]["planets"])
workers() = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(META["stagef"]["workers"])
hd(pl) = haskey(CLI, "horizon-days") ? parse(Float64, CLI["horizon-days"]) : Float64(load_stagef(pl; path=CFG20)["grid"]["horizon_days"])

function cells_of(pl)
    D = CSV.read(joinpath(get(CLI, "cells-dir", SF), pl, "cells_base.csv"), DataFrame)
    haskey(CLI, "only") && (D = D[in.(D.tag, Ref(split(CLI["only"], ","))), :])
    return D
end

function phase_h(pls, split, frac)
    code = lpad(string(round(Int, 100 * frac)), 3, '0'); ctrl = "adaptive_h" * code; phname = "h" * (split == "dev" ? "dev" : "held") * code
    seeds = NamedTuple[]
    for pl in pls
        P = CSV.read(joinpath(SF, pl, "adapt_params.csv"), DataFrame)[1, :]   # development-cell predictor, e_hi and e_lo (Stage F at 20x20, unchanged)
        g = frac * abs(Float64(P.e_lo))                                       # relaxation magnitude, ln-density units
        for r in eachrow(cells_of(pl))
            r.split == split && isfile(joinpath(SF, pl, "cal", r.tag, "guess.csv")) || continue
            push!(seeds, (; config_path=CFG20, planet=pl, tag=r.tag, inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch_idx), arm="pert", ctrl=ctrl, run_id="$(ctrl)_pert_$(r.tag)",
                out=joinpath(SF, pl, "runs"), cal_root=joinpath(SF, pl, "cal"), horizon_days=hd(pl), adapt_pred=String(P.predictor), adapt_e_hi=Float64(P.e_hi), adapt_e_lo=g))
        end
    end
    println("$phname: $(length(seeds)) cases over $(pls), frac=$frac")
    res = run_gram_campaign(stageh_run_task, seeds, workers(); mission_time_s=maximum(hd.(pls)) * 86400.0, worker_files=[joinpath(@__DIR__, "stageh_common.jl")])
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time, not a performance number)")
    foreach(f -> println("FAILED $(f.seed.run_id): $(first(sprint(showerror, f.error), 600))"), res.failed)
    rows = [merge(s.value, (; host=HOST)) for s in res.successful]
    write_rows_m(joinpath(SF, "summary_$(phname)_$(HOST).csv"), rows)
    open(joinpath(SF, "mc_$(phname)_$(HOST).toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_density_family = \"gram\"\nhost = \"$HOST\"")
    end
    res.route === :process || error("route was $(res.route), not :process")
    bad = count(r -> r.status != "ok", rows); println("rows=$(length(rows)) not ok=$bad")
    for r in rows; r.status != "ok" && println("  $(r.planet) $(r.run_id): $(r.error)"); end
end

"""Exact-case check of stageh_edges on an exponential fit: analytic edges, relaxed low edge (g > 0 moves hi UP by H*g), g = 0 equals the Stage G rule, narrow-band anchoring."""
function selftest()
    H = 6000.0; rho0 = 1e-8; h0 = 120000.0; fit = exp_fit(rho0, h0, H, 400000.0, 170.0, "selftest")
    mk(e_hi, g, mb) = (; fit, vrel_first=4000.0, q_low_w_m2=0.14e4, q_high_w_m2=0.32e4, meas=[(; y=0.2, lon_rad=0.0)], pred_name="P1", e_hi, e_lo=g, min_band_m=mb)
    V3 = 4000.0^3; alt(q, off) = h0 - H * log(2q / V3 * exp(-off) / rho0)
    e0 = stageh_edges(mk(0.3, 0.0, 200.0))
    @assert abs(e0.lo - alt(0.32e4, 0.5)) < 1e-6 * H && abs(e0.hi - alt(0.14e4, 0.2)) < 1e-6 * H && !e0.collapsed && e0.margin == 0.0   # g = 0: Stage G edges
    e1 = stageh_edges(mk(0.3, 0.1, 200.0))
    @assert abs(e1.hi - alt(0.14e4, 0.3)) < 1e-6 * H && abs((e1.hi - e0.hi) - H * 0.1) < 1e-6 * H && e1.lo == e0.lo   # hi moves up by H*g, thermal edge untouched
    e2 = stageh_edges(mk(2.0, 0.0, 200.0))
    @assert e2.collapsed && abs(e2.lo - alt(0.32e4, 2.2)) < 1e-6 * H && abs((e2.hi - e2.lo) - 200.0) < 1e-9
    println("stageh selftest ok: g=0 band $(round(e0.hi - e0.lo)) m, g=0.1 band $(round(e1.hi - e1.lo)) m")
end

function main()
    mkpath(SF); phase = get(CLI, "phase", "")
    phase == "selftest" && return selftest()
    pls = String[]
    for pl in PLANETS
        try; sc = load_stagef(pl; path=CFG20); stagef_setup(sc); push!(pls, pl); catch err; println("INFEASIBLE $pl: ", sprint(showerror, err)); end
    end
    isempty(pls) && error("no feasible planet")
    frac = parse(Float64, get(CLI, "frac", "NaN")); (0.0 < frac <= 1.0 || phase == "selftest") || error("--frac in (0,1] required")
    phase == "dev" && phase_h(pls, "dev", frac)
    phase == "held" && phase_h(pls, "heldout", frac)
end
main()
