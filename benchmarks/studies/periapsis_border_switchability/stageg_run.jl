# Stage G driver: asymmetric adaptive band on the held-out cells of the Stage F at 20x20 tree (reads stagef20_config.toml, <root>/stagef20/<planet>/{cells_base.csv,cal,adapt_params.csv}).
# One run_monte_carlo campaign (process route asserted); memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. stageg_run.jl --phase=f2g|selftest --results-root=DIR [--workers=N] [--planets=venus,mars,earth] [--host=NAME]
#        [--horizon-days=D] [--only=tag,tag] [--cells-dir=DIR]
# Output: <root>/stagef20/<planet>/runs/adaptive_asym_pert_<tag>/, rows in <root>/stagef20/summary_f2g_<host>.csv. Plan: tasks/2026-10-15_periapsis-region-border-targeting/stageg_plan.md.
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "stageg_common.jl"))
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

function phase_f2g(pls)
    seeds = NamedTuple[]
    for pl in pls
        P = CSV.read(joinpath(SF, pl, "adapt_params.csv"), DataFrame)[1, :]   # development-cell predictor and e_hi (Stage F at 20x20, unchanged); e_lo is not used
        for r in eachrow(cells_of(pl))
            r.split == "heldout" && isfile(joinpath(SF, pl, "cal", r.tag, "guess.csv")) || continue
            push!(seeds, (; config_path=CFG20, planet=pl, tag=r.tag, inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch_idx), arm="pert", ctrl="adaptive_asym", run_id="adaptive_asym_pert_$(r.tag)",
                out=joinpath(SF, pl, "runs"), cal_root=joinpath(SF, pl, "cal"), horizon_days=hd(pl), adapt_pred=String(P.predictor), adapt_e_hi=Float64(P.e_hi), adapt_e_lo=0.0))
        end
    end
    println("f2g: $(length(seeds)) cases over $(pls)")
    res = run_gram_campaign(stageg_run_task, seeds, workers(); mission_time_s=maximum(hd.(pls)) * 86400.0, worker_files=[joinpath(@__DIR__, "stageg_common.jl")])
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time, not a performance number)")
    foreach(f -> println("FAILED $(f.seed.run_id): $(first(sprint(showerror, f.error), 600))"), res.failed)
    rows = [merge(s.value, (; host=HOST)) for s in res.successful]
    write_rows_m(joinpath(SF, "summary_f2g_$(HOST).csv"), rows)
    open(joinpath(SF, "mc_f2g_$(HOST).toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_density_family = \"gram\"\nhost = \"$HOST\"")
    end
    bad = count(r -> r.status != "ok", rows); println("rows=$(length(rows)) not ok=$bad")
    for r in rows; r.status != "ok" && println("  $(r.planet) $(r.run_id): $(r.error)"); end
end

"""Exact-case check of stageg_edges on an exponential fit (rho = rho0 exp(-(h-h0)/H)): analytic edges, no-buffer low edge, narrow-band anchoring of the thermal edge."""
function selftest()
    H = 6000.0; rho0 = 1e-8; h0 = 120000.0; fit = exp_fit(rho0, h0, H, 400000.0, 170.0, "selftest")
    mk(e_hi, mb) = (; fit, vrel_first=4000.0, q_low_w_m2=0.14e4, q_high_w_m2=0.32e4, meas=[(; y=0.2, lon_rad=0.0)], pred_name="P1", e_hi, min_band_m=mb)
    V3 = 4000.0^3; alt(q, off) = h0 - H * log(2q / V3 * exp(-off) / rho0)   # density 2Q/V^3 e^-(yhat+e) at altitude alt
    e = stageg_edges(mk(0.3, 200.0))
    @assert abs(e.lo - alt(0.32e4, 0.2 + 0.3)) < 1e-6 * H && abs(e.hi - alt(0.14e4, 0.2)) < 1e-6 * H && !e.collapsed && e.margin == 0.0
    e = stageg_edges(mk(2.0, 200.0))   # buffer wider than the band: thermal edge kept, low edge pushed to lo + min band
    @assert e.collapsed && abs(e.lo - alt(0.32e4, 2.2)) < 1e-6 * H && abs((e.hi - e.lo) - 200.0) < 1e-9
    println("stageg selftest ok: lo=$(round(e.lo)) hi=$(round(e.hi))")
end

function main()
    mkpath(SF); phase = get(CLI, "phase", "")
    phase == "selftest" && return selftest()
    pls = String[]
    for pl in PLANETS
        try; sc = load_stagef(pl; path=CFG20); stagef_setup(sc); push!(pls, pl); catch err; println("INFEASIBLE $pl: ", sprint(showerror, err)); end
    end
    isempty(pls) && error("no feasible planet")
    phase == "f2g" && phase_f2g(pls)
end
main()
