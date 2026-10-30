# Stage N-flight cross-check: re-fly ONE held-out cell of the Stage I fixed 1000 m, Stage G (family s = 0, p = 0.95) and Stage F (family s = -1, p = 0.95) arms with the FROZEN Stage L case function on the current
# code (through the compatibility shim), into <root>/stagen_flight/mars/xcheck, for the bitwise comparison with the stored Stage I tables (stagen_flight_analysis.py --phase repro). Not a verdict criterion.
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. stagen_flight_xcheck.jl --results-root=DIR --workers=3 --host=NAME [--dry-run]
include(joinpath(@__DIR__, "stagen_flight_shim.jl"))
include(joinpath(@__DIR__, "stagel_common.jl"))
const CFGL = joinpath(@__DIR__, "stagel_config.toml")
kv = Dict{String, String}(); for a in ARGS; m = match(r"^--([^=]+)(?:=(.*))?$", a); m === nothing || (kv[m.captures[1]] = something(m.captures[2], "1")); end
ROOT = abspath(kv["results-root"]); SFI = joinpath(ROOT, "stagei"); SFL = joinpath(ROOT, "stagel"); SFN = joinpath(ROOT, "stagen_flight"); HOST = get(kv, "host", "xcheck"); PL = "mars"
META = stagef_meta(CFGL); P = CSV.read(joinpath(SFL, PL, "family_params.csv"), DataFrame)[1, :]; hd = Float64(load_stagef(PL; path=CFGL)["grid"]["horizon_days"])
cell = first(CSV.read(joinpath(SFI, PL, "cells_base.csv"), DataFrame) |> d -> d[d.split .== "heldout", :], 1)[1, :]
arm(code, kind, m, p, s) = (; kind, margin_m=m, p_hi=p, s_frac=s, pred=kind == "fixedm" ? "" : String(P.predictor), e_hi=kind == "fixedm" ? 0.0 : Float64(P.e_hi_095), e_lo=kind == "fixedm" ? 0.0 : s * abs(Float64(P.e_lo)), code)
arms = [arm("fm01000", "fixedm", 1000.0, NaN, NaN), arm("fam_p095_sp000", "adaptive_l", NaN, 0.95, 0.0), arm("fam_p095_sm100", "adaptive_l", NaN, 0.95, -1.0)]
@assert all(a.code == stagel_arm_code(a.kind, a.margin_m, a.p_hi, a.s_frac) for a in arms)
outdir = joinpath(SFN, PL, "xcheck")
seeds = [(; config_path=CFGL, planet=PL, tag=cell.tag, inc=cell.inc, raan=cell.raan, omega=cell.omega, epoch_idx=Int(cell.epoch_idx), arm="pert", kind=a.kind, code=a.code, ctrl=a.kind, run_id="$(a.code)_pert_$(cell.tag)",
    out=outdir, cal_root=joinpath(SFI, PL, "cal"), horizon_days=hd, margin_m=a.margin_m, p_hi=a.p_hi, s_frac=a.s_frac, adapt_pred=a.pred, adapt_e_hi=a.e_hi, adapt_e_lo=a.e_lo) for a in arms]
println("xcheck: $(length(seeds)) cases (cell $(cell.tag))")
haskey(kv, "dry-run") && (foreach(s -> println("  ", s.run_id), seeds); println("dry run: nothing flown"); exit(0))
res = run_gram_campaign(stagel_run_task, seeds, parse(Int, get(kv, "workers", "3")); mission_time_s=hd * 86400.0, worker_files=[joinpath(@__DIR__, "stagen_flight_shim.jl"), joinpath(@__DIR__, "stagel_common.jl")])
println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time, not a performance number)")
foreach(f -> println("FAILED $(f.seed.run_id): $(first(sprint(showerror, f.error), 600))"), res.failed)
rows = [merge(s.value, (; host=HOST)) for s in res.successful]; mkpath(SFN); write_rows_m(joinpath(SFN, "summary_xcheck_$(PL)_$(HOST).csv"), rows)
open(joinpath(SFN, "mc_xcheck_$(PL)_$(HOST).toml"), "w") do io
    println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_density_family = \"gram\"\nhost = \"$HOST\"")
end
res.route === :process || error("route was $(res.route), not :process")
bad = count(r -> r.status != "ok", rows); println("rows=$(length(rows)) not ok=$bad"); (isempty(res.failed) && bad == 0) || error("failed or not-ok cases present")
