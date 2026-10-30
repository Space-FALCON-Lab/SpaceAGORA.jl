# Held baselines F-O and N2-I only (seeds 13 to 28, plan R4.5 / R7.3), through the frozen Stage O flight code. The frozen driver stageo_run.jl cannot select arms (its held phase flies N2-1 as well, 16 extra flights), so this
# new driver includes stageo_common.jl (not edited) and rebuilds the flight seeds exactly as stageo_run.jl mkseed does (same config, same arm codes, same run ids, same GRAM tags ODY_s<k>). Usage:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. planner_baselines_run.jl --results-root=DIR --workers=3 --host=NAME [--seeds=13:28] [--dry-run]
# The results root must contain stageo/mars/cal/odyssey (Stage O calibration), stageo/mars/calibrated_arms.csv and stagen_flight/mars/kalman_params.csv (as the Stage O driver). Process route asserted; a summary is never overwritten.
using TOML
const CLI = let kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)(?:=(.*))?$", a); m === nothing || (kv[m.captures[1]] = something(m.captures[2], "1")); end
    kv
end
const DRY = haskey(CLI, "dry-run")
const CFGO = joinpath(@__DIR__, "stageo_config.toml")
const META = TOML.parsefile(CFGO); const SO = META["stageo"]
include(joinpath(@__DIR__, "stageo_common.jl"))
const ROOT = abspath(get(CLI, "results-root", ".")); const SFS = joinpath(ROOT, String(META["stagef"]["results_subdir"])); const SFO = joinpath(SFS, "mars"); const SFK = joinpath(ROOT, "stagen_flight")
const HOST = get(CLI, "host", "local")
workers() = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(META["stagef"]["workers"])
horizon_days(key) = Float64(SO["arms"][key]["horizon_days"])
function kconst()
    P = CSV.read(joinpath(SFK, "mars", "kalman_params.csv"), DataFrame)[1, :]; return (; kr=Float64(P.r), kq=Float64(P.q), km0=Float64(P.m0), kv0=Float64(P.v0))
end
function calibrated_pair(name)
    D = CSV.read(joinpath(SFO, "calibrated_arms.csv"), DataFrame); i = findfirst(==(name), String.(D.arm)); i === nothing && error("arm $name not in calibrated_arms.csv"); return (Float64(D.g_a[i]), Float64(D.t_mult[i]))
end
function mkseed(key, k, dir)   # identical to stageo_run.jl mkseed (default arguments)
    A = SO["arms"][key]; fixed = A["ctrl"] == "fixed"; K = kconst()
    g, t = fixed ? (NaN, NaN) : get(A, "calibrated", false) ? calibrated_pair(A["name"]) : (Float64(A["g_a"]), Float64(A["t_mult"]))
    eps = fixed ? NaN : Float64(A["eps"]); w = fixed ? NaN : Float64(SO["n2"]["w_cap_m"]); margin = Float64(SO["n2"]["fixed_margin_m"])
    code = stageo_arm_code(key, A["ctrl"], margin, eps, t, g, w); tg = String(SO["seed_prefix"]) * string(k)
    return (; config_path=CFGO, band=String(A["band"]), tag=tg, inc=Float64(SO["inclination_deg"]), raan=Float64(SO["raan_deg"]), omega=Float64(SO["omega_deg"]), arm=key, ctrl=String(A["ctrl"]), code,
        eps, w_m=w, aim=g, t_mult=t, kr=K.kr, kq=K.kq, km0=K.km0, kv0=K.kv0, run_id="$(code)_pert_$(tg)", out=dir, cal_root=joinpath(SFO, "cal"),
        horizon_days=horizon_days(key), dt_atm_s=NaN, role="flight", name=String(A["name"]))
end
function cases()
    lo, hi = haskey(CLI, "seeds") ? parse.(Int, split(CLI["seeds"], ":")) : (13, 28)
    return [mkseed(key, k, joinpath(SFO, "runs")) for k in lo:hi, key in ("FO", "N2I")] |> vec
end
ondisk(c) = isfile(joinpath(c.out, c.run_id, "passes.csv"))
function main()
    cs = cases(); todo = filter(c -> !ondisk(c), cs)
    println("baselines: $(length(cs)) flights (F-O and N2-I), $(length(todo)) to fly; arms $(sort(unique(c.name for c in cs))) seeds $(minimum(parse(Int, split(c.tag, "_s")[2]) for c in cs))..$(maximum(parse(Int, split(c.tag, "_s")[2]) for c in cs))")
    DRY && (println("dry run: nothing flown; first run id $(cs[1].run_id), horizon $(cs[1].horizon_days) d"); return)
    isempty(todo) && error("no cases left")
    summary = joinpath(SFS, "summary_held_mars_$(HOST).csv"); isfile(summary) && error("$summary exists: use a distinct --host stamp per call")
    res = run_gram_campaign(stageo_run_task, todo, workers(); mission_time_s=maximum(c.horizon_days for c in todo) * 86400.0, worker_files=[joinpath(@__DIR__, "stageo_common.jl")])
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time, not a performance number)")
    foreach(f -> println("FAILED $(f.seed.run_id): $(first(sprint(showerror, f.error), 600))"), res.failed)
    rows = [merge(s.value, (; host=HOST)) for s in res.successful]; mkpath(SFS); write_rows_m(summary, rows)
    open(joinpath(SFS, "mc_held_mars_$(HOST).toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_density_family = \"gram\"\nhost = \"$HOST\"")
    end
    res.route === :process || error("route was $(res.route), not :process")
    bad = count(r -> r.status != "ok", rows); println("rows=$(length(rows)) not ok=$bad")
    for r in rows; r.status == "ok" && println("  $(r.run_id): passes=$(r.passes) t_end_days=$(round(r.t_end_days; digits=2)) stop_reason=$(r.stop_reason) dv=$(round(r.total_dv_mps; digits=4)) burns=$(r.burns) core_epoch=$(r.gram_core_epoch)"); end
    (isempty(res.failed) && bad == 0) || error("failed or not-ok cases present")
end
main()
