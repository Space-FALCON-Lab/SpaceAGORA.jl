# Driver of the GRAM-on half-lever twin (I1): one run_monte_carlo campaign (process route asserted), memory capped:
#   systemd-run --user --scope -p MemoryMax=24G julia --project=. benchmarks/studies/periapsis_border_switchability/halflever_twin_run.jl --out=DIR --workers=N --host=NAME \
#       [--seeds=1,2,3,4] [--variants=none,p25,m25,p14,f136] [--npass=8] [--dry-run]
# Variants (dv m/s): none 0, p25 +0.25, m25 -0.25, p14 +0.14, f136 +0.136 at apoapsis 2, late25 +0.25 at apoapsis 100 (use --npass=106 and the no-burn twin of the same npass). Runs already on disk (passes.csv present) are skipped. A one-case campaign plans route none: the driver refuses fewer than 2 cases.
using TOML
const TCLI = let kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)(?:=(.*))?$", a); m === nothing || (kv[m.captures[1]] = something(m.captures[2], "1")); end
    kv
end
const TDRY = haskey(TCLI, "dry-run")
const TCFG = joinpath(@__DIR__, "stageo_config.toml")
const TOUT = abspath(get(TCLI, "out", joinpath(@__DIR__, "..", "..", "..", "output", "halflever_twin")))
const TCAL = get(TCLI, "cal-root", joinpath(get(ENV, "PLANNER_TASK_ROOT", "."), "stageo_data/phase_cal_smoke_gate/cal"))
const TVARIANTS = Dict("none" => (0.0, 1), "p25" => (0.25, 1), "m25" => (-0.25, 1), "p14" => (0.14, 1), "f136" => (0.136, 2), "late25" => (0.25, 100))   # name => (dv m/s, apoapsis index of the burn); f136 = the F-O keeper's first burn (0.136 m/s at apoapsis 2, SOURCED decisions.csv)
include(joinpath(@__DIR__, "stageo_common.jl"))
include(joinpath(@__DIR__, "halflever_twin_common.jl"))

function twin_cases()
    seeds = parse.(Int, split(get(TCLI, "seeds", "1,2,3,4"), ",")); vars = split(get(TCLI, "variants", "none,p25,m25,p14,f136"), ",")
    return [(; config_path=TCFG, tag="ODY_s$k", name="s$(k)_$(v)", dv=TVARIANTS[v][1], burn_j=TVARIANTS[v][2], npass=parse(Int, get(TCLI, "npass", "8")), horizon_days=12.0 + 0.5 * parse(Int, get(TCLI, "npass", "8")), out=TOUT, cal_root=TCAL, gate_altitude_m=1.0e6) for k in seeds for v in vars]
end

function main()
    cs = filter(c -> !isfile(joinpath(c.out, c.name, "passes.csv")), twin_cases())
    println("cases to fly: $(length(cs)) -> $(join([c.name for c in cs], ", "))")
    TDRY && return
    length(cs) >= 2 || error("a one-case campaign plans route none; pass at least 2 cases")
    mkpath(TOUT); host = get(TCLI, "host", "local"); w = parse(Int, get(TCLI, "workers", "2"))
    res = run_gram_campaign(twin_task, cs, w; mission_time_s=8.0e5, worker_files=[joinpath(@__DIR__, "stageo_common.jl"), joinpath(@__DIR__, "halflever_twin_common.jl")])
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time, not a performance number)")
    foreach(f -> println("FAILED $(f.seed.name): $(first(sprint(showerror, f.error), 600))"), res.failed)
    rows = [s.value for s in res.successful]; write_rows_m(joinpath(TOUT, "summary_$(host).csv"), rows)
    for r in rows; println("  $(r.name): status=$(r.status) passes=$(r.passes) steps=$(r.n_steps) t_end_s=$(round(r.t_end_s; digits=1)) dv_used=$(round(r.dv_used_mps; digits=5)) run_wall_s=$(round(r.run_wall_s; digits=1)) $(r.error)"); end
    res.route === :process || error("route was $(res.route), not :process")
    (isempty(res.failed) && all(r -> r.status == "ok", rows)) || error("failed or not-ok cases present")
end
main()
