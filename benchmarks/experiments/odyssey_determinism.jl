using SpaceAGORA, DataFrames
const TV = SpaceAGORA.TelemetryVerification
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
cd(ROOT)
scen = TV._load_scenarios_from_manifest(joinpath(ROOT, "test", "telemetry_benchmark_manifest.toml"))
cfg = only(filter(s -> s.name == "odyssey", scen))
orbits = parse(Int, get(ENV, "ODY_ORBITS", "8"))
quick = get(ENV, "ODY_QUICK", "1") == "1"
reps = parse(Int, get(ENV, "ODY_REPS", "2"))
function once()
    args = TV._make_orbit_args(cfg, orbits; cd_scale=1.0, cr_override=cfg.srp_cr)
    args = TV._with_study_settings(args; quick=quick)
    t0 = time()
    r = TV._run_simulation_dataframe(args, cfg.name, cfg.atmosphere_truth, quick ? :quick : :full)
    return r.results_df, time() - t0
end
dfs = Any[]
using CSV
mkpath(joinpath(ROOT, "output", "experiment"))
for i in 1:reps
    df, dt = once()
    push!(dfs, df)
    CSV.write(joinpath(ROOT, "output", "experiment", "odyssey_orbits$(orbits)_rep$(i).csv"), df)
    println("rep$i: rows=", nrow(df), " wall=", round(dt; digits=1), "s")
end
numcols = [c for c in names(dfs[1]) if eltype(dfs[1][!, c]) <: Real]
println("nthreads_default=", Threads.nthreads(:default), " nthreads_interactive=", Threads.nthreads(:interactive), " cpu_threads=", Sys.CPU_THREADS)
println("threads=", Threads.nthreads(), " orbits=", orbits, " quick=", quick, " reps=", reps, " numeric_cols=", length(numcols))
for (i, df) in enumerate(dfs)
    fp = hash(reduce(vcat, [Float64.(df[!, c]) for c in numcols]))
    println("rep$i fingerprint=", string(fp, base=16), " rows=", nrow(df), " last=", [round(df[end, c]; sigdigits=12) for c in numcols[1:min(4, end)]])
end
for i in 2:reps
    a, b = dfs[1], dfs[i]
    same = nrow(a) == nrow(b) && all(c -> isequal(a[!, c], b[!, c]), numcols)
    println("rep$i vs rep1: bitwise_identical=", same)
    if !same && nrow(a) == nrow(b)
        worst = sort([(maximum(abs.(Float64.(a[!, c]) .- Float64.(b[!, c]))), c) for c in numcols]; rev=true)[1:min(5, end)]
        println("   largest column differences: ", worst)
    end
end
