using Dates
using Statistics
using TOML

length(ARGS) == 5 ||
    error("Usage: benchmark_scheduler_runtime.jl prototype|4_v1|8_v5|9_full helpers orbits schedule inclination_deg")
const BENCH_CODE = popfirst!(ARGS)
const RUNNERS = Dict("prototype" => "run_1_prototype.jl",
    "4_v1" => "run_4_v1.jl", "8_v5" => "run_8_v5.jl", "9_full" => "run_9_full.jl")
haskey(RUNNERS, BENCH_CODE) || error("Unknown implementation: $BENCH_CODE")
const LOAD_WALL_S = @elapsed include(joinpath(@__DIR__, RUNNERS[BENCH_CODE]))
const BENCH_ROOT = joinpath(@__DIR__, "runtime_benchmarks")

function benchmark()
    trials = Dict{String, Any}[]
    for trial in 1:4
        GC.gc()
        wrapper_wall = @elapsed result = run_case(;
            output_root=joinpath(BENCH_ROOT, "trial_$trial"))
        @assert result.retcode == "Success"
        @assert result.helper_count == N_HELPERS
        @assert result.requested_target_orbits == ORBITS
        @assert result.schedule == string(SCHEDULE)
        @assert isapprox(result.initial_inclination_deg, TARGET_INC_DEG; atol=1e-10)
        @assert result.duration_s == ORBITS * result.target_period_s
        @assert isfinite(result.runtime_s) && result.runtime_s > 0
        @assert wrapper_wall >= result.runtime_s
        metrics = Dict{String, Any}(string(k) => v for (k, v) in pairs(result))
        metrics["trial"] = trial
        metrics["wrapper_wall_s"] = wrapper_wall
        push!(trials, metrics)
    end
    # Fresh state/model creation must yield the same output on repeated calls.
    for trial in trials[2:end], key in ("final_sma_m", "final_inclination_deg",
        "final_eccentricity", "link_activations", "active_updates", "final_helper")
        @assert trial[key] == trials[1][key] "$BENCH_CODE repeated-run mismatch: $key"
    end
    warm = [t["runtime_s"] for t in trials[2:end]]
    summary = Dict(
        "code" => BENCH_CODE, "helper_count" => N_HELPERS, "orbits" => ORBITS,
        "schedule" => string(SCHEDULE), "target_inclination_deg" => TARGET_INC_DEG,
        "completed_at" => string(now()), "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "cpu_model" => first(Sys.cpu_info()).model,
        "runner_load_wall_s" => LOAD_WALL_S,
        "first_simulation_call_wall_s" => trials[1]["runtime_s"],
        "warm_simulation_median_s" => median(warm),
        "warm_simulation_min_s" => minimum(warm),
        "warm_simulation_max_s" => maximum(warm),
        "warm_wrapper_median_s" => median([t["wrapper_wall_s"] for t in trials[2:end]]),
        "trials" => trials)
    dir = joinpath(BENCH_ROOT, CASE_DIR)
    mkpath(dir)
    open(joinpath(dir, "$BENCH_CODE.toml"), "w") do io
        TOML.print(io, summary; sorted=true)
    end
    println("BENCHMARK: ", BENCH_CODE, " orbits=", ORBITS,
        " first=", summary["first_simulation_call_wall_s"],
        " warm_median=", median(warm), " warm_range=", extrema(warm),
        " wrapper_median=", summary["warm_wrapper_median_s"])
end

benchmark()
