using Dates
using TOML

isempty(ARGS) || error("Usage: benchmark_v9_helper_count_sequence.jl (10, 20, 30 helpers; 100 orbits)")
append!(ARGS, ["10", "100", "gve_sma", "1"])
const LOAD_WALL_S = @elapsed include(joinpath(@__DIR__, "run_9_full.jl"))
const SEQUENCE_ROOT = joinpath(@__DIR__, "runtime_benchmarks", "v9_helper_count_sequence")

function benchmark_sequence()
    trials = Dict{String, Any}[]
    for (order, helpers) in enumerate((10, 20, 30))
        GC.gc()
        wrapper_wall = @elapsed result = run_case(;
            output_root=SEQUENCE_ROOT, helper_count=helpers)
        @assert result.retcode == "Success"
        @assert result.helper_count == helpers
        @assert result.requested_target_orbits == 100
        @assert result.schedule == "gve_sma"
        @assert isapprox(result.initial_inclination_deg, 1.0; atol=1e-10)
        @assert result.duration_s == 100 * result.target_period_s
        @assert isfinite(result.runtime_s) && result.runtime_s > 0
        @assert wrapper_wall >= result.runtime_s
        metrics = Dict{String, Any}(string(k) => v for (k, v) in pairs(result))
        metrics["sequence_order"] = order
        metrics["wrapper_wall_s"] = wrapper_wall
        push!(trials, metrics)
        println("SEQUENCE: order=", order, " helpers=", helpers,
            " simulation_s=", result.runtime_s, " wrapper_s=", wrapper_wall)
    end
    summary = Dict(
        "helper_counts" => [10, 20, 30], "calls_per_count" => 1,
        "same_process" => true, "process_id" => getpid(),
        "prior_warmup_simulations" => 0,
        "runner_load_wall_s" => LOAD_WALL_S,
        "completed_at" => string(now()), "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(), "cpu_model" => first(Sys.cpu_info()).model,
        "trials" => trials)
    open(joinpath(SEQUENCE_ROOT, "sequence.toml"), "w") do io
        TOML.print(io, summary; sorted=true)
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && benchmark_sequence()
