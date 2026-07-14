# Regenerates a reproducible-propagation ("golden") fixture under
# test/golden/<scenario_name>/.
#
# Usage:
#   julia --project=. scripts/regenerate_golden.jl <scenario_name> [spacecraft_count]
#
# Only regenerate a fixture when the underlying physics/model change producing
# a different trajectory is intentional -- treat a diff to a fixture.csv as
# something that needs extra scrutiny in review, not a routine update. If a
# scenario's config.jl didn't change and the fixture still differs, that's a
# regression, not a reason to regenerate.

using Test
using CSV
using DataFrames
using TOML
using Dates

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "test", "helpers", "bootstrap.jl"))

function regenerate_golden(
    scenario_name::AbstractString;
    spacecraft_count::Int=1,
    n_samples::Int=9,
    pos_atol_m::Float64=50.0,
    vel_atol_mps::Float64=0.05
)
    scenario_dir = joinpath(GOLDEN_ROOT, scenario_name)
    config_path = joinpath(scenario_dir, "config.jl")
    isfile(config_path) || error("No config.jl for scenario '$scenario_name' at $config_path")
    include(config_path)

    build_fn_name = Symbol("build_", scenario_name, "_config")
    isdefined(Main, build_fn_name) ||
        error("config.jl for '$scenario_name' must define a function named `$build_fn_name`")
    build_fn = getfield(Main, build_fn_name)

    # config.jl was just `include`d dynamically from inside this function, so
    # its method exists in a newer world than this function was called in.
    args = Base.invokelatest(build_fn)
    df = run_case(args)
    mission_time = args.mission_configuration.mission_time
    sample_times = collect(range(0.0, mission_time, length=n_samples))
    times = Vector{Float64}(df.time)

    ordered_cols = Symbol[:time, :pos_atol_m, :vel_atol_mps]
    columns = Dict{Symbol, Vector{Float64}}(
        :time => sample_times,
        :pos_atol_m => fill(pos_atol_m, n_samples),
        :vel_atol_mps => fill(vel_atol_mps, n_samples)
    )
    for k in 1:spacecraft_count
        for c in 1:3
            pos_col = Symbol("sc$(k)_pos_$(c)")
            push!(ordered_cols, pos_col)
            columns[pos_col] = [interp_linear(times, getproperty(df, pos_col), t) for t in sample_times]
        end
        for c in 1:3
            vel_col = Symbol("sc$(k)_vel_$(c)")
            push!(ordered_cols, vel_col)
            columns[vel_col] = [interp_linear(times, getproperty(df, vel_col), t) for t in sample_times]
        end
    end

    fixture = DataFrame()
    for col in ordered_cols
        fixture[!, col] = columns[col]
    end
    fixture_path = joinpath(scenario_dir, "fixture.csv")
    CSV.write(fixture_path, fixture)

    metadata_path = joinpath(scenario_dir, "metadata.toml")
    metadata = isfile(metadata_path) ? TOML.parsefile(metadata_path) : Dict{String, Any}()
    metadata["git_sha"] = strip(read(`git -C $REPO_ROOT rev-parse HEAD`, String))
    metadata["julia_version"] = string(VERSION)
    metadata["generated_at"] = string(Dates.now())
    get!(metadata, "tier", "pr")
    open(metadata_path, "w") do io
        TOML.print(io, metadata)
    end

    println("Wrote $(nrow(fixture)) rows to $fixture_path")
    println("Updated provenance in $metadata_path")
    return fixture_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error("Usage: julia --project=. scripts/regenerate_golden.jl <scenario_name> [spacecraft_count]")
    scenario_name = ARGS[1]
    spacecraft_count = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
    regenerate_golden(scenario_name; spacecraft_count=spacecraft_count)
end
