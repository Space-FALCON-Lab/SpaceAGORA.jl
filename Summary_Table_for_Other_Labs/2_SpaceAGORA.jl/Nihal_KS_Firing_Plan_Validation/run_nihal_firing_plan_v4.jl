#!/usr/bin/env julia

using CSV
using DataFrames

const V4_INPUT_DIR = joinpath(@__DIR__, "Input", "run_nihal_firing_plan_v4")
const V4_OUTPUT_DIR = joinpath(@__DIR__, "Output", "run_nihal_firing_plan_v4")
const V4_IC_PATH = joinpath(V4_INPUT_DIR, "initial_conditions.csv")
const V4_MISSION_TIME_S = 31535.592545587482
const V4_SCHEDULE_PATTERN = r"^schedule_(.+)\.csv$"

function _validate_v4_schedule(path::String)
    schedule = CSV.read(path, DataFrame)
    names(schedule) == vcat(["start_time_s"], ["$(helper)_to_1" for helper in 2:27]) ||
        error("Unexpected v4 schedule columns in $path.")
    nrow(schedule) >= 2 || error("Expected at least two schedule intervals in $path.")
    times = Float64.(schedule.start_time_s)
    all(isfinite, times) && first(times) == 0.0 && all(diff(times) .> 0.0) ||
        error("Schedule times must be finite, start at zero, and strictly increase in $path.")
    last(times) < V4_MISSION_TIME_S || error("Schedule exceeds the v4 mission duration in $path.")
    flags = Matrix(schedule[:, 2:end])
    all(value -> !ismissing(value) && value in (0, 1), flags) ||
        error("Schedule flags must be binary in $path.")
    all(sum(flags; dims=2) .== 1) || error("Expected exactly one active helper per interval in $path.")
    return (
        schedule=match(V4_SCHEDULE_PATTERN, basename(path)).captures[1],
        intervals=length(times), spacecraft=27,
        minimum_interval_s=minimum(diff(times)), maximum_interval_s=maximum(diff(times)),
        last_start_s=last(times), mission_time_s=V4_MISSION_TIME_S,
    )
end

function _load_v4_firing_plan(ic_df::DataFrame, path::String)
    _validate_v4_schedule(path)
    plan = _load_new_firing_plan(ic_df, path)
    return merge(plan, (mission_time_s=V4_MISSION_TIME_S,))
end

function _validate_v4_initial_conditions(path::String)
    isfile(path) || error("Missing initial conditions: $path. Supply initial_conditions.csv in the v4 input folder or pass its path as the sole argument. No simulation has been run.")
    initial_conditions = CSV.read(path, DataFrame)
    required = ["satellite", "a_km", "e", "i_deg", "raan_deg", "omega_deg", "M_deg"]
    all(column -> column in names(initial_conditions), required) ||
        error("Initial conditions must contain: $(join(required, ", ")).")
    all(value -> value isa Real && isfinite(value), Matrix(initial_conditions[:, required])) ||
        error("Initial conditions must contain finite numeric values.")
    all(isinteger, initial_conditions.satellite) || error("Satellite IDs must be integers.")
    for satellite in 1:27
        count(==(satellite), initial_conditions.satellite) == 1 ||
            error("Expected exactly one initial-condition row for satellite $satellite.")
    end
    all(initial_conditions.a_km .> 0) && all(0 .<= initial_conditions.e .< 1) ||
        error("Initial conditions require positive semimajor axes and eccentricities in [0, 1).")
    return nothing
end

function main_v4(arguments=ARGS)
    length(arguments) <= 1 || error("Usage: run_nihal_firing_plan_v4.jl [--check-inputs | initial_conditions.csv]")
    check_only = arguments == ["--check-inputs"]
    ic_path = isempty(arguments) || check_only ? V4_IC_PATH : abspath(only(arguments))
    schedule_paths = sort(filter(path -> occursin(V4_SCHEDULE_PATTERN, basename(path)),
        readdir(V4_INPUT_DIR; join=true)))
    isempty(schedule_paths) && error("No schedule_*.csv files found in $V4_INPUT_DIR.")
    checks = DataFrame(_validate_v4_schedule.(schedule_paths))
    mkpath(V4_OUTPUT_DIR)
    CSV.write(joinpath(V4_OUTPUT_DIR, "schedule_validation.csv"), checks)
    show(stdout, MIME"text/plain"(), checks)
    println("\nSchedule checks passed. Output directory: $V4_OUTPUT_DIR")
    if check_only
        println("Schedule-only check: orbital simulations were not run. Initial conditions present: $(isfile(ic_path)).")
        return nothing
    end
    _validate_v4_initial_conditions(ic_path)
    println("Using initial conditions: $ic_path")
    main_newnew(; input_dir=V4_INPUT_DIR, output_dir=V4_OUTPUT_DIR,
        variants=("",), ic_path=ic_path, schedule_pattern=V4_SCHEDULE_PATTERN,
        plan_loader=_load_v4_firing_plan)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    if ARGS != ["--check-inputs"]
        _validate_v4_initial_conditions(isempty(ARGS) ? V4_IC_PATH : abspath(first(ARGS)))
        include(joinpath(@__DIR__, "run_nihal_firing_plan_v3.jl"))
    end
    main_v4()
end