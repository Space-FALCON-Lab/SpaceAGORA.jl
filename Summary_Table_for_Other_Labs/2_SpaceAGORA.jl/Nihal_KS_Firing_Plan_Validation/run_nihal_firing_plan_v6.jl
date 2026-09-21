#!/usr/bin/env julia

include(joinpath(@__DIR__, "run_nihal_firing_plan_v4.jl"))

const V6_INPUT_DIR = joinpath(@__DIR__, "Input", "run_nihal_firing_plan_v6")
const V6_OUTPUT_DIR = joinpath(@__DIR__, "Output", "run_nihal_firing_plan_v6")
const V6_IC_PATH = joinpath(@__DIR__, "Input", "initial_conditions.csv")

function main_v6(arguments=ARGS)
    length(arguments) <= 1 || error("Usage: run_nihal_firing_plan_v6.jl [--check-inputs | initial_conditions.csv]")
    check_only = arguments == ["--check-inputs"]
    ic_path = isempty(arguments) || check_only ? V6_IC_PATH : abspath(only(arguments))
    _validate_v4_initial_conditions(ic_path)
    schedule_paths = sort(filter(path -> occursin(V4_SCHEDULE_PATTERN, basename(path)),
        readdir(V6_INPUT_DIR; join=true)))
    isempty(schedule_paths) && error("No schedule_*.csv files found in $V6_INPUT_DIR.")
    checks = DataFrame(_validate_v4_schedule.(schedule_paths))
    mkpath(V6_OUTPUT_DIR)
    CSV.write(joinpath(V6_OUTPUT_DIR, "schedule_validation.csv"), checks)
    show(stdout, MIME"text/plain"(), checks)
    println("\nV6 input checks passed. Initial conditions: $ic_path")
    check_only && return nothing
    main_newnew(; input_dir=V6_INPUT_DIR, output_dir=V6_OUTPUT_DIR,
        variants=("",), ic_path=ic_path, schedule_pattern=V4_SCHEDULE_PATTERN,
        plan_loader=_load_v4_firing_plan)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    if ARGS != ["--check-inputs"]
        include(joinpath(@__DIR__, "run_nihal_firing_plan_v3.jl"))
    end
    main_v6()
end