#!/usr/bin/env julia

if !isdefined(Main, :SpaceAGORACalibration)
    include(normpath(joinpath(@__DIR__, "..", "src", "SpaceAGORACalibration.jl")))
end
using .SpaceAGORACalibration

function main(args::Vector{String})
    spec_path = isempty(args) ? joinpath(@__DIR__, "..", "examples", "demo_spec.toml") : args[1]
    spec = load_spec(spec_path)
    backend = MockBackend()
    result = run_calibration(spec, backend)
    println("run_id=" * result.run_id)
    println("run_dir=" * result.run_dir)
    println("report=" * result.report_path)
    return nothing
end

main(copy(ARGS))
