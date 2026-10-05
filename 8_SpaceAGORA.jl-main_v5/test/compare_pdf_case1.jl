using Test

# Step 1: Load the shared Case 1 example implementation and its run/test helpers.
include(joinpath(@__DIR__, "..", "II_examples", "compare_pdf_case1.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
    # Step 2: Parse optional implementation, orbit-count, and timestep arguments.
    mode = Symbol(isempty(ARGS) ? "v5" : ARGS[1])
    orbits = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 10.0
    dtmax = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 10.0

    # Step 3: Run the laser-on versus laser-off comparison test.
    compare_case1(mode; orbits, dtmax)
end