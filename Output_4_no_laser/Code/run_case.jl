include(normpath(joinpath(@__DIR__, "..", "..",
    "Output_3_rangeOnly_vs_positiveAlongTrack", "Code", "run_case.jl")))

if abspath(PROGRAM_FILE) == @__FILE__
    run_worker(ARGS; allowed_policies=(:no_laser,))
end
