include(joinpath(@__DIR__, "telemetry_odyssey_tuner", "cli_and_plots.jl"))
include(joinpath(@__DIR__, "telemetry_odyssey_tuner", "candidates_and_evaluator.jl"))
include(joinpath(@__DIR__, "telemetry_odyssey_tuner", "reporting.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_tuner()
end
