include(joinpath(@__DIR__, "telemetry_hybrid_tuner", "cli_and_specs.jl"))
include(joinpath(@__DIR__, "telemetry_hybrid_tuner", "search_strategy.jl"))
include(joinpath(@__DIR__, "telemetry_hybrid_tuner", "manifest_and_evaluator.jl"))
include(joinpath(@__DIR__, "telemetry_hybrid_tuner", "reporting.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_hybrid_tuner()
end
