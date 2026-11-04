ENV["SPACEAGORA_SKIP_GMAT_MATRIX"] = "1"
ENV["SPACEAGORA_TELEMETRY_SOLVER_MODE"] = get(ENV, "SPACEAGORA_TELEMETRY_SOLVER_MODE", "rodas5p")

include(joinpath(@__DIR__, "gmat_scenario_matrix.jl"))
