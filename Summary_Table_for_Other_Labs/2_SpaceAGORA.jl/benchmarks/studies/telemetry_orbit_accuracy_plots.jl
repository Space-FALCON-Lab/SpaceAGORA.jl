const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CANONICAL_PLOTS = joinpath(REPO_ROOT, "scripts", "plotting", "telemetry_orbit_accuracy_plots.jl")

isfile(CANONICAL_PLOTS) || throw(ArgumentError("Missing canonical telemetry plotting script: $CANONICAL_PLOTS"))
include(CANONICAL_PLOTS)

main(copy(ARGS))
