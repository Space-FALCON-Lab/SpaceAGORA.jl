# Canonical aggregator: no behavior ownership.
"""Run ownership and accepted-update installation for opt-in RPO planners."""
module RPOPlannerLifecycle
using Random: MersenneTwister
using SHA: sha256
using LinearAlgebra: norm, dot
using ..SimulationModel
import ..RPOPlannerInterfaces as P
import ..SimulationLifecycle as L
import ..DirectRPOPlanning: DirectRPOPlanner
const S = SimulationModel

include("planner_lifecycle.jl")
include("planner_configuration.jl")
end
