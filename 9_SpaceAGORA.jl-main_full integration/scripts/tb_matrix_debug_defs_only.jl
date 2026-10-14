# Compatibility include path. Definitions are maintained in the matrix itself;
# this loader never generates or rewrites source. Access helpers through
# ScenarioMatrixDebugSupport, whose inclusion does not execute any campaign.
if !isdefined(@__MODULE__, :ScenarioMatrixDebugSupport)
    include(joinpath(@__DIR__, "scenario_matrix_debug_support.jl"))
end
