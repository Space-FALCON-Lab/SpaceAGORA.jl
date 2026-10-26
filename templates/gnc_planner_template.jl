# A planning algorithm ("Writing a GNC Algorithm" in the docs, the planning
# profile): it proposes a product for a request, and its family, trusted code in
# SpaceAGORA, validates the proposal and installs only a validated copy. This one
# is an RPO reference planner; replace the body of `plan_rpo!` with your method.
# examples/gnc_algorithm_env/ holds a complete planner of each family, run in CI.
#
# For another family, subtype that family's algorithm type (`AbstractBurnPlanner`,
# `AbstractJointReferencePlanner`), implement `propose!(state, algorithm, request,
# rng)` returning an `AlgorithmResult` with the family's product, and
# run it in a simulation with `AlgorithmGuidance(algorithm=..., spacecraft_index=...)`
# among the guidance effectors.
module ExampleGNCPlannerTemplate

using Random: AbstractRNG
using SpaceAGORA

struct ExamplePlanner <: AbstractRPOPlanner
    # The planner's settings. Run state goes in what `initialize_planner` returns.
end

# What the planner supports: truth states in the target's RTN frame, no
# retiming (declare `retiming=true` and add `retime_rpo!` together).
SpaceAGORA.planner_capabilities(::ExamplePlanner) =
    RPOPlannerCapabilities(state_sources=(:truth,), frames=(:target_rtn,))

# Whether `plan_rpo!` draws from `rng`. Declare `false` for a planner that draws
# nothing; the conformance checks then require the stream to stay untouched.
SpaceAGORA.consumes_rng(::ExamplePlanner) = true

# Fresh run state for each run; the default is `nothing`.
SpaceAGORA.initialize_planner(::ExamplePlanner, run_context) = nothing

# Propose a reference for `request`: an `RPOPlanningResult` with status
# `:candidate` and an `RPOReference` from `request.x_rtn` to `request.goal_rtn_m`
# within `request.constraints`, or a reason it could not (`:infeasible`,
# `:unsupported`, `:failed`). Draw random numbers only from `rng`. Never validate
# or install the reference here: the family does.
function SpaceAGORA.plan_rpo!(state, planner::ExamplePlanner, request::RPOPlanningRequest, rng::AbstractRNG)
    return RPOPlanningResult(request_id=request.request_id, status=:infeasible,
        termination=:not_implemented)
end

end # module ExampleGNCPlannerTemplate
