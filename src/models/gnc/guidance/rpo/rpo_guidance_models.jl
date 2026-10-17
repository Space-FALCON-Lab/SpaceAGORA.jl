"""
Guidance model state for RPO planning, tracking, and optional replanning.

Deprecated: a run that configures this model warns before it starts. Build RPO
runs with `make_rpo_configuration` and a planner instead (`HYPRRPOPlanner` for
HYPR planning), which configures a `PlannerGuidance`. That path also covers the
model's replanning spheres (`geometry_provider=RPOReplanningSpheres(...)`), its
safe hold (`on_rejection=:install_fallback`) and retiming
(`retime_tracking_error_m`). The model will be removed in a later release.

Each run starts the goal, the plan buffer and the replanning record from their
values before the model's first run (`initialize_guidance!`).
"""
Base.@kwdef mutable struct RPOGuidanceModel <: AbstractGuidanceModel
    chaser_idx::Int = 1
    target_idx::Int = 2
    goal_rtn::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    geometry::Any = nothing
    plan_buffer::RPOPlanBuffer = RPOPlanBuffer()
    pso_config::Any = nothing
    safe_distance_m::Float64 = 0.0
    force_replan::Bool = false
    replanning_config::Any = nothing
    replanning_events::Vector{NamedTuple} = NamedTuple[]
    replan_count::Int = 0
    retime_count::Int = 0
    safe_hold_count::Int = 0
    replan_failure_count::Int = 0
    last_replanning_time_s::Float64 = -Inf
    last_replanning_signature::UInt = UInt(0)
    replanning_persistence_count::Int = 0
    replanning_phase::Symbol = :tracking
    pending_replan::Any = nothing
    hold_position_tolerance_m::Float64 = 0.02
    hold_speed_tolerance_mps::Float64 = 0.002
    # The run-written fields and plan buffer before the first run, restored at each run's start.
    run_start::Any = nothing
end
