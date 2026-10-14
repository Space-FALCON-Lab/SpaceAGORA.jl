"""Adaptive collision-sampling controls tied to clearance and curvature."""
Base.@kwdef struct RPOAdaptiveSamplingSettings
    enabled::Bool = true
    max_ds_m::Float64 = 0.50
    far_clearance_m::Float64 = 1.0
    power::Float64 = 1.0
    safe_distance_fraction::Float64 = 0.5
    obstacle_guard_fraction::Float64 = 0.5
end
