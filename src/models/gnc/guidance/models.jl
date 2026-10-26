module GuidanceModels

using ..AbstractTypes: AbstractGuidanceModel
using StaticArrays

export AerobrakingCampaignPropulsiveManeuverGuidanceModel
export ApoapsisTargetPeriapsisRaiseGuidanceModel
export RPOPlan, RPOPlanBuffer, update_rpo_plan_buffer!

include(joinpath(@__DIR__, "maneuvers", "thruster_guidance_models.jl"))
include(joinpath(@__DIR__, "rpo", "rpo_plan_buffer.jl"))

end # module GuidanceModels
