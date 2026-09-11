module BatteryModel

using ..Resources: AbstractResourceModel, ResourceInputs, ResourceOutputs, ResourceState
import ..Resources: initialize_resource_state, step_resource!

export BatteryResourceModel

Base.@kwdef struct BatteryResourceModel <: AbstractResourceModel
    capacity_j::Float64 = 0.0
    min_state_of_charge::Float64 = 0.0
    max_state_of_charge::Float64 = 1.0
end

function initialize_resource_state(::BatteryResourceModel)
    throw(ErrorException("Not implemented: initialize_resource_state(::BatteryResourceModel)"))
end

function step_resource!(::BatteryResourceModel, ::ResourceState, ::ResourceInputs)::ResourceOutputs
    throw(ErrorException("Not implemented: step_resource!(::BatteryResourceModel, ::ResourceState, ::ResourceInputs)"))
end

end # module BatteryModel
