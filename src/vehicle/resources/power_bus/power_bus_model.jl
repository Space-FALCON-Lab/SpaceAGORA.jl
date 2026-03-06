module PowerBusModel

using ..Resources: AbstractResourceModel, ResourceInputs, ResourceOutputs, ResourceState
import ..Resources: initialize_resource_state, step_resource!

export PowerBusResourceModel

Base.@kwdef struct PowerBusResourceModel <: AbstractResourceModel
    nominal_voltage_v::Float64 = 0.0
    max_current_a::Float64 = 0.0
end

function initialize_resource_state(::PowerBusResourceModel)
    throw(ErrorException("Not implemented: initialize_resource_state(::PowerBusResourceModel)"))
end

function step_resource!(::PowerBusResourceModel, ::ResourceState, ::ResourceInputs)::ResourceOutputs
    throw(ErrorException("Not implemented: step_resource!(::PowerBusResourceModel, ::ResourceState, ::ResourceInputs)"))
end

end # module PowerBusModel
