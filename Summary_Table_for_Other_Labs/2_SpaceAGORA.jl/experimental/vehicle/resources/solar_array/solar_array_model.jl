module SolarArrayModel

using ..Resources: AbstractResourceModel, ResourceInputs, ResourceOutputs, ResourceState
import ..Resources: initialize_resource_state, step_resource!

export SolarArrayResourceModel

Base.@kwdef struct SolarArrayResourceModel <: AbstractResourceModel
    max_generation_w::Float64 = 0.0
    efficiency::Float64 = 1.0
end

function initialize_resource_state(::SolarArrayResourceModel)
    throw(ErrorException("Not implemented: initialize_resource_state(::SolarArrayResourceModel)"))
end

function step_resource!(::SolarArrayResourceModel, ::ResourceState, ::ResourceInputs)::ResourceOutputs
    throw(ErrorException("Not implemented: step_resource!(::SolarArrayResourceModel, ::ResourceState, ::ResourceInputs)"))
end

end # module SolarArrayModel
