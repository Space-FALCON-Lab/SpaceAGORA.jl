module LoadModel

using ..Resources: AbstractResourceModel, ResourceInputs, ResourceOutputs, ResourceState
import ..Resources: initialize_resource_state, step_resource!

export LoadResourceModel

Base.@kwdef struct LoadResourceModel <: AbstractResourceModel
    nominal_load_w::Float64 = 0.0
    priority::Int = 0
end

function initialize_resource_state(::LoadResourceModel)
    throw(ErrorException("Not implemented: initialize_resource_state(::LoadResourceModel)"))
end

function step_resource!(::LoadResourceModel, ::ResourceState, ::ResourceInputs)::ResourceOutputs
    throw(ErrorException("Not implemented: step_resource!(::LoadResourceModel, ::ResourceState, ::ResourceInputs)"))
end

end # module LoadModel
