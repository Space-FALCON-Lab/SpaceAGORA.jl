module Resources

export AbstractResourceModel
export ResourceState, ResourceInputs, ResourceOutputs
export initialize_resource_state, step_resource!
export BatteryResourceModel, SolarArrayResourceModel, PowerBusResourceModel, LoadResourceModel

abstract type AbstractResourceModel end

Base.@kwdef mutable struct ResourceState
    state_of_charge::Float64 = 1.0
    available_power_w::Float64 = 0.0
    allocated_power_w::Float64 = 0.0
    power_margin_w::Float64 = 0.0
    terminal_enabled::Bool = true
end

Base.@kwdef struct ResourceInputs
    epoch_s::Float64 = 0.0
    dt_s::Float64 = 0.0
    generation_w::Float64 = 0.0
    requested_load_w::Float64 = 0.0
end

Base.@kwdef struct ResourceOutputs
    delivered_power_w::Float64 = 0.0
    shed_load_w::Float64 = 0.0
    bus_voltage_v::Float64 = 0.0
end

function initialize_resource_state(model::AbstractResourceModel)
    throw(ErrorException("Not implemented: initialize_resource_state(::$(typeof(model)))"))
end

function step_resource!(model::AbstractResourceModel, state::ResourceState, inputs::ResourceInputs)
    throw(ErrorException("Not implemented: step_resource!(::$(typeof(model)), ::ResourceState, ::ResourceInputs)"))
end

include(joinpath(@__DIR__, "battery", "battery_model.jl"))
include(joinpath(@__DIR__, "solar_array", "solar_array_model.jl"))
include(joinpath(@__DIR__, "power_bus", "power_bus_model.jl"))
include(joinpath(@__DIR__, "loads", "load_model.jl"))

using .BatteryModel: BatteryResourceModel
using .SolarArrayModel: SolarArrayResourceModel
using .PowerBusModel: PowerBusResourceModel
using .LoadModel: LoadResourceModel

end # module Resources
