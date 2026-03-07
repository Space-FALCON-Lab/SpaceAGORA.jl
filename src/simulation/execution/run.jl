using .SimulationModel

include(joinpath(@__DIR__, "dispatch.jl"))

function execute_orbital_elements_campaign(args::SimulationConfiguration; isolate_state::Bool=true)
    return execute_campaign(args; isolate_state=isolate_state)
end

function execute_vgamma_campaign(args::SimulationConfiguration; isolate_state::Bool=true)
    return execute_campaign(args; isolate_state=isolate_state)
end

function execute_ae_campaign(args::SimulationConfiguration; isolate_state::Bool=true)
    return execute_campaign(args; isolate_state=isolate_state)
end

function execute_analysis(args::SimulationConfiguration; isolate_state::Bool=true)
    return execute_campaign(args; isolate_state=isolate_state)
end
