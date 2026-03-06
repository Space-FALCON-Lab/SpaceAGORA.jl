using .SimulationModel

@inline function _require_simulation_configuration(args)
    args isa SimulationConfiguration || throw(ArgumentError("execute_campaign expects SimulationConfiguration; got $(typeof(args))."))
    return args
end

"""
    execute_campaign(args::SimulationConfiguration; isolate_state=true, state=nothing)

Canonical campaign execution entrypoint. Delegates to the simulation engine.
"""
function execute_campaign(args::SimulationConfiguration; isolate_state::Bool=true, state=nothing)
    _ = state
    return SimulationEngine.run_simulation(args; isolate_state=isolate_state)
end

# Legacy overload kept for callsites that pass a second positional state argument.
function execute_campaign(args::SimulationConfiguration, state; isolate_state::Bool=true)
    _ = state
    return SimulationEngine.run_simulation(args; isolate_state=isolate_state)
end

function execute_campaign(args; isolate_state::Bool=true, state=nothing)
    _ = isolate_state
    _ = state
    _require_simulation_configuration(args)
    return nothing
end
