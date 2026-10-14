"""
    SimulationResults

What `run_simulation(args; return_results=true)` returns.

- `table::DataFrame`: the results table, same columns and values as the results
  CSV, built in memory even when `simulation_settings.results=false`.
- `configuration`: the configuration that ran. Under the default
  `isolate_state=true` this is the deep copy, so mutated controller state and
  logs (for example the RPO command log) are reachable here and not on the
  caller's object.
- `files::Vector{String}`: paths of existing files the run wrote (results CSV,
  bundle, scene, checkpoint); empty when none.
- `solution`: the ODE solution when `return_solution=true` was also passed,
  otherwise `nothing`.
"""
struct SimulationResults{C, S}
    table::DataFrame
    configuration::C
    files::Vector{String}
    solution::S
end
