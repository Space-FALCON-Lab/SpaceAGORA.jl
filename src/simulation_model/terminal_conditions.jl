module TerminalConditions
    using ..SimulationModel

    export TimeTerminalCondition, OrbitTerminalCondition, check_terminal_condition
    export terminal_condition_value, estimated_terminal_time

    include("../utils/terminal_condition_types.jl")
end # module TerminalConditions
