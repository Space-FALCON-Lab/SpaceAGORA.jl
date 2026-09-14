Base.@kwdef struct AerobrakingPassSchedule
    strategy_by_pass::Dict{Int, AerobrakingStrategyKind} = Dict{Int, AerobrakingStrategyKind}()
end

function strategy_for_pass(schedule::AerobrakingPassSchedule, pass_number::Int, fallback::AerobrakingStrategyKind)
    return get(schedule.strategy_by_pass, pass_number, fallback)
end
