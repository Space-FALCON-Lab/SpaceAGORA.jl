module AerobrakingPolicy

@enum AerobrakingStrategyKind E_EDG T_EDG

abstract type AbstractAerobrakingPolicySelector end

Base.@kwdef struct AerobrakingPolicyConfig
    default_strategy::AerobrakingStrategyKind = E_EDG
end

include(joinpath(@__DIR__, "selector_stub.jl"))
include(joinpath(@__DIR__, "pass_schedule.jl"))

export AerobrakingStrategyKind
export E_EDG, T_EDG
export AbstractAerobrakingPolicySelector
export AerobrakingPolicyConfig
export DRLPolicyAdapterStub
export DefaultAerobrakingPolicySelector
export AerobrakingPassSchedule
export strategy_for_pass
export select_strategy

end # module AerobrakingPolicy
