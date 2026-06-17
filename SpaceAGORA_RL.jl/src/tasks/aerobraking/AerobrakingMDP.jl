struct AerobrakingMDP <: AbstractMDP
    config::AerobrakingScenarioConfig
end

AerobrakingMDP(; kwargs...) = AerobrakingMDP(default_aerobraking_config(; kwargs...))
