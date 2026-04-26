struct DRLPolicyAdapterStub <: AbstractAerobrakingPolicySelector
end

struct DefaultAerobrakingPolicySelector <: AbstractAerobrakingPolicySelector
end

function select_strategy(::DRLPolicyAdapterStub, config::AerobrakingPolicyConfig, input)
    throw(ErrorException("Not implemented: DRLPolicyAdapterStub.select_strategy"))
end

function select_strategy(::DefaultAerobrakingPolicySelector, config::AerobrakingPolicyConfig, input)
    return config.default_strategy
end
