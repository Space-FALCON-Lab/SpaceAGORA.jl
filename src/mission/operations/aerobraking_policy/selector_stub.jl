struct DRLPolicyAdapterStub <: AbstractAerobrakingPolicySelector
end

struct DefaultAerobrakingPolicySelector <: AbstractAerobrakingPolicySelector
end

function select_strategy(::DRLPolicyAdapterStub, config::AerobrakingPolicyConfig, input)
    throw(ErrorException("Not implemented: DRLPolicyAdapterStub.select_strategy"))
end

function select_strategy(::DefaultAerobrakingPolicySelector, config::AerobrakingPolicyConfig, input)
    args = getfield(input, :args)
    if args isa AbstractDict && haskey(args, :aerobraking_strategy)
        token = lowercase(String(args[:aerobraking_strategy]))
        token in ("e-edg", "e_edg", "eedg", "energy") && return E_EDG
        token in ("t-edg", "t_edg", "tedg", "targeting") && return T_EDG
    end
    return config.default_strategy
end
