@inline function strategy_from_kind(kind)
    if kind == E_EDG
        return EEdgStrategy()
    elseif kind == T_EDG
        return TEdgStrategy()
    end
    throw(ArgumentError("Unsupported aerobraking strategy kind: $(kind)"))
end

function dispatch_aerobraking_guidance(kind, input::AerobrakingGuidanceInput)
    return compute_aerobraking_guidance(strategy_from_kind(kind), input)
end

function dispatch_aerobraking_guidance(
    selector::AbstractAerobrakingPolicySelector,
    config::AerobrakingPolicyConfig,
    input::AerobrakingGuidanceInput,
)
    kind = select_strategy(selector, config, input)
    return dispatch_aerobraking_guidance(kind, input)
end
