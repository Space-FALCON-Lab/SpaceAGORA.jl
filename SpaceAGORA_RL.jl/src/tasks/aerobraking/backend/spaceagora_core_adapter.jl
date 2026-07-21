struct SpaceAGORACoreAdapter
    mode::Symbol
end

SpaceAGORACoreAdapter() = SpaceAGORACoreAdapter(:paper_surrogate)

function propagate_pass(adapter::SpaceAGORACoreAdapter, config, state, action::AerobrakingAction, rng::AbstractRNG)
    if adapter.mode == :paper_surrogate
        return paper_surrogate_pass(config, state, action, rng)
    elseif adapter.mode == :spaceagora_marsgram
        return spaceagora_marsgram_pass(config, state, action, rng, adapter.mode)
    elseif adapter.mode == :spaceagora_physics || adapter.mode == :spaceagora_full_physics
        return spaceagora_physics_pass(config, state, action, rng)
    end
    throw(ArgumentError(
        "unsupported aerobraking backend mode: $(adapter.mode). " *
        "Supported modes: paper_surrogate, spaceagora_marsgram, spaceagora_physics"
    ))
end
