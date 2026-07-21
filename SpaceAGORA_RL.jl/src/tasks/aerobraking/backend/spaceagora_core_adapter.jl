struct SpaceAGORACoreAdapter
    mode::Symbol
end

SpaceAGORACoreAdapter() = SpaceAGORACoreAdapter(:paper_surrogate)

function propagate_pass(adapter::SpaceAGORACoreAdapter, config, state, action::AerobrakingAction, rng::AbstractRNG)
    adapter.mode == :paper_surrogate && return paper_surrogate_pass(config, state, action, rng)
    if _is_spaceagora_live_backend(adapter.mode)
        return _spaceagora_physics_propagate_single_pass(config, state, action)
    end
    throw(ArgumentError("unsupported aerobraking backend mode: $(adapter.mode)"))
end
