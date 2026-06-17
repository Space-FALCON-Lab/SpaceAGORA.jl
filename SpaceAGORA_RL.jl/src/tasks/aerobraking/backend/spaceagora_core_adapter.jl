struct SpaceAGORACoreAdapter
    mode::Symbol
end

SpaceAGORACoreAdapter() = SpaceAGORACoreAdapter(:paper_surrogate)

function propagate_pass(adapter::SpaceAGORACoreAdapter, config, state, action::AerobrakingAction, rng::AbstractRNG)
    adapter.mode == :paper_surrogate || throw(ArgumentError("unsupported aerobraking backend mode: $(adapter.mode)"))
    return paper_surrogate_pass(config, state, action, rng)
end
