struct MaskedTransition
    observation::Vector{Float32}
    action_mask::BitVector
    action_index::Int
    reward::Float32
    next_observation::Vector{Float32}
    next_action_mask::BitVector
    terminated::Bool
    truncated::Bool
    info_index::Int
end

mutable struct MaskedReplayBuffer
    observations::Matrix{Float32}
    action_masks::BitMatrix
    actions::Vector{Int}
    rewards::Vector{Float32}
    next_observations::Matrix{Float32}
    next_action_masks::BitMatrix
    terminated::Vector{Bool}
    truncated::Vector{Bool}
    info_index::Vector{Int}
    position::Int
    size::Int
end

function MaskedReplayBuffer(obs_dim::Integer, action_dim::Integer, capacity::Integer)
    obs_dim > 0 || throw(ArgumentError("observation dimension must be positive"))
    action_dim > 0 || throw(ArgumentError("action dimension must be positive"))
    capacity > 0 || throw(ArgumentError("replay capacity must be positive"))
    return MaskedReplayBuffer(
        zeros(Float32, obs_dim, capacity), falses(action_dim, capacity),
        zeros(Int, capacity), zeros(Float32, capacity),
        zeros(Float32, obs_dim, capacity), falses(action_dim, capacity),
        falses(capacity), falses(capacity), zeros(Int, capacity), 1, 0,
    )
end

Base.length(buffer::MaskedReplayBuffer) = buffer.size
capacity(buffer::MaskedReplayBuffer) = length(buffer.actions)

function Base.push!(buffer::MaskedReplayBuffer, transition::MaskedTransition)
    size(buffer.observations, 1) == length(transition.observation) ||
        throw(DimensionMismatch("observation length does not match replay buffer"))
    size(buffer.action_masks, 1) == length(transition.action_mask) ||
        throw(DimensionMismatch("action-mask length does not match replay buffer"))
    transition.action_mask[transition.action_index] ||
        throw(ArgumentError("transition action is invalid under its action mask"))
    idx = buffer.position
    buffer.observations[:, idx] .= transition.observation
    buffer.action_masks[:, idx] .= transition.action_mask
    buffer.actions[idx] = transition.action_index
    buffer.rewards[idx] = transition.reward
    buffer.next_observations[:, idx] .= transition.next_observation
    buffer.next_action_masks[:, idx] .= transition.next_action_mask
    buffer.terminated[idx] = transition.terminated
    buffer.truncated[idx] = transition.truncated
    buffer.info_index[idx] = transition.info_index
    buffer.position = idx == capacity(buffer) ? 1 : idx + 1
    buffer.size = min(buffer.size + 1, capacity(buffer))
    return buffer
end

struct MaskedReplayBatch
    observations::Matrix{Float32}
    action_masks::BitMatrix
    actions::Vector{Int}
    rewards::Vector{Float32}
    next_observations::Matrix{Float32}
    next_action_masks::BitMatrix
    terminated::Vector{Bool}
    truncated::Vector{Bool}
    info_index::Vector{Int}
end

function sample_batch(buffer::MaskedReplayBuffer, batch_size::Integer, rng::AbstractRNG)
    buffer.size >= batch_size || throw(ArgumentError("not enough replay samples"))
    idx = rand(rng, 1:buffer.size, batch_size)
    return MaskedReplayBatch(
        buffer.observations[:, idx], buffer.action_masks[:, idx],
        buffer.actions[idx], buffer.rewards[idx],
        buffer.next_observations[:, idx], buffer.next_action_masks[:, idx],
        buffer.terminated[idx], buffer.truncated[idx], buffer.info_index[idx],
    )
end
