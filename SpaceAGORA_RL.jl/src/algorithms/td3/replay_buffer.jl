struct ContinuousTransition
    observation::Vector{Float32}
    action::Vector{Float32}
    reward::Float32
    next_observation::Vector{Float32}
    terminated::Bool
    truncated::Bool
    info_index::Int
end

ContinuousTransition(observation::Vector{Float32}, action::Real, reward::Real,
                     next_observation::Vector{Float32}, terminated::Bool,
                     truncated::Bool, info_index::Integer) =
    ContinuousTransition(observation, Float32[action], Float32(reward), next_observation,
                         terminated, truncated, Int(info_index))

function continuous_transition(transition::Transition, action::AerobrakingAction)
    return ContinuousTransition(
        transition.observation,
        action.delta_v_mps,
        transition.reward,
        transition.next_observation,
        transition.terminated,
        transition.truncated,
        transition.info_index,
    )
end

mutable struct ContinuousReplayBuffer
    observations::Matrix{Float32}
    actions::Matrix{Float32}
    rewards::Vector{Float32}
    next_observations::Matrix{Float32}
    terminated::Vector{Bool}
    truncated::Vector{Bool}
    info_index::Vector{Int}
    position::Int
    size::Int
end

function ContinuousReplayBuffer(obs_dim::Integer, action_dim::Integer, capacity::Integer)
    capacity > 0 || throw(ArgumentError("continuous replay capacity must be positive"))
    obs_dim > 0 || throw(ArgumentError("continuous replay observation dimension must be positive"))
    action_dim > 0 || throw(ArgumentError("continuous replay action dimension must be positive"))
    return ContinuousReplayBuffer(
        zeros(Float32, obs_dim, capacity),
        zeros(Float32, action_dim, capacity),
        zeros(Float32, capacity),
        zeros(Float32, obs_dim, capacity),
        falses(capacity),
        falses(capacity),
        zeros(Int, capacity),
        1,
        0,
    )
end

Base.length(buffer::ContinuousReplayBuffer) = buffer.size
capacity(buffer::ContinuousReplayBuffer) = size(buffer.actions, 2)

function Base.push!(buffer::ContinuousReplayBuffer, transition::ContinuousTransition)
    size(buffer.observations, 1) == length(transition.observation) ||
        throw(DimensionMismatch("continuous transition observation dimension does not match replay"))
    size(buffer.actions, 1) == length(transition.action) ||
        throw(DimensionMismatch("continuous transition action dimension does not match replay"))
    size(buffer.next_observations, 1) == length(transition.next_observation) ||
        throw(DimensionMismatch("continuous transition next-observation dimension does not match replay"))
    idx = buffer.position
    buffer.observations[:, idx] .= transition.observation
    buffer.actions[:, idx] .= transition.action
    buffer.rewards[idx] = transition.reward
    buffer.next_observations[:, idx] .= transition.next_observation
    buffer.terminated[idx] = transition.terminated
    buffer.truncated[idx] = transition.truncated
    buffer.info_index[idx] = transition.info_index
    buffer.position = idx == capacity(buffer) ? 1 : idx + 1
    buffer.size = min(buffer.size + 1, capacity(buffer))
    return buffer
end

struct ContinuousReplayBatch
    observations::Matrix{Float32}
    actions::Matrix{Float32}
    rewards::Vector{Float32}
    next_observations::Matrix{Float32}
    terminated::Vector{Bool}
    truncated::Vector{Bool}
    info_index::Vector{Int}
end

function sample_batch(buffer::ContinuousReplayBuffer, batch_size::Integer,
                      rng::AbstractRNG)
    buffer.size >= batch_size || throw(ArgumentError("not enough continuous replay samples"))
    indices = rand(rng, 1:buffer.size, batch_size)
    return ContinuousReplayBatch(
        buffer.observations[:, indices],
        buffer.actions[:, indices],
        buffer.rewards[indices],
        buffer.next_observations[:, indices],
        buffer.terminated[indices],
        buffer.truncated[indices],
        buffer.info_index[indices],
    )
end
