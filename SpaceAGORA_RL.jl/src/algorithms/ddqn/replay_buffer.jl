mutable struct ReplayBuffer
    observations::Matrix{Float32}
    actions::Vector{Int}
    rewards::Vector{Float32}
    next_observations::Matrix{Float32}
    terminated::Vector{Bool}
    truncated::Vector{Bool}
    info_index::Vector{Int}
    position::Int
    size::Int
end

function ReplayBuffer(obs_dim::Integer, capacity::Integer)
    capacity > 0 || throw(ArgumentError("replay capacity must be positive"))
    obs_dim > 0 || throw(ArgumentError("observation dimension must be positive"))
    return ReplayBuffer(
        zeros(Float32, obs_dim, capacity),
        zeros(Int, capacity),
        zeros(Float32, capacity),
        zeros(Float32, obs_dim, capacity),
        falses(capacity),
        falses(capacity),
        zeros(Int, capacity),
        1,
        0,
    )
end

Base.length(buffer::ReplayBuffer) = buffer.size
capacity(buffer::ReplayBuffer) = length(buffer.actions)

function Base.push!(buffer::ReplayBuffer, transition::Transition)
    idx = buffer.position
    buffer.observations[:, idx] .= transition.observation
    buffer.actions[idx] = transition.action_index
    buffer.rewards[idx] = transition.reward
    buffer.next_observations[:, idx] .= transition.next_observation
    buffer.terminated[idx] = transition.terminated
    buffer.truncated[idx] = transition.truncated
    buffer.info_index[idx] = transition.info_index
    buffer.position = idx == capacity(buffer) ? 1 : idx + 1
    buffer.size = min(buffer.size + 1, capacity(buffer))
    return buffer
end

struct ReplayBatch
    observations::Matrix{Float32}
    actions::Vector{Int}
    rewards::Vector{Float32}
    next_observations::Matrix{Float32}
    terminated::Vector{Bool}
    truncated::Vector{Bool}
    info_index::Vector{Int}
end

function sample_batch(buffer::ReplayBuffer, batch_size::Integer, rng::AbstractRNG)
    buffer.size >= batch_size || throw(ArgumentError("not enough replay samples"))
    idx = rand(rng, 1:buffer.size, batch_size)
    return ReplayBatch(
        buffer.observations[:, idx],
        buffer.actions[idx],
        buffer.rewards[idx],
        buffer.next_observations[:, idx],
        buffer.terminated[idx],
        buffer.truncated[idx],
        buffer.info_index[idx],
    )
end
