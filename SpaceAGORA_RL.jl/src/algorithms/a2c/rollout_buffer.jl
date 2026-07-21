struct A2CRolloutBatch
    observations::Matrix{Float32}
    actions::Vector{Int}
    returns::Vector{Float32}
end

Base.length(batch::A2CRolloutBatch) = length(batch.actions)

function compute_discounted_returns(rewards::AbstractMatrix{Float32},
                                    done::AbstractMatrix{Bool},
                                    valid::AbstractMatrix{Bool},
                                    bootstrap_values::AbstractVector{Float32},
                                    discount::Real)
    size(rewards) == size(done) == size(valid) ||
        throw(DimensionMismatch("rewards, done, and valid matrices must have the same size"))
    size(rewards, 1) == length(bootstrap_values) ||
        throw(DimensionMismatch("one bootstrap value is required per worker"))
    gamma = Float32(discount)
    returns = zeros(Float32, size(rewards))
    for worker in axes(rewards, 1)
        last_valid = findlast(t -> valid[worker, t], axes(rewards, 2))
        last_valid === nothing && continue
        running = bootstrap_values[worker]
        for t in last_valid:-1:first(axes(rewards, 2))
            valid[worker, t] || continue
            if done[worker, t]
                running = 0f0
            end
            running = rewards[worker, t] + gamma * running
            returns[worker, t] = running
        end
    end
    return returns
end

function flatten_rollout(observations::Array{Float32,3}, actions::AbstractMatrix{Int},
                         returns::AbstractMatrix{Float32}, valid::AbstractMatrix{Bool})
    obs_dim, n_workers, segment_length = size(observations)
    size(actions) == (n_workers, segment_length) ||
        throw(DimensionMismatch("actions must be n_workers by segment_length"))
    size(returns) == (n_workers, segment_length) ||
        throw(DimensionMismatch("returns must be n_workers by segment_length"))
    size(valid) == (n_workers, segment_length) ||
        throw(DimensionMismatch("valid must be n_workers by segment_length"))

    n = count(valid)
    obs_batch = Matrix{Float32}(undef, obs_dim, n)
    action_batch = Vector{Int}(undef, n)
    return_batch = Vector{Float32}(undef, n)
    col = 1
    for t in 1:segment_length, worker in 1:n_workers
        valid[worker, t] || continue
        obs_batch[:, col] .= @view observations[:, worker, t]
        action_batch[col] = actions[worker, t]
        return_batch[col] = returns[worker, t]
        col += 1
    end
    return A2CRolloutBatch(obs_batch, action_batch, return_batch)
end

