"""
    valid_action_mask(config, state)

Return a Boolean vector identifying the discrete actions that are valid in the
current task state. Tasks with a fixed action space may return `trues(n)`.
"""
function valid_action_mask end

function _valid_action_indices(mask::AbstractVector{Bool})
    indices = findall(mask)
    isempty(indices) && throw(ArgumentError("an action mask must permit at least one action"))
    return indices
end

function _masked_argmax(values::AbstractVector{<:Real}, mask::AbstractVector{Bool})
    length(values) == length(mask) ||
        throw(DimensionMismatch("action mask length must match the Q-value vector"))
    valid = _valid_action_indices(mask)
    _, local_index = findmax(view(values, valid))
    return valid[local_index]
end
