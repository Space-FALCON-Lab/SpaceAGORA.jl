mutable struct A3CRollout
    transitions::Vector{Transition}
    actions_mps::Vector{Float32}
end

A3CRollout() = A3CRollout(Transition[], Float32[])
A3CRollout(transitions::Vector{Transition}) = A3CRollout(
    transitions,
    Float32[PAPER_ACTIONS_MPS[transition.action_index] for transition in transitions],
)

Base.length(rollout::A3CRollout) = length(rollout.transitions)
Base.isempty(rollout::A3CRollout) = isempty(rollout.transitions)

function Base.empty!(rollout::A3CRollout)
    empty!(rollout.transitions)
    empty!(rollout.actions_mps)
    return rollout
end

function push_a3c_transition!(rollout::A3CRollout, transition::Transition,
                              action::AerobrakingAction)
    push!(rollout.transitions, transition)
    push!(rollout.actions_mps, Float32(action.delta_v_mps))
    return rollout
end

struct A3CRolloutBatch
    observations::Matrix{Float32}
    actions::Vector{Int}
    returns::Vector{Float32}
    policy_version::Int
end

Base.length(batch::A3CRolloutBatch) = length(batch.actions)

struct ContinuousA3CRolloutBatch
    observations::Matrix{Float32}
    actions_mps::Vector{Float32}
    returns::Vector{Float32}
    policy_version::Int
end

Base.length(batch::ContinuousA3CRolloutBatch) = length(batch.actions_mps)
