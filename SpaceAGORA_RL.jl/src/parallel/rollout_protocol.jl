struct WorkerTransitionMessage
    worker_id::Int
    transition::Transition
    result::AerobrakingStepResult
    summary::Union{Nothing,EpisodeSummary}
end

struct WorkerActionReply
    action_index::Int
    stop::Bool
end
