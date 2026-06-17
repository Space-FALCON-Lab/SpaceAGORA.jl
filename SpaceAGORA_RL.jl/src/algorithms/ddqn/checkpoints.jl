function save_checkpoint(path::AbstractString, learner::DDQNLearner; manifest=nothing)
    mkpath(dirname(path))
    payload = Dict(
        :online => learner.online,
        :target => learner.target,
        :optimizer => learner.optimizer,
        :config => learner.config,
        :schedule => learner.schedule,
        :global_step => learner.global_step,
        :train_steps => learner.train_steps,
        :last_loss => learner.last_loss,
        :action_table => copy(PAPER_ACTIONS_MPS),
        :manifest => manifest,
    )
    serialize(path, payload)
    return path
end

load_checkpoint(path::AbstractString) = deserialize(path)
