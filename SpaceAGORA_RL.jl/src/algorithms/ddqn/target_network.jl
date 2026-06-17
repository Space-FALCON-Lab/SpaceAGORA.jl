function maybe_update_target!(learner)
    learner.config.target_update <= 0 && return false
    if learner.train_steps > 0 && learner.train_steps % learner.config.target_update == 0
        copy_network!(learner.target, learner.online)
        return true
    end
    return false
end
