Base.@kwdef struct RPOHyPRRLTrainingConfig
    episodes::Int = 1_000
    seed::Int = 740
    checkpoint_every_episodes::Int = 100
    checkpoint_directory::String = joinpath(pwd(), "hypr_rl_checkpoints")
end

function rpo_hypr_rl_ddqn_config(config::RPOHyPRRLConfig;
                                 hidden_dim::Int=512,
                                 learning_rate::Real=1.0e-4,
                                 discount::Real=0.99,
                                 batch_size::Int=256,
                                 train_start::Int=5_000,
                                 replay_size::Int=500_000,
                                 target_update::Int=5_000)
    return DDQNConfig(
        learning_rate=Float64(learning_rate), discount=Float64(discount),
        batch_size=batch_size, train_start=train_start, replay_size=replay_size,
        target_update=target_update, hidden_dim=hidden_dim,
        obs_dim=rpo_hypr_rl_observation_dim(config),
        action_dim=action_count(config),
    )
end

function save_hypr_rl_checkpoint(path::AbstractString, learner::DDQNLearner,
                                 config::RPOHyPRRLConfig;
                                 training_metadata=nothing)
    return save_checkpoint(
        path, learner;
        task=:rpo_hypr_rl,
        action_table=rpo_editor_actions(config),
        task_metadata=(config=config, training=training_metadata),
    )
end

"""
Train masked DDQN on a HyPR-RL MDP. Pass a `scenario_sampler` to randomize the
start, goal, geometry, or spacecraft configuration between episodes.
"""
function train_hypr_rl!(mdp::RPOHyPRRLMDP;
                        training::RPOHyPRRLTrainingConfig=RPOHyPRRLTrainingConfig(),
                        ddqn_config::DDQNConfig=rpo_hypr_rl_ddqn_config(mdp.config),
                        schedule::EpsilonSchedule=EpsilonSchedule(
                            decay_start_step=ddqn_config.train_start,
                        ),
                        learner::DDQNLearner=DDQNLearner(
                            MersenneTwister(training.seed), ddqn_config;
                            schedule=schedule,
                        ),
                        scenario_sampler=nothing)
    ddqn_config.obs_dim == rpo_hypr_rl_observation_dim(mdp.config) ||
        throw(DimensionMismatch("DDQN observation dimension does not match HyPR-RL"))
    ddqn_config.action_dim == action_count(mdp.config) ||
        throw(DimensionMismatch("DDQN action dimension does not match HyPR-RL"))
    learner.config.obs_dim == ddqn_config.obs_dim ||
        throw(DimensionMismatch("learner observation dimension does not match training config"))
    learner.config.action_dim == ddqn_config.action_dim ||
        throw(DimensionMismatch("learner action dimension does not match training config"))
    rng = MersenneTwister(training.seed)
    replay = MaskedReplayBuffer(
        ddqn_config.obs_dim, ddqn_config.action_dim, ddqn_config.replay_size,
    )
    episode_returns = zeros(Float64, training.episodes)
    best_objectives = fill(Inf, training.episodes)
    accepted_edits = zeros(Int, training.episodes)
    for episode in 1:training.episodes
        episode_mdp = scenario_sampler === nothing ? mdp :
            RPOHyPRRLMDP(mdp.config, scenario_sampler(rng, episode);
                         evaluator=mdp.evaluator)
        state = reset_scenario(episode_mdp, rng)
        while !state.stopped && state.edit_count < episode_mdp.config.max_edits
            observation = observe_state(episode_mdp, state)
            mask = valid_action_mask(episode_mdp, state)
            action = select_action(
                learner, observation; rng=rng, action_mask=mask,
            )
            result = step_scenario(episode_mdp, state, action, rng)
            next_observation = observe_state(episode_mdp, result.state)
            next_mask = valid_action_mask(episode_mdp, result.state)
            push!(replay, MaskedTransition(
                observation, mask, action, Float32(result.reward),
                next_observation, next_mask, result.terminated, result.truncated,
                episode,
            ))
            maybe_train!(learner, replay, rng)
            episode_returns[episode] += result.reward
            accepted_edits[episode] += result.accepted ? 1 : 0
            state = result.state
            (result.terminated || result.truncated) && break
        end
        best_objectives[episode] = state.best_evaluation.objective
        if training.checkpoint_every_episodes > 0 &&
           episode % training.checkpoint_every_episodes == 0
            checkpoint_path = joinpath(
                training.checkpoint_directory,
                "hypr_rl_episode_$(lpad(episode, 6, '0')).jls",
            )
            save_hypr_rl_checkpoint(
                checkpoint_path, learner, mdp.config;
                training_metadata=(episode=episode,),
            )
        end
    end
    return (
        learner=learner,
        replay=replay,
        episode_returns=episode_returns,
        best_objectives=best_objectives,
        accepted_edits=accepted_edits,
    )
end
