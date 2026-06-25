function learner_policy_action!(learner::DDQNLearner, observation::Vector{Float32}, rng::AbstractRNG; test::Bool=false)
    return select_action(learner, observation; rng=rng, test=test)
end

function snapshot_policy_action(network::QNetwork, schedule::EpsilonSchedule,
                                ddqn_config::DDQNConfig,
                                observation::Vector{Float32}, step::Integer,
                                rng::AbstractRNG; test::Bool=false)
    eps = test ? 0.0 : epsilon_value(schedule, step)
    if !test && rand(rng) < eps
        return rand(rng, 1:ddqn_config.action_dim)
    end
    return argmax(predict_q(network, observation))
end

function run_worker_episode!(session::TrainingSession, episode_index::Int, worker_id::Int, seed::Int;
                             train::Bool=true)
    rng = MersenneTwister(seed)
    config = session.config.scenario
    state = reset_scenario(session.backend, rng)
    obs = observe_state(config, state)
    norm_obs = normalize_observation(obs, config.normalization_bounds)
    summary = empty_episode_summary(episode_index=episode_index, worker_id=worker_id, seed=seed)
    transitions = Transition[]

    while length(transitions) < session.config.training.max_passes_per_campaign
        action_index = learner_policy_action!(session.learner, norm_obs, rng; test=!train)
        result = step_scenario(session.backend, state, action_index, rng)
        transition = transition_from_step(norm_obs, action_index, result, length(transitions) + 1)
        push!(transitions, transition)
        if train
            observe!(session.learner, transition)
            maybe_train!(session.learner, rng)
        end
        summary = update_episode_summary(summary, result)
        state = result.state
        obs = result.raw_observation
        norm_obs = result.normalized_observation
        if result.flags.terminated || result.flags.truncated
            break
        end
    end

    return finalize_episode_summary(summary, config), transitions
end

function run_threaded_worker_episode(config::AerobrakingScenarioConfig,
                                     schedule::EpsilonSchedule,
                                     ddqn_config::DDQNConfig,
                                     policy_snapshot::QNetwork,
                                     episode_index::Int,
                                     worker_id::Int,
                                     seed::Int,
                                     max_passes_per_campaign::Int,
                                     global_step_start::Int;
                                     train::Bool=true)
    rng = MersenneTwister(seed)
    state = reset_scenario(config, rng)
    obs = observe_state(config, state)
    norm_obs = normalize_observation(obs, config.normalization_bounds)
    summary = empty_episode_summary(episode_index=episode_index, worker_id=worker_id, seed=seed)
    transitions = Transition[]

    while length(transitions) < max_passes_per_campaign
        step = global_step_start + length(transitions) + 1
        action_index = snapshot_policy_action(policy_snapshot, schedule, ddqn_config, norm_obs, step, rng; test=!train)
        result = step_scenario(config, state, action_index, rng)
        transition = transition_from_step(norm_obs, action_index, result, length(transitions) + 1)
        push!(transitions, transition)
        summary = update_episode_summary(summary, result)
        state = result.state
        obs = result.raw_observation
        norm_obs = result.normalized_observation
        if result.flags.terminated || result.flags.truncated
            break
        end
    end

    return finalize_episode_summary(summary, config), transitions
end
