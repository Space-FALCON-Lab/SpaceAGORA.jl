function evaluate_policy(policy::AbstractPolicy, config::AerobrakingScenarioConfig;
                         episodes::Int=1,
                         seed::Int=1,
                         policy_name::AbstractString=string(nameof(typeof(policy))))
    summaries = EpisodeSummary[]
    transitions = Transition[]
    pass_rows = NamedTuple[]
    for episode in 1:episodes
        rng = MersenneTwister(seed + episode - 1)
        state = reset_scenario(config, rng)
        obs = observe_state(config, state)
        norm_obs = normalize_observation(obs, config.normalization_bounds)
        summary = empty_episode_summary(episode_index=episode, seed=seed + episode - 1)

        while true
            action_index = policy_action_index(policy, config, state, obs, rng)
            result = step_scenario(config, state, action_index, rng)
            push!(transitions, transition_from_step(norm_obs, action_index, result, length(transitions) + 1))
            summary = update_episode_summary(summary, result)
            state = result.state
            obs = result.raw_observation
            norm_obs = result.normalized_observation
            if result.flags.terminated || result.flags.truncated
                summary = finalize_episode_summary(summary, config)
                append!(pass_rows, pass_log_rows(summary; policy_name=policy_name))
                push!(summaries, summary)
                break
            end
        end
    end
    return (summaries=summaries, transitions=transitions, pass_rows=pass_rows,
            metrics=[episode_metrics(s; policy_name=policy_name) for s in summaries],
            aggregate=aggregate_metrics(summaries; policy_name=policy_name))
end

function evaluate_baselines(config::AerobrakingScenarioConfig; episodes::Int=4, seed::Int=1)
    policies = [
        ("no_maneuver", NoManeuverPolicy()),
        ("random", RandomActionPolicy()),
        ("fixed_corridor", FixedCorridorPolicy()),
        ("aads_heuristic", AADSHeuristicPolicy()),
    ]
    return Dict(name => evaluate_policy(policy, config; episodes=episodes, seed=seed, policy_name=name)
                for (name, policy) in policies)
end
