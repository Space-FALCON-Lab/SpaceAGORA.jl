function evaluate_policy(policy::AbstractPolicy, config::AerobrakingScenarioConfig;
                         episodes::Int=PAPER_IID_EVALUATION_EPISODES,
                         seed::Int=1,
                         policy_name::AbstractString=string(nameof(typeof(policy))),
                         paper_protocol::Bool=true,
                         protected_initialization::ProtectedInitializationConfig=
                             ProtectedInitializationConfig())
    scenario = paper_protocol ? paper_evaluation_scenario(config) : config
    summaries = EpisodeSummary[]
    transitions = Transition[]
    pass_rows = NamedTuple[]
    for episode in 1:episodes
        rng = MersenneTwister(seed + episode - 1)
        state = reset_scenario(scenario, rng)
        summary = empty_episode_summary(episode_index=episode, seed=seed + episode - 1)
        initial = run_protected_initializer(
            scenario,
            state,
            rng,
            summary;
            settings=protected_initialization,
        )
        state = initial.state
        obs = initial.observation
        norm_obs = initial.normalized_observation
        summary = initial.summary

        while !initial.done
            action_index = policy_action_index(policy, scenario, state, obs, rng)
            result = step_scenario(scenario, state, action_index, rng)
            push!(transitions, transition_from_step(norm_obs, action_index, result, length(transitions) + 1))
            summary = update_episode_summary(summary, result)
            state = result.state
            obs = result.raw_observation
            norm_obs = result.normalized_observation
            if result.flags.terminated || result.flags.truncated
                break
            end
        end
        summary = finalize_episode_summary(summary, scenario)
        append!(pass_rows, pass_log_rows(summary; policy_name=policy_name))
        push!(summaries, summary)
    end
    return (summaries=summaries, transitions=transitions, pass_rows=pass_rows,
            metrics=[episode_metrics(s; policy_name=policy_name) for s in summaries],
            aggregate=aggregate_metrics(summaries; policy_name=policy_name))
end

function evaluate_baselines(config::AerobrakingScenarioConfig;
                            episodes::Int=PAPER_IID_EVALUATION_EPISODES,
                            seed::Int=1,
                            paper_protocol::Bool=true,
                            protected_initialization::ProtectedInitializationConfig=
                                ProtectedInitializationConfig())
    policies = [
        ("no_maneuver", NoManeuverPolicy()),
        ("random", RandomActionPolicy()),
        ("fixed_corridor", FixedCorridorPolicy()),
        ("aads_heuristic", AADSHeuristicPolicy()),
    ]
    return Dict(name => evaluate_policy(
                    policy,
                    config;
                    episodes=episodes,
                    seed=seed,
                    policy_name=name,
                    paper_protocol=paper_protocol,
                    protected_initialization=protected_initialization,
                )
                for (name, policy) in policies)
end
