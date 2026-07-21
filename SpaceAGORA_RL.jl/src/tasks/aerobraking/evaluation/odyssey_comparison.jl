function paper_pr_drl_evaluation_config(; training::Bool=false,
                                        process_noise_scale::Real=1.0,
                                        backend_mode::Symbol=:paper_surrogate)
    randomization = AerobrakingRandomizationConfig(
        nominal = false,
        apoapsis_jitter_m = 2.5e3,
        periapsis_jitter_m = 1.0e3,
        nonnominal_inclination_low_deg = 88.6,
        nonnominal_inclination_high_deg = 98.6,
        nonnominal_aop_low_deg = 60.0,
        nonnominal_aop_high_deg = 90.0,
        nonnominal_raan_low_deg = 110.0,
        nonnominal_raan_high_deg = 120.0,
        process_noise = process_noise_scale > 0,
        process_noise_scale = Float64(process_noise_scale),
        aerodynamic_coefficient_dispersion = true,
        aerodynamic_coefficient_span = 0.10,
        marsgram_perturbation_scale = 1.0,
    )
    return default_aerobraking_config(phase="Main",
                                      nominal=false,
                                      max_passes=80,
                                      backend_mode=backend_mode,
                                      training=training,
                                      randomization_config=randomization)
end

function paper_odyssey_flight_evaluation_config(; training::Bool=false,
                                                process_noise_scale::Real=1.0,
                                                backend_mode::Symbol=:paper_surrogate)
    randomization = AerobrakingRandomizationConfig(
        nominal = true,
        apoapsis_jitter_m = 2.5e3,
        periapsis_jitter_m = 1.0e3,
        angle_jitter_deg = 0.0,
        process_noise = process_noise_scale > 0,
        process_noise_scale = Float64(process_noise_scale),
        aerodynamic_coefficient_dispersion = true,
        aerodynamic_coefficient_span = 0.10,
        marsgram_perturbation_scale = 1.0,
    )
    return default_aerobraking_config(phase="Main",
                                      nominal=true,
                                      max_passes=80,
                                      backend_mode=backend_mode,
                                      training=training,
                                      randomization_config=randomization,
                                      nominal_inclination_deg=93.6,
                                      nominal_argument_of_periapsis_deg=115.0,
                                      nominal_raan_deg=89.0)
end

function paper_pr_drl_marsgram_evaluation_config(; training::Bool=false,
                                                 process_noise_scale::Real=1.0)
    return paper_pr_drl_evaluation_config(training=training,
                                          process_noise_scale=process_noise_scale,
                                          backend_mode=:spaceagora_marsgram)
end

function paper_pr_drl_physics_evaluation_config(; training::Bool=false,
                                                process_noise_scale::Real=1.0)
    return paper_pr_drl_evaluation_config(training=training,
                                          process_noise_scale=process_noise_scale,
                                          backend_mode=:spaceagora_physics)
end

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
            transition_action = action_index isa AerobrakingAction ?
                                nearest_action_index(action_index.delta_v_mps) :
                                Int(action_index)
            push!(transitions, transition_from_step(norm_obs, transition_action, result, length(transitions) + 1))
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
