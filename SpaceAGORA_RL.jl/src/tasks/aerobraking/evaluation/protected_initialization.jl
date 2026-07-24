Base.@kwdef struct ProtectedInitializationConfig
    enabled::Bool = true
    corridor_maneuver::Bool = true
    suppress_thermal_terminal::Bool = true
    corridor_low_w_cm2::Float64 = 0.025
    corridor_high_w_cm2::Float64 = 0.40
end

function protected_initialization_config(config)
    return ProtectedInitializationConfig(
        enabled = config.protected_first_pass,
        corridor_maneuver = config.protected_initial_corridor_maneuver,
        suppress_thermal_terminal = config.protected_first_pass_suppress_thermal_terminal,
        corridor_low_w_cm2 = config.protected_corridor_low_w_cm2,
        corridor_high_w_cm2 = config.protected_corridor_high_w_cm2,
    )
end

function protected_corridor_action_index(config::ProtectedInitializationConfig,
                                         result::Union{AerobrakingStepResult,Nothing})
    config.corridor_maneuver || return nothing
    result === nothing && return nothing
    heat_rate = result.raw_observation.max_heat_rate_w_cm2
    if heat_rate > config.corridor_high_w_cm2
        return action_count()
    elseif heat_rate < config.corridor_low_w_cm2
        return 1
    end
    return nothing
end

function protected_initialization_flags(flags::TerminationFlags,
                                        config::ProtectedInitializationConfig)
    config.suppress_thermal_terminal || return flags
    thermal_only_terminal = flags.thermal_violation &&
                            flags.terminated &&
                            !(flags.success ||
                              flags.target_undershoot ||
                              flags.impact ||
                              flags.out_of_drag_passage)
    thermal_only_terminal || return flags
    return TerminationFlags(
        flags.success,
        flags.target_undershoot,
        flags.impact,
        flags.out_of_drag_passage,
        flags.thermal_violation,
        false,
        flags.truncated,
    )
end

function protected_step_result(config::AerobrakingScenarioConfig,
                               result::AerobrakingStepResult,
                               settings::ProtectedInitializationConfig)
    flags = protected_initialization_flags(result.flags, settings)
    flags === result.flags && return result
    reward = paper_reward(
        result.raw_observation,
        config,
        result.action,
        flags,
        config.reward_config,
    )
    return AerobrakingStepResult(
        result.state,
        result.action,
        result.raw_observation,
        result.normalized_observation,
        reward,
        flags,
        result.metrics,
        result.simulation_config,
    )
end

function run_protected_initializer(config::AerobrakingScenarioConfig,
                                   state::AerobrakingDecisionState,
                                   rng::AbstractRNG,
                                   summary::EpisodeSummary;
                                   settings::ProtectedInitializationConfig=
                                       ProtectedInitializationConfig())
    obs = observe_state(config, state)
    norm_obs = normalize_observation(obs, config.normalization_bounds)
    results = AerobrakingStepResult[]
    settings.enabled ||
        return (; state, observation=obs, normalized_observation=norm_obs,
                summary, done=false, results)

    result = protected_step_result(
        config,
        step_scenario(config, state, zero_action_index(), rng),
        settings,
    )
    push!(results, result)
    summary = update_episode_summary(summary, result; protected=true)
    state = result.state
    obs = result.raw_observation
    norm_obs = result.normalized_observation
    done = result.flags.terminated || result.flags.truncated

    corridor_action = done ? nothing : protected_corridor_action_index(settings, result)
    if corridor_action !== nothing
        result = protected_step_result(
            config,
            step_scenario(config, state, corridor_action, rng),
            settings,
        )
        push!(results, result)
        summary = update_episode_summary(summary, result; protected=true)
        state = result.state
        obs = result.raw_observation
        norm_obs = result.normalized_observation
        done = result.flags.terminated || result.flags.truncated
    end

    return (; state, observation=obs, normalized_observation=norm_obs,
            summary, done, results)
end
