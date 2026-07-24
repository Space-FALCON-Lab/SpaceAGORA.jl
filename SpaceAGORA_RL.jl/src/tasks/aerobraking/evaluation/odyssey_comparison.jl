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

function evaluate_policy_modes(
    policy::AbstractPolicy,
    config::AerobrakingScenarioConfig;
    episodes::Int=PAPER_IID_EVALUATION_EPISODES,
    seed::Int=1,
    policy_name::AbstractString=string(nameof(typeof(policy))),
    protected_initialization::ProtectedInitializationConfig=ProtectedInitializationConfig(),
)
    scenarios = paper_evaluation_mode_scenarios(config)
    return Dict(
        mode => evaluate_policy(
            policy,
            scenarios[mode];
            episodes=episodes,
            seed=seed,
            policy_name="$(policy_name)_$(mode)",
            paper_protocol=false,
            protected_initialization=protected_initialization,
        )
        for mode in PAPER_EVALUATION_MODES
    )
end

function evaluate_frozen_checkpoint_modes(
    checkpoint_path::AbstractString,
    config::AerobrakingScenarioConfig;
    episodes::Int=PAPER_IID_EVALUATION_EPISODES,
    seed::Int=1,
    protected_initialization::ProtectedInitializationConfig=ProtectedInitializationConfig(),
)
    policy = load_trained_pr_drl_policy(checkpoint_path)
    checkpoint_name = splitext(basename(checkpoint_path))[1]
    return evaluate_policy_modes(
        policy,
        config;
        episodes=episodes,
        seed=seed,
        policy_name=checkpoint_name,
        protected_initialization=protected_initialization,
    )
end

function frozen_checkpoint_paths(checkpoint_directory::AbstractString)
    isdir(checkpoint_directory) ||
        throw(ArgumentError("checkpoint directory does not exist: $(checkpoint_directory)"))
    paths = filter(readdir(checkpoint_directory; join=true)) do path
        occursin(r"^checkpoint_(?:\d+|final)\.jls$", basename(path))
    end
    sort!(paths; by=path -> begin
        token = replace(splitext(basename(path))[1], "checkpoint_" => "")
        token == "final" ? typemax(Int) : parse(Int, token)
    end)
    return paths
end

function evaluate_frozen_checkpoints(
    checkpoint_paths::AbstractVector{<:AbstractString},
    config::AerobrakingScenarioConfig;
    episodes::Int=PAPER_IID_EVALUATION_EPISODES,
    seed::Int=1,
    output_dir::Union{Nothing,AbstractString}=nothing,
    protected_initialization::ProtectedInitializationConfig=ProtectedInitializationConfig(),
)
    isempty(checkpoint_paths) && throw(ArgumentError("no frozen checkpoints were provided"))
    evaluations = Dict{String,Any}()
    for checkpoint_path in checkpoint_paths
        path = String(checkpoint_path)
        modes = evaluate_frozen_checkpoint_modes(
            path,
            config;
            episodes=episodes,
            seed=seed,
            protected_initialization=protected_initialization,
        )
        evaluations[path] = modes
        if output_dir !== nothing
            checkpoint_name = splitext(basename(path))[1]
            write_evaluation_artifacts(joinpath(output_dir, checkpoint_name), modes)
        end
    end
    return evaluations
end

function evaluate_frozen_checkpoints(
    checkpoint_directory::AbstractString,
    config::AerobrakingScenarioConfig;
    kwargs...,
)
    paths = frozen_checkpoint_paths(checkpoint_directory)
    isempty(paths) &&
        throw(ArgumentError("no checkpoint_*.jls files found in $(checkpoint_directory)"))
    return evaluate_frozen_checkpoints(paths, config; kwargs...)
end
