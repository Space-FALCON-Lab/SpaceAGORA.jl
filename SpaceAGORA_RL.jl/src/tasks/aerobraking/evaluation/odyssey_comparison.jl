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
    metrics = [
        merge(
            episode_metrics(summary; policy_name=policy_name),
            episode_thermal_violation_metrics(summary, scenario),
        )
        for summary in summaries
    ]
    aggregate = merge(
        aggregate_metrics(summaries; policy_name=policy_name),
        aggregate_thermal_violation_metrics(summaries, scenario),
    )
    return (summaries=summaries, transitions=transitions, pass_rows=pass_rows,
            metrics=metrics, aggregate=aggregate)
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

function _checkpoint_validation_record(checkpoint_path::AbstractString,
                                       mode::AbstractString,
                                       aggregate::NamedTuple)
    return merge(
        (
            checkpoint = basename(checkpoint_path),
            checkpoint_path = abspath(checkpoint_path),
            mode = String(mode),
            greedy = true,
        ),
        aggregate,
    )
end

function _validation_rank(record::NamedTuple)
    success_rate = isfinite(record.success_rate) ? record.success_rate : -Inf
    thermal_failure_rate = isfinite(record.thermal_terminal_failure_rate) ?
                           record.thermal_terminal_failure_rate : Inf
    target_error = isfinite(record.mean_target_error_km) ?
                   record.mean_target_error_km : Inf
    thermal_violations = isfinite(record.mean_thermal_violations) ?
                         record.mean_thermal_violations : Inf
    return (
        success_rate,
        -thermal_failure_rate,
        -target_error,
        -thermal_violations,
    )
end

function select_best_validation_checkpoint(records::AbstractVector{<:NamedTuple};
                                           mode::AbstractString="conservative")
    candidates = filter(record -> record.mode == mode, records)
    isempty(candidates) &&
        throw(ArgumentError("no validation records found for mode $(repr(mode))"))
    best = first(candidates)
    best_rank = _validation_rank(best)
    for candidate in Iterators.drop(candidates, 1)
        rank = _validation_rank(candidate)
        if rank > best_rank
            best = candidate
            best_rank = rank
        end
    end
    return best
end

function _format_percentage_delta(value::Real)
    isfinite(value) || return "n/a"
    return @sprintf("%+.1fpp", 100 * Float64(value))
end

function _print_validation_result(record::NamedTuple,
                                  previous::Union{Nothing,NamedTuple})
    success_delta = previous === nothing ? NaN :
                    record.success_rate - previous.success_rate
    @printf(
        "validation checkpoint=%s mode=%s greedy=true episodes=%d success=%.1f%% delta_success=%s mean_target_error_km=%.2f mean_thermal_violations=%.3f\n",
        record.checkpoint,
        record.mode,
        record.episodes,
        100 * record.success_rate,
        _format_percentage_delta(success_delta),
        record.mean_target_error_km,
        record.mean_thermal_violations,
    )
    @printf(
        "validation violations checkpoint=%s mode=%s unprotected_mean_per_episode[low=%.3f soft=%.3f medium=%.3f hard=%.3f] failed_episode_rate[low=%.1f%% soft=%.1f%% medium=%.1f%% hard=%.1f%%]\n",
        record.checkpoint,
        record.mode,
        record.mean_low_thermal_violations,
        record.mean_soft_thermal_violations,
        record.mean_medium_thermal_violations,
        record.mean_hard_thermal_violations,
        100 * record.failed_with_low_violation_rate,
        100 * record.failed_with_soft_violation_rate,
        100 * record.failed_with_medium_violation_rate,
        100 * record.failed_with_hard_violation_rate,
    )
    if record.mode == "conservative"
        previous_low = previous === nothing ? NaN :
                       record.terminal_low_violation_rate -
                       previous.terminal_low_violation_rate
        previous_soft = previous === nothing ? NaN :
                        record.terminal_soft_violation_rate -
                        previous.terminal_soft_violation_rate
        previous_medium = previous === nothing ? NaN :
                          record.terminal_medium_violation_rate -
                          previous.terminal_medium_violation_rate
        previous_hard = previous === nothing ? NaN :
                        record.terminal_hard_violation_rate -
                        previous.terminal_hard_violation_rate
        @printf(
            "validation thermal_failure_attribution checkpoint=%s rate[low=%.1f%% soft=%.1f%% medium=%.1f%% hard=%.1f%%] delta[low=%s soft=%s medium=%s hard=%s]\n",
            record.checkpoint,
            100 * record.terminal_low_violation_rate,
            100 * record.terminal_soft_violation_rate,
            100 * record.terminal_medium_violation_rate,
            100 * record.terminal_hard_violation_rate,
            _format_percentage_delta(previous_low),
            _format_percentage_delta(previous_soft),
            _format_percentage_delta(previous_medium),
            _format_percentage_delta(previous_hard),
        )
    end
    flush(stdout)
    return nothing
end

function validate_frozen_checkpoints(
    checkpoint_directory::AbstractString,
    config::AerobrakingScenarioConfig;
    episodes::Int=PAPER_IID_EVALUATION_EPISODES,
    seed::Int=1,
    output_dir::AbstractString=joinpath(checkpoint_directory, "checkpoint_validation"),
    protected_initialization::ProtectedInitializationConfig=ProtectedInitializationConfig(),
    selection_mode::AbstractString="conservative",
)
    paths = frozen_checkpoint_paths(checkpoint_directory)
    isempty(paths) &&
        throw(ArgumentError("no checkpoint_*.jls files found in $(checkpoint_directory)"))
    episodes > 0 || throw(ArgumentError("validation episodes must be positive"))
    records = NamedTuple[]
    previous_by_mode = Dict{String,NamedTuple}()
    @printf(
        "checkpoint validation starting checkpoints=%d modes=%s episodes_per_mode=%d seed=%d greedy=true output_dir=%s\n",
        length(paths),
        join(PAPER_EVALUATION_MODES, ","),
        episodes,
        seed,
        output_dir,
    )
    flush(stdout)

    for checkpoint_path in paths
        modes = evaluate_frozen_checkpoint_modes(
            checkpoint_path,
            config;
            episodes=episodes,
            seed=seed,
            protected_initialization=protected_initialization,
        )
        checkpoint_name = splitext(basename(checkpoint_path))[1]
        write_evaluation_artifacts(joinpath(output_dir, checkpoint_name), modes)
        for mode in PAPER_EVALUATION_MODES
            record = _checkpoint_validation_record(
                checkpoint_path,
                mode,
                modes[mode].aggregate,
            )
            push!(records, record)
            _print_validation_result(record, get(previous_by_mode, mode, nothing))
            previous_by_mode[mode] = record
        end
    end

    best = select_best_validation_checkpoint(records; mode=selection_mode)
    artifacts = write_checkpoint_validation_artifacts(output_dir, records, best)
    @printf(
        "best validation checkpoint=%s path=%s selection_mode=%s greedy=true criterion=success_rate_then_thermal_failures_then_target_error success=%.1f%% mean_target_error_km=%.2f thermal_terminal_failure_rate=%.1f%%\n",
        best.checkpoint,
        best.checkpoint_path,
        best.mode,
        100 * best.success_rate,
        best.mean_target_error_km,
        100 * best.thermal_terminal_failure_rate,
    )
    @printf(
        "checkpoint validation complete summary=%s best=%s\n",
        artifacts.summary,
        artifacts.best,
    )
    flush(stdout)
    return (records=records, best=best, artifacts=artifacts)
end
