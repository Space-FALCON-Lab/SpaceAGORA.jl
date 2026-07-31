function _evaluate_spaceagora_physics_policy(
    policy::AbstractPolicy,
    scenario::AerobrakingScenarioConfig;
    episodes::Int,
    seed::Int,
    policy_name::AbstractString,
    protected_initialization::ProtectedInitializationConfig,
)
    summaries = EpisodeSummary[]
    transitions = Transition[]
    pass_rows = NamedTuple[]
    for episode in 1:episodes
        summary, episode_transitions =
            run_spaceagora_physics_policy_campaign_episode(
                policy,
                scenario,
                episode,
                Distributed.myid(),
                seed + episode - 1;
                max_passes_per_campaign=scenario.termination_config.max_passes,
                protected_initialization=protected_initialization,
            )
        append!(transitions, episode_transitions)
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
    return (; summaries, transitions, pass_rows, metrics, aggregate)
end

function evaluate_policy(policy::AbstractPolicy, config::AerobrakingScenarioConfig;
                         episodes::Int=PAPER_IID_EVALUATION_EPISODES,
                         seed::Int=1,
                         policy_name::AbstractString=string(nameof(typeof(policy))),
                         paper_protocol::Bool=true,
                         protected_initialization::ProtectedInitializationConfig=
                             ProtectedInitializationConfig())
    scenario = paper_protocol ? paper_evaluation_scenario(config) : config
    if _is_spaceagora_live_backend(scenario.backend_mode)
        return _evaluate_spaceagora_physics_policy(
            policy,
            scenario;
            episodes=episodes,
            seed=seed,
            policy_name=policy_name,
            protected_initialization=protected_initialization,
        )
    end

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

_evaluation_std(values) = length(values) > 1 ? std(values) : 0.0

function _checkpoint_step_from_name(checkpoint_path::AbstractString)
    matched = match(r"^checkpoint_(\d+)\.jls$", basename(checkpoint_path))
    return matched === nothing ? 0 : parse(Int, matched.captures[1])
end

function _checkpoint_validation_record(checkpoint_path::AbstractString,
                                       mode::AbstractString,
                                       result,
                                       payload::AbstractDict)
    aggregate = result.aggregate
    metrics = result.metrics
    values(field::Symbol) = Float64[getfield(metric, field) for metric in metrics]
    statistic(field::Symbol) = let samples = values(field)
        (mean(samples), _evaluation_std(samples))
    end
    reward = statistic(:episode_reward)
    target_distances_km = abs.(values(:target_error_km))
    thermal = statistic(:thermal_violations)
    delta_v = statistic(:total_delta_v_mps)
    maneuvers = statistic(:maneuver_count)
    passages = statistic(:pass_count)
    duration = statistic(:mission_duration_days)
    successes = getfield.(metrics, :success)
    impacts = getfield.(metrics, :impact)
    out_of_passage = getfield.(metrics, :out_of_drag_passage)
    surpassed = [
        !successes[index] &&
        !impacts[index] &&
        !out_of_passage[index] &&
        metrics[index].target_error_km < 0
        for index in eachindex(metrics)
    ]
    unfinished = .!successes .& .!impacts .& .!out_of_passage .& .!surpassed
    episodes = length(metrics)
    global_step = Int(get(
        payload,
        :global_step,
        _checkpoint_step_from_name(checkpoint_path),
    ))
    mean_loss = Float64(get(
        payload,
        :mean_training_loss,
        get(payload, :last_loss, NaN),
    ))
    return merge(
        (
            checkpoint = basename(checkpoint_path),
            checkpoint_path = abspath(checkpoint_path),
            global_step = global_step,
            mode = String(mode),
            greedy = true,
            mean_training_loss = mean_loss,
            training_loss_sum = Float64(get(payload, :training_loss_sum, NaN)),
            training_loss_count = Int(get(payload, :training_loss_count, 0)),
        ),
        aggregate,
        (
            mean_reward = reward[1],
            std_reward = reward[2],
            success_percent = 100 * count(identity, successes) / episodes,
            impact_percent = 100 * count(identity, impacts) / episodes,
            out_of_drag_passage_percent =
                100 * count(identity, out_of_passage) / episodes,
            surpassed_target_percent = 100 * count(identity, surpassed) / episodes,
            unfinished_percent = 100 * count(identity, unfinished) / episodes,
            mean_target_error_km = mean(target_distances_km),
            std_target_error_km = _evaluation_std(target_distances_km),
            mean_thermal_violations = thermal[1],
            std_thermal_violations = thermal[2],
            mean_delta_v_mps = delta_v[1],
            std_delta_v_mps = delta_v[2],
            mean_maneuver_count = maneuvers[1],
            std_maneuver_count = maneuvers[2],
            mean_pass_count = passages[1],
            std_pass_count = passages[2],
            mean_mission_duration_days = duration[1],
            std_mission_duration_days = duration[2],
        ),
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
    checkpoint_stride::Int=1,
    write_plots::Bool=true,
)
    checkpoint_stride > 0 ||
        throw(ArgumentError("checkpoint_stride must be positive"))
    all_paths = frozen_checkpoint_paths(checkpoint_directory)
    isempty(all_paths) &&
        throw(ArgumentError("no checkpoint_*.jls files found in $(checkpoint_directory)"))
    numbered_paths = filter(path -> basename(path) != "checkpoint_final.jls", all_paths)
    paths = numbered_paths[1:checkpoint_stride:end]
    final_path = findfirst(path -> basename(path) == "checkpoint_final.jls", all_paths)
    final_path === nothing || push!(paths, all_paths[final_path])
    episodes > 0 || throw(ArgumentError("validation episodes must be positive"))
    records = NamedTuple[]
    previous_by_mode = Dict{String,NamedTuple}()
    seen_global_steps = Set{Int}()
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
        payload = load_checkpoint(checkpoint_path)
        global_step = Int(get(
            payload,
            :global_step,
            _checkpoint_step_from_name(checkpoint_path),
        ))
        global_step in seen_global_steps && continue
        push!(seen_global_steps, global_step)
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
                modes[mode],
                payload,
            )
            push!(records, record)
            _print_validation_result(record, get(previous_by_mode, mode, nothing))
            previous_by_mode[mode] = record
        end
    end

    best = select_best_validation_checkpoint(records; mode=selection_mode)
    artifacts = write_checkpoint_validation_artifacts(output_dir, records, best)
    plot_paths = write_plots ?
                 write_checkpoint_training_plots(output_dir, records) :
                 String[]
    artifacts = merge(artifacts, (plots=plot_paths,))
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
        "checkpoint validation complete summary=%s best=%s plots=%s\n",
        artifacts.summary,
        artifacts.best,
        isempty(artifacts.plots) ? "disabled" : dirname(first(artifacts.plots)),
    )
    flush(stdout)
    return (records=records, best=best, artifacts=artifacts)
end
