function _td3_snapshot_action(actor::QNetwork, config::TD3Config,
                              observation::Vector{Float32}, step::Integer,
                              rng::AbstractRNG; test::Bool=false)
    config.action_dim == 1 || throw(ArgumentError(
        "the aerobraking TD3 policy requires action_dim=1",
    ))
    normalized = if !test && step <= config.random_steps
        2f0 * rand(rng, Float32) - 1f0
    else
        only(td3_actor_output(actor, observation))
    end
    if !test && step > config.random_steps
        normalized += Float32(config.exploration_noise) * randn(rng, Float32)
    end
    return _continuous_action_from_normalized(clamp(normalized, -1f0, 1f0))
end

function _run_td3_surrogate_episode(config::AerobrakingScenarioConfig,
                                    actor_snapshot::QNetwork,
                                    td3_config::TD3Config,
                                    episode_index::Int,
                                    worker_id::Int,
                                    action_seed::Int,
                                    scenario_seed::Int,
                                    max_passes_per_campaign::Int,
                                    global_step_start::Int,
                                    protected_initialization::ProtectedInitializationConfig)
    action_rng = MersenneTwister(action_seed)
    scenario_rng = MersenneTwister(scenario_seed)
    state = reset_scenario(config, scenario_rng)
    observation = observe_state(config, state)
    normalized_observation = normalize_observation(observation, config.normalization_bounds)
    summary = empty_episode_summary(
        episode_index=episode_index,
        worker_id=worker_id,
        seed=scenario_seed,
    )
    initial = run_protected_initializer(
        config,
        state,
        scenario_rng,
        summary;
        settings=protected_initialization,
    )
    state = initial.state
    normalized_observation = initial.normalized_observation
    summary = initial.summary
    transitions = ContinuousTransition[]
    pass_cap = min(max_passes_per_campaign, config.termination_config.max_passes)

    while !initial.done && summary.pass_count < pass_cap
        action = _td3_snapshot_action(
            actor_snapshot,
            td3_config,
            normalized_observation,
            global_step_start + length(transitions) + 1,
            action_rng,
        )
        result = step_scenario(config, state, action, scenario_rng)
        discrete_transition = transition_from_step(
            normalized_observation,
            action.index,
            result,
            length(transitions) + 1,
        )
        push!(transitions, continuous_transition(discrete_transition, result.action))
        summary = update_episode_summary(summary, result)
        state = result.state
        normalized_observation = result.normalized_observation
        if result.flags.terminated || result.flags.truncated
            break
        end
    end
    return finalize_episode_summary(summary, config), transitions
end

function _launch_spaceagora_physics_streaming_worker!(
    session::TrainingSession{<:TD3Learner},
    event_channel,
    worker_id::Int,
    episode_index::Int;
    simulation_template::Union{Nothing,SpaceAGORAPhysicsSimulationTemplate}=nothing,
    process_id::Union{Nothing,Int}=nothing,
    scenario_seed::Union{Nothing,Int}=nothing,
    successful_case_repeat_index::Int=0,
)
    seed = scenario_seed === nothing ?
           _spaceagora_physics_streaming_worker_seed(
               session.config.training.seed,
               worker_id,
               episode_index,
           ) : Int(scenario_seed)
    scenario_rng = MersenneTwister(seed)
    config = session.config.scenario
    state = reset_scenario(config, scenario_rng)
    normalized_observation = normalize_observation(
        observe_state(config, state),
        config.normalization_bounds,
    )
    protected_first_pass = session.config.training.protected_first_pass
    selected_action = protected_first_pass ? zero_action_index() :
                      select_action(
                          session.learner,
                          normalized_observation;
                          rng=session.rng,
                      )
    _, action = _resolve_policy_action(selected_action)
    summary = empty_episode_summary(
        episode_index=episode_index,
        worker_id=worker_id,
        seed=seed,
    )
    compatibility_config = DDQNConfig(
        obs_dim=session.learner.config.obs_dim,
        action_dim=action_count(),
    )
    template = nothing
    handle = if process_id === nothing
        template = simulation_template === nothing ?
                   _spaceagora_physics_simulation_template(
                       config,
                       state,
                       action;
                       campaign_max_passes=_spaceagora_physics_campaign_mission_pass_cap(
                           config,
                           session.config.training.max_passes_per_campaign,
                       ),
                   ) : simulation_template
        start_spaceagora_physics_campaign_worker!(
            event_channel,
            config,
            EpsilonSchedule(),
            compatibility_config,
            nothing,
            template,
            state,
            normalized_observation,
            selected_action,
            summary,
            episode_index,
            worker_id,
            seed,
            session.config.training.max_passes_per_campaign,
            session.learner.global_step;
            protected_first_pass=protected_first_pass,
            protected_suppress_thermal_terminal=
                session.config.training.protected_first_pass_suppress_thermal_terminal,
        )
    else
        action_channel = RemoteChannel(() -> Channel{Any}(1), process_id)
        future = remotecall(
            run_spaceagora_physics_campaign_process_worker_episode,
            process_id,
            event_channel,
            action_channel,
            config,
            EpsilonSchedule(),
            compatibility_config,
            state,
            normalized_observation,
            selected_action,
            summary,
            episode_index,
            worker_id,
            seed,
            session.config.training.max_passes_per_campaign,
            session.learner.global_step;
            protected_first_pass=protected_first_pass,
            protected_suppress_thermal_terminal=
                session.config.training.protected_first_pass_suppress_thermal_terminal,
        )
        SpaceAGORAPhysicsWorkerHandle(future, action_channel)
    end
    return SpaceAGORAPhysicsActiveWorker(
        worker_id,
        episode_index,
        seed,
        handle,
        0,
        template,
        process_id,
        successful_case_repeat_index,
    )
end

_td3_metric(value::Real) = isfinite(value) ? @sprintf("%.6g", value) : "n/a"

function _print_td3_progress(session::TrainingSession{<:TD3Learner}, summaries,
                             episode_budget::Int, target_global_step::Int,
                             active_workers::Int, start_time::Float64,
                             start_global_step::Int)
    elapsed = time() - start_time
    work_done = session.learner.global_step - start_global_step
    remaining = target_global_step == typemax(Int) ? 0 :
                max(0, target_global_step - session.learner.global_step)
    rate = elapsed > 0 ? work_done / elapsed : 0.0
    eta = rate > 0 && target_global_step != typemax(Int) ? remaining / rate : Inf
    stats = _recent_training_stats(
        summaries,
        100;
        target_tolerance_m=session.config.scenario.reward_config.target_tolerance_m,
    )
    @printf(
        "progress algo=td3 ep=%d/%s steps=%d%s replay=%d train_steps=%d actor_updates=%d actor_loss=%s critic1_loss=%s critic2_loss=%s recent_reward=%.3f recent_mean_thermal_violations=%.2f workers=%d elapsed=%s eta=%s\n",
        length(summaries),
        _budget_label(episode_budget),
        session.learner.global_step,
        target_global_step == typemax(Int) ? "" : "/$(target_global_step)",
        length(session.learner.replay),
        session.learner.train_steps,
        session.learner.actor_updates,
        _td3_metric(session.learner.last_actor_loss),
        _td3_metric(session.learner.last_critic1_loss),
        _td3_metric(session.learner.last_critic2_loss),
        stats.mean_reward,
        stats.mean_thermal_violations,
        active_workers,
        _format_duration(elapsed),
        isfinite(eta) ? _format_duration(eta) : "n/a",
    )
    flush(stdout)
    return nothing
end

function _td3_checkpoint_due!(session::TrainingSession{<:TD3Learner},
                              next_checkpoint_step::Base.RefValue{Int},
                              checkpoint_frequency::Int)
    while session.learner.global_step >= next_checkpoint_step[]
        checkpoint_path = joinpath(
            session.output_dir,
            "checkpoint_$(next_checkpoint_step[]).jls",
        )
        save_checkpoint(checkpoint_path, session.learner; manifest=session.manifest)
        @printf("checkpoint step=%d path=%s\n", next_checkpoint_step[], checkpoint_path)
        flush(stdout)
        next_checkpoint_step[] += checkpoint_frequency
    end
    return nothing
end

function _td3_training_result(session::TrainingSession{<:TD3Learner}, summaries,
                              transitions, aggregate_accumulator,
                              target_global_step::Int, start_time::Float64)
    final_checkpoint_path = joinpath(session.output_dir, "checkpoint_final.jls")
    save_checkpoint(final_checkpoint_path, session.learner; manifest=session.manifest)
    summaries isa DiskBackedHistory && close_history!(summaries)
    transitions isa DiskBackedHistory && close_history!(transitions)
    @printf(
        "training complete final_checkpoint=%s elapsed=%s\n",
        final_checkpoint_path,
        _format_duration(time() - start_time),
    )
    flush(stdout)
    metrics = summaries isa DiskBackedHistory ?
              MappedHistory(summaries, summary -> episode_metrics(summary; policy_name="td3")) :
              [episode_metrics(summary; policy_name="td3") for summary in summaries]
    aggregate = aggregate_accumulator isa EpisodeAggregateAccumulator ?
                aggregate_metrics(aggregate_accumulator; policy_name="td3") :
                aggregate_metrics(summaries; policy_name="td3")
    return (
        summaries=summaries,
        transitions=transitions,
        metrics=metrics,
        aggregate=aggregate,
        global_step=session.learner.global_step,
        target_global_step=target_global_step == typemax(Int) ? nothing : target_global_step,
        output_dir=session.output_dir,
    )
end

function _train_parallel_td3_live!(session::TrainingSession{<:TD3Learner};
                                   global_steps::Int=session.config.training.global_steps,
                                   episodes::Union{Nothing,Int}=nothing,
                                   n_workers::Int=session.config.training.n_workers,
                                   process_ids::Vector{Int}=Int[])
    summaries = DiskBackedHistory(
        EpisodeSummary,
        joinpath(session.output_dir, "training_episode_summaries.jls"),
    )
    transitions = DiskBackedHistory(
        ContinuousTransition,
        joinpath(session.output_dir, "training_transitions.jls"),
    )
    aggregate_accumulator = EpisodeAggregateAccumulator()
    requested_workers = max(1, n_workers)
    worker_backend = isempty(process_ids) ? :threads : :processes
    active_workers = worker_backend == :threads ?
                     min(requested_workers, Threads.nthreads()) :
                     min(requested_workers, length(process_ids))
    target_global_step = global_steps > 0 ? global_steps : typemax(Int)
    episode_budget = episodes === nothing ?
                     (global_steps > 0 ? typemax(Int) : session.config.training.episodes) :
                     max(0, episodes)
    checkpoint_frequency = session.config.training.checkpoint_frequency
    next_checkpoint_step = Ref(checkpoint_frequency > 0 ? checkpoint_frequency : typemax(Int))
    progress_frequency = max(0, session.config.training.progress_frequency)
    next_progress_step = Ref(progress_frequency > 0 ? progress_frequency : typemax(Int))
    start_time = time()
    start_global_step = session.learner.global_step

    @printf(
        "starting TD3 training global_steps=%s episode_cap=%s n_workers=%d active_workers=%d worker_backend=%s train_start=%d random_steps=%d batch_size=%d device=%s architecture=off_policy_pass_streaming\n",
        target_global_step == typemax(Int) ? "none" : string(target_global_step),
        _budget_label(episode_budget),
        requested_workers,
        active_workers,
        string(worker_backend),
        session.learner.config.train_start,
        session.learner.config.random_steps,
        session.learner.config.batch_size,
        training_device_name(session.learner.device),
    )
    @printf("output_dir=%s\n", session.output_dir)
    flush(stdout)

    event_channel = worker_backend == :processes ?
                    RemoteChannel(() -> Channel{Any}(max(32, 2 * active_workers)), myid()) :
                    Channel{Any}(max(32, 2 * active_workers))
    active = Dict{Int,SpaceAGORAPhysicsActiveWorker}()
    next_episode = 1
    for worker_id in 1:min(active_workers, episode_budget)
        active[worker_id] = _launch_spaceagora_physics_streaming_worker!(
            session,
            event_channel,
            worker_id,
            next_episode;
            process_id=worker_backend == :processes ? process_ids[worker_id] : nothing,
        )
        next_episode += 1
    end

    try
        while session.learner.global_step < target_global_step &&
              length(summaries) < episode_budget && !isempty(active)
            event = take!(event_channel)
            worker = get(active, event.worker_id, nothing)
            worker === nothing && continue
            event.error === nothing || throw(ErrorException(event.error))
            event.protected && (worker.protected_events_seen += 1)

            if event.transition !== nothing && !event.protected &&
               session.learner.global_step < target_global_step
                event.result === nothing && throw(ErrorException(
                    "TD3 live transition is missing its exact executed action",
                ))
                transition = continuous_transition(event.transition, event.result.action)
                push!(transitions, transition)
                session.learner.global_step += 1
                observe!(session.learner, transition)
                maybe_train!(session.learner, session.rng)
                _td3_checkpoint_due!(session, next_checkpoint_step, checkpoint_frequency)
            end

            reached_limits = session.learner.global_step >= target_global_step ||
                             length(summaries) >= episode_budget
            if event.done
                final_summary = finalize_episode_summary(event.summary, session.config.scenario)
                push!(summaries, final_summary)
                accumulate_episode!(aggregate_accumulator, final_summary)
                repeat_plan = _next_successful_case_repeat(
                    worker.seed,
                    worker.successful_case_repeat_index,
                    final_summary.success,
                    session.config.training.successful_case_repetitions,
                )
                _finish_spaceagora_physics_streaming_worker!(worker)
                delete!(active, event.worker_id)
                if !reached_limits && next_episode <= episode_budget
                    active[event.worker_id] = _launch_spaceagora_physics_streaming_worker!(
                        session,
                        event_channel,
                        event.worker_id,
                        next_episode;
                        simulation_template=worker.simulation_template,
                        process_id=worker.process_id,
                        scenario_seed=repeat_plan === nothing ? nothing : repeat_plan.seed,
                        successful_case_repeat_index=
                            repeat_plan === nothing ? 0 : repeat_plan.repeat_index,
                    )
                    next_episode += 1
                end
            elseif reached_limits
                put!(worker.handle.action_channel, nothing)
                _finish_spaceagora_physics_streaming_worker!(worker)
                delete!(active, event.worker_id)
            else
                corridor_action = event.protected && worker.protected_events_seen == 1 ?
                                  _protected_corridor_action_index(
                                      session.config.training,
                                      event.result,
                                  ) : nothing
                command = if corridor_action === nothing
                    event.transition === nothing && throw(ErrorException(
                        "TD3 live worker event is missing its next observation",
                    ))
                    select_action(
                        session.learner,
                        event.transition.next_observation;
                        rng=session.rng,
                    )
                else
                    (action_index=corridor_action, protected=true)
                end
                put!(worker.handle.action_channel, command)
            end

            if progress_frequency > 0 &&
               session.learner.global_step >= next_progress_step[]
                _print_td3_progress(
                    session,
                    summaries,
                    episode_budget,
                    target_global_step,
                    active_workers,
                    start_time,
                    start_global_step,
                )
                while next_progress_step[] <= session.learner.global_step
                    next_progress_step[] += progress_frequency
                end
            end
        end
    catch
        _stop_spaceagora_physics_streaming_workers!(active, event_channel)
        close_history!(summaries)
        close_history!(transitions)
        rethrow()
    end
    isempty(active) || _stop_spaceagora_physics_streaming_workers!(active, event_channel)
    _print_td3_progress(session, summaries, episode_budget, target_global_step,
                        active_workers, start_time, start_global_step)
    return _td3_training_result(
        session,
        summaries,
        transitions,
        aggregate_accumulator,
        target_global_step,
        start_time,
    )
end

function _train_parallel_td3_surrogate!(session::TrainingSession{<:TD3Learner};
                                        global_steps::Int=session.config.training.global_steps,
                                        episodes::Union{Nothing,Int}=nothing,
                                        n_workers::Int=session.config.training.n_workers)
    summaries = EpisodeSummary[]
    transitions = ContinuousTransition[]
    requested_workers = max(1, n_workers)
    active_workers = min(requested_workers, Threads.nthreads())
    target_global_step = global_steps > 0 ? global_steps : typemax(Int)
    episode_budget = episodes === nothing ?
                     (global_steps > 0 ? typemax(Int) : session.config.training.episodes) :
                     max(0, episodes)
    checkpoint_frequency = session.config.training.checkpoint_frequency
    next_checkpoint_step = Ref(checkpoint_frequency > 0 ? checkpoint_frequency : typemax(Int))
    progress_frequency = max(0, session.config.training.progress_frequency)
    next_progress_step = Ref(progress_frequency > 0 ? progress_frequency : typemax(Int))
    start_time = time()
    start_global_step = session.learner.global_step
    episode = 1
    repeat_seeds = Union{Nothing,Int}[nothing for _ in 1:active_workers]
    repeat_indices = zeros(Int, active_workers)

    @printf(
        "starting TD3 training global_steps=%s episode_cap=%s n_workers=%d active_workers=%d train_start=%d random_steps=%d batch_size=%d device=%s architecture=off_policy_surrogate_batches\n",
        target_global_step == typemax(Int) ? "none" : string(target_global_step),
        _budget_label(episode_budget),
        requested_workers,
        active_workers,
        session.learner.config.train_start,
        session.learner.config.random_steps,
        session.learner.config.batch_size,
        training_device_name(session.learner.device),
    )
    @printf("output_dir=%s\n", session.output_dir)
    flush(stdout)

    while session.learner.global_step < target_global_step && episode <= episode_budget
        batch_episodes = episode:min(episode_budget, episode + active_workers - 1)
        actor_snapshot = cpu_network(session.learner.actor)
        global_step_start = session.learner.global_step
        jobs = map(enumerate(batch_episodes)) do (worker_id, episode_index)
            action_seed = session.config.training.seed + 10_000 * worker_id + episode_index
            scenario_seed = repeat_indices[worker_id] == 0 ? action_seed :
                            (repeat_seeds[worker_id]::Int)
            return (
                worker_id=worker_id,
                episode_index=episode_index,
                action_seed=action_seed,
                scenario_seed=scenario_seed,
                repeat_index=repeat_indices[worker_id],
            )
        end
        tasks = map(jobs) do job
            Threads.@spawn _run_td3_surrogate_episode(
                session.config.scenario,
                $actor_snapshot,
                session.learner.config,
                $(job.episode_index),
                $(job.worker_id),
                $(job.action_seed),
                $(job.scenario_seed),
                session.config.training.max_passes_per_campaign,
                global_step_start,
                protected_initialization_config(session.config.training),
            )
        end

        for (job, task) in zip(jobs, tasks)
            summary, episode_transitions = fetch(task)
            repeat_plan = _next_successful_case_repeat(
                job.scenario_seed,
                job.repeat_index,
                summary.success,
                session.config.training.successful_case_repetitions,
            )
            repeat_seeds[job.worker_id] = repeat_plan === nothing ? nothing : repeat_plan.seed
            repeat_indices[job.worker_id] = repeat_plan === nothing ? 0 : repeat_plan.repeat_index
            remaining_steps = target_global_step - session.learner.global_step
            remaining_steps <= 0 && continue
            ingest_count = min(length(episode_transitions), remaining_steps)
            ingest_count == length(episode_transitions) && push!(summaries, summary)
            for transition in @view episode_transitions[1:ingest_count]
                push!(transitions, transition)
                session.learner.global_step += 1
                observe!(session.learner, transition)
                maybe_train!(session.learner, session.rng)
                _td3_checkpoint_due!(session, next_checkpoint_step, checkpoint_frequency)
            end
            if progress_frequency > 0 &&
               session.learner.global_step >= next_progress_step[]
                _print_td3_progress(
                    session,
                    summaries,
                    episode_budget,
                    target_global_step,
                    active_workers,
                    start_time,
                    start_global_step,
                )
                while next_progress_step[] <= session.learner.global_step
                    next_progress_step[] += progress_frequency
                end
            end
        end
        episode += length(batch_episodes)
    end
    _print_td3_progress(session, summaries, episode_budget, target_global_step,
                        active_workers, start_time, start_global_step)
    return _td3_training_result(
        session,
        summaries,
        transitions,
        summaries,
        target_global_step,
        start_time,
    )
end

function train_parallel!(session::TrainingSession{<:TD3Learner};
                         global_steps::Int=session.config.training.global_steps,
                         episodes::Union{Nothing,Int}=nothing,
                         n_workers::Int=session.config.training.n_workers)
    if _is_spaceagora_live_backend(session.config.scenario.backend_mode)
        worker_backend = session.config.training.worker_backend
        active_workers = worker_backend == :processes ? max(1, n_workers) :
                         min(max(1, n_workers), Threads.nthreads())
        return _with_spaceagora_physics_outer_parallelism(
            active_workers,
            session.config.scenario,
            worker_backend,
        ) do
            if worker_backend == :processes
                process_ids = setup_isolated_process_workers(active_workers)
                try
                    _remotecall_wait_all(
                        _prewarm_spaceagora_rl_shared_ephemeris_cache!,
                        process_ids,
                        session.config.scenario,
                        session.config.training.max_passes_per_campaign,
                    )
                    _train_parallel_td3_live!(
                        session;
                        global_steps=global_steps,
                        episodes=episodes,
                        n_workers=n_workers,
                        process_ids=process_ids,
                    )
                finally
                    rmprocs(process_ids)
                end
            else
                _prewarm_spaceagora_rl_shared_ephemeris_cache!(
                    session.config.scenario,
                    session.config.training.max_passes_per_campaign,
                )
                _train_parallel_td3_live!(
                    session;
                    global_steps=global_steps,
                    episodes=episodes,
                    n_workers=n_workers,
                )
            end
        end
    end
    return _train_parallel_td3_surrogate!(
        session;
        global_steps=global_steps,
        episodes=episodes,
        n_workers=n_workers,
    )
end
