function _format_duration(seconds::Real)
    total = max(0, round(Int, Float64(seconds)))
    hours = total ÷ 3600
    minutes = (total % 3600) ÷ 60
    secs = total % 60
    return @sprintf("%02d:%02d:%02d", hours, minutes, secs)
end

function _recent_training_stats(summaries::AbstractVector{<:EpisodeSummary}, window::Int;
                                target_tolerance_m::Real=10e3)
    isempty(summaries) && return (
        mean_reward=NaN,
        reached_goal_percent=NaN,
        mean_thermal_violations=NaN,
        mean_passes_to_end=NaN,
        mean_end_distance_km=NaN,
    )
    first_idx = max(1, length(summaries) - max(1, window) + 1)
    recent = summaries[first_idx:end]
    thermal_free = filter(summary -> summary.thermal_violations == 0, recent)
    final_errors_m = getfield.(thermal_free, :target_error_m)
    finite_final_errors_m = filter(isfinite, final_errors_m)
    return (
        mean_reward = mean(getfield.(recent, :episode_reward)),
        reached_goal_percent = isempty(thermal_free) ?
                               NaN :
                               100 *
            count(error_m -> isfinite(error_m) &&
                             abs(error_m) <= target_tolerance_m, final_errors_m) /
            length(thermal_free),
        mean_thermal_violations = mean(getfield.(recent, :thermal_violations)),
        mean_passes_to_end = mean(getfield.(recent, :pass_count)),
        mean_end_distance_km = isempty(finite_final_errors_m) ?
                               NaN :
                               mean(abs, finite_final_errors_m) / 1000,
    )
end

function _budget_label(budget::Int)
    return budget == typemax(Int) ? "none" : string(budget)
end

function _print_training_progress(session::TrainingSession, summaries::AbstractVector{<:EpisodeSummary},
                                  episode_budget::Int, target_global_step::Int,
                                  active_workers::Int, start_time::Float64,
                                  start_global_step::Int)
    completed_episodes = length(summaries)
    elapsed = time() - start_time
    step_limited = target_global_step != typemax(Int)
    work_done = step_limited ? session.learner.global_step - start_global_step : completed_episodes
    work_remaining = if step_limited
        max(0, target_global_step - session.learner.global_step)
    elseif episode_budget == typemax(Int)
        0
    else
        max(0, episode_budget - completed_episodes)
    end
    work_rate = elapsed > 0 ? work_done / elapsed : 0.0
    eta = work_rate > 0 && (step_limited || episode_budget != typemax(Int)) ? work_remaining / work_rate : Inf
    stats = _recent_training_stats(
        summaries,
        100;
        target_tolerance_m=session.config.scenario.reward_config.target_tolerance_m,
    )
    eps = epsilon_value(session.learner.schedule, session.learner.global_step)
    loss = isfinite(session.learner.last_loss) ? @sprintf("%.6g", session.learner.last_loss) : "n/a"
    if step_limited
        @printf(
            "progress ep=%d steps=%d/%d replay=%d train_steps=%d loss=%s eps=%.4f recent_reward=%.3f recent_mean_thermal_violations=%.2f recent_mean_passes_to_end=%.1f recent_reached_goal_no_thermal=%.1f%% recent_mean_end_distance_no_thermal_km=%.2f workers=%d/%d elapsed=%s eta=%s\n",
            completed_episodes,
            session.learner.global_step,
            target_global_step,
            length(session.learner.replay),
            session.learner.train_steps,
            loss,
            eps,
            stats.mean_reward,
            stats.mean_thermal_violations,
            stats.mean_passes_to_end,
            stats.reached_goal_percent,
            stats.mean_end_distance_km,
            active_workers,
            Threads.nthreads(),
            _format_duration(elapsed),
            isfinite(eta) ? _format_duration(eta) : "n/a",
        )
    else
        @printf(
            "progress ep=%d/%s steps=%d replay=%d train_steps=%d loss=%s eps=%.4f recent_reward=%.3f recent_mean_thermal_violations=%.2f recent_mean_passes_to_end=%.1f recent_reached_goal_no_thermal=%.1f%% recent_mean_end_distance_no_thermal_km=%.2f workers=%d/%d elapsed=%s eta=%s\n",
            completed_episodes,
            _budget_label(episode_budget),
            session.learner.global_step,
            length(session.learner.replay),
            session.learner.train_steps,
            loss,
            eps,
            stats.mean_reward,
            stats.mean_thermal_violations,
            stats.mean_passes_to_end,
            stats.reached_goal_percent,
            stats.mean_end_distance_km,
            active_workers,
            Threads.nthreads(),
            _format_duration(elapsed),
            isfinite(eta) ? _format_duration(eta) : "n/a",
        )
    end
    flush(stdout)
end

mutable struct SpaceAGORAPhysicsActiveWorker
    worker_id::Int
    episode_index::Int
    seed::Int
    handle::SpaceAGORAPhysicsWorkerHandle
    protected_events_seen::Int
    simulation_template::Union{Nothing,SpaceAGORAPhysicsSimulationTemplate}
    process_id::Union{Nothing,Int}
    successful_case_repeat_index::Int
end

function _next_successful_case_repeat(seed::Int, repeat_index::Int, success::Bool,
                                      max_repetitions::Int)
    max_repetitions <= 0 && return nothing
    if repeat_index == 0
        return success ? (seed=seed, repeat_index=1) : nothing
    end
    repeat_index < max_repetitions || return nothing
    return (seed=seed, repeat_index=repeat_index + 1)
end

function _spaceagora_physics_streaming_worker_seed(base_seed::Int, worker_id::Int,
                                                   episode_index::Int)
    return base_seed + 10_000 * worker_id + episode_index
end

function _ddqn_master_action(learner::DDQNLearner, observation::Vector{Float32},
                             rng::AbstractRNG; test::Bool=false)
    step = learner.global_step + 1
    eps = test ? 0.0 : epsilon_value(learner.schedule, step)
    if !test && rand(rng) < eps
        return rand(rng, 1:learner.config.action_dim)
    end
    return greedy_action_index(learner, observation)
end

function _protected_corridor_action_index(training::TrainingConfig,
                                          result::Union{AerobrakingStepResult,Nothing})
    return protected_corridor_action_index(
        protected_initialization_config(training),
        result,
    )
end

function _launch_spaceagora_physics_streaming_worker!(session::TrainingSession{<:DDQNLearner},
                                                      event_channel,
                                                      worker_id::Int,
                                                      episode_index::Int;
                                                      simulation_template::Union{Nothing,SpaceAGORAPhysicsSimulationTemplate}=nothing,
                                                      process_id::Union{Nothing,Int}=nothing,
                                                      scenario_seed::Union{Nothing,Int}=nothing,
                                                      successful_case_repeat_index::Int=0)
    seed = scenario_seed === nothing ?
           _spaceagora_physics_streaming_worker_seed(session.config.training.seed,
                                                     worker_id,
                                                     episode_index) :
           Int(scenario_seed)
    rng = MersenneTwister(seed)
    config = session.config.scenario
    state = reset_scenario(config, rng)
    norm_obs = normalize_observation(observe_state(config, state), config.normalization_bounds)
    protected_first_pass = session.config.training.protected_first_pass
    action_index = protected_first_pass ?
                   zero_action_index() :
                   _ddqn_master_action(session.learner, norm_obs, session.rng)
    action = action_from_index(action_index)
    summary = empty_episode_summary(episode_index=episode_index, worker_id=worker_id, seed=seed)
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
                   ) :
                   simulation_template
        start_spaceagora_physics_campaign_worker!(
            event_channel,
            config,
            session.learner.schedule,
            session.learner.config,
            nothing,
            template,
            state,
            norm_obs,
            action_index,
            summary,
            episode_index,
            worker_id,
            seed,
            session.config.training.max_passes_per_campaign,
            session.learner.global_step,
            protected_first_pass = protected_first_pass,
            protected_suppress_thermal_terminal =
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
            session.learner.schedule,
            session.learner.config,
            state,
            norm_obs,
            action_index,
            summary,
            episode_index,
            worker_id,
            seed,
            session.config.training.max_passes_per_campaign,
            session.learner.global_step;
            protected_first_pass = protected_first_pass,
            protected_suppress_thermal_terminal =
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

function _launch_spaceagora_physics_streaming_worker!(
    session::TrainingSession{<:A2CLearner},
    event_channel,
    worker_id::Int,
    episode_index::Int,
    actor_snapshot::QNetwork,
    action_rng::AbstractRNG,
    policy_version::Int;
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
           ) :
           Int(scenario_seed)
    scenario_rng = MersenneTwister(seed)
    config = session.config.scenario
    state = reset_scenario(config, scenario_rng)
    norm_obs = normalize_observation(observe_state(config, state), config.normalization_bounds)
    protected_first_pass = session.config.training.protected_first_pass
    action_index = protected_first_pass ?
                   zero_action_index() :
                   actor_action(actor_snapshot, norm_obs, action_rng; test=false)
    action = action_from_index(action_index)
    summary = empty_episode_summary(episode_index=episode_index, worker_id=worker_id, seed=seed)
    compatibility_config = DDQNConfig(
        obs_dim=session.learner.config.obs_dim,
        action_dim=session.learner.config.action_dim,
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
                   ) :
                   simulation_template
        start_spaceagora_physics_campaign_worker!(
            event_channel,
            config,
            EpsilonSchedule(),
            compatibility_config,
            nothing,
            template,
            state,
            norm_obs,
            action_index,
            summary,
            episode_index,
            worker_id,
            seed,
            session.config.training.max_passes_per_campaign,
            session.learner.global_step;
            protected_first_pass=protected_first_pass,
            protected_suppress_thermal_terminal=
                session.config.training.protected_first_pass_suppress_thermal_terminal,
            policy_version=policy_version,
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
            norm_obs,
            action_index,
            summary,
            episode_index,
            worker_id,
            seed,
            session.config.training.max_passes_per_campaign,
            session.learner.global_step;
            protected_first_pass=protected_first_pass,
            protected_suppress_thermal_terminal=
                session.config.training.protected_first_pass_suppress_thermal_terminal,
            policy_version=policy_version,
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

function _finish_spaceagora_physics_streaming_worker!(worker::SpaceAGORAPhysicsActiveWorker)
    fetch(worker.handle.task)
    return nothing
end

function _stop_spaceagora_physics_streaming_workers!(active::Dict{Int,SpaceAGORAPhysicsActiveWorker},
                                                     event_channel)
    while !isempty(active)
        event = take!(event_channel)
        worker = get(active, event.worker_id, nothing)
        worker === nothing && continue
        event.error === nothing || @warn "SpaceAGORA physics worker stopped with an error during shutdown" worker_id=event.worker_id error=event.error
        if !event.done
            put!(worker.handle.action_channel, nothing)
        end
        _finish_spaceagora_physics_streaming_worker!(worker)
        delete!(active, event.worker_id)
    end
    return nothing
end

function _with_spaceagora_physics_outer_parallelism(f::Function,
                                                    active_workers::Int,
                                                    config::AerobrakingScenarioConfig,
                                                    worker_backend::Symbol)
    env_pairs = _spaceagora_rl_core_ephemeris_env_pairs()
    push!(
        env_pairs,
        "SPACEAGORA_GRAM_ONCE_PER_STEP" =>
            (config.spaceagora_gram_once_per_step ? "1" : "0"),
    )
    if active_workers > 1
        push!(env_pairs, "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1")
    end
    if active_workers > 1 && isempty(strip(get(ENV, "SPACEAGORA_INNER_THREAD_BUDGET", "")))
        inner_budget = worker_backend == :processes ?
                       1 :
                       max(1, fld(Threads.nthreads(), active_workers))
        push!(env_pairs, "SPACEAGORA_INNER_THREAD_BUDGET" => string(inner_budget))
    end
    isempty(env_pairs) && return f()
    return withenv(env_pairs...) do
        f()
    end
end

function _with_spaceagora_physics_outer_parallelism(f::Function, active_workers::Int)
    return _with_spaceagora_physics_outer_parallelism(
        f,
        active_workers,
        default_aerobraking_config(backend_mode=:spaceagora_physics),
        :threads,
    )
end

function _train_parallel_spaceagora_physics_streaming!(session::TrainingSession{<:DDQNLearner};
                                                       global_steps::Int=session.config.training.global_steps,
                                                       episodes::Union{Nothing,Int}=nothing,
                                                       n_workers::Int=session.config.training.n_workers,
                                                       process_ids::Vector{Int}=Int[])
    summaries = DiskBackedHistory(
        EpisodeSummary,
        joinpath(session.output_dir, "training_episode_summaries.jls"),
    )
    transitions = DiskBackedHistory(
        Transition,
        joinpath(session.output_dir, "training_transitions.jls"),
    )
    aggregate_accumulator = EpisodeAggregateAccumulator()
    requested_workers = max(1, n_workers)
    worker_backend = isempty(process_ids) ? :threads : :processes
    active_workers = worker_backend == :threads ?
                     min(requested_workers, Threads.nthreads()) :
                     min(requested_workers, length(process_ids))
    if worker_backend == :threads && active_workers < requested_workers
        @warn "n_workers exceeds available Julia threads; using Threads.nthreads()" requested_workers active_workers
    end
    target_global_step = global_steps > 0 ? global_steps : typemax(Int)
    episode_budget = episodes === nothing ? (global_steps > 0 ? typemax(Int) : session.config.training.episodes) : max(0, episodes)
    checkpoint_frequency = session.config.training.checkpoint_frequency
    next_checkpoint_step = checkpoint_frequency > 0 ? checkpoint_frequency : typemax(Int)
    progress_frequency = max(0, session.config.training.progress_frequency)
    next_progress_step = progress_frequency > 0 && target_global_step != typemax(Int) ?
                         min(target_global_step, session.learner.global_step + progress_frequency) :
                         typemax(Int)
    next_progress_episode = progress_frequency > 0 && target_global_step == typemax(Int) ?
                            min(progress_frequency, episode_budget) :
                            typemax(Int)
    start_time = time()
    start_global_step = session.learner.global_step
    last_progress_step = -1
    last_progress_episode = -1
    printed_initial_progress = false

    @printf(
        "starting %s training global_steps=%s episode_cap=%s n_workers=%d active_workers=%d worker_backend=%s julia_threads=%d outer_parallel=%s inner_thread_budget=%s train_start=%d batch_size=%d checkpoint_frequency=%d progress_frequency=%d successful_case_repetitions=%d architecture=paper_pass_streaming\n",
        algorithm_display_name(session.config.training.algorithm),
        target_global_step == typemax(Int) ? "none" : string(target_global_step),
        _budget_label(episode_budget),
        requested_workers,
        active_workers,
        string(worker_backend),
        Threads.nthreads(),
        get(ENV, "SPACEAGORA_OUTER_PARALLEL_ACTIVE", "0"),
        get(ENV, "SPACEAGORA_INNER_THREAD_BUDGET", "auto"),
        session.learner.config.train_start,
        session.learner.config.batch_size,
        checkpoint_frequency,
        progress_frequency,
        session.config.training.successful_case_repetitions,
    )
    @printf("output_dir=%s\n", session.output_dir)
    @printf(
        "epsilon start=%.4f stop=%.4f decay_steps=%d decay_start_step=%d\n",
        session.learner.schedule.start,
        session.learner.schedule.stop,
        session.learner.schedule.decay_steps,
        session.learner.schedule.decay_start_step,
    )
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
            next_episode,
            process_id=worker_backend == :processes ? process_ids[worker_id] : nothing,
        )
        next_episode += 1
    end

    while session.learner.global_step < target_global_step &&
          length(summaries) < episode_budget &&
          !isempty(active)
        event = take!(event_channel)
        worker = get(active, event.worker_id, nothing)
        worker === nothing && continue

        if event.error !== nothing
            delete!(active, event.worker_id)
            _finish_spaceagora_physics_streaming_worker!(worker)
            _stop_spaceagora_physics_streaming_workers!(active, event_channel)
            close_history!(summaries)
            close_history!(transitions)
            throw(ErrorException(event.error))
        end

        if event.protected
            worker.protected_events_seen += 1
        end

        if event.transition !== nothing &&
           !event.protected &&
           session.learner.global_step < target_global_step
            transition = event.transition
            push!(transitions, transition)
            session.learner.global_step += 1
            observe!(session.learner, transition)
            maybe_train!(session.learner, session.rng)
            while session.learner.global_step >= next_checkpoint_step
                checkpoint_path = joinpath(session.output_dir, "checkpoint_$(next_checkpoint_step).jls")
                save_checkpoint(checkpoint_path, session.learner; manifest=session.manifest)
                @printf("checkpoint step=%d path=%s\n", next_checkpoint_step, checkpoint_path)
                flush(stdout)
                next_checkpoint_step += checkpoint_frequency
            end
        end

        reached_limits = session.learner.global_step >= target_global_step ||
                         length(summaries) >= episode_budget
        if event.done
            repeat_plan = nothing
            if length(summaries) < episode_budget
                final_summary = finalize_episode_summary(event.summary, session.config.scenario)
                push!(summaries, final_summary)
                accumulate_episode!(aggregate_accumulator, final_summary)
                repeat_plan = _next_successful_case_repeat(
                    worker.seed,
                    worker.successful_case_repeat_index,
                    final_summary.success,
                    session.config.training.successful_case_repetitions,
                )
            end
            _finish_spaceagora_physics_streaming_worker!(worker)
            delete!(active, event.worker_id)
            if !reached_limits && next_episode <= episode_budget
                active[event.worker_id] = _launch_spaceagora_physics_streaming_worker!(
                    session,
                    event_channel,
                    event.worker_id,
                    next_episode,
                    simulation_template=worker.simulation_template,
                    process_id=worker.process_id,
                    scenario_seed=repeat_plan === nothing ? nothing : repeat_plan.seed,
                    successful_case_repeat_index=
                        repeat_plan === nothing ? 0 : repeat_plan.repeat_index,
                )
                if repeat_plan !== nothing
                    @printf(
                        "successful_case_repeat worker=%d scenario_seed=%d repeat=%d/%d episode=%d\n",
                        event.worker_id,
                        repeat_plan.seed,
                        repeat_plan.repeat_index,
                        session.config.training.successful_case_repetitions,
                        next_episode,
                    )
                    flush(stdout)
                end
                next_episode += 1
            end
        else
            if reached_limits
                put!(worker.handle.action_channel, nothing)
                _finish_spaceagora_physics_streaming_worker!(worker)
                delete!(active, event.worker_id)
            else
                corridor_action = event.protected && worker.protected_events_seen == 1 ?
                                  _protected_corridor_action_index(session.config.training, event.result) :
                                  nothing
                if corridor_action === nothing
                    next_action = _ddqn_master_action(
                        session.learner,
                        (event.transition::Transition).next_observation,
                        session.rng,
                    )
                    put!(worker.handle.action_channel, next_action)
                else
                    put!(worker.handle.action_channel, (action_index=corridor_action, protected=true))
                end
            end
        end

        if progress_frequency > 0 && !printed_initial_progress &&
           session.learner.global_step > start_global_step
            _print_training_progress(session, summaries, episode_budget,
                                     target_global_step, active_workers,
                                     start_time, start_global_step)
            printed_initial_progress = true
            if target_global_step != typemax(Int)
                last_progress_step = session.learner.global_step
                while next_progress_step <= session.learner.global_step &&
                      next_progress_step < target_global_step
                    next_progress_step += progress_frequency
                end
                next_progress_step = min(next_progress_step, target_global_step)
            else
                last_progress_episode = length(summaries)
                while next_progress_episode <= length(summaries)
                    next_progress_episode += progress_frequency
                end
            end
        elseif progress_frequency > 0 && target_global_step != typemax(Int) &&
               session.learner.global_step >= next_progress_step
            _print_training_progress(session, summaries, episode_budget,
                                     target_global_step, active_workers,
                                     start_time, start_global_step)
            last_progress_step = session.learner.global_step
            while next_progress_step <= session.learner.global_step &&
                  next_progress_step < target_global_step
                next_progress_step += progress_frequency
            end
            next_progress_step = min(next_progress_step, target_global_step)
        elseif progress_frequency > 0 && target_global_step == typemax(Int) &&
               length(summaries) >= next_progress_episode
            _print_training_progress(session, summaries, episode_budget,
                                     target_global_step, active_workers,
                                     start_time, start_global_step)
            last_progress_episode = length(summaries)
            while next_progress_episode <= length(summaries)
                next_progress_episode += progress_frequency
            end
        end
    end

    if !isempty(active)
        _stop_spaceagora_physics_streaming_workers!(active, event_channel)
    end

    if target_global_step != typemax(Int)
        if last_progress_step != session.learner.global_step
            _print_training_progress(session, summaries, episode_budget,
                                     target_global_step, active_workers,
                                     start_time, start_global_step)
        end
    elseif last_progress_episode != length(summaries)
        _print_training_progress(session, summaries, episode_budget,
                                 target_global_step, active_workers,
                                 start_time, start_global_step)
    end

    final_checkpoint_path = joinpath(session.output_dir, "checkpoint_final.jls")
    save_checkpoint(final_checkpoint_path, session.learner; manifest=session.manifest)
    close_history!(summaries)
    close_history!(transitions)
    @printf("training complete final_checkpoint=%s elapsed=%s\n",
            final_checkpoint_path, _format_duration(time() - start_time))
    flush(stdout)
    policy_name = algorithm_report_name(session.config.training.algorithm)
    return (summaries=summaries, transitions=transitions,
            metrics=MappedHistory(summaries, summary -> episode_metrics(summary; policy_name=policy_name)),
            aggregate=aggregate_metrics(aggregate_accumulator; policy_name=policy_name),
            global_step=session.learner.global_step,
            target_global_step=target_global_step == typemax(Int) ? nothing : target_global_step,
            output_dir=session.output_dir)
end

function train_parallel!(session::TrainingSession{<:DDQNLearner};
                         global_steps::Int=session.config.training.global_steps,
                         episodes::Union{Nothing,Int}=nothing,
                         n_workers::Int=session.config.training.n_workers)
    if session.config.training.algorithm == :pr_drl &&
       _is_spaceagora_live_backend(session.config.scenario.backend_mode)
        worker_backend = session.config.training.worker_backend
        active_workers = worker_backend == :processes ?
                         max(1, n_workers) :
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
                    return _train_parallel_spaceagora_physics_streaming!(
                        session;
                        global_steps = global_steps,
                        episodes = episodes,
                        n_workers = n_workers,
                        process_ids = process_ids,
                    )
                finally
                    rmprocs(process_ids)
                end
            else
                _prewarm_spaceagora_rl_shared_ephemeris_cache!(
                    session.config.scenario,
                    session.config.training.max_passes_per_campaign,
                )
                return _train_parallel_spaceagora_physics_streaming!(
                    session;
                    global_steps = global_steps,
                    episodes = episodes,
                    n_workers = n_workers,
                )
            end
        end
    end

    summaries = EpisodeSummary[]
    transitions = Transition[]
    requested_workers = max(1, n_workers)
    active_workers = min(requested_workers, Threads.nthreads())
    active_workers < requested_workers && @warn "n_workers exceeds available Julia threads; using Threads.nthreads()" requested_workers active_workers
    target_global_step = global_steps > 0 ? global_steps : typemax(Int)
    episode_budget = episodes === nothing ? (global_steps > 0 ? typemax(Int) : session.config.training.episodes) : max(0, episodes)
    checkpoint_frequency = session.config.training.checkpoint_frequency
    next_checkpoint_step = checkpoint_frequency > 0 ? checkpoint_frequency : typemax(Int)
    progress_frequency = max(0, session.config.training.progress_frequency)
    next_progress_step = progress_frequency > 0 && target_global_step != typemax(Int) ?
                         min(target_global_step, session.learner.global_step + progress_frequency) :
                         typemax(Int)
    next_progress_episode = progress_frequency > 0 && target_global_step == typemax(Int) ?
                            min(progress_frequency, episode_budget) :
                            typemax(Int)
    start_time = time()
    start_global_step = session.learner.global_step
    last_progress_step = -1
    last_progress_episode = -1
    printed_initial_progress = false

    @printf(
        "starting %s training global_steps=%s episode_cap=%s n_workers=%d active_workers=%d julia_threads=%d train_start=%d batch_size=%d checkpoint_frequency=%d progress_frequency=%d successful_case_repetitions=%d\n",
        algorithm_display_name(session.config.training.algorithm),
        target_global_step == typemax(Int) ? "none" : string(target_global_step),
        _budget_label(episode_budget),
        requested_workers,
        active_workers,
        Threads.nthreads(),
        session.learner.config.train_start,
        session.learner.config.batch_size,
        checkpoint_frequency,
        progress_frequency,
        session.config.training.successful_case_repetitions,
    )
    @printf("output_dir=%s\n", session.output_dir)
    @printf(
        "epsilon start=%.4f stop=%.4f decay_steps=%d decay_start_step=%d\n",
        session.learner.schedule.start,
        session.learner.schedule.stop,
        session.learner.schedule.decay_steps,
        session.learner.schedule.decay_start_step,
    )
    flush(stdout)

    episode = 1
    repeat_seeds = Union{Nothing,Int}[nothing for _ in 1:active_workers]
    repeat_indices = zeros(Int, active_workers)
    while session.learner.global_step < target_global_step && episode <= episode_budget
        batch_episodes = episode:min(episode_budget, episode + active_workers - 1)
        policy_snapshot = cpu_network(session.learner.online)
        global_step_start = session.learner.global_step
        jobs = map(enumerate(batch_episodes)) do (local_worker_id, episode_index)
            action_seed = session.config.training.seed + 10_000 * local_worker_id + episode_index
            scenario_seed = repeat_indices[local_worker_id] == 0 ?
                            action_seed :
                            (repeat_seeds[local_worker_id]::Int)
            return (
                worker_id=local_worker_id,
                episode_index=episode_index,
                action_seed=action_seed,
                scenario_seed=scenario_seed,
                repeat_index=repeat_indices[local_worker_id],
            )
        end
        tasks = map(jobs) do job
            Threads.@spawn run_threaded_worker_episode(
                session.config.scenario,
                session.learner.schedule,
                session.learner.config,
                $policy_snapshot,
                $(job.episode_index),
                $(job.worker_id),
                $(job.action_seed),
                session.config.training.max_passes_per_campaign,
                global_step_start;
                train=true,
                scenario_seed=$(job.scenario_seed),
                protected_initialization=
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
            if repeat_plan !== nothing
                @printf(
                    "successful_case_repeat worker=%d scenario_seed=%d next_repeat=%d/%d after_episode=%d\n",
                    job.worker_id,
                    repeat_plan.seed,
                    repeat_plan.repeat_index,
                    session.config.training.successful_case_repetitions,
                    job.episode_index,
                )
                flush(stdout)
            end
            remaining_steps = target_global_step - session.learner.global_step
            remaining_steps <= 0 && continue
            ingest_count = min(length(episode_transitions), remaining_steps)
            if ingest_count == length(episode_transitions)
                push!(summaries, summary)
            end
            for idx in 1:ingest_count
                transition = episode_transitions[idx]
                push!(transitions, transition)
                session.learner.global_step += 1
                observe!(session.learner, transition)
                maybe_train!(session.learner, session.rng)
                while session.learner.global_step >= next_checkpoint_step
                    checkpoint_path = joinpath(session.output_dir, "checkpoint_$(next_checkpoint_step).jls")
                    save_checkpoint(checkpoint_path, session.learner; manifest=session.manifest)
                    @printf("checkpoint step=%d path=%s\n", next_checkpoint_step, checkpoint_path)
                    flush(stdout)
                    next_checkpoint_step += checkpoint_frequency
                end
            end
            if progress_frequency > 0 && !printed_initial_progress &&
               session.learner.global_step > start_global_step
                _print_training_progress(session, summaries, episode_budget,
                                         target_global_step, active_workers,
                                         start_time, start_global_step)
                printed_initial_progress = true
                if target_global_step != typemax(Int)
                    last_progress_step = session.learner.global_step
                    while next_progress_step <= session.learner.global_step &&
                          next_progress_step < target_global_step
                        next_progress_step += progress_frequency
                    end
                    next_progress_step = min(next_progress_step, target_global_step)
                else
                    last_progress_episode = length(summaries)
                    while next_progress_episode <= length(summaries)
                        next_progress_episode += progress_frequency
                    end
                end
            elseif progress_frequency > 0 && target_global_step != typemax(Int) &&
                   session.learner.global_step >= next_progress_step
                _print_training_progress(session, summaries, episode_budget,
                                         target_global_step, active_workers,
                                         start_time, start_global_step)
                last_progress_step = session.learner.global_step
                while next_progress_step <= session.learner.global_step &&
                      next_progress_step < target_global_step
                    next_progress_step += progress_frequency
                end
                next_progress_step = min(next_progress_step, target_global_step)
            end
        end

        if progress_frequency > 0 && target_global_step == typemax(Int) &&
           length(summaries) >= next_progress_episode
            _print_training_progress(session, summaries, episode_budget,
                                     target_global_step, active_workers,
                                     start_time, start_global_step)
            last_progress_episode = length(summaries)
            while next_progress_episode <= length(summaries)
                next_progress_episode += progress_frequency
            end
        end

        episode += length(batch_episodes)
    end

    if target_global_step != typemax(Int)
        if last_progress_step != session.learner.global_step
            _print_training_progress(session, summaries, episode_budget,
                                     target_global_step, active_workers,
                                     start_time, start_global_step)
        end
    elseif last_progress_episode != length(summaries)
        _print_training_progress(session, summaries, episode_budget,
                                 target_global_step, active_workers,
                                 start_time, start_global_step)
    end

    final_checkpoint_path = joinpath(session.output_dir, "checkpoint_final.jls")
    save_checkpoint(final_checkpoint_path, session.learner; manifest=session.manifest)
    @printf("training complete final_checkpoint=%s elapsed=%s\n",
            final_checkpoint_path, _format_duration(time() - start_time))
    flush(stdout)
    return (summaries=summaries, transitions=transitions,
            metrics=[episode_metrics(s; policy_name=algorithm_report_name(session.config.training.algorithm)) for s in summaries],
            aggregate=aggregate_metrics(summaries; policy_name=algorithm_report_name(session.config.training.algorithm)),
            global_step=session.learner.global_step,
            target_global_step=target_global_step == typemax(Int) ? nothing : target_global_step,
            output_dir=session.output_dir)
end

mutable struct A2CWorkerCampaign
    worker_id::Int
    episode_index::Int
    seed::Int
    scenario_rng::MersenneTwister
    policy_rng::MersenneTwister
    state::AerobrakingDecisionState
    observation::Vector{Float32}
    summary::EpisodeSummary
    ready::Bool
    summary_pending::Bool
end

struct A2CCollectedStep
    transition::Transition
    episode_end::Bool
end

function _a2c_worker_policy_seed(base_seed::Int, worker_id::Int)
    return base_seed + 1_000_000_007 + 10_000 * worker_id
end

function _initialize_a2c_worker_campaign!(worker::A2CWorkerCampaign,
                                          session::TrainingSession{<:A2CLearner},
                                          episode_index::Int, seed::Int)
    config = session.config.scenario
    worker.episode_index = episode_index
    worker.seed = seed
    worker.scenario_rng = MersenneTwister(seed)
    worker.state = reset_scenario(config, worker.scenario_rng)
    worker.summary = empty_episode_summary(episode_index=episode_index,
                                           worker_id=worker.worker_id,
                                           seed=seed)
    initial = run_protected_initializer(
        config,
        worker.state,
        worker.scenario_rng,
        worker.summary;
        settings=protected_initialization_config(session.config.training),
    )
    worker.state = initial.state
    worker.observation = initial.normalized_observation
    worker.summary = initial.summary
    pass_cap = min(session.config.training.max_passes_per_campaign,
                   config.termination_config.max_passes)
    worker.ready = !initial.done && worker.summary.pass_count < pass_cap
    worker.summary_pending = !worker.ready
    return worker
end

function _new_a2c_worker_campaign(session::TrainingSession{<:A2CLearner}, worker_id::Int,
                                  episode_index::Int, seed::Int)
    config = session.config.scenario
    scenario_rng = MersenneTwister(seed)
    state = reset_scenario(config, scenario_rng)
    observation = normalize_observation(observe_state(config, state),
                                        config.normalization_bounds)
    summary = empty_episode_summary(episode_index=episode_index, worker_id=worker_id, seed=seed)
    worker = A2CWorkerCampaign(
        worker_id,
        episode_index,
        seed,
        scenario_rng,
        MersenneTwister(_a2c_worker_policy_seed(session.config.training.seed, worker_id)),
        state,
        observation,
        summary,
        true,
        false,
    )
    return _initialize_a2c_worker_campaign!(worker, session, episode_index, seed)
end

function _prepare_a2c_worker!(session::TrainingSession{<:A2CLearner},
                              worker::A2CWorkerCampaign,
                              summaries::AbstractVector{EpisodeSummary},
                              episode_budget::Int,
                              next_episode_ref::Base.RefValue{Int};
                              allow_new_episode::Bool=true)
    while !worker.ready
        if worker.summary_pending
            push!(summaries, finalize_episode_summary(worker.summary, session.config.scenario))
            worker.summary_pending = false
        end
        allow_new_episode || return false
        next_episode_ref[] < episode_budget || return false
        next_episode_ref[] += 1
        seed = _a2c_worker_seed(session.config.training.seed,
                                worker.worker_id,
                                next_episode_ref[])
        _initialize_a2c_worker_campaign!(worker, session, next_episode_ref[], seed)
    end
    return true
end

function _a2c_worker_seed(base_seed::Int, worker_id::Int, episode_index::Int)
    return base_seed + 10_000 * worker_id + episode_index
end

function _a2c_next_values(learner::A2CLearner, next_observations::Array{Float32,3},
                          valid::AbstractMatrix{Bool})
    size(next_observations, 2) == size(valid, 1) ||
        throw(DimensionMismatch("next-observation worker dimension does not match valid mask"))
    size(next_observations, 3) == size(valid, 2) ||
        throw(DimensionMismatch("next-observation time dimension does not match valid mask"))
    observations = Matrix{Float32}(undef, learner.config.obs_dim, count(valid))
    coordinates = Tuple{Int,Int}[]
    column = 1
    for t in axes(valid, 2), worker in axes(valid, 1)
        valid[worker, t] || continue
        observations[:, column] .= @view next_observations[:, worker, t]
        push!(coordinates, (worker, t))
        column += 1
    end
    predictions = value_predictions(learner, observations)
    next_values = zeros(Float32, size(valid))
    for (value, (worker, t)) in zip(predictions, coordinates)
        next_values[worker, t] = value
    end
    return next_values
end

function _a2c_batch_from_worker_rollouts(
    learner::A2CLearner,
    worker_rollouts::AbstractVector{<:AbstractVector{A2CCollectedStep}},
    policy_version::Integer,
)
    n_workers = length(worker_rollouts)
    segment_length = maximum(length, worker_rollouts; init=0)
    segment_length > 0 || return nothing
    observations = zeros(Float32, learner.config.obs_dim, n_workers, segment_length)
    next_observations = zeros(Float32, learner.config.obs_dim, n_workers, segment_length)
    actions = zeros(Int, n_workers, segment_length)
    rewards = zeros(Float32, n_workers, segment_length)
    episode_end = falses(n_workers, segment_length)
    terminated = falses(n_workers, segment_length)
    valid = falses(n_workers, segment_length)
    for (worker, rollout) in pairs(worker_rollouts), (t, collected) in pairs(rollout)
        transition = collected.transition
        observations[:, worker, t] .= transition.observation
        next_observations[:, worker, t] .= transition.next_observation
        actions[worker, t] = transition.action_index
        rewards[worker, t] = transition.reward
        episode_end[worker, t] = collected.episode_end
        terminated[worker, t] = transition.terminated
        valid[worker, t] = true
    end
    next_values = _a2c_next_values(learner, next_observations, valid)
    returns = compute_discounted_returns(
        rewards,
        episode_end,
        terminated,
        valid,
        next_values,
        learner.config.discount,
    )
    return flatten_rollout(
        observations,
        actions,
        returns,
        valid;
        policy_version=policy_version,
    )
end

function _a2c_rollout_quotas(worker_ids::AbstractVector{Int}, remaining_steps::Int,
                             segment_length::Int)
    quotas = Dict(worker_id => 0 for worker_id in worker_ids)
    remaining = remaining_steps
    for _ in 1:segment_length
        for worker_id in worker_ids
            remaining > 0 || return quotas
            quotas[worker_id] += 1
            remaining -= 1
        end
    end
    return quotas
end

function _collect_a2c_segment!(session::TrainingSession{<:A2CLearner},
                               workers::Vector{A2CWorkerCampaign},
                               summaries::AbstractVector{EpisodeSummary},
                               transitions::AbstractVector{Transition},
                               target_global_step::Int,
                               episode_budget::Int,
                               next_episode_ref::Base.RefValue{Int})
    learner = session.learner
    config = session.config.scenario
    n_workers = length(workers)
    segment_length = learner.config.segment_length
    observations = zeros(Float32, learner.config.obs_dim, n_workers, segment_length)
    actions = zeros(Int, n_workers, segment_length)
    rewards = zeros(Float32, n_workers, segment_length)
    next_observations = zeros(Float32, learner.config.obs_dim, n_workers, segment_length)
    episode_end = falses(n_workers, segment_length)
    terminated = falses(n_workers, segment_length)
    valid = falses(n_workers, segment_length)
    actor_snapshot = cpu_network(learner.actor)
    policy_version = learner.policy_version

    for t in 1:segment_length
        learner.global_step < target_global_step || break
        length(summaries) < episode_budget || break

        for worker in workers
            _prepare_a2c_worker!(
                session,
                worker,
                summaries,
                episode_budget,
                next_episode_ref;
                allow_new_episode=learner.global_step < target_global_step,
            )
        end

        step_jobs = Vector{Union{Nothing,Task}}(nothing, n_workers)
        chosen_actions = Vector{Int}(undef, n_workers)
        previous_observations = Vector{Vector{Float32}}(undef, n_workers)
        remaining_steps = target_global_step - learner.global_step
        launched = 0

        for (worker_index, worker) in pairs(workers)
            worker.ready || continue
            launched < remaining_steps || break
            previous_observations[worker_index] = copy(worker.observation)
            chosen_actions[worker_index] = actor_action(actor_snapshot, worker.observation,
                                                        worker.policy_rng; test=false)
            observations[:, worker_index, t] .= worker.observation
            actions[worker_index, t] = chosen_actions[worker_index]
            step_jobs[worker_index] = Threads.@spawn step_scenario(
                $config,
                $(worker.state),
                $(chosen_actions[worker_index]),
                $(worker.scenario_rng),
            )
            launched += 1
        end
        launched == 0 && break

        for (worker_index, worker) in pairs(workers)
            step_jobs[worker_index] === nothing && continue
            result = fetch(step_jobs[worker_index])
            transition = transition_from_step(previous_observations[worker_index],
                                              chosen_actions[worker_index],
                                              result,
                                              length(transitions) + 1)
            push!(transitions, transition)
            valid[worker_index, t] = true
            rewards[worker_index, t] = transition.reward
            next_observations[:, worker_index, t] .= transition.next_observation
            learner.global_step += 1

            worker.summary = update_episode_summary(worker.summary, result)
            worker.state = result.state
            worker.observation = result.normalized_observation
            campaign_done = transition.terminated || transition.truncated ||
                            worker.summary.pass_count >= session.config.training.max_passes_per_campaign
            episode_end[worker_index, t] = campaign_done
            terminated[worker_index, t] = transition.terminated

            if campaign_done
                worker.ready = false
                worker.summary_pending = true
                _prepare_a2c_worker!(
                    session,
                    worker,
                    summaries,
                    episode_budget,
                    next_episode_ref;
                    allow_new_episode=learner.global_step < target_global_step,
                )
            end
        end
    end

    count(valid) == 0 && return nothing
    next_values = _a2c_next_values(learner, next_observations, valid)
    returns = compute_discounted_returns(rewards, episode_end, terminated, valid,
                                         next_values, learner.config.discount)
    return flatten_rollout(observations, actions, returns, valid;
                           policy_version=policy_version)
end

function _print_a2c_progress(session::TrainingSession{<:A2CLearner},
                             summaries::AbstractVector{<:EpisodeSummary},
                             episode_budget::Int,
                             target_global_step::Int,
                             active_workers::Int,
                             start_time::Float64,
                             start_global_step::Int)
    completed_episodes = length(summaries)
    elapsed = time() - start_time
    step_limited = target_global_step != typemax(Int)
    work_done = step_limited ? session.learner.global_step - start_global_step : completed_episodes
    work_remaining = step_limited ? max(0, target_global_step - session.learner.global_step) :
                     max(0, episode_budget - completed_episodes)
    work_rate = elapsed > 0 ? work_done / elapsed : 0.0
    eta = work_rate > 0 ? work_remaining / work_rate : Inf
    stats = _recent_training_stats(
        summaries,
        100;
        target_tolerance_m=session.config.scenario.reward_config.target_tolerance_m,
    )
    loss = isfinite(session.learner.last_loss) ? @sprintf("%.6g", session.learner.last_loss) : "n/a"
    policy_loss = isfinite(session.learner.last_policy_loss) ?
                  @sprintf("%.6g", session.learner.last_policy_loss) : "n/a"
    value_loss = isfinite(session.learner.last_value_loss) ?
                 @sprintf("%.6g", session.learner.last_value_loss) : "n/a"
    entropy = isfinite(session.learner.last_entropy) ?
              @sprintf("%.6g", session.learner.last_entropy) : "n/a"
    explained_variance = isfinite(session.learner.last_explained_variance) ?
                         @sprintf("%.4f", session.learner.last_explained_variance) : "n/a"
    @printf(
        "progress algo=a2c ep=%d/%s steps=%d%s train_steps=%d policy_version=%d loss=%s policy_loss=%s value_loss=%s entropy=%s explained_variance=%s recent_reward=%.3f recent_mean_thermal_violations=%.2f recent_mean_passes_to_end=%.1f recent_reached_goal_no_thermal=%.1f%% recent_mean_end_distance_no_thermal_km=%.2f workers=%d/%d elapsed=%s eta=%s\n",
        completed_episodes,
        _budget_label(episode_budget),
        session.learner.global_step,
        step_limited ? "/$(target_global_step)" : "",
        session.learner.train_steps,
        session.learner.policy_version,
        loss,
        policy_loss,
        value_loss,
        entropy,
        explained_variance,
        stats.mean_reward,
        stats.mean_thermal_violations,
        stats.mean_passes_to_end,
        stats.reached_goal_percent,
        stats.mean_end_distance_km,
        active_workers,
        Threads.nthreads(),
        _format_duration(elapsed),
        isfinite(eta) ? _format_duration(eta) : "n/a",
    )
    flush(stdout)
end

function _shutdown_a2c_live_workers!(active::Dict{Int,SpaceAGORAPhysicsActiveWorker},
                                     parked::Set{Int}, event_channel)
    for worker_id in collect(parked)
        worker = get(active, worker_id, nothing)
        worker === nothing && continue
        put!(worker.handle.action_channel, nothing)
        _finish_spaceagora_physics_streaming_worker!(worker)
        delete!(active, worker_id)
    end
    empty!(parked)
    isempty(active) || _stop_spaceagora_physics_streaming_workers!(active, event_channel)
    return nothing
end

function _train_parallel_a2c_spaceagora_physics_streaming!(
    session::TrainingSession{<:A2CLearner};
    global_steps::Int=session.config.training.global_steps,
    episodes::Union{Nothing,Int}=nothing,
    n_workers::Int=session.config.training.n_workers,
    process_ids::Vector{Int}=Int[],
)
    summaries = DiskBackedHistory(
        EpisodeSummary,
        joinpath(session.output_dir, "training_episode_summaries.jls"),
    )
    transitions = DiskBackedHistory(
        Transition,
        joinpath(session.output_dir, "training_transitions.jls"),
    )
    aggregate_accumulator = EpisodeAggregateAccumulator()
    requested_workers = max(1, n_workers)
    worker_backend = isempty(process_ids) ? :threads : :processes
    active_worker_limit = worker_backend == :threads ?
                          min(requested_workers, Threads.nthreads()) :
                          min(requested_workers, length(process_ids))
    if worker_backend == :threads && active_worker_limit < requested_workers
        @warn "n_workers exceeds available Julia threads; using Threads.nthreads()" requested_workers active_workers=active_worker_limit
    end
    target_global_step = global_steps > 0 ? global_steps : typemax(Int)
    episode_budget = episodes === nothing ?
                     (global_steps > 0 ? typemax(Int) : session.config.training.episodes) :
                     max(0, episodes)
    active_worker_limit = min(active_worker_limit, episode_budget)
    checkpoint_frequency = session.config.training.checkpoint_frequency
    next_checkpoint_step = checkpoint_frequency > 0 ? checkpoint_frequency : typemax(Int)
    progress_frequency = max(0, session.config.training.progress_frequency)
    next_progress_step = progress_frequency > 0 && target_global_step != typemax(Int) ?
                         min(target_global_step, session.learner.global_step + progress_frequency) :
                         typemax(Int)
    next_progress_episode = progress_frequency > 0 && target_global_step == typemax(Int) ?
                            min(progress_frequency, episode_budget) :
                            typemax(Int)
    start_time = time()
    start_global_step = session.learner.global_step
    last_progress_step = -1
    last_progress_episode = -1

    @printf(
        "starting A2C training global_steps=%s episode_cap=%s n_workers=%d active_workers=%d worker_backend=%s julia_threads=%d segment_length=%d train_start=%d checkpoint_frequency=%d progress_frequency=%d device=%s architecture=on_policy_pass_streaming\n",
        target_global_step == typemax(Int) ? "none" : string(target_global_step),
        _budget_label(episode_budget),
        requested_workers,
        active_worker_limit,
        string(worker_backend),
        Threads.nthreads(),
        session.learner.config.segment_length,
        session.learner.config.train_start,
        checkpoint_frequency,
        progress_frequency,
        training_device_name(session.learner.device),
    )
    @printf("output_dir=%s\n", session.output_dir)
    flush(stdout)

    event_channel = worker_backend == :processes ?
                    RemoteChannel(() -> Channel{Any}(max(32, 2 * active_worker_limit)), myid()) :
                    Channel{Any}(max(32, 2 * active_worker_limit))
    active = Dict{Int,SpaceAGORAPhysicsActiveWorker}()
    parked = Set{Int}()
    parked_observations = Dict{Int,Vector{Float32}}()
    simulation_templates = Dict{Int,SpaceAGORAPhysicsSimulationTemplate}()
    pending_repeats = Dict{Int,Any}()
    action_rngs = Dict(
        worker_id => MersenneTwister(
            _a2c_worker_policy_seed(session.config.training.seed, worker_id),
        )
        for worker_id in 1:active_worker_limit
    )
    next_episode = 1

    try
        while session.learner.global_step < target_global_step &&
              length(summaries) < episode_budget
            worker_ids = sort!(collect(keys(active)))
            available_new_episodes = episode_budget == typemax(Int) ?
                                     active_worker_limit :
                                     max(0, episode_budget - next_episode + 1)
            for worker_id in 1:active_worker_limit
                worker_id in worker_ids && continue
                available_new_episodes > 0 || break
                push!(worker_ids, worker_id)
                available_new_episodes -= 1
            end
            sort!(unique!(worker_ids))
            isempty(worker_ids) && break

            remaining_steps = target_global_step - session.learner.global_step
            quotas = _a2c_rollout_quotas(
                worker_ids,
                remaining_steps,
                session.learner.config.segment_length,
            )
            policy_version = session.learner.policy_version
            actor_snapshot = cpu_network(session.learner.actor)
            rollout_by_worker = Dict(worker_id => A2CCollectedStep[] for worker_id in worker_ids)
            collected = Dict(worker_id => 0 for worker_id in worker_ids)

            participating = Int[]
            for worker_id in worker_ids
                quota = quotas[worker_id]
                if quota == 0
                    if worker_id in parked
                        worker = active[worker_id]
                        put!(worker.handle.action_channel, nothing)
                        _finish_spaceagora_physics_streaming_worker!(worker)
                        delete!(active, worker_id)
                        delete!(parked, worker_id)
                        delete!(parked_observations, worker_id)
                    end
                    continue
                end
                if haskey(active, worker_id)
                    worker_id in parked ||
                        throw(ErrorException("A2C live worker $worker_id was not parked at an update boundary"))
                    observation = parked_observations[worker_id]
                    action_index = actor_action(
                        actor_snapshot,
                        observation,
                        action_rngs[worker_id];
                        test=false,
                    )
                    put!(
                        active[worker_id].handle.action_channel,
                        (action_index=action_index, protected=false,
                         policy_version=policy_version),
                    )
                    delete!(parked, worker_id)
                    delete!(parked_observations, worker_id)
                    push!(participating, worker_id)
                else
                    next_episode <= episode_budget || continue
                    repeat_plan = pop!(pending_repeats, worker_id, nothing)
                    active[worker_id] = _launch_spaceagora_physics_streaming_worker!(
                        session,
                        event_channel,
                        worker_id,
                        next_episode,
                        actor_snapshot,
                        action_rngs[worker_id],
                        policy_version;
                        simulation_template=get(simulation_templates, worker_id, nothing),
                        process_id=worker_backend == :processes ? process_ids[worker_id] : nothing,
                        scenario_seed=repeat_plan === nothing ? nothing : repeat_plan.seed,
                        successful_case_repeat_index=
                            repeat_plan === nothing ? 0 : repeat_plan.repeat_index,
                    )
                    next_episode += 1
                    push!(participating, worker_id)
                end
            end
            worker_ids = participating
            isempty(worker_ids) && break

            while any(worker_id -> collected[worker_id] < quotas[worker_id], worker_ids)
                event = take!(event_channel)
                worker = get(active, event.worker_id, nothing)
                worker === nothing && continue
                if event.policy_version != policy_version
                    event.done || put!(worker.handle.action_channel, nothing)
                    _finish_spaceagora_physics_streaming_worker!(worker)
                    delete!(active, event.worker_id)
                    throw(ErrorException(
                        "A2C worker $(event.worker_id) returned policy version $(event.policy_version); expected $policy_version",
                    ))
                end
                if event.error !== nothing
                    event.done || put!(worker.handle.action_channel, nothing)
                    _finish_spaceagora_physics_streaming_worker!(worker)
                    delete!(active, event.worker_id)
                    throw(ErrorException(event.error))
                end

                if event.protected
                    worker.protected_events_seen += 1
                elseif event.transition !== nothing
                    transition = event.transition
                    push!(transitions, transition)
                    push!(
                        rollout_by_worker[event.worker_id],
                        A2CCollectedStep(transition, event.done),
                    )
                    collected[event.worker_id] += 1
                    session.learner.global_step += 1
                end

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
                    worker.simulation_template === nothing ||
                        (simulation_templates[event.worker_id] = worker.simulation_template)
                    _finish_spaceagora_physics_streaming_worker!(worker)
                    delete!(active, event.worker_id)

                    can_continue = collected[event.worker_id] < quotas[event.worker_id] &&
                                   next_episode <= episode_budget &&
                                   session.learner.global_step < target_global_step
                    if can_continue
                        active[event.worker_id] = _launch_spaceagora_physics_streaming_worker!(
                            session,
                            event_channel,
                            event.worker_id,
                            next_episode,
                            actor_snapshot,
                            action_rngs[event.worker_id],
                            policy_version;
                            simulation_template=get(simulation_templates, event.worker_id, nothing),
                            process_id=worker_backend == :processes ?
                                       process_ids[event.worker_id] : nothing,
                            scenario_seed=repeat_plan === nothing ? nothing : repeat_plan.seed,
                            successful_case_repeat_index=
                                repeat_plan === nothing ? 0 : repeat_plan.repeat_index,
                        )
                        next_episode += 1
                    else
                        repeat_plan === nothing ||
                            (pending_repeats[event.worker_id] = repeat_plan)
                        if collected[event.worker_id] < quotas[event.worker_id]
                            quotas[event.worker_id] = collected[event.worker_id]
                        end
                    end
                elseif collected[event.worker_id] >= quotas[event.worker_id]
                    transition = event.transition
                    transition === nothing &&
                        throw(ErrorException("A2C worker parked without a transition observation"))
                    push!(parked, event.worker_id)
                    parked_observations[event.worker_id] = transition.next_observation
                else
                    corridor_action = event.protected && worker.protected_events_seen == 1 ?
                                      _protected_corridor_action_index(
                                          session.config.training,
                                          event.result,
                                      ) :
                                      nothing
                    if corridor_action === nothing
                        transition = event.transition
                        transition === nothing &&
                            throw(ErrorException("A2C worker event is missing its next observation"))
                        next_action = actor_action(
                            actor_snapshot,
                            transition.next_observation,
                            action_rngs[event.worker_id];
                            test=false,
                        )
                        put!(
                            worker.handle.action_channel,
                            (action_index=next_action, protected=false,
                             policy_version=policy_version),
                        )
                    else
                        put!(
                            worker.handle.action_channel,
                            (action_index=corridor_action, protected=true,
                             policy_version=policy_version),
                        )
                    end
                end
            end

            ordered_rollouts = [rollout_by_worker[worker_id] for worker_id in worker_ids]
            batch = _a2c_batch_from_worker_rollouts(
                session.learner,
                ordered_rollouts,
                policy_version,
            )
            batch === nothing || maybe_train!(session.learner, batch)

            while session.learner.global_step >= next_checkpoint_step
                checkpoint_path = joinpath(
                    session.output_dir,
                    "checkpoint_$(next_checkpoint_step).jls",
                )
                save_checkpoint(checkpoint_path, session.learner; manifest=session.manifest)
                @printf("checkpoint step=%d path=%s\n", next_checkpoint_step, checkpoint_path)
                flush(stdout)
                next_checkpoint_step += checkpoint_frequency
            end

            if progress_frequency > 0 && target_global_step != typemax(Int) &&
               session.learner.global_step >= next_progress_step
                _print_a2c_progress(
                    session,
                    summaries,
                    episode_budget,
                    target_global_step,
                    active_worker_limit,
                    start_time,
                    start_global_step,
                )
                last_progress_step = session.learner.global_step
                while next_progress_step <= session.learner.global_step &&
                      next_progress_step < target_global_step
                    next_progress_step += progress_frequency
                end
                next_progress_step = min(next_progress_step, target_global_step)
            elseif progress_frequency > 0 && target_global_step == typemax(Int) &&
                   length(summaries) >= next_progress_episode
                _print_a2c_progress(
                    session,
                    summaries,
                    episode_budget,
                    target_global_step,
                    active_worker_limit,
                    start_time,
                    start_global_step,
                )
                last_progress_episode = length(summaries)
                while next_progress_episode <= length(summaries)
                    next_progress_episode += progress_frequency
                end
            end
        end
    catch
        _shutdown_a2c_live_workers!(active, parked, event_channel)
        close_history!(summaries)
        close_history!(transitions)
        rethrow()
    end

    _shutdown_a2c_live_workers!(active, parked, event_channel)
    if target_global_step != typemax(Int)
        if last_progress_step != session.learner.global_step
            _print_a2c_progress(session, summaries, episode_budget, target_global_step,
                                active_worker_limit, start_time, start_global_step)
        end
    elseif last_progress_episode != length(summaries)
        _print_a2c_progress(session, summaries, episode_budget, target_global_step,
                            active_worker_limit, start_time, start_global_step)
    end

    final_checkpoint_path = joinpath(session.output_dir, "checkpoint_final.jls")
    save_checkpoint(final_checkpoint_path, session.learner; manifest=session.manifest)
    close_history!(summaries)
    close_history!(transitions)
    @printf("training complete final_checkpoint=%s elapsed=%s\n",
            final_checkpoint_path, _format_duration(time() - start_time))
    flush(stdout)
    return (
        summaries=summaries,
        transitions=transitions,
        metrics=MappedHistory(summaries, summary -> episode_metrics(summary; policy_name="a2c")),
        aggregate=aggregate_metrics(aggregate_accumulator; policy_name="a2c"),
        global_step=session.learner.global_step,
        target_global_step=target_global_step == typemax(Int) ? nothing : target_global_step,
        output_dir=session.output_dir,
    )
end

function train_parallel!(session::TrainingSession{<:A2CLearner};
                         global_steps::Int=session.config.training.global_steps,
                         episodes::Union{Nothing,Int}=nothing,
                         n_workers::Int=session.config.training.n_workers)
    if _is_spaceagora_live_backend(session.config.scenario.backend_mode)
        worker_backend = session.config.training.worker_backend
        active_workers = worker_backend == :processes ?
                         max(1, n_workers) :
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
                    return _train_parallel_a2c_spaceagora_physics_streaming!(
                        session;
                        global_steps=global_steps,
                        episodes=episodes,
                        n_workers=n_workers,
                        process_ids=process_ids,
                    )
                finally
                    rmprocs(process_ids)
                end
            end
            _prewarm_spaceagora_rl_shared_ephemeris_cache!(
                session.config.scenario,
                session.config.training.max_passes_per_campaign,
            )
            return _train_parallel_a2c_spaceagora_physics_streaming!(
                session;
                global_steps=global_steps,
                episodes=episodes,
                n_workers=n_workers,
            )
        end
    end

    summaries = DiskBackedHistory(
        EpisodeSummary,
        joinpath(session.output_dir, "training_episode_summaries.jls"),
    )
    transitions = DiskBackedHistory(
        Transition,
        joinpath(session.output_dir, "training_transitions.jls"),
    )
    requested_workers = max(1, n_workers)
    target_global_step = global_steps > 0 ? global_steps : typemax(Int)
    episode_budget = episodes === nothing ? (global_steps > 0 ? typemax(Int) : session.config.training.episodes) : max(0, episodes)
    thread_workers = min(requested_workers, Threads.nthreads())
    thread_workers < requested_workers && @warn "n_workers exceeds available Julia threads; using Threads.nthreads()" requested_workers active_workers=thread_workers
    active_workers = min(thread_workers, episode_budget)
    checkpoint_frequency = session.config.training.checkpoint_frequency
    next_checkpoint_step = checkpoint_frequency > 0 ? checkpoint_frequency : typemax(Int)
    progress_frequency = max(0, session.config.training.progress_frequency)
    next_progress_step = progress_frequency > 0 && target_global_step != typemax(Int) ?
                         min(target_global_step, session.learner.global_step + progress_frequency) :
                         typemax(Int)
    next_progress_episode = progress_frequency > 0 && target_global_step == typemax(Int) ?
                            min(progress_frequency, episode_budget) :
                            typemax(Int)
    start_time = time()
    start_global_step = session.learner.global_step
    last_progress_step = -1
    last_progress_episode = -1
    printed_initial_progress = false

    @printf(
        "starting A2C training global_steps=%s episode_cap=%s n_workers=%d active_workers=%d julia_threads=%d segment_length=%d train_start=%d checkpoint_frequency=%d progress_frequency=%d device=%s\n",
        target_global_step == typemax(Int) ? "none" : string(target_global_step),
        _budget_label(episode_budget),
        requested_workers,
        active_workers,
        Threads.nthreads(),
        session.learner.config.segment_length,
        session.learner.config.train_start,
        checkpoint_frequency,
        progress_frequency,
        training_device_name(session.learner.device),
    )
    @printf("output_dir=%s\n", session.output_dir)
    flush(stdout)

    next_episode_ref = Ref(active_workers)
    workers = [
        _new_a2c_worker_campaign(
            session,
            worker_id,
            worker_id,
            _a2c_worker_seed(session.config.training.seed, worker_id, worker_id),
        )
        for worker_id in 1:active_workers
    ]

    while session.learner.global_step < target_global_step && length(summaries) < episode_budget
        batch = _collect_a2c_segment!(session, workers, summaries, transitions,
                                      target_global_step, episode_budget, next_episode_ref)
        batch === nothing && break
        maybe_train!(session.learner, batch)

        while session.learner.global_step >= next_checkpoint_step
            checkpoint_path = joinpath(session.output_dir, "checkpoint_$(next_checkpoint_step).jls")
            save_checkpoint(checkpoint_path, session.learner; manifest=session.manifest)
            @printf("checkpoint step=%d path=%s\n", next_checkpoint_step, checkpoint_path)
            flush(stdout)
            next_checkpoint_step += checkpoint_frequency
        end

        if progress_frequency > 0 && !printed_initial_progress &&
           session.learner.global_step > start_global_step
            _print_a2c_progress(session, summaries, episode_budget, target_global_step,
                                active_workers, start_time, start_global_step)
            printed_initial_progress = true
            if target_global_step != typemax(Int)
                last_progress_step = session.learner.global_step
                while next_progress_step <= session.learner.global_step &&
                      next_progress_step < target_global_step
                    next_progress_step += progress_frequency
                end
                next_progress_step = min(next_progress_step, target_global_step)
            else
                last_progress_episode = length(summaries)
                while next_progress_episode <= length(summaries)
                    next_progress_episode += progress_frequency
                end
            end
        elseif progress_frequency > 0 && target_global_step != typemax(Int) &&
               session.learner.global_step >= next_progress_step
            _print_a2c_progress(session, summaries, episode_budget, target_global_step,
                                active_workers, start_time, start_global_step)
            last_progress_step = session.learner.global_step
            while next_progress_step <= session.learner.global_step &&
                  next_progress_step < target_global_step
                next_progress_step += progress_frequency
            end
            next_progress_step = min(next_progress_step, target_global_step)
        elseif progress_frequency > 0 && target_global_step == typemax(Int) &&
               length(summaries) >= next_progress_episode
            _print_a2c_progress(session, summaries, episode_budget, target_global_step,
                                active_workers, start_time, start_global_step)
            last_progress_episode = length(summaries)
            while next_progress_episode <= length(summaries)
                next_progress_episode += progress_frequency
            end
        end
    end

    if target_global_step != typemax(Int)
        if last_progress_step != session.learner.global_step
            _print_a2c_progress(session, summaries, episode_budget, target_global_step,
                                active_workers, start_time, start_global_step)
        end
    elseif last_progress_episode != length(summaries)
        _print_a2c_progress(session, summaries, episode_budget, target_global_step,
                            active_workers, start_time, start_global_step)
    end

    final_checkpoint_path = joinpath(session.output_dir, "checkpoint_final.jls")
    save_checkpoint(final_checkpoint_path, session.learner; manifest=session.manifest)
    close_history!(summaries)
    close_history!(transitions)
    @printf("training complete final_checkpoint=%s elapsed=%s\n",
            final_checkpoint_path, _format_duration(time() - start_time))
    flush(stdout)
    return (summaries=summaries, transitions=transitions,
            metrics=MappedHistory(summaries, summary -> episode_metrics(summary; policy_name="a2c")),
            aggregate=aggregate_metrics(summaries; policy_name="a2c"),
            global_step=session.learner.global_step,
            target_global_step=target_global_step == typemax(Int) ? nothing : target_global_step,
            output_dir=session.output_dir)
end
