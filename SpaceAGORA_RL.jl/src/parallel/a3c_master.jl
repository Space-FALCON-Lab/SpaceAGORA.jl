mutable struct A3CSurrogateWorker
    worker_id::Int
    episode_index::Int
    seed::Int
    successful_case_repeat_index::Int
    scenario_rng::MersenneTwister
    policy_rng::MersenneTwister
    state::AerobrakingDecisionState
    observation::Vector{Float32}
    summary::EpisodeSummary
    local_model::A3CLocalModel
    rollout::A3CRollout
    ready::Bool
end

mutable struct A3CLiveWorkerState
    local_model::A3CLocalModel
    rollout::A3CRollout
    policy_rng::MersenneTwister
end

_a3c_worker_policy_seed(base_seed::Int, worker_id::Int) =
    base_seed + 1_500_000_017 + 10_000 * worker_id

_a3c_worker_seed(base_seed::Int, worker_id::Int, episode_index::Int) =
    base_seed + 10_000 * worker_id + episode_index

_a3c_format_metric(value::Real) = isfinite(value) ? @sprintf("%.6g", value) : "n/a"

function _a3c_prewarm_spaceagora_cache!(config::AerobrakingScenarioConfig,
                                        max_passes_per_campaign::Integer)
    _prewarm_spaceagora_rl_shared_ephemeris_cache!(config, max_passes_per_campaign)
    return nothing
end

function _initialize_a3c_surrogate_worker!(worker::A3CSurrogateWorker,
                                            session::TrainingSession{<:A3CLearner},
                                            episode_index::Int, seed::Int;
                                            successful_case_repeat_index::Int=0)
    config = session.config.scenario
    worker.episode_index = episode_index
    worker.seed = seed
    worker.successful_case_repeat_index = successful_case_repeat_index
    worker.scenario_rng = MersenneTwister(seed)
    worker.state = reset_scenario(config, worker.scenario_rng)
    worker.summary = empty_episode_summary(
        episode_index=episode_index,
        worker_id=worker.worker_id,
        seed=seed,
    )
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
    pass_cap = min(
        session.config.training.max_passes_per_campaign,
        config.termination_config.max_passes,
    )
    worker.ready = !initial.done && worker.summary.pass_count < pass_cap
    return worker
end

function _new_a3c_surrogate_worker(session::TrainingSession{<:A3CLearner},
                                   worker_id::Int, episode_index::Int, seed::Int)
    config = session.config.scenario
    scenario_rng = MersenneTwister(seed)
    state = reset_scenario(config, scenario_rng)
    worker = A3CSurrogateWorker(
        worker_id,
        episode_index,
        seed,
        0,
        scenario_rng,
        MersenneTwister(_a3c_worker_policy_seed(session.config.training.seed, worker_id)),
        state,
        normalize_observation(observe_state(config, state), config.normalization_bounds),
        empty_episode_summary(episode_index=episode_index, worker_id=worker_id, seed=seed),
        A3CLocalModel(session.learner),
        A3CRollout(),
        true,
    )
    return _initialize_a3c_surrogate_worker!(worker, session, episode_index, seed)
end

function _a3c_update_rollout!(learner::A3CLearner, local_model::A3CLocalModel,
                              rollout::A3CRollout)
    isempty(rollout) && return nothing
    batch = a3c_rollout_batch(rollout, local_model, learner.config)
    result = a3c_update!(learner, local_model, batch)
    empty!(rollout)
    return result
end

function _a3c_prepare_surrogate_worker!(session::TrainingSession{<:A3CLearner},
                                        worker::A3CSurrogateWorker,
                                        summaries,
                                        aggregate_accumulator::EpisodeAggregateAccumulator,
                                        episode_budget::Int,
                                        next_episode::Base.RefValue{Int})
    while !worker.ready
        final_summary = finalize_episode_summary(worker.summary, session.config.scenario)
        push!(summaries, final_summary)
        accumulate_episode!(aggregate_accumulator, final_summary)
        length(summaries) >= episode_budget && return false
        next_episode[] <= episode_budget || return false
        repeat_plan = _next_successful_case_repeat(
            worker.seed,
            worker.successful_case_repeat_index,
            final_summary.success,
            session.config.training.successful_case_repetitions,
        )
        seed = repeat_plan === nothing ?
               _a3c_worker_seed(session.config.training.seed, worker.worker_id, next_episode[]) :
               repeat_plan.seed
        repeat_index = repeat_plan === nothing ? 0 : repeat_plan.repeat_index
        _initialize_a3c_surrogate_worker!(
            worker,
            session,
            next_episode[],
            seed;
            successful_case_repeat_index=repeat_index,
        )
        next_episode[] += 1
    end
    return true
end

function _print_a3c_progress(session::TrainingSession{<:A3CLearner}, summaries,
                             episode_budget::Int, target_global_step::Int,
                             active_workers::Int, start_time::Float64,
                             start_global_step::Int)
    stats = _recent_training_stats(summaries, 100)
    elapsed = time() - start_time
    completed = session.learner.global_step - start_global_step
    eta = if target_global_step == typemax(Int) || completed <= 0
        NaN
    else
        elapsed * max(0, target_global_step - session.learner.global_step) / completed
    end
    @printf(
        "progress algo=a3c ep=%d/%s steps=%d%s train_steps=%d policy_version=%d loss=%s policy_loss=%s value_loss=%s entropy=%s explained_variance=%s mean_policy_lag=%s max_policy_lag=%d dropped_stale_updates=%d recent_reward=%.3f recent_mean_thermal_violations=%.2f workers=%d elapsed=%s eta=%s\n",
        length(summaries),
        _budget_label(episode_budget),
        session.learner.global_step,
        target_global_step == typemax(Int) ? "" : "/$(target_global_step)",
        session.learner.train_steps,
        session.learner.policy_version,
        _a3c_format_metric(session.learner.last_loss),
        _a3c_format_metric(session.learner.last_policy_loss),
        _a3c_format_metric(session.learner.last_value_loss),
        _a3c_format_metric(session.learner.last_entropy),
        _a3c_format_metric(session.learner.last_explained_variance),
        _a3c_format_metric(mean_policy_lag(session.learner)),
        session.learner.max_observed_policy_lag,
        session.learner.dropped_stale_updates,
        stats.mean_reward,
        stats.mean_thermal_violations,
        active_workers,
        _format_duration(elapsed),
        isfinite(eta) ? _format_duration(eta) : "n/a",
    )
    flush(stdout)
    return nothing
end

function _a3c_checkpoint_due!(session::TrainingSession{<:A3CLearner},
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

function _a3c_training_result(session::TrainingSession{<:A3CLearner}, summaries,
                              transitions, aggregate_accumulator,
                              target_global_step::Int, start_time::Float64)
    final_checkpoint_path = joinpath(session.output_dir, "checkpoint_final.jls")
    save_checkpoint(final_checkpoint_path, session.learner; manifest=session.manifest)
    close_history!(summaries)
    close_history!(transitions)
    @printf(
        "training complete final_checkpoint=%s elapsed=%s\n",
        final_checkpoint_path,
        _format_duration(time() - start_time),
    )
    flush(stdout)
    return (
        summaries=summaries,
        transitions=transitions,
        metrics=MappedHistory(summaries, summary -> episode_metrics(summary; policy_name="a3c")),
        aggregate=aggregate_metrics(aggregate_accumulator; policy_name="a3c"),
        global_step=session.learner.global_step,
        target_global_step=target_global_step == typemax(Int) ? nothing : target_global_step,
        output_dir=session.output_dir,
    )
end

function _train_parallel_a3c_surrogate!(session::TrainingSession{<:A3CLearner};
                                        global_steps::Int=session.config.training.global_steps,
                                        episodes::Union{Nothing,Int}=nothing,
                                        n_workers::Int=session.config.training.n_workers)
    summaries = DiskBackedHistory(
        EpisodeSummary,
        joinpath(session.output_dir, "training_episode_summaries.jls"),
    )
    transitions = DiskBackedHistory(
        Transition,
        joinpath(session.output_dir, "training_transitions.jls"),
    )
    aggregate_accumulator = EpisodeAggregateAccumulator()
    target_global_step = global_steps > 0 ? global_steps : typemax(Int)
    episode_budget = episodes === nothing ?
                     (global_steps > 0 ? typemax(Int) : session.config.training.episodes) :
                     max(0, episodes)
    active_workers = min(max(1, n_workers), episode_budget)
    checkpoint_frequency = session.config.training.checkpoint_frequency
    next_checkpoint_step = Ref(checkpoint_frequency > 0 ? checkpoint_frequency : typemax(Int))
    progress_frequency = max(0, session.config.training.progress_frequency)
    next_progress_step = Ref(progress_frequency > 0 ? progress_frequency : typemax(Int))
    start_time = time()
    start_global_step = session.learner.global_step

    @printf(
        "starting A3C training global_steps=%s episode_cap=%s n_workers=%d active_workers=%d t_max=%d device=%s architecture=asynchronous_parameter_server\n",
        target_global_step == typemax(Int) ? "none" : string(target_global_step),
        _budget_label(episode_budget),
        n_workers,
        active_workers,
        session.learner.config.t_max,
        training_device_name(session.learner.device),
    )
    @printf("output_dir=%s\n", session.output_dir)
    flush(stdout)

    workers = Dict{Int,A3CSurrogateWorker}()
    next_episode = Ref(1)
    for worker_id in 1:active_workers
        workers[worker_id] = _new_a3c_surrogate_worker(
            session,
            worker_id,
            next_episode[],
            _a3c_worker_seed(session.config.training.seed, worker_id, next_episode[]),
        )
        next_episode[] += 1
    end
    event_channel = Channel{Any}(max(32, 2 * active_workers))
    inflight = Dict{Int,Task}()

    try
        while session.learner.global_step < target_global_step &&
              length(summaries) < episode_budget &&
              !isempty(workers)
            for worker_id in sort!(collect(keys(workers)))
                haskey(inflight, worker_id) && continue
                session.learner.global_step + length(inflight) < target_global_step || break
                worker = workers[worker_id]
                if !_a3c_prepare_surrogate_worker!(
                    session,
                    worker,
                    summaries,
                    aggregate_accumulator,
                    episode_budget,
                    next_episode,
                )
                    delete!(workers, worker_id)
                    continue
                end
                previous_observation = copy(worker.observation)
                selected_action = select_action(
                    worker.local_model,
                    previous_observation;
                    rng=worker.policy_rng,
                )
                action_index, _ = _resolve_policy_action(selected_action)
                task = Threads.@spawn begin
                    try
                        result = step_scenario(
                            $(session.config.scenario),
                            $(worker.state),
                            $selected_action,
                            $(worker.scenario_rng),
                        )
                        put!($event_channel, (
                            worker_id=$worker_id,
                            observation=$previous_observation,
                            action_index=$action_index,
                            result=result,
                            error=nothing,
                        ))
                    catch err
                        put!($event_channel, (
                            worker_id=$worker_id,
                            observation=$previous_observation,
                            action_index=$action_index,
                            result=nothing,
                            error=sprint(showerror, err, catch_backtrace()),
                        ))
                    end
                end
                inflight[worker_id] = task
            end
            isempty(inflight) && break

            event = take!(event_channel)
            fetch(pop!(inflight, event.worker_id))
            event.error === nothing || throw(ErrorException(event.error))
            worker = workers[event.worker_id]
            result = event.result::AerobrakingStepResult
            transition = transition_from_step(
                event.observation,
                event.action_index,
                result,
                length(transitions) + 1,
            )
            push!(transitions, transition)
            push_a3c_transition!(worker.rollout, transition, result.action)
            session.learner.global_step += 1
            worker.summary = update_episode_summary(worker.summary, result)
            worker.state = result.state
            worker.observation = result.normalized_observation
            pass_cap = min(
                session.config.training.max_passes_per_campaign,
                session.config.scenario.termination_config.max_passes,
            )
            campaign_done = transition.terminated || transition.truncated ||
                            worker.summary.pass_count >= pass_cap
            if length(worker.rollout) >= session.learner.config.t_max || campaign_done
                _a3c_update_rollout!(session.learner, worker.local_model, worker.rollout)
            end
            worker.ready = !campaign_done

            checkpoint_frequency > 0 && _a3c_checkpoint_due!(
                session,
                next_checkpoint_step,
                checkpoint_frequency,
            )
            if progress_frequency > 0 && session.learner.global_step >= next_progress_step[]
                _print_a3c_progress(
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
    finally
        foreach(wait, values(inflight))
    end

    for worker in values(workers)
        _a3c_update_rollout!(session.learner, worker.local_model, worker.rollout)
        if !worker.ready && length(summaries) < episode_budget
            final_summary = finalize_episode_summary(worker.summary, session.config.scenario)
            push!(summaries, final_summary)
            accumulate_episode!(aggregate_accumulator, final_summary)
        end
    end
    _print_a3c_progress(
        session,
        summaries,
        episode_budget,
        target_global_step,
        active_workers,
        start_time,
        start_global_step,
    )
    return _a3c_training_result(
        session,
        summaries,
        transitions,
        aggregate_accumulator,
        target_global_step,
        start_time,
    )
end

function _launch_a3c_live_worker!(session::TrainingSession{<:A3CLearner}, event_channel,
                                  worker_id::Int, episode_index::Int,
                                  worker_state::A3CLiveWorkerState;
                                  simulation_template::Union{Nothing,SpaceAGORAPhysicsSimulationTemplate}=nothing,
                                  process_id::Union{Nothing,Int}=nothing,
                                  scenario_seed::Union{Nothing,Int}=nothing,
                                  successful_case_repeat_index::Int=0)
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
    selected_action = protected_first_pass ?
                      zero_action_index() :
                      select_action(
                          worker_state.local_model,
                          norm_obs;
                          rng=worker_state.policy_rng,
                      )
    action_index, action = _resolve_policy_action(selected_action)
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
            policy_version=worker_state.local_model.policy_version,
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
            policy_version=worker_state.local_model.policy_version,
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

function _train_parallel_a3c_live!(session::TrainingSession{<:A3CLearner};
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
    active_worker_limit = worker_backend == :threads ?
                          min(requested_workers, Threads.nthreads()) :
                          min(requested_workers, length(process_ids))
    target_global_step = global_steps > 0 ? global_steps : typemax(Int)
    episode_budget = episodes === nothing ?
                     (global_steps > 0 ? typemax(Int) : session.config.training.episodes) :
                     max(0, episodes)
    active_worker_limit = min(active_worker_limit, episode_budget)
    checkpoint_frequency = session.config.training.checkpoint_frequency
    next_checkpoint_step = Ref(checkpoint_frequency > 0 ? checkpoint_frequency : typemax(Int))
    progress_frequency = max(0, session.config.training.progress_frequency)
    next_progress_step = Ref(progress_frequency > 0 ? progress_frequency : typemax(Int))
    start_time = time()
    start_global_step = session.learner.global_step

    @printf(
        "starting A3C training global_steps=%s episode_cap=%s n_workers=%d active_workers=%d worker_backend=%s t_max=%d device=%s architecture=asynchronous_parameter_server\n",
        target_global_step == typemax(Int) ? "none" : string(target_global_step),
        _budget_label(episode_budget),
        requested_workers,
        active_worker_limit,
        string(worker_backend),
        session.learner.config.t_max,
        training_device_name(session.learner.device),
    )
    @printf("output_dir=%s\n", session.output_dir)
    flush(stdout)

    event_channel = worker_backend == :processes ?
                    RemoteChannel(() -> Channel{Any}(max(32, 2 * active_worker_limit)), myid()) :
                    Channel{Any}(max(32, 2 * active_worker_limit))
    active = Dict{Int,SpaceAGORAPhysicsActiveWorker}()
    worker_states = Dict(
        worker_id => A3CLiveWorkerState(
            A3CLocalModel(session.learner),
            A3CRollout(),
            MersenneTwister(_a3c_worker_policy_seed(session.config.training.seed, worker_id)),
        )
        for worker_id in 1:active_worker_limit
    )
    next_episode = 1
    for worker_id in 1:active_worker_limit
        active[worker_id] = _launch_a3c_live_worker!(
            session,
            event_channel,
            worker_id,
            next_episode,
            worker_states[worker_id];
            process_id=worker_backend == :processes ? process_ids[worker_id] : nothing,
        )
        next_episode += 1
    end

    try
        while session.learner.global_step < target_global_step &&
              length(summaries) < episode_budget &&
              !isempty(active)
            event = take!(event_channel)
            worker = get(active, event.worker_id, nothing)
            worker === nothing && continue
            state = worker_states[event.worker_id]
            if event.policy_version != state.local_model.policy_version
                event.done || put!(worker.handle.action_channel, nothing)
                _finish_spaceagora_physics_streaming_worker!(worker)
                delete!(active, event.worker_id)
                throw(ErrorException(
                    "A3C worker $(event.worker_id) returned policy version $(event.policy_version); expected $(state.local_model.policy_version)",
                ))
            end
            if event.error !== nothing
                event.done || put!(worker.handle.action_channel, nothing)
                _finish_spaceagora_physics_streaming_worker!(worker)
                delete!(active, event.worker_id)
                throw(ErrorException(event.error))
            end

            event.protected && (worker.protected_events_seen += 1)
            if event.transition !== nothing && !event.protected &&
               session.learner.global_step < target_global_step
                transition = event.transition
                push!(transitions, transition)
                push_a3c_transition!(state.rollout, transition, event.result.action)
                session.learner.global_step += 1
                if length(state.rollout) >= session.learner.config.t_max || event.done
                    _a3c_update_rollout!(session.learner, state.local_model, state.rollout)
                end
                checkpoint_frequency > 0 && _a3c_checkpoint_due!(
                    session,
                    next_checkpoint_step,
                    checkpoint_frequency,
                )
            end

            reached_limits = session.learner.global_step >= target_global_step ||
                             length(summaries) >= episode_budget
            if event.done
                event.protected || _a3c_update_rollout!(session.learner, state.local_model, state.rollout)
                final_summary = finalize_episode_summary(event.summary, session.config.scenario)
                repeat_plan = nothing
                if length(summaries) < episode_budget
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
                if !reached_limits && length(summaries) < episode_budget &&
                   next_episode <= episode_budget
                    active[event.worker_id] = _launch_a3c_live_worker!(
                        session,
                        event_channel,
                        event.worker_id,
                        next_episode,
                        state;
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
                                  ) :
                                  nothing
                if corridor_action === nothing
                    transition = event.transition
                    transition === nothing && throw(ErrorException(
                        "A3C worker event is missing its next observation",
                    ))
                    next_action = select_action(
                        state.local_model,
                        transition.next_observation;
                        rng=state.policy_rng,
                    )
                    put!(
                        worker.handle.action_channel,
                        _actor_critic_action_command(
                            next_action;
                            protected=false,
                            policy_version=state.local_model.policy_version,
                        ),
                    )
                else
                    put!(
                        worker.handle.action_channel,
                        (action_index=corridor_action, protected=true,
                         policy_version=state.local_model.policy_version),
                    )
                end
            end

            if progress_frequency > 0 && session.learner.global_step >= next_progress_step[]
                _print_a3c_progress(
                    session,
                    summaries,
                    episode_budget,
                    target_global_step,
                    active_worker_limit,
                    start_time,
                    start_global_step,
                )
                while next_progress_step[] <= session.learner.global_step
                    next_progress_step[] += progress_frequency
                end
            end
        end
    finally
        isempty(active) || _stop_spaceagora_physics_streaming_workers!(active, event_channel)
    end

    for state in values(worker_states)
        _a3c_update_rollout!(session.learner, state.local_model, state.rollout)
    end
    _print_a3c_progress(
        session,
        summaries,
        episode_budget,
        target_global_step,
        active_worker_limit,
        start_time,
        start_global_step,
    )
    return _a3c_training_result(
        session,
        summaries,
        transitions,
        aggregate_accumulator,
        target_global_step,
        start_time,
    )
end

function train_parallel!(session::TrainingSession{<:A3CLearner};
                         global_steps::Int=session.config.training.global_steps,
                         episodes::Union{Nothing,Int}=nothing,
                         n_workers::Int=session.config.training.n_workers)
    if !_is_spaceagora_live_backend(session.config.scenario.backend_mode)
        return _train_parallel_a3c_surrogate!(
            session;
            global_steps=global_steps,
            episodes=episodes,
            n_workers=n_workers,
        )
    end

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
                    _a3c_prewarm_spaceagora_cache!,
                    process_ids,
                    session.config.scenario,
                    session.config.training.max_passes_per_campaign,
                )
                return _train_parallel_a3c_live!(
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
        return _train_parallel_a3c_live!(
            session;
            global_steps=global_steps,
            episodes=episodes,
            n_workers=n_workers,
        )
    end
end
