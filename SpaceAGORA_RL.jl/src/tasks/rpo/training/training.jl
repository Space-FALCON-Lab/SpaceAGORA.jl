Base.@kwdef struct RPOHyPRRLTrainingConfig
    episodes::Int = 100_000
    seed::Int = 740
    n_workers::Int = 16
    worker_backend::Symbol = :processes
    epsilon_start::Float64 = 1.0
    epsilon_stop::Float64 = 0.01
    epsilon_decay_start_episode::Int = 10_000
    epsilon_decay_end_episode::Int = 20_000
    successful_case_repetitions::Int = 3
    progress_every_episodes::Int = 10
    checkpoint_every_episodes::Int = 1_000
    checkpoint_directory::String = joinpath(pwd(), "hypr_rl_checkpoints")
end

struct RPOHyPRRLWorkerEpisode
    episode_index::Int
    worker_id::Int
    seed::Int
    episode_return::Float64
    best_objective::Float64
    terminal_feasible::Bool
    terminal_reward_correction::Float64
    accepted_edits::Int
    successful_case_repeat_index::Int
    transitions::Vector{MaskedTransition}
end

function rpo_hypr_rl_epsilon_schedule(training::RPOHyPRRLTrainingConfig)
    return EpsilonSchedule(
        start=training.epsilon_start,
        stop=training.epsilon_stop,
        decay_start_step=training.epsilon_decay_start_episode,
        decay_steps=training.epsilon_decay_end_episode -
                    training.epsilon_decay_start_episode,
    )
end

function rpo_hypr_rl_ddqn_config(config::RPOHyPRRLConfig;
                                 hidden_dim::Int=512,
                                 learning_rate::Real=1.0e-4,
                                 discount::Real=0.99,
                                 batch_size::Int=256,
                                 train_start::Int=5_000,
                                 replay_size::Int=500_000,
                                 target_update::Int=5_000)
    return DDQNConfig(
        learning_rate=Float64(learning_rate), discount=Float64(discount),
        batch_size=batch_size, train_start=train_start, replay_size=replay_size,
        target_update=target_update, hidden_dim=hidden_dim,
        obs_dim=rpo_hypr_rl_observation_dim(config),
        action_dim=action_count(config),
    )
end

function save_hypr_rl_checkpoint(path::AbstractString, learner::DDQNLearner,
                                 config::RPOHyPRRLConfig;
                                 training_metadata=nothing)
    return save_checkpoint(
        path, learner;
        task=:rpo_hypr_rl,
        algorithm=:pr_drl,
        action_table=rpo_editor_actions(config),
        task_metadata=(config=config, training=training_metadata),
    )
end

function _hypr_rl_snapshot_action(network::QNetwork, schedule::EpsilonSchedule,
                                  ddqn_config::DDQNConfig, observation,
                                  action_mask::AbstractVector{Bool}, episode::Int,
                                  rng::AbstractRNG)
    epsilon = epsilon_value(schedule, episode)
    if rand(rng) < epsilon
        return rand(rng, _valid_action_indices(action_mask))
    end
    return _masked_argmax(predict_q(network, observation), action_mask)
end

function _run_hypr_rl_worker_episode(mdp::RPOHyPRRLMDP,
                                     schedule::EpsilonSchedule,
                                     ddqn_config::DDQNConfig,
                                     policy_snapshot::QNetwork,
                                     episode_index::Int,
                                     worker_id::Int,
                                     scenario_seed::Int,
                                     action_seed::Int,
                                     successful_case_repeat_index::Int,
                                     terminal_evaluator)
    scenario_rng = MersenneTwister(scenario_seed)
    action_rng = MersenneTwister(action_seed)
    state = reset_scenario(mdp, scenario_rng)
    transitions = MaskedTransition[]
    episode_return = 0.0
    accepted_edits = 0
    while !state.stopped && state.edit_count < mdp.config.max_edits
        observation = observe_state(mdp, state)
        mask = valid_action_mask(mdp, state)
        action = _hypr_rl_snapshot_action(
            policy_snapshot, schedule, ddqn_config, observation, mask,
            episode_index,
            action_rng,
        )
        result = step_scenario(mdp, state, action, scenario_rng)
        next_observation = observe_state(mdp, result.state)
        next_mask = valid_action_mask(mdp, result.state)
        push!(transitions, MaskedTransition(
            observation, mask, action, Float32(result.reward),
            next_observation, next_mask, result.terminated, result.truncated,
            episode_index,
        ))
        episode_return += result.reward
        accepted_edits += result.accepted ? 1 : 0
        state = result.state
        (result.terminated || result.truncated) && break
    end
    terminal_evaluation = terminal_evaluator(
        mdp.scenario, mdp.config, state.best_control_points_rtn,
        state.best_attitude_progress, state.best_attitude_quaternions,
    )
    full_terminal_reward = _rpo_terminal_reward(
        state.seed_evaluation, terminal_evaluation, mdp.config,
    )
    stopped_by_policy = !isempty(transitions) &&
        rpo_editor_actions(mdp.config)[transitions[end].action_index].kind == :stop
    edit_terminal_reward = stopped_by_policy ? _rpo_terminal_reward(
        state.seed_evaluation, state.best_evaluation, mdp.config,
    ) : 0.0
    terminal_correction = full_terminal_reward - edit_terminal_reward
    if !isempty(transitions)
        final_transition = transitions[end]
        transitions[end] = MaskedTransition(
            final_transition.observation,
            final_transition.action_mask,
            final_transition.action_index,
            Float32(final_transition.reward + terminal_correction),
            final_transition.next_observation,
            final_transition.next_action_mask,
            true,
            final_transition.truncated,
            final_transition.info_index,
        )
    end
    episode_return += terminal_correction
    return RPOHyPRRLWorkerEpisode(
        episode_index, worker_id, scenario_seed, episode_return,
        terminal_evaluation.objective, terminal_evaluation.feasible,
        terminal_correction, accepted_edits, successful_case_repeat_index,
        transitions,
    )
end

function _run_hypr_rl_worker_episode(mdp::RPOHyPRRLMDP,
                                     schedule::EpsilonSchedule,
                                     ddqn_config::DDQNConfig,
                                     policy_snapshot::QNetwork,
                                     episode_index::Int,
                                     worker_id::Int,
                                     scenario_seed::Int,
                                     action_seed::Int,
                                     successful_case_repeat_index::Int)
    return _run_hypr_rl_worker_episode(
        mdp, schedule, ddqn_config, policy_snapshot, episode_index, worker_id,
        scenario_seed, action_seed, successful_case_repeat_index, mdp.evaluator,
    )
end

function _prepare_hypr_rl_process_worker!()
    LinearAlgebra.BLAS.set_num_threads(1)
    _spaceagora_rpo_modules()
    return nothing
end

function _setup_hypr_rl_process_workers(n_workers::Int)
    project_file = Base.active_project()
    project_directory = project_file === nothing ? package_root() : dirname(project_file)
    process_ids = addprocs(
        n_workers;
        exeflags=Cmd([
            "--project=$(project_directory)",
            "--threads=1",
            "--compiled-modules=existing",
        ]),
    )
    try
        _remotecall_wait_all(Base.eval, process_ids, Main, :(using SpaceAGORA_RL))
    catch
        rmprocs(process_ids)
        rethrow()
    end
    return process_ids
end

function _hypr_rl_active_worker_count(training::RPOHyPRRLTrainingConfig)
    requested = max(1, training.n_workers)
    if training.worker_backend == :threads
        active = min(requested, Threads.nthreads())
        active < requested && @warn(
            "HyPR-RL n_workers exceeds available Julia threads; using Threads.nthreads()",
            requested_workers=requested,
            active_workers=active,
        )
        return active
    elseif training.worker_backend == :processes
        return requested
    end
    throw(ArgumentError("HyPR-RL worker_backend must be :threads or :processes"))
end

function _hypr_rl_episode_mdp(mdp::RPOHyPRRLMDP, scenario_sampler,
                              scenario_seed::Int, episode::Int)
    scenario_sampler === nothing && return mdp
    scenario_rng = MersenneTwister(scenario_seed)
    return RPOHyPRRLMDP(
        mdp.config, scenario_sampler(scenario_rng, episode);
        evaluator=mdp.evaluator,
    )
end

function _launch_hypr_rl_rollout_batch(mdp::RPOHyPRRLMDP,
                                       scenario_sampler,
                                       training::RPOHyPRRLTrainingConfig,
                                       learner::DDQNLearner,
                                       job_specs,
                                       process_ids::Vector{Int},
                                       terminal_evaluator)
    policy_snapshot = cpu_network(learner.online)
    jobs = map(enumerate(job_specs)) do (local_worker_id, spec)
        episode = spec.episode
        scenario_seed = spec.scenario_seed
        action_seed = training.seed + 1_000_000_007 +
                      10_000 * local_worker_id + episode
        episode_mdp = _hypr_rl_episode_mdp(
            mdp, scenario_sampler, scenario_seed, episode,
        )
        return (
            mdp=episode_mdp,
            episode=episode,
            worker_id=local_worker_id,
            scenario_seed=scenario_seed,
            action_seed=action_seed,
            successful_case_repeat_index=spec.repeat_index,
        )
    end
    if training.worker_backend == :threads
        tasks = map(jobs) do job
            Threads.@spawn _run_hypr_rl_worker_episode(
                $(job.mdp), learner.schedule, learner.config, $policy_snapshot,
                $(job.episode), $(job.worker_id), $(job.scenario_seed),
                $(job.action_seed), $(job.successful_case_repeat_index),
                terminal_evaluator,
            )
        end
        return [fetch(task) for task in tasks]
    end
    futures = map(jobs) do job
        process_id = process_ids[job.worker_id]
        remotecall(
            _run_hypr_rl_worker_episode, process_id,
            job.mdp, learner.schedule, learner.config, policy_snapshot,
            job.episode, job.worker_id, job.scenario_seed, job.action_seed,
            job.successful_case_repeat_index, terminal_evaluator,
        )
    end
    return [fetch(future) for future in futures]
end

function _hypr_rl_eta_seconds(episode::Int, total_episodes::Int, elapsed_s::Real)
    episode <= 0 && return Inf
    rate = episode / max(Float64(elapsed_s), eps(Float64))
    return rate > 0.0 ? max(total_episodes - episode, 0) / rate : Inf
end

function _hypr_rl_training_progress(episode::Int, training::RPOHyPRRLTrainingConfig,
                                    learner::DDQNLearner, replay::MaskedReplayBuffer,
                                    returns, objectives, terminal_feasible,
                                    repeat_indices, active_workers::Int,
                                    start_time::Float64)
    first_index = max(1, episode - 99)
    mean_return = mean(view(returns, first_index:episode))
    finite_objectives = filter(isfinite, view(objectives, first_index:episode))
    mean_objective = isempty(finite_objectives) ? Inf : mean(finite_objectives)
    success_rate = mean(view(terminal_feasible, first_index:episode))
    elapsed = time() - start_time
    eta = _hypr_rl_eta_seconds(episode, training.episodes, elapsed)
    repeat_count = count(>(0), view(repeat_indices, 1:episode))
    @printf(
        "HyPR-RL progress episodes=%d/%d steps=%d replay=%d train_steps=%d loss=%s eps=%.4f mean_return_100=%.5g mean_objective_100=%.5g terminal_success_100=%.3f successful_case_replays=%d workers=%d backend=%s elapsed=%s eta=%s\n",
        episode, training.episodes, learner.global_step, length(replay),
        learner.train_steps,
        isfinite(learner.last_loss) ? @sprintf("%.6g", learner.last_loss) : "n/a",
        epsilon_value(learner.schedule, episode), mean_return,
        mean_objective, success_rate, repeat_count, active_workers,
        string(training.worker_backend), _format_duration(elapsed),
        isfinite(eta) ? _format_duration(eta) : "n/a",
    )
    flush(stdout)
    return nothing
end

"""
Train HyPR-RL with the PR-DRL central-learner/parallel-rollout strategy.

Every rollout batch uses one frozen policy snapshot. Workers own their scenario
and RNG state and return masked transitions; only the central task mutates the
learner, replay buffer, optimizer, and checkpoint state.
"""
function train_hypr_rl!(mdp::RPOHyPRRLMDP;
                        training::RPOHyPRRLTrainingConfig=RPOHyPRRLTrainingConfig(),
                        ddqn_config::DDQNConfig=rpo_hypr_rl_ddqn_config(mdp.config),
                        schedule::EpsilonSchedule=
                            rpo_hypr_rl_epsilon_schedule(training),
                        learner::DDQNLearner=DDQNLearner(
                            MersenneTwister(training.seed), ddqn_config;
                            schedule=schedule,
                        ),
                        scenario_sampler=nothing,
                        terminal_evaluator=mdp.evaluator)
    training.episodes >= 0 || throw(ArgumentError("training episodes must be nonnegative"))
    training.n_workers > 0 || throw(ArgumentError("n_workers must be positive"))
    0.0 <= training.epsilon_stop <= training.epsilon_start <= 1.0 ||
        throw(ArgumentError("epsilon values must satisfy 0 <= stop <= start <= 1"))
    training.epsilon_decay_start_episode >= 0 || throw(ArgumentError(
        "epsilon_decay_start_episode must be nonnegative",
    ))
    training.epsilon_decay_end_episode >= training.epsilon_decay_start_episode ||
        throw(ArgumentError(
            "epsilon_decay_end_episode must not precede epsilon_decay_start_episode",
        ))
    training.successful_case_repetitions >= 0 || throw(ArgumentError(
        "successful_case_repetitions must be nonnegative",
    ))
    training.worker_backend in (:threads, :processes) ||
        throw(ArgumentError("worker_backend must be :threads or :processes"))
    ddqn_config.obs_dim == rpo_hypr_rl_observation_dim(mdp.config) ||
        throw(DimensionMismatch("DDQN observation dimension does not match HyPR-RL"))
    ddqn_config.action_dim == action_count(mdp.config) ||
        throw(DimensionMismatch("DDQN action dimension does not match HyPR-RL"))
    learner.config.obs_dim == ddqn_config.obs_dim ||
        throw(DimensionMismatch("learner observation dimension does not match training config"))
    learner.config.action_dim == ddqn_config.action_dim ||
        throw(DimensionMismatch("learner action dimension does not match training config"))

    active_workers = _hypr_rl_active_worker_count(training)
    replay = MaskedReplayBuffer(
        ddqn_config.obs_dim, ddqn_config.action_dim, ddqn_config.replay_size,
    )
    episode_returns = zeros(Float64, training.episodes)
    best_objectives = fill(Inf, training.episodes)
    terminal_feasible = falses(training.episodes)
    terminal_reward_corrections = zeros(Float64, training.episodes)
    accepted_edits = zeros(Int, training.episodes)
    scenario_seeds = zeros(Int, training.episodes)
    successful_case_repeat_indices = zeros(Int, training.episodes)
    central_rng = MersenneTwister(training.seed)
    process_ids = Int[]
    start_time = time()
    @printf(
        "starting HyPR-RL PR-DRL training episodes=%d requested_workers=%d active_workers=%d worker_backend=%s julia_threads=%d curve=bezier max_bezier_waypoints=%d obs_dim=%d action_dim=%d train_start=%d batch_size=%d epsilon=episode_linear:%.3f@%d->%.3f@%d successful_case_repetitions=%d eta=pending\n",
        training.episodes, training.n_workers, active_workers,
        string(training.worker_backend), Threads.nthreads(),
        mdp.config.max_translation_waypoints, ddqn_config.obs_dim,
        ddqn_config.action_dim, ddqn_config.train_start, ddqn_config.batch_size,
        training.epsilon_start, training.epsilon_decay_start_episode,
        training.epsilon_stop, training.epsilon_decay_end_episode,
        training.successful_case_repetitions,
    )
    flush(stdout)

    try
        if training.worker_backend == :processes && training.episodes > 0
            process_ids = _setup_hypr_rl_process_workers(active_workers)
            _remotecall_wait_all(_prepare_hypr_rl_process_worker!, process_ids)
        end
        next_checkpoint_episode = training.checkpoint_every_episodes > 0 ?
            training.checkpoint_every_episodes : typemax(Int)
        next_progress_episode = training.progress_every_episodes > 0 ?
            training.progress_every_episodes : typemax(Int)
        next_episode = 1
        next_base_case = 1
        pending_replays = NamedTuple{(:seed, :repeat_index),
                                     Tuple{Int, Int}}[]
        while next_episode <= training.episodes
            jobs_in_batch = min(
                active_workers, training.episodes - next_episode + 1,
            )
            job_specs = NamedTuple[]
            for _ in 1:jobs_in_batch
                scenario_spec = if isempty(pending_replays)
                    seed = training.seed + next_base_case
                    next_base_case += 1
                    (seed=seed, repeat_index=0)
                else
                    popfirst!(pending_replays)
                end
                push!(job_specs, (
                    episode=next_episode,
                    scenario_seed=scenario_spec.seed,
                    repeat_index=scenario_spec.repeat_index,
                ))
                next_episode += 1
            end
            worker_results = _launch_hypr_rl_rollout_batch(
                mdp, scenario_sampler, training, learner, job_specs,
                process_ids, terminal_evaluator,
            )
            sort!(worker_results; by=result -> result.episode_index)
            for worker_result in worker_results
                episode = worker_result.episode_index
                for transition in worker_result.transitions
                    learner.global_step += 1
                    push!(replay, transition)
                    maybe_train!(learner, replay, central_rng)
                end
                episode_returns[episode] = worker_result.episode_return
                best_objectives[episode] = worker_result.best_objective
                terminal_feasible[episode] = worker_result.terminal_feasible
                terminal_reward_corrections[episode] =
                    worker_result.terminal_reward_correction
                accepted_edits[episode] = worker_result.accepted_edits
                scenario_seeds[episode] = worker_result.seed
                successful_case_repeat_indices[episode] =
                    worker_result.successful_case_repeat_index
                repeat_plan = _next_successful_case_repeat(
                    worker_result.seed,
                    worker_result.successful_case_repeat_index,
                    worker_result.terminal_feasible,
                    training.successful_case_repetitions,
                )
                repeat_plan === nothing || push!(pending_replays, repeat_plan)
                if episode >= next_checkpoint_episode
                    checkpoint_path = joinpath(
                        training.checkpoint_directory,
                        "hypr_rl_episode_$(lpad(next_checkpoint_episode, 6, '0')).jls",
                    )
                    save_hypr_rl_checkpoint(
                        checkpoint_path, learner, mdp.config;
                        training_metadata=(
                            episode=episode,
                            worker_backend=training.worker_backend,
                            n_workers=active_workers,
                            successful_case_repetitions=
                                training.successful_case_repetitions,
                            epsilon_decay_start_episode=
                                training.epsilon_decay_start_episode,
                            epsilon_decay_end_episode=
                                training.epsilon_decay_end_episode,
                        ),
                    )
                    next_checkpoint_episode += training.checkpoint_every_episodes
                end
            end
            completed_episode = worker_results[end].episode_index
            if completed_episode >= next_progress_episode
                _hypr_rl_training_progress(
                    completed_episode, training, learner, replay, episode_returns,
                    best_objectives, terminal_feasible,
                    successful_case_repeat_indices, active_workers, start_time,
                )
                while next_progress_episode <= completed_episode
                    next_progress_episode += training.progress_every_episodes
                end
            end
        end
        if training.episodes > 0 &&
           (training.progress_every_episodes <= 0 ||
            mod(training.episodes, training.progress_every_episodes) != 0)
            _hypr_rl_training_progress(
                training.episodes, training, learner, replay, episode_returns,
                best_objectives, terminal_feasible,
                successful_case_repeat_indices, active_workers, start_time,
            )
        end
    finally
        isempty(process_ids) || rmprocs(process_ids)
    end
    return (
        learner=learner,
        replay=replay,
        episode_returns=episode_returns,
        best_objectives=best_objectives,
        terminal_feasible=terminal_feasible,
        terminal_reward_corrections=terminal_reward_corrections,
        accepted_edits=accepted_edits,
        scenario_seeds=scenario_seeds,
        successful_case_repeat_indices=successful_case_repeat_indices,
        active_workers=active_workers,
        worker_backend=training.worker_backend,
    )
end

train_hypr_rl_parallel!(args...; kwargs...) = train_hypr_rl!(args...; kwargs...)
