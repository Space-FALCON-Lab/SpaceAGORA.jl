function _format_duration(seconds::Real)
    total = max(0, round(Int, Float64(seconds)))
    hours = total ÷ 3600
    minutes = (total % 3600) ÷ 60
    secs = total % 60
    return @sprintf("%02d:%02d:%02d", hours, minutes, secs)
end

function _recent_training_stats(summaries::Vector{EpisodeSummary}, window::Int)
    isempty(summaries) && return (mean_reward=NaN, success_rate=NaN, mean_passes=NaN)
    first_idx = max(1, length(summaries) - max(1, window) + 1)
    recent = summaries[first_idx:end]
    return (
        mean_reward = mean(getfield.(recent, :episode_reward)),
        success_rate = count(summary -> summary.success, recent) / length(recent),
        mean_passes = mean(getfield.(recent, :pass_count)),
    )
end

function _print_training_progress(session::TrainingSession, summaries::Vector{EpisodeSummary},
                                  episodes::Int, active_workers::Int,
                                  start_time::Float64, progress_window::Int)
    done = length(summaries)
    elapsed = time() - start_time
    episode_rate = elapsed > 0 ? done / elapsed : 0.0
    eta = episode_rate > 0 ? (episodes - done) / episode_rate : Inf
    stats = _recent_training_stats(summaries, progress_window)
    eps = epsilon_value(session.learner.schedule, session.learner.global_step)
    loss = isfinite(session.learner.last_loss) ? @sprintf("%.6g", session.learner.last_loss) : "n/a"
    @printf(
        "progress ep=%d/%d steps=%d replay=%d train_steps=%d loss=%s eps=%.4f recent_reward=%.3f recent_success=%.1f%% recent_passes=%.1f workers=%d/%d elapsed=%s eta=%s\n",
        done,
        episodes,
        session.learner.global_step,
        length(session.learner.replay),
        session.learner.train_steps,
        loss,
        eps,
        stats.mean_reward,
        100 * stats.success_rate,
        stats.mean_passes,
        active_workers,
        Threads.nthreads(),
        _format_duration(elapsed),
        isfinite(eta) ? _format_duration(eta) : "n/a",
    )
    flush(stdout)
end

function train_parallel!(session::TrainingSession;
                         episodes::Int=session.config.training.episodes,
                         n_workers::Int=session.config.training.n_workers)
    summaries = EpisodeSummary[]
    transitions = Transition[]
    requested_workers = max(1, n_workers)
    active_workers = min(requested_workers, Threads.nthreads())
    active_workers < requested_workers && @warn "n_workers exceeds available Julia threads; using Threads.nthreads()" requested_workers active_workers
    checkpoint_frequency = session.config.training.checkpoint_frequency
    next_checkpoint_step = checkpoint_frequency > 0 ? checkpoint_frequency : typemax(Int)
    progress_frequency = max(0, session.config.training.progress_frequency)
    next_progress_episode = progress_frequency > 0 ? min(progress_frequency, episodes) : typemax(Int)
    start_time = time()

    @printf(
        "starting DDQN training episodes=%d n_workers=%d active_workers=%d julia_threads=%d train_start=%d batch_size=%d checkpoint_frequency=%d progress_frequency=%d\n",
        episodes,
        requested_workers,
        active_workers,
        Threads.nthreads(),
        session.learner.config.train_start,
        session.learner.config.batch_size,
        checkpoint_frequency,
        progress_frequency,
    )
    @printf("output_dir=%s\n", session.output_dir)
    flush(stdout)

    episode = 1
    while episode <= episodes
        batch_episodes = episode:min(episodes, episode + active_workers - 1)
        policy_snapshot = copy(session.learner.online)
        global_step_start = session.learner.global_step
        tasks = map(enumerate(batch_episodes)) do (local_worker_id, episode_index)
            seed = session.config.training.seed + 10_000 * local_worker_id + episode_index
            Threads.@spawn run_threaded_worker_episode(
                session.config.scenario,
                session.learner.schedule,
                session.learner.config,
                $policy_snapshot,
                $episode_index,
                $local_worker_id,
                $seed,
                session.config.training.max_steps,
                global_step_start;
                train=true,
            )
        end

        for task in tasks
            summary, episode_transitions = fetch(task)
            push!(summaries, summary)
            append!(transitions, episode_transitions)
            for transition in episode_transitions
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
        end

        if progress_frequency > 0 && length(summaries) >= next_progress_episode
            _print_training_progress(session, summaries, episodes, active_workers,
                                     start_time, progress_frequency)
            while next_progress_episode <= length(summaries)
                next_progress_episode += progress_frequency
            end
        end

        episode += length(batch_episodes)
    end

    if isempty(summaries) || progress_frequency == 0 || length(summaries) < next_progress_episode - progress_frequency
        _print_training_progress(session, summaries, episodes, active_workers,
                                 start_time, max(1, progress_frequency))
    end

    final_checkpoint_path = joinpath(session.output_dir, "checkpoint_final.jls")
    save_checkpoint(final_checkpoint_path, session.learner; manifest=session.manifest)
    @printf("training complete final_checkpoint=%s elapsed=%s\n",
            final_checkpoint_path, _format_duration(time() - start_time))
    flush(stdout)
    return (summaries=summaries, transitions=transitions,
            metrics=[episode_metrics(s; policy_name="ddqn") for s in summaries],
            aggregate=aggregate_metrics(summaries; policy_name="ddqn"),
            output_dir=session.output_dir)
end
