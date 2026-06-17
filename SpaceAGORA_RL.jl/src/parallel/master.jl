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

function _budget_label(budget::Int)
    return budget == typemax(Int) ? "none" : string(budget)
end

function _print_training_progress(session::TrainingSession, summaries::Vector{EpisodeSummary},
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
    stats = _recent_training_stats(summaries, 100)
    eps = epsilon_value(session.learner.schedule, session.learner.global_step)
    loss = isfinite(session.learner.last_loss) ? @sprintf("%.6g", session.learner.last_loss) : "n/a"
    if step_limited
        @printf(
            "progress ep=%d steps=%d/%d replay=%d train_steps=%d loss=%s eps=%.4f recent_reward=%.3f recent_success=%.1f%% recent_passes=%.1f workers=%d/%d elapsed=%s eta=%s\n",
            completed_episodes,
            session.learner.global_step,
            target_global_step,
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
    else
        @printf(
            "progress ep=%d/%s steps=%d replay=%d train_steps=%d loss=%s eps=%.4f recent_reward=%.3f recent_success=%.1f%% recent_passes=%.1f workers=%d/%d elapsed=%s eta=%s\n",
            completed_episodes,
            _budget_label(episode_budget),
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
    end
    flush(stdout)
end

function train_parallel!(session::TrainingSession;
                         global_steps::Int=session.config.training.global_steps,
                         episodes::Union{Nothing,Int}=nothing,
                         n_workers::Int=session.config.training.n_workers)
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

    @printf(
        "starting DDQN training global_steps=%s episode_cap=%s n_workers=%d active_workers=%d julia_threads=%d train_start=%d batch_size=%d checkpoint_frequency=%d progress_frequency=%d\n",
        target_global_step == typemax(Int) ? "none" : string(target_global_step),
        _budget_label(episode_budget),
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
    while session.learner.global_step < target_global_step && episode <= episode_budget
        batch_episodes = episode:min(episode_budget, episode + active_workers - 1)
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
            if progress_frequency > 0 && target_global_step != typemax(Int) &&
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
            metrics=[episode_metrics(s; policy_name="ddqn") for s in summaries],
            aggregate=aggregate_metrics(summaries; policy_name="ddqn"),
            global_step=session.learner.global_step,
            target_global_step=target_global_step == typemax(Int) ? nothing : target_global_step,
            output_dir=session.output_dir)
end
