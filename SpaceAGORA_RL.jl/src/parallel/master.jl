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
                    save_checkpoint(joinpath(session.output_dir, "checkpoint_$(next_checkpoint_step).jls"),
                                    session.learner; manifest=session.manifest)
                    next_checkpoint_step += checkpoint_frequency
                end
            end
        end

        episode += length(batch_episodes)
    end

    save_checkpoint(joinpath(session.output_dir, "checkpoint_final.jls"), session.learner;
                    manifest=session.manifest)
    return (summaries=summaries, transitions=transitions,
            metrics=[episode_metrics(s; policy_name="ddqn") for s in summaries],
            aggregate=aggregate_metrics(summaries; policy_name="ddqn"),
            output_dir=session.output_dir)
end
