#!/usr/bin/env julia

"""
Evaluate a completed SpaceAGORA_RL run using the protocol and artifacts reported
in Falcone and Putnam (2023), "Autonomous Decision-Making for Aerobraking via
Parallel Randomized Deep Reinforcement Learning."

The script accepts a run directory, discovers its manifest/config/checkpoints,
and writes raw episode/pass data, paper tables, and paper-style figures.
Run with `--help` for the command-line interface.
"""

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using CSV
using DataFrames
using Distributed
using Plots
using Printf
using SHA
using Statistics
using TOML
using SpaceAGORA_RL

const ODYSSEY_REFERENCE = (
    maneuver_count = 15.0,
    mission_duration_days = 23.3,
    total_mission_delta_v_mps = 31.78,
    thermal_violations = 9.0,
    pass_count = 210.0,
)
const PAPER_TARGET_TOLERANCE_KM = 10.0

const PAPER_GENERALIZATION_CASES = (
    "nominal",
    "exponential_density",
    "aggressive_atmosphere",
    "short_campaign",
    "long_campaign",
    "accurate_simulator",
)

function _usage(io::IO=stdout)
    println(io, """
Usage:
  julia --project=SpaceAGORA_RL.jl \\
    SpaceAGORA_RL.jl/scripts/evaluate_rl_run.jl RUN_DIR [options]

Required:
  RUN_DIR                     RL run containing manifest.toml and checkpoints.

Options:
  --output DIR                Output directory (default: RUN_DIR/paper_evaluation,
                              or RUN_DIR/final_flight_comparison in comparison mode).
  --config PATH               Override the config recorded in the run manifest.
  --checkpoint PATH|final     Policy checkpoint (default: checkpoint_final.jls,
                              otherwise the largest numbered checkpoint).
  --final-flight-comparison   Only evaluate checkpoint_final.jls against AADS and
                              the Mars Odyssey mission reference. Uses the
                              thermal-tolerant protocol and --episodes campaigns
                              (default: 40), including the percentage that finish
                              within ±10 km of the target apoapsis radius.
  --episodes N                IID PR-DRL/AADS episodes for Table V/Fig. 10
                              (default: 40).
  --flight-episodes N         Odyssey-geometry episodes for Figs. 11-12
                              (default: 40).
  --generalization-episodes N Episodes per Table VI case (default: 100).
  --checkpoint-episodes N     Episodes per mode and checkpoint for Figs. 7-8
                              (default: 40).
  --checkpoint-stride N       Evaluate every Nth numbered checkpoint (default: 1).
  --processes N               Worker processes for final-flight comparison
                              (default: min(4, available CPU threads)).
  --threads-per-process N     Julia threads available inside each worker process
                              (default: 1).
  --progress-every N          Print progress every N completed campaigns
                              (default: 1).
  --skip-checkpoint-sweep     Skip Figs. 7-8 checkpoint evaluation.
  --skip-generalization       Skip Table VI generalization evaluation.
  --help                      Show this message.

Outputs include the raw episode/pass CSVs, checkpoint_metrics.csv,
paper_table_v_pr_drl_vs_aads.csv, paper_table_vi_generalization.csv,
paper_fig07_*.png through paper_fig12_*.png, the dedicated episode-completion
and final-target-distance chart, and evaluation_manifest.toml.
Final-flight-comparison mode writes its four-panel PNG, comparison CSV, raw
episode/pass CSVs, and evaluation_manifest.toml only.
""")
end

function _parse_cli(args)
    isempty(args) && throw(ArgumentError("RUN_DIR is required; use --help for usage"))
    any(==("--help"), args) && return (help=true,)
    startswith(first(args), "--") &&
        throw(ArgumentError("RUN_DIR must be the first argument"))

    options = Dict{Symbol,Any}(
        :run_dir => first(args),
        :output => nothing,
        :config => nothing,
        :checkpoint => nothing,
        :final_flight_comparison => false,
        :episodes => PAPER_IID_EVALUATION_EPISODES,
        :flight_episodes => PAPER_IID_EVALUATION_EPISODES,
        :generalization_episodes => PAPER_GENERALIZATION_EVALUATION_EPISODES,
        :checkpoint_episodes => PAPER_IID_EVALUATION_EPISODES,
        :checkpoint_stride => 1,
        :processes => max(1, min(4, Sys.CPU_THREADS)),
        :threads_per_process => 1,
        :progress_every => 1,
        :checkpoint_sweep => true,
        :generalization => true,
    )
    value_options = Dict(
        "--output" => :output,
        "--config" => :config,
        "--checkpoint" => :checkpoint,
        "--episodes" => :episodes,
        "--flight-episodes" => :flight_episodes,
        "--generalization-episodes" => :generalization_episodes,
        "--checkpoint-episodes" => :checkpoint_episodes,
        "--checkpoint-stride" => :checkpoint_stride,
        "--processes" => :processes,
        "--threads-per-process" => :threads_per_process,
        "--progress-every" => :progress_every,
    )
    integer_options = Set([
        :episodes,
        :flight_episodes,
        :generalization_episodes,
        :checkpoint_episodes,
        :checkpoint_stride,
        :processes,
        :threads_per_process,
        :progress_every,
    ])

    index = 2
    while index <= length(args)
        arg = args[index]
        if arg == "--skip-checkpoint-sweep"
            options[:checkpoint_sweep] = false
            index += 1
        elseif arg == "--final-flight-comparison"
            options[:final_flight_comparison] = true
            index += 1
        elseif arg == "--skip-generalization"
            options[:generalization] = false
            index += 1
        elseif haskey(value_options, arg)
            index == length(args) &&
                throw(ArgumentError("missing value for $arg"))
            key = value_options[arg]
            value = args[index + 1]
            options[key] = key in integer_options ? parse(Int, value) : value
            index += 2
        else
            throw(ArgumentError("unknown option: $arg"))
        end
    end
    for key in integer_options
        options[key] > 0 || throw(ArgumentError("$key must be positive"))
    end
    return (; help=false, options...)
end

function _resolve_existing_path(path::AbstractString, run_dir::AbstractString)
    candidates = unique([
        abspath(path),
        abspath(joinpath(run_dir, path)),
        abspath(joinpath(dirname(run_dir), path)),
        abspath(joinpath(dirname(dirname(run_dir)), path)),
        abspath(joinpath(dirname(dirname(dirname(run_dir))), path)),
    ])
    found = findfirst(isfile, candidates)
    return found === nothing ? nothing : candidates[found]
end

function _config_path(run_dir::AbstractString, manifest::Dict, override)
    if override !== nothing
        path = _resolve_existing_path(String(override), run_dir)
        path === nothing && throw(ArgumentError("config does not exist: $override"))
        return path
    end
    recorded = strip(String(get(manifest, "config_path", "")))
    if !isempty(recorded)
        path = _resolve_existing_path(recorded, run_dir)
        path !== nothing && return path
        @warn "manifest config_path was not found; using paper replication defaults" recorded
    end
    return default_config_path()
end

function _numeric_checkpoint_step(path::AbstractString)
    match_result = match(r"^checkpoint_(\d+)\.jls$", basename(path))
    return match_result === nothing ? nothing : parse(Int, match_result.captures[1])
end

function _checkpoint_path(run_dir::AbstractString, requested)
    if requested !== nothing && requested != "final"
        path = _resolve_existing_path(String(requested), run_dir)
        path === nothing && throw(ArgumentError("checkpoint does not exist: $requested"))
        return path
    end
    final_path = joinpath(run_dir, "checkpoint_final.jls")
    if (requested === nothing || requested == "final") && isfile(final_path)
        return final_path
    end
    numbered = filter(path -> _numeric_checkpoint_step(path) !== nothing,
                      frozen_checkpoint_paths(run_dir))
    isempty(numbered) && throw(ArgumentError("no usable checkpoint found in $run_dir"))
    return last(numbered)
end

function _final_checkpoint_path(run_dir::AbstractString)
    path = joinpath(run_dir, "checkpoint_final.jls")
    isfile(path) ||
        throw(ArgumentError("final-policy comparison requires $path"))
    return path
end

function _load_policy(checkpoint_path::AbstractString)
    payload = load_checkpoint(checkpoint_path)
    algorithm = Symbol(get(payload, :algorithm, haskey(payload, :actor) ? :a2c : :pr_drl))
    policy = algorithm == :a2c ?
             load_trained_a2c_policy(checkpoint_path) :
             load_trained_pr_drl_policy(checkpoint_path)
    return policy, algorithm, payload
end

_std(values) = length(values) > 1 ? std(values) : 0.0

function _format_wall_time(seconds::Real)
    total = max(0, round(Int, seconds))
    hours, remainder = divrem(total, 3600)
    minutes, secs = divrem(remainder, 60)
    return hours > 0 ? @sprintf("%02d:%02d:%02d", hours, minutes, secs) :
           @sprintf("%02d:%02d", minutes, secs)
end

function _start_evaluation_workers(
    processes::Int,
    threads_per_process::Int,
    scenario::AerobrakingScenarioConfig,
)
    project_file = Base.active_project()
    project_file === nothing &&
        error("cannot start evaluation workers without an active Julia project")
    project_dir = dirname(project_file)
    @printf(
        "starting evaluation workers processes=%d threads_per_process=%d project=%s\n",
        processes,
        threads_per_process,
        project_dir,
    )
    flush(stdout)
    gram_once_per_step =
        SpaceAGORA_RL._spaceagora_physics_gram_once_per_step(scenario)
    process_ids = addprocs(
        processes;
        exeflags=Cmd([
            "--project=$(project_dir)",
            "--threads=$(threads_per_process)",
        ]),
    )
    try
        @sync for process_id in process_ids
            @async begin
                @printf("evaluation worker initializing pid=%d\n", process_id)
                flush(stdout)
                remotecall_wait(
                    process_id,
                    threads_per_process,
                    length(process_ids),
                    gram_once_per_step,
                ) do thread_count, process_count, use_gram_once_per_step
                    ENV["SPACEAGORA_OUTER_PARALLEL_ACTIVE"] =
                        process_count > 1 ? "1" : "0"
                    ENV["SPACEAGORA_INNER_THREAD_BUDGET"] = string(thread_count)
                    ENV["SPACEAGORA_GRAM_ONCE_PER_STEP"] =
                        use_gram_once_per_step ? "1" : "0"
                    return nothing
                end
                remotecall_wait(
                    Base.eval,
                    process_id,
                    Main,
                    :(using Distributed; using SpaceAGORA_RL),
                )
                actual_threads = remotecall_fetch(process_id) do
                    Threads.nthreads()
                end
                @printf(
                    "evaluation worker ready pid=%d threads=%d gram_once_per_step=%s\n",
                    process_id,
                    actual_threads,
                    gram_once_per_step,
                )
                flush(stdout)
            end
        end
    catch
        rmprocs(process_ids)
        rethrow()
    end
    return process_ids
end

function _result_from_summaries(summaries, scenario, policy_name::AbstractString)
    sort!(summaries; by=summary -> summary.episode_index)
    pass_rows = NamedTuple[]
    for summary in summaries
        append!(
            pass_rows,
            SpaceAGORA_RL.pass_log_rows(summary; policy_name=policy_name),
        )
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
    return (
        summaries=summaries,
        transitions=Transition[],
        pass_rows=pass_rows,
        metrics=metrics,
        aggregate=aggregate,
    )
end

function _evaluate_pair_parallel(
    policy,
    scenario,
    episodes::Int,
    seed::Int,
    protected_initialization,
    process_ids::Vector{Int};
    paper_protocol::Bool=true,
    progress_every::Int=1,
)
    isempty(process_ids) &&
        throw(ArgumentError("parallel evaluation requires at least one worker process"))
    policies = Dict{String,Any}(
        "aads_heuristic" => AADSHeuristicPolicy(),
        "trained_pr_drl" => policy,
    )
    jobs = [
        (
            policy_name=policy_name,
            episode=episode,
            campaign_seed=seed + episode - 1,
        )
        for policy_name in ("aads_heuristic", "trained_pr_drl")
        for episode in 1:episodes
    ]
    total = length(jobs)
    progress = RemoteChannel(() -> Channel{Any}(total), myid())
    pool = CachingPool(process_ids)
    start_time = time()
    @printf(
        "campaign evaluation starting total=%d policies=2 campaigns_per_policy=%d workers=%d\n",
        total,
        episodes,
        length(process_ids),
    )
    flush(stdout)

    evaluation_task = @async pmap(pool, jobs; batch_size=1) do job
        evaluation_policy = policies[job.policy_name]
        summary = if SpaceAGORA_RL._is_spaceagora_live_backend(scenario.backend_mode)
            campaign_summary, _ =
                SpaceAGORA_RL.run_spaceagora_physics_policy_campaign_episode(
                    evaluation_policy,
                    scenario,
                    job.episode,
                    Distributed.myid(),
                    job.campaign_seed;
                    max_passes_per_campaign=scenario.termination_config.max_passes,
                    protected_initialization=protected_initialization,
                )
            campaign_summary
        else
            episode_result = SpaceAGORA_RL.evaluate_policy(
                evaluation_policy,
                scenario;
                episodes=1,
                seed=job.campaign_seed,
                policy_name=job.policy_name,
                paper_protocol=paper_protocol,
                protected_initialization=protected_initialization,
            )
            only(episode_result.summaries)
        end
        summary.episode_index = job.episode
        summary.worker_id = Distributed.myid()
        put!(
            progress,
            (
                policy_name=job.policy_name,
                episode=job.episode,
                worker_id=summary.worker_id,
                pass_count=summary.pass_count,
                thermal_violations=summary.thermal_violations,
                target_error_km=abs(summary.target_error_m) / 1000,
            ),
        )
        return (policy_name=job.policy_name, summary=summary)
    end

    completed = 0
    last_status_time = time()
    while completed < total
        if isready(progress)
            message = take!(progress)
            completed += 1
            last_status_time = time()
            if completed % progress_every == 0 || completed == total
                elapsed = time() - start_time
                eta = completed == 0 ? Inf :
                      elapsed * (total - completed) / completed
                @printf(
                    "evaluation progress=%d/%d (%.1f%%) policy=%s campaign=%d/%d worker=%d passes=%d thermal_violations=%d end_distance_km=%.3f elapsed=%s eta=%s\n",
                    completed,
                    total,
                    100 * completed / total,
                    message.policy_name,
                    message.episode,
                    episodes,
                    message.worker_id,
                    message.pass_count,
                    message.thermal_violations,
                    message.target_error_km,
                    _format_wall_time(elapsed),
                    _format_wall_time(eta),
                )
                flush(stdout)
            end
        elseif istaskdone(evaluation_task)
            break
        else
            if time() - last_status_time >= 30
                @printf(
                    "evaluation active completed=%d/%d workers=%d elapsed=%s waiting_for_campaigns=true\n",
                    completed,
                    total,
                    length(process_ids),
                    _format_wall_time(time() - start_time),
                )
                flush(stdout)
                last_status_time = time()
            end
            sleep(0.1)
        end
    end
    evaluated = fetch(evaluation_task)
    summaries = Dict(
        "aads_heuristic" => EpisodeSummary[],
        "trained_pr_drl" => EpisodeSummary[],
    )
    for item in evaluated
        push!(summaries[item.policy_name], item.summary)
    end
    return Dict(
        policy_name =>
            _result_from_summaries(summaries[policy_name], scenario, policy_name)
        for policy_name in ("aads_heuristic", "trained_pr_drl")
    )
end

function _evaluate_pair(policy, scenario, episodes, seed, protected_initialization;
                        paper_protocol::Bool=true)
    return Dict(
        "aads_heuristic" => evaluate_policy(
            AADSHeuristicPolicy(),
            scenario;
            episodes=episodes,
            seed=seed,
            policy_name="aads_heuristic",
            paper_protocol=paper_protocol,
            protected_initialization=protected_initialization,
        ),
        "trained_pr_drl" => evaluate_policy(
            policy,
            scenario;
            episodes=episodes,
            seed=seed,
            policy_name="trained_pr_drl",
            paper_protocol=paper_protocol,
            protected_initialization=protected_initialization,
        ),
    )
end

function _result_frames(results::Dict)
    metrics = DataFrame(vcat([result.metrics for result in values(results)]...))
    pass_logs = DataFrame(vcat([result.pass_rows for result in values(results)]...))
    return metrics, pass_logs
end

function _table_v(metrics::DataFrame)
    rows = NamedTuple[]
    for policy in ("trained_pr_drl", "aads_heuristic")
        data = metrics[metrics.policy .== policy, :]
        nrow(data) == 0 && continue
        errors = abs.(Float64.(data.target_error_km))
        push!(rows, (
            policy = policy == "trained_pr_drl" ? "PR-DRL" : "AADS",
            episodes = nrow(data),
            mean_reward = mean(data.episode_reward),
            std_reward = _std(data.episode_reward),
            mean_low_thermal_violations = mean(data.low_thermal_violations),
            mean_soft_thermal_violations = mean(data.soft_thermal_violations),
            mean_medium_thermal_violations = mean(data.medium_thermal_violations),
            mean_hard_thermal_violations = mean(data.hard_thermal_violations),
            mean_total_thermal_violations = mean(data.thermal_violations),
            goal_within_10km_percent = 100 * count(<=(10.0), errors) / nrow(data),
            goal_within_5km_percent = 100 * count(<=(5.0), errors) / nrow(data),
            mean_abs_target_error_km = mean(errors),
            std_abs_target_error_km = _std(errors),
        ))
    end
    return DataFrame(rows)
end

function _policy_label(policy)
    policy == "trained_pr_drl" && return "PR-DRL"
    policy == "aads_heuristic" && return "AADS"
    return replace(String(policy), "_" => " ")
end

function _paper_figure_10(pass_logs::DataFrame, scenario, path::AbstractString)
    heat = plot(
        xlabel="Atmospheric passage",
        ylabel="Heat rate (W/cm^2)",
        title="(a) Heat rate",
        legend=:outerright,
    )
    action = plot(
        xlabel="Atmospheric passage",
        ylabel="ABM Delta-V (m/s)",
        title="(b) Apoapsis maneuver",
        legend=:outerright,
    )
    hline!(heat, [scenario.reward_config.heat_low_w_cm2];
           label="low corridor", linestyle=:dash, color=:gray)
    hline!(heat, [scenario.reward_config.heat_high_w_cm2];
           label="high corridor", linestyle=:dash, color=:red)
    colors = Dict("trained_pr_drl" => :steelblue, "aads_heuristic" => :darkorange)
    for policy in ("trained_pr_drl", "aads_heuristic")
        policy_rows = pass_logs[pass_logs.policy .== policy, :]
        episodes = first(sort(unique(policy_rows.episode)),
                         min(10, length(unique(policy_rows.episode))))
        first_trace = true
        for episode in episodes
            trace = policy_rows[policy_rows.episode .== episode, :]
            label = first_trace ? _policy_label(policy) : ""
            plot!(heat, trace.pass, trace.heat_rate_w_cm2;
                  label=label, color=colors[policy], alpha=0.45)
            plot!(action, trace.pass, trace.action_delta_v_mps;
                  label=label, color=colors[policy], alpha=0.45, seriestype=:steppre)
            first_trace = false
        end
    end
    savefig(plot(heat, action; layout=(2, 1), size=(960, 760)), path)
    return path
end

function _metric_stats(metrics::DataFrame, policy::String, field::Symbol)
    rows = metrics[metrics.policy .== policy, :]
    values = Float64.(rows[!, field])
    return mean(values), _std(values)
end

function _target_reached_stats(metrics::DataFrame, policy::String;
                               tolerance_km::Real=PAPER_TARGET_TOLERANCE_KM)
    rows = metrics[metrics.policy .== policy, :]
    episodes = nrow(rows)
    episodes > 0 || return (
        count=0,
        percent=NaN,
        mean_abs_error_km=NaN,
        std_abs_error_km=NaN,
    )
    target_errors_km = abs.(Float64.(rows.target_error_km))
    reached = count(<=(Float64(tolerance_km)), target_errors_km)
    return (
        count=reached,
        percent=100 * reached / episodes,
        mean_abs_error_km=mean(target_errors_km),
        std_abs_error_km=_std(target_errors_km),
    )
end

function _flight_performance_table(metrics::DataFrame)
    rows = NamedTuple[]
    for (policy, label) in (
        ("trained_pr_drl", "PR-DRL"),
        ("aads_heuristic", "AADS"),
    )
        maneuvers = _metric_stats(metrics, policy, :maneuver_count)
        duration = _metric_stats(metrics, policy, :mission_duration_days)
        delta_v = _metric_stats(metrics, policy, :total_mission_delta_v_mps)
        thermal = _metric_stats(metrics, policy, :thermal_violations)
        target_reached = _target_reached_stats(metrics, policy)
        push!(rows, (
            policy=label,
            episodes=count(==(policy), metrics.policy),
            mean_maneuver_count=maneuvers[1],
            std_maneuver_count=maneuvers[2],
            mean_duration_days=duration[1],
            std_duration_days=duration[2],
            mean_total_delta_v_mps=delta_v[1],
            std_total_delta_v_mps=delta_v[2],
            mean_thermal_violations=thermal[1],
            std_thermal_violations=thermal[2],
            target_reached_10km_count=target_reached.count,
            target_reached_10km_percent=target_reached.percent,
            mean_abs_final_target_error_km=target_reached.mean_abs_error_km,
            std_abs_final_target_error_km=target_reached.std_abs_error_km,
        ))
    end
    push!(rows, (
        policy="Mars Odyssey",
        episodes=1,
        mean_maneuver_count=ODYSSEY_REFERENCE.maneuver_count,
        std_maneuver_count=0.0,
        mean_duration_days=ODYSSEY_REFERENCE.mission_duration_days,
        std_duration_days=0.0,
        mean_total_delta_v_mps=ODYSSEY_REFERENCE.total_mission_delta_v_mps,
        std_total_delta_v_mps=0.0,
        mean_thermal_violations=ODYSSEY_REFERENCE.thermal_violations,
        std_thermal_violations=0.0,
        target_reached_10km_count=missing,
        target_reached_10km_percent=NaN,
        mean_abs_final_target_error_km=NaN,
        std_abs_final_target_error_km=NaN,
    ))
    return DataFrame(rows)
end

function _paper_figure_11(metrics::DataFrame, path::AbstractString)
    policies = ["trained_pr_drl", "aads_heuristic", "odyssey_flight"]
    labels = ["PR-DRL", "AADS", "Mars Odyssey"]
    specs = [
        (
            :maneuver_count,
            "Number of maneuvers",
            ODYSSEY_REFERENCE.maneuver_count,
            "(a) Number of maneuvers",
        ),
        (
            :mission_duration_days,
            "Duration (days)",
            ODYSSEY_REFERENCE.mission_duration_days,
            "(b) Campaign duration",
        ),
        (
            :total_mission_delta_v_mps,
            "Total ΔV (m/s)",
            ODYSSEY_REFERENCE.total_mission_delta_v_mps,
            "(c) Total ΔV",
        ),
        (
            :thermal_violations,
            "Thermal violations",
            ODYSSEY_REFERENCE.thermal_violations,
            "(d) Thermal violations",
        ),
    ]
    panels = Plots.Plot[]
    for (field, ylabel, reference, title) in specs
        means = Float64[]
        errors = Float64[]
        for policy in policies[1:2]
            μ, σ = _metric_stats(metrics, policy, field)
            push!(means, μ)
            push!(errors, min(σ, μ))
        end
        push!(means, reference)
        push!(errors, 0.0)
        push!(panels, bar(
            labels,
            means;
            yerror=errors,
            ylabel=ylabel,
            legend=false,
            xrotation=15,
            title=title,
        ))
    end
    target_reached_percent = [
        _target_reached_stats(metrics, policy).percent
        for policy in policies[1:2]
    ]
    target_panel = bar(
        labels[1:2],
        target_reached_percent;
        ylabel="Campaigns (%)",
        legend=false,
        xrotation=15,
        ylims=(0, 105),
        yticks=0:20:100,
        title="(e) Target reached (±10 km)",
    )
    for (index, percentage) in pairs(target_reached_percent)
        annotate!(target_panel, index, min(102.0, percentage + 4),
                  text(@sprintf("%.1f%%", percentage), 9))
    end
    push!(panels, target_panel)
    campaigns = count(==("trained_pr_drl"), metrics.policy)
    savefig(plot(
        panels...;
        layout=(2, 3),
        size=(1300, 820),
        plot_title="$(campaigns)-campaign thermal-tolerant flight comparison",
        left_margin=8Plots.mm,
        bottom_margin=5Plots.mm,
    ), path)
    return path
end

function _paper_figure_12(pass_logs::DataFrame, scenario, path::AbstractString)
    policy_rows = pass_logs[pass_logs.policy .== "trained_pr_drl", :]
    episode = minimum(policy_rows.episode)
    trace = policy_rows[policy_rows.episode .== episode, :]
    odyssey_heat = vcat(fill(0.20, 45), fill(0.13, 77), fill(0.04, 88))
    p = plot(trace.pass, trace.heat_rate_w_cm2;
             label="PR-DRL example", linewidth=2, color=:steelblue,
             xlabel="Orbit / atmospheric passage", ylabel="Heat rate (W/cm^2)")
    plot!(p, 1:length(odyssey_heat), odyssey_heat;
          label="Mars Odyssey phase medians (Table III)",
          linewidth=2, linestyle=:dash, color=:black)
    hline!(p, [scenario.reward_config.heat_low_w_cm2];
           label="low corridor", linestyle=:dot, color=:gray)
    hline!(p, [scenario.reward_config.heat_high_w_cm2];
           label="high corridor", linestyle=:dot, color=:red)
    savefig(p, path)
    return path
end

function _checkpoint_record(checkpoint_path, mode, result, payload)
    metrics = DataFrame(result.metrics)
    surpassed_target =
        .!metrics.success .&
        .!metrics.impact .&
        .!metrics.out_of_drag_passage .&
        (metrics.target_error_km .< 0)
    unfinished =
        .!metrics.success .&
        .!metrics.impact .&
        .!metrics.out_of_drag_passage .&
        .!surpassed_target
    step = Int(get(payload, :global_step,
                   something(_numeric_checkpoint_step(checkpoint_path), 0)))
    stat(field) = (mean(Float64.(metrics[!, field])),
                   _std(Float64.(metrics[!, field])))
    reward = stat(:episode_reward)
    target_errors_km = abs.(Float64.(metrics.target_error_km))
    target = (mean(target_errors_km), _std(target_errors_km))
    thermal = stat(:thermal_violations)
    delta_v = stat(:total_delta_v_mps)
    maneuvers = stat(:maneuver_count)
    passes = stat(:pass_count)
    duration = stat(:mission_duration_days)
    return (
        checkpoint=basename(checkpoint_path),
        global_step=step,
        mode=mode,
        episodes=nrow(metrics),
        checkpoint_loss=Float64(get(
            payload,
            :mean_training_loss,
            get(payload, :last_loss, NaN),
        )),
        training_loss_sum=Float64(get(payload, :training_loss_sum, NaN)),
        training_loss_count=Int(get(payload, :training_loss_count, 0)),
        mean_reward=reward[1], std_reward=reward[2],
        success_percent=100 * mean(metrics.success),
        surpassed_target_percent=100 * mean(surpassed_target),
        impact_percent=100 * mean(metrics.impact),
        out_of_drag_passage_percent=100 * mean(metrics.out_of_drag_passage),
        unfinished_percent=100 * mean(unfinished),
        not_reached_goal_percent=100 * mean(surpassed_target .| unfinished),
        mean_target_error_km=target[1], std_target_error_km=target[2],
        mean_thermal_violations=thermal[1], std_thermal_violations=thermal[2],
        mean_delta_v_mps=delta_v[1], std_delta_v_mps=delta_v[2],
        mean_maneuver_count=maneuvers[1], std_maneuver_count=maneuvers[2],
        mean_pass_count=passes[1], std_pass_count=passes[2],
        mean_mission_duration_days=duration[1], std_mission_duration_days=duration[2],
    )
end

function _checkpoint_sweep(run_dir, scenario, episodes, seed, stride,
                           protected_initialization, output_dir)
    paths = filter(path -> _numeric_checkpoint_step(path) !== nothing,
                   frozen_checkpoint_paths(run_dir))
    paths = paths[1:stride:end]
    final_path = joinpath(run_dir, "checkpoint_final.jls")
    isfile(final_path) && push!(paths, final_path)
    isempty(paths) && return DataFrame()

    records = NamedTuple[]
    seen_steps = Set{Int}()
    for path in paths
        policy, _, payload = _load_policy(path)
        step = Int(get(payload, :global_step,
                       something(_numeric_checkpoint_step(path), 0)))
        step in seen_steps && continue
        push!(seen_steps, step)
        @printf("checkpoint evaluation step=%d checkpoint=%s\n", step, path)
        modes = evaluate_policy_modes(
            policy,
            scenario;
            episodes=episodes,
            seed=seed,
            policy_name=splitext(basename(path))[1],
            protected_initialization=protected_initialization,
        )
        for mode in PAPER_EVALUATION_MODES
            push!(records, _checkpoint_record(path, mode, modes[mode], payload))
        end
    end
    data = sort(DataFrame(records), [:global_step, :mode])
    CSV.write(joinpath(output_dir, "checkpoint_metrics.csv"), data)
    return data
end

function _checkpoint_figures(data::DataFrame, output_dir::AbstractString)
    nrow(data) == 0 && return String[]
    paths = String[]
    colors = Dict("conservative" => :seagreen, "tolerant" => :darkorange)

    loss_rows = unique(
        data[:, [
            :global_step,
            :checkpoint_loss,
            :training_loss_sum,
            :training_loss_count,
        ]],
        :global_step,
    )
    sort!(loss_rows, :global_step)
    interval_loss = copy(loss_rows.checkpoint_loss)
    previous_sum = 0.0
    previous_count = 0
    for index in eachindex(interval_loss)
        current_sum = loss_rows.training_loss_sum[index]
        current_count = loss_rows.training_loss_count[index]
        if isfinite(current_sum) && current_count > previous_count
            interval_loss[index] =
                (current_sum - previous_sum) / (current_count - previous_count)
            previous_sum = current_sum
            previous_count = current_count
        end
    end
    p_loss = plot(loss_rows.global_step, interval_loss;
                  xlabel="Training step", ylabel="Average optimizer loss",
                  label="training loss", color=:steelblue, marker=:circle)
    loss_path = joinpath(output_dir, "paper_fig07a_checkpoint_loss.png")
    savefig(p_loss, loss_path)
    push!(paths, loss_path)

    p_reward = plot(xlabel="Training step", ylabel="Evaluation reward")
    for mode in PAPER_EVALUATION_MODES
        rows = data[data.mode .== mode, :]
        plot!(p_reward, rows.global_step, rows.mean_reward;
              ribbon=rows.std_reward, label=mode, color=colors[mode], marker=:circle)
    end
    reward_path = joinpath(output_dir, "paper_fig07b_evaluation_reward.png")
    savefig(p_reward, reward_path)
    push!(paths, reward_path)

    rows = data[data.mode .== "tolerant", :]
    completion_values = hcat(
        rows.success_percent,
        rows.surpassed_target_percent,
        rows.impact_percent,
        rows.out_of_drag_passage_percent,
        rows.unfinished_percent,
    )
    completion = plot(rows.global_step, completion_values;
                      xlabel="Training step", ylabel="Episodes (%)",
                      label=permutedims([
                          "reached goal",
                          "surpassed target",
                          "impact",
                          "out of passage",
                          "unfinished",
                      ]),
                      title="(a) Episode completion", marker=:circle,
                      linewidth=2, ylims=(0, 100))
    target = plot(rows.global_step, rows.mean_target_error_km;
                  ribbon=(
                      min.(rows.std_target_error_km, rows.mean_target_error_km),
                      rows.std_target_error_km,
                  ), label="mean ± standard deviation",
                  xlabel="Training step", ylabel="Absolute target distance (km)",
                  title="(b) Final distance to target apoapsis radius",
                  color=:steelblue, linewidth=2, marker=:circle)
    thermal = plot(rows.global_step, rows.mean_thermal_violations;
                   ribbon=rows.std_thermal_violations, label=false,
                   xlabel="Training step", ylabel="Violations",
                   title="(c) Thermal violations", marker=:circle)
    performance = plot(rows.global_step, rows.mean_delta_v_mps;
                       ribbon=rows.std_delta_v_mps, label="Delta-V",
                       xlabel="Training step", ylabel="Value",
                       title="(d) Delta-V and ABMs", marker=:circle)
    plot!(performance, rows.global_step, rows.mean_maneuver_count;
          ribbon=rows.std_maneuver_count, label="ABMs", marker=:circle)
    duration = plot(rows.global_step, rows.mean_pass_count;
                    ribbon=rows.std_pass_count, label="passages",
                    xlabel="Training step", ylabel="Value",
                    title="(e) Passages and duration", marker=:circle)
    plot!(duration, rows.global_step, rows.mean_mission_duration_days;
          ribbon=rows.std_mission_duration_days, label="days", marker=:circle)
    figure8_path = joinpath(output_dir, "paper_fig08_checkpoint_performance.png")
    savefig(plot(completion, target, thermal, performance, duration;
                 layout=(3, 2), size=(1100, 1050)), figure8_path)
    push!(paths, figure8_path)

    completion_target_path = joinpath(
        output_dir,
        "episode_completion_and_final_target_distance.png",
    )
    savefig(
        plot(
            completion,
            target;
            layout=(1, 2),
            size=(1200, 480),
            plot_title="Greedy thermal-tolerant checkpoint evaluation",
            left_margin=8Plots.mm,
            bottom_margin=5Plots.mm,
        ),
        completion_target_path,
    )
    push!(paths, completion_target_path)
    return paths
end

function _repeated_action_stats(pass_logs::DataFrame)
    rows = NamedTuple[]
    policy_rows = pass_logs[pass_logs.policy .== "trained_pr_drl", :]
    by_action = Dict{Float64,Vector{Int}}()
    for episode_rows in groupby(policy_rows, :episode)
        actions = Float64.(episode_rows.action_delta_v_mps)
        isempty(actions) && continue
        run_action = first(actions)
        run_length = 1
        for action in Iterators.drop(actions, 1)
            if action == run_action
                run_length += 1
            else
                push!(get!(by_action, run_action, Int[]), run_length)
                run_action = action
                run_length = 1
            end
        end
        push!(get!(by_action, run_action, Int[]), run_length)
    end
    for action in sort(collect(keys(by_action)))
        lengths = by_action[action]
        push!(rows, (
            action_delta_v_mps=action,
            mean_consecutive_repetitions=mean(lengths),
            std_consecutive_repetitions=_std(lengths),
            action_runs=length(lengths),
        ))
    end
    return DataFrame(rows)
end

function _paper_figure_9(pass_logs::DataFrame, output_dir::AbstractString)
    stats = _repeated_action_stats(pass_logs)
    csv_path = joinpath(output_dir, "repeated_action_metrics.csv")
    CSV.write(csv_path, stats)
    plot_path = joinpath(output_dir, "paper_fig09_repeated_actions.png")
    p = bar(string.(stats.action_delta_v_mps), stats.mean_consecutive_repetitions;
            yerror=stats.std_consecutive_repetitions,
            xlabel="Action Delta-V (m/s)",
            ylabel="Mean consecutive selections",
            legend=false)
    savefig(p, plot_path)
    return (metrics=csv_path, figure=plot_path)
end

function _randomization_variant(config; process_noise_scale=config.process_noise_scale,
                                marsgram_perturbation_scale=config.marsgram_perturbation_scale)
    return AerobrakingRandomizationConfig(
        nominal=config.nominal,
        apoapsis_jitter_m=config.apoapsis_jitter_m,
        periapsis_jitter_m=config.periapsis_jitter_m,
        angle_jitter_deg=config.angle_jitter_deg,
        nonnominal_inclination_low_deg=config.nonnominal_inclination_low_deg,
        nonnominal_inclination_high_deg=config.nonnominal_inclination_high_deg,
        nonnominal_aop_low_deg=config.nonnominal_aop_low_deg,
        nonnominal_aop_high_deg=config.nonnominal_aop_high_deg,
        nonnominal_raan_low_deg=config.nonnominal_raan_low_deg,
        nonnominal_raan_high_deg=config.nonnominal_raan_high_deg,
        initial_date_start=config.initial_date_start,
        initial_date_days=config.initial_date_days,
        randomize_initial_time_of_day=config.randomize_initial_time_of_day,
        initial_true_anomaly_jitter_deg=config.initial_true_anomaly_jitter_deg,
        process_noise=process_noise_scale != 0,
        process_noise_scale=process_noise_scale,
        aerodynamic_coefficient_dispersion=config.aerodynamic_coefficient_dispersion,
        aerodynamic_coefficient_span=config.aerodynamic_coefficient_span,
        aerodynamic_cd_span=config.aerodynamic_cd_span,
        aerodynamic_cl_span=config.aerodynamic_cl_span,
        marsgram_perturbation_scale=marsgram_perturbation_scale,
        marsgram_seed_base=config.marsgram_seed_base,
    )
end

function _scenario_variant(base; phase=base.phase, backend_mode=base.backend_mode,
                           randomization=base.randomization_config)
    return default_aerobraking_config(
        phase=phase,
        nominal=randomization.nominal,
        max_passes=base.termination_config.max_passes,
        backend_mode=backend_mode,
        training=false,
        spaceagora_atmosphere_model=base.spaceagora_atmosphere_model,
        spaceagora_gram_once_per_step=base.spaceagora_gram_once_per_step,
        spaceagora_tabulated_flight_file=base.spaceagora_tabulated_flight_file,
        spaceagora_tabulated_flight_sigma=base.spaceagora_tabulated_flight_sigma,
        spaceagora_gravity_harmonics_degree=base.spaceagora_gravity_harmonics_degree,
        spaceagora_gravity_harmonics_order=base.spaceagora_gravity_harmonics_order,
        spaceagora_gravity_harmonics_file=base.spaceagora_gravity_harmonics_file,
        reward_config=base.reward_config,
        termination_config=base.termination_config,
        randomization_config=randomization,
    )
end

function _generalization_scenarios(base)
    aggressive = _randomization_variant(
        base.randomization_config;
        process_noise_scale=max(0.8, 2 * base.randomization_config.process_noise_scale),
        marsgram_perturbation_scale=2.0,
    )
    return Dict(
        "nominal" => _scenario_variant(base),
        "exponential_density" => _scenario_variant(base; backend_mode=:paper_surrogate),
        "aggressive_atmosphere" => _scenario_variant(base; randomization=aggressive),
        "short_campaign" => _scenario_variant(base; phase="Walkout"),
        "long_campaign" => _scenario_variant(base; phase="Campaign"),
        "accurate_simulator" => _scenario_variant(base),
    )
end

function _generalization_record(case_name, result, nominal_reward)
    summaries = result.summaries
    field_values(field) = Float64[getfield(summary, field) for summary in summaries]
    rewards = field_values(:episode_reward)
    thermal = field_values(:thermal_violations)
    distance = field_values(:target_error_m) ./ 1000
    delta_v = field_values(:total_delta_v_mps)
    duration = field_values(:mission_duration_days)
    length_values = field_values(:pass_count)
    maneuvers = field_values(:maneuver_count)
    return (
        case=case_name,
        episodes=length(summaries),
        generalization_gap=nominal_reward - mean(rewards),
        mean_reward=mean(rewards), std_reward=_std(rewards),
        mean_thermal_violations=mean(thermal), std_thermal_violations=_std(thermal),
        reached_goal_percent=100 * mean(getfield.(summaries, :success)),
        mean_goal_distance_km=mean(distance), std_goal_distance_km=_std(distance),
        mean_delta_v_mps=mean(delta_v), std_delta_v_mps=_std(delta_v),
        mean_mission_duration_days=mean(duration), std_mission_duration_days=_std(duration),
        mean_episode_length=mean(length_values), std_episode_length=_std(length_values),
        mean_maneuver_count=mean(maneuvers), std_maneuver_count=_std(maneuvers),
    )
end

function _generalization_table(policy, scenario, episodes, seed,
                               protected_initialization, output_dir)
    scenarios = _generalization_scenarios(scenario)
    results = Dict{String,Any}()
    for case_name in PAPER_GENERALIZATION_CASES
        @printf("generalization evaluation case=%s episodes=%d\n", case_name, episodes)
        results[case_name] = evaluate_policy(
            policy,
            scenarios[case_name];
            episodes=episodes,
            seed=seed,
            policy_name="trained_pr_drl_$case_name",
            paper_protocol=false,
            protected_initialization=protected_initialization,
        )
    end
    nominal_reward = mean(getfield.(results["nominal"].summaries, :episode_reward))
    table = DataFrame([
        _generalization_record(case_name, results[case_name], nominal_reward)
        for case_name in PAPER_GENERALIZATION_CASES
    ])
    path = joinpath(output_dir, "paper_table_vi_generalization.csv")
    CSV.write(path, table)
    return path
end

function _write_manifest(output_dir; run_dir, config_path, config_sha256,
                         checkpoint_path, algorithm,
                         options, artifacts)
    artifact_dict = Dict(String(key) => String(value) for (key, value) in pairs(artifacts))
    manifest = Dict{String,Any}(
        "run_dir" => abspath(run_dir),
        "config_path" => abspath(config_path),
        "config_sha256" => config_sha256,
        "checkpoint_path" => abspath(checkpoint_path),
        "algorithm" => String(algorithm),
        "seed" => options.seed,
        "iid_episodes" => options.episodes,
        "flight_episodes" => options.flight_episodes,
        "generalization_episodes" => options.generalization_episodes,
        "checkpoint_episodes" => options.checkpoint_episodes,
        "paper_reference" => Dict(
            "title" => "Autonomous Decision-Making for Aerobraking via Parallel Randomized Deep Reinforcement Learning",
            "doi" => "10.1109/TAES.2022.3221697",
            "odyssey_reference_source" => "Paper Table III, Main Phase II + Endgame + Walkout",
            "figure_12_reference" => "Piecewise phase median heat rates from Paper Table III",
        ),
        "generalization_case_mapping" => Dict(
            "nominal" => "run evaluation scenario",
            "exponential_density" => "paper_surrogate exponential atmosphere",
            "aggressive_atmosphere" => "2x process-noise scale (minimum 0.8) and MarsGRAM perturbation scale 2",
            "short_campaign" => "Walkout phase initial conditions",
            "long_campaign" => "full Campaign initial conditions",
            "accurate_simulator" => "run backend with SpaceAGORA_RL evaluation solver settings",
        ),
        "artifacts" => artifact_dict,
    )
    path = joinpath(output_dir, "evaluation_manifest.toml")
    open(path, "w") do io
        TOML.print(io, manifest)
    end
    return path
end

function _final_flight_scenario(resolved)
    max_passes = max(1000, resolved.scenario.termination_config.max_passes)
    odyssey_defaults = paper_odyssey_flight_evaluation_config(
        backend_mode=resolved.scenario.backend_mode,
        max_passes=max_passes,
    )
    return paper_evaluation_scenario(
        resolved.scenario;
        max_passes=max_passes,
        randomization_config=odyssey_defaults.randomization_config,
        terminal_on_thermal_violation=false,
    )
end

function evaluate_final_flight_comparison(options)
    run_dir = abspath(options.run_dir)
    isdir(run_dir) || throw(ArgumentError("run directory does not exist: $run_dir"))
    options.checkpoint in (nothing, "final") ||
        throw(ArgumentError("--final-flight-comparison always uses checkpoint_final.jls"))

    manifest_path = joinpath(run_dir, "manifest.toml")
    isfile(manifest_path) || throw(ArgumentError("run has no manifest.toml: $run_dir"))
    run_manifest = TOML.parsefile(manifest_path)
    config_path = _config_path(run_dir, run_manifest, options.config)
    config_digest = bytes2hex(sha256(read(config_path)))
    recorded_digest = strip(String(get(run_manifest, "config_sha256", "")))
    if !isempty(recorded_digest) && recorded_digest != config_digest
        @warn "evaluation config differs from the config recorded by the training run" config_path recorded_digest config_digest
    end

    resolved = resolve_config(config_path)
    checkpoint_path = _final_checkpoint_path(run_dir)
    policy, algorithm, _ = _load_policy(checkpoint_path)
    output_dir = abspath(options.output === nothing ?
                        joinpath(run_dir, "final_flight_comparison") :
                        options.output)
    mkpath(output_dir)
    scenario = _final_flight_scenario(resolved)
    scenario.termination_config.terminal_on_thermal_violation &&
        error("final flight comparison must use thermal-tolerant termination")
    seed = resolved.training.validation_seed
    protected = protected_initialization_config(resolved.training)

    @printf(
        "final flight comparison run=%s algorithm=%s checkpoint=%s campaigns_per_policy=%d thermal_terminal=false processes=%d threads_per_process=%d progress_every=%d config=%s output=%s\n",
        run_dir,
        algorithm,
        checkpoint_path,
        options.episodes,
        options.processes,
        options.threads_per_process,
        options.progress_every,
        config_path,
        output_dir,
    )
    flush(stdout)
    process_ids = _start_evaluation_workers(
        options.processes,
        options.threads_per_process,
        scenario,
    )
    results = try
        _evaluate_pair_parallel(
            policy,
            scenario,
            options.episodes,
            seed,
            protected,
            process_ids;
            paper_protocol=false,
            progress_every=options.progress_every,
        )
    finally
        println("stopping evaluation workers")
        flush(stdout)
        rmprocs(process_ids)
    end
    raw_paths = write_evaluation_artifacts(
        joinpath(output_dir, "final_policy_vs_aads"),
        results,
    )
    metrics, _ = _result_frames(results)
    comparison_path = joinpath(
        output_dir,
        "final_flight_performance_comparison.csv",
    )
    comparison = _flight_performance_table(metrics)
    CSV.write(comparison_path, comparison)
    for policy in ("PR-DRL", "AADS")
        row = only(eachrow(comparison[comparison.policy .== policy, :]))
        @printf(
            "final checkpoint target reached policy=%s tolerance_km=%.1f campaigns=%d/%d percent=%.1f%%\n",
            policy,
            PAPER_TARGET_TOLERANCE_KM,
            row.target_reached_10km_count,
            row.episodes,
            row.target_reached_10km_percent,
        )
    end
    flush(stdout)
    figure_path = _paper_figure_11(
        metrics,
        joinpath(output_dir, "final_flight_performance_comparison.png"),
    )

    evaluation_manifest = joinpath(output_dir, "evaluation_manifest.toml")
    manifest = Dict{String,Any}(
        "evaluation_mode" => "final_flight_comparison",
        "run_dir" => run_dir,
        "checkpoint_path" => abspath(checkpoint_path),
        "config_path" => abspath(config_path),
        "config_sha256" => config_digest,
        "algorithm" => String(algorithm),
        "campaigns_per_policy" => options.episodes,
        "worker_processes" => options.processes,
        "threads_per_process" => options.threads_per_process,
        "progress_every_campaigns" => options.progress_every,
        "seed" => seed,
        "terminal_on_thermal_violation" => false,
        "target_reached_tolerance_km" => PAPER_TARGET_TOLERANCE_KM,
        "policy_action_selection" => "greedy",
        "odyssey_reference" => Dict(
            "maneuver_count" => ODYSSEY_REFERENCE.maneuver_count,
            "mission_duration_days" => ODYSSEY_REFERENCE.mission_duration_days,
            "total_mission_delta_v_mps" =>
                ODYSSEY_REFERENCE.total_mission_delta_v_mps,
            "thermal_violations" => ODYSSEY_REFERENCE.thermal_violations,
        ),
        "artifacts" => Dict(
            "episode_metrics" => raw_paths.metrics,
            "summary_metrics" => raw_paths.aggregate,
            "pass_logs" => raw_paths.pass_logs,
            "comparison_csv" => comparison_path,
            "comparison_figure" => figure_path,
        ),
    )
    open(evaluation_manifest, "w") do io
        TOML.print(io, manifest)
    end
    println("final flight comparison complete: ", output_dir)
    println("comparison figure: ", figure_path)
    println("comparison metrics: ", comparison_path)
    return (
        output_dir=output_dir,
        checkpoint=checkpoint_path,
        metrics=comparison_path,
        figure=figure_path,
        manifest=evaluation_manifest,
    )
end

function evaluate_run(options)
    options.final_flight_comparison &&
        return evaluate_final_flight_comparison(options)
    run_dir = abspath(options.run_dir)
    isdir(run_dir) || throw(ArgumentError("run directory does not exist: $run_dir"))
    manifest_path = joinpath(run_dir, "manifest.toml")
    isfile(manifest_path) || throw(ArgumentError("run has no manifest.toml: $run_dir"))
    run_manifest = TOML.parsefile(manifest_path)
    config_path = _config_path(run_dir, run_manifest, options.config)
    config_digest = bytes2hex(sha256(read(config_path)))
    recorded_digest = strip(String(get(run_manifest, "config_sha256", "")))
    if !isempty(recorded_digest) && recorded_digest != config_digest
        @warn "evaluation config differs from the config recorded by the training run" config_path recorded_digest config_digest
    end
    resolved = resolve_config(config_path)
    checkpoint_path = _checkpoint_path(run_dir, options.checkpoint)
    policy, algorithm, _ = _load_policy(checkpoint_path)
    output_dir = abspath(options.output === nothing ?
                        joinpath(run_dir, "paper_evaluation") : options.output)
    mkpath(output_dir)
    seed = resolved.training.validation_seed
    protected = protected_initialization_config(resolved.training)
    scenario = paper_evaluation_scenario(
        resolved.scenario;
        max_passes=max(1000, resolved.scenario.termination_config.max_passes),
    )
    runtime_options = merge(options, (; seed=seed))

    @printf("evaluating run=%s algorithm=%s checkpoint=%s config=%s output=%s\n",
            run_dir, algorithm, checkpoint_path, config_path, output_dir)

    iid_dir = joinpath(output_dir, "iid_pr_drl_vs_aads")
    iid = _evaluate_pair(policy, scenario, options.episodes, seed, protected)
    iid_paths = write_evaluation_artifacts(iid_dir, iid)
    iid_metrics, iid_pass_logs = _result_frames(iid)
    table_v_path = joinpath(output_dir, "paper_table_v_pr_drl_vs_aads.csv")
    CSV.write(table_v_path, _table_v(iid_metrics))
    fig10_path = _paper_figure_10(
        iid_pass_logs,
        scenario,
        joinpath(output_dir, "paper_fig10_heat_rate_and_delta_v.png"),
    )
    fig9 = _paper_figure_9(iid_pass_logs, output_dir)

    flight_scenario = _final_flight_scenario(resolved)
    flight_dir = joinpath(output_dir, "odyssey_flight_comparison")
    flight = _evaluate_pair(
        policy,
        flight_scenario,
        options.flight_episodes,
        seed,
        protected;
        paper_protocol=false,
    )
    flight_paths = write_evaluation_artifacts(flight_dir, flight)
    flight_metrics, flight_pass_logs = _result_frames(flight)
    fig11_path = _paper_figure_11(
        flight_metrics,
        joinpath(output_dir, "paper_fig11_flight_performance_summary.png"),
    )
    fig12_path = _paper_figure_12(
        flight_pass_logs,
        flight_scenario,
        joinpath(output_dir, "paper_fig12_heat_rate_odyssey_comparison.png"),
    )

    checkpoint_csv = ""
    checkpoint_figure_paths = String[]
    if options.checkpoint_sweep
        checkpoint_data = _checkpoint_sweep(
            run_dir,
            scenario,
            options.checkpoint_episodes,
            seed,
            options.checkpoint_stride,
            protected,
            output_dir,
        )
        checkpoint_figure_paths = _checkpoint_figures(checkpoint_data, output_dir)
        checkpoint_csv = nrow(checkpoint_data) == 0 ? "" :
                         joinpath(output_dir, "checkpoint_metrics.csv")
    end

    table_vi_path = ""
    if options.generalization
        table_vi_path = _generalization_table(
            policy,
            scenario,
            options.generalization_episodes,
            seed,
            protected,
            output_dir,
        )
    end

    artifacts = (
        iid_episode_metrics=iid_paths.metrics,
        iid_summary_metrics=iid_paths.aggregate,
        iid_pass_logs=iid_paths.pass_logs,
        flight_episode_metrics=flight_paths.metrics,
        flight_summary_metrics=flight_paths.aggregate,
        flight_pass_logs=flight_paths.pass_logs,
        checkpoint_metrics=checkpoint_csv,
        repeated_action_metrics=fig9.metrics,
        paper_table_v=table_v_path,
        paper_table_vi=table_vi_path,
        paper_fig07a=length(checkpoint_figure_paths) >= 1 ? checkpoint_figure_paths[1] : "",
        paper_fig07b=length(checkpoint_figure_paths) >= 2 ? checkpoint_figure_paths[2] : "",
        paper_fig08=length(checkpoint_figure_paths) >= 3 ? checkpoint_figure_paths[3] : "",
        episode_completion_and_final_target_distance=
            length(checkpoint_figure_paths) >= 4 ? checkpoint_figure_paths[4] : "",
        paper_fig09=fig9.figure,
        paper_fig10=fig10_path,
        paper_fig11=fig11_path,
        paper_fig12=fig12_path,
    )
    evaluation_manifest = _write_manifest(
        output_dir;
        run_dir=run_dir,
        config_path=config_path,
        config_sha256=config_digest,
        checkpoint_path=checkpoint_path,
        algorithm=algorithm,
        options=runtime_options,
        artifacts=artifacts,
    )
    println("evaluation complete: ", output_dir)
    println("evaluation manifest: ", evaluation_manifest)
    return (output_dir=output_dir, artifacts=artifacts, manifest=evaluation_manifest)
end

function main(args=ARGS)
    options = _parse_cli(args)
    if options.help
        _usage()
        return nothing
    end
    return evaluate_run(options)
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main()
    catch error
        showerror(stderr, error)
        println(stderr)
        println(stderr, "Use --help for usage.")
        exit(1)
    end
end
