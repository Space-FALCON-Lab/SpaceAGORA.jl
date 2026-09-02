#!/usr/bin/env julia

"""
Evaluate a completed SpaceAGORA_RL run using the protocol and artifacts reported
in Falcone and Putnam (2023), "Autonomous Decision-Making for Aerobraking via
Parallel Randomized Deep Reinforcement Learning."

The script accepts one or more run directories, discovers their
manifest/config/checkpoints, and writes raw episode/pass data, paper tables,
and paper-style figures.
Run with `--help` for the command-line interface.
"""

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using CSV
using DataFrames
using Dates
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
  --checkpoint PATH|best|final
                              Policy checkpoint (default: best validation
                              checkpoint, then final/latest as a fallback).
  --compare-run RUN_DIR       Add a trained run to a combined flight comparison.
                              May be repeated. Additional runs use their own
                              recorded config and best/final checkpoint. All
                              policies are evaluated on the primary RUN_DIR's
                              scenario and campaign seeds. This mode adds AADS
                              and the Mars Odyssey reference automatically.
  --wind-mode MODE            MarsGRAM wind mode: zero, nominal, or perturbed.
                              Overrides spaceagora_physics.gram_wind_mode.
  --final-flight-comparison   Only evaluate the selected checkpoint against AADS and
                              the Mars Odyssey mission reference. Uses the
                              thermal-tolerant protocol and --episodes campaigns
                              (default: 40), including the percentage that finish
                              within ±10 km of the target apoapsis radius.
  --episodes N                IID PR-DRL/AADS episodes for Table V/Fig. 10,
                              or campaigns per policy in flight comparisons
                              (default: 40).
  --flight-episodes N         Odyssey-geometry episodes for Figs. 11-12
                              (default: 40).
  --generalization-episodes N Episodes per IID reference and Table VI-inspired
                              SpaceAGORA case (default: 100).
  --checkpoint-episodes N     Episodes per mode and checkpoint for Figs. 7-8
                              (default: 40).
  --checkpoint-stride N       Evaluate every Nth numbered checkpoint (default: 1).
  --processes N               Worker processes for parallel campaign evaluation,
                              including generalization (default: 16).
  --threads-per-process N     Julia threads available inside each worker process
                              (default: 1).
  --progress-every N          Print progress every N completed campaigns
                              (default: 1).
  --generalization-only       Run only the frozen-policy SpaceAGORA
                              generalization evaluation suite.
  --skip-checkpoint-sweep     Skip Figs. 7-8 checkpoint evaluation.
  --skip-generalization       Skip Table VI generalization evaluation.
  --help                      Show this message.

Outputs include the raw episode/pass CSVs, checkpoint_metrics.csv,
paper_table_v_pr_drl_vs_aads.csv, the generalization_evaluation_suite directory,
generalization_results_table.pdf,
paper_fig07_*.png through paper_fig12_*.png, the dedicated episode-completion
and final-target-distance chart, and evaluation_manifest.toml.
Final-flight-comparison mode writes its four-panel PNG, comparison CSV, raw
episode/pass CSVs, and evaluation_manifest.toml only.
With --compare-run, the equivalent combined artifacts are written outside the
training runs under outputs/comparisons by default; use --output to choose an
explicit comparison directory.
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
        :comparison_runs => String[],
        :wind_mode => nothing,
        :final_flight_comparison => false,
        :episodes => PAPER_IID_EVALUATION_EPISODES,
        :flight_episodes => PAPER_IID_EVALUATION_EPISODES,
        :generalization_episodes => PAPER_GENERALIZATION_EVALUATION_EPISODES,
        :checkpoint_episodes => PAPER_IID_EVALUATION_EPISODES,
        :checkpoint_stride => 1,
        :processes => 16,
        :threads_per_process => 1,
        :progress_every => 1,
        :generalization_only => false,
        :checkpoint_sweep => true,
        :generalization => true,
    )
    value_options = Dict(
        "--output" => :output,
        "--config" => :config,
        "--checkpoint" => :checkpoint,
        "--wind-mode" => :wind_mode,
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
        elseif arg == "--generalization-only"
            options[:generalization_only] = true
            index += 1
        elseif arg == "--compare-run"
            index == length(args) &&
                throw(ArgumentError("missing value for $arg"))
            push!(options[:comparison_runs], args[index + 1])
            index += 2
        elseif haskey(value_options, arg)
            index == length(args) &&
                throw(ArgumentError("missing value for $arg"))
            key = value_options[arg]
            value = args[index + 1]
            options[key] = if key in integer_options
                parse(Int, value)
            elseif key == :wind_mode
                canonical_gram_wind_mode(value)
            else
                value
            end
            index += 2
        else
            throw(ArgumentError("unknown option: $arg"))
        end
    end
    for key in integer_options
        options[key] > 0 || throw(ArgumentError("$key must be positive"))
    end
    options[:generalization_only] && !options[:generalization] &&
        throw(ArgumentError("--generalization-only cannot be combined with --skip-generalization"))
    options[:generalization_only] && options[:final_flight_comparison] &&
        throw(ArgumentError("--generalization-only cannot be combined with --final-flight-comparison"))
    options[:generalization_only] && !isempty(options[:comparison_runs]) &&
        throw(ArgumentError("--generalization-only cannot be combined with --compare-run"))
    return (; help=false, options...)
end

function _resolve_existing_path(path::AbstractString, run_dir::AbstractString)
    candidates = unique([
        abspath(path),
        abspath(joinpath(run_dir, path)),
        abspath(joinpath(dirname(run_dir), path)),
        abspath(joinpath(dirname(dirname(run_dir)), path)),
        abspath(joinpath(dirname(dirname(dirname(run_dir))), path)),
        abspath(joinpath(dirname(dirname(dirname(dirname(run_dir)))), path)),
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

function _best_checkpoint_path(run_dir::AbstractString)
    record_path = joinpath(
        run_dir,
        "checkpoint_validation",
        "best_validation_checkpoint.txt",
    )
    isfile(record_path) || return nothing

    checkpoint_name = nothing
    for line in eachline(record_path)
        key_value = split(line, "="; limit=2)
        length(key_value) == 2 || continue
        if strip(first(key_value)) == "checkpoint"
            checkpoint_name = strip(last(key_value))
            break
        end
    end
    checkpoint_name === nothing && throw(ArgumentError(
        "best-checkpoint record has no checkpoint entry: $record_path",
    ))
    isempty(checkpoint_name) && throw(ArgumentError(
        "best-checkpoint record has an empty checkpoint entry: $record_path",
    ))
    basename(checkpoint_name) == checkpoint_name || throw(ArgumentError(
        "best-checkpoint record must name a run-local checkpoint: $record_path",
    ))

    checkpoint_path = joinpath(run_dir, checkpoint_name)
    isfile(checkpoint_path) || throw(ArgumentError(
        "best validation checkpoint does not exist: $checkpoint_path",
    ))
    return checkpoint_path
end

function _checkpoint_path(run_dir::AbstractString, requested)
    if requested !== nothing && requested != "best" && requested != "final"
        path = _resolve_existing_path(String(requested), run_dir)
        path === nothing && throw(ArgumentError("checkpoint does not exist: $requested"))
        return path
    end
    if requested === nothing || requested == "best"
        best_path = _best_checkpoint_path(run_dir)
        best_path !== nothing && return best_path
        requested == "best" && throw(ArgumentError(
            "run has no checkpoint_validation/best_validation_checkpoint.txt: $run_dir",
        ))
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

function _load_policy(checkpoint_path::AbstractString)
    payload = load_checkpoint(checkpoint_path)
    algorithm = Symbol(get(payload, :algorithm, haskey(payload, :actor) ? :a2c : :pr_drl))
    policy = if algorithm == :a2c
        load_trained_a2c_policy(checkpoint_path)
    elseif algorithm == :a3c
        load_trained_a3c_policy(checkpoint_path)
    elseif algorithm == :td3
        load_trained_td3_policy(checkpoint_path)
    else
        load_trained_pr_drl_policy(checkpoint_path)
    end
    return policy, algorithm, payload
end

function _load_evaluation_source(run_directory::AbstractString;
                                 config_override=nothing,
                                 checkpoint_requested=nothing,
                                 wind_mode=nothing)
    run_dir = abspath(run_directory)
    isdir(run_dir) || throw(ArgumentError("run directory does not exist: $run_dir"))
    manifest_path = joinpath(run_dir, "manifest.toml")
    isfile(manifest_path) || throw(ArgumentError("run has no manifest.toml: $run_dir"))
    run_manifest = TOML.parsefile(manifest_path)
    config_path = _config_path(run_dir, run_manifest, config_override)
    config_digest = bytes2hex(sha256(read(config_path)))
    recorded_digest = strip(String(get(run_manifest, "config_sha256", "")))
    if !isempty(recorded_digest) && recorded_digest != config_digest
        @warn "evaluation config differs from the config recorded by the training run" run_dir config_path recorded_digest config_digest
    end
    resolved = resolve_config(config_path; gram_wind_mode=wind_mode)
    checkpoint_path = _checkpoint_path(run_dir, checkpoint_requested)
    policy, algorithm, payload = _load_policy(checkpoint_path)
    run_id = String(get(run_manifest, "run_id", basename(run_dir)))
    return (
        run_dir=run_dir,
        run_id=run_id,
        manifest_path=manifest_path,
        config_path=config_path,
        config_sha256=config_digest,
        checkpoint_path=checkpoint_path,
        policy=policy,
        algorithm=algorithm,
        checkpoint_payload=payload,
        resolved=resolved,
    )
end

function _comparison_policy_specs(sources)
    algorithm_counts = Dict{Symbol,Int}()
    for source in sources
        algorithm_counts[source.algorithm] = get(algorithm_counts, source.algorithm, 0) + 1
    end
    occurrences = Dict{Symbol,Int}()
    return map(sources) do source
        occurrence = get(occurrences, source.algorithm, 0) + 1
        occurrences[source.algorithm] = occurrence
        base_key = "trained_$(source.algorithm)"
        duplicated = algorithm_counts[source.algorithm] > 1
        policy_key = duplicated ? "$(base_key)_$(occurrence)" : base_key
        algorithm_label = SpaceAGORA_RL.algorithm_display_name(source.algorithm)
        policy_label = duplicated ? "$(algorithm_label) ($(source.run_id))" : algorithm_label
        return (key=policy_key, label=policy_label, source=source)
    end
end

_field_values(value) = Tuple(getfield(value, index) for index in 1:fieldcount(typeof(value)))

function _comparison_environment_signature(resolved)
    scenario = resolved.scenario
    training = resolved.training
    return (
        phase=scenario.phase,
        backend_mode=scenario.backend_mode,
        atmosphere_model=scenario.spaceagora_atmosphere_model,
        gram_wind_mode=scenario.spaceagora_gram_wind_mode,
        gram_once_per_step=scenario.spaceagora_gram_once_per_step,
        mars_mgcm_dust_levels=scenario.spaceagora_mars_mgcm_dust_levels,
        mars_dust_storm=scenario.spaceagora_mars_dust_storm,
        integration=_field_values(scenario.spaceagora_integration_config),
        gravity_degree=scenario.spaceagora_gravity_harmonics_degree,
        gravity_order=scenario.spaceagora_gravity_harmonics_order,
        gravity_file=scenario.spaceagora_gravity_harmonics_file,
        reward=_field_values(scenario.reward_config),
        termination=_field_values(scenario.termination_config),
        randomization=_field_values(scenario.randomization_config),
        normalization=_field_values(scenario.normalization_bounds),
        protected_initialization=(
            training.protected_first_pass,
            training.protected_initial_corridor_maneuver,
            training.protected_first_pass_suppress_thermal_terminal,
            training.protected_corridor_low_w_cm2,
            training.protected_corridor_high_w_cm2,
        ),
    )
end

function _default_comparison_output_dir(sources)
    primary_parent = dirname(first(sources).run_dir)
    output_root = basename(primary_parent) == "runs" ? dirname(primary_parent) : primary_parent
    comparison_name = join(
        replace.(basename.(getfield.(sources, :run_dir)), r"[^A-Za-z0-9._-]" => "_"),
        "_vs_",
    )
    return joinpath(output_root, "comparisons", comparison_name)
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

function _result_from_parallel_episodes(episodes, scenario, policy_name::AbstractString)
    sort!(episodes; by=episode -> episode.summary.episode_index)
    summaries = EpisodeSummary[episode.summary for episode in episodes]
    result = _result_from_summaries(summaries, scenario, policy_name)
    transitions = Transition[]
    for episode in episodes
        append!(transitions, episode.transitions)
    end
    return merge(result, (; transitions))
end

function _evaluate_policies_parallel(
    policies::Dict{String,Any},
    policy_order::Vector{String},
    scenario,
    episodes::Int,
    seed::Int,
    protected_initialization,
    process_ids::Vector{Int};
    paper_protocol::Bool=true,
    progress_every::Int=1,
    collect_transitions::Bool=false,
    campaign_callback=nothing,
)
    isempty(process_ids) &&
        throw(ArgumentError("parallel evaluation requires at least one worker process"))
    isempty(policy_order) &&
        throw(ArgumentError("parallel evaluation requires at least one policy"))
    all(haskey(policies, policy_name) for policy_name in policy_order) ||
        throw(ArgumentError("policy_order contains a policy missing from policies"))
    jobs = [
        (
            policy_name=policy_name,
            episode=episode,
            campaign_seed=seed + episode - 1,
        )
        for policy_name in policy_order
        for episode in 1:episodes
    ]
    total = length(jobs)
    progress = RemoteChannel(() -> Channel{Any}(total), myid())
    pool = CachingPool(process_ids)
    start_time = time()
    @printf(
        "campaign evaluation starting total=%d policies=%d campaigns_per_policy=%d workers=%d\n",
        total,
        length(policy_order),
        episodes,
        length(process_ids),
    )
    flush(stdout)

    evaluation_task = @async pmap(pool, jobs; batch_size=1) do job
        evaluation_policy = policies[job.policy_name]
        summary, episode_transitions = if SpaceAGORA_RL._is_spaceagora_live_backend(scenario.backend_mode)
            campaign_summary, campaign_transitions =
                SpaceAGORA_RL.run_spaceagora_physics_policy_campaign_episode(
                    evaluation_policy,
                    scenario,
                    job.episode,
                    Distributed.myid(),
                    job.campaign_seed;
                    max_passes_per_campaign=scenario.termination_config.max_passes,
                    protected_initialization=protected_initialization,
                )
            campaign_summary, campaign_transitions
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
            only(episode_result.summaries), episode_result.transitions
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
        return (
            policy_name=job.policy_name,
            summary=summary,
            transitions=collect_transitions ? episode_transitions : Transition[],
        )
    end

    completed = 0
    last_status_time = time()
    while completed < total
        if isready(progress)
            message = take!(progress)
            completed += 1
            last_status_time = time()
            campaign_callback === nothing || campaign_callback(completed, message)
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
    episode_results = Dict(policy_name => Any[] for policy_name in policy_order)
    for item in evaluated
        push!(episode_results[item.policy_name], item)
    end
    return Dict(
        policy_name =>
            _result_from_parallel_episodes(
                episode_results[policy_name],
                scenario,
                policy_name,
            )
        for policy_name in policy_order
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
    policies = Dict{String,Any}(
        "aads_heuristic" => AADSHeuristicPolicy(),
        "trained_pr_drl" => policy,
    )
    return _evaluate_policies_parallel(
        policies,
        ["aads_heuristic", "trained_pr_drl"],
        scenario,
        episodes,
        seed,
        protected_initialization,
        process_ids;
        paper_protocol=paper_protocol,
        progress_every=progress_every,
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

function _table_v(metrics::DataFrame; trained_label::AbstractString="PR-DRL")
    rows = NamedTuple[]
    for policy in ("trained_pr_drl", "aads_heuristic")
        data = metrics[metrics.policy .== policy, :]
        nrow(data) == 0 && continue
        errors = abs.(Float64.(data.target_error_km))
        push!(rows, (
            policy = policy == "trained_pr_drl" ? String(trained_label) : "AADS",
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

function _policy_label(policy; trained_label::AbstractString="PR-DRL")
    policy == "trained_pr_drl" && return String(trained_label)
    policy == "aads_heuristic" && return "AADS"
    return replace(String(policy), "_" => " ")
end

function _paper_figure_10(pass_logs::DataFrame, scenario, path::AbstractString;
                          trained_label::AbstractString="PR-DRL")
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
            label = first_trace ? _policy_label(policy; trained_label=trained_label) : ""
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

_default_flight_policy_specs(trained_label::AbstractString="PR-DRL") = [
    (key="trained_pr_drl", label=String(trained_label)),
    (key="aads_heuristic", label="AADS"),
]

function _flight_performance_table(metrics::DataFrame,
                                   policy_specs=_default_flight_policy_specs())
    rows = NamedTuple[]
    for spec in policy_specs
        policy = spec.key
        label = spec.label
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

function _paper_figure_11(metrics::DataFrame, path::AbstractString,
                          policy_specs=_default_flight_policy_specs())
    policies = getfield.(policy_specs, :key)
    policy_labels = getfield.(policy_specs, :label)
    labels = vcat(policy_labels, ["Mars Odyssey"])
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
        for policy in policies
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
    target_stats = [
        _target_reached_stats(metrics, policy)
        for policy in policies
    ]
    target_reached_percent = getfield.(target_stats, :percent)
    target_panel = bar(
        policy_labels,
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
    target_distance_means = getfield.(target_stats, :mean_abs_error_km)
    target_distance_errors = getfield.(target_stats, :std_abs_error_km)
    target_distance_upper = maximum(target_distance_means .+ target_distance_errors)
    push!(panels, bar(
        policy_labels,
        target_distance_means;
        yerror=target_distance_errors,
        ylabel="Absolute final distance (km)",
        legend=false,
        xrotation=15,
        ylims=(0, max(1.0, 1.1 * target_distance_upper)),
        title="(f) Final target distance",
    ))
    campaigns = minimum(count(==(policy), metrics.policy) for policy in policies)
    savefig(plot(
        panels...;
        layout=(2, 3),
        size=(max(1300, 260 * length(labels)), 820),
        plot_title="$(campaigns) campaigns per policy: thermal-tolerant flight comparison",
        left_margin=8Plots.mm,
        bottom_margin=5Plots.mm,
    ), path)
    return path
end

function _paper_figure_12(pass_logs::DataFrame, scenario, path::AbstractString;
                          trained_label::AbstractString="PR-DRL")
    policy_rows = pass_logs[pass_logs.policy .== "trained_pr_drl", :]
    episode = minimum(policy_rows.episode)
    trace = policy_rows[policy_rows.episode .== episode, :]
    odyssey_heat = vcat(fill(0.20, 45), fill(0.13, 77), fill(0.04, 88))
    p = plot(trace.pass, trace.heat_rate_w_cm2;
             label="$(trained_label) example", linewidth=2, color=:steelblue,
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

function _transition_batch(transitions)
    isempty(transitions) && return nothing
    return (
        observations=hcat(getfield.(transitions, :observation)...),
        next_observations=hcat(getfield.(transitions, :next_observation)...),
        actions=Int[getfield(transition, :action_index) for transition in transitions],
        rewards=Float32[getfield(transition, :reward) for transition in transitions],
        terminated=Bool[getfield(transition, :terminated) for transition in transitions],
        truncated=Bool[getfield(transition, :truncated) for transition in transitions],
    )
end

function _ddqn_evaluation_td_loss(payload, transitions)
    batch = _transition_batch(transitions)
    batch === nothing && return NaN
    online = payload[:online]
    target = get(payload, :target, online)
    config = payload[:config]
    q_values = predict_q(online, batch.observations)
    online_next = predict_q(online, batch.next_observations)
    target_next = predict_q(target, batch.next_observations)
    targets = compute_ddqn_targets(
        online_next,
        target_next,
        batch.rewards,
        batch.terminated,
        batch.truncated,
        config.discount;
        bootstrap_truncated=config.bootstrap_truncated,
    )
    selected = Float32[q_values[action, index]
                       for (index, action) in pairs(batch.actions)]
    return mean(abs2, selected .- targets)
end

function _a2c_evaluation_td_loss(payload, transitions)
    batch = _transition_batch(transitions)
    batch === nothing && return NaN
    critic = payload[:critic]
    config = payload[:config]
    values = vec(predict_q(critic, batch.observations))
    next_values = vec(predict_q(critic, batch.next_observations))
    targets = similar(batch.rewards)
    for index in eachindex(targets)
        bootstrap = batch.terminated[index] ? 0f0 : next_values[index]
        targets[index] = batch.rewards[index] + Float32(config.discount) * bootstrap
    end
    return mean(abs2, values .- targets)
end

function _td3_evaluation_td_loss(payload, transitions, actions_mps)
    batch = _transition_batch(transitions)
    batch === nothing && return NaN
    actions_mps === nothing && throw(ArgumentError(
        "TD3 generalization loss requires exact continuous actions",
    ))
    length(actions_mps) == length(transitions) || throw(DimensionMismatch(
        "TD3 exact-action count must match the transition count",
    ))
    config = payload[:config]
    action_scale = Float32((CONTINUOUS_ACTION_HIGH_MPS - CONTINUOUS_ACTION_LOW_MPS) / 2)
    action_midpoint = Float32((CONTINUOUS_ACTION_HIGH_MPS + CONTINUOUS_ACTION_LOW_MPS) / 2)
    normalized_actions = reshape(
        clamp.((Float32.(actions_mps) .- action_midpoint) ./ action_scale, -1f0, 1f0),
        config.action_dim,
        :,
    )
    target_actor = get(payload, :target_actor, payload[:actor])
    target_critic1 = get(payload, :target_critic1, payload[:critic1])
    target_critic2 = get(payload, :target_critic2, payload[:critic2])
    target_actions = td3_actor_output(target_actor, batch.next_observations)
    target_inputs = vcat(batch.next_observations, target_actions)
    next_q1 = vec(predict_q(target_critic1, target_inputs))
    next_q2 = vec(predict_q(target_critic2, target_inputs))
    targets = similar(batch.rewards)
    for index in eachindex(targets)
        done = batch.terminated[index] ||
               (!config.bootstrap_truncated && batch.truncated[index])
        bootstrap = done ? 0f0 : min(next_q1[index], next_q2[index])
        targets[index] = batch.rewards[index] + Float32(config.discount) * bootstrap
    end
    inputs = vcat(batch.observations, normalized_actions)
    q1 = vec(predict_q(payload[:critic1], inputs))
    q2 = vec(predict_q(payload[:critic2], inputs))
    return (mean(abs2, q1 .- targets) + mean(abs2, q2 .- targets)) / 2
end

function _generalization_evaluation_loss(payload, algorithm::Symbol, transitions;
                                         actions_mps=nothing)
    algorithm in (:ddqn, :pr_drl) &&
        return _ddqn_evaluation_td_loss(payload, transitions)
    algorithm in (:a2c, :a3c) && return _a2c_evaluation_td_loss(payload, transitions)
    algorithm == :td3 && return _td3_evaluation_td_loss(payload, transitions, actions_mps)
    throw(ArgumentError("unsupported generalization-loss algorithm: $algorithm"))
end

_generalization_loss_definition(algorithm::Symbol) =
    algorithm in (:ddqn, :pr_drl) ? "ddqn_action_value_td_mse" :
    algorithm in (:a2c, :a3c) ? "$(algorithm)_critic_value_td_mse" :
    algorithm == :td3 ? "td3_twin_critic_action_value_td_mse" : "unsupported"

function _generalization_record(case_name, result, evaluation_loss, reference_loss)
    summaries = result.summaries
    field_values(field) = Float64[getfield(summary, field) for summary in summaries]
    rewards = field_values(:episode_reward)
    thermal = field_values(:thermal_violations)
    distance = abs.(field_values(:target_error_m)) ./ 1000
    delta_v = field_values(:total_delta_v_mps)
    duration = field_values(:mission_duration_days)
    length_values = field_values(:pass_count)
    maneuvers = field_values(:maneuver_count)
    successes = getfield.(summaries, :success)
    return (
        case=case_name,
        episodes=length(summaries),
        evaluation_td_loss=evaluation_loss,
        generalization_gap=evaluation_loss - reference_loss,
        mean_reward=mean(rewards), std_reward=_std(rewards),
        mean_thermal_violations=mean(thermal), std_thermal_violations=_std(thermal),
        reached_goal_fraction=mean(successes),
        reached_goal_percent=100 * mean(successes),
        mean_goal_distance_km=mean(distance), std_goal_distance_km=_std(distance),
        mean_delta_v_mps=mean(delta_v), std_delta_v_mps=_std(delta_v),
        mean_mission_duration_days=mean(duration), std_mission_duration_days=_std(duration),
        mean_episode_length=mean(length_values), std_episode_length=_std(length_values),
        mean_maneuver_count=mean(maneuvers), std_maneuver_count=_std(maneuvers),
    )
end

function _generalization_case_config_record(case_name, scenario)
    integration = scenario.spaceagora_integration_config
    return (
        case=case_name,
        reference_case=case_name == GENERALIZATION_EVALUATION_REFERENCE_CASE,
        phase=scenario.phase,
        backend_mode=String(scenario.backend_mode),
        atmosphere_model=String(scenario.spaceagora_atmosphere_model),
        gram_wind_mode=String(scenario.spaceagora_gram_wind_mode),
        initial_apoapsis_radius_km=scenario.initial_apoapsis_radius_m / 1000,
        final_apoapsis_radius_km=scenario.final_apoapsis_radius_m / 1000,
        nominal_periapsis_altitude_km=scenario.nominal_periapsis_altitude_m / 1000,
        nominal_argument_of_periapsis_deg=
            rad2deg(scenario.nominal_argument_of_periapsis_rad),
        nominal_epoch=string(scenario.nominal_epoch),
        nominal_initial_conditions=scenario.randomization_config.nominal,
        terminal_on_thermal_violation=
            scenario.termination_config.terminal_on_thermal_violation,
        marsgram_perturbation_scale=
            scenario.randomization_config.marsgram_perturbation_scale,
        mars_mgcm_dust_levels=scenario.spaceagora_mars_mgcm_dust_levels === nothing ?
                              "" : join(scenario.spaceagora_mars_mgcm_dust_levels, ","),
        mars_dust_storm=scenario.spaceagora_mars_dust_storm === nothing ?
                        "" : join(scenario.spaceagora_mars_dust_storm, ","),
        solver_mode=String(integration.solver_mode),
        split_imex_solver=String(integration.split_imex_solver),
        reltol_orbit=integration.reltol_orbit,
        abstol_orbit=integration.abstol_orbit,
        dt_max_orbit_s=integration.dt_max_orbit_s,
        reltol_atmosphere=integration.reltol_atmosphere,
        abstol_atmosphere=integration.abstol_atmosphere,
        dt_max_atmosphere_s=integration.dt_max_atmosphere_s,
    )
end

function _write_generalization_progress(
    path;
    status,
    current_case,
    case_index,
    case_count,
    case_episode,
    episodes_per_case,
    completed_episodes,
    total_episodes,
    started_at,
    worker_processes=1,
    threads_per_process=1,
)
    progress = Dict(
        "status" => String(status),
        "current_case" => String(current_case),
        "case_index" => case_index,
        "case_count" => case_count,
        "case_episode" => case_episode,
        "episodes_per_case" => episodes_per_case,
        "completed_episodes" => completed_episodes,
        "total_episodes" => total_episodes,
        "percent_complete" => total_episodes == 0 ? 100.0 :
                              100 * completed_episodes / total_episodes,
        "elapsed_seconds" => time() - started_at,
        "process_id" => getpid(),
        "worker_processes" => worker_processes,
        "threads_per_process" => threads_per_process,
        "updated_utc" => string(now(UTC)),
    )
    temporary_path = path * ".tmp"
    open(temporary_path, "w") do io
        TOML.print(io, progress)
    end
    mv(temporary_path, path; force=true)
    return path
end

function _render_generalization_report(csv_path::AbstractString,
                                       output_path::AbstractString)
    project_file = Base.active_project()
    project_file === nothing &&
        error("cannot render the generalization report without an active Julia project")
    renderer_path = joinpath(@__DIR__, "render_generalization_table.jl")
    isfile(renderer_path) || error("generalization report renderer is missing: $renderer_path")
    project_dir = dirname(project_file)
    @printf("rendering generalization PDF: %s\n", output_path)
    flush(stdout)
    run(`$(Base.julia_cmd()) --project=$project_dir $renderer_path $csv_path $output_path`)
    isfile(output_path) ||
        error("generalization report renderer did not create: $output_path")
    return String(output_path)
end

function _generalization_table(policy, checkpoint_payload, algorithm, scenario,
                               episodes, seed, protected_initialization, output_dir;
                               progress_every=1,
                               processes=16,
                               threads_per_process=1)
    progress_every > 0 || throw(ArgumentError("progress_every must be positive"))
    processes > 0 || throw(ArgumentError("processes must be positive"))
    threads_per_process > 0 ||
        throw(ArgumentError("threads_per_process must be positive"))
    suite = generalization_evaluation_suite(scenario)
    results = Dict{String,Any}()
    losses = Dict{String,Float64}()
    case_configs = Dict(suite)
    suite_dir = joinpath(output_dir, "generalization_evaluation_suite")
    mkpath(suite_dir)
    progress_path = joinpath(suite_dir, "progress.toml")
    case_count = length(suite)
    total_episodes = case_count * episodes
    started_at = time()
    _write_generalization_progress(
        progress_path;
        status="running",
        current_case="",
        case_index=0,
        case_count=case_count,
        case_episode=0,
        episodes_per_case=episodes,
        completed_episodes=0,
        total_episodes=total_episodes,
        started_at=started_at,
        worker_processes=processes,
        threads_per_process=threads_per_process,
    )
    process_ids = Int[]
    try
        process_ids = _start_evaluation_workers(
            processes,
            threads_per_process,
            scenario,
        )
        for (case_index, (case_name, case_scenario)) in enumerate(suite)
            completed_before_case = (case_index - 1) * episodes
            _write_generalization_progress(
                progress_path;
                status="running",
                current_case=case_name,
                case_index=case_index,
                case_count=case_count,
                case_episode=0,
                episodes_per_case=episodes,
                completed_episodes=completed_before_case,
                total_episodes=total_episodes,
                started_at=started_at,
                worker_processes=length(process_ids),
                threads_per_process=threads_per_process,
            )
            @printf(
                "generalization evaluation case=%s episodes=%d workers=%d\n",
                case_name,
                episodes,
                length(process_ids),
            )
            flush(stdout)
            completed_in_case = Ref(0)
            campaign_callback = function (case_completed, _)
                completed_in_case[] = case_completed
                completed_episodes = completed_before_case + case_completed
                _write_generalization_progress(
                    progress_path;
                    status="running",
                    current_case=case_name,
                    case_index=case_index,
                    case_count=case_count,
                    case_episode=case_completed,
                    episodes_per_case=episodes,
                    completed_episodes=completed_episodes,
                    total_episodes=total_episodes,
                    started_at=started_at,
                    worker_processes=length(process_ids),
                    threads_per_process=threads_per_process,
                )
                if completed_episodes % progress_every == 0 ||
                   case_completed == episodes
                    @printf(
                        "generalization progress case=%s completed=%d/%d overall=%d/%d (%.1f%%)\n",
                        case_name,
                        case_completed,
                        episodes,
                        completed_episodes,
                        total_episodes,
                        100 * completed_episodes / total_episodes,
                    )
                    flush(stdout)
                end
            end
            policy_name = "frozen_policy_$case_name"
            results[case_name] = try
                only_policy = Dict{String,Any}(policy_name => policy)
                parallel_results = _evaluate_policies_parallel(
                    only_policy,
                    [policy_name],
                    case_scenario,
                    episodes,
                    seed,
                    protected_initialization,
                    process_ids;
                    paper_protocol=false,
                    progress_every=progress_every,
                    collect_transitions=true,
                    campaign_callback=campaign_callback,
                )
                parallel_results[policy_name]
            catch
                _write_generalization_progress(
                    progress_path;
                    status="failed",
                    current_case=case_name,
                    case_index=case_index,
                    case_count=case_count,
                    case_episode=completed_in_case[],
                    episodes_per_case=episodes,
                    completed_episodes=completed_before_case + completed_in_case[],
                    total_episodes=total_episodes,
                    started_at=started_at,
                    worker_processes=length(process_ids),
                    threads_per_process=threads_per_process,
                )
                rethrow()
            end
            losses[case_name] = _generalization_evaluation_loss(
                checkpoint_payload,
                algorithm,
                results[case_name].transitions,
                actions_mps=algorithm == :td3 ? Float32[
                    action
                    for summary in results[case_name].summaries
                    for (action, protected) in zip(
                        summary.action_trace,
                        summary.protected_trace,
                    )
                    if !protected
                ] : nothing,
            )
        end
    finally
        if !isempty(process_ids)
            println("stopping generalization evaluation workers")
            flush(stdout)
            rmprocs(process_ids)
        end
    end
    reference_loss = losses[GENERALIZATION_EVALUATION_REFERENCE_CASE]
    all_case_names = first.(suite)
    all_cases = DataFrame([
        _generalization_record(
            case_name,
            results[case_name],
            losses[case_name],
            reference_loss,
        )
        for case_name in all_case_names
    ])
    table_vi = all_cases[in.(all_cases.case, Ref(GENERALIZATION_EVALUATION_CASES)), :]

    episode_rows = NamedTuple[]
    pass_rows = NamedTuple[]
    for case_name in all_case_names
        append!(episode_rows, [
            merge((case=case_name,), episode_metrics(summary; policy_name="frozen_policy"))
            for summary in results[case_name].summaries
        ])
        append!(pass_rows, [
            merge((case=case_name,), row)
            for row in results[case_name].pass_rows
        ])
    end

    table_path = joinpath(suite_dir, "table_vi_metrics.csv")
    all_cases_path = joinpath(suite_dir, "all_cases_with_iid_reference.csv")
    episodes_path = joinpath(suite_dir, "episode_metrics.csv")
    passes_path = joinpath(suite_dir, "pass_logs.csv")
    configs_path = joinpath(suite_dir, "case_configurations.csv")
    report_pdf_path = joinpath(suite_dir, "generalization_results_table.pdf")
    CSV.write(table_path, table_vi)
    CSV.write(all_cases_path, all_cases)
    CSV.write(episodes_path, DataFrame(episode_rows))
    CSV.write(passes_path, DataFrame(pass_rows))
    CSV.write(configs_path, DataFrame([
        _generalization_case_config_record(case_name, case_configs[case_name])
        for case_name in all_case_names
    ]))
    _write_generalization_progress(
        progress_path;
        status="complete",
        current_case=last(all_case_names),
        case_index=case_count,
        case_count=case_count,
        case_episode=episodes,
        episodes_per_case=episodes,
        completed_episodes=total_episodes,
        total_episodes=total_episodes,
        started_at=started_at,
        worker_processes=processes,
        threads_per_process=threads_per_process,
    )
    return (
        table=table_path,
        all_cases=all_cases_path,
        episode_metrics=episodes_path,
        pass_logs=passes_path,
        case_configurations=configs_path,
        progress=progress_path,
        report_pdf=report_pdf_path,
        loss_definition=_generalization_loss_definition(algorithm),
    )
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
        "policy_action_selection" => "greedy",
        "policy_updates_during_evaluation" => false,
        "generalization_protocol" => "paper-inspired SpaceAGORA-native",
        "generalization_loss_definition" => _generalization_loss_definition(algorithm),
        "generalization_reference_case" => GENERALIZATION_EVALUATION_REFERENCE_CASE,
        "seed" => options.seed,
        "iid_episodes" => options.episodes,
        "flight_episodes" => options.flight_episodes,
        "generalization_episodes" => options.generalization_episodes,
        "checkpoint_episodes" => options.checkpoint_episodes,
        "worker_processes" => options.processes,
        "threads_per_process" => options.threads_per_process,
        "gram_wind_mode" => String(options.gram_wind_mode),
        "paper_reference" => Dict(
            "title" => "Autonomous Decision-Making for Aerobraking via Parallel Randomized Deep Reinforcement Learning",
            "doi" => "10.1109/TAES.2022.3221697",
            "odyssey_reference_source" => "Paper Table III, Main Phase II + Endgame + Walkout",
            "figure_12_reference" => "Piecewise phase median heat rates from Paper Table III",
        ),
        "generalization_case_mapping" => Dict(
            "iid_reference" => "held-out draws from the policy run's training distribution",
            "nominal" => "nominal SpaceAGORA initial conditions with the run's native GRAM setup",
            "exponential_density" => "nominal case with only the SpaceAGORA density model changed to exponential",
            "aggressive_atmosphere" => "nominal case with 2x native GRAM perturbations, DUSTTAU 0.3, and a global maximum-intensity storm",
            "short_campaign" => "nominal SpaceAGORA Walkout phase initial conditions",
            "long_campaign" => "nominal SpaceAGORA full-Campaign initial conditions",
            "high_accuracy_spaceagora" => "nominal case with the same SpaceAGORA solver and tighter tolerances and step limits",
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

function evaluate_multi_run_comparison(options)
    requested_run_dirs = [String(options.run_dir); String.(options.comparison_runs)]
    absolute_run_dirs = abspath.(requested_run_dirs)
    length(unique(absolute_run_dirs)) == length(absolute_run_dirs) ||
        throw(ArgumentError("multi-run comparison contains a duplicate run directory"))

    sources = map(eachindex(absolute_run_dirs)) do index
        _load_evaluation_source(
            absolute_run_dirs[index];
            config_override=index == 1 ? options.config : nothing,
            checkpoint_requested=index == 1 ? options.checkpoint : nothing,
            wind_mode=options.wind_mode,
        )
    end
    specs = _comparison_policy_specs(sources)
    primary = first(sources)
    primary_signature = _comparison_environment_signature(primary.resolved)
    environment_matches = Bool[]
    for source in sources
        matches = _comparison_environment_signature(source.resolved) == primary_signature
        push!(environment_matches, matches)
        if !matches
            @warn "comparison run was trained with a different environment; evaluating it on the primary run's scenario for a common test" primary_run=primary.run_dir comparison_run=source.run_dir
        end
    end

    output_dir = abspath(options.output === nothing ?
                        _default_comparison_output_dir(sources) : options.output)
    mkpath(output_dir)
    scenario = _final_flight_scenario(primary.resolved)
    scenario.termination_config.terminal_on_thermal_violation &&
        error("multi-run flight comparison must use thermal-tolerant termination")
    seed = primary.resolved.training.validation_seed
    protected = protected_initialization_config(primary.resolved.training)

    policies = Dict{String,Any}("aads_heuristic" => AADSHeuristicPolicy())
    for spec in specs
        policies[spec.key] = spec.source.policy
        @printf(
            "comparison source policy=%s run=%s algorithm=%s checkpoint=%s config=%s environment_matches_primary=%s\n",
            spec.label,
            spec.source.run_dir,
            spec.source.algorithm,
            spec.source.checkpoint_path,
            spec.source.config_path,
            _comparison_environment_signature(spec.source.resolved) == primary_signature,
        )
    end
    policy_order = vcat(getfield.(specs, :key), ["aads_heuristic"])
    @printf(
        "multi-run flight comparison policies=%d trained_runs=%d campaigns_per_policy=%d scenario_source=%s wind_mode=%s processes=%d threads_per_process=%d output=%s\n",
        length(policy_order),
        length(specs),
        options.episodes,
        primary.run_dir,
        String(primary.resolved.scenario.spaceagora_gram_wind_mode),
        options.processes,
        options.threads_per_process,
        output_dir,
    )
    flush(stdout)

    process_ids = _start_evaluation_workers(
        options.processes,
        options.threads_per_process,
        scenario,
    )
    results = try
        _evaluate_policies_parallel(
            policies,
            policy_order,
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
        joinpath(output_dir, "trained_policies_vs_aads"),
        results,
    )
    metrics, _ = _result_frames(results)
    flight_policy_specs = vcat(
        [(key=spec.key, label=spec.label) for spec in specs],
        [(key="aads_heuristic", label="AADS")],
    )
    comparison = _flight_performance_table(metrics, flight_policy_specs)
    comparison_path = joinpath(
        output_dir,
        "multi_run_flight_performance_comparison.csv",
    )
    CSV.write(comparison_path, comparison)
    for spec in flight_policy_specs
        row = only(eachrow(comparison[comparison.policy .== spec.label, :]))
        @printf(
            "flight comparison target reached policy=%s tolerance_km=%.1f campaigns=%d/%d percent=%.1f%%\n",
            spec.label,
            PAPER_TARGET_TOLERANCE_KM,
            row.target_reached_10km_count,
            row.episodes,
            row.target_reached_10km_percent,
        )
    end
    flush(stdout)
    figure_path = _paper_figure_11(
        metrics,
        joinpath(output_dir, "multi_run_flight_performance_comparison.png"),
        flight_policy_specs,
    )

    source_records = [
        Dict{String,Any}(
            "run_dir" => source.run_dir,
            "run_id" => source.run_id,
            "algorithm" => String(source.algorithm),
            "policy_key" => spec.key,
            "policy_label" => spec.label,
            "checkpoint_path" => abspath(source.checkpoint_path),
            "checkpoint_global_step" => Int(get(source.checkpoint_payload, :global_step, 0)),
            "config_path" => abspath(source.config_path),
            "config_sha256" => source.config_sha256,
            "environment_matches_primary" => environment_matches[index],
        )
        for (index, (source, spec)) in enumerate(zip(sources, specs))
    ]
    evaluation_manifest = joinpath(output_dir, "evaluation_manifest.toml")
    manifest = Dict{String,Any}(
        "evaluation_mode" => "multi_run_flight_comparison",
        "scenario_source_run_dir" => primary.run_dir,
        "campaigns_per_policy" => options.episodes,
        "worker_processes" => options.processes,
        "threads_per_process" => options.threads_per_process,
        "progress_every_campaigns" => options.progress_every,
        "seed" => seed,
        "terminal_on_thermal_violation" => false,
        "target_reached_tolerance_km" => PAPER_TARGET_TOLERANCE_KM,
        "policy_action_selection" => "greedy",
        "gram_wind_mode" => String(primary.resolved.scenario.spaceagora_gram_wind_mode),
        "source_runs" => source_records,
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
    println("multi-run flight comparison complete: ", output_dir)
    println("comparison figure: ", figure_path)
    println("comparison metrics: ", comparison_path)
    return (
        output_dir=output_dir,
        checkpoints=getfield.(sources, :checkpoint_path),
        metrics=comparison_path,
        figure=figure_path,
        manifest=evaluation_manifest,
    )
end

function evaluate_final_flight_comparison(options)
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

    resolved = resolve_config(config_path; gram_wind_mode=options.wind_mode)
    checkpoint_path = _checkpoint_path(run_dir, options.checkpoint)
    policy, algorithm, _ = _load_policy(checkpoint_path)
    trained_label = SpaceAGORA_RL.algorithm_display_name(algorithm)
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
        "final flight comparison run=%s algorithm=%s checkpoint=%s campaigns_per_policy=%d thermal_terminal=false wind_mode=%s processes=%d threads_per_process=%d progress_every=%d config=%s output=%s\n",
        run_dir,
        algorithm,
        checkpoint_path,
        options.episodes,
        String(resolved.scenario.spaceagora_gram_wind_mode),
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
    flight_policy_specs = _default_flight_policy_specs(trained_label)
    comparison = _flight_performance_table(metrics, flight_policy_specs)
    CSV.write(comparison_path, comparison)
    for spec in flight_policy_specs
        row = only(eachrow(comparison[comparison.policy .== spec.label, :]))
        @printf(
            "flight comparison target reached policy=%s tolerance_km=%.1f campaigns=%d/%d percent=%.1f%%\n",
            spec.label,
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
        flight_policy_specs,
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
        "gram_wind_mode" => String(resolved.scenario.spaceagora_gram_wind_mode),
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
    !isempty(options.comparison_runs) &&
        return evaluate_multi_run_comparison(options)
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
    resolved = resolve_config(config_path; gram_wind_mode=options.wind_mode)
    checkpoint_path = _checkpoint_path(run_dir, options.checkpoint)
    policy, algorithm, checkpoint_payload = _load_policy(checkpoint_path)
    trained_label = SpaceAGORA_RL.algorithm_display_name(algorithm)
    default_output_name = options.generalization_only ?
                          "generalization_evaluation" : "paper_evaluation"
    output_dir = abspath(options.output === nothing ?
                        joinpath(run_dir, default_output_name) : options.output)
    mkpath(output_dir)
    seed = resolved.training.validation_seed
    protected = protected_initialization_config(resolved.training)
    scenario = paper_evaluation_scenario(
        resolved.scenario;
        max_passes=max(1000, resolved.scenario.termination_config.max_passes),
    )
    runtime_options = merge(
        options,
        (; seed=seed, gram_wind_mode=resolved.scenario.spaceagora_gram_wind_mode),
    )

    @printf("evaluating run=%s algorithm=%s checkpoint=%s wind_mode=%s config=%s output=%s\n",
            run_dir, algorithm, checkpoint_path,
            String(resolved.scenario.spaceagora_gram_wind_mode), config_path, output_dir)

    if options.generalization_only
        paths = _generalization_table(
            policy,
            checkpoint_payload,
            algorithm,
            scenario,
            options.generalization_episodes,
            seed,
            protected,
            output_dir,
            progress_every=options.progress_every,
            processes=options.processes,
            threads_per_process=options.threads_per_process,
        )
        artifacts = (
            generalization_table_vi=paths.table,
            paper_table_vi=paths.table,
            generalization_all_cases=paths.all_cases,
            generalization_episode_metrics=paths.episode_metrics,
            generalization_pass_logs=paths.pass_logs,
            generalization_case_configurations=paths.case_configurations,
            generalization_progress=paths.progress,
            generalization_report_pdf=paths.report_pdf,
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
        _render_generalization_report(paths.all_cases, paths.report_pdf)
        println("generalization evaluation complete: ", paths.table)
        println("generalization PDF: ", paths.report_pdf)
        return (output_dir=output_dir, artifacts=artifacts, manifest=evaluation_manifest)
    end

    iid_dir = joinpath(output_dir, "iid_pr_drl_vs_aads")
    iid = _evaluate_pair(policy, scenario, options.episodes, seed, protected)
    iid_paths = write_evaluation_artifacts(iid_dir, iid)
    iid_metrics, iid_pass_logs = _result_frames(iid)
    table_v_path = joinpath(output_dir, "paper_table_v_pr_drl_vs_aads.csv")
    CSV.write(table_v_path, _table_v(iid_metrics; trained_label=trained_label))
    fig10_path = _paper_figure_10(
        iid_pass_logs,
        scenario,
        joinpath(output_dir, "paper_fig10_heat_rate_and_delta_v.png"),
        trained_label=trained_label,
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
    flight_policy_specs = _default_flight_policy_specs(trained_label)
    fig11_path = _paper_figure_11(
        flight_metrics,
        joinpath(output_dir, "paper_fig11_flight_performance_summary.png"),
        flight_policy_specs,
    )
    fig12_path = _paper_figure_12(
        flight_pass_logs,
        flight_scenario,
        joinpath(output_dir, "paper_fig12_heat_rate_odyssey_comparison.png"),
        trained_label=trained_label,
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

    generalization_paths = (
        table="",
        all_cases="",
        episode_metrics="",
        pass_logs="",
        case_configurations="",
        progress="",
        report_pdf="",
    )
    if options.generalization
        generalization_paths = _generalization_table(
            policy,
            checkpoint_payload,
            algorithm,
            scenario,
            options.generalization_episodes,
            seed,
            protected,
            output_dir,
            progress_every=options.progress_every,
            processes=options.processes,
            threads_per_process=options.threads_per_process,
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
        generalization_table_vi=generalization_paths.table,
        paper_table_vi=generalization_paths.table,
        generalization_all_cases=generalization_paths.all_cases,
        generalization_episode_metrics=generalization_paths.episode_metrics,
        generalization_pass_logs=generalization_paths.pass_logs,
        generalization_case_configurations=generalization_paths.case_configurations,
        generalization_progress=generalization_paths.progress,
        generalization_report_pdf=generalization_paths.report_pdf,
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
    if options.generalization
        _render_generalization_report(
            generalization_paths.all_cases,
            generalization_paths.report_pdf,
        )
        println("generalization PDF: ", generalization_paths.report_pdf)
    end
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
