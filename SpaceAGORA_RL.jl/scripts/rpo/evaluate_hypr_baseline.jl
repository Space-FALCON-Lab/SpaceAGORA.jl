#!/usr/bin/env julia

include(joinpath(@__DIR__, "evaluate_hypr_rl.jl"))

function _hypr_baseline_task_config(task)
    return RPOHyPRRLConfig(
        safe_distance_m=task["safe_distance_m"],
        max_translation_waypoints=task["max_translation_waypoints"],
        max_attitude_waypoints=task["max_attitude_waypoints"],
        max_edits=task["max_edits"],
        coordinate_scale_m=task["coordinate_scale_m"],
        fuel_weight=task["fuel_weight"],
        duration_weight=task["duration_weight"],
        allocation_error_weight=task["allocation_error_weight"],
        wheel_weight=task["wheel_weight"],
    )
end

function _run_hypr_baseline_cases(config, scenarios, seeds, n_workers::Int)
    n_cases = length(scenarios)
    results = Vector{Any}(undef, n_cases)
    active_workers = min(max(n_workers, 1), n_cases)
    println(
        "evaluating baseline HYPR with PSO retimed-fuel objective " *
        "cases=$n_cases active_workers=$active_workers terminal=full_lqmpc",
    )
    if active_workers == 1
        for index in eachindex(scenarios)
            results[index] = try
                evaluate_hypr_pso_baseline_case(
                    config, scenarios[index], seeds[index],
                )
            catch error
                (plan=nothing, runtime_s=NaN, error=sprint(showerror, error))
            end
            println("baseline evaluation progress $index/$n_cases")
        end
        return results
    end

    process_ids = SpaceAGORA_RL._setup_hypr_rl_process_workers(active_workers)
    try
        foreach(fetch, [remotecall(
            SpaceAGORA_RL._prepare_hypr_rl_process_worker!, process_id,
        ) for process_id in process_ids])
        for first_index in 1:active_workers:n_cases
            last_index = min(first_index + active_workers - 1, n_cases)
            indices = first_index:last_index
            futures = [remotecall(
                evaluate_hypr_pso_baseline_case,
                process_ids[mod1(offset, active_workers)],
                config, scenarios[index], seeds[index],
            ) for (offset, index) in enumerate(indices)]
            for (index, future) in zip(indices, futures)
                results[index] = try
                    fetch(future)
                catch error
                    (plan=nothing, runtime_s=NaN, error=sprint(showerror, error))
                end
                println("baseline evaluation progress $index/$n_cases")
            end
        end
    finally
        rmprocs(process_ids)
    end
    return results
end

function _hypr_baseline_case_metrics(case_index, scenario, result, config)
    metrics = _case_metrics(case_index, scenario, result, config)
    plan = _ntget(result, :plan, nothing)
    diagnostics = plan === nothing ? nothing : plan.diagnostics
    return merge(
        metrics,
        (
            planner_runtime_s=Float64(_ntget(
                result, :planner_runtime_s, NaN,
            )),
            terminal_runtime_s=Float64(_ntget(
                result, :terminal_runtime_s, NaN,
            )),
            pso_iterations=Int(_ntget(
                diagnostics, :pso_iterations, 0,
            )),
            pso_refinement_improved=Bool(_ntget(
                diagnostics, :pso_refinement_improved, false,
            )),
            pso_early_stopped=Bool(_ntget(
                diagnostics, :pso_early_stopped, false,
            )),
        ),
    )
end

function baseline_main(args=ARGS)
    config_path = isempty(args) ?
        joinpath(@__DIR__, "..", "..", "configs", "rpo", "hypr_rl.toml") :
        args[1]
    raw = TOML.parsefile(config_path)
    task = raw["task"]
    scenario_config = raw["scenario"]
    training_config = raw["training"]
    evaluation_config = raw["evaluation"]
    config = _hypr_baseline_task_config(task)
    n_cases = length(args) >= 3 ? parse(Int, args[3]) :
        Int(get(evaluation_config, "cases", 100))
    n_cases > 0 || throw(ArgumentError("number of baseline cases must be positive"))
    output_directory = length(args) >= 2 ? args[2] : String(get(
        evaluation_config,
        "hypr_baseline_output_directory",
        "SpaceAGORA_RL.jl/outputs/hypr_rl/evaluation_hypr_pso_retimed_fuel",
    ))
    n_workers = length(args) >= 4 ? parse(Int, args[4]) :
        Int(get(evaluation_config, "n_workers", training_config["n_workers"]))
    evaluation_seed = Int(get(
        evaluation_config, "seed", training_config["seed"] + 10_000_000,
    ))
    station_asset = Symbol(scenario_config["station_asset"])
    modules = SpaceAGORA_RL._spaceagora_rpo_modules()
    pso_config = Base.invokelatest(
        getproperty(modules.guidance, :rpo_740_mpc_final_pso_config);
        safe_distance_m=config.safe_distance_m,
        n_particles=Int(evaluation_config["hypr_baseline_particles"]),
        n_iters=Int(evaluation_config["hypr_baseline_iterations"]),
        n_waypoints=Int(evaluation_config["hypr_baseline_waypoints"]),
        refinement_enable=true,
    )
    base_scenario = build_rpo_hypr_rl_scenario(
        station_asset=station_asset,
        station_points=scenario_config["station_points"],
        station_seed=scenario_config["station_seed"],
        station_keepout_radius_m=scenario_config["station_keepout_radius_m"],
        pso_config=pso_config,
    )
    scenario_sampler = build_rpo_hypr_rl_endpoint_sampler(
        base_scenario;
        station_asset=station_asset,
        safe_distance_m=config.safe_distance_m,
        endpoint_clearance_margin_m=
            scenario_config["endpoint_clearance_margin_m"],
        endpoint_max_clearance_m=scenario_config["endpoint_max_clearance_m"],
        min_separation_m=scenario_config["min_separation_m"],
        surrounded_max_distance_m=scenario_config["surrounded_max_distance_m"],
        max_sampling_tries=scenario_config["max_sampling_tries"],
    )
    scenario_seeds = [evaluation_seed + 2 * index for index in 1:n_cases]
    planner_seeds = [evaluation_seed + 2 * index + 1 for index in 1:n_cases]
    scenarios = [sample_rpo_hypr_rl_scenario(
        scenario_sampler, MersenneTwister(scenario_seeds[index]),
    ) for index in 1:n_cases]

    mkpath(output_directory)
    trajectory_directory = joinpath(output_directory, "trajectories")
    mkpath(trajectory_directory)
    station_mesh = _load_station_mesh(station_asset)
    mesh_path = _write_station_mesh_asset(
        joinpath(output_directory, "station_mesh.js"), station_mesh, station_asset,
    )
    results = _run_hypr_baseline_cases(
        config, scenarios, planner_seeds, n_workers,
    )
    metrics = [_hypr_baseline_case_metrics(
        index, scenarios[index], results[index], config,
    ) for index in 1:n_cases]
    rows = DataFrame(metrics)
    summary = merge(
        _evaluation_summary(rows, nothing, Dict{Symbol, Any}()),
        (
            mean_planner_runtime_s=_finite_mean(rows.planner_runtime_s),
            mean_terminal_runtime_s=_finite_mean(rows.terminal_runtime_s),
            mean_pso_iterations=_finite_mean(rows.pso_iterations),
            pso_refinement_improved_count=count(rows.pso_refinement_improved),
            pso_early_stopped_count=count(rows.pso_early_stopped),
        ),
    )
    CSV.write(joinpath(output_directory, "cases.csv"), rows)
    CSV.write(joinpath(output_directory, "summary.csv"), DataFrame([summary]))
    bezier_plot_spacing_m = Float64(get(
        evaluation_config, "bezier_plot_spacing_m", 0.05,
    ))
    for index in 1:n_cases
        _write_case_plot(
            joinpath(trajectory_directory, @sprintf("case_%03d.html", index)),
            index, scenarios[index], results[index], metrics[index], config;
            station_mesh_script="../$(basename(mesh_path))",
            bezier_plot_spacing_m=bezier_plot_spacing_m,
            planner_label="HYPR PSO (retimed fuel)",
        )
    end
    index_path = _write_index(
        joinpath(output_directory, "index.html"), rows, summary;
        planner_label="HYPR PSO (retimed fuel)",
    )
    manifest = Dict(
        "config" => abspath(config_path),
        "cases" => n_cases,
        "evaluation_seed" => evaluation_seed,
        "scenario_seeds" => scenario_seeds,
        "planner_seeds" => planner_seeds,
        "workers" => min(max(n_workers, 1), n_cases),
        "safe_distance_m" => config.safe_distance_m,
        "station_points" => scenario_config["station_points"],
        "pso_particles" => pso_config.n_particles,
        "pso_iterations" => pso_config.n_iters,
        "pso_waypoints" => pso_config.n_waypoints,
        "pso_objective" => "retimed_feedforward_six_thruster_fuel_plus_wheel",
        "terminal_evaluator" => "full_lqmpc_six_thruster_attitude",
        "index_html" => abspath(index_path),
    )
    open(joinpath(output_directory, "evaluation_manifest.toml"), "w") do io
        TOML.print(io, manifest; sorted=true)
    end
    println("wrote retimed-fuel HYPR baseline to ", abspath(output_directory))
    println("open ", abspath(index_path))
    return (rows=rows, summary=summary, results=results)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    baseline_main()
end
