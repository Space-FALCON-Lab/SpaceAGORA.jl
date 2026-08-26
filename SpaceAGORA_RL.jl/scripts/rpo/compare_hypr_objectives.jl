#!/usr/bin/env julia

include(joinpath(@__DIR__, "evaluate_hypr_baseline.jl"))

function _run_hypr_comparison_cases(config, scenarios, seeds, n_workers::Int)
    n_cases = length(scenarios)
    results = Vector{Any}(undef, n_cases)
    active_workers = min(max(n_workers, 1), n_cases)
    println(
        "comparing original and retimed-fuel HYPR " *
        "cases=$n_cases active_workers=$active_workers terminal=full_lqmpc",
    )
    if active_workers == 1
        for index in eachindex(scenarios)
            results[index] = evaluate_hypr_pso_comparison_case(
                config, scenarios[index], seeds[index],
            )
            println("HYPR comparison progress $index/$n_cases")
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
                evaluate_hypr_pso_comparison_case,
                process_ids[mod1(offset, active_workers)],
                config, scenarios[index], seeds[index],
            ) for (offset, index) in enumerate(indices)]
            for (index, future) in zip(indices, futures)
                results[index] = try
                    fetch(future)
                catch error
                    failure = (
                        plan=nothing, runtime_s=NaN, planner_runtime_s=NaN,
                        terminal_runtime_s=NaN, error=sprint(showerror, error),
                    )
                    (original=failure, retimed_fuel=failure)
                end
                println("HYPR comparison progress $index/$n_cases")
            end
        end
    finally
        rmprocs(process_ids)
    end
    return results
end

function _paired_hypr_metrics(original_rows, retimed_rows)
    return [begin
        original = original_rows[index, :]
        retimed = retimed_rows[index, :]
        both_successful = original.success && retimed.success
        finite_fuel = isfinite(original.propellant_used_g) &&
            isfinite(retimed.propellant_used_g)
        fuel_savings_g = finite_fuel ?
            original.propellant_used_g - retimed.propellant_used_g : NaN
        fuel_savings_percent = finite_fuel && original.propellant_used_g > 0.0 ?
            100.0 * fuel_savings_g / original.propellant_used_g : NaN
        runtime_ratio = isfinite(original.runtime_s) && original.runtime_s > 0.0 &&
            isfinite(retimed.runtime_s) ? retimed.runtime_s / original.runtime_s : NaN
        planner_runtime_ratio = isfinite(original.planner_runtime_s) &&
            original.planner_runtime_s > 0.0 &&
            isfinite(retimed.planner_runtime_s) ?
            retimed.planner_runtime_s / original.planner_runtime_s : NaN
        (
            case=original.case,
            both_successful=both_successful,
            original_success=original.success,
            retimed_fuel_success=retimed.success,
            original_propellant_used_g=original.propellant_used_g,
            retimed_fuel_propellant_used_g=retimed.propellant_used_g,
            fuel_savings_g=fuel_savings_g,
            fuel_savings_percent=fuel_savings_percent,
            original_minimum_clearance_m=original.minimum_clearance_m,
            retimed_fuel_minimum_clearance_m=retimed.minimum_clearance_m,
            original_safety_violations=original.safety_violation_count,
            retimed_fuel_safety_violations=retimed.safety_violation_count,
            original_reference_path_length_m=original.reference_path_length_m,
            retimed_fuel_reference_path_length_m=retimed.reference_path_length_m,
            original_duration_s=original.duration_s,
            retimed_fuel_duration_s=retimed.duration_s,
            original_runtime_s=original.runtime_s,
            retimed_fuel_runtime_s=retimed.runtime_s,
            runtime_ratio=runtime_ratio,
            original_planner_runtime_s=original.planner_runtime_s,
            retimed_fuel_planner_runtime_s=retimed.planner_runtime_s,
            planner_runtime_ratio=planner_runtime_ratio,
            original_terminal_runtime_s=original.terminal_runtime_s,
            retimed_fuel_terminal_runtime_s=retimed.terminal_runtime_s,
            original_pso_iterations=original.pso_iterations,
            retimed_fuel_pso_iterations=retimed.pso_iterations,
            original_pso_objective=original.edit_objective,
            retimed_fuel_pso_objective=retimed.edit_objective,
        )
    end for index in 1:nrow(original_rows)]
end

function _finite_median(values)
    finite = filter(isfinite, Float64.(collect(skipmissing(values))))
    return isempty(finite) ? NaN : median(finite)
end

function _paired_hypr_summary(rows::DataFrame)
    paired = rows[rows.both_successful, :]
    savings = paired.fuel_savings_g
    tolerance_g = 1.0e-9
    return (
        cases=nrow(rows),
        paired_success_count=nrow(paired),
        original_success_count=count(rows.original_success),
        retimed_fuel_success_count=count(rows.retimed_fuel_success),
        mean_paired_fuel_savings_g=_finite_mean(savings),
        median_paired_fuel_savings_g=_finite_median(savings),
        mean_paired_fuel_savings_percent=
            _finite_mean(paired.fuel_savings_percent),
        retimed_fuel_lower_count=count(>(tolerance_g), savings),
        original_lower_count=count(<(-tolerance_g), savings),
        equal_fuel_count=count(value -> abs(value) <= tolerance_g, savings),
        mean_runtime_ratio=_finite_mean(rows.runtime_ratio),
        mean_planner_runtime_ratio=_finite_mean(rows.planner_runtime_ratio),
    )
end

function _write_paired_case(path, case_index, original, retimed)
    title = @sprintf("HYPR objective comparison case %03d", case_index)
    open(path, "w") do io
        print(io, """<!doctype html><html><head><meta charset=\"utf-8\"><title>$(title)</title>
<style>html,body{height:100%;margin:0;font:14px system-ui;color:#18212b}.header{height:4rem;box-sizing:border-box;padding:.6rem 1rem;background:#f2f5f7}.plots{display:grid;grid-template-columns:1fr 1fr;height:calc(100% - 4rem);gap:2px;background:#bcc5ce}.panel{background:white;display:flex;flex-direction:column}.label{text-align:center;padding:.35rem;font-weight:600}.panel iframe{border:0;flex:1;width:100%}</style>
</head><body><div class=\"header\"><b>$(title)</b><br>Original: $(@sprintf("%.5f g", original.propellant_used_g)) fuel, $(@sprintf("%.3f s", original.runtime_s)) runtime · Retimed fuel: $(@sprintf("%.5f g", retimed.propellant_used_g)) fuel, $(@sprintf("%.3f s", retimed.runtime_s)) runtime</div>
<div class=\"plots\"><div class=\"panel\"><div class=\"label\">Original HyPR proxy objective</div><iframe src=\"../original/trajectories/case_$(@sprintf("%03d", case_index)).html\"></iframe></div>
<div class=\"panel\"><div class=\"label\">HyPR retimed-fuel objective</div><iframe src=\"../retimed_fuel/trajectories/case_$(@sprintf("%03d", case_index)).html\"></iframe></div></div></body></html>""")
    end
    return path
end

function _write_hypr_comparison_index(path, rows::DataFrame, summary)
    open(path, "w") do io
        print(io, """<!doctype html><html><head><meta charset=\"utf-8\"><title>HYPR objective comparison</title>
<style>body{font:14px system-ui;margin:2rem;color:#18212b}table{border-collapse:collapse;width:100%}th,td{padding:.45rem;border-bottom:1px solid #d9dfe5;text-align:right}th{background:#f2f5f7;position:sticky;top:0}td:first-child,th:first-child{text-align:left}.positive{color:#16834a}.negative{color:#c0392b}.summary{columns:2;margin-bottom:2rem}</style>
</head><body><h1>Original versus retimed-fuel HYPR</h1><div class=\"summary\">""")
        for name in propertynames(summary)
            print(io, "<div><b>$(_html_escape(name)):</b> $(_html_escape(getproperty(summary, name)))</div>")
        end
        print(io, """</div><p><a href=\"paired_cases.csv\">Paired CSV</a> · <a href=\"comparison_summary.csv\">Comparison summary</a> · <a href=\"planner_summaries.csv\">Planner summaries</a></p>
<table><thead><tr><th>case</th><th>original success</th><th>retimed success</th><th>original fuel (g)</th><th>retimed fuel (g)</th><th>fuel savings (g)</th><th>original runtime (s)</th><th>retimed runtime (s)</th><th>runtime ratio</th></tr></thead><tbody>""")
        for row in eachrow(rows)
            savings_class = row.fuel_savings_g >= 0.0 ? "positive" : "negative"
            print(io, @sprintf(
                "<tr><td><a href=\"case_comparisons/case_%03d.html\">case %03d</a></td><td>%s</td><td>%s</td><td>%.5f</td><td>%.5f</td><td class=\"%s\">%.5f</td><td>%.3f</td><td>%.3f</td><td>%.2f</td></tr>",
                row.case, row.case, row.original_success,
                row.retimed_fuel_success, row.original_propellant_used_g,
                row.retimed_fuel_propellant_used_g, savings_class,
                row.fuel_savings_g, row.original_runtime_s,
                row.retimed_fuel_runtime_s, row.runtime_ratio,
            ))
        end
        print(io, "</tbody></table></body></html>")
    end
    return path
end

function comparison_main(args=ARGS)
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
    n_cases > 0 || throw(ArgumentError("number of comparison cases must be positive"))
    output_directory = length(args) >= 2 ? args[2] : String(get(
        evaluation_config,
        "hypr_comparison_output_directory",
        "SpaceAGORA_RL.jl/outputs/hypr_rl/evaluation_hypr_objective_comparison",
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
    sampler = build_rpo_hypr_rl_endpoint_sampler(
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
        sampler, MersenneTwister(scenario_seeds[index]),
    ) for index in 1:n_cases]
    results = _run_hypr_comparison_cases(
        config, scenarios, planner_seeds, n_workers,
    )
    original_results = [result.original for result in results]
    retimed_results = [result.retimed_fuel for result in results]
    original_metrics = [_hypr_baseline_case_metrics(
        index, scenarios[index], original_results[index], config,
    ) for index in 1:n_cases]
    retimed_metrics = [_hypr_baseline_case_metrics(
        index, scenarios[index], retimed_results[index], config,
    ) for index in 1:n_cases]
    original_rows = DataFrame(original_metrics)
    retimed_rows = DataFrame(retimed_metrics)
    original_summary = _evaluation_summary(
        original_rows, nothing, Dict{Symbol, Any}(),
    )
    retimed_summary = _evaluation_summary(
        retimed_rows, nothing, Dict{Symbol, Any}(),
    )
    paired_rows = DataFrame(_paired_hypr_metrics(original_rows, retimed_rows))
    comparison_summary = _paired_hypr_summary(paired_rows)

    original_directory = joinpath(output_directory, "original")
    retimed_directory = joinpath(output_directory, "retimed_fuel")
    comparison_directory = joinpath(output_directory, "case_comparisons")
    for directory in (
        joinpath(original_directory, "trajectories"),
        joinpath(retimed_directory, "trajectories"),
        comparison_directory,
    )
        mkpath(directory)
    end
    mesh_path = _write_station_mesh_asset(
        joinpath(output_directory, "station_mesh.js"),
        _load_station_mesh(station_asset),
        station_asset,
    )
    CSV.write(joinpath(original_directory, "cases.csv"), original_rows)
    CSV.write(
        joinpath(original_directory, "summary.csv"),
        DataFrame([original_summary]),
    )
    CSV.write(joinpath(retimed_directory, "cases.csv"), retimed_rows)
    CSV.write(
        joinpath(retimed_directory, "summary.csv"),
        DataFrame([retimed_summary]),
    )
    CSV.write(joinpath(output_directory, "paired_cases.csv"), paired_rows)
    CSV.write(
        joinpath(output_directory, "comparison_summary.csv"),
        DataFrame([comparison_summary]),
    )
    CSV.write(
        joinpath(output_directory, "planner_summaries.csv"),
        DataFrame([
            merge((planner="original_hypr_proxy",), original_summary),
            merge((planner="hypr_retimed_fuel",), retimed_summary),
        ]),
    )
    bezier_plot_spacing_m = Float64(get(
        evaluation_config, "bezier_plot_spacing_m", 0.05,
    ))
    for index in 1:n_cases
        _write_case_plot(
            joinpath(original_directory, "trajectories", @sprintf("case_%03d.html", index)),
            index, scenarios[index], original_results[index], original_metrics[index], config;
            station_mesh_script="../../$(basename(mesh_path))",
            bezier_plot_spacing_m=bezier_plot_spacing_m,
            planner_label="Original HYPR (proxy objective)",
        )
        _write_case_plot(
            joinpath(retimed_directory, "trajectories", @sprintf("case_%03d.html", index)),
            index, scenarios[index], retimed_results[index], retimed_metrics[index], config;
            station_mesh_script="../../$(basename(mesh_path))",
            bezier_plot_spacing_m=bezier_plot_spacing_m,
            planner_label="HYPR (retimed-fuel objective)",
        )
        _write_paired_case(
            joinpath(comparison_directory, @sprintf("case_%03d.html", index)),
            index, original_rows[index, :], retimed_rows[index, :],
        )
    end
    _write_index(
        joinpath(original_directory, "index.html"), original_rows, original_summary;
        planner_label="Original HYPR (proxy objective)",
    )
    _write_index(
        joinpath(retimed_directory, "index.html"), retimed_rows, retimed_summary;
        planner_label="HYPR (retimed-fuel objective)",
    )
    index_path = _write_hypr_comparison_index(
        joinpath(output_directory, "index.html"), paired_rows, comparison_summary,
    )
    manifest = Dict(
        "config" => abspath(config_path),
        "cases" => n_cases,
        "evaluation_seed" => evaluation_seed,
        "scenario_seeds" => scenario_seeds,
        "planner_seeds" => planner_seeds,
        "workers" => min(max(n_workers, 1), n_cases),
        "pso_particles" => pso_config.n_particles,
        "pso_iterations" => pso_config.n_iters,
        "pso_waypoints" => pso_config.n_waypoints,
        "original_objective" => "hypr_length_obstacle_proxy_fuel",
        "current_objective" => "retimed_feedforward_six_thruster_fuel_plus_wheel",
        "terminal_evaluator" => "full_lqmpc_six_thruster_attitude_for_both",
        "index_html" => abspath(index_path),
    )
    open(joinpath(output_directory, "evaluation_manifest.toml"), "w") do io
        TOML.print(io, manifest; sorted=true)
    end
    println("wrote paired HYPR comparison to ", abspath(output_directory))
    println("open ", abspath(index_path))
    return (
        original_rows=original_rows,
        retimed_rows=retimed_rows,
        paired_rows=paired_rows,
        comparison_summary=comparison_summary,
        results=results,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    comparison_main()
end
