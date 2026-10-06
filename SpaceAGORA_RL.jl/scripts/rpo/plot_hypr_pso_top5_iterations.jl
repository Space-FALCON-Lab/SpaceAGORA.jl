#!/usr/bin/env julia

include(joinpath(@__DIR__, "compare_hypr_objectives.jl"))

function _hypr_top5_recorder(position_to_path)
    history = Any[]
    callback = function (iteration, pbest, pbest_cost, gbest, gbest_cost,
                         start, goal, n_waypoints, _)
        order = sortperm(pbest_cost)
        positions = Vector{Vector{Float64}}()
        costs = Float64[]
        if isfinite(gbest_cost)
            push!(positions, Vector{Float64}(gbest))
            push!(costs, Float64(gbest_cost))
        end
        for index in order
            isfinite(pbest_cost[index]) || continue
            candidate = Vector{Float64}(view(pbest, :, index))
            any(existing -> existing == candidate, positions) && continue
            push!(positions, candidate)
            push!(costs, Float64(pbest_cost[index]))
            length(positions) >= 5 && break
        end
        paths = [Matrix{Float64}(Base.invokelatest(
            position_to_path, position, start, goal, n_waypoints,
        )) for position in positions]
        push!(history, (
            iteration=Int(iteration),
            costs=costs,
            paths=paths,
        ))
        return nothing
    end
    return history, callback
end

function _run_hypr_top5_history(scenario, config, pso_config, planner_seed;
                                 retimed::Bool)
    modules = SpaceAGORA_RL._spaceagora_rpo_modules()
    position_to_path = getproperty(modules.guidance, :rpo_position_to_path)
    history, callback = _hypr_top5_recorder(position_to_path)
    objective_evaluator = retimed ?
        RPOHyPRRLPSOObjectiveEvaluator(config, scenario) : nothing
    active_config = if retimed
        _, retimed_config, _, _ = SpaceAGORA_RL._rpo_spaceagora_settings(
            scenario, config,
        )
        retimed_config
    else
        Base.invokelatest(
            getproperty(modules.guidance, :rpo_pso_config), pso_config;
            safe_distance_m=config.safe_distance_m,
        )
    end
    elapsed_s = @elapsed result = Base.invokelatest(
        getproperty(modules.guidance, :rpo_pso_plan_path),
        scenario.start_rtn,
        scenario.goal_rtn,
        scenario.geometry,
        active_config;
        safe_distance_m=config.safe_distance_m,
        rng=MersenneTwister(planner_seed),
        swarm_callback=callback,
        objective_evaluator=objective_evaluator,
    )
    isempty(history) && throw(ErrorException("PSO produced no iteration history"))
    return history, result, elapsed_s
end

function _hypr_history_frame(history, iteration::Int)
    index = findlast(frame -> frame.iteration <= iteration, history)
    return history[index === nothing ? 1 : index]
end

function _hypr_sample_top_path(path, sample_path, spacing_m::Real;
                               curve_type::Symbol=:bezier)
    sampled = Matrix{Float64}(Base.invokelatest(
        sample_path, path, Float64(spacing_m); curve_type=curve_type,
    ))
    size(sampled, 2) <= 240 && return sampled
    indices = unique(round.(Int, range(1, size(sampled, 2); length=240)))
    return sampled[:, indices]
end

function _hypr_reference_minimum_clearance(path, scenario;
                                           curve_type::Symbol)
    path === nothing && return NaN
    modules = SpaceAGORA_RL._spaceagora_rpo_modules()
    samples = Base.invokelatest(
        getproperty(modules.guidance, :rpo_sample_path), path, 0.02;
        curve_type=curve_type,
    )
    stats = Base.invokelatest(
        getproperty(modules.guidance, :rpo_clearance_stats_from_samples),
        samples, scenario.geometry, 0.0,
    )
    return stats.min_clearance
end

function _hypr_top5_trace(path, cost, rank::Int, scene, sample_path,
                          spacing_m::Real, color)
    sampled = _hypr_sample_top_path(path, sample_path, spacing_m)
    return Dict(
        "type" => "scatter3d",
        "mode" => "lines",
        "scene" => scene,
        "name" => @sprintf("rank %d · %.6g", rank, cost),
        "x" => vec(sampled[1, :]),
        "y" => vec(sampled[2, :]),
        "z" => vec(sampled[3, :]),
        "line" => Dict("color" => color, "width" => rank == 1 ? 8 : 4),
        "opacity" => rank == 1 ? 1.0 : 0.72,
        "hovertemplate" => @sprintf(
            "rank %d<br>objective %.8g<br>x %%{x:.3f}<br>y %%{y:.3f}<br>z %%{z:.3f}<extra></extra>",
            rank, cost,
        ),
        "showlegend" => true,
    )
end

function _hypr_empty_top5_trace(rank::Int, scene, color)
    return Dict(
        "type" => "scatter3d", "mode" => "lines", "scene" => scene,
        "name" => "rank $rank · unavailable", "x" => Float64[],
        "y" => Float64[], "z" => Float64[],
        "line" => Dict("color" => color, "width" => rank == 1 ? 8 : 4),
        "visible" => false,
    )
end

function _hypr_reference_trace(path, name, color, dash, scene, sample_path,
                               spacing_m::Real; cost=nothing, width=6,
                               curve_type::Symbol=:bezier)
    if path === nothing
        return Dict(
            "type" => "scatter3d", "mode" => "lines", "scene" => scene,
            "name" => "$name · unavailable", "x" => Float64[],
            "y" => Float64[], "z" => Float64[], "visible" => false,
        )
    end
    sampled = _hypr_sample_top_path(
        path, sample_path, spacing_m; curve_type=curve_type,
    )
    label = cost === nothing ? name : @sprintf("%s · %.6g", name, cost)
    hover_cost = cost === nothing ? "" : @sprintf("<br>objective %.8g", cost)
    return Dict(
        "type" => "scatter3d", "mode" => "lines", "scene" => scene,
        "name" => label,
        "x" => vec(sampled[1, :]), "y" => vec(sampled[2, :]),
        "z" => vec(sampled[3, :]),
        "line" => Dict("color" => color, "width" => width, "dash" => dash),
        "opacity" => 0.95,
        "hovertemplate" => "$name$hover_cost<br>x %{x:.3f}<br>y %{y:.3f}<br>z %{z:.3f}<extra></extra>",
        "showlegend" => true,
    )
end

function _hypr_top5_frame_traces(frame, scene, sample_path, spacing_m, colors)
    traces = Any[]
    for rank in 1:5
        if rank <= length(frame.paths)
            push!(traces, _hypr_top5_trace(
                frame.paths[rank], frame.costs[rank], rank, scene,
                sample_path, spacing_m, colors[rank],
            ))
        else
            push!(traces, _hypr_empty_top5_trace(rank, scene, colors[rank]))
        end
    end
    return traces
end

function _hypr_endpoint_trace(point, name, color, symbol, scene)
    return Dict(
        "type" => "scatter3d", "mode" => "markers", "scene" => scene,
        "name" => name, "x" => [point[1]], "y" => [point[2]],
        "z" => [point[3]], "marker" => Dict(
            "color" => color, "size" => 7, "symbol" => symbol,
        ), "showlegend" => false,
    )
end

function _hypr_top5_title(case_index, iteration, original_frame, retimed_frame)
    original_cost = isempty(original_frame.costs) ? Inf : original_frame.costs[1]
    retimed_cost = isempty(retimed_frame.costs) ? Inf : retimed_frame.costs[1]
    return @sprintf(
        "Case %03d · PSO iteration %d · original best %.6g · retimed best %.6g",
        case_index, iteration, original_cost, retimed_cost,
    )
end

function _write_hypr_top5_html(path, case_index, scenario, original_history,
                               retimed_history, original_result, retimed_result,
                               station_asset::Symbol;
                               spacing_m::Real=0.08)
    modules = SpaceAGORA_RL._spaceagora_rpo_modules()
    sample_path = getproperty(modules.guidance, :rpo_sample_path)
    mesh = _station_mesh_trace(_load_station_mesh(station_asset), station_asset)
    mesh["opacity"] = 0.42
    original_mesh = copy(mesh)
    original_mesh["scene"] = "scene"
    original_mesh["showlegend"] = false
    retimed_mesh = copy(mesh)
    retimed_mesh["scene"] = "scene2"
    retimed_mesh["showlegend"] = false
    colors = [
        "rgb(25,92,190)", "rgb(238,126,34)", "rgb(38,150,88)",
        "rgb(146,78,180)", "rgb(210,55,70)",
    ]
    maximum_iteration = max(
        maximum(frame.iteration for frame in original_history),
        maximum(frame.iteration for frame in retimed_history),
    )
    iterations = collect(1:maximum_iteration)
    first_iteration = first(iterations)
    original_first = _hypr_history_frame(original_history, first_iteration)
    retimed_first = _hypr_history_frame(retimed_history, first_iteration)
    data = Any[
        original_mesh,
        _hypr_endpoint_trace(
            scenario.start_rtn, "start", "rgb(35,155,85)", "circle", "scene",
        ),
        _hypr_endpoint_trace(
            scenario.goal_rtn, "goal", "rgb(195,50,55)", "diamond", "scene",
        ),
    ]
    append!(data, _hypr_top5_frame_traces(
        original_first, "scene", sample_path, spacing_m, colors,
    ))
    push!(data, _hypr_reference_trace(
        original_result.warmstart_path, "RRT-Connect warm start",
        "rgb(55,55,55)", "dash", "scene", sample_path, spacing_m;
        width=5, curve_type=:polyline,
    ))
    push!(data, _hypr_reference_trace(
        original_result.initial_seed_path, "RRT → Bézier PSO seed",
        "rgb(245,190,35)", "dot", "scene", sample_path, spacing_m;
        width=7,
    ))
    push!(data, _hypr_reference_trace(
        original_result.path, "final post-refined", "rgb(0,180,215)",
        "solid", "scene", sample_path, spacing_m;
        cost=original_result.cost, width=9,
    ))
    append!(data, Any[
        retimed_mesh,
        _hypr_endpoint_trace(
            scenario.start_rtn, "start", "rgb(35,155,85)", "circle", "scene2",
        ),
        _hypr_endpoint_trace(
            scenario.goal_rtn, "goal", "rgb(195,50,55)", "diamond", "scene2",
        ),
    ])
    append!(data, _hypr_top5_frame_traces(
        retimed_first, "scene2", sample_path, spacing_m, colors,
    ))
    push!(data, _hypr_reference_trace(
        retimed_result.warmstart_path, "RRT-Connect warm start",
        "rgb(55,55,55)", "dash", "scene2", sample_path, spacing_m;
        width=5, curve_type=:polyline,
    ))
    push!(data, _hypr_reference_trace(
        retimed_result.initial_seed_path, "RRT → Bézier PSO seed",
        "rgb(245,190,35)", "dot", "scene2", sample_path, spacing_m;
        width=7,
    ))
    push!(data, _hypr_reference_trace(
        retimed_result.path, "final post-refined", "rgb(225,20,205)",
        "solid", "scene2", sample_path, spacing_m;
        cost=retimed_result.cost, width=9,
    ))

    animated_trace_indices = vcat(collect(3:7), collect(14:18))
    frames = Any[]
    for iteration in iterations
        original_frame = _hypr_history_frame(original_history, iteration)
        retimed_frame = _hypr_history_frame(retimed_history, iteration)
        frame_traces = _hypr_top5_frame_traces(
            original_frame, "scene", sample_path, spacing_m, colors,
        )
        append!(frame_traces, _hypr_top5_frame_traces(
            retimed_frame, "scene2", sample_path, spacing_m, colors,
        ))
        push!(frames, Dict(
            "name" => string(iteration),
            "data" => frame_traces,
            "traces" => animated_trace_indices,
            "layout" => Dict("title" => Dict("text" => _hypr_top5_title(
                case_index, iteration, original_frame, retimed_frame,
            ))),
        ))
    end

    scene_axis = Dict("title" => "m", "showbackground" => true)
    scene = Dict(
        "domain" => Dict("x" => [0.0, 0.485], "y" => [0.0, 1.0]),
        "aspectmode" => "data",
        "xaxis" => merge(copy(scene_axis), Dict("title" => "radial (m)")),
        "yaxis" => merge(copy(scene_axis), Dict("title" => "along-track (m)")),
        "zaxis" => merge(copy(scene_axis), Dict("title" => "cross-track (m)")),
    )
    scene2 = deepcopy(scene)
    scene2["domain"] = Dict("x" => [0.515, 1.0], "y" => [0.0, 1.0])
    slider_steps = [Dict(
        "label" => string(iteration), "method" => "animate",
        "args" => [[string(iteration)], Dict(
            "mode" => "immediate",
            "frame" => Dict("duration" => 0, "redraw" => true),
            "transition" => Dict("duration" => 0),
        )],
    ) for iteration in iterations]
    layout = Dict(
        "title" => Dict("text" => _hypr_top5_title(
            case_index, first_iteration, original_first, retimed_first,
        )),
        "margin" => Dict("l" => 5, "r" => 5, "b" => 70, "t" => 95),
        "scene" => scene, "scene2" => scene2,
        "legend" => Dict("orientation" => "h", "y" => -0.08),
        "annotations" => [
            Dict(
                "text" => "Original HYPR · proxy objective", "x" => 0.2425,
                "y" => 1.04, "xref" => "paper", "yref" => "paper",
                "showarrow" => false, "font" => Dict("size" => 15),
            ),
            Dict(
                "text" => "New HYPR · retimed fuel + optimized pointing",
                "x" => 0.7575, "y" => 1.04, "xref" => "paper",
                "yref" => "paper", "showarrow" => false,
                "font" => Dict("size" => 15),
            ),
        ],
        "sliders" => [Dict(
            "active" => 0, "currentvalue" => Dict("prefix" => "iteration "),
            "pad" => Dict("t" => 35), "steps" => slider_steps,
        )],
        "updatemenus" => [Dict(
            "type" => "buttons", "direction" => "left", "x" => 0.0,
            "y" => -0.01, "showactive" => false, "buttons" => [
                Dict("label" => "Play", "method" => "animate", "args" => [
                    nothing, Dict(
                        "fromcurrent" => true,
                        "frame" => Dict("duration" => 450, "redraw" => true),
                        "transition" => Dict("duration" => 0),
                    ),
                ]),
                Dict("label" => "Pause", "method" => "animate", "args" => [
                    [nothing], Dict(
                        "mode" => "immediate",
                        "frame" => Dict("duration" => 0, "redraw" => false),
                        "transition" => Dict("duration" => 0),
                    ),
                ]),
            ],
        )],
    )
    open(path, "w") do io
        print(io, """<!doctype html><html><head><meta charset=\"utf-8\">
<title>Case $(@sprintf("%03d", case_index)) PSO top-five history</title>
<script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script></head>
<body style=\"margin:0\"><div id=\"comparison\" style=\"width:100vw;height:100vh\"></div>
<script>
const data = $(JSON.json(data));
const layout = $(JSON.json(layout));
const frames = $(JSON.json(frames));
Plotly.newPlot('comparison', data, layout, {responsive:true}).then(function() {
  return Plotly.addFrames('comparison', frames);
});
</script></body></html>""")
    end
    return path
end

function _write_hypr_top5_costs(path, original_history, retimed_history)
    rows = NamedTuple[]
    for (method, history) in (
        ("original_proxy", original_history),
        ("retimed_pointing_optimized", retimed_history),
    )
        for frame in history, rank in eachindex(frame.costs)
            push!(rows, (
                method=method, iteration=frame.iteration, rank=rank,
                objective_cost=frame.costs[rank],
            ))
        end
    end
    CSV.write(path, DataFrame(rows))
    return path
end

function top5_main(args=ARGS)
    config_path = isempty(args) ?
        joinpath(@__DIR__, "..", "..", "configs", "rpo", "hypr_rl.toml") :
        args[1]
    output_directory = length(args) >= 2 ? args[2] : joinpath(
        @__DIR__, "..", "..", "outputs", "hypr_rl",
        "case_024_pso_top5_bezier_comparison",
    )
    case_index = length(args) >= 3 ? parse(Int, args[3]) : 24
    raw = TOML.parsefile(config_path)
    task = raw["task"]
    scenario_config = raw["scenario"]
    evaluation_config = raw["evaluation"]
    config = _hypr_baseline_task_config(task)
    evaluation_seed = Int(evaluation_config["seed"])
    scenario_seed = evaluation_seed + 2 * case_index
    planner_seed = evaluation_seed + 2 * case_index + 1
    station_asset = Symbol(scenario_config["station_asset"])
    station_points = 10_000
    sample_ds_m = 0.5
    modules = SpaceAGORA_RL._spaceagora_rpo_modules()
    pso_config = Base.invokelatest(
        getproperty(modules.guidance, :rpo_740_mpc_final_pso_config);
        safe_distance_m=config.safe_distance_m,
        n_particles=Int(evaluation_config["hypr_baseline_particles"]),
        n_iters=Int(evaluation_config["hypr_baseline_iterations"]),
        n_waypoints=Int(evaluation_config["hypr_baseline_waypoints"]),
        sample_ds_m=sample_ds_m,
        curve_type=:bezier,
        refinement_enable=true,
    )
    base_scenario = build_rpo_hypr_rl_scenario(
        station_asset=station_asset,
        station_points=station_points,
        station_seed=scenario_config["station_seed"],
        station_keepout_radius_m=scenario_config["station_keepout_radius_m"],
        pso_config=pso_config,
    )
    sampler = build_rpo_hypr_rl_endpoint_sampler(
        base_scenario;
        station_asset=station_asset,
        safe_distance_m=config.safe_distance_m,
        endpoint_clearance_margin_m=scenario_config["endpoint_clearance_margin_m"],
        endpoint_max_clearance_m=scenario_config["endpoint_max_clearance_m"],
        min_separation_m=scenario_config["min_separation_m"],
        surrounded_max_distance_m=scenario_config["surrounded_max_distance_m"],
        max_sampling_tries=scenario_config["max_sampling_tries"],
    )
    scenario = sample_rpo_hypr_rl_scenario(
        sampler, MersenneTwister(scenario_seed),
    )
    mkpath(output_directory)
    println("case=$case_index scenario_seed=$scenario_seed planner_seed=$planner_seed")
    println("capturing original proxy PSO history")
    original_history, original_result, original_runtime_s =
        _run_hypr_top5_history(
            scenario, config, pso_config, planner_seed; retimed=false,
        )
    println("capturing retimed pointing-optimized PSO history")
    retimed_history, retimed_result, retimed_runtime_s =
        _run_hypr_top5_history(
            scenario, config, pso_config, planner_seed; retimed=true,
        )
    html_path = _write_hypr_top5_html(
        joinpath(output_directory, @sprintf(
            "case_%03d_pso_top5_iterations.html", case_index,
        )),
        case_index, scenario, original_history, retimed_history,
        original_result, retimed_result, station_asset,
    )
    costs_path = _write_hypr_top5_costs(
        joinpath(output_directory, "top5_costs.csv"),
        original_history, retimed_history,
    )
    warmstart_paths_match = original_result.warmstart_path !== nothing &&
        retimed_result.warmstart_path !== nothing &&
        original_result.warmstart_path == retimed_result.warmstart_path
    original_warmstart_polyline_clearance_m = _hypr_reference_minimum_clearance(
        original_result.warmstart_path, scenario; curve_type=:polyline,
    )
    original_warmstart_bezier_render_clearance_m = _hypr_reference_minimum_clearance(
        original_result.warmstart_path, scenario; curve_type=:bezier,
    )
    retimed_warmstart_polyline_clearance_m = _hypr_reference_minimum_clearance(
        retimed_result.warmstart_path, scenario; curve_type=:polyline,
    )
    retimed_warmstart_bezier_render_clearance_m = _hypr_reference_minimum_clearance(
        retimed_result.warmstart_path, scenario; curve_type=:bezier,
    )
    original_initial_bezier_seed_clearance_m = _hypr_reference_minimum_clearance(
        original_result.initial_seed_path, scenario; curve_type=:bezier,
    )
    retimed_initial_bezier_seed_clearance_m = _hypr_reference_minimum_clearance(
        retimed_result.initial_seed_path, scenario; curve_type=:bezier,
    )
    manifest = Dict(
        "case" => case_index,
        "config" => abspath(config_path),
        "scenario_seed" => scenario_seed,
        "planner_seed" => planner_seed,
        "station_asset" => String(station_asset),
        "station_points" => station_points,
        "pso_particles" => pso_config.n_particles,
        "pso_iterations_requested" => pso_config.n_iters,
        "pso_waypoints_requested" => pso_config.n_waypoints,
        "pso_sample_ds_m" => sample_ds_m,
        "pso_curve_type" => "bezier",
        "archive" => "global best plus four lowest distinct personal-best trajectories after each recorded iteration",
        "reference_trajectories" => "RRT-Connect warm start, its fitted initial Bezier PSO seed, and the final post-refined path are overlaid in each panel",
        "warmstart_display_curve_type" => "polyline",
        "warmstart_paths_match" => warmstart_paths_match,
        "original_warmstart_polyline_min_clearance_m" => original_warmstart_polyline_clearance_m,
        "original_warmstart_incorrect_bezier_render_min_clearance_m" => original_warmstart_bezier_render_clearance_m,
        "retimed_warmstart_polyline_min_clearance_m" => retimed_warmstart_polyline_clearance_m,
        "retimed_warmstart_incorrect_bezier_render_min_clearance_m" => retimed_warmstart_bezier_render_clearance_m,
        "original_initial_bezier_seed_min_clearance_m" => original_initial_bezier_seed_clearance_m,
        "retimed_initial_bezier_seed_min_clearance_m" => retimed_initial_bezier_seed_clearance_m,
        "original_iterations_recorded" => length(original_history),
        "retimed_iterations_recorded" => length(retimed_history),
        "original_runtime_s" => original_runtime_s,
        "retimed_runtime_s" => retimed_runtime_s,
        "original_final_pso_objective" => original_history[end].costs[1],
        "retimed_final_pso_objective" => retimed_history[end].costs[1],
        "original_post_refinement_objective" => original_result.cost,
        "retimed_post_refinement_objective" => retimed_result.cost,
        "html" => abspath(html_path),
        "costs_csv" => abspath(costs_path),
    )
    open(joinpath(output_directory, "run_manifest.toml"), "w") do io
        TOML.print(io, manifest; sorted=true)
    end
    println("wrote ", abspath(html_path))
    return (
        html_path=html_path, costs_path=costs_path,
        original_history=original_history, retimed_history=retimed_history,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    top5_main()
end
