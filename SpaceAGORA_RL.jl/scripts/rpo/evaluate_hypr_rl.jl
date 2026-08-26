#!/usr/bin/env julia

using CSV
using DataFrames
using Distributed
using JSON
using LinearAlgebra
using Printf
using Random
using Statistics
using TOML
using SpaceAGORA_RL

_ntget(value, name::Symbol, default) =
    value !== nothing && hasproperty(value, name) ? getproperty(value, name) : default

function _latest_hypr_rl_checkpoint(directory::AbstractString)
    isdir(directory) || throw(ArgumentError("checkpoint directory does not exist: $directory"))
    candidates = filter(readdir(directory; join=true)) do path
        basename(path) == "hypr_rl_final.jls" ||
            occursin(r"^hypr_rl_episode_[0-9]+\.jls$", basename(path))
    end
    isempty(candidates) && throw(ArgumentError("no HyPR-RL checkpoints found in $directory"))
    return argmax(path -> stat(path).mtime, candidates)
end

function _resolve_checkpoint(requested::AbstractString, directory::AbstractString)
    path = lowercase(requested) == "latest" ?
        _latest_hypr_rl_checkpoint(directory) : requested
    isfile(path) || throw(ArgumentError("checkpoint does not exist: $path"))
    return abspath(path)
end

function _validate_checkpoint(path::AbstractString)
    payload = load_checkpoint(path)
    get(payload, :task, nothing) == :rpo_hypr_rl ||
        throw(ArgumentError("checkpoint is not a HyPR-RL checkpoint: $path"))
    network = payload[:online]
    all(field -> all(isfinite, getfield(network, field)), fieldnames(typeof(network))) ||
        throw(ArgumentError("checkpoint policy contains non-finite weights: $path"))
    metadata = get(payload, :task_metadata, nothing)
    config = _ntget(metadata, :config, nothing)
    config isa RPOHyPRRLConfig ||
        throw(ArgumentError("checkpoint does not contain its RPOHyPRRLConfig: $path"))
    return payload, config
end

function _path_length(path::AbstractMatrix)
    size(path, 2) < 2 && return 0.0
    return sum(norm(view(path, :, index + 1) - view(path, :, index))
               for index in 1:(size(path, 2) - 1))
end

function _tracked_path(plan)
    plan === nothing && return zeros(3, 0)
    evaluator = _ntget(plan.diagnostics, :evaluator, nothing)
    history = _ntget(evaluator, :state_history, zeros(6, 0))
    return size(history, 1) >= 3 ? Matrix{Float64}(history[1:3, :]) : zeros(3, 0)
end

function _clearance_stats(path::AbstractMatrix, scenario, safe_distance_m::Real)
    isempty(path) && return (minimum_m=NaN, violations=-1, clearances=Float64[])
    modules = SpaceAGORA_RL._spaceagora_rpo_modules()
    clearance = getproperty(modules.navigation, :rpo_clearance_distance_to_station)
    values = [Float64(Base.invokelatest(
        clearance, view(path, :, index), scenario.geometry,
    )) for index in axes(path, 2)]
    return (
        minimum_m=minimum(values),
        violations=count(value -> value + 1.0e-9 < safe_distance_m, values),
        clearances=values,
    )
end

function _case_metrics(case_index::Int, scenario, result, config)
    plan = _ntget(result, :plan, nothing)
    runtime_s = Float64(_ntget(result, :runtime_s, NaN))
    error_message = String(_ntget(result, :error, ""))
    if plan === nothing
        return (
            case=case_index, success=false, stopped=false, edits=0,
            runtime_s=runtime_s, minimum_clearance_m=NaN,
            reference_minimum_clearance_m=NaN, tracked_minimum_clearance_m=NaN,
            safety_violation_count=-1, reference_safety_violation_count=-1,
            propellant_used_kg=NaN, propellant_used_g=NaN, duration_s=NaN,
            reference_path_length_m=NaN, tracked_path_length_m=NaN,
            final_position_error_m=NaN, objective=NaN, edit_objective=NaN,
            seed_objective=NaN, allocation_error_impulse_mps=NaN,
            thruster_saturation_fraction=NaN, wheel_energy_j=NaN,
            wheel_peak_momentum_nms=NaN,
            thruster_1_impulse_ns=NaN, thruster_2_impulse_ns=NaN,
            thruster_3_impulse_ns=NaN, thruster_4_impulse_ns=NaN,
            thruster_5_impulse_ns=NaN, thruster_6_impulse_ns=NaN,
            start_r_m=scenario.start_rtn[1], start_t_m=scenario.start_rtn[2],
            start_n_m=scenario.start_rtn[3], goal_r_m=scenario.goal_rtn[1],
            goal_t_m=scenario.goal_rtn[2], goal_n_m=scenario.goal_rtn[3],
            evaluation_status="worker_error", error_message=error_message,
        )
    end

    reference_stats = _clearance_stats(
        plan.r_ref_rtn, scenario, config.safe_distance_m,
    )
    tracked = _tracked_path(plan)
    tracked_stats = _clearance_stats(tracked, scenario, config.safe_distance_m)
    evaluator = _ntget(plan.diagnostics, :evaluator, nothing)
    isempty(error_message) &&
        (error_message = String(_ntget(evaluator, :error_message, "")))
    impulses = Vector{Float64}(_ntget(
        plan.diagnostics, :thruster_impulse_ns, zeros(6),
    ))
    length(impulses) == 6 || (impulses = fill(NaN, 6))
    terminal_minimum = Float64(_ntget(plan.diagnostics, :min_clearance_m, NaN))
    minimum_clearance = isfinite(tracked_stats.minimum_m) ?
        tracked_stats.minimum_m : terminal_minimum
    violation_count = Int(_ntget(
        evaluator, :keepout_violations, tracked_stats.violations,
    ))
    status = hasproperty(evaluator, :evaluator_mode) ?
        String(evaluator.evaluator_mode) : String(_ntget(evaluator, :reason, :unknown))
    propellant = plan.propellant_used_kg
    return (
        case=case_index, success=plan.valid, stopped=plan.stopped,
        edits=plan.edits, runtime_s=runtime_s,
        minimum_clearance_m=minimum_clearance,
        reference_minimum_clearance_m=reference_stats.minimum_m,
        tracked_minimum_clearance_m=tracked_stats.minimum_m,
        safety_violation_count=violation_count,
        reference_safety_violation_count=reference_stats.violations,
        propellant_used_kg=propellant, propellant_used_g=1_000.0 * propellant,
        duration_s=Float64(_ntget(plan.diagnostics, :duration_s, NaN)),
        reference_path_length_m=_path_length(plan.r_ref_rtn),
        tracked_path_length_m=_path_length(tracked),
        final_position_error_m=Float64(_ntget(
            plan.diagnostics, :final_position_error_m, NaN,
        )),
        objective=plan.cost,
        edit_objective=Float64(_ntget(plan.diagnostics, :edit_cost, NaN)),
        seed_objective=plan.seed_cost,
        allocation_error_impulse_mps=Float64(_ntget(
            plan.diagnostics, :allocation_error_impulse_mps, NaN,
        )),
        thruster_saturation_fraction=Float64(_ntget(
            plan.diagnostics, :thruster_saturation_fraction, NaN,
        )),
        wheel_energy_j=Float64(_ntget(plan.diagnostics, :wheel_energy_j, NaN)),
        wheel_peak_momentum_nms=Float64(_ntget(
            plan.diagnostics, :wheel_peak_momentum_nms, NaN,
        )),
        thruster_1_impulse_ns=impulses[1], thruster_2_impulse_ns=impulses[2],
        thruster_3_impulse_ns=impulses[3], thruster_4_impulse_ns=impulses[4],
        thruster_5_impulse_ns=impulses[5], thruster_6_impulse_ns=impulses[6],
        start_r_m=scenario.start_rtn[1], start_t_m=scenario.start_rtn[2],
        start_n_m=scenario.start_rtn[3], goal_r_m=scenario.goal_rtn[1],
        goal_t_m=scenario.goal_rtn[2], goal_n_m=scenario.goal_rtn[3],
        evaluation_status=status, error_message=error_message,
    )
end

function _finite_mean(values)
    finite = filter(isfinite, Float64.(collect(skipmissing(values))))
    return isempty(finite) ? NaN : mean(finite)
end

function _finite_minimum(values)
    finite = filter(isfinite, Float64.(collect(skipmissing(values))))
    return isempty(finite) ? NaN : minimum(finite)
end

function _evaluation_summary(rows::DataFrame, checkpoint, payload)
    violation_counts = filter(>=(0), rows.safety_violation_count)
    successful = rows.success
    metadata = get(payload, :task_metadata, nothing)
    training = _ntget(metadata, :training, nothing)
    checkpoint_episode = Int(_ntget(
        training, :episodes, _ntget(training, :episode, 0),
    ))
    successful_rows = rows[rows.success, :]
    return (
        cases=nrow(rows), success_count=count(successful),
        failure_count=count(!, successful),
        success_rate=count(successful) / max(nrow(rows), 1),
        mean_minimum_clearance_m=_finite_mean(rows.minimum_clearance_m),
        worst_minimum_clearance_m=_finite_minimum(rows.minimum_clearance_m),
        mean_successful_minimum_clearance_m=
            _finite_mean(successful_rows.minimum_clearance_m),
        total_safety_violations=sum(violation_counts; init=0),
        mean_safety_violations=isempty(violation_counts) ? NaN : mean(violation_counts),
        cases_with_safety_violations=count(>(0), violation_counts),
        cases_without_full_trajectory_metrics=count(<(0), rows.safety_violation_count),
        mean_propellant_used_kg=_finite_mean(rows.propellant_used_kg),
        mean_propellant_used_g=_finite_mean(rows.propellant_used_g),
        mean_successful_propellant_used_g=
            _finite_mean(successful_rows.propellant_used_g),
        mean_duration_s=_finite_mean(rows.duration_s),
        mean_reference_path_length_m=_finite_mean(rows.reference_path_length_m),
        mean_tracked_path_length_m=_finite_mean(rows.tracked_path_length_m),
        mean_final_position_error_m=_finite_mean(rows.final_position_error_m),
        mean_objective=_finite_mean(rows.objective),
        mean_successful_final_position_error_m=
            _finite_mean(successful_rows.final_position_error_m),
        mean_successful_objective=_finite_mean(successful_rows.objective),
        mean_edits=mean(rows.edits), mean_runtime_s=_finite_mean(rows.runtime_s),
        mean_allocation_error_impulse_mps=_finite_mean(rows.allocation_error_impulse_mps),
        mean_thruster_saturation_fraction=_finite_mean(rows.thruster_saturation_fraction),
        mean_wheel_energy_j=_finite_mean(rows.wheel_energy_j),
        mean_wheel_peak_momentum_nms=_finite_mean(rows.wheel_peak_momentum_nms),
        checkpoint=checkpoint === nothing ? "" : abspath(checkpoint),
        checkpoint_episode=checkpoint_episode,
        checkpoint_global_step=Int(get(payload, :global_step, 0)),
        checkpoint_train_steps=Int(get(payload, :train_steps, 0)),
        checkpoint_last_loss=Float64(get(payload, :last_loss, NaN)),
        checkpoint_mean_training_loss=Float64(get(payload, :mean_training_loss, NaN)),
    )
end

function _downsample(path::AbstractMatrix, maximum_points::Int)
    size(path, 2) <= maximum_points && return path
    indices = unique(round.(Int, range(1, size(path, 2); length=maximum_points)))
    return path[:, indices]
end

function _line_trace(path, name, color; dash="solid", markers=false, width=5)
    sampled = _downsample(path, 2_000)
    return Dict(
        "type" => "scatter3d", "name" => name,
        "x" => vec(sampled[1, :]), "y" => vec(sampled[2, :]),
        "z" => vec(sampled[3, :]),
        "mode" => markers ? "lines+markers" : "lines",
        "line" => Dict("color" => color, "width" => width, "dash" => dash),
        "marker" => Dict("color" => color, "size" => 3),
    )
end

function _marker_trace(path, name, color; size=5, visible=true)
    return Dict(
        "type" => "scatter3d", "name" => name,
        "x" => vec(path[1, :]), "y" => vec(path[2, :]),
        "z" => vec(path[3, :]), "mode" => "markers",
        "marker" => Dict("color" => color, "size" => size),
        "visible" => visible,
    )
end

function _sample_planned_bezier(control_points; spacing_m::Real=0.05)
    modules = SpaceAGORA_RL._spaceagora_rpo_modules()
    sampler = Base.invokelatest(
        getproperty, modules.guidance, :rpo_sample_path_bezier,
    )
    return Matrix{Float64}(Base.invokelatest(
        sampler, control_points, Float64(spacing_m),
    ))
end

function _load_station_mesh(station_asset::Symbol)
    spaceagora = SpaceAGORA_RL._load_spaceagora!(load_gramsuite=false)
    loader = Base.invokelatest(
        getproperty, spaceagora, :load_rpo_station_cad_triangles,
    )
    return Matrix{Float64}(Base.invokelatest(loader, station_asset))
end

function _station_mesh_trace(triangles::AbstractMatrix, station_asset::Symbol)
    size(triangles, 1) == 3 ||
        throw(ArgumentError("station mesh triangles must have three rows"))
    size(triangles, 2) % 3 == 0 ||
        throw(ArgumentError("station mesh must contain complete triangles"))

    vertices = NTuple{3, Float64}[]
    vertex_indices = Dict{NTuple{3, Float64}, Int}()
    triangle_indices = Vector{Int}(undef, size(triangles, 2))
    for column in axes(triangles, 2)
        vertex = (
            Float64(triangles[1, column]), Float64(triangles[2, column]),
            Float64(triangles[3, column]),
        )
        triangle_indices[column] = get!(vertex_indices, vertex) do
            push!(vertices, vertex)
            length(vertices) - 1 # Plotly mesh indices are zero based.
        end
    end

    return Dict(
        "type" => "mesh3d", "name" => "$(station_asset) station mesh",
        "x" => [vertex[1] for vertex in vertices],
        "y" => [vertex[2] for vertex in vertices],
        "z" => [vertex[3] for vertex in vertices],
        "i" => triangle_indices[1:3:end],
        "j" => triangle_indices[2:3:end],
        "k" => triangle_indices[3:3:end],
        "color" => "rgb(145,155,165)", "opacity" => 0.70,
        "flatshading" => true, "hoverinfo" => "skip",
        "lighting" => Dict(
            "ambient" => 0.55, "diffuse" => 0.75, "specular" => 0.20,
            "roughness" => 0.75, "fresnel" => 0.10,
        ),
        "lightposition" => Dict("x" => 100, "y" => 200, "z" => 300),
    )
end

function _write_station_mesh_asset(path, triangles, station_asset::Symbol)
    trace = _station_mesh_trace(triangles, station_asset)
    open(path, "w") do io
        print(io, "window.HYPR_RL_STATION_MESH = ", JSON.json(trace), ";\n")
    end
    return path
end

function _write_case_plot(path, case_index, scenario, result, metrics, config;
                          station_mesh_script::Union{Nothing, String}=nothing,
                          station_max_points::Int=2_000,
                          bezier_plot_spacing_m::Real=0.05,
                          planner_label::AbstractString="HyPR-RL")
    plan = _ntget(result, :plan, nothing)
    traces = Any[]
    if station_mesh_script === nothing
        station = getproperty(getproperty(scenario.geometry, :station), :points_body)
        station = _downsample(station, station_max_points)
        push!(traces, Dict(
            "type" => "scatter3d", "name" => "station point cloud (mesh unavailable)",
            "x" => vec(station[1, :]), "y" => vec(station[2, :]),
            "z" => vec(station[3, :]), "mode" => "markers",
            "marker" => Dict("color" => "rgb(145,155,165)", "size" => 2,
                             "opacity" => 0.30), "hoverinfo" => "skip",
        ))
    end
    if plan !== nothing
        if size(plan.path_rtn, 2) > 1
            planned_bezier = _sample_planned_bezier(
                plan.path_rtn; spacing_m=bezier_plot_spacing_m,
            )
            push!(traces, _line_trace(
                planned_bezier, "planned Bezier trajectory", "rgb(35,90,205)";
                width=7,
            ))
            push!(traces, _marker_trace(
                plan.path_rtn, "Bezier control waypoints", "rgb(70,75,85)";
                size=5,
            ))
            sampled = _downsample(planned_bezier, 2_000)
            unsafe = _clearance_stats(
                sampled, scenario, config.safe_distance_m,
            ).clearances .+ 1.0e-9 .< config.safe_distance_m
            if any(unsafe)
                push!(traces, Dict(
                    "type" => "scatter3d", "name" => "planned-path violations",
                    "x" => vec(sampled[1, unsafe]), "y" => vec(sampled[2, unsafe]),
                    "z" => vec(sampled[3, unsafe]), "mode" => "markers",
                    "marker" => Dict("color" => "rgb(215,45,45)", "size" => 5),
                ))
            end
        end
        if size(plan.r_ref_rtn, 2) > 1
            push!(traces, _marker_trace(
                _downsample(plan.r_ref_rtn, 500),
                "retimed reference samples", "rgb(80,160,235)";
                size=2, visible="legendonly",
            ))
        end
        tracked = _tracked_path(plan)
        if size(tracked, 2) > 1
            push!(traces, _line_trace(
                tracked, "executed LQ-MPC trajectory", "rgb(235,135,35)";
                width=6,
            ))
            sampled = _downsample(tracked, 2_000)
            unsafe = _clearance_stats(
                sampled, scenario, config.safe_distance_m,
            ).clearances .+ 1.0e-9 .< config.safe_distance_m
            if any(unsafe)
                push!(traces, Dict(
                    "type" => "scatter3d", "name" => "executed-path violations",
                    "x" => vec(sampled[1, unsafe]), "y" => vec(sampled[2, unsafe]),
                    "z" => vec(sampled[3, unsafe]), "mode" => "markers",
                    "marker" => Dict("color" => "rgb(215,45,45)", "size" => 5),
                ))
            end
        end
    end
    for (name, point, color, symbol) in (
        ("start", scenario.start_rtn, "rgb(35,155,85)", "circle"),
        ("goal", scenario.goal_rtn, "rgb(195,50,55)", "diamond"),
    )
        push!(traces, Dict(
            "type" => "scatter3d", "name" => name, "x" => [point[1]],
            "y" => [point[2]], "z" => [point[3]], "mode" => "markers",
            "marker" => Dict("color" => color, "size" => 7, "symbol" => symbol),
        ))
    end
    title = @sprintf(
        "%s case %03d | success=%s | min clearance=%.3f m | violations=%d | fuel=%.4f g",
        planner_label, case_index, metrics.success, metrics.minimum_clearance_m,
        metrics.safety_violation_count, metrics.propellant_used_g,
    )
    layout = Dict(
        "title" => title, "margin" => Dict("l" => 0, "r" => 0, "b" => 0, "t" => 55),
        "legend" => Dict("orientation" => "h"),
        "scene" => Dict(
            "aspectmode" => "data", "xaxis" => Dict("title" => "radial (m)"),
            "yaxis" => Dict("title" => "along-track (m)"),
            "zaxis" => Dict("title" => "cross-track (m)"),
        ),
    )
    mesh_script_tag = station_mesh_script === nothing ? "" :
        "<script src=\"$(station_mesh_script)\"></script>"
    plot_data = station_mesh_script === nothing ? JSON.json(traces) :
        "[window.HYPR_RL_STATION_MESH].concat($(JSON.json(traces)))"
    open(path, "w") do io
        print(io, """<!doctype html><html><head><meta charset=\"utf-8\"><title>$(title)</title>
<script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>$(mesh_script_tag)</head>
<body style=\"margin:0\"><div id=\"trajectory\" style=\"width:100vw;height:100vh\"></div>
<script>Plotly.newPlot('trajectory', $(plot_data), $(JSON.json(layout)), {responsive:true});</script>
</body></html>""")
    end
    return path
end

_html_escape(value) = replace(string(value), '&' => "&amp;", '<' => "&lt;",
                              '>' => "&gt;", '"' => "&quot;")

function _write_index(path, rows::DataFrame, summary;
                      planner_label::AbstractString="HyPR-RL")
    open(path, "w") do io
        print(io, """<!doctype html><html><head><meta charset=\"utf-8\"><title>$(planner_label) evaluation</title>
<style>body{font:14px system-ui;margin:2rem;color:#18212b}table{border-collapse:collapse;width:100%}
th,td{padding:.45rem;border-bottom:1px solid #d9dfe5;text-align:right}th{background:#f2f5f7;position:sticky;top:0}
td:first-child,th:first-child{text-align:left}.ok{color:#16834a}.bad{color:#c0392b}.summary{columns:2;margin-bottom:2rem}</style>
</head><body><h1>$(planner_label) $(nrow(rows))-case evaluation</h1><div class=\"summary\">""")
        for name in propertynames(summary)
            print(io, "<div><b>$(_html_escape(name)):</b> $(_html_escape(getproperty(summary, name)))</div>")
        end
        print(io, """</div><p><a href=\"cases.csv\">Per-case CSV</a> · <a href=\"summary.csv\">Summary CSV</a></p>
<table><thead><tr><th>case</th><th>success</th><th>min clearance (m)</th><th>violations</th>
<th>fuel (g)</th><th>duration (s)</th><th>final error (m)</th><th>runtime (s)</th></tr></thead><tbody>""")
        for row in eachrow(rows)
            class_name = row.success ? "ok" : "bad"
            print(io, @sprintf(
                "<tr><td><a href=\"trajectories/case_%03d.html\">case %03d</a></td><td class=\"%s\">%s</td><td>%.4f</td><td>%d</td><td>%.5f</td><td>%.3f</td><td>%.4f</td><td>%.3f</td></tr>",
                row.case, row.case, class_name, row.success,
                row.minimum_clearance_m, row.safety_violation_count,
                row.propellant_used_g, row.duration_s,
                row.final_position_error_m, row.runtime_s,
            ))
        end
        print(io, "</tbody></table></body></html>")
    end
    return path
end

function _run_cases(config, scenarios, policy, policy_seeds, n_workers::Int)
    n_cases = length(scenarios)
    results = Vector{Any}(undef, n_cases)
    active_workers = min(max(n_workers, 1), n_cases)
    println("evaluating HyPR-RL cases=$n_cases active_workers=$active_workers terminal=full_lqmpc")
    if active_workers == 1
        for index in eachindex(scenarios)
            results[index] = try
                evaluate_hypr_rl_policy_case(
                    config, scenarios[index], policy, policy_seeds[index],
                )
            catch error
                (plan=nothing, runtime_s=NaN, error=sprint(showerror, error))
            end
            println("evaluation progress $index/$n_cases")
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
                evaluate_hypr_rl_policy_case,
                process_ids[mod1(offset, active_workers)],
                config, scenarios[index], policy, policy_seeds[index],
            ) for (offset, index) in enumerate(indices)]
            for (index, future) in zip(indices, futures)
                results[index] = try
                    fetch(future)
                catch error
                    (plan=nothing, runtime_s=NaN, error=sprint(showerror, error))
                end
                println("evaluation progress $index/$n_cases")
            end
        end
    finally
        rmprocs(process_ids)
    end
    return results
end

function main(args=ARGS)
    config_path = length(args) >= 2 ? args[2] :
        joinpath(@__DIR__, "..", "..", "configs", "rpo", "hypr_rl.toml")
    raw = TOML.parsefile(config_path)
    scenario_config = raw["scenario"]
    training_config = raw["training"]
    evaluation_config = get(raw, "evaluation", Dict{String, Any}())
    checkpoint = _resolve_checkpoint(
        isempty(args) ? "latest" : args[1], training_config["checkpoint_directory"],
    )
    payload, config = _validate_checkpoint(checkpoint)
    policy = load_hypr_rl_policy(checkpoint)
    n_cases = length(args) >= 4 ? parse(Int, args[4]) :
        Int(get(evaluation_config, "cases", 100))
    n_cases > 0 || throw(ArgumentError("number of evaluation cases must be positive"))
    output_directory = length(args) >= 3 ? args[3] : String(get(
        evaluation_config, "output_directory",
        "SpaceAGORA_RL.jl/outputs/hypr_rl/evaluation_100_cases",
    ))
    n_workers = length(args) >= 5 ? parse(Int, args[5]) :
        Int(get(evaluation_config, "n_workers", training_config["n_workers"]))
    evaluation_seed = Int(get(
        evaluation_config, "seed", training_config["seed"] + 10_000_000,
    ))
    station_max_points = Int(get(evaluation_config, "station_max_plot_points", 2_000))
    bezier_plot_spacing_m = Float64(get(
        evaluation_config, "bezier_plot_spacing_m", 0.05,
    ))
    bezier_plot_spacing_m > 0.0 ||
        throw(ArgumentError("evaluation Bezier plot spacing must be positive"))
    station_asset = Symbol(scenario_config["station_asset"])

    base_scenario = build_rpo_hypr_rl_scenario(
        station_asset=station_asset,
        station_points=scenario_config["station_points"],
        station_seed=scenario_config["station_seed"],
        station_keepout_radius_m=scenario_config["station_keepout_radius_m"],
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
    policy_seeds = [evaluation_seed + 2 * index + 1 for index in 1:n_cases]
    scenarios = [sample_rpo_hypr_rl_scenario(
        scenario_sampler, MersenneTwister(scenario_seeds[index]),
    ) for index in 1:n_cases]

    mkpath(output_directory)
    trajectory_directory = joinpath(output_directory, "trajectories")
    mkpath(trajectory_directory)
    station_mesh = try
        _load_station_mesh(station_asset)
    catch error
        @warn "station mesh unavailable; trajectory plots will use the point-cloud fallback" station_asset exception=(error, catch_backtrace())
        nothing
    end
    station_mesh_script = if station_mesh === nothing
        nothing
    else
        mesh_path = _write_station_mesh_asset(
            joinpath(output_directory, "station_mesh.js"), station_mesh, station_asset,
        )
        println(
            "wrote shared station mesh triangles=", size(station_mesh, 2) ÷ 3,
            " path=", abspath(mesh_path),
        )
        "../station_mesh.js"
    end
    results = _run_cases(config, scenarios, policy, policy_seeds, n_workers)
    metrics = [_case_metrics(index, scenarios[index], results[index], config)
               for index in 1:n_cases]
    rows = DataFrame(metrics)
    summary = _evaluation_summary(rows, checkpoint, payload)
    CSV.write(joinpath(output_directory, "cases.csv"), rows)
    CSV.write(joinpath(output_directory, "summary.csv"), DataFrame([summary]))
    for index in 1:n_cases
        _write_case_plot(
            joinpath(trajectory_directory, @sprintf("case_%03d.html", index)),
            index, scenarios[index], results[index], metrics[index], config;
            station_mesh_script=station_mesh_script,
            station_max_points=station_max_points,
            bezier_plot_spacing_m=bezier_plot_spacing_m,
        )
    end
    index_path = _write_index(joinpath(output_directory, "index.html"), rows, summary)
    manifest = Dict(
        "checkpoint" => checkpoint, "config" => abspath(config_path),
        "cases" => n_cases, "evaluation_seed" => evaluation_seed,
        "scenario_seeds" => scenario_seeds, "policy_seeds" => policy_seeds,
        "workers" => min(max(n_workers, 1), n_cases),
        "safe_distance_m" => config.safe_distance_m,
        "endpoint_distribution" => "canonical_hypr_surface_shell",
        "endpoint_min_clearance_m" =>
            scenario_sampler.endpoint_min_clearance_m,
        "endpoint_max_clearance_m" =>
            scenario_sampler.endpoint_max_clearance_m,
        "endpoint_min_separation_m" => scenario_sampler.min_separation_m,
        "endpoint_surrounded_max_distance_m" =>
            scenario_sampler.surrounded_max_distance_m,
        "edit_evaluator" => "retimed_feedforward",
        "terminal_evaluator" => "full_lqmpc",
        "station_visualization" => station_mesh === nothing ? "point_cloud" : "cad_mesh",
        "station_mesh_triangles" => station_mesh === nothing ? 0 : size(station_mesh, 2) ÷ 3,
        "index_html" => abspath(index_path),
    )
    open(joinpath(output_directory, "evaluation_manifest.toml"), "w") do io
        TOML.print(io, manifest; sorted=true)
    end
    println("wrote HyPR-RL evaluation to ", abspath(output_directory))
    println("success_rate=", summary.success_rate,
            " mean_minimum_clearance_m=", summary.mean_minimum_clearance_m,
            " total_safety_violations=", summary.total_safety_violations,
            " mean_propellant_used_g=", summary.mean_propellant_used_g)
    println("open ", abspath(index_path))
    return (rows=rows, summary=summary, results=results)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
