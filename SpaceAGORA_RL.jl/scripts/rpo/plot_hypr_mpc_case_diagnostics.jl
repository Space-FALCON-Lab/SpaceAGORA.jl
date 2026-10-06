#!/usr/bin/env julia

include(joinpath(@__DIR__, "compare_hypr_objectives.jl"))

function _hypr_mpc_diagnostic_config(raw)
    task = raw["task"]
    scenario = raw["scenario"]
    evaluation = raw["evaluation"]
    config = _hypr_baseline_task_config(task)
    modules = SpaceAGORA_RL._spaceagora_rpo_modules()
    pso_config = Base.invokelatest(
        getproperty(modules.guidance, :rpo_740_mpc_final_pso_config);
        safe_distance_m=config.safe_distance_m,
        n_particles=Int(evaluation["hypr_baseline_particles"]),
        n_iters=Int(evaluation["hypr_baseline_iterations"]),
        n_waypoints=Int(evaluation["hypr_baseline_waypoints"]),
        refinement_enable=true,
    )
    base_scenario = build_rpo_hypr_rl_scenario(
        station_asset=Symbol(scenario["station_asset"]),
        station_points=Int(scenario["station_points"]),
        station_seed=Int(scenario["station_seed"]),
        station_keepout_radius_m=Float64(scenario["station_keepout_radius_m"]),
        pso_config=pso_config,
    )
    sampler = build_rpo_hypr_rl_endpoint_sampler(
        base_scenario;
        station_asset=Symbol(scenario["station_asset"]),
        safe_distance_m=config.safe_distance_m,
        endpoint_clearance_margin_m=Float64(scenario["endpoint_clearance_margin_m"]),
        endpoint_max_clearance_m=Float64(scenario["endpoint_max_clearance_m"]),
        min_separation_m=Float64(scenario["min_separation_m"]),
        surrounded_max_distance_m=Float64(scenario["surrounded_max_distance_m"]),
        max_sampling_tries=Int(scenario["max_sampling_tries"]),
    )
    return config, sampler, Int(evaluation["seed"])
end

function _hypr_mpc_tracking_error(plan, dt_s::Real)
    evaluator = plan.diagnostics.evaluator
    state_history = evaluator.state_history
    n_states = size(state_history, 2)
    n_reference = size(plan.r_ref_rtn, 2)
    time_s = collect(0:(n_states - 1)) .* Float64(dt_s)
    error_m = [norm(
        view(state_history, 1:3, index) -
        view(plan.r_ref_rtn, :, min(index, n_reference)),
    ) for index in 1:n_states]
    return time_s, error_m
end

function _hypr_mpc_thruster_forces(plan, config, dt_s::Real)
    evaluator = plan.diagnostics.evaluator
    commands = evaluator.command_history
    attitudes = evaluator.actual_attitude_history
    n_steps = size(commands, 2)
    dt = Float64(dt_s)
    allocator = SpaceAGORA_RL._rpo_thruster_allocator(config)
    scheduler = SpaceAGORA_RL._rpo_thruster_pulse_scheduler(config)
    force_n = zeros(6, n_steps)
    for index in 1:n_steps
        actuation = SpaceAGORA_RL._thruster_translation_step(
            view(commands, :, index), view(attitudes, :, index + 1),
            dt, config, allocator, scheduler,
        )
        force_n[:, index] .= actuation.thruster_impulse_ns ./ dt
    end
    return collect(0:(n_steps - 1)) .* dt, force_n
end

function _write_hypr_mpc_diagnostic_html(path, case_index::Int, plan,
                                          tracking_time_s, tracking_error_m,
                                          thrust_time_s, thrust_force_n)
    traces = Any[
        Dict(
            "type" => "scatter", "mode" => "lines", "name" => "position error",
            "x" => tracking_time_s, "y" => tracking_error_m,
            "line" => Dict("color" => "rgb(200, 65, 55)", "width" => 2),
            "xaxis" => "x", "yaxis" => "y",
        ),
    ]
    colors = [
        "rgb(34, 94, 168)", "rgb(92, 148, 216)", "rgb(41, 136, 89)",
        "rgb(117, 193, 125)", "rgb(166, 96, 195)", "rgb(218, 154, 219)",
    ]
    for thruster_index in 1:6
        push!(traces, Dict(
            "type" => "scatter", "mode" => "lines",
            "name" => "thruster $(thruster_index)",
            "x" => thrust_time_s, "y" => vec(thrust_force_n[thruster_index, :]),
            "line" => Dict("color" => colors[thruster_index], "width" => 1.5),
            "xaxis" => "x2", "yaxis" => "y2",
        ))
    end
    title = @sprintf(
        "HYPR original/proxy case %03d | fuel %.4f g | final position error %.3e m",
        case_index, 1_000.0 * plan.propellant_used_kg,
        plan.diagnostics.final_position_error_m,
    )
    layout = Dict(
        "title" => title,
        "height" => 760,
        "margin" => Dict("l" => 75, "r" => 35, "t" => 70, "b" => 60),
        "legend" => Dict("orientation" => "h", "y" => -0.16),
        "grid" => Dict("rows" => 2, "columns" => 1, "pattern" => "independent"),
        "xaxis" => Dict("title" => "time (s)"),
        "yaxis" => Dict("title" => "position tracking error (m)"),
        "xaxis2" => Dict("title" => "time (s)"),
        "yaxis2" => Dict("title" => "average thruster force (N)"),
    )
    open(path, "w") do io
        print(io, """<!doctype html><html><head><meta charset=\"utf-8\"><title>$(title)</title>
<script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script></head>
<body style=\"margin:0\"><div id=\"diagnostic\" style=\"width:100vw;height:100vh\"></div>
<script>Plotly.newPlot('diagnostic', $(JSON.json(traces)), $(JSON.json(layout)), {responsive:true});</script>
</body></html>""")
    end
    return path
end

function plot_hypr_mpc_case_diagnostics(
    config_path::AbstractString=joinpath(@__DIR__, "..", "..", "configs", "rpo", "hypr_rl.toml");
    case_index::Integer=1,
    output_directory::AbstractString=joinpath(
        @__DIR__, "..", "..", "outputs", "hypr_rl", "mpc_diagnostics",
    ),
)
    case_index > 0 || throw(ArgumentError("case_index must be positive"))
    raw = TOML.parsefile(config_path)
    config, sampler, evaluation_seed = _hypr_mpc_diagnostic_config(raw)
    scenario_seed = evaluation_seed + 2 * Int(case_index)
    planner_seed = scenario_seed + 1
    scenario = sample_rpo_hypr_rl_scenario(sampler, MersenneTwister(scenario_seed))
    result = evaluate_hypr_original_baseline_case(config, scenario, planner_seed)
    plan = result.plan
    plan.valid || throw(ErrorException("case $(case_index) did not produce a valid terminal trajectory"))

    _, _, tracking, _ = SpaceAGORA_RL._rpo_spaceagora_settings(scenario, config)
    tracking_time_s, tracking_error_m = _hypr_mpc_tracking_error(plan, tracking.dt_s)
    thrust_time_s, thrust_force_n = _hypr_mpc_thruster_forces(plan, config, tracking.dt_s)

    mkpath(output_directory)
    prefix = joinpath(output_directory, @sprintf("case_%03d", case_index))
    error_csv = "$(prefix)_mpc_tracking_error.csv"
    thrust_csv = "$(prefix)_thruster_force.csv"
    html_path = "$(prefix)_mpc_diagnostics.html"
    CSV.write(error_csv, DataFrame(time_s=tracking_time_s, position_tracking_error_m=tracking_error_m))
    CSV.write(thrust_csv, DataFrame(
        time_s=thrust_time_s,
        thruster_1_n=vec(thrust_force_n[1, :]),
        thruster_2_n=vec(thrust_force_n[2, :]),
        thruster_3_n=vec(thrust_force_n[3, :]),
        thruster_4_n=vec(thrust_force_n[4, :]),
        thruster_5_n=vec(thrust_force_n[5, :]),
        thruster_6_n=vec(thrust_force_n[6, :]),
    ))
    _write_hypr_mpc_diagnostic_html(
        html_path, Int(case_index), plan,
        tracking_time_s, tracking_error_m, thrust_time_s, thrust_force_n,
    )
    println("HYPR case $(case_index) diagnostic")
    println("  scenario_seed=$(scenario_seed), planner_seed=$(planner_seed)")
    println("  fuel_g=$(1_000.0 * plan.propellant_used_kg)")
    println("  tracking_error_csv=$(abspath(error_csv))")
    println("  thrust_csv=$(abspath(thrust_csv))")
    println("  html=$(abspath(html_path))")
    return (html=html_path, tracking_error_csv=error_csv, thrust_csv=thrust_csv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    config_path = isempty(ARGS) ?
        joinpath(@__DIR__, "..", "..", "configs", "rpo", "hypr_rl.toml") : ARGS[1]
    case_index = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
    output_directory = length(ARGS) >= 3 ? ARGS[3] : joinpath(
        @__DIR__, "..", "..", "outputs", "hypr_rl", "mpc_diagnostics",
    )
    plot_hypr_mpc_case_diagnostics(
        config_path; case_index=case_index, output_directory=output_directory,
    )
end
