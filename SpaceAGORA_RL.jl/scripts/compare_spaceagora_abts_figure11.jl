#!/usr/bin/env julia

ENV["GKSwstype"] = "100"

using CSV
using DataFrames
using Plots
using Printf
using Statistics

const PACKAGE_DIR = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_INPUT = joinpath(
    PACKAGE_DIR,
    "outputs",
    "hyperparameter_sweeps",
    "a2c_one_factor_20260814",
    "common_evaluation",
    "multi_run_flight_performance_comparison.csv",
)
const DEFAULT_OUTPUT = joinpath(
    PACKAGE_DIR,
    "outputs",
    "comparisons",
    "spaceagora_vs_abts_figure11_20260817",
)
const DEFAULT_EPISODE_METRICS = joinpath(
    dirname(DEFAULT_INPUT),
    "trained_policies_vs_aads",
    "episode_metrics.csv",
)
const PAPER_PATH = joinpath(
    PACKAGE_DIR,
    "Reference Code (ABTS)",
    "Ref_Paper_Autonomous_Decision-Making_for_Aerobraking_via_Parallel_Rand.pdf",
)

function usage(io::IO=stdout)
    println(io, """
Usage:
  julia --project=SpaceAGORA_RL.jl \\
    SpaceAGORA_RL.jl/scripts/compare_spaceagora_abts_figure11.jl \\
    [SPACEAGORA_COMPARISON_CSV] [OUTPUT_DIRECTORY]

Defaults:
  input:  $DEFAULT_INPUT
  output: $DEFAULT_OUTPUT
  signed-distance episodes: $DEFAULT_EPISODE_METRICS
""")
end

function find_row(frame::DataFrame, predicate, description::AbstractString)
    rows = frame[[predicate(String(policy)) for policy in frame.policy], :]
    nrow(rows) == 1 || error("expected one $description row, found $(nrow(rows))")
    return only(eachrow(rows))
end

function signed_goal_distance(frame::DataFrame, policy::AbstractString)
    rows = frame[frame.policy .== policy, :]
    isempty(rows) && error("no episode metrics found for policy $policy")
    values = Float64.(rows.target_error_km)
    return (mean(values), std(values), length(values))
end

function spaceagora_record(
    row,
    method::AbstractString,
    goal_distance;
    notes::AbstractString="",
)
    target_percent = ismissing(row.target_reached_10km_percent) ?
        missing : Float64(row.target_reached_10km_percent)
    target_error = ismissing(row.mean_abs_final_target_error_km) ?
        missing : Float64(row.mean_abs_final_target_error_km)
    return (
        source="SpaceAGORA",
        method=String(method),
        episodes=Int(row.episodes),
        maneuver_count_mean=Float64(row.mean_maneuver_count),
        maneuver_count_std=Float64(row.std_maneuver_count),
        duration_days_mean=Float64(row.mean_duration_days),
        duration_days_std=Float64(row.std_duration_days),
        total_delta_v_mps_mean=Float64(row.mean_total_delta_v_mps),
        total_delta_v_mps_std=Float64(row.std_total_delta_v_mps),
        thermal_violations_mean=Float64(row.mean_thermal_violations),
        thermal_violations_std=Float64(row.std_thermal_violations),
        target_reached_10km_percent=ismissing(target_percent) || !isfinite(target_percent) ?
            missing : target_percent,
        mean_abs_final_target_error_km=ismissing(target_error) || !isfinite(target_error) ?
            missing : target_error,
        std_abs_final_target_error_km=ismissing(row.std_abs_final_target_error_km) ||
            !isfinite(Float64(row.std_abs_final_target_error_km)) ?
            missing : Float64(row.std_abs_final_target_error_km),
        mean_signed_goal_distance_km=goal_distance[1],
        std_signed_goal_distance_km=goal_distance[2],
        goal_distance_episodes=goal_distance[3],
        goal_distance_source=ismissing(goal_distance[1]) ? missing :
            "SpaceAGORA episode metrics",
        notes=String(notes),
    )
end

function abts_record(
    method::AbstractString,
    episodes::Integer,
    maneuvers,
    duration,
    delta_v,
    thermal;
    goal_distance=(missing, missing, missing),
    notes::AbstractString="Digitized from Figure 11; values are approximate.",
)
    return (
        source="ABTS paper",
        method=String(method),
        episodes=Int(episodes),
        maneuver_count_mean=Float64(maneuvers[1]),
        maneuver_count_std=Float64(maneuvers[2]),
        duration_days_mean=Float64(duration[1]),
        duration_days_std=Float64(duration[2]),
        total_delta_v_mps_mean=Float64(delta_v[1]),
        total_delta_v_mps_std=Float64(delta_v[2]),
        thermal_violations_mean=Float64(thermal[1]),
        thermal_violations_std=Float64(thermal[2]),
        target_reached_10km_percent=missing,
        mean_abs_final_target_error_km=missing,
        std_abs_final_target_error_km=missing,
        mean_signed_goal_distance_km=goal_distance[1],
        std_signed_goal_distance_km=goal_distance[2],
        goal_distance_episodes=goal_distance[3],
        goal_distance_source=ismissing(goal_distance[1]) ? missing :
            "Figure 8(b), policy at 1.1M training steps",
        notes=String(notes),
    )
end

function comparison_data(input_path::AbstractString, episode_metrics_path::AbstractString)
    frame = CSV.read(input_path, DataFrame)
    episodes = CSV.read(episode_metrics_path, DataFrame)
    a2c = find_row(
        frame,
        policy -> occursin("a2c_sweep_rollout_length_16", policy),
        "rollout-16 A2C",
    )
    pr_drl = find_row(frame, policy -> policy == "PR-DRL", "PR-DRL")
    aads = find_row(frame, policy -> policy == "AADS", "AADS")
    odyssey = find_row(frame, policy -> policy == "Mars Odyssey", "Mars Odyssey")

    rows = [
        spaceagora_record(
            a2c,
            "A2C",
            signed_goal_distance(episodes, "trained_a2c_5");
            notes="Current best completed tuning: rollout 16, entropy 0.01, learning rate 1e-4, successful repetitions 9; checkpoint step 205056.",
        ),
        spaceagora_record(
            pr_drl,
            "PR-DRL",
            signed_goal_distance(episodes, "trained_pr_drl");
            notes="Same comparison checkpoint used previously; checkpoint step 505000.",
        ),
        spaceagora_record(aads, "AADS", signed_goal_distance(episodes, "aads_heuristic")),
        spaceagora_record(odyssey, "Odyssey", (missing, missing, missing)),
        abts_record(
            "PR-DRL",
            100,
            (22.7, 4.4),
            (24.3, 1.8),
            (28.0, 1.1),
            (0.24, 0.53);
            goal_distance=(4.8, 3.3, 40),
            notes="Flight metrics digitized from Figure 11. Signed goal distance digitized separately from Figure 8(b) at the selected 1.1M-step policy.",
        ),
        abts_record(
            "AADS",
            100,
            (11.9, 1.5),
            (17.8, 3.5),
            (31.5, 2.0),
            (13.0, 9.3);
            notes="Digitized from Figure 11. The plotted total delta-V is about 31.5 m/s, although the adjacent prose implies 34.2 m/s.",
        ),
        abts_record(
            "Odyssey",
            1,
            (15.0, 0.0),
            (23.3, 0.0),
            (31.8, 0.0),
            (18.0, 0.0);
            notes="Fixed Mars Odyssey flight-data reference plotted in Figure 11.",
        ),
    ]
    return DataFrame(rows)
end

function metric_values(
    frame::DataFrame,
    source::AbstractString,
    mean_field::Symbol,
    std_field::Symbol,
)
    methods = ["A2C", "PR-DRL", "AADS", "Odyssey"]
    values = Float64[]
    errors = Float64[]
    positions = Float64[]
    for (position, method) in enumerate(methods)
        rows = frame[(frame.source .== source) .& (frame.method .== method), :]
        isempty(rows) && continue
        row = only(eachrow(rows))
        (ismissing(row[mean_field]) || ismissing(row[std_field])) && continue
        value = Float64(row[mean_field])
        error = Float64(row[std_field])
        (!isfinite(value) || !isfinite(error)) && continue
        push!(positions, position)
        push!(values, value)
        push!(errors, error)
    end
    return positions, values, errors
end

function value_label(value::Real)
    return abs(value) < 1 ? @sprintf("%.2f", value) : @sprintf("%.1f", value)
end

function metric_panel(frame::DataFrame, field::Symbol, title::AbstractString, ylabel::AbstractString;
                      show_legend::Bool=false)
    methods = ["A2C", "PR-DRL", "AADS", "Odyssey"]
    std_field = Symbol(replace(String(field), "_mean" => "_std"))
    space_x, space_y, space_error = metric_values(frame, "SpaceAGORA", field, std_field)
    paper_x, paper_y, paper_error = metric_values(frame, "ABTS paper", field, std_field)
    space_x .-= 0.19
    paper_x .+= 0.19
    upper_bound = maximum(vcat(space_y .+ space_error, paper_y .+ paper_error))
    ymax = upper_bound == 0 ? 1.0 : 1.19 * upper_bound

    panel = plot(
        ;
        title=title,
        ylabel=ylabel,
        xticks=(1:4, methods),
        xlims=(0.55, 4.45),
        ylims=(0, ymax),
        legend=show_legend ? :topright : false,
        grid=:y,
        gridalpha=0.22,
        framestyle=:box,
        bottom_margin=6Plots.mm,
        left_margin=5Plots.mm,
    )
    bar!(
        panel,
        space_x,
        space_y;
        bar_width=0.34,
        color="#0072B2",
        linecolor="#0072B2",
        label=show_legend ? "SpaceAGORA" : false,
    )
    bar!(
        panel,
        paper_x,
        paper_y;
        bar_width=0.34,
        color="#D55E00",
        linecolor="#D55E00",
        label=show_legend ? "ABTS paper" : false,
    )

    for (xs, ys, errors) in (
        (space_x, space_y, space_error),
        (paper_x, paper_y, paper_error),
    )
        lower_errors = min.(errors, ys)
        scatter!(
            panel,
            xs,
            ys;
            yerror=(lower_errors, errors),
            markershape=:none,
            color=:black,
            linewidth=1.5,
            label=false,
        )
        for (x, y, error) in zip(xs, ys, errors)
            annotate!(
                panel,
                x,
                min(y + error + 0.025 * ymax, 0.965 * ymax),
                text(value_label(y), 8, :black, :center),
            )
        end
    end
    annotate!(panel, 1.19, 0.035 * ymax, text("not in paper", 7, :gray40, :center))
    return panel
end

function signed_goal_distance_panel(frame::DataFrame)
    methods = ["A2C", "PR-DRL", "AADS", "Odyssey"]
    mean_field = :mean_signed_goal_distance_km
    std_field = :std_signed_goal_distance_km
    space_x, space_y, space_error = metric_values(
        frame,
        "SpaceAGORA",
        mean_field,
        std_field,
    )
    paper_x, paper_y, paper_error = metric_values(frame, "ABTS paper", mean_field, std_field)
    space_x .-= 0.19
    paper_x .+= 0.19
    lower = minimum(vcat(space_y .- space_error, paper_y .- paper_error, [0.0]))
    upper = maximum(vcat(space_y .+ space_error, paper_y .+ paper_error, [0.0]))
    padding = 0.12 * max(upper - lower, 1.0)
    ymin = lower - padding
    ymax = upper + padding

    panel = plot(
        ;
        title="(e) Signed final distance to goal",
        ylabel="km (+ above target)",
        xticks=(1:4, methods),
        xlims=(0.55, 4.45),
        ylims=(ymin, ymax),
        legend=false,
        grid=:y,
        gridalpha=0.22,
        framestyle=:box,
        bottom_margin=6Plots.mm,
        left_margin=5Plots.mm,
    )
    hline!(panel, [0.0]; color=:gray35, linewidth=1.2, label=false)
    for (xs, ys, errors, color) in (
        (space_x, space_y, space_error, "#0072B2"),
        (paper_x, paper_y, paper_error, "#D55E00"),
    )
        bar!(panel, xs, ys; bar_width=0.34, color=color, linecolor=color, label=false)
        scatter!(
            panel,
            xs,
            ys;
            yerror=errors,
            markershape=:none,
            color=:black,
            linewidth=1.5,
            label=false,
        )
        for (x, y) in zip(xs, ys)
            offset = 0.025 * (ymax - ymin)
            label_y = y >= 0 ? y + offset : y - offset
            annotate!(panel, x, label_y, text(value_label(y), 8, :black, :center))
        end
    end
    return panel
end

function source_note_panel()
    panel = plot(
        ;
        framestyle=:none,
        xlims=(0, 1),
        ylims=(0, 1),
        xticks=false,
        yticks=false,
        legend=false,
    )
    annotate!(
        panel,
        0.5,
        0.58,
        text(
            "Distance-panel source\n\nABTS PR-DRL: Figure 8(b),\nselected 1.1M-step policy,\n40 testing episodes.\nDifferent test set from Figure 11.\n\nABTS AADS and Odyssey means\nare not reported in the paper.",
            10,
            :gray25,
            :center,
        ),
    )
    return panel
end

function write_report(
    path::AbstractString,
    frame::DataFrame,
    input_path::AbstractString,
    episode_metrics_path::AbstractString,
)
    spaceagora = frame[frame.source .== "SpaceAGORA", :]
    a2c = only(eachrow(spaceagora[spaceagora.method .== "A2C", :]))
    pr_drl = only(eachrow(spaceagora[spaceagora.method .== "PR-DRL", :]))
    aads = only(eachrow(spaceagora[spaceagora.method .== "AADS", :]))
    open(path, "w") do io
        println(io, "# SpaceAGORA versus ABTS paper results")
        println(io)
        println(io, "SpaceAGORA source: `$(abspath(input_path))`.")
        println(io)
        println(io, "SpaceAGORA signed-distance source: `$(abspath(episode_metrics_path))`.")
        println(io)
        println(io, "ABTS Figure 11 source: `$(abspath(PAPER_PATH))`, PDF page 13 (journal page 3067).")
        println(io)
        println(io, "ABTS signed-distance source: Figure 8(b) in the same PDF, PDF page 11 (journal page 3065).")
        println(io)
        println(io, "Current best completed A2C tuning: rollout length 16, entropy coefficient 0.01, learning rate 1e-4, successful-case repetitions 9; checkpoint step 205056.")
        println(io)
        println(io, "| Source | Method | Maneuvers | Duration (days) | Total ΔV (m/s) | Thermal violations | Target ≤10 km | Mean absolute target error (km) | Mean signed goal distance (km) |")
        println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|")
        for row in eachrow(frame)
            reached = ismissing(row.target_reached_10km_percent) ? "—" : @sprintf("%.1f%%", row.target_reached_10km_percent)
            target_error = ismissing(row.mean_abs_final_target_error_km) ? "—" :
                @sprintf("%.2f ± %.2f", row.mean_abs_final_target_error_km, row.std_abs_final_target_error_km)
            signed_distance = ismissing(row.mean_signed_goal_distance_km) ? "—" :
                @sprintf("%+.2f ± %.2f", row.mean_signed_goal_distance_km, row.std_signed_goal_distance_km)
            @printf(
                io,
                "| %s | %s | %.2f ± %.2f | %.2f ± %.2f | %.2f ± %.2f | %.2f ± %.2f | %s | %s | %s |\n",
                row.source,
                row.method,
                row.maneuver_count_mean,
                row.maneuver_count_std,
                row.duration_days_mean,
                row.duration_days_std,
                row.total_delta_v_mps_mean,
                row.total_delta_v_mps_std,
                row.thermal_violations_mean,
                row.thermal_violations_std,
                reached,
                target_error,
                signed_distance,
            )
        end
        println(io)
        println(io, "ABTS flight-performance values are approximate visual digitizations from Figure 11 in Falcone and Putnam (2023). Error bars are standard deviations. Figure 11 truncates thermal lower error bars at zero.")
        println(io)
        println(io, "The signed distance is final apoapsis minus target apoapsis, so a negative value indicates overshooting below the target. The ABTS PR-DRL value (approximately +4.8 ± 3.3 km) is digitized from Figure 8(b) at the selected 1.1-million-step policy and averages 40 test episodes. The paper does not report corresponding mean values for AADS or Odyssey.")
        println(io)
        println(io, "Figure 8(b) and Figure 11 do not use the same evaluation set: Figure 8(b) is a 40-episode policy test during training, whereas Figure 11 is the 100-episode Mars Odyssey flight-condition comparison. The distance panel is therefore contextual rather than a controlled Figure 11 comparison.")
        println(io)
        println(io, "The paper's AADS ΔV bar is approximately 31.5 m/s, while its adjacent prose implies 34.2 m/s. This report keeps the plotted value because the requested source is Figure 11.")
        println(io)
        println(io, "Odyssey thermal counts are not directly equivalent across sources: SpaceAGORA's reference is 9, while Figure 11 plots 18. The paper also states that Odyssey used a more lenient heat-rate corridor than PR-DRL and AADS.")
        println(io)
        println(io, "## SpaceAGORA A2C summary")
        println(io)
        @printf(
            io,
            "Against SpaceAGORA PR-DRL, A2C uses %.1f%% less ΔV and has %.1f%% fewer thermal violations while finishing %.1f%% sooner, but reaches the ±10 km target %.1f percentage points less often and has %.1f%% greater mean absolute target error.\n",
            100 * (1 - a2c.total_delta_v_mps_mean / pr_drl.total_delta_v_mps_mean),
            100 * (1 - a2c.thermal_violations_mean / pr_drl.thermal_violations_mean),
            100 * (1 - a2c.duration_days_mean / pr_drl.duration_days_mean),
            pr_drl.target_reached_10km_percent - a2c.target_reached_10km_percent,
            100 * (a2c.mean_abs_final_target_error_km / pr_drl.mean_abs_final_target_error_km - 1),
        )
        println(io)
        @printf(
            io,
            "Against SpaceAGORA AADS, A2C uses %.1f%% less ΔV and has %.1f%% fewer thermal violations, reaches the ±10 km target %.1f percentage points more often, and has %.1f%% lower mean absolute target error; its campaign is %.1f%% longer.\n",
            100 * (1 - a2c.total_delta_v_mps_mean / aads.total_delta_v_mps_mean),
            100 * (1 - a2c.thermal_violations_mean / aads.thermal_violations_mean),
            a2c.target_reached_10km_percent - aads.target_reached_10km_percent,
            100 * (1 - a2c.mean_abs_final_target_error_km / aads.mean_abs_final_target_error_km),
            100 * (a2c.duration_days_mean / aads.duration_days_mean - 1),
        )
        println(io)
        println(io, "Cross-source differences should not be interpreted as algorithm-only effects: the SpaceAGORA evaluation and the ABTS paper runs do not share identical dynamics, action accounting, atmosphere realizations, or thermal-count definitions.")
    end
end

function main(args=ARGS)
    if any(arg -> arg in ("-h", "--help"), args)
        usage()
        return
    end
    length(args) <= 2 || error("expected at most two arguments; use --help")
    input_path = abspath(length(args) >= 1 ? args[1] : DEFAULT_INPUT)
    output_dir = abspath(length(args) >= 2 ? args[2] : DEFAULT_OUTPUT)
    episode_metrics_path = abspath(joinpath(
        dirname(input_path),
        "trained_policies_vs_aads",
        "episode_metrics.csv",
    ))
    isfile(input_path) || error("comparison CSV not found: $input_path")
    isfile(episode_metrics_path) || error("episode metrics not found: $episode_metrics_path")
    mkpath(output_dir)

    frame = comparison_data(input_path, episode_metrics_path)
    csv_path = joinpath(output_dir, "spaceagora_vs_abts_figure11.csv")
    CSV.write(csv_path, frame)

    panels = [
        metric_panel(frame, :maneuver_count_mean, "(a) Maneuvers", "Count"; show_legend=true),
        metric_panel(frame, :duration_days_mean, "(b) Campaign duration", "Days"),
        metric_panel(frame, :total_delta_v_mps_mean, "(c) Total ΔV", "m/s"),
        metric_panel(frame, :thermal_violations_mean, "(d) Thermal violations", "Count"),
        signed_goal_distance_panel(frame),
        source_note_panel(),
    ]
    figure = plot(
        panels...;
        layout=(2, 3),
        size=(1800, 1000),
        dpi=180,
        plot_title="SpaceAGORA best A2C comparison with ABTS paper results",
        plot_titlefontsize=15,
        fontfamily="sans-serif",
        guidefontsize=11,
        tickfontsize=9,
        titlefontsize=11,
        legendfontsize=9,
    )
    png_path = joinpath(output_dir, "spaceagora_vs_abts_figure11_comparison.png")
    pdf_path = joinpath(output_dir, "spaceagora_vs_abts_figure11_comparison.pdf")
    savefig(figure, png_path)
    savefig(figure, pdf_path)

    report_path = joinpath(output_dir, "README.md")
    write_report(report_path, frame, input_path, episode_metrics_path)
    println("wrote $csv_path")
    println("wrote $png_path")
    println("wrote $pdf_path")
    println("wrote $report_path")
end

main()
