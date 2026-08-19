#!/usr/bin/env julia

"""Render SpaceAGORA generalization metrics as a one-page PDF report."""

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using CSV
using DataFrames
using Plots
using Printf
using TOML

const CASE_LABELS = Dict(
    "iid_reference" => "IID reference",
    "nominal" => "Nominal",
    "exponential_density" => "Exponential density",
    "aggressive_atmosphere" => "Aggressive atmosphere",
    "short_campaign" => "Short campaign",
    "long_campaign" => "Long campaign",
    "high_accuracy_spaceagora" => "High-accuracy SpaceAGORA",
)

function _usage(io::IO=stdout)
    println(io, """
Usage:
  julia --project=SpaceAGORA_RL.jl \\
    SpaceAGORA_RL.jl/scripts/render_generalization_table.jl \\
    ALL_CASES_CSV [OUTPUT_PDF]

ALL_CASES_CSV should be the generalization suite's
all_cases_with_iid_reference.csv. OUTPUT_PDF defaults to
generalization_results_table.pdf beside the CSV.
""")
end

function _mean_std(row, mean_field, std_field; digits=2)
    format = Printf.Format("%.$(digits)f")
    return Printf.format(format, row[mean_field]) * " +/- " *
           Printf.format(format, row[std_field])
end

function _cell_text(row, field; digits=3)
    return Printf.format(Printf.Format("%.$(digits)f"), row[field])
end

function _draw_table!(report, headers, rows;
                      left, right, top, row_height, widths,
                      header_size=9, cell_size=9,
                      highlight_column=nothing)
    length(headers) == length(widths) ||
        throw(ArgumentError("table headers and widths must have equal length"))
    isapprox(sum(widths), 1.0; atol=1e-8) ||
        throw(ArgumentError("table widths must sum to one"))
    edges = left .+ (right - left) .* cumsum(vcat(0.0, widths))
    border = :gray70

    for column in eachindex(headers)
        x0, x1 = edges[column], edges[column + 1]
        y0, y1 = top - row_height, top
        plot!(
            report,
            Shape([x0, x1, x1, x0], [y0, y0, y1, y1]);
            fillcolor="#243B53",
            linecolor=border,
            linewidth=0.5,
            label=false,
        )
        annotate!(
            report,
            (x0 + x1) / 2,
            (y0 + y1) / 2,
            text(headers[column], header_size, :white, :center),
        )
    end

    for (row_index, row) in enumerate(rows)
        y1 = top - row_index * row_height
        y0 = y1 - row_height
        row_fill = iseven(row_index) ? "#F5F7FA" : :white
        for column in eachindex(row)
            x0, x1 = edges[column], edges[column + 1]
            fill = column == highlight_column ? "#D9EAF7" : row_fill
            plot!(
                report,
                Shape([x0, x1, x1, x0], [y0, y0, y1, y1]);
                fillcolor=fill,
                linecolor=border,
                linewidth=0.5,
                label=false,
            )
            alignment = column == 1 ? :left : :center
            x = column == 1 ? x0 + 0.008 * (right - left) : (x0 + x1) / 2
            annotate!(
                report,
                x,
                (y0 + y1) / 2,
                text(row[column], cell_size, :black, alignment),
            )
        end
    end
    return report
end

function _metadata(csv_path)
    manifest_path = joinpath(dirname(dirname(csv_path)), "evaluation_manifest.toml")
    isfile(manifest_path) || return nothing
    return TOML.parsefile(manifest_path)
end

function render_generalization_table(csv_path::AbstractString,
                                     output_path::AbstractString)
    isfile(csv_path) || throw(ArgumentError("metrics CSV does not exist: $csv_path"))
    frame = CSV.read(csv_path, DataFrame)
    isempty(frame) && throw(ArgumentError("metrics CSV contains no rows: $csv_path"))

    labels = [get(CASE_LABELS, String(row.case), replace(String(row.case), "_" => " "))
              for row in eachrow(frame)]
    rows = collect(eachrow(frame))
    case_headers = replace.(labels, " " => "\n")
    iid_index = findfirst(==("iid_reference"), String.(frame.case))
    iid_column = iid_index === nothing ? nothing : iid_index + 1
    primary_rows = [
        [
            "Episodes (N)",
            [string(row.episodes) for row in rows]...,
        ],
        ["TD loss", [_cell_text(row, :evaluation_td_loss; digits=3) for row in rows]...],
        ["Generalization gap", [_cell_text(row, :generalization_gap; digits=3) for row in rows]...],
        ["Reward\n(mean +/- SD)", [_mean_std(row, :mean_reward, :std_reward; digits=2) for row in rows]...],
        [
            "Thermal violations\n(mean +/- SD)",
            [_mean_std(row, :mean_thermal_violations,
                       :std_thermal_violations; digits=2) for row in rows]...,
        ],
        ["Success", [@sprintf("%.1f%%", row.reached_goal_percent) for row in rows]...],
        [
            "Goal distance, km\n(mean +/- SD)",
            [_mean_std(row, :mean_goal_distance_km,
                       :std_goal_distance_km; digits=2) for row in rows]...,
        ],
    ]
    campaign_rows = [
        [
            "Delta-v, m/s\n(mean +/- SD)",
            [_mean_std(row, :mean_delta_v_mps,
                       :std_delta_v_mps; digits=2) for row in rows]...,
        ],
        [
            "Duration, days\n(mean +/- SD)",
            [_mean_std(row, :mean_mission_duration_days,
                       :std_mission_duration_days; digits=2) for row in rows]...,
        ],
        [
            "Passes / episode\n(mean +/- SD)",
            [_mean_std(row, :mean_episode_length,
                       :std_episode_length; digits=2) for row in rows]...,
        ],
        [
            "Maneuvers\n(mean +/- SD)",
            [_mean_std(row, :mean_maneuver_count,
                       :std_maneuver_count; digits=2) for row in rows]...,
        ],
    ]

    metadata = _metadata(csv_path)
    subtitle = if metadata === nothing
        "Frozen-policy SpaceAGORA evaluation"
    else
        checkpoint = basename(String(metadata["checkpoint_path"]))
        wind = String(metadata["gram_wind_mode"])
        episodes = Int(metadata["generalization_episodes"])
        "Frozen PR-DRL policy: $checkpoint  |  $episodes episodes/case  |  " *
        "native SpaceAGORA physics  |  MarsGRAM winds: $wind"
    end

    gr()
    report = plot(
        ;
        size=(1800, 1080),
        xlims=(0, 1),
        ylims=(0, 1),
        framestyle=:none,
        xticks=false,
        yticks=false,
        legend=false,
        background_color=:white,
        foreground_color=:black,
        fontfamily="DejaVu Sans",
    )
    annotate!(report, 0.04, 0.955,
              text("PR-DRL Generalization Evaluation", 22, :navy, :left))
    annotate!(report, 0.04, 0.915, text(subtitle, 10, :gray30, :left))

    annotate!(report, 0.04, 0.865,
              text("Primary generalization metrics", 13, :navy, :left))
    _draw_table!(
        report,
        ["Metric", case_headers...],
        primary_rows;
        left=0.04,
        right=0.96,
        top=0.83,
        row_height=0.055,
        widths=[0.20, fill(0.80 / length(rows), length(rows))...],
        header_size=7,
        cell_size=7,
        highlight_column=iid_column,
    )

    annotate!(report, 0.04, 0.335,
              text("Campaign and control statistics", 13, :navy, :left))
    _draw_table!(
        report,
        ["Metric", case_headers...],
        campaign_rows;
        left=0.04,
        right=0.96,
        top=0.305,
        row_height=0.045,
        widths=[0.20, fill(0.80 / length(rows), length(rows))...],
        header_size=7,
        cell_size=7,
        highlight_column=iid_column,
    )

    annotate!(
        report,
        0.04,
        0.018,
        text(
            "Notes: Generalization gap = evaluation TD loss - IID-reference TD loss. " *
            "Thermal violations are totals per episode.\n" *
            "Goal distance is absolute final " *
            "apoapsis-target error. The blue column is the held-out IID reference.",
            7,
            :gray35,
            :left,
        ),
    )

    mkpath(dirname(output_path))
    savefig(report, output_path)
    return output_path
end

function main(args=ARGS)
    any(==("--help"), args) && return _usage()
    length(args) in (1, 2) || throw(ArgumentError("expected ALL_CASES_CSV [OUTPUT_PDF]"))
    csv_path = abspath(args[1])
    output_path = length(args) == 2 ? abspath(args[2]) :
                  joinpath(dirname(csv_path), "generalization_results_table.pdf")
    path = render_generalization_table(csv_path, output_path)
    println("wrote generalization table PDF: ", path)
    return path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
