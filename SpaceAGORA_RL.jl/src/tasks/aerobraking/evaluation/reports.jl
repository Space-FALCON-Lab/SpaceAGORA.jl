ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using CSV
using DataFrames
using Plots

function write_evaluation_artifacts(output_dir::AbstractString, results::Dict)
    mkpath(output_dir)

    metric_rows = NamedTuple[]
    aggregate_rows = NamedTuple[]
    pass_rows = NamedTuple[]
    for (policy, result) in sort(collect(results); by=first)
        append!(metric_rows, result.metrics)
        push!(aggregate_rows, result.aggregate)
        append!(pass_rows, result.pass_rows)
    end

    metrics_path = joinpath(output_dir, "episode_metrics.csv")
    aggregate_path = joinpath(output_dir, "summary_metrics.csv")
    pass_path = joinpath(output_dir, "pass_logs.csv")
    CSV.write(metrics_path, DataFrame(metric_rows))
    CSV.write(aggregate_path, DataFrame(aggregate_rows))
    CSV.write(pass_path, DataFrame(pass_rows))
    return (metrics=metrics_path, aggregate=aggregate_path, pass_logs=pass_path)
end

function write_checkpoint_validation_artifacts(output_dir::AbstractString,
                                               records::AbstractVector{<:NamedTuple},
                                               best::NamedTuple)
    mkpath(output_dir)
    summary_path = joinpath(output_dir, "checkpoint_validation_summary.csv")
    best_path = joinpath(output_dir, "best_validation_checkpoint.txt")
    CSV.write(summary_path, DataFrame(records))
    open(best_path, "w") do io
        println(io, "checkpoint=", best.checkpoint)
        println(io, "checkpoint_path=", best.checkpoint_path)
        println(io, "selection_mode=", best.mode)
        println(io, "greedy=", best.greedy)
        println(io, "episodes=", best.episodes)
        println(io, "success_rate=", best.success_rate)
        println(io, "mean_target_error_km=", best.mean_target_error_km)
        println(io, "thermal_terminal_failure_rate=", best.thermal_terminal_failure_rate)
    end
    return (summary=summary_path, best=best_path)
end

function _save_training_metric_plot(plot_object, output_dir::AbstractString,
                                    filename::AbstractString,
                                    paths::Vector{String})
    path = joinpath(output_dir, filename)
    savefig(plot_object, path)
    push!(paths, path)
    return path
end

function _mean_std_plot(rows::DataFrame, mean_field::Symbol, std_field::Symbol;
                        ylabel::AbstractString, title::AbstractString,
                        label=false, color=:steelblue)
    return plot(
        rows.global_step,
        rows[!, mean_field];
        ribbon=rows[!, std_field],
        xlabel="Training step",
        ylabel=ylabel,
        title=title,
        label=label,
        color=color,
        linewidth=2,
        marker=:circle,
    )
end

"""
Write the checkpoint-evaluation plots used to assess policy learning over time.

The policy-performance panels use the thermal-tolerant evaluation mode: thermal
violations are counted, but do not stop the campaign before its final outcome is
known. Loss is averaged over optimizer updates since the preceding evaluated
checkpoint, and reward is averaged over that checkpoint's evaluation episodes.
"""
function write_checkpoint_training_plots(
    output_dir::AbstractString,
    records::AbstractVector{<:NamedTuple},
)
    isempty(records) && return String[]
    mkpath(output_dir)
    data = sort(DataFrame(records), [:global_step, :mode])
    tolerant = data[data.mode .== "tolerant", :]
    nrow(tolerant) > 0 ||
        throw(ArgumentError("training plots require tolerant checkpoint evaluations"))
    paths = String[]

    completion_values = hcat(
        tolerant.success_percent,
        tolerant.surpassed_target_percent,
        tolerant.impact_percent,
        tolerant.out_of_drag_passage_percent,
        tolerant.unfinished_percent,
    )
    completion = plot(
        tolerant.global_step,
        completion_values;
        xlabel="Training step",
        ylabel="Episodes (%)",
        title="(a) Episode completion",
        label=permutedims([
            "reached goal",
            "surpassed target",
            "impact",
            "out of passage",
            "unfinished",
        ]),
        linewidth=2,
        marker=:circle,
        ylims=(0, 100),
    )
    _save_training_metric_plot(
        completion,
        output_dir,
        "training_a_episode_completion.png",
        paths,
    )

    target = _mean_std_plot(
        tolerant,
        :mean_target_error_km,
        :std_target_error_km;
        ylabel="Absolute final distance (km)",
        title="(b) Final distance to target apoapsis radius",
    )
    _save_training_metric_plot(
        target,
        output_dir,
        "training_b_final_target_distance.png",
        paths,
    )

    thermal = _mean_std_plot(
        tolerant,
        :mean_thermal_violations,
        :std_thermal_violations;
        ylabel="Violations per episode",
        title="(c) Thermal violations",
        color=:firebrick,
    )
    _save_training_metric_plot(
        thermal,
        output_dir,
        "training_c_thermal_violations.png",
        paths,
    )

    delta_v = _mean_std_plot(
        tolerant,
        :mean_delta_v_mps,
        :std_delta_v_mps;
        ylabel="ΔV (m/s)",
        title="(d1) Average ΔV",
        color=:darkorange,
    )
    maneuvers = _mean_std_plot(
        tolerant,
        :mean_maneuver_count,
        :std_maneuver_count;
        ylabel="Maneuvers per episode",
        title="(d2) Number of maneuvers",
        color=:mediumpurple,
    )
    delta_v_and_maneuvers = plot(
        delta_v,
        maneuvers;
        layout=(2, 1),
        size=(900, 800),
    )
    _save_training_metric_plot(
        delta_v_and_maneuvers,
        output_dir,
        "training_d_delta_v_and_maneuvers.png",
        paths,
    )

    passages = _mean_std_plot(
        tolerant,
        :mean_pass_count,
        :std_pass_count;
        ylabel="Passages per episode",
        title="(e1) Number of atmospheric passages",
        color=:seagreen,
    )
    duration = _mean_std_plot(
        tolerant,
        :mean_mission_duration_days,
        :std_mission_duration_days;
        ylabel="Episode duration (days)",
        title="(e2) Episode length",
        color=:sienna,
    )
    passages_and_duration = plot(
        passages,
        duration;
        layout=(2, 1),
        size=(900, 800),
    )
    _save_training_metric_plot(
        passages_and_duration,
        output_dir,
        "training_e_passages_and_episode_length.png",
        paths,
    )

    loss_rows = unique(
        data[:, [
            :global_step,
            :mean_training_loss,
            :training_loss_sum,
            :training_loss_count,
        ]],
        :global_step,
    )
    sort!(loss_rows, :global_step)
    interval_loss = copy(loss_rows.mean_training_loss)
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
    finite_loss = isfinite.(interval_loss)
    averaged_loss = if any(finite_loss)
        plot(
            loss_rows.global_step[finite_loss],
            interval_loss[finite_loss];
            xlabel="Training step",
            ylabel="Mean optimizer loss",
            title="Averaged training loss",
            label=false,
            color=:steelblue,
            linewidth=2,
            marker=:circle,
        )
    else
        plot(
            xlabel="Training step",
            ylabel="Average optimizer loss",
            title="Averaged training loss (not recorded)",
            label=false,
        )
    end
    _save_training_metric_plot(
        averaged_loss,
        output_dir,
        "training_average_loss.png",
        paths,
    )

    averaged_reward = plot(
        xlabel="Training step",
        ylabel="Mean episode reward",
        title="Averaged policy reward",
    )
    colors = Dict("conservative" => :seagreen, "tolerant" => :darkorange)
    for mode in PAPER_EVALUATION_MODES
        mode_rows = data[data.mode .== mode, :]
        nrow(mode_rows) == 0 && continue
        plot!(
            averaged_reward,
            mode_rows.global_step,
            mode_rows.mean_reward;
            ribbon=mode_rows.std_reward,
            label=mode,
            color=colors[mode],
            linewidth=2,
            marker=:circle,
        )
    end
    _save_training_metric_plot(
        averaged_reward,
        output_dir,
        "training_average_reward.png",
        paths,
    )
    return paths
end
