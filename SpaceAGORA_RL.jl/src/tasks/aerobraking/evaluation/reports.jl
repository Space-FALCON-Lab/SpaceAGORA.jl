using CSV
using DataFrames

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
