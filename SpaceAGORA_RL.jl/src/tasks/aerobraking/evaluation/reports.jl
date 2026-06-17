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
