ENV["GKSwstype"] = "100"
using CSV
using DataFrames
using Plots

runs_dir = normpath(joinpath(@__DIR__, "..", "outputs", "runs"))
output_dir = normpath(joinpath(runs_dir, "..", "evaluation_reward_charts"))
mkpath(output_dir)
index = ["# Non-exploratory reward versus iteration", "",
    "Saved checkpoint validations with `greedy=true`. Each point is mean total episode reward; shading is ± one episode standard deviation, not a confidence interval. Conservative and tolerant evaluation modes are shown separately. Iteration is `training_loss_count`, which increments alongside `train_steps` in the DDQN, A2C, and A3C learners. Curves connect evaluated checkpoints only; no smoothing or new simulations are applied. CSVs retain environment steps, episode counts, and checkpoint names. Runs may use different reward configurations.", "",
    "| Run | Checkpoints | PNG | PDF | Data |", "|---|---:|---|---|---|"]
skipped = String[]
for run_dir in sort(readdir(runs_dir; join=true))
    isdir(run_dir) || continue
    run = basename(run_dir)
    source = joinpath(run_dir, "checkpoint_validation", "checkpoint_validation_summary.csv")
    if !isfile(source)
        push!(skipped, "$run: no saved checkpoint validation summary")
        continue
    end
    data = CSV.read(source, DataFrame)
    filter!(row -> row.greedy === true && isfinite(row.mean_reward), data)
    if isempty(data)
        push!(skipped, "$run: no finite greedy evaluation rewards")
        continue
    end
    sort!(data, [:mode, :training_loss_count])
    panels = []
    for mode in sort(unique(data.mode))
        rows = data[data.mode .== mode, :]
        @assert all(>=(0), rows.training_loss_count) "Invalid iteration in $run"
        @assert all(>=(0), rows.std_reward) "Invalid reward deviation in $run"
        color = mode == "conservative" ? :seagreen : :darkorange
        panel = plot(rows.training_loss_count, rows.mean_reward;
            ribbon=rows.std_reward, fillalpha=0.18, color=color,
            marker=:circle, markersize=3, linewidth=1.5,
            label="Mean ± episode SD", title="Greedy — $mode",
            xlabel="Training iteration (optimizer updates)", ylabel="Evaluation episode reward",
            legend=:best, gridalpha=0.15, framestyle=:box, margin=6Plots.mm)
        push!(panels, panel)
    end
    algorithm = uppercase(replace(last(split(run, '-')), '_' => '-'))
    chart = plot(panels...; layout=(length(panels), 1),
        size=(1100, 450 * length(panels)), dpi=180,
        plot_title="$algorithm — $(first(split(run, '_')))")
    stem = joinpath(output_dir, run)
    savefig(chart, stem * ".png")
    savefig(chart, stem * ".pdf")
    exported = select(data, :training_loss_count => :iteration,
        :global_step => :environment_steps, :checkpoint, :mode, :greedy,
        :episodes, :mean_reward, :std_reward)
    CSV.write(stem * ".csv", exported)
    count = length(unique(data.checkpoint))
    push!(index, "| $run | $count | [PNG]($run.png) | [PDF]($run.pdf) | [CSV]($run.csv) |")
    println("Plotted $run: $count checkpoints")
end
append!(index, ["", "## Runs without saved greedy reward data", ""])
append!(index, ["- $item" for item in skipped])
write(joinpath(output_dir, "README.md"), join(index, '\n') * "\n")
println("Charts and index: $output_dir")
