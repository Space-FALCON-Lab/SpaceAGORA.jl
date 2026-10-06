ENV["GKSwstype"] = "100"
using Plots
using Statistics

# Plot the recent mean episode reward recorded at each terminal progress report.
runs_dir = normpath(joinpath(@__DIR__, "..", "outputs", "runs"))
output_dir = normpath(joinpath(runs_dir, "..", "training_reward_charts"))
mkpath(output_dir)
index = ["# Training reward versus iteration", "",
    "Iteration is the logged `train_steps` optimizer-update counter. Reward is the logged `recent_reward`: mean total episode reward over the configured recent-episode window. Pale curves show this logged mean; dark curves show an additional trailing mean over up to 100 logged samples. Logs sample training periodically, so these are not rewards for individual optimizer updates. Repeated iteration counters retain the latest sample, including iteration zero during warmup. Runs may use different reward configurations.", "",
    "| Run | Samples | Last iteration | PNG | PDF | Data |", "|---|---:|---:|---|---|---|"]
skipped = String[]
for run_dir in sort(readdir(runs_dir; join=true))
    isdir(run_dir) || continue
    run = basename(run_dir)
    log_path = joinpath(run_dir, "terminal_output.txt")
    if !isfile(log_path)
        push!(skipped, "$run: no terminal log")
        continue
    end
    iterations, steps, rewards = Int[], Int[], Float64[]
    for line in eachline(log_path)
        startswith(line, "progress ") || continue
        fields = Dict(m.captures[1] => m.captures[2] for m in eachmatch(r"(\w+)=([^\s]+)", line))
        reward = tryparse(Float64, get(fields, "recent_reward", "n/a"))
        (reward === nothing || !isfinite(reward)) && continue
        iteration = parse(Int, fields["train_steps"])
        step = parse(Int, first(split(fields["steps"], '/')))
        if !isempty(iterations) && iteration == last(iterations)
            rewards[end], steps[end] = reward, step
        else
            push!(iterations, iteration)
            push!(steps, step)
            push!(rewards, reward)
        end
    end
    if isempty(rewards)
        push!(skipped, "$run: no finite reward measurements in terminal log")
        continue
    end
    @assert issorted(iterations) "Optimizer counter resets in $run"
    window = min(100, length(rewards))
    smooth = [mean(@view rewards[max(1, i-window+1):i]) for i in eachindex(rewards)]
    algorithm = uppercase(replace(last(split(run, '-')), '_' => '-'))
    title = "$algorithm — training reward\n$(first(split(run, '_')))"
    ylabel = "Recent mean episode reward"
    chart = plot(iterations, rewards; label="Logged recent mean reward", color=:steelblue,
        alpha=0.3, linewidth=0.7, xlabel="Training iteration (optimizer updates)",
        ylabel=ylabel, title=title, size=(1100, 650), dpi=180,
        legend=:topright, gridalpha=0.15, framestyle=:box,
        margin=6Plots.mm, titlefontsize=12)
    plot!(chart, iterations, smooth; label="Trailing mean ($window logged samples)",
        color=:navy, linewidth=2)
    stem = joinpath(output_dir, run)
    savefig(chart, stem * ".png")
    savefig(chart, stem * ".pdf")
    open(stem * ".csv", "w") do io
        println(io, "iteration,environment_steps,recent_mean_episode_reward,trailing_mean")
        for i in eachindex(rewards)
            println(io, "$(iterations[i]),$(steps[i]),$(rewards[i]),$(smooth[i])")
        end
    end
    push!(index, "| $run | $(length(rewards)) | $(last(iterations)) | [PNG]($run.png) | [PDF]($run.pdf) | [CSV]($run.csv) |")
    println("Plotted $run: $(length(rewards)) samples")
end
append!(index, ["", "## Runs without logged reward data", ""])
append!(index, ["- $item" for item in skipped])
write(joinpath(output_dir, "README.md"), join(index, '\n') * "\n")
println("Charts and index: $output_dir")
