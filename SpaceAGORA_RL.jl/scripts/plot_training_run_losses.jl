ENV["GKSwstype"] = "100"
using Plots
using Statistics

# Plot the latest optimizer loss recorded at each terminal progress report.
runs_dir = normpath(joinpath(@__DIR__, "..", "outputs", "runs"))
output_dir = normpath(joinpath(runs_dir, "..", "training_loss_charts"))
mkpath(output_dir)
index = ["# Training loss versus iteration", "",
    "Iteration is the logged `train_steps` optimizer-update counter. Pale curves show logged loss; dark curves show a trailing mean over up to 100 logged samples. Logs sample training periodically, so these are not every optimizer update. Linear loss axes preserve negative actor–critic losses. Different algorithms use different objectives; their loss magnitudes are not directly comparable.", "",
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
    iterations, steps, losses = Int[], Int[], Float64[]
    for line in eachline(log_path)
        startswith(line, "progress ") || continue
        fields = Dict(m.captures[1] => m.captures[2] for m in eachmatch(r"(\w+)=([^\s]+)", line))
        loss = tryparse(Float64, get(fields, "loss", "n/a"))
        (loss === nothing || !isfinite(loss)) && continue
        iteration = parse(Int, fields["train_steps"])
        step = parse(Int, first(split(fields["steps"], '/')))
        if !isempty(iterations) && iteration == last(iterations)
            losses[end], steps[end] = loss, step
        else
            push!(iterations, iteration)
            push!(steps, step)
            push!(losses, loss)
        end
    end
    if isempty(losses)
        push!(skipped, "$run: no finite loss measurements in terminal log")
        continue
    end
    @assert issorted(iterations) "Optimizer counter resets in $run"
    window = min(100, length(losses))
    smooth = [mean(@view losses[max(1, i-window+1):i]) for i in eachindex(losses)]
    algorithm = uppercase(replace(last(split(run, '-')), '_' => '-'))
    title = "$algorithm — training loss\n$(first(split(run, '_')))"
    ylabel = occursin("pr_drl", run) ? "Q-learning loss" : "Total actor–critic loss"
    chart = plot(iterations, losses; label="Logged loss", color=:steelblue,
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
        println(io, "iteration,environment_steps,logged_loss,trailing_mean")
        for i in eachindex(losses)
            println(io, "$(iterations[i]),$(steps[i]),$(losses[i]),$(smooth[i])")
        end
    end
    push!(index, "| $run | $(length(losses)) | $(last(iterations)) | [PNG]($run.png) | [PDF]($run.pdf) | [CSV]($run.csv) |")
    println("Plotted $run: $(length(losses)) samples")
end
append!(index, ["", "## Runs without logged loss data", ""])
append!(index, ["- $item" for item in skipped])
write(joinpath(output_dir, "README.md"), join(index, '\n') * "\n")
println("Charts and index: $output_dir")
