# Writes the two Odyssey figures from a fresh run of the tutorial's default comparison
# (panel caps 90 and 30 degrees). Terminal, from the repository root:
#   julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl --output=output/odyssey_figures
#   julia --project=examples scripts/docs_figures/odyssey_figures.jl output/odyssey_figures
using CSV, DataFrames, Plots

length(ARGS) == 1 || error("Usage: odyssey_figures.jl <comparison directory with cap_90_deg/ and cap_30_deg/>")
dir = ARGS[1]
out = normpath(joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "figures"))
mkpath(out)
baseline = CSV.read(joinpath(dir, "cap_90_deg", "trajectory.csv"), DataFrame)
variant = CSV.read(joinpath(dir, "cap_30_deg", "trajectory.csv"), DataFrame)

# 1. The plotting snippet from the Odyssey tutorial, unchanged except for the file it writes.
drag = plot(baseline.elapsed_s, baseline.drag_N; label="cap 90 deg",
    xlabel="Elapsed time (s)", ylabel="Drag (N)")
plot!(drag, variant.elapsed_s, variant.drag_N; label="cap 30 deg")
height = plot(baseline.elapsed_s, baseline.height_km; label="cap 90 deg",
    xlabel="Elapsed time (s)", ylabel="Ellipsoidal height (km)")
plot!(height, variant.elapsed_s, variant.height_km; label="cap 30 deg")
savefig(plot(drag, height; layout=(2, 1), size=(800, 700)),
    joinpath(out, "odyssey_panel_cap_comparison.png"))

# 2. One pass, annotated: height, drag, panel heat rate and panel angle for both caps.
accent, warm, muted = colorant"#2f5f7f", colorant"#d97706", colorant"#6b7280"
runs = ((90, baseline, accent), (30, variant, warm))
t, h = baseline.elapsed_s, baseline.height_km
i_min = argmin(h)
ei = round(Int, h[1])
common = (; grid=true, gridalpha=0.15, framestyle=:axes, legend=false, tickfontsize=8, guidefontsize=9,
          left_margin=6Plots.mm, xlims=(-15, t[end] + 25))

p_h = plot(t, h; color=accent, lw=2, ylabel="height, km", ylims=(100, 270), common...)
hline!(p_h, [ei]; color=muted, ls=:dash, lw=1)
scatter!(p_h, [t[i_min]], [h[i_min]]; color=accent, ms=4, msw=0)
annotate!(p_h, 18, 222, text("entry interface, EI_km = $ei:\nthe run starts here", 8, :left, :black))
annotate!(p_h, t[end] - 160, 222, text("exit through $ei km\nends the run", 8, :left, :black))
annotate!(p_h, t[i_min], h[i_min] + 38,
          text("periapsis: $(round(Int, h[i_min])) km at $(round(Int, t[i_min])) s", 8, :center, :black))

p_d = plot(; ylabel="drag, N", common..., legend=:topright, legendfontsize=8, foreground_color_legend=nothing)
p_q = plot(; ylabel="panel heat rate,\nW/cm²", common...)
p_a = plot(; ylabel="panel angle, deg", ylims=(0, 105), yticks=[0, 30, 60, 90],
           xlabel="elapsed time since the entry interface, s", common...)
for (cap, df, c) in runs
    plot!(p_d, df.elapsed_s, df.drag_N; color=c, lw=1.8, label="panel cap $(cap)°")
    plot!(p_q, df.elapsed_s, df.heat_left_W_cm2; color=c, lw=1.8)
    annotate!(p_q, 470, 0.62 * maximum(df.heat_left_W_cm2),
              text("heat load $(round(Int, df.heat_left_J_cm2[end])) J/cm²", 8, :left, c))
    plot!(p_a, df.elapsed_s, df.panel_left_deg; seriestype=:steppost, color=c, lw=1.8)
end
annotate!(p_a, 40, 12, text("the first command moves the panels to the cap and holds them there", 8, :left, warm))
annotate!(p_a, 420, 58, text("mode: max energy depletion after the first update\n" *
                             "guidance every 3 s, control every 0.1 s (example settings)", 8, :center, :black))
for p in (p_h, p_d, p_q, p_a)
    vline!(p, [t[i_min]]; color=colorant"#d1d5db", lw=0.8, label="")
end
fig = plot(p_h, p_d, p_q, p_a; layout=grid(4, 1; heights=[0.31, 0.23, 0.23, 0.23]), size=(900, 960),
           plot_title="Anatomy of one aerobraking pass: the Odyssey tutorial with two panel caps",
           plot_titlefontsize=11, link=:x, bottom_margin=2Plots.mm)
savefig(fig, joinpath(out, "aerobraking_pass_anatomy.png"))
println("wrote odyssey_panel_cap_comparison.png and aerobraking_pass_anatomy.png in ", out)
