using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

const MU = 3.986004418e14; const C = 3.0e8; const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include(joinpath(@__DIR__, "..", "..", "functions", "4_Diagnostics.jl"))
include(joinpath(@__DIR__, "..", "..", "functions", "12_CSV_Write_Read.jl"))
const ż = SVector(0.0, 0.0, 1.0)
include(joinpath(@__DIR__, "..", "..", "functions", "5_OE_Converters.jl"))

# ── Tuning Knobs ──────────────────────────────────────────────────────────────
csv_file = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", 
    "target_h1010km_i0.0deg", 
    "timeseries_N2_T7200s_h1000km_t1010km_ih0.0deg_it0.0deg_B100_Pin1e+04_rmin0m_rmax2e+05_J2T.csv"))
title_fontsize   = 12
xlabel_fontsize  = 10
ylabel_fontsize  = 10
zlabel_fontsize  = 10
grid_alpha       = 0.4          # grid line darkness (0=invisible, 1=solid)
# ─────────────────────────────────────────────────────────────────────────────
Plots.default(gridalpha=grid_alpha)

# Load CSV data
sol, _, p = load_timeseries_csv(csv_file)
println("Loaded CSV from: $csv_file")
println("  Satellites: $(length(sol.u[1]) ÷ 6)")
println("  Time points: $(length(sol.t))")

N_sat = length(sol.u[1]) ÷ 6  # total satellites; target = last
_N_label  = "_N$(N_sat)"
_J2_label = occursin("_J2T", basename(csv_file)) ? "_J2T" : occursin("_J2F", basename(csv_file)) ? "_J2F" : ""
target_idx = N_sat
helper_idx = 1

# Extract relative positions: helper position - target position
rel_x = Float64[]
rel_y = Float64[]
rel_z = Float64[]
rel_dist = Float64[]

# Extract target absolute positions
target_x = Float64[]
target_y = Float64[]
target_z = Float64[]
target_dist = Float64[]

# Extract helper absolute positions
helper_x = Float64[]
helper_y = Float64[]
helper_z = Float64[]

for u in sol.u
    r_target = @SVector [u[idx(target_idx,1)], u[idx(target_idx,2)], u[idx(target_idx,3)]]
    r_helper = @SVector [u[idx(helper_idx,1)], u[idx(helper_idx,2)], u[idx(helper_idx,3)]]
    rel_pos = r_helper - r_target
    
    push!(rel_x, rel_pos[1])
    push!(rel_y, rel_pos[2])
    push!(rel_z, rel_pos[3])
    push!(rel_dist, norm(rel_pos))
    
    push!(target_x, r_target[1])
    push!(target_y, r_target[2])
    push!(target_z, r_target[3])
    push!(target_dist, norm(r_target))

    push!(helper_x, r_helper[1])
    push!(helper_y, r_helper[2])
    push!(helper_z, r_helper[3])
end

# Convert time to hours for plotting
t_hours = sol.t ./ 3600

# ── Panel 1: Relative Position X over time ────────────────────────────────────
p1 = plot(t_hours, rel_x ./ 1000, label="ΔX", color=:blue, lw=2,
          xlabel="Time, hours", ylabel="ΔX, km",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          grid=true, legend=:best)

# ── Panel 2: Relative Position Y over time ────────────────────────────────────
p2 = plot(t_hours, rel_y ./ 1000, label="ΔY", color=:blue, lw=2,
          xlabel="Time, hours", ylabel="ΔY, km",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          grid=true, legend=:best)

# ── Panel 3: Relative Position Z over time ────────────────────────────────────
p3 = plot(t_hours, rel_z ./ 1000, label="ΔZ", color=:blue, lw=2,
          xlabel="Time, hours", ylabel="ΔZ, km",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          grid=true, legend=:best)

# ── Panel 4: Relative Distance over time ──────────────────────────────────────
p4 = plot(t_hours, rel_dist ./ 1000, label="Distance", color=:blue, lw=2,
          xlabel="Time, hours", ylabel="Distance, km",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          grid=true, legend=:best)

# Combine time series plots
fig_ts = plot(p1, p2, p3, p4, layout=grid(2, 2), size=(1200, 800))

# ── 3D Trajectory Plot ────────────────────────────────────────────────────────
p_3d = plot3d(rel_x ./ 1000, rel_y ./ 1000, rel_z ./ 1000, 
              label="Helper Relative Trajectory", 
              color=:blue, lw=2, 
              xlabel="ΔX, km", ylabel="ΔY, km", zlabel="ΔZ, km",
              xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize, 
              zguidefontsize=zlabel_fontsize,
              size=(900, 800),
              margin=5Plots.mm)

# Mark start and end points
scatter3d!([rel_x[1]/1000], [rel_y[1]/1000], [rel_z[1]/1000], 
           label="Start", markercolor=:white, markerstrokecolor=:blue,
           markerstrokewidth=2, markersize=8)
scatter3d!([rel_x[end]/1000], [rel_y[end]/1000], [rel_z[end]/1000], 
           label="End", color=:blue, markersize=8)

# Mark the target at the origin (relative frame)
scatter3d!([0.0], [0.0], [0.0],
           label="Target", color=:red, marker=:circle, markersize=8)

# Save figures
output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "relative_position"))
mkpath(output_dir)

ts_path = joinpath(output_dir, "relative_position_timeseries$(_N_label)$(_J2_label).png")
ts_path_pdf = joinpath(output_dir, "relative_position_timeseries$(_N_label)$(_J2_label).pdf")
savefig(fig_ts, ts_path)
savefig(fig_ts, ts_path_pdf)
println("Saved to: $ts_path")
println("Saved to: $ts_path_pdf")

traj_path = joinpath(output_dir, "relative_position_3d_trajectory$(_N_label)$(_J2_label).png")
traj_path_pdf = joinpath(output_dir, "relative_position_3d_trajectory$(_N_label)$(_J2_label).pdf")
savefig(p_3d, traj_path)
savefig(p_3d, traj_path_pdf)
println("Saved to: $traj_path")
println("Saved to: $traj_path_pdf")

display(fig_ts)
display(p_3d)

# ── Absolute Position X over time (helper vs target) ──────────────────────────
p_ax = plot(t_hours, helper_x ./ 1000, label="Helper", color=:blue, lw=2,
            xlabel="Time, hours", ylabel="X, km",
            xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
            grid=true, legend=:best)
plot!(p_ax, t_hours, target_x ./ 1000, label="Target", color=:red, lw=2)

# ── Absolute Position Y over time (helper vs target) ──────────────────────────
p_ay = plot(t_hours, helper_y ./ 1000, label="Helper", color=:blue, lw=2,
            xlabel="Time, hours", ylabel="Y, km",
            xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
            grid=true, legend=:best)
plot!(p_ay, t_hours, target_y ./ 1000, label="Target", color=:red, lw=2)

# ── Absolute Position Z over time (helper vs target) ──────────────────────────
p_az = plot(t_hours, helper_z ./ 1000, label="Helper", color=:blue, lw=2,
            xlabel="Time, hours", ylabel="Z, km",
            xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
            grid=true, legend=:best)
plot!(p_az, t_hours, target_z ./ 1000, label="Target", color=:red, lw=2)

# Combine into 3x1 layout
fig_abs = plot(p_ax, p_ay, p_az, layout=grid(3, 1), size=(900, 1000))

abs_path = joinpath(output_dir, "absolute_position_xyz_timeseries$(_N_label)$(_J2_label).png")
abs_path_pdf = joinpath(output_dir, "absolute_position_xyz_timeseries$(_N_label)$(_J2_label).pdf")
savefig(fig_abs, abs_path)
savefig(fig_abs, abs_path_pdf)
println("Saved to: $abs_path")
println("Saved to: $abs_path_pdf")

display(fig_abs)

# ── Relative Position X over time (helper relative to target) ─────────────────
p_rx = plot(t_hours, zeros(length(t_hours)), label="Target", color=:red, lw=2, linestyle=:dash,
            xlabel="Time, hours", ylabel="X, km",
            xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
            grid=true, legend=:best)
plot!(p_rx, t_hours, rel_x ./ 1000, label="ΔX (Helper − Target)", color=:blue, lw=2)

# ── Relative Position Y over time (helper relative to target) ─────────────────
p_ry = plot(t_hours, zeros(length(t_hours)), label="Target", color=:red, lw=2, linestyle=:dash,
            xlabel="Time, hours", ylabel="Y, km",
            xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
            grid=true, legend=:best)
plot!(p_ry, t_hours, rel_y ./ 1000, label="ΔY (Helper − Target)", color=:blue, lw=2)

# ── Relative Position Z over time (helper relative to target) ─────────────────
p_rz = plot(t_hours, zeros(length(t_hours)), label="Target", color=:red, lw=2, linestyle=:dash,
            xlabel="Time, hours", ylabel="Z, km",
            xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
            grid=true, legend=:best)
plot!(p_rz, t_hours, rel_z ./ 1000, label="ΔZ (Helper − Target)", color=:blue, lw=2)

# Combine into 3x1 layout
fig_rel = plot(p_rx, p_ry, p_rz, layout=grid(3, 1), size=(900, 1000))

rel_xyz_path = joinpath(output_dir, "relative_position_xyz_timeseries$(_N_label)$(_J2_label).png")
rel_xyz_path_pdf = joinpath(output_dir, "relative_position_xyz_timeseries$(_N_label)$(_J2_label).pdf")
savefig(fig_rel, rel_xyz_path)
savefig(fig_rel, rel_xyz_path_pdf)
println("Saved to: $rel_xyz_path")
println("Saved to: $rel_xyz_path_pdf")

display(fig_rel)

# ══ Combined Figure Tuning Knobs ══════════════════════════════════════════════
combined_fig_size          = (1600, 950)  # overall (width, height) in pixels
combined_3d_width_frac     = 0.45         # fraction of width for the 3D panel (0–1)
combined_3x1_panel_heights = [1/3, 1/3, 1/3]  # relative heights of ΔX, ΔY, ΔZ panels (auto-normalised)
# 3D panel font/marker sizes
combined_3d_xlabel_fontsize   = 13         # x-axis label font size (3D panel)
combined_3d_ylabel_fontsize   = 13         # y-axis label font size (3D panel)
combined_3d_zlabel_fontsize   = 13         # z-axis label font size (3D panel)
combined_3d_tick_fontsize     = 13         # tick label font size (3D panel)
combined_3d_legend_fontsize    = 13         # legend font size (3D panel)
combined_3d_legend_markersize = 15         # legend marker size (3D panel)
combined_3d_legend_position   = (0.1, 0.1)  # legend anchor: symbol (:topright, :topleft, :bottomright, :best, …)
                                              #   OR (x, y) tuple — fractions of panel width/height, e.g.
                                              #   (0.0, 0.0)=bottom-left  (0.5, 0.5)=centre  (1.0, 1.0)=top-right
# 3×1 time-series panels font/marker sizes
combined_ts_xlabel_fontsize   = 13         # x-axis label font size (3×1 panels)
combined_ts_ylabel_fontsize   = 13         # y-axis label font size (3×1 panels)
combined_ts_tick_fontsize     = 13         # tick label font size (3×1 panels)
combined_ts_legend_fontsize   = 13         # legend font size (3×1 panels)
combined_ts_legend_markersize = 13         # legend marker size (3×1 panels)
# ═════════════════════════════════════════════════════════════════════════════

# ── Build combined-specific subplots ─────────────────────────────────────────
pc_3d = plot3d(rel_x ./ 1000, rel_y ./ 1000, rel_z ./ 1000,
               label="Helper trajectory", color=:blue, lw=2,
               xlabel="ΔX, km", ylabel="ΔY, km", zlabel="ΔZ, km",
               xguidefontsize=combined_3d_xlabel_fontsize,
               yguidefontsize=combined_3d_ylabel_fontsize,
               zguidefontsize=combined_3d_zlabel_fontsize,
               xtickfontsize=combined_3d_tick_fontsize,
               ytickfontsize=combined_3d_tick_fontsize,
               ztickfontsize=combined_3d_tick_fontsize,
               legendfontsize=combined_3d_legend_fontsize,
               legend_markersize=combined_3d_legend_markersize,
               margin=4Plots.mm, legend=combined_3d_legend_position)
scatter3d!(pc_3d, [rel_x[1]/1000],   [rel_y[1]/1000],   [rel_z[1]/1000],
           label="Start", markercolor=:white, markerstrokecolor=:blue,
           markerstrokewidth=1, markersize=8)
scatter3d!(pc_3d, [rel_x[end]/1000], [rel_y[end]/1000], [rel_z[end]/1000],
           label="End",    color=:blue,    markersize=8)
scatter3d!(pc_3d, [0.0], [0.0], [0.0],
           label="Target", color=:red, marker=:circle, markersize=8)

pc_rx = plot(t_hours, zeros(length(t_hours)), label="Target", color=:red, lw=2, linestyle=:dash,
             xlabel="Time, hours", ylabel="X, km",
             xguidefontsize=combined_ts_xlabel_fontsize,
             yguidefontsize=combined_ts_ylabel_fontsize,
             xtickfontsize=combined_ts_tick_fontsize,
             ytickfontsize=combined_ts_tick_fontsize,
             legendfontsize=combined_ts_legend_fontsize,
             legend_markersize=combined_ts_legend_markersize,
             grid=true, legend=:best)
plot!(pc_rx, t_hours, rel_x ./ 1000, label="ΔX (Helper − Target)", color=:blue, lw=2)

pc_ry = plot(t_hours, zeros(length(t_hours)), label="Target", color=:red, lw=2, linestyle=:dash,
             xlabel="Time, hours", ylabel="Y, km",
             xguidefontsize=combined_ts_xlabel_fontsize,
             yguidefontsize=combined_ts_ylabel_fontsize,
             xtickfontsize=combined_ts_tick_fontsize,
             ytickfontsize=combined_ts_tick_fontsize,
             legendfontsize=combined_ts_legend_fontsize,
             legend_markersize=combined_ts_legend_markersize,
             grid=true, legend=:best)
plot!(pc_ry, t_hours, rel_y ./ 1000, label="ΔY (Helper − Target)", color=:blue, lw=2)

pc_rz = plot(t_hours, zeros(length(t_hours)), label="Target", color=:red, lw=2, linestyle=:dash,
             xlabel="Time, hours", ylabel="Z, km",
             xguidefontsize=combined_ts_xlabel_fontsize,
             yguidefontsize=combined_ts_ylabel_fontsize,
             xtickfontsize=combined_ts_tick_fontsize,
             ytickfontsize=combined_ts_tick_fontsize,
             legendfontsize=combined_ts_legend_fontsize,
             legend_markersize=combined_ts_legend_markersize,
             grid=true, legend=:best)
plot!(pc_rz, t_hours, rel_z ./ 1000, label="ΔZ (Helper − Target)", color=:blue, lw=2)

# ── Combined: 3D trajectory + 3×1 relative position time series ──────────────
_ph = combined_3x1_panel_heights ./ sum(combined_3x1_panel_heights)
l_combined = eval(Meta.parse(
    "@layout [a{$(combined_3d_width_frac)w} " *
    "[b{$(_ph[1])h}; c{$(_ph[2])h}; d{$(_ph[3])h}]]"))
fig_combined = plot(pc_3d, pc_rx, pc_ry, pc_rz,
                    layout=l_combined, size=combined_fig_size,
                    margin=4Plots.mm)

combined_path     = joinpath(output_dir, "combined_3d_and_relative_xyz$(_N_label)$(_J2_label).png")
combined_path_pdf = joinpath(output_dir, "combined_3d_and_relative_xyz$(_N_label)$(_J2_label).pdf")
savefig(fig_combined, combined_path)
savefig(fig_combined, combined_path_pdf)
println("Saved to: $combined_path")
println("Saved to: $combined_path_pdf")

display(fig_combined)