using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

const MU = 3.986004418e14; const C = 3.0e8; const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include(joinpath(@__DIR__, "..", "functions", "4_Diagnostics.jl"))
include(joinpath(@__DIR__, "..", "functions", "12_CSV_Write_Read.jl"))
const Ẑ = SVector(0.0, 0.0, 1.0)
include(joinpath(@__DIR__, "..", "functions", "5_OE_Converters.jl"))

# ── Tuning Knobs ──────────────────────────────────────────────────────────────
csv_file = normpath(joinpath(@__DIR__, "..", "output", "CSV", 
    "target_h1010km_i0.0deg", 
    "timeseries_N2_T7200s_h1000km_t1010km_ih0.0deg_it0.0deg_B100_Pin1e+04_rmin0m_rmax2e+05_J2T.csv"))
title_fontsize   = 12
xlabel_fontsize  = 10
ylabel_fontsize  = 10
zlabel_fontsize  = 10
# ─────────────────────────────────────────────────────────────────────────────

# Load CSV data
sol, _, p = load_timeseries_csv(csv_file)
println("Loaded CSV from: $csv_file")
println("  Satellites: $(length(sol.u[1]) ÷ 6)")
println("  Time points: $(length(sol.t))")

N_sat = length(sol.u[1]) ÷ 6  # total satellites; target = last
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
end

# Convert time to hours for plotting
t_hours = sol.t ./ 3600

# ── Panel 1: Relative Position X over time ────────────────────────────────────
p1 = plot(t_hours, rel_x ./ 1000, label="ΔX", color=:red, lw=2,
          xlabel="Time [hours]", ylabel="ΔX [km]",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          title="Relative Position X Component",
          titlefontsize=title_fontsize, grid=true, legend=:best)

# ── Panel 2: Relative Position Y over time ────────────────────────────────────
p2 = plot(t_hours, rel_y ./ 1000, label="ΔY", color=:green, lw=2,
          xlabel="Time [hours]", ylabel="ΔY [km]",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          title="Relative Position Y Component",
          titlefontsize=title_fontsize, grid=true, legend=:best)

# ── Panel 3: Relative Position Z over time ────────────────────────────────────
p3 = plot(t_hours, rel_z ./ 1000, label="ΔZ", color=:blue, lw=2,
          xlabel="Time [hours]", ylabel="ΔZ [km]",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          title="Relative Position Z Component",
          titlefontsize=title_fontsize, grid=true, legend=:best)

# ── Panel 4: Relative Distance over time ──────────────────────────────────────
p4 = plot(t_hours, rel_dist ./ 1000, label="Distance", color=:purple, lw=2,
          xlabel="Time [hours]", ylabel="Distance [km]",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          title="Relative Distance between Satellites",
          titlefontsize=title_fontsize, grid=true, legend=:best)

# Combine time series plots
fig_ts = plot(p1, p2, p3, p4, layout=grid(2, 2), size=(1200, 800))

# ── 3D Trajectory Plot ────────────────────────────────────────────────────────
p_3d = plot3d(rel_x ./ 1000, rel_y ./ 1000, rel_z ./ 1000, 
              label="Helper Relative Trajectory", 
              color=:viridis, lw=2, 
              xlabel="ΔX [km]", ylabel="ΔY [km]", zlabel="ΔZ [km]",
              xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize, 
              zguidefontsize=zlabel_fontsize,
              title="3D Relative Position Trajectory",
              titlefontsize=title_fontsize, 
              size=(900, 800),
              margin=5Plots.mm)

# Mark start and end points
scatter3d!([rel_x[1]/1000], [rel_y[1]/1000], [rel_z[1]/1000], 
           label="Start", color=:green, markersize=8)
scatter3d!([rel_x[end]/1000], [rel_y[end]/1000], [rel_z[end]/1000], 
           label="End", color=:red, markersize=8)

# Mark the target at the origin (relative frame)
scatter3d!([0.0], [0.0], [0.0],
           label="Target", color=:orange, marker=:star5, markersize=12)

# Save figures
output_dir = normpath(joinpath(@__DIR__, "..", "output", "images_reconstructed_from_csv", "relative_position"))
mkpath(output_dir)

ts_path = joinpath(output_dir, "relative_position_timeseries.png")
ts_path_pdf = joinpath(output_dir, "relative_position_timeseries.pdf")
savefig(fig_ts, ts_path)
savefig(fig_ts, ts_path_pdf)
println("Saved to: $ts_path")
println("Saved to: $ts_path_pdf")

traj_path = joinpath(output_dir, "relative_position_3d_trajectory.png")
traj_path_pdf = joinpath(output_dir, "relative_position_3d_trajectory.pdf")
savefig(p_3d, traj_path)
savefig(p_3d, traj_path_pdf)
println("Saved to: $traj_path")
println("Saved to: $traj_path_pdf")

display(fig_ts)
display(p_3d)