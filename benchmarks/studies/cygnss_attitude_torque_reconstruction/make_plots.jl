##
# Generates every figure for the CYGNSS attitude/torque reconstruction study
# from the data written by the run_*.jl scripts to data/plot_data/*.arrow.
# Doesn't re-run any simulation, so it's fast to iterate on plot styling --
# run the relevant run_*.jl script first if the underlying data is missing
# or stale.
#
# Usage: julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/make_plots.jl
##

using Plots
using Arrow
using DataFrames
using Statistics

const STUDY_DIR = @__DIR__
const PLOTS_DIR = joinpath(STUDY_DIR, "plots")
const PLOT_DATA_DIR = joinpath(STUDY_DIR, "data", "plot_data")
mkpath(PLOTS_DIR)

Plots.default(
    size=(850, 850),
    xtickfontsize=14, ytickfontsize=14, guidefontsize=14, legendfontsize=12,
    left_margin=2Plots.mm, bottom_margin=2Plots.mm, top_margin=2Plots.mm, right_margin=2Plots.mm,
)

_load(name) = DataFrame(Arrow.Table(joinpath(PLOT_DATA_DIR, name)))

# ==============================================================================
# Kinematic torque back-out
# ==============================================================================

kb = _load("kinematic_backout.arrow")

torque_fig = plot(
    kb.t_s, kb.tau_x;
    label="tau_x", lw=1.4, color=:steelblue,
    xlabel="Time since t=0s (s)", ylabel="Torque (N*m)", legend=:topright,
)
plot!(torque_fig, kb.t_s, kb.tau_y; label="tau_y", lw=1.4, color=:darkorange)
plot!(torque_fig, kb.t_s, kb.tau_z; label="tau_z", lw=1.4, color=:seagreen)
savefig(torque_fig, joinpath(PLOTS_DIR, "cygnss_kinematic_torque_components.pdf"))

torque_mag_fig = plot(
    kb.t_s, kb.tau_mag;
    label="|tau|", lw=1.6, color=:purple,
    xlabel="Time since t=0s (s)", ylabel="Torque magnitude (N*m)", legend=:topright,
)
savefig(torque_mag_fig, joinpath(PLOTS_DIR, "cygnss_kinematic_torque_magnitude.pdf"))

for (i, name) in enumerate(("q1", "q2", "q3", "q4")) # scalar-first: q1 is the scalar component
    sim_col, gt_col = Symbol("q$(i)_sim"), Symbol("q$(i)_gt")
    sp = plot(
        kb.t_s, kb[!, sim_col];
        label="kinematic-backout replay", lw=1.4, color=:steelblue,
        xlabel="Time since t=0s (s)", ylabel="value", legend=:best,
    )
    plot!(sp, kb.t_s, kb[!, gt_col]; label="telemetry", lw=1.0, ls=:dash, color=:darkorange)
    savefig(sp, joinpath(PLOTS_DIR, "cygnss_kinematic_backout_quaternion_$(name)_timeseries.pdf"))
end

err_fig = plot(
    kb.t_s, kb.angle_err_deg;
    label="angle error", lw=1.6, color=:steelblue,
    xlabel="Time since t=0s (s)", ylabel="Quaternion angle error (deg)", legend=:topleft,
)
plot!(err_fig, kb.t_s, kb.running_rms_deg; label="running RMS", lw=1.6, ls=:dash, color=:black)
hline!(err_fig, [sqrt(mean(kb.angle_err_deg .^ 2))]; label="total RMS", lw=1.0, ls=:dot, color=:gray)
savefig(err_fig, joinpath(PLOTS_DIR, "cygnss_kinematic_backout_attitude_error_timeseries.pdf"))

for name in ("omega_x", "omega_y", "omega_z")
    sim_col, gt_col = Symbol("$(name)_sim"), Symbol("$(name)_gt")
    sp = plot(
        kb.t_s, kb[!, sim_col];
        label="replay", lw=1.4, color=:steelblue,
        xlabel="Time since t=0s (s)", ylabel="rad/s", legend=:best,
    )
    plot!(sp, kb.t_s, kb[!, gt_col]; label="telemetry", lw=1.0, ls=:dash, color=:darkorange)
    savefig(sp, joinpath(PLOTS_DIR, "cygnss_kinematic_backout_angular_rate_$(name)_timeseries.pdf"))
end

for name in ("x", "y", "z")
    sim_col, gt_col = Symbol("pos_$(name)_sim"), Symbol("pos_$(name)_gt")
    sp = plot(
        kb.t_s, kb[!, sim_col];
        label="simulated", lw=1.4, color=:steelblue,
        xlabel="Time since t=0s (s)", ylabel="km", legend=:best,
    )
    plot!(sp, kb.t_s, kb[!, gt_col]; label="telemetry", lw=1.0, ls=:dash, color=:darkorange)
    savefig(sp, joinpath(PLOTS_DIR, "cygnss_kinematic_backout_position_$(name)_timeseries.pdf"))
end

pos_err_fig = plot(
    kb.t_s, kb.pos_err_km;
    label="position error", lw=1.6, color=:steelblue,
    xlabel="Time since t=0s (s)", ylabel="Position error (km)", legend=:topleft,
)
savefig(pos_err_fig, joinpath(PLOTS_DIR, "cygnss_kinematic_backout_position_error_timeseries.pdf"))

vel_err_fig = plot(
    kb.t_s, kb.vel_err_km_s;
    label="velocity error", lw=1.6, color=:steelblue,
    xlabel="Time since t=0s (s)", ylabel="Velocity error (km/s)", legend=:topleft,
)
savefig(vel_err_fig, joinpath(PLOTS_DIR, "cygnss_kinematic_backout_velocity_error_timeseries.pdf"))

println("kinematic torque back-out plots written")

# ==============================================================================
# Position accuracy: IC-correction x gravity-fidelity ablation + best case
# ==============================================================================

ablation = _load("position_accuracy_ablation.arrow")
xs = 1:nrow(ablation)
bar_fig = bar(
    xs .- 0.15, ablation.pos_mean_km;
    bar_width=0.3, label="mean error", color=:steelblue,
    ylabel="Position error (km)", xlabel="",
    xticks=(xs, ablation.label), xtickfontsize=10, legend=:topright,
)
bar!(bar_fig, xs .+ 0.15, ablation.pos_max_km; bar_width=0.3, label="max error", color=:darkorange)
savefig(bar_fig, joinpath(PLOTS_DIR, "cygnss_position_accuracy_ablation.pdf"))

vel_bar_fig = bar(
    xs .- 0.15, ablation.vel_mean_km_s;
    bar_width=0.3, label="mean error", color=:steelblue,
    ylabel="Velocity error (km/s)", xlabel="",
    xticks=(xs, ablation.label), xtickfontsize=10, legend=:topright,
)
bar!(vel_bar_fig, xs .+ 0.15, ablation.vel_max_km_s; bar_width=0.3, label="max error", color=:darkorange)
savefig(vel_bar_fig, joinpath(PLOTS_DIR, "cygnss_velocity_accuracy_ablation.pdf"))

best_case = _load("position_accuracy_best_case.arrow")
for name in ("x", "y", "z")
    sim_col, gt_col = Symbol("pos_$(name)_sim"), Symbol("pos_$(name)_gt")
    sp = plot(
        best_case.t_s, best_case[!, sim_col];
        label="simulated (SMA-corrected IC, 50x50 EGM)", lw=1.4, color=:steelblue,
        xlabel="Time since t=0s (s)", ylabel="km", legend=:best,
    )
    plot!(sp, best_case.t_s, best_case[!, gt_col]; label="telemetry", lw=1.0, ls=:dash, color=:darkorange)
    savefig(sp, joinpath(PLOTS_DIR, "cygnss_position_accuracy_best_case_$(name)_timeseries.pdf"))
end

pos_err_fig2 = plot(
    best_case.t_s, best_case.pos_err_km;
    label="position error (SMA-corrected IC, 50x50 EGM)", lw=1.6, color=:steelblue,
    xlabel="Time since t=0s (s)", ylabel="Position error (km)", legend=:topleft,
)
savefig(pos_err_fig2, joinpath(PLOTS_DIR, "cygnss_position_accuracy_best_case_error_timeseries.pdf"))

vel_err_fig2 = plot(
    best_case.t_s, best_case.vel_err_km_s;
    label="velocity error (SMA-corrected IC, 50x50 EGM)", lw=1.6, color=:steelblue,
    xlabel="Time since t=0s (s)", ylabel="Velocity error (km/s)", legend=:topleft,
)
savefig(vel_err_fig2, joinpath(PLOTS_DIR, "cygnss_position_accuracy_best_case_velocity_error_timeseries.pdf"))

println("position accuracy plots written")

# ==============================================================================
# Independent-effect quantification
# ==============================================================================

summary = _load("independent_effects_summary.arrow")
xs2 = 1:nrow(summary)
bar_fig2 = bar(
    xs2 .- 0.15, summary.attitude_mean_deg;
    bar_width=0.3, label="mean error", color=:steelblue,
    ylabel="Quaternion angle error (deg)", xlabel="",
    xticks=(xs2, summary.label), xtickfontsize=10, xrotation=25, legend=:topright,
)
bar!(bar_fig2, xs2 .+ 0.15, summary.attitude_max_deg; bar_width=0.3, label="max error", color=:darkorange)
savefig(bar_fig2, joinpath(PLOTS_DIR, "cygnss_independent_effects_summary.pdf"))

pos_bar_fig2 = bar(
    xs2 .- 0.15, summary.pos_mean_km;
    bar_width=0.3, label="mean error", color=:steelblue,
    ylabel="Position error (km)", xlabel="",
    xticks=(xs2, summary.label), xtickfontsize=10, xrotation=25, legend=:topright,
)
bar!(pos_bar_fig2, xs2 .+ 0.15, summary.pos_max_km; bar_width=0.3, label="max error", color=:darkorange)
savefig(pos_bar_fig2, joinpath(PLOTS_DIR, "cygnss_independent_effects_position_summary.pdf"))

vel_bar_fig2 = bar(
    xs2 .- 0.15, summary.vel_mean_km_s;
    bar_width=0.3, label="mean error", color=:steelblue,
    ylabel="Velocity error (km/s)", xlabel="",
    xticks=(xs2, summary.label), xtickfontsize=10, xrotation=25, legend=:topright,
)
bar!(vel_bar_fig2, xs2 .+ 0.15, summary.vel_max_km_s; bar_width=0.3, label="max error", color=:darkorange)
savefig(vel_bar_fig2, joinpath(PLOTS_DIR, "cygnss_independent_effects_velocity_summary.pdf"))

ts = _load("independent_effects_position_timeseries.arrow")
pos_ts_fig = plot(xlabel="Time since t=0s (s)", ylabel="Position error (km)", legend=:topleft)
colors = (:steelblue, :darkorange, :seagreen, :purple, :firebrick)
for (row, c) in zip(eachrow(summary), colors)
    plot!(pos_ts_fig, ts.t_s, ts[!, Symbol(row.key)]; label=row.label, lw=1.4, color=c)
end
savefig(pos_ts_fig, joinpath(PLOTS_DIR, "cygnss_independent_effects_position_timeseries.pdf"))

vel_ts = _load("independent_effects_velocity_timeseries.arrow")
vel_ts_fig = plot(xlabel="Time since t=0s (s)", ylabel="Velocity error (km/s)", legend=:topleft)
for (row, c) in zip(eachrow(summary), colors)
    plot!(vel_ts_fig, vel_ts.t_s, vel_ts[!, Symbol(row.key)]; label=row.label, lw=1.4, color=c)
end
savefig(vel_ts_fig, joinpath(PLOTS_DIR, "cygnss_independent_effects_velocity_timeseries.pdf"))

println("independent-effects plots written")

# ==============================================================================
# 48hr reference replication
# ==============================================================================

r48 = _load("position_accuracy_48hr.arrow")
err_fig48 = plot(
    r48.t_hr, r48.ex;
    label="x error", lw=1.2, alpha=0.9, color=:steelblue,
    xlabel="Time (hr)", ylabel="Error (km)", legend=:topleft,
)
plot!(err_fig48, r48.t_hr, r48.ey; label="y error", lw=1.2, alpha=0.9, color=:darkorange)
plot!(err_fig48, r48.t_hr, r48.ez; label="z error", lw=1.2, alpha=0.9, color=:seagreen)
plot!(err_fig48, r48.t_hr, r48.running_total_rmse; label="total running RMSE", lw=2.0, color=:black)
savefig(err_fig48, joinpath(PLOTS_DIR, "cygnss_48hr_reference_replication_error_timeseries.pdf"))

println("48hr reference replication plot written")

println()
println("All plots written to: $(PLOTS_DIR)")
