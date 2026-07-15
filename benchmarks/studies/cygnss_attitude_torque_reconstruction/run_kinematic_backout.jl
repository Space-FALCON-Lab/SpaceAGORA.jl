##
# Kinematic torque back-out: the reaction-wheel torque required to explain
# the measured CYGNSS attitude trajectory in a fully torque-free (no
# gravity-gradient/aero/magnetic) environment, computed directly from the
# measured q(t)/w(t) kinematics via the modified Euler equation (derivation
# in common.jl), then replayed through SpaceAGORA's engine to check the
# trajectory match. See README.md ("Kinematic torque back-out") for results
# and interpretation.
#
# Usage: julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_kinematic_backout.jl
##

include(joinpath(@__DIR__, "common.jl"))
using Plots
using Printf

Plots.default(left_margin=10Plots.mm, bottom_margin=8Plots.mm)

const T_CAL = 0.0
const T_END = 3600.0

println("=== Running kinematic torque back-out case over the full hour ===")
result = run_cygnss_case(
    t_cal=T_CAL, t_end=T_END,
    dynamic_effectors_factory=(t_offset) -> (InverseSquaredGravityModel(gravity_gradient=false),),
    label="kinematic_backout_full_hour", dt_max=1.0,
    control_effectors_factory=(t_offset) -> (KinematicTorqueReplayControlModel(t_offset),),
)
sol = result.sol

# ==============================================================================
# Sample everything on a common time grid
# ==============================================================================

n_samples = 600
sample_t = range(T_CAL, T_END, length=n_samples)
t_s = collect(sample_t) .- T_CAL

torque_kin = [τ_kinematic_backout(t) for t in sample_t]
tau_x = [τ[1] for τ in torque_kin]
tau_y = [τ[2] for τ in torque_kin]
tau_z = [τ[3] for τ in torque_kin]
tau_mag = [norm(τ) for τ in torque_kin]

q_sim_sl = [SVector{4, Float64}(sol(t - T_CAL).sc[1].q) for t in sample_t] # scalar-last
q_sim_sf = [SVector{4, Float64}(q[4], q[1], q[2], q[3]) for q in q_sim_sl] # scalar-first, for comparison
q_gt_sf = [q_truth_scalarfirst(t) for t in sample_t]

angle_err_deg = [
    rad2deg(2.0 * acos(clamp(abs(dot(q_sim_sf[i], q_gt_sf[i])), 0.0, 1.0)))
    for i in eachindex(sample_t)
]

ω_sim = [SVector{3, Float64}(sol(t - T_CAL).sc[1].ω) for t in sample_t]
ω_gt = [ω_meas(t) for t in sample_t]

pos_sim = [SVector{3, Float64}(sol(t - T_CAL).sc[1].pos) for t in sample_t]
pos_gt = [r_truth(t) for t in sample_t]
pos_err_km = [norm(pos_sim[i] - pos_gt[i]) / 1000.0 for i in eachindex(sample_t)]

println("Kinematic-derived torque magnitude: mean=$(mean(tau_mag)) N*m, max=$(maximum(tau_mag)) N*m")
println("Reconstruction error vs telemetry (full hour): mean=$(mean(angle_err_deg)) deg, max=$(maximum(angle_err_deg)) deg")
println("Position error vs telemetry (full hour, gravity-only propagation): mean=$(mean(pos_err_km)) km, max=$(maximum(pos_err_km)) km")

# ==============================================================================
# Plot 1: kinematics-derived torque components
# ==============================================================================

torque_fig = plot(
    t_s, tau_x;
    label="tau_x", lw=1.4, color=:steelblue,
    title="Kinematics-Derived Reaction-Wheel Torque (Torque-Free Back-Out)",
    xlabel="Time since t=0s (s)", ylabel="Torque (N*m)",
    legend=:topright, size=(1000, 550), titlefontsize=11,
)
plot!(torque_fig, t_s, tau_y; label="tau_y", lw=1.4, color=:darkorange)
plot!(torque_fig, t_s, tau_z; label="tau_z", lw=1.4, color=:seagreen)
torque_path = joinpath(PLOTS_DIR, "cygnss_kinematic_torque_components.png")
savefig(torque_fig, torque_path)
println("torque components plot: $(torque_path)")

torque_mag_fig = plot(
    t_s, tau_mag;
    label="|tau|", lw=1.6, color=:purple,
    title="Kinematics-Derived Reaction-Wheel Torque Magnitude",
    xlabel="Time since t=0s (s)", ylabel="Torque magnitude (N*m)",
    legend=:topright, size=(1000, 450), titlefontsize=11,
)
torque_mag_path = joinpath(PLOTS_DIR, "cygnss_kinematic_torque_magnitude.png")
savefig(torque_mag_fig, torque_mag_path)
println("torque magnitude plot: $(torque_mag_path)")

# ==============================================================================
# Plot 2: quaternion time series, kinematic-backout replay vs telemetry
# ==============================================================================

component_names = ("q1 (scalar)", "q2", "q3", "q4")
quat_subplots = map(1:4) do i
    sp = plot(
        t_s, [q[i] for q in q_sim_sf];
        label="kinematic-backout replay", lw=1.4, color=:steelblue,
        title=component_names[i], xlabel="Time since t=0s (s)", ylabel="value",
        legend=(i == 1 ? :best : false),
    )
    plot!(sp, t_s, [q[i] for q in q_gt_sf]; label="telemetry", lw=1.0, ls=:dash, color=:darkorange)
    return sp
end
quat_fig = plot(quat_subplots...; layout=(2, 2), size=(1000, 700), plot_title="CYGNSS Full Hour: Kinematic Torque Back-Out Replay vs. Telemetry")
quat_path = joinpath(PLOTS_DIR, "cygnss_kinematic_backout_quaternion_timeseries.png")
savefig(quat_fig, quat_path)
println("quaternion time series plot: $(quat_path)")

# ==============================================================================
# Plot 3: attitude error time series (full hour)
# ==============================================================================

running_rms_deg = sqrt.(cumsum(angle_err_deg .^ 2) ./ collect(1:length(angle_err_deg)))
err_fig = plot(
    t_s, angle_err_deg;
    label="angle error", lw=1.6, color=:steelblue,
    title="CYGNSS Full-Hour Reconstruction Error: Kinematic Torque Back-Out vs. Telemetry",
    xlabel="Time since t=0s (s)", ylabel="Quaternion angle error (deg)",
    legend=:topleft, size=(1000, 500), titlefontsize=11,
)
plot!(err_fig, t_s, running_rms_deg; label="running RMS", lw=1.6, ls=:dash, color=:black)
hline!(err_fig, [sqrt(mean(angle_err_deg .^ 2))]; label="total RMS", lw=1.0, ls=:dot, color=:gray)
err_path = joinpath(PLOTS_DIR, "cygnss_kinematic_backout_attitude_error_timeseries.png")
savefig(err_fig, err_path)
println("attitude error time series plot: $(err_path)")

# ==============================================================================
# Plot 4: angular velocity comparison (sanity check on the back-out itself)
# ==============================================================================

ω_names = ("omega_x", "omega_y", "omega_z")
ω_subplots = map(1:3) do i
    sp = plot(
        t_s, [w[i] for w in ω_sim];
        label="replay", lw=1.4, color=:steelblue,
        title=ω_names[i], xlabel="Time since t=0s (s)", ylabel="rad/s",
        legend=(i == 1 ? :best : false),
    )
    plot!(sp, t_s, [w[i] for w in ω_gt]; label="telemetry", lw=1.0, ls=:dash, color=:darkorange)
    return sp
end
ω_fig = plot(ω_subplots...; layout=(3, 1), size=(1000, 800), plot_title="CYGNSS Full Hour: Body Rate, Kinematic Torque Back-Out Replay vs. Telemetry")
ω_path = joinpath(PLOTS_DIR, "cygnss_kinematic_backout_angular_rate_timeseries.png")
savefig(ω_fig, ω_path)
println("angular rate time series plot: $(ω_path)")

# ==============================================================================
# Plot 5: position comparison (translational dynamics is gravity-only here --
# this checks the orbit propagation independent of the torque back-out, which
# only affects attitude)
# ==============================================================================

pos_names = ("x", "y", "z")
pos_subplots = map(1:3) do i
    sp = plot(
        t_s, [p[i] / 1000.0 for p in pos_sim];
        label="simulated", lw=1.4, color=:steelblue,
        title=pos_names[i], xlabel="Time since t=0s (s)", ylabel="km",
        legend=(i == 1 ? :best : false),
    )
    plot!(sp, t_s, [p[i] / 1000.0 for p in pos_gt]; label="telemetry", lw=1.0, ls=:dash, color=:darkorange)
    return sp
end
pos_fig = plot(pos_subplots...; layout=(3, 1), size=(1000, 800), plot_title="CYGNSS Full Hour: ECI Position, Simulated vs. Telemetry")
pos_path = joinpath(PLOTS_DIR, "cygnss_kinematic_backout_position_timeseries.png")
savefig(pos_fig, pos_path)
println("position time series plot: $(pos_path)")

pos_err_fig = plot(
    t_s, pos_err_km;
    label="position error", lw=1.6, color=:steelblue,
    title="CYGNSS Full-Hour Position Error: Simulated vs. Telemetry",
    xlabel="Time since t=0s (s)", ylabel="Position error (km)",
    legend=:topleft, size=(1000, 500), titlefontsize=11,
)
pos_err_path = joinpath(PLOTS_DIR, "cygnss_kinematic_backout_position_error_timeseries.png")
savefig(pos_err_fig, pos_err_path)
println("position error time series plot: $(pos_err_path)")

println()
println("All plots written to: $(PLOTS_DIR)")
