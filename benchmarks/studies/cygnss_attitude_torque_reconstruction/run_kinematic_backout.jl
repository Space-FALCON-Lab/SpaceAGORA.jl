##
# Kinematic torque back-out: the reaction-wheel torque required to explain
# the measured CYGNSS attitude trajectory in a fully torque-free (no
# gravity-gradient/aero/magnetic) environment, computed directly from the
# measured q(t)/w(t) kinematics via the modified Euler equation (derivation
# in common.jl), then replayed through SpaceAGORA's engine to check the
# trajectory match. See README.md ("Kinematic torque back-out") for results
# and interpretation.
#
# This script only runs the simulation and writes the plotted quantities to
# data/plot_data/kinematic_backout.arrow -- run make_plots.jl separately to
# generate the figures without re-running the simulation.
#
# Usage: julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_kinematic_backout.jl
##

include(joinpath(@__DIR__, "common.jl"))

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

vel_sim = [SVector{3, Float64}(sol(t - T_CAL).sc[1].vel) for t in sample_t]
vel_gt = [v_truth(t) for t in sample_t]
vel_err_km_s = [norm(vel_sim[i] - vel_gt[i]) / 1000.0 for i in eachindex(sample_t)]

println("Kinematic-derived torque magnitude: mean=$(mean(tau_mag)) N*m, max=$(maximum(tau_mag)) N*m")
println("Reconstruction error vs telemetry (full hour): mean=$(mean(angle_err_deg)) deg, max=$(maximum(angle_err_deg)) deg")
println("Position error vs telemetry (full hour, gravity-only propagation): mean=$(mean(pos_err_km)) km, max=$(maximum(pos_err_km)) km")
println("Velocity error vs telemetry (full hour, gravity-only propagation): mean=$(mean(vel_err_km_s)) km/s, max=$(maximum(vel_err_km_s)) km/s")

# ==============================================================================
# Save plotted quantities (see make_plots.jl for the figures)
# ==============================================================================

running_rms_deg = sqrt.(cumsum(angle_err_deg .^ 2) ./ collect(1:length(angle_err_deg)))

data_out = DataFrame(
    t_s=t_s,
    tau_x=tau_x, tau_y=tau_y, tau_z=tau_z, tau_mag=tau_mag,
    q1_sim=[q[1] for q in q_sim_sf], q2_sim=[q[2] for q in q_sim_sf], q3_sim=[q[3] for q in q_sim_sf], q4_sim=[q[4] for q in q_sim_sf],
    q1_gt=[q[1] for q in q_gt_sf], q2_gt=[q[2] for q in q_gt_sf], q3_gt=[q[3] for q in q_gt_sf], q4_gt=[q[4] for q in q_gt_sf],
    angle_err_deg=angle_err_deg, running_rms_deg=running_rms_deg,
    omega_x_sim=[w[1] for w in ω_sim], omega_y_sim=[w[2] for w in ω_sim], omega_z_sim=[w[3] for w in ω_sim],
    omega_x_gt=[w[1] for w in ω_gt], omega_y_gt=[w[2] for w in ω_gt], omega_z_gt=[w[3] for w in ω_gt],
    pos_x_sim=[p[1] / 1000.0 for p in pos_sim], pos_y_sim=[p[2] / 1000.0 for p in pos_sim], pos_z_sim=[p[3] / 1000.0 for p in pos_sim],
    pos_x_gt=[p[1] / 1000.0 for p in pos_gt], pos_y_gt=[p[2] / 1000.0 for p in pos_gt], pos_z_gt=[p[3] / 1000.0 for p in pos_gt],
    pos_err_km=pos_err_km,
    vel_x_sim=[v[1] / 1000.0 for v in vel_sim], vel_y_sim=[v[2] / 1000.0 for v in vel_sim], vel_z_sim=[v[3] / 1000.0 for v in vel_sim],
    vel_x_gt=[v[1] / 1000.0 for v in vel_gt], vel_y_gt=[v[2] / 1000.0 for v in vel_gt], vel_z_gt=[v[3] / 1000.0 for v in vel_gt],
    vel_err_km_s=vel_err_km_s,
)
data_path = joinpath(PLOT_DATA_DIR, "kinematic_backout.arrow")
Arrow.write(data_path, data_out)
println("plot data written to: $(data_path)")

rmse_out = DataFrame(
    key=["kinematic_backout"], label=["Kinematic torque back-out"],
    rmse_pos_km=[result.rmse_pos_km], rmse_vel_km_s=[result.rmse_vel_km_s],
    rmse_angle_deg=[result.rmse_angle_deg], rmse_ang_vel_deg_s=[result.rmse_ang_vel_deg_s],
)
rmse_path = joinpath(PLOT_DATA_DIR, "kinematic_backout_rmse.arrow")
Arrow.write(rmse_path, rmse_out)
println("RMSE table data written to: $(rmse_path)")
