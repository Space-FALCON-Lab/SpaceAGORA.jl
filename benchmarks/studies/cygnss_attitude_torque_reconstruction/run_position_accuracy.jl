##
# Position-accuracy replication: applies the exact IC-noise fix from
# test/gmat_scenario_matrix.jl's _build_cygnss_48hr_reference (N=5 vis-viva
# SMA averaging + velocity-magnitude rescale, no direction change) and that
# reference's full 50x50 EGM gravity field (not just J2) to this study's
# telemetry window, to see how close a full-hour CYGNSS orbit propagation
# gets to that reference's ~1.2-1.6 km/48hr position RMSE. See README.md
# ("Position accuracy replication") for the ablation table and result.
#
# This script only runs the simulations and writes the plotted quantities to
# data/plot_data/position_accuracy_ablation.arrow and
# data/plot_data/position_accuracy_best_case.arrow -- run make_plots.jl
# separately to generate the figures without re-running the simulations.
#
# Usage: julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_position_accuracy.jl
##

include(joinpath(@__DIR__, "common.jl"))

const T0 = 0.0
const TEND = 3600.0
const HARMONICS_FILE = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
isfile(HARMONICS_FILE) || error("Missing harmonics file: $HARMONICS_FILE")

make_harmonics_gravity() = GravitationalHarmonicsModel(
    50, 50, HARMONICS_FILE, planet;
    coefficients_normalized=true, j2_source=:file_c20,
)

results = Dict{String, Any}()

println("=== raw IC, two-body only (existing study baseline) ===")
results["raw_twobody"] = run_cygnss_case(
    t_cal=T0, t_end=TEND,
    dynamic_effectors_factory=(t_offset) -> (InverseSquaredGravityModel(gravity_gradient=false),),
    label="raw_twobody", dt_max=1.0, sma_correct_ic=false,
)

println("=== raw IC, + J2 ===")
results["raw_j2"] = run_cygnss_case(
    t_cal=T0, t_end=TEND,
    dynamic_effectors_factory=(t_offset) -> (InverseSquaredJ2GravityModel(gravity_gradient=false),),
    label="raw_j2", dt_max=1.0, sma_correct_ic=false,
)

println("=== raw IC, + 50x50 EGM harmonics ===")
results["raw_harmonics"] = run_cygnss_case(
    t_cal=T0, t_end=TEND,
    dynamic_effectors_factory=(t_offset) -> (make_harmonics_gravity(),),
    label="raw_harmonics", dt_max=1.0, sma_correct_ic=false,
)

println("=== SMA-corrected IC, two-body only ===")
results["sma_twobody"] = run_cygnss_case(
    t_cal=T0, t_end=TEND,
    dynamic_effectors_factory=(t_offset) -> (InverseSquaredGravityModel(gravity_gradient=false),),
    label="sma_twobody", dt_max=1.0, sma_correct_ic=true,
)

println("=== SMA-corrected IC, + J2 ===")
results["sma_j2"] = run_cygnss_case(
    t_cal=T0, t_end=TEND,
    dynamic_effectors_factory=(t_offset) -> (InverseSquaredJ2GravityModel(gravity_gradient=false),),
    label="sma_j2", dt_max=1.0, sma_correct_ic=true,
)

println("=== SMA-corrected IC, + 50x50 EGM harmonics (exact replication of the 48hr reference setup) ===")
results["sma_harmonics"] = run_cygnss_case(
    t_cal=T0, t_end=TEND,
    dynamic_effectors_factory=(t_offset) -> (make_harmonics_gravity(),),
    label="sma_harmonics", dt_max=1.0, sma_correct_ic=true,
)

println()
println("=== SUMMARY (position error only -- these all use the wheel-only attitude control channel; attitude numbers are the same 87.56 deg baseline throughout and not the point of this comparison) ===")
order = ["raw_twobody", "raw_j2", "raw_harmonics", "sma_twobody", "sma_j2", "sma_harmonics"]
labels = ["Raw IC\ntwo-body", "Raw IC\n+J2", "Raw IC\n+50x50 EGM", "SMA IC\ntwo-body", "SMA IC\n+J2", "SMA IC\n+50x50 EGM"]
for k in order
    r = results[k]
    println(rpad(k, 16), " position mean=", round(r.pos_mean_km, digits=4), " km  max=", round(r.pos_max_km, digits=4), " km  final=", round(r.pos_final_km, digits=4), " km",
        "  |  velocity mean=", round(r.vel_mean_km_s, digits=6), " km/s  max=", round(r.vel_max_km_s, digits=6), " km/s  final=", round(r.vel_final_km_s, digits=6), " km/s")
end

pos_means_km = [results[k].pos_mean_km for k in order]
pos_maxes_km = [results[k].pos_max_km for k in order]
vel_means_km_s = [results[k].vel_mean_km_s for k in order]
vel_maxes_km_s = [results[k].vel_max_km_s for k in order]

ablation_out = DataFrame(
    key=order, label=labels,
    pos_mean_km=pos_means_km, pos_max_km=pos_maxes_km,
    vel_mean_km_s=vel_means_km_s, vel_max_km_s=vel_maxes_km_s,
)
ablation_path = joinpath(PLOT_DATA_DIR, "position_accuracy_ablation.arrow")
Arrow.write(ablation_path, ablation_out)
println("ablation plot data written to: $(ablation_path)")

# ==============================================================================
# Best-case position time series and component comparison
# ==============================================================================

best = results["sma_harmonics"]
sol = best.sol
n_samples = 600
sample_t = range(T0, TEND, length=n_samples)
t_s = collect(sample_t) .- T0
pos_sim = [SVector{3, Float64}(sol(t - T0).sc[1].pos) for t in sample_t]
pos_gt = [r_truth(t) for t in sample_t]
pos_err_km = [norm(pos_sim[i] - pos_gt[i]) / 1000.0 for i in eachindex(sample_t)]

vel_sim = [SVector{3, Float64}(sol(t - T0).sc[1].vel) for t in sample_t]
vel_gt = [v_truth(t) for t in sample_t]
vel_err_km_s = [norm(vel_sim[i] - vel_gt[i]) / 1000.0 for i in eachindex(sample_t)]

best_case_out = DataFrame(
    t_s=t_s,
    pos_x_sim=[p[1] / 1000.0 for p in pos_sim], pos_y_sim=[p[2] / 1000.0 for p in pos_sim], pos_z_sim=[p[3] / 1000.0 for p in pos_sim],
    pos_x_gt=[p[1] / 1000.0 for p in pos_gt], pos_y_gt=[p[2] / 1000.0 for p in pos_gt], pos_z_gt=[p[3] / 1000.0 for p in pos_gt],
    pos_err_km=pos_err_km,
    vel_x_sim=[v[1] / 1000.0 for v in vel_sim], vel_y_sim=[v[2] / 1000.0 for v in vel_sim], vel_z_sim=[v[3] / 1000.0 for v in vel_sim],
    vel_x_gt=[v[1] / 1000.0 for v in vel_gt], vel_y_gt=[v[2] / 1000.0 for v in vel_gt], vel_z_gt=[v[3] / 1000.0 for v in vel_gt],
    vel_err_km_s=vel_err_km_s,
)
best_case_path = joinpath(PLOT_DATA_DIR, "position_accuracy_best_case.arrow")
Arrow.write(best_case_path, best_case_out)
println("best-case plot data written to: $(best_case_path)")
