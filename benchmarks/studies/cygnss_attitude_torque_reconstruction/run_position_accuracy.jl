##
# Position-accuracy replication: applies the exact IC-noise fix from
# test/gmat_scenario_matrix.jl's _build_cygnss_48hr_reference (N=5 vis-viva
# SMA averaging + velocity-magnitude rescale, no direction change) and that
# reference's full 50x50 EGM gravity field (not just J2) to this study's
# telemetry window, to see how close a full-hour CYGNSS orbit propagation
# gets to that reference's ~1.2-1.6 km/48hr position RMSE. See README.md
# ("Position accuracy replication") for the ablation table and result.
#
# Usage: julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_position_accuracy.jl
##

include(joinpath(@__DIR__, "common.jl"))
using Plots

Plots.default(left_margin=10Plots.mm, bottom_margin=8Plots.mm)

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
labels = ["Raw IC\ntwo-body", "Raw IC\n+J2", "Raw IC\n+50x50 EGM", "SMA-corrected IC\ntwo-body", "SMA-corrected IC\n+J2", "SMA-corrected IC\n+50x50 EGM"]
for k in order
    r = results[k]
    println(rpad(k, 16), " position mean=", round(r.pos_mean_km, digits=4), " km  max=", round(r.pos_max_km, digits=4), " km  final=", round(r.pos_final_km, digits=4), " km")
end

pos_means_km = [results[k].pos_mean_km for k in order]
pos_maxes_km = [results[k].pos_max_km for k in order]
xs = 1:length(order)
bar_fig = bar(
    xs .- 0.15, pos_means_km;
    bar_width=0.3, label="mean error", color=:steelblue,
    title="CYGNSS Full-Hour Position Error: IC Correction x Gravity Fidelity Ablation",
    ylabel="Position error (km)", xlabel="",
    xticks=(xs, labels), size=(1100, 550), titlefontsize=11, legend=:topright,
    xtickfontsize=8,
)
bar!(bar_fig, xs .+ 0.15, pos_maxes_km; bar_width=0.3, label="max error", color=:darkorange)
bar_path = joinpath(PLOTS_DIR, "cygnss_position_accuracy_ablation.png")
savefig(bar_fig, bar_path)
println("ablation bar chart: $(bar_path)")

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

pos_names = ("x", "y", "z")
pos_subplots = map(1:3) do i
    sp = plot(
        t_s, [p[i] / 1000.0 for p in pos_sim];
        label="simulated (SMA-corrected IC, 50x50 EGM)", lw=1.4, color=:steelblue,
        title=pos_names[i], xlabel="Time since t=0s (s)", ylabel="km",
        legend=(i == 1 ? :best : false),
    )
    plot!(sp, t_s, [p[i] / 1000.0 for p in pos_gt]; label="telemetry", lw=1.0, ls=:dash, color=:darkorange)
    return sp
end
pos_fig = plot(pos_subplots...; layout=(3, 1), size=(1000, 800), plot_title="CYGNSS Full Hour: ECI Position, Best-Case Replication vs. Telemetry")
pos_path = joinpath(PLOTS_DIR, "cygnss_position_accuracy_best_case_timeseries.png")
savefig(pos_fig, pos_path)
println("best-case position time series plot: $(pos_path)")

pos_err_fig = plot(
    t_s, pos_err_km;
    label="position error (SMA-corrected IC, 50x50 EGM)", lw=1.6, color=:steelblue,
    title="CYGNSS Full-Hour Position Error: Best-Case Replication",
    xlabel="Time since t=0s (s)", ylabel="Position error (km)",
    legend=:topleft, size=(1000, 500), titlefontsize=11,
)
pos_err_path = joinpath(PLOTS_DIR, "cygnss_position_accuracy_best_case_error_timeseries.png")
savefig(pos_err_fig, pos_err_path)
println("best-case position error time series plot: $(pos_err_path)")
