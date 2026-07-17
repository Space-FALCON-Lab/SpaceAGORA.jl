##
# Full-hour, independent-effect quantification: gravity-gradient, GRAM
# aerodynamic torque, and real dplCmd-telemetry-driven magnetic torque, each
# added on top of the Omega_rw/J_RW wheel-telemetry replay baseline, then all
# combined. See README.md ("Independent-effect results") for the numbers and
# interpretation.
#
# Usage: julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_independent_effects.jl
##

include(joinpath(@__DIR__, "common.jl"))
using Plots

Plots.default(left_margin=10Plots.mm, bottom_margin=14Plots.mm)

const T0 = 0.0
const TEND = 3600.0

results = Dict{String, Any}()

println("=== FULL HOUR: wheel-only baseline ===")
results["wheel_only"] = run_cygnss_case(
    t_cal=T0, t_end=TEND,
    dynamic_effectors_factory=(t_offset) -> (InverseSquaredGravityModel(gravity_gradient=false),),
    label="wheel_only", dt_max=1.0,
)

println("=== FULL HOUR: + gravity-gradient only ===")
results["gg"] = run_cygnss_case(
    t_cal=T0, t_end=TEND,
    dynamic_effectors_factory=(t_offset) -> (InverseSquaredGravityModel(gravity_gradient=true),),
    label="gravity_gradient", dt_max=1.0,
)

println("=== FULL HOUR: + aero (GRAM, free-molecular, per-link density) only ===")
results["aero"] = run_cygnss_case(
    t_cal=T0, t_end=TEND,
    dynamic_effectors_factory=(t_offset) -> (InverseSquaredGravityModel(gravity_gradient=false), AerodynamicCoefficientfM()),
    label="aero_gram", dt_max=1.0,
    density_model=GRAMAtmosphereModel(planet_name="earth"),
)

println("=== FULL HOUR: + magnetic (telemetry dplCmd) only ===")
results["magnetic"] = run_cygnss_case(
    t_cal=T0, t_end=TEND,
    dynamic_effectors_factory=(t_offset) -> (InverseSquaredGravityModel(gravity_gradient=false), TelemetryMagnetorquerModel(m_body_itp, t_offset)),
    label="magnetic_telemetry", dt_max=1.0,
)

println("=== FULL HOUR: full combination (GG + aero + magnetic) ===")
results["full"] = run_cygnss_case(
    t_cal=T0, t_end=TEND,
    dynamic_effectors_factory=(t_offset) -> (InverseSquaredGravityModel(gravity_gradient=true), AerodynamicCoefficientfM(), TelemetryMagnetorquerModel(m_body_itp, t_offset)),
    label="full_combo", dt_max=1.0,
    density_model=GRAMAtmosphereModel(planet_name="earth"),
)

println()
println("=== SUMMARY ===")
order = ["wheel_only", "gg", "aero", "magnetic", "full"]
labels = ["Wheel-only", "+ Gravity-gradient", "+ Aero (GRAM)", "+ Magnetic (dplCmd)", "All combined"]
for k in order
    r = results[k]
    println(
        rpad(k, 15),
        " attitude mean=", round(r.mean, digits=3), " deg  max=", round(r.max, digits=3), " deg  final=", round(r.final, digits=3), " deg",
        "  |  position mean=", round(r.pos_mean_km, digits=4), " km  max=", round(r.pos_max_km, digits=4), " km  final=", round(r.pos_final_km, digits=4), " km",
    )
end

means = [results[k].mean for k in order]
maxes = [results[k].max for k in order]

# Plain Plots.jl has no native `groupedbar` (that's a StatsPlots recipe, not a
# dependency here) -- dodge two bar series manually with offset x positions.
xs = 1:length(order)
bar_fig = bar(
    xs .- 0.15, means;
    bar_width=0.3, label="mean error", color=:steelblue,
    title="CYGNSS Full-Hour Attitude Reconstruction Error by Effect (Omega_rw/J_RW Wheel Baseline)",
    ylabel="Quaternion angle error (deg)", xlabel="",
    xticks=(xs, labels), size=(1000, 550), titlefontsize=11, legend=:topright,
    xrotation=15,
)
bar!(bar_fig, xs .+ 0.15, maxes; bar_width=0.3, label="max error", color=:darkorange)
bar_path = joinpath(PLOTS_DIR, "cygnss_independent_effects_summary.png")
savefig(bar_fig, bar_path)
println("attitude summary bar chart: $(bar_path)")

# ==============================================================================
# Position comparison
#
# Different effect combinations propagate different translational dynamics
# too, not just torque: the aero cases add real drag force
# (AerodynamicCoefficientfM's force_ii), so position drift differs between
# the wheel-only/gravity-gradient-only cases (gravity force only, no drag)
# and the aero/full-combo cases.
# ==============================================================================

pos_means_km = [results[k].pos_mean_km for k in order]
pos_maxes_km = [results[k].pos_max_km for k in order]

pos_bar_fig = bar(
    xs .- 0.15, pos_means_km;
    bar_width=0.3, label="mean error", color=:steelblue,
    title="CYGNSS Full-Hour Position Reconstruction Error by Effect",
    ylabel="Position error (km)", xlabel="",
    xticks=(xs, labels), size=(1000, 550), titlefontsize=11, legend=:topright,
    xrotation=15,
)
bar!(pos_bar_fig, xs .+ 0.15, pos_maxes_km; bar_width=0.3, label="max error", color=:darkorange)
pos_bar_path = joinpath(PLOTS_DIR, "cygnss_independent_effects_position_summary.png")
savefig(pos_bar_fig, pos_bar_path)
println("position summary bar chart: $(pos_bar_path)")

sample_t = range(T0, TEND, length=300)
t_s = collect(sample_t) .- T0
pos_ts_fig = plot(
    title="CYGNSS Full-Hour Position Error vs. Time by Effect",
    xlabel="Time since t=0s (s)", ylabel="Position error (km)",
    legend=:topleft, size=(1000, 550), titlefontsize=11,
)
colors = (:steelblue, :darkorange, :seagreen, :purple, :firebrick)
for (k, lbl, c) in zip(order, labels, colors)
    plot!(pos_ts_fig, t_s, results[k].pos_errs_km; label=lbl, lw=1.4, color=c)
end
pos_ts_path = joinpath(PLOTS_DIR, "cygnss_independent_effects_position_timeseries.png")
savefig(pos_ts_fig, pos_ts_path)
println("position error time series plot: $(pos_ts_path)")
