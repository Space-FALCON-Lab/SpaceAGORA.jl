##
# Full-hour, independent-effect quantification: gravity-gradient, GRAM
# aerodynamic torque, and real dplCmd-telemetry-driven magnetic torque, each
# added on top of the Omega_rw/J_RW wheel-telemetry replay baseline, then all
# combined. See README.md ("Independent-effect results") for the numbers and
# interpretation.
#
# This script only runs the simulations and writes the plotted quantities to
# data/plot_data/independent_effects_summary.arrow and
# data/plot_data/independent_effects_position_timeseries.arrow -- run
# make_plots.jl separately to generate the figures without re-running the
# simulations.
#
# Usage: julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_independent_effects.jl
##

include(joinpath(@__DIR__, "common.jl"))

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

# Per-link atmosphere sampling is opt-in (defaults off so baseline simulation
# physics matches the single-sample behavior); the aero cases here want it.
SimulationModel.DynamicEffectors.AerodynamicEffectors.set_per_link_atmosphere!(true)

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

SimulationModel.DynamicEffectors.AerodynamicEffectors.set_per_link_atmosphere!(false)

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
        "  |  velocity mean=", round(r.vel_mean_km_s, digits=6), " km/s  max=", round(r.vel_max_km_s, digits=6), " km/s  final=", round(r.vel_final_km_s, digits=6), " km/s",
    )
end

means = [results[k].mean for k in order]
maxes = [results[k].max for k in order]

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
vel_means_km_s = [results[k].vel_mean_km_s for k in order]
vel_maxes_km_s = [results[k].vel_max_km_s for k in order]

summary_out = DataFrame(
    key=order, label=labels,
    attitude_mean_deg=means, attitude_max_deg=maxes,
    pos_mean_km=pos_means_km, pos_max_km=pos_maxes_km,
    vel_mean_km_s=vel_means_km_s, vel_max_km_s=vel_maxes_km_s,
)
summary_path = joinpath(PLOT_DATA_DIR, "independent_effects_summary.arrow")
Arrow.write(summary_path, summary_out)
println("summary plot data written to: $(summary_path)")

sample_t = range(T0, TEND, length=300)
t_s = collect(sample_t) .- T0
ts_out = DataFrame(t_s=t_s)
for k in order
    ts_out[!, Symbol(k)] = results[k].pos_errs_km
end
ts_path = joinpath(PLOT_DATA_DIR, "independent_effects_position_timeseries.arrow")
Arrow.write(ts_path, ts_out)
println("position error time series plot data written to: $(ts_path)")

vel_ts_out = DataFrame(t_s=t_s)
for k in order
    vel_ts_out[!, Symbol(k)] = results[k].vel_errs_km_s
end
vel_ts_path = joinpath(PLOT_DATA_DIR, "independent_effects_velocity_timeseries.arrow")
Arrow.write(vel_ts_path, vel_ts_out)
println("velocity error time series plot data written to: $(vel_ts_path)")
