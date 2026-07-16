##
# Generates a single combined LaTeX table of RMSE metrics for the CYGNSS
# study, from the data written by run_position_accuracy_48hr.jl and
# run_kinematic_backout.jl to data/plot_data/*_rmse.arrow. Doesn't re-run
# any simulation.
#
# Position and velocity RMSE come from the 48-hour reference-replication
# case (the validated TV.run_verification pipeline, run over the full
# 48hr PVT telemetry file) -- the longer window is the meaningful test of
# translational accuracy. Quaternion and angular-velocity RMSE come from
# the kinematic torque back-out case (the full-hour slew maneuver window),
# since that is where attitude reconstruction is actually exercised --
# the 48hr case and the other 1hr studies use a wheel-only attitude
# baseline that isn't representative (see run_independent_effects.jl and
# run_position_accuracy.jl headers).
#
# All four metrics are magnitude (norm) errors, not per-axis/elementwise:
#   - position/velocity RMSE: mean of the per-axis (x/y/z) RMSEs from the
#     TV.run_verification pipeline (same methodology as the ~1.58 km
#     reference figure quoted throughout this study's README)
#   - quaternion RMSE: sqrt(mean(angle^2)), angle = 2*acos(|q_sim . q_gt|)
#   - angular velocity RMSE: sqrt(mean(norm(w_sim - w_gt)^2))
#
# Usage: julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/make_rmse_table.jl
##

using Arrow
using DataFrames
using Printf

const STUDY_DIR = @__DIR__
const PLOT_DATA_DIR = joinpath(STUDY_DIR, "data", "plot_data")
const TABLES_DIR = joinpath(STUDY_DIR, "tables")
mkpath(TABLES_DIR)

_load(name) = DataFrame(Arrow.Table(joinpath(PLOT_DATA_DIR, name)))

position_velocity = _load("position_accuracy_48hr_rmse.arrow")
attitude = _load("kinematic_backout_rmse.arrow")

# LaTeX-formatted scientific notation, e.g. 1.489e-06 -> "1.489 \times 10^{-6}"
function _sci_tex(x::Float64; digits::Int=3)::String
    s = @sprintf("%.*e", digits, x)
    mantissa, exponent = split(s, "e")
    return "\$$(mantissa) \\times 10^{$(parse(Int, exponent))}\$"
end

io = IOBuffer()
println(io, "\\begin{table}[htbp]")
println(io, "\\centering")
println(io, "\\caption{CYGNSS reconstruction RMSE vs.\\ telemetry. Position and velocity are from the 48-hour reference-replication case; quaternion and angular velocity are from the full-hour slew-maneuver kinematic torque back-out case.}")
println(io, "\\label{tab:cygnss_rmse_summary}")
println(io, "\\begin{tabular}{rrrr}")
println(io, "\\toprule")
println(io, "Position RMSE (km) & Velocity RMSE (km/s) & Quaternion RMSE (deg) & Angular velocity RMSE (deg/s) \\\\")
println(io, "\\midrule")
@printf(
    io, "%.3f & %.5f & %.3f & %s \\\\\n",
    position_velocity.rmse_pos_km[1], position_velocity.rmse_vel_km_s[1],
    attitude.rmse_angle_deg[1], _sci_tex(attitude.rmse_ang_vel_deg_s[1]),
)
println(io, "\\bottomrule")
println(io, "\\end{tabular}")
println(io, "\\end{table}")

tex_path = joinpath(TABLES_DIR, "cygnss_rmse_tables.tex")
write(tex_path, String(take!(io)))
println("RMSE LaTeX table written to: $(tex_path)")
