#!/usr/bin/env julia
# Reads one scenario's feather file and plots the time history of cumulative
# ΔvR, ΔvT, ΔvN for that specific case.
#
# Run from the repository root:
#   julia --project=. ORACLE/post_processing/plot_dv_timeseries.jl [options]
#
# Options (all optional — defaults shown below):
#   --helper-altitude-km   1000.0
#   --target-altitude-km   1000.0
#   --helper-inclination   0.0
#   --target-inclination   0.0
#   --helpers              200
#   --orbits               80.0      
#   --schedule             naive_next_entering
#   --output-dir           output/paper_plot_mode

using Arrow
using DataFrames
get!(ENV, "GKSwstype", "100")   # headless rendering (remove for interactive window)
using Plots
using LaTeXStrings
using Printf

gr()

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# ── Argument parsing ──────────────────────────────────────────────────────────
function _parse(argv)
    opts = Dict{String,Any}(
        "helper-altitude-km"  => 1000.0,
        "target-altitude-km"  => 1000.0,
        "helper-inclination"  => 0.0,
        "target-inclination"  => 0.0,
        "helpers"             => 200,
        "orbits"              => 80.0,
        "schedule"            => "naive_next_entering",
        "output-dir"          => joinpath(REPO_ROOT, "output", "paper_plot_mode"),
    )
    i = 1
    while i <= length(argv)
        arg = argv[i]
        (arg == "--help" || arg == "-h") && (println("Pass --key value pairs matching the options above."); exit(0))
        startswith(arg, "--") || (i += 1; continue)
        key = arg[3:end]
        haskey(opts, key) || error("Unknown option: $arg")
        i < length(argv) || error("Missing value for $arg")
        v = argv[i+1]
        opts[key] = opts[key] isa Int    ? parse(Int,     v) :
                    opts[key] isa Float64 ? parse(Float64, v) : v
        i += 2
    end
    return opts
end

opts = _parse(ARGS)

helper_alt_km = opts["helper-altitude-km"]
target_alt_km = opts["target-altitude-km"]
helper_inc    = opts["helper-inclination"]
target_inc    = opts["target-inclination"]
n_helpers     = opts["helpers"]
orbits        = opts["orbits"]
schedule      = opts["schedule"]
output_dir    = opts["output-dir"]

# ── Locate the feather file ───────────────────────────────────────────────────
alt_folder = @sprintf("h%.0fkm_t%.0fkm",    helper_alt_km, target_alt_km)
inc_folder = @sprintf("ih%.1fdeg_it%.1fdeg", helper_inc,    target_inc)
n_folder   = @sprintf("N%d",                 n_helpers)

base_path = joinpath(output_dir, alt_folder, inc_folder, n_folder)
isdir(base_path) || error("Scenario folder not found:\n  $base_path")

# Find the T*s folder dynamically (period depends on target altitude + orbits)
t_folders = filter(d -> startswith(d, "T") && endswith(d, "s"), readdir(base_path))
isempty(t_folders) && error("No T*s folder found inside:\n  $base_path")

feather_path = joinpath(base_path, first(t_folders), schedule, "simulation_results.feather")
isfile(feather_path) || error(
    "Feather file not found:\n  $feather_path\n" *
    "Check --schedule, --helper-altitude-km, --target-altitude-km, --helpers, etc."
)

println("Reading: $feather_path")
df = DataFrame(Arrow.Table(feather_path))

# ── Extract ΔV columns ───────────────────────────────────────────────────────
t_s  = Float64.(df.time)
dv_r = Float64.(df.dv_r_accumulated)
dv_t = Float64.(df.dv_t_accumulated)
dv_n = Float64.(df.dv_n_accumulated)

# Active helper index (if saved)
active_col = hasproperty(df, :laser_active_helper) ? Float64.(df.laser_active_helper) : nothing

# ── Compute cumulative orbit count from target's instantaneous semi-major axis ─
const MU_EARTH = 3.986004418e14   # [m³/s²]

function _semi_major(px, py, pz, vx, vy, vz)
    r = sqrt(px^2 + py^2 + pz^2)
    v2 = vx^2 + vy^2 + vz^2
    return -MU_EARTH / (v2 - 2*MU_EARTH/r)   # vis-viva: a = -μ/(v²-2μ/r)
end

T_series = [
    2pi * sqrt(_semi_major(
        df.sc1_pos_1[k], df.sc1_pos_2[k], df.sc1_pos_3[k],
        df.sc1_vel_1[k], df.sc1_vel_2[k], df.sc1_vel_3[k],
    )^3 / MU_EARTH)
    for k in eachindex(t_s)
]
dt       = diff(t_s)
T_mid    = (T_series[1:end-1] .+ T_series[2:end]) ./ 2
n_orbits = cumsum([0.0; dt ./ T_mid])   # cumulative orbit count

# ── Build the title string ────────────────────────────────────────────────────
title_str = @sprintf(
    "h_helper=%.0f km  h_target=%.0f km  Δi_h=%.1f°  Δi_t=%.1f°  N=%d  sched=%s",
    helper_alt_km, target_alt_km, helper_inc, target_inc, n_helpers, schedule
)

# ── Plot ──────────────────────────────────────────────────────────────────────
n_panels = active_col !== nothing ? 4 : 3

p1 = plot(n_orbits, dv_r;
    label=false, color=:blue, lw=1.5,
    ylabel=L"\Delta v_R\ \mathrm{(m/s)}",
    title=title_str, titlefontsize=9,
    xlabel="", xformatter=_->""
)
hline!(p1, [0.0]; color=:black, lw=0.6, linestyle=:dash, label=false)

p2 = plot(n_orbits, dv_t;
    label=false, color=:red, lw=1.5,
    ylabel=L"\Delta v_T\ \mathrm{(m/s)}",
    xlabel="", xformatter=_->""
)
hline!(p2, [0.0]; color=:black, lw=0.6, linestyle=:dash, label=false)

p3 = plot(n_orbits, dv_n;
    label=false, color=:green, lw=1.5,
    ylabel=L"\Delta v_N\ \mathrm{(m/s)}",
    xlabel=active_col !== nothing ? "" : L"\mathrm{Orbits}",
    xformatter=active_col !== nothing ? (_->"") : :auto
)
hline!(p3, [0.0]; color=:black, lw=0.6, linestyle=:dash, label=false)

panels = [p1, p2, p3]

if active_col !== nothing
    p4 = plot(n_orbits, active_col;
        label=false, color=:purple, lw=1.0,
        ylabel="Active helper",
        xlabel=L"\mathrm{Orbits}"
    )
    push!(panels, p4)
end

fig = plot(panels...;
    layout=(n_panels, 1),
    size=(800, 200 * n_panels + 60),
    left_margin=8Plots.mm,
    bottom_margin=4Plots.mm,
    link=:x,
)

# ── Save ──────────────────────────────────────────────────────────────────────
img_dir = joinpath(REPO_ROOT, "output", "images", "delta_v_timeseries")
mkpath(img_dir)

stem = @sprintf("dv_ts_h%.0f_t%.0f_ih%.1f_it%.1f_N%d_%s",
    helper_alt_km, target_alt_km, helper_inc, target_inc, n_helpers, schedule)

out_png = joinpath(img_dir, stem * ".png")
out_pdf = joinpath(img_dir, stem * ".pdf")
savefig(fig, out_png)
savefig(fig, out_pdf)
display(fig)
println("Saved → ", out_png)
println("Saved → ", out_pdf)
