# Plotting helpers for oracle_laser_links.jl results.
# All plots are saved as PNG to IMG_DIR.
# Requires: Plots.jl, Arrow.jl, DataFrames.jl, LinearAlgebra.

import Plots
using Arrow, DataFrames
using LinearAlgebra
using Printf

# Gravitational parameter for semi-major axis calculation
const _MU_EARTH_PLOT = 3.986004418e14

function _sma_from_rv(r1, r2, r3, v1, v2, v3)
    r = sqrt(r1^2 + r2^2 + r3^2)
    v2 = v1^2 + v2^2 + v3^2
    return -_MU_EARTH_PLOT / (v2 - 2.0 * _MU_EARTH_PLOT / r)
end

"""
    plot_oracle_laser_results(feather_path, tracker, results_dir; n_helpers)

Generates 4 diagnostic plots from a completed oracle_laser_links run and saves them
to `results_dir/images/`. Requires the feather path, the LaserImpulseTracker, and
the number of helper satellites.
"""
function plot_oracle_laser_results(
    feather_path::String,
    tracker,
    results_dir::String;
    n_helpers::Int = 10,
)
    img_dir = joinpath(results_dir, "images")
    mkpath(img_dir)

    df = DataFrame(Arrow.Table(feather_path))
    t  = df.time

    # ── Plot 1: Cumulative ΔV (R, T, N) vs time, one series per active link ─
    p1 = Plots.plot(
        xlabel="time (s)", ylabel="cumulative Δv (m/s)",
        title="Laser Δv in RTN frame (per link)"
    )
    for link in keys(tracker.dv_R_hist)
        label_suffix = " $(link)"
        Plots.plot!(p1, tracker.t_hist, tracker.dv_R_hist[link]; label="ΔvR$label_suffix", lw=2)
        Plots.plot!(p1, tracker.t_hist, tracker.dv_T_hist[link]; label="ΔvT$label_suffix", lw=2)
        Plots.plot!(p1, tracker.t_hist, tracker.dv_N_hist[link]; label="ΔvN$label_suffix", lw=2)
    end
    Plots.savefig(p1, joinpath(img_dir, "dv_RTN.png"))

    # ── Plot 2: Target altitude vs time ────────────────────────────────────
    p2 = Plots.plot(
        t, df.sc1_altitude ./ 1e3;
        label="target", xlabel="time (s)", ylabel="altitude (km)",
        title="Target altitude", lw=2, color=:blue
    )
    Plots.savefig(p2, joinpath(img_dir, "altitude_vs_time.png"))

    # ── Plot 3: Target semi-major axis change Δa vs time ───────────────────
    a0 = _sma_from_rv(df.sc1_pos_1[1], df.sc1_pos_2[1], df.sc1_pos_3[1],
                      df.sc1_vel_1[1], df.sc1_vel_2[1], df.sc1_vel_3[1])
    Δa = [_sma_from_rv(df.sc1_pos_1[k], df.sc1_pos_2[k], df.sc1_pos_3[k],
                       df.sc1_vel_1[k], df.sc1_vel_2[k], df.sc1_vel_3[k]) - a0
          for k in eachindex(t)]
    p3 = Plots.plot(
        t, Δa;
        label="Δa", xlabel="time (s)", ylabel="Δa (m)",
        title="Target semi-major axis change", lw=2, color=:red
    )
    Plots.savefig(p3, joinpath(img_dir, "delta_sma.png"))

    # ── Plot 4: Satellite orbits in ECI XY plane ────────────────────────────
    θ = range(0, 2π, 200)
    R_e = 6_378_137.0
    p4 = Plots.plot(
        R_e .* cos.(θ) ./ 1e6, R_e .* sin.(θ) ./ 1e6;
        label="Earth", color=:blue, lw=1,
        xlabel="ECI x (Mm)", ylabel="ECI y (Mm)",
        title="Satellite orbits (ECI XY)", aspect_ratio=1
    )
    Plots.plot!(p4, df.sc1_pos_1 ./ 1e6, df.sc1_pos_2 ./ 1e6;
                label="target", lw=2, color=:red)
    for k in 2:(n_helpers + 1)
        xcol = Symbol("sc$(k)_pos_1")
        ycol = Symbol("sc$(k)_pos_2")
        Plots.plot!(p4, df[!, xcol] ./ 1e6, df[!, ycol] ./ 1e6;
                    label=(k == 2 ? "helpers" : ""), lw=1, color=:gray, alpha=0.5)
    end
    Plots.savefig(p4, joinpath(img_dir, "orbits_XY.png"))

    println("  Plots saved to: $img_dir")
    return (dv_rtn=p1, altitude=p2, delta_sma=p3, orbits=p4)
end
