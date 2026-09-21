#!/usr/bin/env julia

const _PP_REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
import Pkg; Pkg.activate(_PP_REPO_ROOT; io=devnull)
# Reads a simulation_results.feather file and prints an orbital-element summary
# for the target satellite (sc1) in the same format as run_case2_laser_links.jl.
#
# Usage:
#   julia --project=. ORACLE/post_processing/print_oe_summary.jl PATH/TO/simulation_results.feather

using Arrow
using DataFrames
using LinearAlgebra
using Printf

const MU_EARTH = 3.986004418e14   # [m³/s²]

# ── Argument ──────────────────────────────────────────────────────────────────
isempty(ARGS) && error(
    "Usage: julia --project=. ORACLE/post_processing/print_oe_summary.jl <feather_path>"
)
feather_path = abspath(ARGS[1])
isfile(feather_path) || error("File not found: $feather_path")

# ── Read ──────────────────────────────────────────────────────────────────────
df = DataFrame(Arrow.Table(feather_path))
n  = nrow(df)

# ── Orbital elements from pos+vel ─────────────────────────────────────────────
function orbital_elements(px, py, pz, vx, vy, vz)
    r_vec = [px, py, pz]
    v_vec = [vx, vy, vz]
    r = norm(r_vec)
    v2 = dot(v_vec, v_vec)

    a    = -MU_EARTH / (v2 - 2.0 * MU_EARTH / r)          # semi-major axis [m]
    h    = cross(r_vec, v_vec)                              # angular momentum
    hmag = norm(h)
    e_vec = cross(v_vec, h) / MU_EARTH - r_vec / r         # eccentricity vector
    e    = norm(e_vec)                                      # eccentricity
    i    = acos(clamp(h[3] / hmag, -1.0, 1.0))             # inclination [rad]
    n_asc = cross([0.0, 0.0, 1.0], h)
    n_mag = norm(n_asc)
    raan  = n_mag < 1e-12 ? 0.0 : begin
        Ω = acos(clamp(n_asc[1] / n_mag, -1.0, 1.0))
        n_asc[2] < 0.0 ? 2π - Ω : Ω
    end
    return (a=a, e=e, i_deg=rad2deg(i), raan_deg=rad2deg(raan))
end

oe0 = orbital_elements(
    df.sc1_pos_1[1], df.sc1_pos_2[1], df.sc1_pos_3[1],
    df.sc1_vel_1[1], df.sc1_vel_2[1], df.sc1_vel_3[1],
)
oef = orbital_elements(
    df.sc1_pos_1[n], df.sc1_pos_2[n], df.sc1_pos_3[n],
    df.sc1_vel_1[n], df.sc1_vel_2[n], df.sc1_vel_3[n],
)

# ── Final accumulated ΔV ──────────────────────────────────────────────────────
dv_r = Float64(df.dv_r_accumulated[n])
dv_t = Float64(df.dv_t_accumulated[n])
dv_n = Float64(df.dv_n_accumulated[n])

# ── Active helper stats (if saved) ────────────────────────────────────────────
active_helper_col = hasproperty(df, :laser_active_helper) ?
    Float64.(df.laser_active_helper) : nothing

active_steps = active_helper_col !== nothing ?
    count(x -> x > 0, active_helper_col) : "N/A"

# Count switches (activations) — consecutive changes in active helper (ignoring 0→0)
activations = if active_helper_col !== nothing
    count(i -> active_helper_col[i] > 0 && active_helper_col[i] != active_helper_col[i-1],
          2:length(active_helper_col))
else
    "N/A"
end

# ── Orbit count ───────────────────────────────────────────────────────────────
t_s = Float64.(df.time)
T_series = [begin
    r = sqrt(df.sc1_pos_1[k]^2 + df.sc1_pos_2[k]^2 + df.sc1_pos_3[k]^2)
    v2 = df.sc1_vel_1[k]^2 + df.sc1_vel_2[k]^2 + df.sc1_vel_3[k]^2
    a_k = -MU_EARTH / (v2 - 2.0 * MU_EARTH / r)
    2π * sqrt(a_k^3 / MU_EARTH)
end for k in 1:n]
dt = diff(t_s)
T_mid = (T_series[1:end-1] .+ T_series[2:end]) ./ 2
orbits_elapsed = sum(dt ./ T_mid)

# ── Print ─────────────────────────────────────────────────────────────────────
println("ORACLE Open Cavity Case — feather summary (target satellite only)")
println("  file                  : $feather_path")
println("  saved timesteps       : $n")
@printf("  mission time          : %.2f s  (%.4f orbits elapsed)\n", t_s[end], orbits_elapsed)
println("  link activations      : $activations")
println("  accepted active steps : $active_steps")
println()
println("  ── Initial orbital elements ──────────────────────")
@printf("  a   [m]   : % .6e\n", oe0.a)
@printf("  e         : % .6e\n", oe0.e)
@printf("  i   [deg] : % .6e\n", oe0.i_deg)
@printf("  RAAN[deg] : % .6e\n", oe0.raan_deg)
println()
println("  ── Final orbital elements ────────────────────────")
@printf("  a   [m]   : % .6e\n", oef.a)
@printf("  e         : % .6e\n", oef.e)
@printf("  i   [deg] : % .6e\n", oef.i_deg)
@printf("  RAAN[deg] : % .6e\n", oef.raan_deg)
println()
println("  ── Changes (final − initial) ─────────────────────")
@printf("  da  [m]   : % .6e\n", oef.a       - oe0.a)
@printf("  de        : % .6e\n", oef.e       - oe0.e)
@printf("  di  [deg] : % .6e\n", oef.i_deg   - oe0.i_deg)
@printf("  dRAAN[deg]: % .6e\n", oef.raan_deg - oe0.raan_deg)
println()
println("  ── Accumulated laser ΔV [m/s] ───────────────────")
@printf("  dV_RTN    : R=% .6e  T=% .6e  N=% .6e\n", dv_r, dv_t, dv_n)
