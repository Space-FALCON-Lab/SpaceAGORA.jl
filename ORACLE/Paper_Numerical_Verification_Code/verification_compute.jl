# =============================================================================
# Section IV.F — Validation Diagnostics (Compute Layer)
# =============================================================================
# Functions that read a Feather file produced by run_verification.jl and
# compute all conservation-law residuals described in Sec. IV.F:
#
#   1. Net OCL force:                F_net = Σ F_i^(OCL)                 [N]
#   2. Net OCL torque:               τ_net = Σ r_i × F_i^(OCL)          [N·m]
#   3. OCL angular-momentum exchange:ε_H   = ∫(r_t×F_{t←h} + r_h×F_{h←t})dt [kg·m²/s]
#   4. ΔE_orb:                       total orbital mechanical energy change [J]
#   5. W_OCL:                        accumulated OCL mechanical work       [J]
#   6. ε_E = ΔE_orb − W_OCL:        energy residual (should ≈ 0)          [J]
#   7. H_z:                          z-component of total angular momentum  [kg·m²/s]
#   8. Orbit count (non-dimensional time axis)
#
# Equations referenced: paper Eqs. (13)–(23)
# =============================================================================

using Arrow
using DataFrames
using LinearAlgebra
using StaticArrays
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# Physical constants
# ─────────────────────────────────────────────────────────────────────────────
const _C_LIGHT = 299_792_458.0   # speed of light [m/s]

# ─────────────────────────────────────────────────────────────────────────────
# Laser parameters container (mirrors Table 2 defaults)
# ─────────────────────────────────────────────────────────────────────────────
struct LaserParams
    range_m::Float64
    power_w::Float64
    magnification::Float64
    beta::Float64
    eta::Float64
end

LaserParams(;
    range_km=200.0, power_w=10_000.0,
    magnification=100.0, beta=1.0, eta=1.0
) = LaserParams(range_km * 1e3, power_w, magnification, beta, eta)

# ─────────────────────────────────────────────────────────────────────────────
# Gravitational potential consistent with SpaceAGORA's gravity model
# ─────────────────────────────────────────────────────────────────────────────
"""
    grav_potential(r::SVector{3}, mu, RE, J2, use_j2) → Φ [m²/s²]

Gravitational potential per unit mass at position `r` [m].
When `use_j2=true`, includes the J2 zonal harmonic term:
    Φ = −μ/r + μ J2 R_E² (3z² − r²) / (2r⁵)
This matches the `InverseSquaredJ2GravityModel` used in SpaceAGORA (Eq. 2).
"""
function grav_potential(
    r::SVector{3,Float64},
    mu::Float64,
    RE::Float64,
    J2::Float64,
    use_j2::Bool,
)::Float64
    rn = norm(r)
    phi = -mu / rn
    if use_j2
        z   = r[3]
        phi += mu * J2 * RE^2 * (3z^2 - rn^2) / (2 * rn^5)
    end
    return phi
end

# ─────────────────────────────────────────────────────────────────────────────
# OCL force recomputed from saved positions
# ─────────────────────────────────────────────────────────────────────────────
"""
    ocl_force_on_target(r_t, r_h, lp, is_active) → SVector{3} [N]

OCL force on the target due to the helper, using Eq. (4).
Returns the zero vector when `is_active = false` or when the pair is out of range.
"""
function ocl_force_on_target(
    r_t::SVector{3,Float64},
    r_h::SVector{3,Float64},
    lp::LaserParams,
    is_active::Bool,
)::SVector{3,Float64}
    is_active || return @SVector[0.0, 0.0, 0.0]
    rel = r_t - r_h             # from helper to target
    rho = norm(rel)
    (rho <= 0.0 || rho > lp.range_m) && return @SVector[0.0, 0.0, 0.0]
    F_mag = lp.eta * lp.beta * lp.magnification * lp.power_w / _C_LIGHT   # [N]
    return F_mag * rel / rho
end

# ─────────────────────────────────────────────────────────────────────────────
# Orbit count from flat position/time arrays  (mirrors _orbit_count_from_flat_sol)
# ─────────────────────────────────────────────────────────────────────────────
"""
    compute_orbit_counts(t, r_target_cols, mu) → Vector{Float64}

Compute non-dimensional orbit count for the target satellite at each saved
time step. Uses instantaneous semi-major axis at each step.
"""
function compute_orbit_counts(
    t::AbstractVector{Float64},
    r_cols::NTuple{3, <:AbstractVector},   # (x, y, z) column arrays
    v_cols::NTuple{3, <:AbstractVector},
    mu::Float64,
)::Vector{Float64}
    n = length(t)
    T_series = Vector{Float64}(undef, n)
    for k in 1:n
        r = SVector{3,Float64}(r_cols[1][k], r_cols[2][k], r_cols[3][k])
        v = SVector{3,Float64}(v_cols[1][k], v_cols[2][k], v_cols[3][k])
        a = norm(r)^2 / (2norm(r) - dot(v, v) * norm(r) / mu)  # vis-viva: a = μr/(2μ - v²r)
        # Simpler: a = −μ / (2ε),  ε = v²/2 − μ/r
        eps = 0.5 * dot(v, v) - mu / norm(r)
        a_safe = -mu / (2 * eps)
        T_series[k] = 2π * sqrt(abs(a_safe)^3 / mu)
    end
    dt     = diff(t)
    T_mid  = (T_series[1:end-1] .+ T_series[2:end]) ./ 2
    return cumsum([0.0; dt ./ T_mid])
end

# ─────────────────────────────────────────────────────────────────────────────
# Main diagnostic computation
# ─────────────────────────────────────────────────────────────────────────────
"""
    compute_validation_diagnostics(df, lp, mu, RE, J2, mass_kg, use_j2)
        → NamedTuple

Compute all Section IV.F validation residuals from a Feather DataFrame `df`.

**Returns** a NamedTuple with fields:
- `t`            : time vector [s]
- `orbit_counts` : non-dimensional orbit counter
- `F_net_mag`    : magnitude of net OCL force |Σ F_i^(OCL)| [N]   (Eq. 21)
- `tau_net_mag`  : magnitude of net OCL torque |Σ r_i×F_i^(OCL)| [N·m] (Eq. 21)
- `delta_Eorb`   : ΔE_orb = E_orb(t) − E_orb(0) [J]              (Eq. 13)
- `W_ocl`        : accumulated OCL work W_OCL(t) [J]               (Eq. 14)
- `eps_E`        : energy residual ε_E = ΔE_orb − W_OCL [J]       (Eq. 22)
- `eps_H_mag`    : |ε_H^(OCL)| = |∫(r_t×F_{t←h} + r_h×F_{h←t})dt| [kg·m²/s] (Eq. 23)
- `H_z`          : z-component of total angular momentum [kg·m²/s] (Eq. 17)
- `laser_active` : 0/1 link-active indicator
"""
function compute_validation_diagnostics(
    df::DataFrame,
    lp::LaserParams,
    mu::Float64,
    RE::Float64,
    J2::Float64,
    mass_kg::Float64,
    use_j2::Bool,
)
    t     = Float64.(df.time)
    npts  = length(t)

    # ── Auto-detect total spacecraft count from feather column names ──────────
    # SpaceAGORA writes sc1_pos_1 … sc{N}_pos_1 for each spacecraft.
    n_sc = count(c -> occursin(r"^sc\d+_pos_1$", c), names(df))
    n_sc >= 2 || error("Expected at least 2 spacecraft columns in feather (found $n_sc).")
    n_helpers = n_sc - 1   # spacecraft 1 = target, 2..N = helpers

    # ── Extract ECI states for ALL spacecraft [target + helpers] ──────────────
    # all_r[i][k] = SVector{3} position of spacecraft i at saved time k
    # all_v[i][k] = SVector{3} velocity of spacecraft i at saved time k
    all_r_x = [Float64.(df[!, "sc$(i)_pos_1"]) for i in 1:n_sc]
    all_r_y = [Float64.(df[!, "sc$(i)_pos_2"]) for i in 1:n_sc]
    all_r_z = [Float64.(df[!, "sc$(i)_pos_3"]) for i in 1:n_sc]
    all_v_x = [Float64.(df[!, "sc$(i)_vel_1"]) for i in 1:n_sc]
    all_v_y = [Float64.(df[!, "sc$(i)_vel_2"]) for i in 1:n_sc]
    all_v_z = [Float64.(df[!, "sc$(i)_vel_3"]) for i in 1:n_sc]

    # Convenience accessors for target (sc1)
    r_t_x, r_t_y, r_t_z = all_r_x[1], all_r_y[1], all_r_z[1]
    v_t_x, v_t_y, v_t_z = all_v_x[1], all_v_y[1], all_v_z[1]

    # ── Link-active indicator ─────────────────────────────────────────────────
    # laser_active_helper: 0 = no active link, >0 = index of active helper sc
    laser_active_col = ("laser_active_helper" ∈ names(df)) ?
                       Int.(df[!, "laser_active_helper"]) :
                       zeros(Int, npts)
    laser_active_bool = laser_active_col .> 0

    # ── Orbit count ───────────────────────────────────────────────────────────
    orbit_counts = compute_orbit_counts(
        t,
        (r_t_x, r_t_y, r_t_z),
        (v_t_x, v_t_y, v_t_z),
        mu,
    )

    # ── Per-timestep quantities ───────────────────────────────────────────────
    F_net_mag    = Vector{Float64}(undef, npts)
    tau_net_mag  = Vector{Float64}(undef, npts)
    delta_Eorb   = Vector{Float64}(undef, npts)
    H_z          = Vector{Float64}(undef, npts)
    # Per-satellite H_z: rows = spacecraft index, cols = time steps
    H_z_per_sat  = Matrix{Float64}(undef, n_sc, npts)

    E0 = 0.0   # will be set at k=1

    for k in 1:npts
        r_t = SVector{3,Float64}(r_t_x[k], r_t_y[k], r_t_z[k])
        v_t = SVector{3,Float64}(v_t_x[k], v_t_y[k], v_t_z[k])

        # ── Active helper position (for OCL force / torque) ───────────────────
        active_idx = laser_active_col[k]   # sc index of active helper (0 = none)
        is_active  = active_idx > 0
        r_active_h = is_active ?
            SVector{3,Float64}(all_r_x[active_idx][k], all_r_y[active_idx][k], all_r_z[active_idx][k]) :
            @SVector[0.0, 0.0, 0.0]

        # ── Recomputed OCL force (target ↔ active helper only) ───────────────
        F_t = ocl_force_on_target(r_t, r_active_h, lp, is_active)
        F_h = -F_t                                # Newton's 3rd law (Eq. 5)

        # Net OCL force (Eq. 21) — should be ≈ 0 at machine precision
        F_net = F_t + F_h
        F_net_mag[k] = norm(F_net)

        # Net OCL torque (Eq. 21) — should be ≈ 0 (forces collinear with displacement)
        tau = cross(r_t, F_t) + cross(r_active_h, F_h)
        tau_net_mag[k] = norm(tau)

        # ── Total orbital mechanical energy over ALL spacecraft (Eq. 13) ──────
        E_tot = 0.0
        for i in 1:n_sc
            r_i = SVector{3,Float64}(all_r_x[i][k], all_r_y[i][k], all_r_z[i][k])
            v_i = SVector{3,Float64}(all_v_x[i][k], all_v_y[i][k], all_v_z[i][k])
            phi_i = grav_potential(r_i, mu, RE, J2, use_j2)
            E_tot += 0.5 * mass_kg * dot(v_i, v_i) + mass_kg * phi_i
        end
        k == 1 && (E0 = E_tot)
        delta_Eorb[k] = E_tot - E0

        # ── z-component of total angular momentum over ALL spacecraft (Eq. 17)
        H_z_k = 0.0
        for i in 1:n_sc
            r_i = SVector{3,Float64}(all_r_x[i][k], all_r_y[i][k], all_r_z[i][k])
            v_i = SVector{3,Float64}(all_v_x[i][k], all_v_y[i][k], all_v_z[i][k])
            H_z_i = (r_i[1] * v_i[2] - r_i[2] * v_i[1]) * mass_kg
            H_z_per_sat[i, k] = H_z_i
            H_z_k += H_z_i
        end
        H_z[k] = H_z_k
    end

    # ── Accumulated OCL work from tracker (already integrated at all ODE steps)
    W_ocl = ("W_ocl_accumulated" ∈ names(df)) ?
            Float64.(df[!, "W_ocl_accumulated"]) :
            zeros(npts)

    # ── Energy residual (Eq. 22) ──────────────────────────────────────────────
    eps_E = delta_Eorb .- W_ocl

    # ── Angular-momentum exchange residual from tracker (Eq. 23) ─────────────
    eps_H_x = ("eps_H_x" ∈ names(df)) ? Float64.(df[!, "eps_H_x"]) : zeros(npts)
    eps_H_y = ("eps_H_y" ∈ names(df)) ? Float64.(df[!, "eps_H_y"]) : zeros(npts)
    eps_H_z = ("eps_H_z" ∈ names(df)) ? Float64.(df[!, "eps_H_z"]) : zeros(npts)
    eps_H_mag = sqrt.(eps_H_x.^2 .+ eps_H_y.^2 .+ eps_H_z.^2)

    # ── Shift H_z so it starts at zero for plotting ───────────────────────────
    H_z_change = H_z .- H_z[1]
    # Per-satellite ΔH_z(t) = H_z_i(t) - H_z_i(0)
    H_z_per_sat_change = H_z_per_sat .- H_z_per_sat[:, 1]

    return (
        t                   = t,
        orbit_counts        = orbit_counts,
        F_net_mag           = F_net_mag,
        tau_net_mag         = tau_net_mag,
        delta_Eorb          = delta_Eorb,
        W_ocl               = W_ocl,
        eps_E               = eps_E,
        eps_H_mag           = eps_H_mag,
        H_z_change          = H_z_change,
        H_z_per_sat_change  = H_z_per_sat_change,   # Matrix: n_sc × npts
        n_sc                = n_sc,
        laser_active        = laser_active_bool,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Convenience loader: reads a feather + meta TOML and returns diagnostics
# ─────────────────────────────────────────────────────────────────────────────
"""
    load_and_compute(results_dir, lp; mu, RE, J2, mass_kg) → NamedTuple

Load the Feather written by `run_verification.jl`, read the metadata TOML to
determine whether J2 was active, and return the full diagnostics NamedTuple.

IMPORTANT — the constants here must match what SpaceAGORA actually used:
  • μ  = 3.98600436233e14  (DE421 value, from `Earth().μ`)
  • RE = 6.371008366666667e6  (volumetric MEAN radius from SPICE pck00011,
                               stored in `planet.Rp_m` and used by
                               `_inverse_squared_j2_gravity_accel`).
                               NOT the WGS-84 equatorial radius 6,378,137 m —
                               using that value causes a 0.22 % mismatch in the
                               J₂ potential which appears as a ~150 J sinusoidal
                               energy residual at the orbital period.
  • J2 = 1.08263e-3  (from `Earth().J2`)
"""
function load_and_compute(
    results_dir::String,
    lp::LaserParams;
    mu::Float64      = 3.98600436233e14,         # SpaceAGORA DE421 value
    RE::Float64      = 6.371008366666667e6,       # planet.Rp_m (mean radius, SPICE pck00011)
    J2::Float64      = 1.08263e-3,
    mass_kg::Float64 = 227.0,
)
    feather_path = joinpath(results_dir, "simulation_results.feather")
    isfile(feather_path) || error("Feather file not found: $feather_path")
    df = DataFrame(Arrow.Table(feather_path))

    # Read use_j2 from metadata (default: true if meta not present)
    use_j2 = true
    meta_path = joinpath(results_dir, "verification_meta.toml")
    if isfile(meta_path)
        for line in eachline(meta_path)
            m = match(r"^use_j2\s*=\s*(true|false)", line)
            if m !== nothing
                use_j2 = m.captures[1] == "true"
                break
            end
        end
    end

    n_sc_detected = count(c -> occursin(r"^sc\d+_pos_1$", c), names(df))
    println("  Loading: $results_dir  (J2=$(use_j2), $(nrow(df)) rows, $n_sc_detected spacecraft)")
    return compute_validation_diagnostics(df, lp, mu, RE, J2, mass_kg, use_j2)
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary printer
# ─────────────────────────────────────────────────────────────────────────────
"""
    print_diagnostics_summary(diag; label="")

Print a concise table of peak residuals for quick validation.
"""
function print_diagnostics_summary(diag; label::String="")
    tag = label == "" ? "" : "[$label] "
    println("$(tag)Conservation-law residual summary:")
    @printf("  Peak |F_net^OCL|    = %.3e N\n",    maximum(diag.F_net_mag))
    @printf("  Peak |τ_net^OCL|    = %.3e N·m\n",  maximum(diag.tau_net_mag))
    @printf("  Peak |ε_E|          = %.3e J\n",    maximum(abs.(diag.eps_E)))
    @printf("  Peak |ε_H^OCL|      = %.3e kg·m²/s\n", maximum(diag.eps_H_mag))
    @printf("  ΔE_orb (final)      = %.4e J\n",    diag.delta_Eorb[end])
    @printf("  W_OCL (final)       = %.4e J\n",    diag.W_ocl[end])
    @printf("  ΔH_z   (final)      = %.4e kg·m²/s\n", diag.H_z_change[end])
    @printf("  Orbits elapsed      = %.4f\n",      diag.orbit_counts[end])
    println()
end
