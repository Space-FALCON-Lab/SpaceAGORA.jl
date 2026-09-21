"""
SectionF_validation_main.jl

Section F — Validation Through Numerical Simulation
Paper: "Assessing Small Satellite Maneuverability in Proliferated LEO
        Constellations with Open Cavity Laser-Interlinks"

Runs a two-satellite OCL simulation and verifies the four conservation-law
residuals defined in Eqs. (21)–(23) of the paper:

  • Net OCL force  F_net^(OCL)   (should be ≈ 0 at machine precision)
  • Net OCL torque τ_net^(OCL)   (should be ≈ 0 at machine precision)
  • Energy residual ε_E(t)       (should stay at numerical integration error)
  • Angular-momentum residual ε_H^(OCL)(t) (should stay ≈ 0)

OCL parameters (Table 2):
  Laser range     L_max  = 200 km
  Laser power     P_ij   = 10 000 W
  OCL magnific.   B      = 100
  Geometric loss  β_ij   = 1.0  (absorbed into B · P / c formula)

Satellite configuration:
  Helper: circular orbit, altitude = 1000 km, inclination = 0°
  Target: circular orbit, altitude = 1100 km, inclination = 0°
  Altitude separation = 100 km  (→ intermittent link case)

Two runs are performed:
  1. No J2 perturbation
  2. With J2 perturbation

Outputs:
  Figures saved to  post_processing/output/
  Console summary of residual magnitudes
"""

# ── Activate GR backend (no GUI window needed, saves PNGs) ────────────────────
using Plots
gr()

using OrdinaryDiffEq
using LinearAlgebra, StaticArrays
using Printf

# ── Physical constants (must match functions/) ────────────────────────────────
const MU      = 3.986004418e14          # Earth μ  [m³/s²]
const C       = 3.0e8                   # speed of light [m/s]
const R_EARTH = 6_378_137.0            # Earth mean radius [m]
const R_ATMDEF = R_EARTH + 100_000.0   # Kármán line [m]
const Ẑ       = SVector(0.0, 0.0, 1.0)

@inline idx(i, off) = 6*(i-1) + off   # state-vector index helper

# ── Include shared function libraries ─────────────────────────────────────────
const FUNC_DIR = normpath(joinpath(@__DIR__, "..", "functions"))
include(joinpath(FUNC_DIR, "1_LOS_Metrics.jl"))
include(joinpath(FUNC_DIR, "2_Laser_Forces_ver2.jl"))
include(joinpath(FUNC_DIR, "3_Dynamics.jl"))
include(joinpath(FUNC_DIR, "4_Diagnostics.jl"))
include(joinpath(FUNC_DIR, "5_OE_Converters.jl"))

# ── Include Section-F specific functions ──────────────────────────────────────
include(joinpath(@__DIR__, "SectionF_diagnostics.jl"))
include(joinpath(@__DIR__, "SectionF_plots.jl"))

# ── Output directory ──────────────────────────────────────────────────────────
const OUT_DIR = normpath(joinpath(@__DIR__, "output", "SectionF"))
mkpath(OUT_DIR)

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  Simulation parameters
# ═══════════════════════════════════════════════════════════════════════════════
# OCL parameters — Table 2
const L_MAX    = 200_000.0   # maximum laser-link range [m]
const P_LASER  = 1.0e4       # transmitted laser power  [W]
const B_FACTOR = 100.0       # OCL thrust magnification factor
const BETA     = 1.0         # geometric transmission factor (already = 1)
# The OCL force magnitude is  f = B * P / c  (Eq. 4 with β=1, K=1)
# which equals 100 × 10000 / 3e8 ≈ 3.33 × 10⁻³ N

# Satellite masses
const MASS_KG  = 227.0       # [kg]  (same as in Table 3 / paper default)

# Orbital configuration
const H_HELPER_KM = 1000.0   # helper altitude  [km]
const H_TARGET_KM = 1100.0   # target altitude  [km]

# Simulation duration: 1.2 target orbital periods
# T_orbit(1100 km) ≈ 2π √( (R_E + 1100e3)³ / μ )
let a_target = R_EARTH + H_TARGET_KM*1e3
    global T_ORBIT_TARGET = 2π * sqrt(a_target^3 / MU)
end
const N_ORBITS_SIM = 1.2
const T_SIM        = N_ORBITS_SIM * T_ORBIT_TARGET

@printf("\n── Simulation parameters ──────────────────────────────────────────\n")
@printf("  Helper altitude       : %.0f km\n",  H_HELPER_KM)
@printf("  Target altitude       : %.0f km\n",  H_TARGET_KM)
@printf("  Altitude separation   : %.0f km\n",  H_TARGET_KM - H_HELPER_KM)
@printf("  Target orbital period : %.1f s  (%.3f h)\n", T_ORBIT_TARGET, T_ORBIT_TARGET/3600)
@printf("  Simulation time       : %.1f orbits  = %.1f s\n", N_ORBITS_SIM, T_SIM)
@printf("  OCL range             : %.0f km\n",  L_MAX/1e3)
@printf("  OCL power             : %.0f W\n",   P_LASER)
@printf("  OCL magnification B   : %.0f\n",     B_FACTOR)
@printf("  Satellite mass        : %.0f kg\n",  MASS_KG)
@printf("──────────────────────────────────────────────────────────────────\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  Build initial state vector
# ═══════════════════════════════════════════════════════════════════════════════
# Satellite 1 = helper  (1000 km, i=0°)
# Satellite 2 = target  (1100 km, i=0°)
# Both start at argument of latitude u = 0° (same ascending node).

oe_helper = (a_m=R_EARTH + H_HELPER_KM*1e3, e=0.0, i_deg=0.0,
             Ω_deg=0.0, ω_deg=0.0, ν_deg=0.0)
oe_target = (a_m=R_EARTH + H_TARGET_KM*1e3, e=0.0, i_deg=0.0,
             Ω_deg=0.0, ω_deg=0.0, ν_deg=0.0)

s_helper = state_from_OE(oe_helper.a_m; e=oe_helper.e, i_deg=oe_helper.i_deg,
                          Ω_deg=oe_helper.Ω_deg, ω_deg=oe_helper.ω_deg,
                          ν_deg=oe_helper.ν_deg)
s_target = state_from_OE(oe_target.a_m; e=oe_target.e, i_deg=oe_target.i_deg,
                          Ω_deg=oe_target.Ω_deg, ω_deg=oe_target.ω_deg,
                          ν_deg=oe_target.ν_deg)

u0 = vcat(collect(s_helper), collect(s_target))

# ═══════════════════════════════════════════════════════════════════════════════
# 3.  Build parameter dictionary
# ═══════════════════════════════════════════════════════════════════════════════
N_SATS = 2
masses = fill(MASS_KG, N_SATS)

Pm = zeros(N_SATS, N_SATS)   # no single-pass beams; all thrust via cavity

# OCL cavity: helper (index 1) ↔ target (index 2)
# No leakage terms — pure cavity validation case
cavity = Dict{Tuple{Int,Int}, Dict{Symbol,Any}}()
cavity[(1, 2)] = Dict(:B => B_FACTOR, :Pin => P_LASER)

function make_p(; use_J2::Bool)
    return Dict(
        :mu           => MU,
        :c            => C,
        :N            => N_SATS,
        :masses       => masses,
        :Pmatrix      => Pm,
        :cavity       => cavity,
        :use_los      => true,
        :R_atm        => R_ATMDEF,
        :atm_clearance => 0.0,
        :min_range    => 0.0,
        :max_range    => L_MAX,
        :use_J2       => use_J2,
        :J2           => 1.08262668e-3,
        :Re           => R_EARTH,
        # helper_ids / target_ids used by CSV writer only
        :helper_ids   => [1],
        :target_ids   => [2],
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# 4.  Run simulations
# ═══════════════════════════════════════════════════════════════════════════════

function run_case(; use_J2::Bool)
    label = use_J2 ? "J2 active" : "J2 inactive"
    println("\n════════════════════════════════════════")
    println("Running: ", label)
    println("════════════════════════════════════════")

    p = make_p(use_J2=use_J2)
    prob = ODEProblem(nbody_photon!, u0, (0.0, T_SIM), p)

    elapsed = @elapsed begin
        sol = solve(prob, Vern9();
                    reltol        = 1e-12,
                    abstol        = 1e-12,
                    saveat        = T_SIM / 2000,  # ~2000 output points for smooth plots
                    save_everystep = false)
    end
    @printf("  ODE solve: %.2f s\n", elapsed)
    @printf("  Steps saved: %d\n", length(sol.t))

    return sol, p
end

sol_no_j2, p_no_j2 = run_case(use_J2=false)
sol_j2,    p_j2    = run_case(use_J2=true)

# ═══════════════════════════════════════════════════════════════════════════════
# 5.  Compute diagnostics
# ═══════════════════════════════════════════════════════════════════════════════
println("\nComputing conservation-law diagnostics (no J2) …")
diag_no_j2 = compute_sectionF_diagnostics(sol_no_j2, p_no_j2;
                                           helper_idx=1, target_idx=2)

println("Computing conservation-law diagnostics (J2 active) …")
diag_j2    = compute_sectionF_diagnostics(sol_j2, p_j2;
                                           helper_idx=1, target_idx=2)

# ── Console summary ────────────────────────────────────────────────────────────
print_diagnostic_summary(diag_no_j2; label="J2 inactive")
print_diagnostic_summary(diag_j2;    label="J2 active")

# ═══════════════════════════════════════════════════════════════════════════════
# 6.  Generate plots  (Figures 2 and 3 of the paper)
# ═══════════════════════════════════════════════════════════════════════════════
println("\nGenerating plots …")

fig2 = plot_figure2(diag_no_j2, diag_j2;
                    Δh_label = "Δh = 100 km",
                    save_dir  = OUT_DIR,
                    filename  = "SectionF_Figure2_force_torque_angmom.png")

fig3 = plot_figure3(diag_no_j2, diag_j2;
                    Δh_label = "Δh = 100 km",
                    save_dir  = OUT_DIR,
                    filename  = "SectionF_Figure3_energy_residual.png")

println("\nDone. Output saved to: ", OUT_DIR)

# ═══════════════════════════════════════════════════════════════════════════════
# 7.  Optional: display figures interactively
# ═══════════════════════════════════════════════════════════════════════════════
display(fig2)
display(fig3)
