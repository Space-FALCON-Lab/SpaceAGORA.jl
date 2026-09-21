# Test script for saturation_number() in functions/4_Diagnostics.jl
#
# Formula:  α* = acos( (rh² + rt² − L²) / (2·rh·rt·cos(Δi)) )
#           N_sat = ⌈π / α*⌉
#
# Run from the REPL:  include("test_saturation_number.jl")

using LinearAlgebra, StaticArrays, Printf

# Constants (mirror those in the main codebase)
const R_EARTH  = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include("functions/4_Diagnostics.jl")

# ── Helper ────────────────────────────────────────────────────────────────────
function check(label, N_got, N_expected; tol=0)
    ok = (N_got == N_expected)
    status = ok ? "PASS" : "FAIL"
    println("[$status]  $label")
    println("       got=$N_got  expected=$N_expected")
    !ok && println("       *** mismatch! ***")
    println()
end

# ── Manual reference for h_helper=1000 km, h_target=1150 km, Δi=0 ──────────
# r_h = 6 378 137 + 1 000 000 = 7 378 137 m
# r_t = 6 378 137 + 1 150 000 = 7 528 137 m
# arg  = (r_h²+r_t²−L²)/(2·r_h·r_t) ≈ 0.99982
# α*   = acos(0.99982)              ≈ 1.08°
# N    = ⌈π / 0.01885⌉ = ⌈166.6⌉   = 167
r_h_ref  = R_EARTH + 1000e3
r_t_ref  = R_EARTH + 1150e3
arg_ref  = (r_h_ref^2 + r_t_ref^2 - (200e3)^2) / (2 * r_h_ref * r_t_ref)
α_ref    = acos(arg_ref)
N_ref    = ceil(Int, π / α_ref)

println("="^60)
println("  saturation_number() — unit tests")
println("="^60)
println()

# ── Test 1: h_helper=1000 km, h_target=1150 km, Δi=0° ──────────────────────
N1, α1 = saturation_number(1000e3, 0.0, 1150e3, 0.0, 200e3)
println("Test 1 — h_helper=1000 km, h_target=1150 km, Δi=0°, L=200 km")
@printf("  α* = %.4f°,  N_sat = %d\n", α1, N1)
check("N_sat matches manual reference (expected $N_ref)", N1, N_ref)

# # ── Test 2: same as Test 1 (repeat baseline) ────────────────────────────────
# N2, α2 = saturation_number(1000e3, 0.0, 1150e3, 0.0, 200e3)
# println("Test 2 — h_helper=1000 km, h_target=1150 km, Δi=0° (repeat of Test 1)")
# @printf("  α* = %.4f°,  N_sat = %d\n", α2, N2)
# check("N_sat(Test 2) == N_sat(Test 1)", N2, N1)
# println("  Match check (N2 == N1): ", N2 == N1 ? "PASS" : "FAIL", "\n")

# # ── Test 3: same as Test 1 (consistency check) ──────────────────────────────
# N3, α3 = saturation_number(1000e3, 0.0, 1150e3, 0.0, 200e3)
# println("Test 3 — h_helper=1000 km, h_target=1150 km, Δi=0° (consistency check)")
# @printf("  α* = %.4f°,  N_sat = %d\n", α3, N3)
# println("  Match check (N3 == N1): ", N3 == N1 ? "PASS" : "FAIL", "\n")

# # ── Test 4: small inclination difference ──────────────────────────────────────
# N4, α4 = saturation_number(1000e3, 0.0, 1150e3, 0.5, 200e3)
# println("Test 4 — h_helper=1000 km, h_target=1150 km, Δi=0.5°, L=200 km")
# @printf("  α* = %.4f°,  N_sat = %d\n", α4, N4)
# println("  N4 >= N1 check: ", N4 >= N1 ? "PASS" : "FAIL", "\n")

# # ── Test 5: target below helper orbit ───────────────────────────────────────
# # radial gap = 1000 - 950 = 50 km < L=200 km → geometry is feasible
# N5, α5 = saturation_number(1000e3, 0.0, 950e3, 0.0, 200e3)
# println("Test 5 — h_helper=1000 km, h_target=950 km, Δi=0°, L=200 km")
# @printf("  α* = %.4f°,  N_sat = %d\n", α5, N5)
# println()

# # ── Test 6: helper and target on same orbit ──────────────────────────────────
# N6, α6 = saturation_number(1000e3, 0.0, 1000e3, 0.0, 200e3)
# println("Test 6 — h_helper=h_target=1000 km, Δi=0° (same orbit)")
# @printf("  α* = %.4f°,  N_sat = %d\n", α6, N6)
# println()

# ── Test 7: simulation grid — tabulate N_sat for all (h, i) pairs ────────────
println("─"^60)
println("  Saturation numbers for simulation parameter grid")
println("  L_max = 200 km,  h_helper = 1000 km (fixed),  h_target varies")
println("  N/A = geometry infeasible (max orbit separation > L_max)")
println("─"^60)
altitudes_m  = [850e3, 950e3, 1000e3, 1050e3, 1150e3]
incl_diffs   = [0.0, 0.5, 1.0, 1.5]
L_m          = 200e3
h_helper_m   = 1000e3

print("  Δi →      ")
for di in incl_diffs; @printf("%8.1f°", di); end; println()
println("  h_target ↓")
for h_t in altitudes_m
    @printf("  %4.0f km   ", h_t/1e3)
    for di in incl_diffs
        result = try
            N, _ = saturation_number(h_helper_m, 0.0, h_t, di, L_m)
            @sprintf("%9d", N)
        catch
            @sprintf("%9s", "N/A")
        end
        print(result)
    end
    println()
end
println()

# # ── Test 8: error on impossible geometry ──────────────────────────────────────
# println("Test 8 — Δi = 90°, L=1 m  (should throw an error)")
# try
#     saturation_number(1000e3, 0.0, 1000e3, 90.0, 1.0)
#     println("[FAIL]  No error thrown — expected an error!\n")
# catch e
#     println("[PASS]  Caught expected error:")
#     println("        ", e.msg, "\n")
# end

# println("="^60)
# println("  All tests complete.")
# println("="^60)
