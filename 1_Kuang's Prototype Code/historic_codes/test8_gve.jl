import Makie
import WGLMakie
WGLMakie.activate!()
const GLMakie = Makie
import GeometryBasics: Point3f, Vec3f
using OrdinaryDiffEq, DiffEqCallbacks
using LinearAlgebra, StaticArrays
using Plots, Printf
using Statistics
using Revise
using FileIO
using GeometryBasics
using Observables
using DelimitedFiles

print("\033c")  # Clears the terminal on Windows

#############
# Constants #
#############
const MU       = 3.986004418e14         # Earth μ [m^3/s^2]
const C        = 3.0e8                  # speed of light [m/s]
const R_EARTH  = 6_378_137.0           # Earth mean radius [m]
const R_ATMDEF = R_EARTH + 100_000.0    # default "atmosphere radius" (Kármán line) [m]
const Ẑ       = SVector(0.0, 0.0, 1.0)

@inline idx(i, off) = 6*(i-1) + off  # state indexing helper

#############
# Functions #
#############
include("functions/1_LOS_Metrics.jl")
include("functions/2_Laser_Forces_ver2.jl")
include("functions/3_Dynamics.jl")
include("functions/4_Diagnostics.jl")
include("functions/5_OE_Converters.jl")
include("functions/6_OE_and_dv_in_RTN.jl")
include("functions/7_Plots.jl")
include("functions/8_LoS_time_series.jl")
include("functions/9_Runners.jl")
include("functions/10_Animation_ver2.jl")
include("functions/12_CSV_Write_Read.jl")

########
# Main #
########
println("\nSTART GVE SCHEDULER TEST:")

# -------------------------------------------------------------------------
# Simulation setup (identical geometry to test8.jl)
# -------------------------------------------------------------------------
helper_num             = 50
target_altitude_km     = 1000.0
target_inclination_deg = 0.0
target_nu_deg          = 0.0

min_range = 0.0
max_range = 200e3          # laser range limit [m]

# Cavity parameters shared by every helper link
CAVITY_B   = 100.0         # power-buildup factor
CAVITY_PIN = 1e4           # input power [W]

# Simulation duration
T_seconds = 1 * 3600     # 150 hours (same as test8.jl)

# -------------------------------------------------------------------------
# Schedules to compare.
# :none        — original behavior: activate every in-range helper simultaneously
# :gve_sma     — activate only the helper that maximises da/dt (semi-major axis)
# :gve_ecc     — activate only the helper that maximises de/dt (eccentricity)
# :gve_inc     — activate only the helper that maximises di/dt (inclination)
# :gve_raan    — activate only the helper that maximises dΩ/dt (RAAN)
# :gve_argp    — activate only the helper that maximises dω/dt (arg. of periapsis)
# -------------------------------------------------------------------------
SCHEDULES_TO_RUN = [:gve_sma] #[:none, :gve_sma, :gve_ecc, :gve_inc]

println("\nHelper count = $helper_num")
println("Schedules    = $SCHEDULES_TO_RUN")
println("Mission time = $(T_seconds/3600) hr")
println()

results_summary = []

for schedule in SCHEDULES_TO_RUN
    println("─────────────────────────────────────────────────────")
    @printf("Schedule: %s\n", schedule)

    # --- Build orbital elements ---
    helper_oe = [
        (a_m=R_EARTH+1000e3, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0,
         ν_deg=(360.0 / helper_num) * (j - 1)) for j in 1:helper_num
    ]
    target_orbit = (a_m=R_EARTH + target_altitude_km*1e3,
                    e=0.0, i_deg=target_inclination_deg,
                    Ω_deg=0.0, ω_deg=0.0, ν_deg=target_nu_deg)
    oe = vcat(helper_oe, [target_orbit])
    N  = length(oe)            # total satellites
    target_sat_idx = N         # target is always last

    # --- Cavity map: all helpers linked to the target ---
    Pm  = zeros(N, N)
    cav = Dict{Tuple{Int,Int},Dict{Symbol,Any}}()
    for j in 1:helper_num
        cav[(j, target_sat_idx)] = Dict(:B => CAVITY_B, :Pin => CAVITY_PIN)
    end

    # --- Run simulation ---
    elapsed = @elapsed begin
        sol, p, logs, masses = run_open_cavity_multi(
            oe;
            mass_kg      = 227.0,
            Pm           = Pm,
            cavity       = cav,
            use_los      = true,
            min_range    = min_range,
            max_range    = max_range,
            stop_on_dv   = false,
            T_seconds    = T_seconds,
            verbose      = false,
            result_plots = true,    # each schedule writes to its own IMG_DIR
            target_only  = true,
            IMG_DIR      = "Kuang's Prototype Code/output/images_gve_$(schedule)/",
            helper_num   = helper_num,
            # ── GVE scheduler parameters ─────────────────────────────
            # :gve_schedule  — which element to maximise (:none = default behavior)
            # :gve_target_idx — index of the target satellite (N = last in oe list)
            gve_schedule    = schedule,
            gve_target_idx  = target_sat_idx,
        )
    end

    # --- Extract final orbital elements of the target ---
    u_final   = sol.u[end]
    r_f = SVector{3,Float64}(u_final[idx(target_sat_idx,1)],
                              u_final[idx(target_sat_idx,2)],
                              u_final[idx(target_sat_idx,3)])
    v_f = SVector{3,Float64}(u_final[idx(target_sat_idx,4)],
                              u_final[idx(target_sat_idx,5)],
                              u_final[idx(target_sat_idx,6)])
    u_init    = sol.u[1]
    r_0 = SVector{3,Float64}(u_init[idx(target_sat_idx,1)],
                              u_init[idx(target_sat_idx,2)],
                              u_init[idx(target_sat_idx,3)])
    v_0 = SVector{3,Float64}(u_init[idx(target_sat_idx,4)],
                              u_init[idx(target_sat_idx,5)],
                              u_init[idx(target_sat_idx,6)])
    oe0 = rv2coe(r_0, v_0, MU)
    oef = rv2coe(r_f, v_f, MU)

    da_m    = oef.a - oe0.a
    de      = oef.e - oe0.e
    di_deg  = rad2deg(oef.i - oe0.i)
    dΩ_deg  = rad2deg(oef.Ω - oe0.Ω)

    @printf("  Elapsed: %.2f s (%.2f min)\n", elapsed, elapsed/60)
    @printf("  Δa  = %+.3f m\n",   da_m)
    @printf("  Δe  = %+.6f\n",     de)
    @printf("  Δi  = %+.6f deg\n", di_deg)
    @printf("  ΔΩ  = %+.6f deg\n", dΩ_deg)

    push!(results_summary, (
        schedule = schedule,
        da_m     = da_m,
        de       = de,
        di_deg   = di_deg,
        dΩ_deg   = dΩ_deg,
        elapsed_s = elapsed,
    ))
end

# -------------------------------------------------------------------------
# Summary table
# -------------------------------------------------------------------------
println("\n╔══════════════════════════════════════════════════════════════╗")
println("║                  GVE Scheduler Comparison                   ║")
println("╠═══════════════╦══════════════╦══════════╦══════════╦═════════╣")
println("║ Schedule      ║   Δa [m]     ║  Δe      ║ Δi [deg] ║ ΔΩ[deg] ║")
println("╠═══════════════╬══════════════╬══════════╬══════════╬═════════╣")
for r in results_summary
    @printf("║ %-13s ║ %+12.3f ║ %+8.5f ║ %+8.5f ║ %+7.4f ║\n",
            r.schedule, r.da_m, r.de, r.di_deg, r.dΩ_deg)
end
println("╚═══════════════╩══════════════╩══════════╩══════════╩═════════╝")
println("\nDone.")
