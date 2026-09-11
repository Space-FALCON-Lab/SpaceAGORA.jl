#!/usr/bin/env julia
# Runs Kuang's prototype code for 30 target orbits with a GVE helper-selection scheduler
# (gve_sma, gve_ecc, or gve_inc, passed as ARGS[1], default gve_sma) and a target
# inclination (deg, ARGS[2], default 0.0), then produces the comparison-report
# deliverables for every helper<->target pair. Does not modify the original
# test16_options.jl scenario file (its two dead-wiring quirks — T_seconds hardcoded to
# 63071s, gve_schedule hardcoded to :none — are worked around here rather than "fixed in
# place").

import GeometryBasics: Point3f, Vec3f
using OrdinaryDiffEq, DiffEqCallbacks
using LinearAlgebra, StaticArrays
using Plots, Printf
using Statistics
using DelimitedFiles

const PROTOTYPE_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "1_Kuang's Prototype Code"))
const OUT_DIR        = normpath(joinpath(@__DIR__, ".."))
const SCHEDULE       = length(ARGS) >= 1 ? Symbol(ARGS[1]) : :gve_sma
const TARGET_INC_DEG_ARG = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.0
const SCENARIO_KEY   = TARGET_INC_DEG_ARG == 0.0 ? String(SCHEDULE) : @sprintf("%s_it%.1fdeg", SCHEDULE, TARGET_INC_DEG_ARG)

#############
# Constants #
#############
const MU       = 3.986004418e14
const C        = 3.0e8
const R_EARTH  = 6_378_137.0
const R_ATMDEF = R_EARTH + 100_000.0
const Ẑ       = SVector(0.0, 0.0, 1.0)

@inline idx(i, off) = 6*(i-1) + off

include(joinpath(PROTOTYPE_DIR, "functions", "1_LOS_Metrics.jl"))
include(joinpath(PROTOTYPE_DIR, "functions", "2_Laser_Forces_ver2.jl"))
include(joinpath(PROTOTYPE_DIR, "functions", "3_Dynamics.jl"))
include(joinpath(PROTOTYPE_DIR, "functions", "4_Diagnostics.jl"))
include(joinpath(PROTOTYPE_DIR, "functions", "5_OE_Converters.jl"))
include(joinpath(PROTOTYPE_DIR, "functions", "6_OE_and_dv_in_RTN.jl"))
include(joinpath(PROTOTYPE_DIR, "functions", "7_Plots.jl"))
include(joinpath(PROTOTYPE_DIR, "functions", "8_LoS_time_series.jl"))
include(joinpath(PROTOTYPE_DIR, "functions", "9_Runners.jl"))
include(joinpath(PROTOTYPE_DIR, "functions", "12_CSV_Write_Read.jl"))

include(joinpath(@__DIR__, "encounter_analysis.jl"))

############
# Settings #
############
const HELPERS               = 10
const HELPER_ALT_KM         = 1050.0
const TARGET_ALT_KM         = 1000.0
const TARGET_INC_DEG        = TARGET_INC_DEG_ARG
const HELPER_INC_DEG        = 0.0
const TARGET_NU_DEG         = 0.0
const ORBITS                = 30.0
const LASER_RANGE_KM        = 200.0
const LASER_POWER_W         = 10_000.0
const MAGNIFICATION         = 100.0
const MASS_KG                = 227.0

println("\nSTART SIMULATION (prototype, $(SCHEDULE) schedule, target inc=$(TARGET_INC_DEG)°, 30 orbits):")

helper_num = HELPERS
helper_oe = [
    (a_m = R_EARTH + HELPER_ALT_KM*1e3, e = 0.0,
     i_deg = HELPER_INC_DEG, Ω_deg = 0.0, ω_deg = 0.0,
     ν_deg = (360.0 / helper_num) * (j - 1)) for j in 1:helper_num
]
target_orbit = (a_m = R_EARTH + TARGET_ALT_KM*1e3, e = 0.0,
                i_deg = TARGET_INC_DEG, Ω_deg = 0.0, ω_deg = 0.0,
                ν_deg = TARGET_NU_DEG)
oe = vcat(helper_oe, [target_orbit])
target_idx = helper_num + 1   # satellite N = target (matches test16_options.jl layout)

Pm  = zeros(length(oe), length(oe))
cav = Dict{Tuple{Int,Int},Dict{Symbol,Any}}()
for j in 1:helper_num
    cav[(j, target_idx)] = Dict(:B => MAGNIFICATION, :Pin => LASER_POWER_W)
end

a_target  = R_EARTH + TARGET_ALT_KM*1e3
T_orbit   = 2π * sqrt(a_target^3 / MU)
T_seconds = ORBITS * T_orbit
max_range_m = LASER_RANGE_KM * 1e3

elapsed = @elapsed begin
    global sol, p, logs, masses = run_open_cavity_multi(oe;
        mass_kg      = MASS_KG,
        Pm           = Pm,
        cavity       = cav,
        use_los      = true,
        min_range    = 0.0,
        max_range    = max_range_m,
        stop_on_dv   = false,
        T_seconds    = T_seconds,        # << fixed: was hardcoded 63071 (=10 orbits) in test16_options.jl
        verbose      = true,
        result_plots = false,            # skip prototype's own plot suite; we generate our own below
        target_only  = true,
        IMG_DIR      = joinpath(OUT_DIR, SCENARIO_KEY, "prototype", "_unused_images") * "/",
        helper_num   = helper_num,
        gve_schedule = SCHEDULE,         # << fixed: was hardcoded :none in test16_options.jl (opts.schedule was dead)
    )
end
println(@sprintf("Simulation runtime: %.2f s (%.2f min), %d saved steps", elapsed, elapsed/60, length(sol.t)))

# Save raw timeseries CSV for the record
case_csv_dir = joinpath(OUT_DIR, SCENARIO_KEY, "prototype", "csv")
mkpath(case_csv_dir)
save_timeseries_csv(sol, p, helper_oe, target_orbit, csv_dir=case_csv_dir)

###################
# Post-processing #
###################
t = sol.t
n = length(t)

rt = [SVector(sol.u[k][idx(target_idx,1)], sol.u[k][idx(target_idx,2)], sol.u[k][idx(target_idx,3)]) for k in 1:n]
vt = [SVector(sol.u[k][idx(target_idx,4)], sol.u[k][idx(target_idx,5)], sol.u[k][idx(target_idx,6)]) for k in 1:n]

# Which helper (if any) has an active open-cavity link to the target at each saved step,
# recomputed exactly as the GVE scheduler selected it during simulation. NOTE: laser_forces'
# `current_helpers` busy-matrix marks entire rows/columns (not just the single active cell),
# so it can't be used to identify which specific helper is active — call the scheduler's own
# selection function directly instead (deterministic given state + p[:gve_target_idx]).
N_total = helper_num + 1
active_helper = zeros(Int, n)
for k in 1:n
    u = sol.u[k]
    r = Array{Float64}(undef, 3, N_total)
    for i in 1:N_total
        r[1,i] = u[idx(i,1)]; r[2,i] = u[idx(i,2)]; r[3,i] = u[idx(i,3)]
    end
    best = _gve_select_best_cavity(cav, target_idx, u, r, SCHEDULE,
                                    true, R_ATMDEF, 5_000.0, 0.0, max_range_m, MU)
    if !isempty(best)
        (i, j) = first(keys(best))
        active_helper[k] = i == target_idx ? j : i
    end
end

outdir = joinpath(OUT_DIR, SCENARIO_KEY, "prototype")
const CODE_LABEL = "ORACLE prototype code"
const SCENARIO_LABEL = @sprintf("Helpers @ %.0f km, target @ %.0f km, %.1f° inclination, %.0f km laser range, B=%.0f, P=%.0f kW, %s",
                                 HELPER_ALT_KM, TARGET_ALT_KM, TARGET_INC_DEG, LASER_RANGE_KM, MAGNIFICATION, LASER_POWER_W/1000, SCHEDULE)
all_rows = NamedTuple[]
pair_results = NamedTuple[]
for j in 1:helper_num
    rh = [SVector(sol.u[k][idx(j,1)], sol.u[k][idx(j,2)], sol.u[k][idx(j,3)]) for k in 1:n]
    vh = [SVector(sol.u[k][idx(j,4)], sol.u[k][idx(j,5)], sol.u[k][idx(j,6)]) for k in 1:n]
    active_mask = active_helper .== j
    pair_label = @sprintf("helper%02d_vs_target", j)
    res = run_pair_analysis(t, rt, vt, rh, vh, max_range_m, pair_label, active_mask, outdir, CODE_LABEL, SCENARIO_LABEL)
    append!(all_rows, res.rows)
    push!(pair_results, res)
    println(@sprintf("  %s: %d encounter(s)", pair_label, length(res.rows)))
end

write_contact_window_table(all_rows,
    joinpath(outdir, "plot3_contact_windows", "contact_windows.csv"),
    joinpath(outdir, "plot3_contact_windows", "contact_windows.md"),
    joinpath(outdir, "plot3_contact_windows", "contact_windows_bar.png"), CODE_LABEL, SCENARIO_LABEL)

encountered = filter(r -> !isempty(r.encounters), pair_results)
plot_relspeed_combined(encountered, CODE_LABEL, SCENARIO_LABEL, joinpath(outdir, "plot1_relative_speed"))
plot_range_rangerate_combined(encountered, max_range_m, CODE_LABEL, SCENARIO_LABEL, joinpath(outdir, "plot2_range_and_rangerate"))

println("\nPrototype analysis complete. Outputs in: ", outdir)

