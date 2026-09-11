#!/usr/bin/env julia
# Post-processes SpaceAGORA's feather output (30 orbits, gve_sma/gve_ecc/gve_inc
# schedule) into the same deliverables produced for the prototype code. Pass the
# schedule name as ARGS[1] (default gve_sma) and the target inclination in degrees as
# ARGS[2] (default 0.0).

using Arrow, DataFrames, StaticArrays, Printf, LinearAlgebra

const OUT_DIR = normpath(joinpath(@__DIR__, ".."))
const SCHEDULE = length(ARGS) >= 1 ? ARGS[1] : "gve_sma"
const TARGET_INC_DEG = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.0
const SCENARIO_KEY = TARGET_INC_DEG == 0.0 ? SCHEDULE : @sprintf("%s_it%.1fdeg", SCHEDULE, TARGET_INC_DEG)
const FEATHER_PATH = joinpath(OUT_DIR, SCENARIO_KEY, "spaceagora_raw", "single_case_mode",
    "h1050km_t1000km", @sprintf("ih0.0deg_it%.1fdeg", TARGET_INC_DEG), "N10", "T189214s",
    "$(SCHEDULE)_e0.0000_nu0.0000", "simulation_results.feather")

include(joinpath(@__DIR__, "encounter_analysis.jl"))

isfile(FEATHER_PATH) || error("Feather file not found: $FEATHER_PATH")
df = DataFrame(Arrow.Table(FEATHER_PATH))
t = df.time
n = length(t)
println("Loaded SpaceAGORA feather: $n saved steps, t=[$(t[1]), $(t[end])] s")

const HELPERS = 10
const MAX_RANGE_M = 200e3   # matches --laser-range-km 200 (default)
const HELPER_ALT_KM = 1050.0
const TARGET_ALT_KM = 1000.0
const MAGNIFICATION = 100.0   # default --magnification
const LASER_POWER_W = 10_000.0   # default --laser-power-w

rt = [SVector(df[k, :sc1_pos_1], df[k, :sc1_pos_2], df[k, :sc1_pos_3]) for k in 1:n]
vt = [SVector(df[k, :sc1_vel_1], df[k, :sc1_vel_2], df[k, :sc1_vel_3]) for k in 1:n]

laser_active_helper = df.laser_active_helper   # sc index of active helper (0 = no link)

outdir = joinpath(OUT_DIR, SCENARIO_KEY, "spaceagora")
const CODE_LABEL = "ORACLE with SpaceAGORA engine"
const SCENARIO_LABEL = @sprintf("Helpers @ %.0f km, target @ %.0f km, %.1f° inclination, %.0f km laser range, B=%.0f, P=%.0f kW, %s",
                                 HELPER_ALT_KM, TARGET_ALT_KM, TARGET_INC_DEG, MAX_RANGE_M/1000, MAGNIFICATION, LASER_POWER_W/1000, SCHEDULE)
all_rows = NamedTuple[]
pair_results = NamedTuple[]
for hj in 1:HELPERS
    sc_idx = hj + 1   # helper spacecraft index 2..11 (target = sc1)
    rh = [SVector(df[k, Symbol("sc$(sc_idx)_pos_1")], df[k, Symbol("sc$(sc_idx)_pos_2")], df[k, Symbol("sc$(sc_idx)_pos_3")]) for k in 1:n]
    vh = [SVector(df[k, Symbol("sc$(sc_idx)_vel_1")], df[k, Symbol("sc$(sc_idx)_vel_2")], df[k, Symbol("sc$(sc_idx)_vel_3")]) for k in 1:n]
    active_mask = laser_active_helper .== sc_idx
    pair_label = @sprintf("helper%02d_vs_target", hj)
    res = run_pair_analysis(t, rt, vt, rh, vh, MAX_RANGE_M, pair_label, active_mask, outdir, CODE_LABEL, SCENARIO_LABEL)
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
plot_range_rangerate_combined(encountered, MAX_RANGE_M, CODE_LABEL, SCENARIO_LABEL, joinpath(outdir, "plot2_range_and_rangerate"))

println("\nSpaceAGORA analysis complete. Outputs in: ", outdir)


println("\nSpaceAGORA analysis complete. Outputs in: ", outdir)
