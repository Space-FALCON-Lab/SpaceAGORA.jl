using LinearAlgebra, StaticArrays, Plots, LaTeXStrings, Printf, DelimitedFiles

if !@isdefined(MU);      const MU      = 3.986004418e14; end
if !@isdefined(C);       const C       = 3.0e8;          end
if !@isdefined(R_EARTH); const R_EARTH = 6_378_137.0;    end
if !@isdefined(R_ATMDEF)
    const R_ATMDEF = 6_478_137.0
end
@inline idx(i, off) = 6*(i-1) + off

include("../functions/1_LOS_Metrics.jl")
include("../functions/2_Laser_Forces_ver2.jl")
include("../functions/4_Diagnostics.jl")
include("../functions/5_OE_Converters.jl")
include("../functions/6_OE_and_dv_in_RTN.jl")
include("../functions/7_Plots.jl")
include("../functions/12_CSV_Write_Read.jl")

# ── Target CSV ───────────────────────────────────────────────────────────────
# Produced by Main31_coupling_v3.jl with:
#   helper_counts=1, alt_km=1000, inc_deg=0.0, nu_deg=1.0
#   T_seconds=150*3600, use_J2=true
csv_path = normpath(joinpath(@__DIR__, "..", "output", "CSV",
    "target_h1000km_i0.0deg_nu1.00deg",
    "timeseries_N2_T540000s_h1000km_t1000km_ih0.0deg_it0.0deg_B100_Pin1e+04_rmin0m_rmax2e+05_J2T_toyT.csv"))

output_dir = normpath(joinpath(@__DIR__, "..", "output", "images_reconstructed_from_csv", "OE_diff"))
mkpath(output_dir)

# ── Load solution ─────────────────────────────────────────────────────────────
sol, metadata, p = load_timeseries_csv(csv_path)
println("Loaded CSV: $(basename(csv_path))")
println("  N satellites : $(p[:N])")
println("  Time span    : $(sol.t[1]) – $(sol.t[end]) s  ($(length(sol.t)) steps)")

# ── Plot OE differences: Sat 2 − Sat 1 ───────────────────────────────────────
# sat1 = 1 (helper), sat2 = 2 (target)
plts = report_and_plot_OE_diff(sol, p[:mu];
                                sat1=1, sat2=2,
                                IMG_DIR=output_dir,
                                fn_prefix="OE_diff_sat1_sat2_toy")

plt_a, plt_e, plt_i, plt_Ω, plt_ω, plt_v = plts

for sp in plts
    display(sp)
end

println("Saved OE diff plots to: $output_dir/OE_diff_sat1_sat2_toy/")
