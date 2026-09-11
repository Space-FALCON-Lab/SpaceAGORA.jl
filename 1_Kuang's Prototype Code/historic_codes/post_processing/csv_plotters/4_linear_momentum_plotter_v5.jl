using LinearAlgebra, Plots, Printf, StaticArrays, DelimitedFiles

const MU      = 3.986004418e14
const C       = 3.0e8
const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/12_CSV_Write_Read.jl")

# --- Load all CSVs and plot ---
data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h850km_i0.0deg"))
data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h850km_i0.0deg"))
dir_tag   = basename(data_dir)
csv_files = filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true))

plt1 = scatter(title="Net Satellite Momentum ($dir_tag)", xlabel="t (s)",
               ylabel="Total Momentum (kg·m/s)", formatter=:scientific, legend=:outertopright)
plt2 = scatter(title="Linear Momentum Fractional Drift ($dir_tag)", xlabel="t (s)",
               ylabel="ΔP/P(0) = (P(t) - P(0)) / P(0)", formatter=:scientific, legend=:outertopright)

for csv_path in csv_files
    sol, _, p = load_timeseries_csv(csv_path)
    masses = p[:masses]
    Pmag = [total_momentum(u, masses)[2] for u in sol.u]
    ΔP   = (Pmag .- Pmag[1]) ./ Pmag[1]
    N = p[:N]
    scatter!(plt1, sol.t, Pmag, label="N=$N", ms=2, markerstrokewidth=0)
    scatter!(plt2, sol.t, ΔP,   label="N=$N", ms=2, markerstrokewidth=0)
end

# --- Save ---
output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "linear_momentum"))
mkpath(output_dir)
savefig(plt1, joinpath(output_dir, "linear_momentum_$(dir_tag).png"))
savefig(plt2, joinpath(output_dir, "linear_momentum_fractional_drift_$(dir_tag).png"))
display(plt1); display(plt2)
println("Saved to: ", output_dir)


