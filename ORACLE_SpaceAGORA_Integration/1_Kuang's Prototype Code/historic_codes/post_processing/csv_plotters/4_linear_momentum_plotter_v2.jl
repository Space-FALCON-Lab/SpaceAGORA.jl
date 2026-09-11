using LinearAlgebra, Plots, Printf, StaticArrays, DelimitedFiles

const MU      = 3.986004418e14
const C       = 3.0e8
const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/12_CSV_Write_Read.jl")

# --- Load CSV ---
csv_path = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV",
                             "target_h850km_i0.0deg",
                             "timeseries_N301_T180000s_h1000km_t850km_ih0.0deg_it0.0deg.csv"))
csv_tag  = replace(basename(csv_path), ".csv" => "")
sol, _, p = load_timeseries_csv(csv_path)
masses = p[:masses]
N = length(masses)

# --- Compute momentum ---
Pvals = [total_momentum(u, masses) for u in sol.u] # compute total momentum at each time step
# p[k] = vector total linear momentum P, magnitude total linear momentum norm(P), magnitude linear momentum for each satellitepmag
Pmag  = [Pvals[k][2] for k in eachindex(sol.t)] # extract the magnitude of the total momentum for plotting

# --- Plot ---
plt = plot(sol.t, [Pvals[k][3][1] for k in eachindex(sol.t)],
           title="Net Satellite Momentum ($csv_tag)", xlabel="t (s)", ylabel="Momentum (kg·m/s)",
           label="Sat 1", legend=:topleft, lw=1.5)
ax2 = twinx(plt)
scatter!(ax2, sol.t, Pmag, label="Total", ms=2, color=:red,
         ylabel="Total Momentum (kg·m/s)", formatter=:scientific, legend=:topright)

# --- Save ---
output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "linear_momentum"))
mkpath(output_dir)
fn = "linear_momentum_$(csv_tag).png"
savefig(plt, joinpath(output_dir, fn))
display(plt)
println("Saved: ", joinpath(output_dir, fn))


