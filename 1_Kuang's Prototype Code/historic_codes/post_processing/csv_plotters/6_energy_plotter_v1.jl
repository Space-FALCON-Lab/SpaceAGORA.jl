using LinearAlgebra, Plots, Printf, StaticArrays, DelimitedFiles

const MU      = 3.986004418e14
const C       = 3.0e8
const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/12_CSV_Write_Read.jl")


# --- Load all CSVs and plot ---
data_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h850km_i0.0deg"))
dir_tag  = basename(data_dir)
csv_files = filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true))

plt = scatter(title="Satellite Orbital Energy ($dir_tag)", xlabel="t (s)",
              ylabel="Total Orbital Energy (J)", formatter=:scientific, legend=:outertopright)

for csv_path in csv_files
    sol, _, p = load_timeseries_csv(csv_path)
    masses = p[:masses]
    Etot = [sum(orbital_energy(u, masses, MU)) for u in sol.u]
    N = p[:N]
    scatter!(plt, sol.t, Etot, label="N=$N", ms=2)
end

# --- Save ---
output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "orbital_energy"))
mkpath(output_dir)
fn = "orbital_energy_$(dir_tag).png"
savefig(plt, joinpath(output_dir, fn))
display(plt)
println("Saved: ", joinpath(output_dir, fn))


