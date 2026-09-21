using LinearAlgebra, Plots, Printf, StaticArrays, DelimitedFiles

const MU      = 3.986004418e14
const C       = 3.0e8
const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/12_CSV_Write_Read.jl")


# --- Load all CSVs and plot ---
data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h850km_i0.0deg"))
dir_tag   = basename(data_dir)
csv_files = filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true))

plt1 = scatter(title="Satellite Orbital Energy Variation ($dir_tag)", xlabel="t (s)",
              ylabel="ΔE = E(t) - E(0) (J)", formatter=:scientific, legend=:outertopright)
plt2 = scatter(title="Orbital Energy Fractional Drift ($dir_tag)", xlabel="t (s)",
              ylabel="ΔE/|E(0)| = (E(t) - E(0)) / |E(0)|", formatter=:scientific, legend=:outertopright)

for csv_path in csv_files
    sol, _, p = load_timeseries_csv(csv_path)
    masses = p[:masses]
    Etot = [sum(orbital_energy(u, masses, MU)) for u in sol.u]
    ΔE      = Etot .- Etot[1]
    ΔE_frac = ΔE ./ abs(Etot[1])
    N = p[:N]

    # 1000 points evenly spaced in time via linear interpolation
    t_uniform = collect(LinRange(sol.t[1], sol.t[end], 2000))
    function interp(vals, t_query)
        i = clamp(searchsortedlast(sol.t, t_query), 1, length(sol.t)-1)
        α = (t_query - sol.t[i]) / (sol.t[i+1] - sol.t[i])
        return vals[i] * (1-α) + vals[i+1] * α
    end
    ΔE_u      = [interp(ΔE,      t) for t in t_uniform]
    ΔE_frac_u = [interp(ΔE_frac, t) for t in t_uniform]

    scatter!(plt1, t_uniform, ΔE_u,      label="N=$N", ms=3, markerstrokewidth=0)
    scatter!(plt2, t_uniform, ΔE_frac_u, label="N=$N", ms=3, markerstrokewidth=0)
end

# --- Save ---
output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "orbital_energy"))
mkpath(output_dir)
savefig(plt1, joinpath(output_dir, "orbital_energy_variation_$(dir_tag).png"))
savefig(plt2, joinpath(output_dir, "orbital_energy_fractional_drift_$(dir_tag).png"))
display(plt1); display(plt2)
println("Saved to: ", output_dir)


