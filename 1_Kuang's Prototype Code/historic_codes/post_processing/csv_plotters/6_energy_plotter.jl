using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

const MU      = 3.986004418e14
const C       = 3.0e8
const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/12_CSV_Write_Read.jl")


# --- Load all CSVs and plot ---
# data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h850km_i0.0deg"))
data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h1000km_i0.0deg_nu-0.75deg"))
# data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h1000km_i0.0deg_nu1.50deg"))
dir_tag   = basename(data_dir)
csv_files = filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true))

show_variance = false   # ← set to false to hide the ±1σ shaded band

plt1 = plot(title="Satellite Orbital Energy Variation ($dir_tag)", xlabel="t (s)",
            ylabel="ΔE = E(t) - E(0) (J)", formatter=:scientific, legend=:outertopright)
plt2 = plot(title="Orbital Energy Fractional Drift ($dir_tag)", xlabel="t (s)",
            ylabel="ΔE/|E(0)| = (E(t) - E(0)) / |E(0)|", formatter=:scientific, legend=:outertopright)

# collect interpolated curves for mean computation
const N_PTS_E = 2000
all_t_E   = Vector{Vector{Float64}}()
all_ΔE    = Vector{Vector{Float64}}()
all_ΔE_fr = Vector{Vector{Float64}}()

for csv_path in csv_files
    sol, _, p = load_timeseries_csv(csv_path)
    masses = p[:masses]
    Etot = [sum(orbital_energy(u, masses, MU)) for u in sol.u]
    ΔE      = Etot .- Etot[1]
    ΔE_frac = ΔE ./ abs(Etot[1])
    N = p[:N]

    t_uniform = collect(LinRange(sol.t[1], sol.t[end], N_PTS_E))
    function interp(vals, t_query)
        i = clamp(searchsortedlast(sol.t, t_query), 1, length(sol.t)-1)
        α = (t_query - sol.t[i]) / (sol.t[i+1] - sol.t[i])
        return vals[i] * (1-α) + vals[i+1] * α
    end
    ΔE_u      = [interp(ΔE,      t) for t in t_uniform]
    ΔE_frac_u = [interp(ΔE_frac, t) for t in t_uniform]

    push!(all_t_E,   t_uniform)
    push!(all_ΔE,    ΔE_u)
    push!(all_ΔE_fr, ΔE_frac_u)

    plot!(plt1, t_uniform, ΔE_u,      label="N=$N", color=:grey, alpha=0.5, lw=1)
    plot!(plt2, t_uniform, ΔE_frac_u, label="N=$N", color=:grey, alpha=0.5, lw=1)
end

# mean and variance across all N cases
mean_ΔE    = mean(stack(all_ΔE),    dims=2) |> vec
mean_ΔE_fr = mean(stack(all_ΔE_fr), dims=2) |> vec
std_ΔE     = std(stack(all_ΔE),     dims=2) |> vec
std_ΔE_fr  = std(stack(all_ΔE_fr),  dims=2) |> vec

if show_variance
    plot!(plt1, all_t_E[1], mean_ΔE    .+ std_ΔE,    fillrange=mean_ΔE    .- std_ΔE,
          fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
    plot!(plt2, all_t_E[1], mean_ΔE_fr .+ std_ΔE_fr, fillrange=mean_ΔE_fr .- std_ΔE_fr,
          fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
end
plot!(plt1, all_t_E[1], mean_ΔE,    label="Mean", color=:blue, lw=2)
plot!(plt2, all_t_E[1], mean_ΔE_fr, label="Mean", color=:blue, lw=2)

# --- Save ---
output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "orbital_energy"))
mkpath(output_dir)
savefig(plt1, joinpath(output_dir, "orbital_energy_variation_$(dir_tag).png"))
savefig(plt2, joinpath(output_dir, "orbital_energy_fractional_drift_$(dir_tag).png"))
display(plt1); display(plt2)
println("Saved to: ", output_dir)


