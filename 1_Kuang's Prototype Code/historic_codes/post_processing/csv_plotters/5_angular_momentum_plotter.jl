using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

const MU      = 3.986004418e14
const C       = 3.0e8
const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/12_CSV_Write_Read.jl")

# --- Load all CSVs and plot ---
# data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h850km_i0.0deg"))
# data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h1000km_i0.0deg_nu-0.75deg"))
# data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h1000km_i0.0deg_nu1.50deg"))
dir_tag   = basename(data_dir)
csv_files = filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true))

show_variance = false   # ← set to false to hide the ±1σ shaded band

plt1 = plot(title="Net Satellite Angular Momentum ($dir_tag)", xlabel="t (s)",
            ylabel="Total Angular Momentum (kg·m²/s)", formatter=:scientific, legend=:outertopright)
plt2 = plot(title="Angular Momentum Fractional Drift ($dir_tag)", xlabel="t (s)",
            ylabel="ΔH/H(0) = (H(t) - H(0)) / H(0)", formatter=:scientific, legend=:outertopright)

# collect interpolated curves for mean computation
const N_PTS_H = 2000
all_t_H    = Vector{Vector{Float64}}()
all_Hmag   = Vector{Vector{Float64}}()
all_ΔH     = Vector{Vector{Float64}}()

for csv_path in csv_files
    sol, _, p = load_timeseries_csv(csv_path)
    masses = p[:masses]
    Hmag = [angular_momentum(u, masses)[2] for u in sol.u]
    ΔH   = (Hmag .- Hmag[1]) ./ Hmag[1]
    N = p[:N]

    t_uniform = collect(LinRange(sol.t[1], sol.t[end], N_PTS_H))
    function interp(vals, t_query)
        i = clamp(searchsortedlast(sol.t, t_query), 1, length(sol.t)-1)
        α = (t_query - sol.t[i]) / (sol.t[i+1] - sol.t[i])
        return vals[i] * (1-α) + vals[i+1] * α
    end
    Hmag_u = [interp(Hmag, t) for t in t_uniform]
    ΔH_u   = [interp(ΔH,   t) for t in t_uniform]

    push!(all_t_H,  t_uniform)
    push!(all_Hmag, Hmag_u)
    push!(all_ΔH,   ΔH_u)

    plot!(plt1, t_uniform, Hmag_u, label="N=$N", color=:grey, alpha=0.5, lw=1)
    plot!(plt2, t_uniform, ΔH_u,   label="N=$N", color=:grey, alpha=0.5, lw=1)
end

# mean line across all N cases
mean_Hmag = mean(stack(all_Hmag), dims=2) |> vec
mean_ΔH   = mean(stack(all_ΔH),   dims=2) |> vec
std_Hmag  = std(stack(all_Hmag),  dims=2) |> vec
std_ΔH    = std(stack(all_ΔH),    dims=2) |> vec

if show_variance
    plot!(plt1, all_t_H[1], mean_Hmag .+ std_Hmag, fillrange=mean_Hmag .- std_Hmag,
          fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
    plot!(plt2, all_t_H[1], mean_ΔH   .+ std_ΔH,   fillrange=mean_ΔH   .- std_ΔH,
          fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
end
plot!(plt1, all_t_H[1], mean_Hmag, label="Mean", color=:blue, lw=2)
plot!(plt2, all_t_H[1], mean_ΔH,   label="Mean", color=:blue, lw=2)

# --- Save ---
output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "angular_momentum"))
mkpath(output_dir)
savefig(plt1, joinpath(output_dir, "angular_momentum_$(dir_tag).png"))
savefig(plt2, joinpath(output_dir, "angular_momentum_fractional_drift_$(dir_tag).png"))
display(plt1); display(plt2)
println("Saved to: ", output_dir)


