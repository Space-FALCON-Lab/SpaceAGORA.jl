using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

const MU      = 3.986004418e14
const C       = 3.0e8
const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/12_CSV_Write_Read.jl")

# --- Load all CSVs and plot ---
#data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h850km_i0.0deg"))
# data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h1000km_i0.0deg_nu-0.75deg"))
data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h1000km_i0.0deg_nu1.5deg"))
#data_dir  = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "loop_for_num_of_helpers_from_main_29_test_850km_set/target_h850km_i0.0deg"))
dir_tag   = basename(data_dir)
csv_files = filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true))

show_variance = true   # ← set to false to hide the ±1σ shaded band

plt1 = plot(title="Net Satellite Momentum ($dir_tag)", xlabel="t (s)",
            ylabel="Total Momentum (kg·m/s)", formatter=:scientific, legend=:outertopright)
plt2 = plot(title="Linear Momentum Fractional Drift ($dir_tag)", xlabel="t (s)",
            ylabel="ΔP/P(0) = (P(t) - P(0)) / P(0)", formatter=:scientific, legend=:outertopright)

# collect interpolated curves for mean computation
const N_PTS  = 2000
all_t    = Vector{Vector{Float64}}()
all_Pmag = Vector{Vector{Float64}}()
all_ΔP   = Vector{Vector{Float64}}()

for csv_path in csv_files
    sol, _, p = load_timeseries_csv(csv_path)
    masses = p[:masses]
    Pmag = [total_momentum(u, masses)[2] for u in sol.u]
    ΔP   = (Pmag .- Pmag[1]) ./ Pmag[1]
    N = p[:N]

    t_uniform = collect(LinRange(sol.t[1], sol.t[end], N_PTS))
    function interp(vals, t_query)
        i = clamp(searchsortedlast(sol.t, t_query), 1, length(sol.t)-1)
        α = (t_query - sol.t[i]) / (sol.t[i+1] - sol.t[i])
        return vals[i] * (1-α) + vals[i+1] * α
    end
    Pmag_u = [interp(Pmag, t) for t in t_uniform]
    ΔP_u   = [interp(ΔP,   t) for t in t_uniform]

    push!(all_t,    t_uniform)
    push!(all_Pmag, Pmag_u)
    push!(all_ΔP,   ΔP_u)

    plot!(plt1, t_uniform, Pmag_u, label="N=$N", color=:grey, alpha=0.5, lw=1)
    plot!(plt2, t_uniform, ΔP_u,   label="N=$N", color=:grey, alpha=0.5, lw=1)
end

# mean line across all N cases
mean_Pmag = mean(stack(all_Pmag), dims=2) |> vec
mean_ΔP   = mean(stack(all_ΔP),   dims=2) |> vec
std_Pmag  = std(stack(all_Pmag),  dims=2) |> vec
std_ΔP    = std(stack(all_ΔP),    dims=2) |> vec

if show_variance
    plot!(plt1, all_t[1], mean_Pmag .+ std_Pmag, fillrange=mean_Pmag .- std_Pmag,
          fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
    plot!(plt2, all_t[1], mean_ΔP   .+ std_ΔP,   fillrange=mean_ΔP   .- std_ΔP,
          fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
end
plot!(plt1, all_t[1], mean_Pmag, label="Mean", color=:blue, lw=2)
plot!(plt2, all_t[1], mean_ΔP,   label="Mean", color=:blue, lw=2)

# --- Save ---
output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "linear_momentum"))
mkpath(output_dir)
savefig(plt1, joinpath(output_dir, "linear_momentum_$(dir_tag).png"))
savefig(plt2, joinpath(output_dir, "linear_momentum_fractional_drift_$(dir_tag).png"))
display(plt1); display(plt2)
println("Saved to: ", output_dir)


