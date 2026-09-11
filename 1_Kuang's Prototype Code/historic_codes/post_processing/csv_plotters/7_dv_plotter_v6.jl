using LinearAlgebra, Plots, Printf, StaticArrays, DelimitedFiles

const MU      = 3.986004418e14
const C       = 3.0e8
const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/5_OE_Converters.jl")
include("../../functions/6_OE_and_dv_in_RTN.jl")
include("../../functions/12_CSV_Write_Read.jl")

# --- Config ---
data_dir   = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", "target_h850km_i0.0deg"))
dir_tag    = basename(data_dir)
csv_files  = sort(filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true)),
                  by = f -> parse(Int, match(r"_N(\d+)_", basename(f))[1]))
output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "delta_v"))
mkpath(output_dir)

# --- Compute dΔv/dt slope for each CSV ---
linfit(x, y) = (hcat(x, ones(length(x))) \ y)[1]   # returns slope only

cache_path    = joinpath(output_dir, "dvdt_slopes_$(dir_tag).csv")
N_helpers_vec = Int[]
slope_rows    = Vector{Vector{Float64}}()

if isfile(cache_path)
    raw = readdlm(cache_path, ',', Float64; skipstart=1)
    N_helpers_vec = Int.(raw[:, 1])
    slope_rows    = [raw[i, 2:4] for i in 1:size(raw, 1)]
    println("Loaded slopes from cache: $cache_path")
else
    for csv_path in csv_files
        sol, _, p  = load_timeseries_csv(csv_path)
        t, Δv_hist = delta_v_RTN_time_series(sol, p)
        dv         = Δv_hist[p[:N]]          # target = last satellite
        N_h        = p[:N] - 1               # number of helpers
        slopes     = [linfit(t, dv[row, :]) for row in 1:3]
        push!(N_helpers_vec, N_h)
        push!(slope_rows, slopes)
        @printf("  N_helpers=%d  dΔv/dt: R=%+.4e  T=%+.4e  N=%+.4e  m/s²\n",
                N_h, slopes...)
    end
    open(cache_path, "w") do io
        writedlm(io, ["N_helpers" "slope_R" "slope_T" "slope_N"], ',')
        writedlm(io, hcat(N_helpers_vec, stack(slope_rows)'), ',')
    end
    println("Saved slopes to cache: $cache_path")
end

slope_mat = stack(slope_rows)'   # (n_files × 3)

# --- Plot ---
# order: R (blue dash), T (green solid), N (orange dot) — N last so it renders on top
styles = [(:blue, :dash), (:green, :solid), (:orange, :dot)]
labels = ["Δv_R", "Δv_T", "Δv_N"]
cols   = [1, 2, 3]

plt = plot(title="dΔv/dt (RTN) vs Number of Helpers ($dir_tag)",
           xlabel="Number of helpers", ylabel="dΔv/dt (m/s²)",
           legend=:outertopright, formatter=:scientific)

for (col, lbl, (clr, ls)) in zip(cols, labels, styles)
    plot!(plt, N_helpers_vec, slope_mat[:, col]; label=lbl, color=clr, linestyle=ls, lw=2)
end

savefig(plt, joinpath(output_dir, "dvdt_RTN_vs_N_helpers_$(dir_tag).png"))
display(plt)
println("Saved to: ", output_dir)

