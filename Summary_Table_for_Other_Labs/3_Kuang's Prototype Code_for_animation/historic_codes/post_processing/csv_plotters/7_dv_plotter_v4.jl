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

N_helpers_vec = Int[]
slope_rows    = Vector{Vector{Float64}}()

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

slope_mat = stack(slope_rows)'   # (n_files × 3)

# --- Plot ---
labels = ["Δv_R", "Δv_T", "Δv_N"]
colors = [:blue, :orange, :green]

plts = [plot(title="d$(lbl)/dt vs N helpers ($dir_tag)",
             xlabel="Number of helpers", ylabel="dΔv/dt (m/s²)",
             legend=false, formatter=:scientific)
        for lbl in labels]

for (col, clr, plt) in zip(1:3, colors, plts)
    y = slope_mat[:, col]
    plot!(plt, N_helpers_vec, y; color=clr, marker=:circle, ms=5, lw=2)
end

fig = plot(plts...; layout=(3,1), size=(700, 900))
savefig(fig, joinpath(output_dir, "dvdt_RTN_vs_N_helpers_$(dir_tag).png"))
display(fig)
println("Saved to: ", output_dir)

