using LinearAlgebra, Plots, Printf, StaticArrays, DelimitedFiles

const MU      = 3.986004418e14
const C       = 3.0e8
const R_EARTH = 6_378_137.0
const R_ATMDEF = R_EARTH + 100_000.0    # default "atmosphere radius" (Kármán line) [m]
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/5_OE_Converters.jl")
include("../../functions/2_Laser_Forces_ver2.jl")
include("../../functions/6_OE_and_dv_in_RTN.jl")
include("../../functions/12_CSV_Write_Read.jl")

# --- Config ---
csv_path   = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV",
                               "target_h850km_i0.0deg",
                               "timeseries_N101_T540000s_h1000km_t850km_ih0.0deg_it0.0deg.csv"))
output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "delta_v"))
mkpath(output_dir)

# --- Load & compute ---
sol, _, p  = load_timeseries_csv(csv_path)
csv_tag    = replace(basename(csv_path), ".csv" => "")
t, Δv_hist = delta_v_RTN_time_series(sol, p)
dv         = Δv_hist[p[:N]]   # target = last satellite

# --- Linear fit: dv ≈ slope*t + intercept ---
labels   = ["Δv_R", "Δv_T", "Δv_N"]
suffixes = ["_R_only", "_T_only", "_N_only"]

A    = hcat(t, ones(length(t)))              # [t | 1]
fits = [A \ dv[row, :] for row in 1:3]      # fits[row] = [slope, intercept]

println("\nLinear fit rates of change (Δv = slope·t + intercept):")
for (lbl, f) in zip(labels, fits)
    @printf("  %-6s  slope = %+.6e m/s²   intercept = %+.6e m/s\n", lbl, f[1], f[2])
end

# --- Plot ---

# individual components (data + linear fit overlay)
for (row, lbl, suf, f) in zip(1:3, labels, suffixes, fits)
    fit_line = f[1] .* t .+ f[2]
    plt = plot(t, dv[row, :]; label=lbl, xlabel="t (s)", ylabel="m/s",
               title="$lbl component (RTN)", color=:grey, alpha=0.7, lw=1)
    plot!(plt, t, fit_line; label="fit ($(Printf.@sprintf("%+.3e", f[1])) m/s²)",
          color=:red, lw=2, linestyle=:dash)
    savefig(plt, joinpath(output_dir, "dv_RTN_$(csv_tag)$(suf).png"))
end

# all components together
plt_all = plot(t, dv[1, :]; label="Δv_R", xlabel="t (s)", ylabel="m/s",
               title="Δv components (RTN)")
plot!(plt_all, t, dv[2, :]; label="Δv_T")
plot!(plt_all, t, dv[3, :]; label="Δv_N")
for (row, lbl, f) in zip(1:3, labels, fits)
    plot!(plt_all, t, f[1] .* t .+ f[2];
          label="$(lbl) fit", lw=2, linestyle=:dash)
end
savefig(plt_all, joinpath(output_dir, "dv_RTN_$(csv_tag).png"))

display(plt_all)
println("Saved to: ", output_dir)

