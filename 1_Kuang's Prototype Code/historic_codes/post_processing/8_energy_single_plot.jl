# Single plot — Orbital Energy Fractional Drift ΔE/|E₀|
# Set the tuning knobs below and run to plot one figure.

using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

const MU = 3.986004418e14; const C = 3.0e8; const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include(joinpath(@__DIR__, "..", "functions", "4_Diagnostics.jl"))
include(joinpath(@__DIR__, "..", "functions", "12_CSV_Write_Read.jl"))

# ── Tuning Knobs ──────────────────────────────────────────────────────────────
target_altitudes_km     = 950 #[850, 950, 1000, 1050, 1150]
target_inclinations_deg = 1.0 #[0.0, 0.5, 1.0, 1.5]
filter_N         = 100   # helper count, e.g. 101  (matches _N101_ in filename)
filter_B         = 2000   # e.g. 100.0
filter_Pin       = 1e5   # e.g. 1e4
filter_min_range = 0.0   # set both or neither, e.g. 0.0
filter_max_range = 3000e3   # e.g. Inf

# ── Layout Knobs ──────────────────────────────────────────────────────────────
gap_inner           = 0;  margin_left = 0;  margin_bottom = 0
margin_right        = 0.5; margin_top = 0
tick_fontsize       = 8;  maintitle_fontsize = 12;  main_label_fontsize = 9
title_fontsize      = 12
subtitle_fontsize   = 9
title_top_margin_mm = 6     # mm of space above the axes (pushes title up, subtitle down)
subtitle_drop_frac  = -0.08  # subtitle drops this fraction of the y-range below the top edge (0 = flush top, 0.1 = 10% down)
row_label_fontsize  = 8;  col_label_fontsize = 8
main_xlabel         = "t (hr)";  main_ylabel = "ΔE/|E₀|  [×10⁻⁶]"
ylabel_gap          = 0.5;  xlabel_gap = 0.5
# ─────────────────────────────────────────────────────────────────────────────

# Build the list of filename fragments that every accepted CSV must contain
csv_filters = String[]
filter_N         !== nothing && push!(csv_filters, @sprintf("_N%d_",         filter_N+1))
filter_B         !== nothing && push!(csv_filters, @sprintf("_B%.4g",        filter_B))
filter_Pin       !== nothing && push!(csv_filters, @sprintf("_Pin%.4g",      filter_Pin))
(filter_min_range !== nothing && filter_max_range !== nothing) &&
    push!(csv_filters, @sprintf("_rmin%.0fm_rmax%.4g", filter_min_range, filter_max_range))

# Linear interpolation helper
function interp1(xs, ys, xq)
    ii = clamp(searchsortedlast(xs, xq), 1, length(xs)-1)
    α  = (xq - xs[ii]) / (xs[ii+1] - xs[ii])
    ys[ii] * (1-α) + ys[ii+1] * α
end


h     = target_altitudes_km
i_deg = target_inclinations_deg
data_dir = normpath(joinpath(@__DIR__, "..", "output", "CSV",
                             @sprintf("target_h%dkm_i%.1fdeg", h, i_deg)))

csv_files = isdir(data_dir) ?
    filter(f -> endswith(f,".csv") && all(fr->occursin(fr,basename(f)), csv_filters),
           readdir(data_dir; join=true)) : String[]

isempty(csv_files) && error("No matching CSVs found in: $data_dir  (filters: $csv_filters)")

all_ΔE   = Vector{Vector{Float64}}()
t_common = nothing
for csv_path in csv_files
    sol, _, p = try load_timeseries_csv(csv_path) catch e; @warn e; continue end
    Etot  = [sum(orbital_energy(u, p[:masses], MU)) for u in sol.u]
    ΔE    = (Etot .- Etot[1]) ./ abs(Etot[1])
    t_uni = collect(LinRange(sol.t[1], sol.t[end], 2000))
    push!(all_ΔE, [interp1(sol.t, ΔE, t) for t in t_uni])
    isnothing(t_common) && (global t_common = t_uni)
end

isempty(all_ΔE) && error("Failed to load any CSVs.")

t_hr    = t_common ./ 3600.0
mat     = hcat(all_ΔE...)
mean_ΔE = vec(mean(mat, dims=2))
std_ΔE  = vec(std(mat,  dims=2))

_sub_parts = [@sprintf("h=%dkm, i=%.1f°", h, i_deg)]
filter_N         !== nothing && push!(_sub_parts, "N=$(filter_N+1)")
filter_B         !== nothing && push!(_sub_parts, "B=$(filter_B)")
filter_Pin       !== nothing && push!(_sub_parts, "Pin=$(filter_Pin)")
filter_max_range !== nothing &&
    push!(_sub_parts, "Range=$(filter_max_range/1000) km")
_subtitle_str = join(_sub_parts, "  |  ")

fig = plot(legend=false, grid=true, tick_direction=:out,
           xlabel="t (hr)", ylabel="ΔE/|E₀|  [×10⁻⁶]",
           title="Orbital Energy Fractional Drift",
           titlefontsize=title_fontsize,
           top_margin=title_top_margin_mm*Plots.mm)
for c in all_ΔE
    plot!(fig, t_hr, c .* 1e6, color=:grey, alpha=0.4, lw=0.5, label="")
end
plot!(fig, t_hr, (mean_ΔE .+ std_ΔE) .* 1e6,
      fillrange=(mean_ΔE .- std_ΔE) .* 1e6,
      fillalpha=0.2, fillcolor=:blue, linealpha=0, label="")
plot!(fig, t_hr, mean_ΔE .* 1e6, color=:blue, lw=1.5, label="")
_yl = ylims(fig)
annotate!(fig, mean(xlims(fig)), _yl[2] - subtitle_drop_frac*(_yl[2]-_yl[1]),
          Plots.text(_subtitle_str, subtitle_fontsize, :center, :top))

output_dir = normpath(joinpath(@__DIR__, "..", "output", "images_reconstructed_from_csv", "orbital_energy"))
mkpath(output_dir)
out_path = joinpath(output_dir, @sprintf("energy_fractional_h%dkm_i%.1fdeg.png", h, i_deg))
savefig(fig, out_path)
display(fig)
println("Saved to: ", out_path)
