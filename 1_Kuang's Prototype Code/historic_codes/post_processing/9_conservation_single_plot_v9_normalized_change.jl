# 3×1 combined plot — Linear Momentum / Angular Momentum / Orbital Energy
# Fractional drift ΔP/P₀, ΔH/H₀ [×10⁻¹¹], ΔE/|E₀| [×10⁻⁶]
# Set the tuning knobs below and run.

using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

const MU = 3.986004418e14; const C = 3.0e8; const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include(joinpath(@__DIR__, "..", "functions", "4_Diagnostics.jl"))
include(joinpath(@__DIR__, "..", "functions", "12_CSV_Write_Read.jl"))
const Ẑ = SVector(0.0, 0.0, 1.0)   # ECI z-axis (North Pole), required by rv2coe
include(joinpath(@__DIR__, "..", "functions", "5_OE_Converters.jl"))

# ── Tuning Knobs ──────────────────────────────────────────────────────────────
target_altitudes_km     = 1100 #[850, 950, 1000, 1050, 1150]
target_inclinations_deg = 0.0 #[0.0, 0.5, 1.0, 1.5]
filter_N         = 1   # helper count (nothing = all)
filter_B         = 100.0   # e.g. 100.0
filter_Pin       = 1e4   # e.g. 1e4
filter_min_range = 0.0   # set both or neither
filter_max_range = 2e5
filter_J2        = false  # true = only _J2T files, false = only _J2F files, nothing = all
subtitle_fontsize   = 10
xlabel_fontsize     = 9   # x-axis label font size
ylabel_fontsize     = 9   # y-axis label font size
title_top_margin_mm = 5     # mm of space above axes of top panel
subtitle_drop_frac  = -0.1 # subtitle y-offset from top edge (0 = flush, 0.1 = 10% down)
panel_top_margin_mm = 3     # mm top margin for middle and bottom panels
# ─────────────────────────────────────────────────────────────────────────────

# Build filename filter fragments
csv_filters = String[]
filter_N         !== nothing && push!(csv_filters, @sprintf("_N%d_",         filter_N+1))
filter_B         !== nothing && push!(csv_filters, @sprintf("_B%.4g",        filter_B))
filter_Pin       !== nothing && push!(csv_filters, @sprintf("_Pin%.4g",      filter_Pin))
(filter_min_range !== nothing && filter_max_range !== nothing) &&
    push!(csv_filters, @sprintf("_rmin%.0fm_rmax%.4g", filter_min_range, filter_max_range))
filter_J2 === true  && push!(csv_filters, "_J2T")
filter_J2 === false && push!(csv_filters, "_J2F")

# Linear interpolation helper
function interp1(xs, ys, xq)
    ii = clamp(searchsortedlast(xs, xq), 1, length(xs)-1)
    α  = (xq - xs[ii]) / (xs[ii+1] - xs[ii])
    ys[ii] * (1-α) + ys[ii+1] * α
end

# Compute cumulative orbit count for the TARGET satellite (last satellite in state vector).
# Uses the instantaneous semi-major axis a(t) at each saved time step, so it correctly
# accounts for orbital period changes caused by the laser.
# orbit_count(t) = ∫₀ᵗ dt'/T(t')  ≈  cumsum(Δt / T(tₖ))
function orbit_count_from_sol(sol)
    N_sat = length(sol.u[1]) ÷ 6          # total satellites; target = last
    # compute instantaneous period T(t) at each saved time point
    T_series = [begin
        u = sol.u[k]
        r = @SVector [u[idx(N_sat,1)], u[idx(N_sat,2)], u[idx(N_sat,3)]]
        v = @SVector [u[idx(N_sat,4)], u[idx(N_sat,5)], u[idx(N_sat,6)]]
        a = rv2coe(r, v, MU).a
        2π * sqrt(a^3 / MU)
    end for k in eachindex(sol.t)]
    # integrate: use mid-point rule over each interval
    dt = diff(sol.t)
    T_mid = (T_series[1:end-1] .+ T_series[2:end]) ./ 2
    orbit_increments = dt ./ T_mid
    return cumsum([0.0; orbit_increments])
end

h     = target_altitudes_km
i_deg = target_inclinations_deg
data_dir = normpath(joinpath(@__DIR__, "..", "output", "CSV",
                             @sprintf("target_h%dkm_i%.1fdeg", h, i_deg)))

csv_files = isdir(data_dir) ?
    filter(f -> endswith(f,".csv") && all(fr->occursin(fr,basename(f)), csv_filters) &&
               (filter_J2 !== nothing || (!occursin("_J2T", basename(f)) && !occursin("_J2F", basename(f)))),
           readdir(data_dir; join=true)) : String[]

isempty(csv_files) && error("No matching CSVs found in: $data_dir  (filters: $csv_filters)")

# ── Load all CSVs once, compute all three quantities ─────────────────────────
all_ΔP   = Vector{Vector{Float64}}()
all_ΔH   = Vector{Vector{Float64}}()
all_ΔE   = Vector{Vector{Float64}}()
max_orbits_all = Float64[]   # track max orbit count per CSV for common grid
_loaded = []  # (sol.t, orbit_count, ΔP, ΔH, ΔE) per CSV
for csv_path in csv_files
    sol, _, p = try load_timeseries_csv(csv_path) catch e; @warn e; continue end

    Pmag = [total_momentum(u, p[:masses])[2]       for u in sol.u]
    Hmag = [angular_momentum(u, p[:masses])[2]     for u in sol.u]
    Etot = [sum(orbital_energy(u, p[:masses], MU)) for u in sol.u]

    ΔP = (Pmag .- Pmag[1]) ./ Pmag[1]
    ΔH = (Hmag .- Hmag[1]) ./ Hmag[1]
    ΔE = (Etot .- Etot[1]) ./ abs(Etot[1])

    oc = orbit_count_from_sol(sol)   # orbit count at each sol.t sample
    push!(max_orbits_all, oc[end])
    push!(_loaded, (t=sol.t, oc=oc, ΔP=ΔP, ΔH=ΔH, ΔE=ΔE))
end

isempty(_loaded) && error("Failed to load any CSVs.")

# Common orbit-count grid: 0 → minimum of all max orbit counts (no extrapolation)
orbit_common = collect(LinRange(0.0, minimum(max_orbits_all), 2000))

for d in _loaded
    push!(all_ΔP, [interp1(d.oc, d.ΔP, o) for o in orbit_common])
    push!(all_ΔH, [interp1(d.oc, d.ΔH, o) for o in orbit_common])
    push!(all_ΔE, [interp1(d.oc, d.ΔE, o) for o in orbit_common])
end

isempty(all_ΔP) && error("Failed to load any CSVs.")

function _stats(all_curves)
    mat = hcat(all_curves...)
    vec(mean(mat, dims=2)), vec(std(mat, dims=2))
end

# Auto-scale: returns (scale, ylabel_suffix) so axis ticks sit in a tidy range.
# base_scale=100 converts fractional → percentage before computing the exponent.
function _auto_scale(mean_curve; base_scale=100.0)
    max_val = maximum(abs, mean_curve) * base_scale
    max_val == 0 && return base_scale, " [%]"
    n = floor(Int, log10(max_val))
    scale = base_scale * 10.0^(-n)
    _sup = Dict('0'=>'⁰','1'=>'¹','2'=>'²','3'=>'³','4'=>'⁴',
                '5'=>'⁵','6'=>'⁶','7'=>'⁷','8'=>'⁸','9'=>'⁹','-'=>'⁻')
    exp_str = join([_sup[c] for c in string(n)])
    return scale, " [×10$(exp_str) %]"
end

mean_ΔP, std_ΔP = _stats(all_ΔP)
mean_ΔH, std_ΔH = _stats(all_ΔH)
mean_ΔE, std_ΔE = _stats(all_ΔE)

scale_P, suffix_P = _auto_scale(mean_ΔP)
scale_H, suffix_H = _auto_scale(mean_ΔH)
scale_E, suffix_E = _auto_scale(mean_ΔE)

x_axis = orbit_common   # use orbit count as x-axis everywhere

# ── Build subtitle string ─────────────────────────────────────────────────────
_sub_parts = [@sprintf("h=%dkm, i=%.1f°", h, i_deg)]
filter_N         !== nothing && push!(_sub_parts, "N=$(filter_N+1)")
filter_B         !== nothing && push!(_sub_parts, "B=$(filter_B)")
filter_Pin       !== nothing && push!(_sub_parts, "Pin=$(filter_Pin)")
filter_max_range !== nothing && push!(_sub_parts, "Range=$(filter_max_range/1000) km")
filter_J2 === true  && push!(_sub_parts, "J2=T")
filter_J2 === false && push!(_sub_parts, "J2=F")
_subtitle_str = join(_sub_parts, "  |  ")

# ── Helper: draw grey curves + mean ± std band on a subplot ──────────────────
function _draw_panel!(plt, t_hr, all_c, mean_c, std_c; scale=1.0)
    for c in all_c
        plot!(plt, t_hr, c .* scale, color=:grey, alpha=0.4, lw=0.5, label="")
    end
    plot!(plt, t_hr, (mean_c .+ std_c) .* scale,
          fillrange=(mean_c .- std_c) .* scale,
          fillalpha=0.2, fillcolor=:blue, linealpha=0, label="")
    plot!(plt, t_hr, mean_c .* scale, color=:blue, lw=1.5, label="")
end

# ── Panel 1: Linear Momentum ──────────────────────────────────────────────────
p1 = plot(legend=false, grid=true, tick_direction=:out,
          xlabel="Number of orbits (target satellite)", ylabel="ΔP/P₀$(suffix_P)",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          # title="Linear Momentum, Angular Momentum and Orbital Energy",
          title="Linear Momentum Change (Normalized)",
          titlefontsize=title_fontsize,
          top_margin=title_top_margin_mm*Plots.mm)
_draw_panel!(p1, x_axis, all_ΔP, mean_ΔP, std_ΔP; scale=scale_P)
_yl = ylims(p1)
annotate!(p1, mean(xlims(p1)), _yl[2] - subtitle_drop_frac*(_yl[2]-_yl[1]),
          Plots.text(_subtitle_str, subtitle_fontsize, :center, :top))

# ── Panel 2: Angular Momentum ─────────────────────────────────────────────────
p2 = plot(legend=false, grid=true, tick_direction=:out,
          xlabel="Number of orbits (target satellite)", ylabel="ΔH/H₀$(suffix_H)",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          title="Angular Momentum Change (Normalized)",
          titlefontsize=title_fontsize,
          top_margin=panel_top_margin_mm*Plots.mm)
_draw_panel!(p2, x_axis, all_ΔH, mean_ΔH, std_ΔH; scale=scale_H)
_yl2 = ylims(p2)
annotate!(p2, mean(xlims(p2)), _yl2[2] - subtitle_drop_frac*(_yl2[2]-_yl2[1]),
          Plots.text(_subtitle_str, subtitle_fontsize, :center, :top))

# ── Panel 3: Orbital Energy ───────────────────────────────────────────────────
p3 = plot(legend=false, grid=true, tick_direction=:out,
           xlabel="Number of orbits (target satellite)", ylabel="ΔE/|E₀|$(suffix_E)",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          title="Orbital Energy Change (Normalized)",
          titlefontsize=title_fontsize,
          top_margin=panel_top_margin_mm*Plots.mm)
_draw_panel!(p3, x_axis, all_ΔE, mean_ΔE, std_ΔE; scale=scale_E)
_yl3 = ylims(p3)
annotate!(p3, mean(xlims(p3)), _yl3[2] - subtitle_drop_frac*(_yl3[2]-_yl3[1]),
          Plots.text(_subtitle_str, subtitle_fontsize, :center, :top))

# ── Combine and save ──────────────────────────────────────────────────────────
fig = plot(p1, p2, p3, layout=grid(3, 1, heights=panel_heights), size=fig_size)

output_dir = normpath(joinpath(@__DIR__, "..", "output", "images_reconstructed_from_csv", "conservation"))
mkpath(output_dir)
out_path     = joinpath(output_dir, @sprintf("conservation_h%dkm_i%.1fdeg.png", h, i_deg))
out_path_pdf = joinpath(output_dir, @sprintf("conservation_h%dkm_i%.1fdeg.pdf", h, i_deg))
savefig(fig, out_path)
savefig(fig, out_path_pdf)
display(fig)
println("Saved to: ", out_path)
println("Saved to: ", out_path_pdf)
