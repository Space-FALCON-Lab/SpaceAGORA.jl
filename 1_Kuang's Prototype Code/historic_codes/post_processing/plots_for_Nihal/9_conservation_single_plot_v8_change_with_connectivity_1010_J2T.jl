# 3×1 combined plot — Linear Momentum / Angular Momentum / Orbital Energy
# Fractional drift ΔP/P₀, ΔH/H₀ [×10⁻¹¹], ΔE/|E₀| [×10⁻⁶]
# Set the tuning knobs below and run.

using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

const MU = 3.986004418e14; const C = 3.0e8; const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include(joinpath(@__DIR__, "..", "..", "functions", "4_Diagnostics.jl"))
include(joinpath(@__DIR__, "..", "..", "functions", "12_CSV_Write_Read.jl"))
const Ẑ = SVector(0.0, 0.0, 1.0)   # ECI z-axis (North Pole), required by rv2coe
include(joinpath(@__DIR__, "..", "..", "functions", "5_OE_Converters.jl"))

# ── Tuning Knobs ──────────────────────────────────────────────────────────────
target_altitudes_km     = 1010 #[850, 950, 1000, 1050, 1150]
target_inclinations_deg = 0.0 #[0.0, 0.5, 1.0, 1.5]
filter_N         = 1   # helper count (nothing = all)
filter_B         = 100.0   # e.g. 100.0
filter_Pin       = 1e4   # e.g. 1e4
filter_min_range = 0.0   # set both or neither
filter_max_range = 2e5
filter_J2        = true  # true = only _J2T files, false = only _J2F files, nothing = all
range_threshold_m = 200000  # threshold for "in range" (meters)
title_fontsize   = 10
subtitle_fontsize   = 9
xlabel_fontsize     = 9   # x-axis label font size
ylabel_fontsize     = 9   # y-axis label font size
title_top_margin_mm = 3     # mm of space above axes of top panel
subtitle_drop_frac  = -0.17 # subtitle y-offset from top edge (0 = flush, 0.1 = 10% down)
panel_top_margin_mm = 3     # mm top margin for middle and bottom panels
panel_heights = [1/4, 1/4, 1/4, 1/4]  # relative heights for the four panels
fig_size = (850, 850)      # figure size in pixels
annot_dx_p1 =  -0.03;  annot_dy_p1 =  0.05   # ΔP panel label offset (fraction of axis range)
annot_dx_p2 =  0.01;  annot_dy_p2 =  0.04   # ΔH panel label offset
annot_dx_p3 =  -0.03;  annot_dy_p3 =  0.05   # Δε panel label offset
grid_alpha   = 0.4   # grid line darkness (0=invisible, 1=solid)
# ─────────────────────────────────────────────────────────────────────────────
Plots.default(gridalpha=grid_alpha)

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
data_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV",
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
all_in_range = Vector{Vector{Float64}}()
max_orbits_all = Float64[]   # track max orbit count per CSV for common grid
_loaded = []  # (sol.t, orbit_count, ΔP, ΔH, ΔE, in_range) per CSV
for csv_path in csv_files
    sol, _, p = try load_timeseries_csv(csv_path) catch e; @warn e; continue end

    Pmag = [total_momentum(u, p[:masses])[2]       for u in sol.u]
    Hmag = [angular_momentum(u, p[:masses])[2]     for u in sol.u]
    Etot = [sum(orbital_energy(u, p[:masses], MU)) for u in sol.u]

    ΔP = (Pmag .- Pmag[1])
    ΔH = (Hmag .- Hmag[1])
    ΔE = (Etot .- Etot[1])

    N_sat = length(sol.u[1]) ÷ 6
    in_range_series = [begin
        u = sol.u[k]
        r_target = @SVector [u[idx(N_sat,1)], u[idx(N_sat,2)], u[idx(N_sat,3)]]
        min_dist = minimum(1:N_sat-1) do i
            r_helper = @SVector [u[idx(i,1)], u[idx(i,2)], u[idx(i,3)]]
            norm(r_target - r_helper)
        end
        min_dist < range_threshold_m ? 1.0 : 0.0
    end for k in eachindex(sol.t)]

    oc = orbit_count_from_sol(sol)   # orbit count at each sol.t sample
    push!(max_orbits_all, oc[end])
    push!(_loaded, (t=sol.t, oc=oc, ΔP=ΔP, ΔH=ΔH, ΔE=ΔE, in_range=in_range_series))
end

isempty(_loaded) && error("Failed to load any CSVs.")

# Common orbit-count grid: 0 → minimum of all max orbit counts (no extrapolation)
orbit_common = collect(LinRange(0.0, minimum(max_orbits_all), 2000))

for d in _loaded
    push!(all_ΔP, [interp1(d.oc, d.ΔP, o) for o in orbit_common])
    push!(all_ΔH, [interp1(d.oc, d.ΔH, o) for o in orbit_common])
    push!(all_ΔE, [interp1(d.oc, d.ΔE, o) for o in orbit_common])
    push!(all_in_range, [interp1(d.oc, d.in_range, o) for o in orbit_common])
end

isempty(all_ΔP) && error("Failed to load any CSVs.")

function _stats(all_curves)
    mat = hcat(all_curves...)
    vec(mean(mat, dims=2)), vec(std(mat, dims=2))
end

mean_ΔP, std_ΔP = _stats(all_ΔP)
mean_ΔH, std_ΔH = _stats(all_ΔH)
mean_ΔE, std_ΔE = _stats(all_ΔE)
mean_in_range, std_in_range = _stats(all_in_range)

x_axis = orbit_common   # use orbit count as x-axis everywhere

# Parse the CSV filename to extract parameters for subtitle
function parse_csv_filename(filename)
    # filename like "timeseries_N2_T7200s_h1000km_t1010km_ih0.0deg_it0.0deg.csv"
    parts = split(basename(filename), '_')
    # Remove .csv from last part
    parts[end] = replace(parts[end], ".csv" => "")
    N = parse(Int, replace(parts[2], "N" => ""))
    T = parse(Float64, replace(parts[3], "T" => "", "s" => ""))
    h_helper = parse(Float64, replace(parts[4], "h" => "", "km" => ""))
    h_target = parse(Float64, replace(parts[5], "t" => "", "km" => ""))
    i_helper = parse(Float64, replace(parts[6], "ih" => "", "deg" => ""))
    i_target = parse(Float64, replace(parts[7], "it" => "", "deg" => ""))
    return N, T, h_helper, h_target, i_helper, i_target
end

# Assume single CSV for subtitle
if length(csv_files) == 1
    N_csv, T_csv, h_h, h_t, i_h, i_t = parse_csv_filename(csv_files[1])
else
    # Fallback to tuning knobs if multiple CSVs
    N_csv = filter_N !== nothing ? filter_N + 1 : nothing
    h_t = h
    i_t = i_deg
    h_h = nothing
    i_h = nothing
    T_csv = nothing
end

# ── Build subtitle string ─────────────────────────────────────────────────────
_sub_parts = [@sprintf("h=%dkm, i=%.1f°", h_t, i_t)]
N_csv !== nothing && push!(_sub_parts, "N=$N_csv")
h_h !== nothing && push!(_sub_parts, @sprintf("h_helper=%dkm", h_h))
i_h !== nothing && push!(_sub_parts, @sprintf("i_helper=%.1f°", i_h))
T_csv !== nothing && push!(_sub_parts, @sprintf("T=%.0fs", T_csv))
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

# Fixed-point formatter: 3 sig figs, never scientific notation, Unicode minus
_fmt(x) = x == 0.0 ? "0" :
    replace(@sprintf("%.*f", max(0, 3 - floor(Int, log10(abs(x))) - 1), x), "-" => "\u2212")

# ── Panel 1: Linear Momentum ──────────────────────────────────────────────────
p1 = plot(legend=false, grid=true, tick_direction=:out,
          xlabel="Number of orbits (target satellite)", ylabel="ΔP, kg⋅m/s",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          top_margin=title_top_margin_mm*Plots.mm)
_draw_panel!(p1, x_axis, all_ΔP, mean_ΔP, std_ΔP; scale=1.0)
if orbit_common[end] >= 1.0
    _v1  = interp1(orbit_common, mean_ΔP, 1.0)
    _yl1 = ylims(p1);  _xl1 = xlims(p1)
    plot!(p1, [1.0, 1.0], [_yl1[1], _v1], color=:gray, linestyle=:dash, lw=0.8, label="")
    scatter!(p1, [1.0], [_v1], markersize=5, color=:black, markerstrokewidth=0.5,
             markerstrokecolor=:black, label="")
    _ax1 = 1.0 + annot_dx_p1 * (_xl1[2] - _xl1[1])
    _ay1 = _v1 + annot_dy_p1 * (_yl1[2] - _yl1[1])
    annotate!(p1, _ax1, _ay1,
              Plots.text("(1, $(_fmt(_v1)))", Plots.default(:tickfontsize), :left, :bottom, :black))
    ylims!(p1, _yl1);  xlims!(p1, _xl1)
end

# ── Panel 2: Angular Momentum ─────────────────────────────────────────────────
p2 = plot(legend=false, grid=true, tick_direction=:out,
          xlabel="Number of orbits (target satellite)", ylabel="ΔH, kg⋅m²/s",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          top_margin=panel_top_margin_mm*Plots.mm)
_draw_panel!(p2, x_axis, all_ΔH, mean_ΔH, std_ΔH; scale=1.0)
if orbit_common[end] >= 1.0
    _v1  = interp1(orbit_common, mean_ΔH, 1.0)
    _yl2 = ylims(p2);  _xl2 = xlims(p2)
    plot!(p2, [1.0, 1.0], [_yl2[1], _v1], color=:gray, linestyle=:dash, lw=0.8, label="")
    scatter!(p2, [1.0], [_v1], markersize=5, color=:black, markerstrokewidth=0.5,
             markerstrokecolor=:black, label="")
    _ax2 = 1.0 + annot_dx_p2 * (_xl2[2] - _xl2[1])
    _ay2 = _v1 + annot_dy_p2 * (_yl2[2] - _yl2[1])
    annotate!(p2, _ax2, _ay2,
              Plots.text("(1, $(_fmt(_v1)))", Plots.default(:tickfontsize), :left, :bottom, :black))
    ylims!(p2, _yl2);  xlims!(p2, _xl2)
end

# ── Panel 3: Orbital Energy ───────────────────────────────────────────────────
p3 = plot(legend=false, grid=true, tick_direction=:out,
          xlabel="Number of orbits (target satellite)", ylabel="Δε, J",
          xguidefontsize=xlabel_fontsize, yguidefontsize=ylabel_fontsize,
          top_margin=panel_top_margin_mm*Plots.mm)
_draw_panel!(p3, x_axis, all_ΔE, mean_ΔE, std_ΔE; scale=1.0)
if orbit_common[end] >= 1.0
    _v1  = interp1(orbit_common, mean_ΔE, 1.0)
    _yl3 = ylims(p3);  _xl3 = xlims(p3)
    plot!(p3, [1.0, 1.0], [_yl3[1], _v1], color=:gray, linestyle=:dash, lw=0.8, label="")
    scatter!(p3, [1.0], [_v1], markersize=5, color=:black, markerstrokewidth=0.5,
             markerstrokecolor=:black, label="")
    _ax3 = 1.0 + annot_dx_p3 * (_xl3[2] - _xl3[1])
    _ay3 = _v1 + annot_dy_p3 * (_yl3[2] - _yl3[1])
    annotate!(p3, _ax3, _ay3,
              Plots.text("(1, $(_fmt(_v1)))", Plots.default(:tickfontsize), :left, :bottom, :black))
    ylims!(p3, _yl3);  xlims!(p3, _xl3)
end

# ── Combine and save ──────────────────────────────────────────────────────────
fig = plot(p1, p2, p3, layout=grid(3, 1, heights=panel_heights[1:3] ./ sum(panel_heights[1:3])), size=fig_size)

output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "conservation"))
mkpath(output_dir)
_N_label  = N_csv !== nothing ? @sprintf("_N%d", N_csv) : ""
_J2_label = filter_J2 === true ? "_J2T" : filter_J2 === false ? "_J2F" : ""
out_path     = joinpath(output_dir, @sprintf("conservation_h%dkm_i%.1fdeg%s%s.png", h, i_deg, _N_label, _J2_label))
out_path_pdf = joinpath(output_dir, @sprintf("conservation_h%dkm_i%.1fdeg%s%s.pdf", h, i_deg, _N_label, _J2_label))
savefig(fig, out_path)
savefig(fig, out_path_pdf)
display(fig)
println("Saved to: ", out_path)
println("Saved to: ", out_path_pdf)
