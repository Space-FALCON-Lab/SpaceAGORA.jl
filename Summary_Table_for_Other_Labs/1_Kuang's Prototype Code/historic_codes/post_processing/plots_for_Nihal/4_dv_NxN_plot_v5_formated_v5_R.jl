# NxN grid of final Δv vs Number of Helpers subplots.
# Rows = altitudes, Columns = inclinations.
# Each subplot shows Δv_R, Δv_T, Δv_N at end of simulation vs number of helper satellites.
# Calls plot_dv(h, i) for data; blank with "No Data" if unavailable.

include(joinpath(@__DIR__, "..", "csv_plotter_functions", "7_dv_function.jl"))
using LaTeXStrings
using Printf

# ── Knobs ──────────────────────────────────────────────────────────────────────
# Negative values pull subplots together; positive push them apart (in mm).
gap_inner   = 0   # inner subplot margin on all four sides (drives inter-subplot gap)
# Outer figure margins — increase if axis labels/ticks get clipped
margin_left   = 5
margin_bottom = -2
margin_right  =  0
margin_top    =  0
tick_fontsize        = 10   # x and y tick number size
subtitle_fontsize    = 10   # per-subplot title font size
maintitle_fontsize   = 10  # main figure title font size
main_label_fontsize  = 10   # shared x/y axis label font size
legend_fontsize      = 8  # legend font size — tune this knob
main_xlabel          = L"N_{\mathrm{helpers}}"                          # shared x-axis label
main_ylabel          = L"\Delta v,\ \mathrm{m/s}"  # shared y-axis label
row_label_fontsize   = 10   # font size for h=Xkm row labels on first-column subplots
col_label_fontsize   = 10   # font size for i=X.X° col labels on bottom-row subplots
ylabel_gap           = 0.5 # x-position of y-label within its strip (0=far left, 1=touching grid)
xlabel_gap           = 0.5 # y-position of x-label within its strip (0=far from grid, 1=touching grid)
# filter_helpers defined per-altitude inside the loop (h=1000 uses 249 in place of 250)
h_helper_km          = 1000    # helper orbit altitude [km] — used for saturation number
L_max_m              = 200e3   # maximum laser range [m]    — used for saturation number
fig_width            = 160*n_cols + 120   # figure width  [px] — override freely
fig_height           = 160*n_rows + 30   # figure height [px] — override freely
# ───────────────────────────────────────────────────────────────────────────────

target_altitudes_km     = [850, 950, 1000, 1050, 1150]
target_inclinations_deg = [0.0, 0.5, 1.0]

n_rows = length(target_altitudes_km)

inc_colors = [:blue, :red, :green]
inc_labels = [latexstring("i=$(round(i, digits=1))^{\\circ}")
              for i in target_inclinations_deg]

fig_width  = 500
fig_height = 200*n_rows + 80

subplots = []

for (row_idx, h) in enumerate(target_altitudes_km)
    local sp
    filter_helpers = h == 1000 ? sort(unique([1, 50, 150, 200, 249, 300])) :
                                 sort(unique([1, 50, 150, 200, 250, 300]))

    is_bottom = row_idx == n_rows
    row_ylab  = latexstring("h=$(h)\\,\\mathrm{km}")

    all_xs   = Vector{Vector{Float64}}()
    all_dvR  = Vector{Vector{Float64}}()
    all_Nsat = Vector{Union{Float64,Nothing}}()

    for i_deg in target_inclinations_deg
        local res = try
            plot_dv(h, i_deg)
        catch
            nothing
        end

        if res !== nothing && !isempty(filter_helpers)
            mask = [n in filter_helpers for n in res.N_helpers_vec]
            res  = (N_helpers_vec = res.N_helpers_vec[mask],
                    dv_timeseries = res.dv_timeseries[mask])
        end
        if h == 1000 && res !== nothing
            res = (N_helpers_vec = [n == 249 ? 250 : n for n in res.N_helpers_vec],
                   dv_timeseries = res.dv_timeseries)
        end

        if res === nothing
            push!(all_xs,   Float64[])
            push!(all_dvR,  Float64[])
            push!(all_Nsat, nothing)
        else
            N_vec = Float64.(res.N_helpers_vec)
            dvR   = [ts.R[end] for ts in res.dv_timeseries]
            push!(all_xs,  N_vec)
            push!(all_dvR, dvR)

            local N_sat_val = try
                N_tmp, _ = saturation_number(h_helper_km * 1e3, 0.0, h * 1e3, i_deg, L_max_m)
                N_tmp
            catch
                nothing
            end
            push!(all_Nsat,
                  (N_sat_val !== nothing &&
                   N_sat_val <= res.N_helpers_vec[end]) ? Float64(N_sat_val) : nothing)
        end
    end

    all_vals_flat = vcat(all_dvR...)
    if isempty(all_vals_flat)
        sp = plot(legend=false, grid=false, framestyle=:box,
                  ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                  left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                  top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)
        annotate!(sp, 0.5, 0.5, text("No Data", :grey, :center, 8))
        push!(subplots, sp)
        continue
    end

    dv_max = maximum(all_vals_flat); dv_min = minimum(all_vals_flat)
    half_span = max((dv_max - dv_min) / 2, abs(dv_max) * 0.5, 1e-9)
    y_lo = dv_min - half_span
    y_hi = dv_max + half_span
    ytick_vals   = collect(range(y_lo, y_hi, length=5))
    ytick_labels = [latexstring(@sprintf("%.2f", v)) for v in ytick_vals]

    ref_xs = first(filter(!isempty, all_xs))

    sp = plot(legend= row_idx == 1 ? :topright : false,
              grid=true,
              ylabel=row_ylab,  yguidefontsize=row_label_fontsize,
              xlabel= is_bottom ? main_xlabel : "",
              xguidefontsize=main_label_fontsize,
              xtickfontsize=tick_fontsize, ytickfontsize=tick_fontsize,
              tick_direction=:out,
              xticks=(ref_xs, latexstring.(string.(round.(Int, ref_xs)))),
              xrotation=45,
              ylims=(y_lo, y_hi),
              yticks=(ytick_vals, ytick_labels),
              legendfontsize=legend_fontsize,
              left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
              top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)

    for (k, (xs, dvR, N_sat_val)) in enumerate(zip(all_xs, all_dvR, all_Nsat))
        isempty(xs) && continue
        plot!(sp, xs, dvR;
              color=inc_colors[k], linestyle=:solid, lw=1.5,
              marker=:circle, ms=3, label=inc_labels[k])
        if N_sat_val !== nothing
            vline!(sp, [N_sat_val]; color=inc_colors[k], linestyle=:dash, lw=1.2, label=false)
        end
    end

    push!(subplots, sp)
end

# ── Shared axis label strips ──────────────────────────────────────────────────
yl_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false,
             xlims=(0,1), ylims=(0,1))
annotate!(yl_sp, ylabel_gap, 0.5,
          text(main_ylabel, :black, :center, main_label_fontsize, rotation=90))

corner_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false)

l = @layout [a{0.06w} grid(n_rows, 1)]

fig = plot(yl_sp, subplots...,
           layout=l,
           size=(fig_width, fig_height),
           left_margin=margin_left*Plots.mm,   right_margin=margin_right*Plots.mm,
           top_margin=margin_top*Plots.mm,      bottom_margin=margin_bottom*Plots.mm)

output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "delta_v"))
mkpath(output_dir)
out_path     = joinpath(output_dir, "dv_final_NxN_plot.png")
out_path_pdf = joinpath(output_dir, "dv_final_NxN_plot.pdf")
savefig(fig, out_path)
savefig(fig, out_path_pdf)
display(fig)
println("Saved to: ", out_path)
println("Saved to: ", out_path_pdf)

