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
margin_bottom = 0
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
filter_helpers       = sort(unique([1, 50, 200, 150, 200, 250, 300]))  # helper counts to include (N_csv = helper+1); empty = all
h_helper_km          = 1000    # helper orbit altitude [km] — used for saturation number
L_max_m              = 200e3   # maximum laser range [m]    — used for saturation number
fig_width            = 160*n_cols + 120   # figure width  [px] — override freely
fig_height           = 160*n_rows + 30   # figure height [px] — override freely
# ───────────────────────────────────────────────────────────────────────────────

target_altitudes_km     = [850, 950, 1000, 1050, 1150]
target_inclinations_deg = [0.0, 0.5, 1.0]

n_cols = length(target_inclinations_deg)
n_rows = length(target_altitudes_km)

# Build one Plots.jl layout with n_rows × n_cols subplots
subplots = []

for (row_idx, h) in enumerate(target_altitudes_km)
    for (col_idx, i_deg) in enumerate(target_inclinations_deg)
        local sp
        local res = try
            plot_dv(h, i_deg)
        catch
            nothing
        end

        # Filter to the requested helper counts (N in CSV = helper + 1)
        if res !== nothing && !isempty(filter_helpers)
            mask = [n in filter_helpers for n in res.N_helpers_vec]
            res  = (N_helpers_vec = res.N_helpers_vec[mask],
                    dv_timeseries = res.dv_timeseries[mask])
        end

        is_left   = col_idx == 1
        is_bottom = row_idx == n_rows
        row_ylab = is_left   ? latexstring("h=$(h)\\,\\mathrm{km}")                    : ""
        col_xlab = is_bottom ? latexstring("i=$(round(i_deg, digits=1))^{\\circ}") : ""

        if res === nothing
            # blank subplot with "No Data" annotation
            sp = plot(legend=false, grid=false, framestyle=:box,
                      ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                      xlabel=col_xlab, xguidefontsize=col_label_fontsize,
                      left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                      top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)
            annotate!(sp, 0.5, 0.5, text("No Data", :grey, :center, 8))
        else
            dv_final_R = [ts.R[end] for ts in res.dv_timeseries]
            xs = 1:length(res.N_helpers_vec)  # evenly spaced indices
            let pad = max((maximum(dv_final_R) - minimum(dv_final_R)) * 0.1, 1e-9)
                global y_lo = min(0.0, minimum(dv_final_R) - pad)
                global y_hi = maximum(dv_final_R) + pad
            end
            ytick_vals   = collect(range(y_lo, y_hi, length=5))
            ytick_labels = [latexstring(@sprintf("%.2f", v)) for v in ytick_vals]
            sp = plot(legend=false,
                      grid=true,
                      ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                      xlabel=col_xlab, xguidefontsize=col_label_fontsize,
                      xtickfontsize=tick_fontsize, ytickfontsize=tick_fontsize,
                      tick_direction=:out,
                      xticks=(xs, latexstring.(string.(res.N_helpers_vec))),
                      xrotation=45,
                      ylims=(y_lo, y_hi),
                      yticks=(ytick_vals, ytick_labels),
                      left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                      top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)
            plot!(sp, xs, dv_final_R;
                  label=L"\Delta v_R", color=:blue,   linestyle=:solid,   lw=1.5, marker=:circle, ms=3)

            # Saturation number — vertical red dashed line
            local N_sat_val = try
                N_tmp, _ = saturation_number(h_helper_km * 1e3, 0.0, h * 1e3, i_deg, L_max_m)
                N_tmp
            catch
                nothing
            end
            if h != 1000 && N_sat_val !== nothing && N_sat_val <= res.N_helpers_vec[end]
                nv = res.N_helpers_vec
                x_sat = if N_sat_val <= nv[1]
                    1.0
                else
                    k = findfirst(j -> nv[j+1] >= N_sat_val, 1:length(nv)-1)
                    k + (N_sat_val - nv[k]) / (nv[k+1] - nv[k])
                end
                vline!(sp, [x_sat]; color=:red, linestyle=:dash, lw=1.5, label=latexstring("N_{\\mathrm{sat}}=$(N_sat_val)"))
            end
        end

        push!(subplots, sp)
    end
end

# ── Shared axis label strips ─────────────────────────────────────────────────
# Left strip: rotated y-label; Bottom strip: x-label; Corner: blank filler
yl_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false,
             xlims=(0,1), ylims=(0,1))
annotate!(yl_sp, ylabel_gap, 0.5, text(main_ylabel, :black, :center, main_label_fontsize, rotation=90))

xl_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false,
             xlims=(0,1), ylims=(0,1))
annotate!(xl_sp, 0.5, xlabel_gap, text(main_xlabel, :black, :center, main_label_fontsize))

corner_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false)

# Legend strip — horizontal legend row between grid and x-label
legend_sp = plot(framestyle=:none, ticks=nothing, grid=false,
                 xlims=(0,1), ylims=(0,1), legend=:top,
                 legendfontsize=legend_fontsize, legend_columns=-1,
                 left_margin=20*Plots.mm, right_margin=20*Plots.mm)
plot!(legend_sp, [NaN], [NaN]; label=L"\Delta v_R",         color=:blue, linestyle=:solid, lw=1.5)
plot!(legend_sp, [NaN], [NaN]; label=L"N_{\mathrm{sat}}",   color=:red,  linestyle=:dash,  lw=1.5)
corner_leg_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false)

# 3-row macro layout: [ y-strip | NxN grid ]
#                     [ corner  | x-strip  ]
#                     [ corner  | legend   ]
l = @layout [a{0.06w} grid(n_rows, n_cols) ; b{0.001w} c{0.001h} ; d{0.001w} e{0.05h}] # change gap around main xlabel and main ylabel here !!!!!

fig = plot(yl_sp, subplots..., corner_sp, xl_sp, corner_leg_sp, legend_sp,
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

