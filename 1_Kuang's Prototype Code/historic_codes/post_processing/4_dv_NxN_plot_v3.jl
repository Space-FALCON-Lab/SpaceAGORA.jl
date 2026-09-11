# NxN grid of Δv change vs Number of Helpers subplots.
# Rows = altitudes, Columns = inclinations.
# Each subplot shows ΔΔv_R, ΔΔv_T, ΔΔv_N (end−start) vs number of helper satellites.
# Calls plot_dv(h, i) for data; blank with "No Data" if unavailable.

include(joinpath(@__DIR__, "csv_plotter_functions", "7_dv_function.jl"))

# ── Knobs ──────────────────────────────────────────────────────────────────────
# Negative values pull subplots together; positive push them apart (in mm).
gap_inner   = 0   # inner subplot margin on all four sides (drives inter-subplot gap)
# Outer figure margins — increase if axis labels/ticks get clipped
margin_left   = 0
margin_bottom = 0
margin_right  =  0.5
margin_top    =  0
tick_fontsize        = 8   # x and y tick number size
subtitle_fontsize    = 7   # per-subplot title font size
maintitle_fontsize   = 12  # main figure title font size
main_label_fontsize  = 9   # shared x/y axis label font size
legend_fontsize      = 8  # legend font size — tune this knob
main_xlabel          = "Number of helpers"   # shared x-axis label
 main_ylabel          = "ΔΔv  [m/s]"           # shared y-axis label
row_label_fontsize   = 8   # font size for h=Xkm row labels on first-column subplots
col_label_fontsize   = 8   # font size for i=X.X° col labels on bottom-row subplots
ylabel_gap           = 0.5 # x-position of y-label within its strip (0=far left, 1=touching grid)
xlabel_gap           = 0.5 # y-position of x-label within its strip (0=far from grid, 1=touching grid)
filter_helpers       = sort(unique([1, 50, 200, 150, 200, 250, 300]))  # helper counts to include (N_csv = helper+1); empty = all
# ───────────────────────────────────────────────────────────────────────────────

target_altitudes_km     = [850, 950, 1000, 1050, 1150]
target_inclinations_deg = [0.0, 0.5, 1.0, 1.5]

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
                    dv_change_R   = res.dv_change_R[mask],
                    dv_change_T   = res.dv_change_T[mask],
                    dv_change_N   = res.dv_change_N[mask])
        end

        is_left   = col_idx == 1
        is_bottom = row_idx == n_rows
        row_ylab = is_left   ? @sprintf("h=%dkm", h)          : ""
        col_xlab = is_bottom ? @sprintf("i=%.1f\u00b0", i_deg) : ""

        if res === nothing
            # blank subplot with "No Data" annotation
            sp = plot(legend=false, grid=false, framestyle=:box,
                      ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                      xlabel=col_xlab, xguidefontsize=col_label_fontsize,
                      left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                      top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)
            annotate!(sp, 0.5, 0.5, text("No Data", :grey, :center, 8))
        else
            xs = 1:length(res.N_helpers_vec)  # evenly spaced indices
            sp = plot(legend=false,
                      grid=true,
                      ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                      xlabel=col_xlab, xguidefontsize=col_label_fontsize,
                      xtickfontsize=tick_fontsize, ytickfontsize=tick_fontsize,
                      tick_direction=:out,
                      xticks=(xs, string.(res.N_helpers_vec)),
                      left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                      top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)
            plot!(sp, xs, res.dv_change_R;
                  label="ΔΔv_R", color=:blue,   linestyle=:solid,   lw=1.5, marker=:circle, ms=3)
            plot!(sp, xs, res.dv_change_T;
                  label="ΔΔv_T", color=:orange, linestyle=:dash,    lw=1.5, marker=:circle, ms=3)
            plot!(sp, xs, res.dv_change_N;
                  label="ΔΔv_N", color=:green,  linestyle=:dashdot, lw=1.5, marker=:circle, ms=3)
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
                 legendfontsize=legend_fontsize, legend_columns=-1)
plot!(legend_sp, [0,0], [2,2]; label="ΔΔv_R", color=:blue,   linestyle=:solid,   lw=1.5)
plot!(legend_sp, [0,0], [2,2]; label="ΔΔv_T", color=:orange, linestyle=:dash,    lw=1.5)
plot!(legend_sp, [0,0], [2,2]; label="ΔΔv_N", color=:green,  linestyle=:dashdot, lw=1.5)
corner_leg_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false)

# 3-row macro layout: [ y-strip | NxN grid ]
#                     [ corner  | x-strip  ]
#                     [ corner  | legend   ]
l = @layout [a{0.06w} grid(n_rows, n_cols) ; b{0.001w} c{0.001h} ; d{0.001w} e{0.05h}] # change gap around main xlabel and main ylabel here !!!!!

fig = plot(yl_sp, subplots..., corner_sp, xl_sp, corner_leg_sp, legend_sp,
           layout=l,
           size=(200*n_cols + 70, 160*n_rows + 50),
           plot_title="Change of Δv in RTN vs Number of Helper Satellites",
           plot_titlefontsize=maintitle_fontsize,
           left_margin=margin_left*Plots.mm,   right_margin=margin_right*Plots.mm,
           top_margin=margin_top*Plots.mm,      bottom_margin=margin_bottom*Plots.mm)

output_dir = normpath(joinpath(@__DIR__, "..", "output", "images_reconstructed_from_csv", "delta_v"))
mkpath(output_dir)
out_path = joinpath(output_dir, "dv_change_NxN_plot.png")
savefig(fig, out_path)
display(fig)
println("Saved to: ", out_path)

