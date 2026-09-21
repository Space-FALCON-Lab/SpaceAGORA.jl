# NxN grid of final semimajor axis vs number of helper satellites.
# Rows = altitudes, Columns = inclinations.
# Each subplot shows a(end) for every helper count in filter_helpers.

include(joinpath(@__DIR__, "..", "csv_plotter_functions", "9_OE_function.jl"))
using LaTeXStrings
using Printf

# ── Knobs ──────────────────────────────────────────────────────────────────────
gap_inner            = 0
margin_left          = 5
margin_bottom        = 0
margin_right         = 0
margin_top           = 0
tick_fontsize        = 10
main_label_fontsize  = 10
legend_fontsize      = 8
main_xlabel          = L"N_{\mathrm{helpers}}"
main_ylabel          = L"\Delta a,\ \mathrm{m}"
row_label_fontsize   = 10
col_label_fontsize   = 10
ylabel_gap           = 0.5
xlabel_gap           = 0.5
filter_helpers       = sort(unique([1, 50, 100, 150, 200, 250, 300]))
# ───────────────────────────────────────────────────────────────────────────────

target_altitudes_km     = [850, 950, 1000, 1050, 1150]
target_inclinations_deg = [0.0, 0.5, 1.0]

n_cols = length(target_inclinations_deg)
n_rows = length(target_altitudes_km)

fig_width  = 160*n_cols + 120
fig_height = 160*n_rows + 30

subplots = []

for (row_idx, h) in enumerate(target_altitudes_km)
    for (col_idx, i_deg) in enumerate(target_inclinations_deg)
        local sp
        local res = try
            plot_OE(h, i_deg; element=:a, show_variance=false,
                    filter_helpers=filter_helpers)
        catch
            nothing
        end

        is_left   = col_idx == 1
        is_bottom = row_idx == n_rows
        row_ylab  = is_left   ? latexstring("h=$(h)\\,\\mathrm{km}")               : ""
        col_xlab  = is_bottom ? latexstring("i=$(round(i_deg, digits=1))^{\\circ}") : ""

        if res === nothing
            sp = plot(legend=false, grid=false, framestyle=:box,
                      ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                      xlabel=col_xlab, xguidefontsize=col_label_fontsize,
                      left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                      top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)
            annotate!(sp, 0.5, 0.5, text("No Data", :grey, :center, 8))
        else
            # Final semimajor axis for each helper count — sort ascending
            N_helpers_vec = res.N_vals .- 1
            a_final       = [v[end] - v[1] for v in res.all_vals]   # Δa = final − initial
            ord           = sortperm(N_helpers_vec)
            N_helpers_vec = N_helpers_vec[ord]
            a_final       = a_final[ord]
            xs            = N_helpers_vec             # actual values → proper spacing

            pad       = max((maximum(a_final) - minimum(a_final)) * 0.1, 1.0)
            y_lo      = minimum(a_final) - pad
            y_hi      = maximum(a_final) + pad
            ytick_vals   = collect(range(y_lo, y_hi, length=5))
            step_est     = (y_hi - y_lo) / 4
            ndec         = step_est >= 1 ? 0 : max(0, -floor(Int, log10(abs(step_est))))
            ytick_labels = step_est >= 1 ? [latexstring(string(round(Int, v))) for v in ytick_vals] :
                                           [latexstring(@sprintf("%.*f", ndec, v)) for v in ytick_vals]

            sp = plot(legend=false, grid=true,
                      ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                      xlabel=col_xlab, xguidefontsize=col_label_fontsize,
                      xtickfontsize=tick_fontsize, ytickfontsize=tick_fontsize,
                      tick_direction=:out,
                      xticks=(xs, latexstring.(string.(N_helpers_vec))),
                      xrotation=45,
                      ylims=(y_lo, y_hi),
                      yticks=(ytick_vals, ytick_labels),
                      left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                      top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)
            plot!(sp, xs, a_final;
                  color=:blue, linestyle=:solid, lw=1.5, marker=:circle, ms=3, label=false)
        end

        push!(subplots, sp)
    end
end

# ── Shared axis label strips ──────────────────────────────────────────────────
yl_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false,
             xlims=(0,1), ylims=(0,1))
annotate!(yl_sp, ylabel_gap, 0.5,
          text(main_ylabel, :black, :center, main_label_fontsize, rotation=90))

xl_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false,
             xlims=(0,1), ylims=(0,1))
annotate!(xl_sp, 0.5, xlabel_gap, text(main_xlabel, :black, :center, main_label_fontsize))

corner_sp     = plot(framestyle=:none, ticks=nothing, legend=false, grid=false)

# Legend strip
legend_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false,
                 xlims=(0,1), ylims=(0,1))
plot!(legend_sp, [0.35, 0.43], [0.5, 0.5]; color=:blue, linestyle=:solid, lw=1.5)
scatter!(legend_sp, [0.39], [0.5]; color=:blue, marker=:circle, ms=3, label=false)
annotate!(legend_sp, 0.45, 0.5,
          text(L"\Delta a", :black, :left, legend_fontsize))
corner_leg_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false)

l = @layout [a{0.06w} grid(n_rows, n_cols) ; b{0.001w} c{0.001h} ; d{0.001w} e{0.05h}]

fig = plot(yl_sp, subplots..., corner_sp, xl_sp, corner_leg_sp, legend_sp,
           layout=l,
           size=(fig_width, fig_height),
           left_margin=margin_left*Plots.mm,   right_margin=margin_right*Plots.mm,
           top_margin=margin_top*Plots.mm,      bottom_margin=margin_bottom*Plots.mm)

output_dir = normpath(joinpath(@__DIR__, "..", "..", "output",
                                "images_reconstructed_from_csv", "orbital_elements"))
mkpath(output_dir)
out_path     = joinpath(output_dir, "OE_a_final_NxN_plot.png")
out_path_pdf = joinpath(output_dir, "OE_a_final_NxN_plot.pdf")
savefig(fig, out_path)
savefig(fig, out_path_pdf)
display(fig)
println("Saved to: ", out_path)
println("Saved to: ", out_path_pdf)
