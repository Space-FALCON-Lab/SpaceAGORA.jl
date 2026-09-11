# NxN grid of final eccentricity vs number of helper satellites.
# Rows = altitudes, Columns = inclinations.
# Each subplot shows e(end) for every helper count in filter_helpers.

include(joinpath(@__DIR__, "..", "csv_plotter_functions", "9_OE_function.jl"))
using LaTeXStrings
using Printf

# ── Exactly n ticks: ceil(hi) as last tick, extend down n-1 steps ────────────
function nice_ticks(lo, hi; n=5)
    span = hi - lo
    span == 0 && (span = abs(lo) > 0 ? abs(lo) * 0.2 : 1.0)
    raw  = span / (n - 1)
    mag  = 10.0^floor(log10(raw))
    for mult in (1.0, 2.0, 5.0, 10.0, 20.0, 50.0)
        step = mult * mag
        hi_t = ceil(hi  / step) * step   # round UP from max data (closest upper limit)
        lo_t = hi_t - (n - 1) * step     # extend down exactly n-1 steps
        lo_t <= lo || continue            # must cover min data too
        vals = [lo_t + k * step for k in 0:(n-1)]
        ndec = max(0, -floor(Int, log10(step)))
        labs = step >= 1 ? [latexstring(string(round(Int, v))) for v in vals] :
                           [latexstring(@sprintf("%.*f", ndec, v)) for v in vals]
        return vals, labs
    end
    vals = collect(range(lo, hi, length=n))
    return vals, [latexstring(@sprintf("%.2g", v)) for v in vals]
end

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
main_ylabel          = L"\Delta e"
row_label_fontsize   = 10
col_label_fontsize   = 10
ylabel_gap           = 0.5
xlabel_gap           = 0.5
# filter_helpers defined per-altitude inside the loop (h=1000 uses 249 in place of 250)
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
        filter_helpers = h == 1000 ? sort(unique([1, 50, 100, 150, 200, 249, 300])) :
                                     sort(unique([1, 50, 100, 150, 200, 250, 300]))
        local res = try
            plot_OE(h, i_deg; element=:e, show_variance=false,
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
            # Final eccentricity for each helper count — sort ascending
            N_helpers_vec = res.N_vals .- 1
            if h == 1000; N_helpers_vec = [n == 249 ? 250 : n for n in N_helpers_vec]; end
            e_final       = [v[end] - v[1] for v in res.all_vals]   # Δe = final − initial
            ord           = sortperm(N_helpers_vec)
            N_helpers_vec = N_helpers_vec[ord]
            e_final       = e_final[ord]
            xs            = N_helpers_vec             # actual values → proper spacing

            _emin = minimum(e_final); _emax = maximum(e_final)
            _ectr  = (_emin + _emax) / 2
            _ehalf = max((_emax - _emin) / 2 * 1.1, abs(_ectr) * 0.5, 1e-9)
            ytk_vals, ytk_labs = nice_ticks(_ectr - _ehalf, _ectr + _ehalf; n=5)

            sp = plot(legend=false, grid=true,
                      ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                      xlabel=col_xlab, xguidefontsize=col_label_fontsize,
                      xtickfontsize=tick_fontsize, ytickfontsize=tick_fontsize,
                      tick_direction=:out,
                      xticks=(xs, latexstring.(string.(N_helpers_vec))),
                      xrotation=45,
                      ylims=(ytk_vals[1], ytk_vals[end]),
                      yticks=(ytk_vals, ytk_labs),
                      left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                      top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)
                plot!(sp, xs, e_final;
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
          text(L"\Delta e", :black, :left, legend_fontsize))
corner_leg_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false)

l = @layout [a{0.06w} grid(n_rows, n_cols) ; b{0.001w} c{0.05h}]

fig = plot(yl_sp, subplots..., corner_sp, xl_sp,
           layout=l,
           size=(fig_width, fig_height),
           left_margin=margin_left*Plots.mm,   right_margin=margin_right*Plots.mm,
           top_margin=margin_top*Plots.mm,      bottom_margin=margin_bottom*Plots.mm)

output_dir = normpath(joinpath(@__DIR__, "..", "..", "output",
                                "images_reconstructed_from_csv", "orbital_elements"))
mkpath(output_dir)
out_path     = joinpath(output_dir, "OE_e_final_NxN_plot.png")
out_path_pdf = joinpath(output_dir, "OE_e_final_NxN_plot.pdf")
savefig(fig, out_path)
savefig(fig, out_path_pdf)
display(fig)
println("Saved to: ", out_path)
println("Saved to: ", out_path_pdf)
