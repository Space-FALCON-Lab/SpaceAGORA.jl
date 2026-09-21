# 5x4 combined OE grid: rows = altitudes, columns = [Δa, Δe, Δi, ΔΩ].
# Each subplot shows 3 curves (one per target inclination: 0°, 0.5°, 1°).
# Horizontal legend strip across the bottom of the figure.

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
        hi_t = ceil(hi  / step) * step
        lo_t = hi_t - (n - 1) * step
        lo_t <= lo || continue
        vals = [lo_t + k * step for k in 0:(n-1)]
        ndec = max(0, -floor(Int, log10(step)))
        labs = step >= 1 ? [latexstring(string(round(Int, v))) for v in vals] :
                           [latexstring(@sprintf("%.*f", ndec, v)) for v in vals]
        return vals, labs
    end
    vals = collect(range(lo, hi, length=n))
    return vals, [latexstring(@sprintf("%.2g", v)) for v in vals]
end

# ── Knobs ─────────────────────────────────────────────────────────────────────
gap_inner           = 0
margin_left         = 0
margin_bottom       = 0
margin_right        = 2
margin_top          = 2
tick_fontsize       = 12
main_label_fontsize = 12
col_title_fontsize  = 12
row_label_fontsize  = 12
legend_fontsize     = 12
marker_stroke_width = 0.1   # black outline width on circle markers (0 = no outline)
main_xlabel         = L"N_{\mathrm{helpers}}"
h_helper_km         = 1000    # helper orbit altitude [km] — used for saturation number
L_max_m             = 200e3   # maximum laser range [m]    — used for saturation number

target_altitudes_km     = [850, 950, 1000, 1050, 1150]
target_inclinations_deg = [0.0, 0.5, 1.0]
elements                = [:a, :e, :i, :Ω]
col_titles              = [L"\Delta a,\ \mathrm{m}", L"\Delta e",
                           L"\Delta i,\ \mathrm{deg}", L"\Delta\Omega,\ \mathrm{deg}"]

inc_colors = [:blue, :red, :green]
inc_labels = [latexstring("i=$(round(i, digits=1))^{\\circ}")
              for i in target_inclinations_deg]

n_rows = length(target_altitudes_km)
n_cols = length(elements)

fig_width  = 220*n_cols + 80   # ~960
fig_height = 180*n_rows + 50   # ~950

# ── Build subplots (row-major: for each altitude, loop over elements) ─────────
subplots = []

for (row_idx, h) in enumerate(target_altitudes_km)
    filter_helpers = h == 1000 ? sort(unique([1, 50, 100, 150, 200, 249, 300])) :
                                 sort(unique([1, 50, 100, 150, 200, 250, 300]))

    is_top    = row_idx == 1
    is_bottom = row_idx == n_rows

    for (col_idx, elem) in enumerate(elements)
        is_left = col_idx == 1

        # ── Collect data for all three inclinations ──────────────────────────
        all_xs   = Vector{Vector{Float64}}()
        all_yf   = Vector{Vector{Float64}}()
        all_Nsat = Vector{Union{Float64,Nothing}}()

        for i_deg in target_inclinations_deg
            local res = try
                plot_OE(h, i_deg; element=elem, show_variance=false,
                        filter_helpers=filter_helpers)
            catch
                nothing
            end
            if res === nothing
                push!(all_xs,   Float64[])
                push!(all_yf,   Float64[])
                push!(all_Nsat, nothing)
            else
                N_vec = res.N_vals .- 1
                if h == 1000
                    N_vec = [n == 249 ? 250 : n for n in N_vec]
                end
                vals = [v[end] - v[1] for v in res.all_vals]
                # Snap near-zero noise to 0 for i and Ω
                if elem in (:i, :Ω) && maximum(abs.(vals)) < 1e-6
                    fill!(vals, 0.0)
                end
                ord = sortperm(N_vec)
                push!(all_xs, Float64.(N_vec[ord]))
                push!(all_yf,  vals[ord])

                local N_sat_val = try
                    N_tmp, _ = saturation_number(h_helper_km * 1e3, 0.0, h * 1e3, i_deg, L_max_m)
                    N_tmp
                catch
                    nothing
                end
                push!(all_Nsat,
                      (N_sat_val !== nothing &&
                       N_sat_val <= N_vec[ord][end]) ? Float64(N_sat_val) : nothing)
            end
        end

        all_vals_flat = vcat(all_yf...)

        row_ylab  = is_left  ? latexstring("h=$(h)\\,\\mathrm{km}") : ""
        col_title = is_top   ? col_titles[col_idx] : ""
        x_lab     = is_bottom ? main_xlabel : ""

        # ── Empty data guard ─────────────────────────────────────────────────
        if isempty(all_vals_flat)
            sp = plot(legend=false, grid=false, framestyle=:box,
                      title=col_title,  titlefontsize=col_title_fontsize,
                      ylabel=row_ylab,  yguidefontsize=row_label_fontsize,
                      xlabel=x_lab,     xguidefontsize=main_label_fontsize,
                      left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                      top_margin=gap_inner*Plots.mm,  bottom_margin=gap_inner*Plots.mm)
            annotate!(sp, 0.5, 0.5, text("No Data", :grey, :center, 8))
            push!(subplots, sp)
            continue
        end

        # ── y-axis range ─────────────────────────────────────────────────────
        _vmin  = minimum(all_vals_flat); _vmax = maximum(all_vals_flat)
        _vctr  = (_vmin + _vmax) / 2
        _vhalf = max((_vmax - _vmin) / 2 * 1.1, abs(_vctr) * 0.5, 1e-9)

        ytk_vals, ytk_labs =
            if elem == :i && maximum(abs.(all_vals_flat)) < 1e-6
                nice_ticks(-0.0020, 0.0020; n=5)
            elseif elem == :Ω && maximum(abs.(all_vals_flat)) < 1e-6
                nice_ticks(-20.0, 20.0; n=5)
            else
                nice_ticks(_vctr - _vhalf, _vctr + _vhalf; n=5)
            end

        ytk_step = ytk_vals[2] - ytk_vals[1]
        ylim_top = elem in (:i, :Ω) ? ytk_vals[end] + ytk_step : ytk_vals[end]

        ref_xs = first(filter(!isempty, all_xs))

        # ── Create subplot ───────────────────────────────────────────────────
        sp = plot(legend=false, grid=true,
                  title=col_title,  titlefontsize=col_title_fontsize,
                  ylabel=row_ylab,  yguidefontsize=row_label_fontsize,
                  xlabel=x_lab,     xguidefontsize=main_label_fontsize,
                  xtickfontsize=tick_fontsize, ytickfontsize=tick_fontsize,
                  tick_direction=:out,
                  xticks=(ref_xs, latexstring.(string.(round.(Int, ref_xs)))),
                  xrotation=45,
                  ylims=(ytk_vals[1], ylim_top),
                  yticks=(ytk_vals, ytk_labs),
                  left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                  top_margin=gap_inner*Plots.mm,  bottom_margin=gap_inner*Plots.mm)

        for (k, (xs, yf, N_sat_val)) in enumerate(zip(all_xs, all_yf, all_Nsat))
            isempty(xs) && continue
            plot!(sp, xs, yf;
                  color=inc_colors[k], linestyle=:solid, lw=1.5,
                  marker=:circle, ms=3, markerstrokewidth=marker_stroke_width, label=false)
            if N_sat_val !== nothing
                vline!(sp, [N_sat_val]; color=inc_colors[k], linestyle=:dash, lw=1.2, label=false)
            end
        end

        push!(subplots, sp)
    end
end

# ── Horizontal legend strip across the bottom ─────────────────────────────────
legend_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false,
                 xlims=(0,1), ylims=(0,1))

# Three columns (one per inclination colour); solid row on top, dashed row below.
legend_xs = [0.15, 0.45, 0.75]   # [i=0°, i=0.5°, i=1°]

# Row 1 – solid lines with circle markers
for (xc, col, lab) in zip(legend_xs, inc_colors, inc_labels)
    plot!(legend_sp, [xc, xc + 0.06], [0.72, 0.72]; color=col, lw=1.5, linestyle=:solid)
    scatter!(legend_sp, [xc + 0.03], [0.72]; color=col, marker=:circle, ms=3,
             markerstrokewidth=marker_stroke_width, label=false)
    annotate!(legend_sp, xc + 0.07, 0.72, text(lab, :black, :left, legend_fontsize))
end

# Row 2 – dashed lines (N_sat), same colours, aligned under row 1
for (xc, col) in zip(legend_xs, inc_colors)
    plot!(legend_sp, [xc, xc + 0.06], [0.28, 0.28]; color=col, lw=1.2, linestyle=:dash)
    annotate!(legend_sp, xc + 0.07, 0.28, text(L"N_{\mathrm{sat}}", :black, :left, legend_fontsize))
end

# ── Assemble and save ─────────────────────────────────────────────────────────
l = @layout [grid(n_rows, n_cols)
             b{0.07h}]

fig = plot(subplots..., legend_sp,
           layout=l,
           size=(fig_width, fig_height),
           left_margin=margin_left*Plots.mm,   right_margin=margin_right*Plots.mm,
           top_margin=margin_top*Plots.mm,      bottom_margin=margin_bottom*Plots.mm)

output_dir = normpath(joinpath(@__DIR__, "..", "..", "output",
                                "images_reconstructed_from_csv", "orbital_elements"))
mkpath(output_dir)
out_path     = joinpath(output_dir, "OE_all_final_NxN_plot.png")
out_path_pdf = joinpath(output_dir, "OE_all_final_NxN_plot.pdf")
savefig(fig, out_path)
savefig(fig, out_path_pdf)
display(fig)
println("Saved to: ", out_path)
println("Saved to: ", out_path_pdf)
