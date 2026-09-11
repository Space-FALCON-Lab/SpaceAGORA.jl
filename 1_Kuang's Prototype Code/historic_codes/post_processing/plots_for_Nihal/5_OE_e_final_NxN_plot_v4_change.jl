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

inc_colors  = [:blue, :red, :green]
inc_labels  = [latexstring("i=$(round(i, digits=1))^{\\circ}") for i in target_inclinations_deg]

n_rows = length(target_altitudes_km)

fig_width  = 500
fig_height = 200*n_rows + 60

subplots = []

for (row_idx, h) in enumerate(target_altitudes_km)
    local sp
    filter_helpers = h == 1000 ? sort(unique([1, 50, 100, 150, 200, 249, 300])) :
                                 sort(unique([1, 50, 100, 150, 200, 250, 300]))

    is_bottom = row_idx == n_rows
    row_ylab  = latexstring("h=$(h)\\,\\mathrm{km}")

    all_xs      = Vector{Vector{Float64}}()
    all_efinals = Vector{Vector{Float64}}()

    for i_deg in target_inclinations_deg
        local res = try
            plot_OE(h, i_deg; element=:e, show_variance=false,
                    filter_helpers=filter_helpers)
        catch
            nothing
        end
        if res === nothing
            push!(all_xs, Float64[])
            push!(all_efinals, Float64[])
        else
            N_helpers_vec = res.N_vals .- 1
            if h == 1000; N_helpers_vec = [n == 249 ? 250 : n for n in N_helpers_vec]; end
            e_final = [v[end] - v[1] for v in res.all_vals]
            ord = sortperm(N_helpers_vec)
            push!(all_xs,      Float64.(N_helpers_vec[ord]))
            push!(all_efinals, e_final[ord])
        end
    end

    all_vals_flat = vcat(all_efinals...)
    if isempty(all_vals_flat)
        sp = plot(legend=false, grid=false, framestyle=:box,
                  ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                  left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                  top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)
        annotate!(sp, 0.5, 0.5, text("No Data", :grey, :center, 8))
        push!(subplots, sp)
        continue
    end

    _emin = minimum(all_vals_flat); _emax = maximum(all_vals_flat)
    _ectr  = (_emin + _emax) / 2
    _ehalf = max((_emax - _emin) / 2 * 1.1, abs(_ectr) * 0.5, 1e-9)
    ytk_vals, ytk_labs = nice_ticks(_ectr - _ehalf, _ectr + _ehalf; n=5)

    ref_xs = first(filter(!isempty, all_xs))

    sp = plot(legend= row_idx == 1 ? :topright : false,
              grid=true,
              ylabel=row_ylab,    yguidefontsize=row_label_fontsize,
              xlabel= is_bottom ? main_xlabel : "",
              xguidefontsize=main_label_fontsize,
              xtickfontsize=tick_fontsize, ytickfontsize=tick_fontsize,
              tick_direction=:out,
              xticks=(ref_xs, latexstring.(string.(round.(Int, ref_xs)))),
              xrotation=45,
              ylims=(ytk_vals[1], ytk_vals[end]),
              yticks=(ytk_vals, ytk_labs),
              legendfontsize=legend_fontsize,
              left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
              top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)

    for (k, (xs, yf)) in enumerate(zip(all_xs, all_efinals))
        isempty(xs) && continue
        plot!(sp, xs, yf;
              color=inc_colors[k], linestyle=:solid, lw=1.5,
              marker=:circle, ms=3, label=inc_labels[k])
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
