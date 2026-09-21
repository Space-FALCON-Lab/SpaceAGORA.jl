# 1×N grid of Linear Momentum Fractional Drift subplots — varying ν (true anomaly offset).
# Fixed altitude and inclination; each column = one ν value.
# Each subplot calls plot_linear_momentum(h, i, nu=ν) for data; blank with "No Data" if unavailable.

include(joinpath(@__DIR__, "csv_plotter_functions", "4_linear_momentum_function.jl"))

# ── Knobs ──────────────────────────────────────────────────────────────────────
gap_inner   = 0      # inner subplot margin on all four sides (mm)
margin_left   = 0
margin_bottom = 0
margin_right  =  0.5
margin_top    =  0
tick_fontsize        = 8
subtitle_fontsize    = 7
maintitle_fontsize   = 12
main_label_fontsize  = 9
main_xlabel          = "t (hr)"
main_ylabel          = "ΔP/P₀  [×10⁻³]"
row_label_fontsize   = 8   # font size for h=Xkm label on the first-column subplot
col_label_fontsize   = 8   # font size for ν=X.XX° labels on bottom-row subplots
ylabel_gap           = 0.5 # x-position of y-label within its strip (0=far left, 1=touching grid)
xlabel_gap           = 0.5 # y-position of x-label within its strip (0=far from grid, 1=touching grid)
# ───────────────────────────────────────────────────────────────────────────────

target_h_km  = 1000
target_i_deg = 0.0
target_nu_deg = [-1.5, -0.75, 0.75, 1.5]   # true anomaly offsets (degrees)

n_rows = 1
n_cols = length(target_nu_deg)

# Build subplots — single row, one column per ν value
subplots = []

for (col_idx, nu) in enumerate(target_nu_deg)
    local sp
    local res = try
        plot_linear_momentum(target_h_km, target_i_deg; nu=nu, show_variance=true)
    catch
        nothing
    end

    is_left   = col_idx == 1
    is_bottom = true   # only one row, so always show bottom label
    row_ylab = is_left   ? @sprintf("h=%dkm", target_h_km) : ""
    col_xlab = @sprintf("\u03bd=%.2f\u00b0", nu)            # ν=X.XX° on every subplot

    if res === nothing
        sp = plot(legend=false, grid=false, framestyle=:box,
                  ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                  xlabel=col_xlab, xguidefontsize=col_label_fontsize,
                  left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                  top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)
        annotate!(sp, 0.5, 0.5, text("No Data", :grey, :center, 8))
    else
        t_hr = res.t ./ 3600.0
        sp = plot(legend=false, grid=true,
                  ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                  xlabel=col_xlab, xguidefontsize=col_label_fontsize,
                  xtickfontsize=tick_fontsize, ytickfontsize=tick_fontsize,
                  tick_direction=:out,
                  left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                  top_margin=gap_inner*Plots.mm, bottom_margin=gap_inner*Plots.mm)
        # individual grey curves
        for ΔP_curve in res.all_ΔP
            plot!(sp, t_hr, ΔP_curve .* 1e3, color=:grey, alpha=0.4, lw=0.5, label="")
        end
        # ±1σ band
        plot!(sp, t_hr, (res.mean_ΔP .+ res.std_ΔP) .* 1e3,
              fillrange=(res.mean_ΔP .- res.std_ΔP) .* 1e3,
              fillalpha=0.2, fillcolor=:blue, linealpha=0, label="")
        # mean line
        plot!(sp, t_hr, res.mean_ΔP .* 1e3, color=:blue, lw=1.2, label="")
    end

    push!(subplots, sp)
end

# ── Shared axis label strips ─────────────────────────────────────────────────
yl_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false,
             xlims=(0,1), ylims=(0,1))
annotate!(yl_sp, ylabel_gap, 0.5, text(main_ylabel, :black, :center, main_label_fontsize, rotation=90))

xl_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false,
             xlims=(0,1), ylims=(0,1))
annotate!(xl_sp, 0.5, xlabel_gap, text(main_xlabel, :black, :center, main_label_fontsize))

corner_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false)

# 2×2 macro layout: [ y-strip | 1×N grid ]
#                   [ corner  | x-strip  ]
l = @layout [a{0.06w} grid(n_rows, n_cols) ; b{0.001w} c{0.001h}]

fig = plot(yl_sp, subplots..., corner_sp, xl_sp,
           layout=l,
           size=(200*n_cols + 70, 160*n_rows + 50),
           plot_title="Linear Momentum Fractional Drift  ΔP/P₀  (h=$(target_h_km)km, i=$(target_i_deg)°)",
           plot_titlefontsize=maintitle_fontsize,
           left_margin=margin_left*Plots.mm,   right_margin=margin_right*Plots.mm,
           top_margin=margin_top*Plots.mm,      bottom_margin=margin_bottom*Plots.mm)

output_dir = normpath(joinpath(@__DIR__, "..", "output", "images_reconstructed_from_csv", "linear_momentum"))
mkpath(output_dir)
out_path = joinpath(output_dir, "linear_momentum_nu_plot_h$(target_h_km)km_i$(target_i_deg)deg.png")
savefig(fig, out_path)
display(fig)
println("Saved to: ", out_path)

