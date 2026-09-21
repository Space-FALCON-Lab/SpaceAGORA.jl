# NxN grid of orbital-element time series subplots.
# Rows = altitudes, Columns = inclinations.
# One figure is produced per element in `elements_to_plot`.
# All six OE (a, e, i, Ω, ω, ν) are cached together on the first call for each
# (h, i_deg) combination, so subsequent element plots load instantly from cache.

include(joinpath(@__DIR__, "csv_plotter_functions", "9_OE_function.jl"))

# ── Which elements to produce NxN figures for ─────────────────────────────────
elements_to_plot = [:a, :e, :i, :Ω]   # add :ω or :ν here as needed

# ── Layout knobs ──────────────────────────────────────────────────────────────
gap_inner        = 0     # inner subplot margin on all four sides (mm)
margin_left      = 0
margin_bottom    = 0
margin_right     = 0.5
margin_top       = 0
tick_fontsize        = 8
maintitle_fontsize   = 12
main_label_fontsize  = 9
main_xlabel          = "Number of orbits [target satellite]"
row_label_fontsize   = 8
col_label_fontsize   = 8
ylabel_gap           = 0.5
xlabel_gap           = 0.5
# ─────────────────────────────────────────────────────────────────────────────

target_altitudes_km     = [850, 950, 1000, 1050, 1150]
target_inclinations_deg = [0.0, 0.5, 1.0]

n_cols = length(target_inclinations_deg)
n_rows = length(target_altitudes_km)

output_dir = normpath(joinpath(@__DIR__, "..", "output",
                                "images_reconstructed_from_csv", "orbital_elements"))
mkpath(output_dir)

# ── Per-element metadata ──────────────────────────────────────────────────────
element_meta = Dict(
    :a => (ylabel="a, m",       title="Semi-major axis",          fname="OE_a_NxN_plot.png"),
    :e => (ylabel="e",           title="Eccentricity",             fname="OE_e_NxN_plot.png"),
    :i => (ylabel="i, deg",      title="Inclination",              fname="OE_i_NxN_plot.png"),
    :Ω => (ylabel="Ω, deg",      title="RAAN",                     fname="OE_Omega_NxN_plot.png"),
    :ω => (ylabel="ω/u, deg",    title="ω (non-circ) / u (circ)", fname="OE_omega_NxN_plot.png"),
    :ν => (ylabel="ν, deg",      title="True anomaly",             fname="OE_nu_NxN_plot.png"),
)

# ── Main loop: one NxN figure per element ─────────────────────────────────────
for element in elements_to_plot
    meta = element_meta[element]
    println("\n══ Plotting element :$element ══")

    subplots = []

    for (row_idx, h) in enumerate(target_altitudes_km)
        for (col_idx, i_deg) in enumerate(target_inclinations_deg)
            local sp
            local res
            local elapsed = @elapsed begin
                res = try
                    plot_OE(h, i_deg; element=element, show_variance=false,
                            filter_helpers=[300])
                catch e
                    println("\n=== ERROR h=$h i=$i_deg element=$element ===")
                    println(typeof(e), ": ", e)
                    Base.show_backtrace(stdout, catch_backtrace())
                    println()
                    nothing
                end
            end
            @printf("  h=%4d km  i=%.1f°  →  %.2f s\n", h, i_deg, elapsed)

            is_left   = col_idx == 1
            is_bottom = row_idx == n_rows
            row_ylab  = is_left   ? @sprintf("h=%dkm", h)           : ""
            col_xlab  = is_bottom ? @sprintf("i=%.1f\u00b0", i_deg) : ""

            if res === nothing
                sp = plot(legend=false, grid=false, framestyle=:box,
                          ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                          xlabel=col_xlab, xguidefontsize=col_label_fontsize,
                          left_margin=gap_inner*Plots.mm,  right_margin=gap_inner*Plots.mm,
                          top_margin=gap_inner*Plots.mm,   bottom_margin=gap_inner*Plots.mm)
                annotate!(sp, 0.5, 0.5, text("No Data", :grey, :center, 8))
            else
                orbits = res.orbits
                sp = plot(legend=false, grid=true,
                          ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                          xlabel=col_xlab, xguidefontsize=col_label_fontsize,
                          xtickfontsize=tick_fontsize, ytickfontsize=tick_fontsize,
                          tick_direction=:out,
                          left_margin=gap_inner*Plots.mm,  right_margin=gap_inner*Plots.mm,
                          top_margin=gap_inner*Plots.mm,   bottom_margin=gap_inner*Plots.mm)
                # single curve for N_helpers=300
                plot!(sp, orbits, res.all_vals[1], color=:blue, lw=1.2, label="")
            end

            push!(subplots, sp)
        end
    end

    # ── Shared axis label strips ──────────────────────────────────────────────
    yl_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false,
                 xlims=(0,1), ylims=(0,1))
    annotate!(yl_sp, ylabel_gap, 0.5,
              text(meta.ylabel, :black, :center, main_label_fontsize, rotation=90))

    xl_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false,
                 xlims=(0,1), ylims=(0,1))
    annotate!(xl_sp, 0.5, xlabel_gap,
              text(main_xlabel, :black, :center, main_label_fontsize))

    corner_sp = plot(framestyle=:none, ticks=nothing, legend=false, grid=false)

    l = @layout [a{0.06w} grid(n_rows, n_cols) ; b{0.001w} c{0.001h}]

    fig = plot(yl_sp, subplots..., corner_sp, xl_sp,
               layout=l,
               size=(200*n_cols + 70, 160*n_rows + 50),
               plot_title="$(meta.title) — target satellite",
               plot_titlefontsize=maintitle_fontsize,
               left_margin=margin_left*Plots.mm,    right_margin=margin_right*Plots.mm,
               top_margin=margin_top*Plots.mm,       bottom_margin=margin_bottom*Plots.mm)

    out_path = joinpath(output_dir, meta.fname)
    savefig(fig, out_path)
    display(fig)
    println("Saved: ", out_path)
end

println("\nAll OE NxN plots complete.")

