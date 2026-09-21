#!/usr/bin/env julia

const _PP_REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
import Pkg; Pkg.activate(_PP_REPO_ROOT; io=devnull)
# Reads simulation_results.feather files produced by run_case2_laser_links.jl
# (--paper-grid --feather-only) and plots a 3×5 ΔvRTN grid.
#
# Rows    = [ΔvR, ΔvT, ΔvN]
# Columns = target altitude h_t  [km]
# X-axis  = N_helpers  (one curve per inclination delta)
#
# Run from the repository root:
#   julia --project=. ORACLE/post_processing/plot_dv_from_feather.jl

using Arrow
using DataFrames
get!(ENV, "GKSwstype", "100")   # headless rendering (remove if you want a window)
using Plots
using LaTeXStrings
using Printf

gr()

const REPO_ROOT  = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_DIR = joinpath(REPO_ROOT, "output", "paper_plot_mode")

# ── Schedule argument: --schedule positive_along_track (default: naive_next_entering) ─
_sched_idx = findfirst(==("--schedule"), ARGS)
const SCHEDULE = _sched_idx !== nothing ? ARGS[_sched_idx + 1] : "naive_next_entering"

# ── Target eccentricity argument: --target-ecc 0.01 (default: 0.0) ───────────
_ecc_idx = findfirst(==("--target-ecc"), ARGS)
const TARGET_ECC = _ecc_idx !== nothing ? parse(Float64, ARGS[_ecc_idx + 1]) : 0.0

# ── Target true anomaly argument: --target-nu-deg 1.0 (default: 0.0) ────────
_nu_idx = findfirst(==("--target-nu-deg"), ARGS)
const TARGET_NU_DEG = _nu_idx !== nothing ? parse(Float64, ARGS[_nu_idx + 1]) : 0.0

# ── Grid parameters — must match what was used in the simulation sweep ────────
const HELPER_ALT_KM          = 1000.0
const TARGET_INCLINATION_DEG = 0.0
const TARGET_ALTITUDES_KM    = [1150.0, 1050.0, 1000.0, 950.0, 850.0]
const INCLINATION_DELTAS_DEG = [0.0, 0.5, 1.0]
const HELPER_COUNTS          = [1, 50, 100, 150, 200, 250, 300]   # must match PAPER_HELPER_COUNTS

# ── Read the final cumulative ΔV for one scenario from its feather file ───────
# Returns (dv_r, dv_t, dv_n) or nothing if the file is not found.
function read_final_dv(output_dir, helper_alt_km, target_alt_km,
                        helper_inc_deg, target_inc_deg, n_helpers;
                        schedule::String="", target_ecc::Float64=0.0,
                        target_nu_deg::Float64=0.0)
    alt_folder = @sprintf("h%.0fkm_t%.0fkm",     helper_alt_km,  target_alt_km)
    inc_folder = @sprintf("ih%.1fdeg_it%.1fdeg",  helper_inc_deg, target_inc_deg)
    n_folder   = @sprintf("N%d",                  n_helpers)
    base_path  = joinpath(output_dir, alt_folder, inc_folder, n_folder)
    isdir(base_path) || return nothing

    # T_s folder name depends on orbits × target period — find it dynamically
    t_folders = filter(d -> startswith(d, "T") && endswith(d, "s"),
                       readdir(base_path))
    isempty(t_folders) && return nothing

    # Schedule folder: always use _e{ecc}_nu{nu} suffix (matches renamed folder convention)
    sched_folder = if isempty(schedule)
        ""
    else
        @sprintf("%s_e%.4f_nu%.4f", schedule, target_ecc, target_nu_deg)
    end
    t_root = joinpath(base_path, first(t_folders))
    feather_path = isempty(sched_folder) ?
        joinpath(t_root, "simulation_results.feather") :
        joinpath(t_root, sched_folder, "simulation_results.feather")
    isfile(feather_path) || return nothing

    tbl = Arrow.Table(feather_path)
    n   = length(tbl.time)
    return (
        dv_r = Float64(tbl.dv_r_accumulated[n]),
        dv_t = Float64(tbl.dv_t_accumulated[n]),
        dv_n = Float64(tbl.dv_n_accumulated[n]),
    )
end

# ── Nice axis ticks (ported from the prototype) ───────────────────────────────
function nice_ticks(lo, hi; n=5)
    span = hi - lo
    span == 0 && (span = abs(lo) > 0 ? abs(lo) * 0.2 : 1.0)
    raw  = span / (n - 1)
    mag  = 10.0^floor(log10(raw))
    for mult in (1.0, 2.0, 5.0, 10.0, 20.0, 50.0)
        step  = mult * mag
        hi_t  = ceil(hi  / step) * step
        lo_t  = hi_t - (n - 1) * step
        lo_t <= lo || continue
        vals  = [lo_t + k * step for k in 0:(n-1)]
        ndec  = max(0, -floor(Int, log10(step)))
        labs  = step >= 1 ? [latexstring(string(round(Int, v))) for v in vals] :
                             [latexstring(@sprintf("%.*f", ndec, v)) for v in vals]
        return vals, labs
    end
    vals = collect(range(lo, hi, length=n))
    return vals, [latexstring(@sprintf("%.2g", v)) for v in vals]
end

# ── Styling knobs (match the prototype) ──────────────────────────────────────
gap_inner           = 0
margin_left         = 5
margin_bottom       = 0
margin_right        = 0
margin_top          = 0
tick_fontsize       = 12
main_label_fontsize = 12
col_title_fontsize  = 12
row_label_fontsize  = 12
legend_fontsize     = 12
marker_stroke_width = 0.1
main_xlabel         = L"N_{\mathrm{helpers}}"

components  = [:R, :T, :N]
comp_labels = [L"\Delta v_R\ \mathrm{(m/s)}",
               L"\Delta v_T\ \mathrm{(m/s)}",
               L"\Delta v_N\ \mathrm{(m/s)}"]

inc_colors = [:blue, :red, :green]
inc_labels = [latexstring("\\Delta i=$(round(d, digits=1))^{\\circ}")
              for d in INCLINATION_DELTAS_DEG]

# columns correspond to target altitudes
n_rows     = length(components)           # 3  (ΔvR, ΔvT, ΔvN)
n_cols     = length(TARGET_ALTITUDES_KM)  # 5  (target altitudes)
fig_width  = 220 * n_cols + 80
fig_height = 180 * n_rows + 50

# ── Build subplots (row-major: component × h_t) ───────────────────────────────
subplots = []

for (row_idx, comp) in enumerate(components)
    is_top    = row_idx == 1
    is_bottom = row_idx == n_rows

    for (col_idx, t_alt) in enumerate(TARGET_ALTITUDES_KM)
        is_left = col_idx == 1

        # collect one (xs, dvs) vector per inclination curve
        all_xs  = Vector{Vector{Float64}}()
        all_dvs = Vector{Vector{Float64}}()

        for delta_i in INCLINATION_DELTAS_DEG
            t_inc = TARGET_INCLINATION_DEG + delta_i
            xs    = Float64[]
            dvs   = Float64[]
            for n in HELPER_COUNTS
                res = read_final_dv(OUTPUT_DIR, HELPER_ALT_KM, t_alt,
                                    TARGET_INCLINATION_DEG, t_inc, n;
                                    schedule=SCHEDULE, target_ecc=TARGET_ECC,
                                    target_nu_deg=TARGET_NU_DEG)
                res === nothing && continue
                push!(xs,  Float64(n))
                push!(dvs, comp == :R ? res.dv_r :
                            comp == :T ? res.dv_t : res.dv_n)
            end
            push!(all_xs,  xs)
            push!(all_dvs, dvs)
        end

        all_vals  = vcat(all_dvs...)
        col_title = is_top    ? latexstring("h_{\\mathrm{t}}=$(round(Int,t_alt))\\,\\mathrm{km}") : ""
        row_ylab  = is_left   ? comp_labels[row_idx] : ""
        x_lab     = is_bottom ? main_xlabel : ""

        # ── No-data guard ────────────────────────────────────────────────────
        if isempty(all_vals)
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
        _vmin  = minimum(all_vals); _vmax = maximum(all_vals)
        _vctr  = (_vmin + _vmax) / 2
        _vhalf = max((_vmax - _vmin) / 2 * 1.1, abs(_vctr) * 0.5, 1e-9)
        ytk_vals, ytk_labs = nice_ticks(_vctr - _vhalf, _vctr + _vhalf; n=5)

        ref_xs = first(filter(!isempty, all_xs))

        sp = plot(legend=false, grid=true,
                  title=col_title,  titlefontsize=col_title_fontsize,
                  ylabel=row_ylab,  yguidefontsize=row_label_fontsize,
                  xlabel=x_lab,     xguidefontsize=main_label_fontsize,
                  xtickfontsize=tick_fontsize, ytickfontsize=tick_fontsize,
                  tick_direction=:out,
                  xticks=(ref_xs, latexstring.(string.(round.(Int, ref_xs)))),
                  xrotation=45,
                  ylims=(ytk_vals[1], ytk_vals[end]),
                  yticks=(ytk_vals, ytk_labs),
                  left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                  top_margin=gap_inner*Plots.mm,  bottom_margin=gap_inner*Plots.mm)

        for (k, (xs, dvs)) in enumerate(zip(all_xs, all_dvs))
            isempty(xs) && continue
            plot!(sp, xs, dvs;
                  color=inc_colors[k], linestyle=:solid, lw=1.5,
                  marker=:circle, ms=3,
                  markerstrokewidth=marker_stroke_width,
                  label=false)
        end

        push!(subplots, sp)
    end
end

# ── Horizontal legend strip ───────────────────────────────────────────────────
legend_sp = plot(framestyle=:none, ticks=nothing, legend=false,
                 grid=false, xlims=(0,1), ylims=(0,1))

legend_xs = [0.15, 0.45, 0.75]
for (xc, col, lab) in zip(legend_xs, inc_colors, inc_labels)
    plot!(legend_sp, [xc, xc+0.06], [0.6, 0.6];
          color=col, lw=1.5, linestyle=:solid)
    scatter!(legend_sp, [xc+0.03], [0.6];
             color=col, marker=:circle, ms=3,
             markerstrokewidth=marker_stroke_width, label=false)
    annotate!(legend_sp, xc+0.07, 0.6,
              text(lab, :black, :left, legend_fontsize))
end

# ── Assemble and save ─────────────────────────────────────────────────────────
l = @layout [grid(n_rows, n_cols)
             b{0.07h}]

fig = plot(subplots..., legend_sp,
           layout=l,
           size=(fig_width, fig_height),
           left_margin=margin_left*Plots.mm,  right_margin=margin_right*Plots.mm,
           top_margin=margin_top*Plots.mm,    bottom_margin=margin_bottom*Plots.mm)

out_dir  = joinpath(REPO_ROOT, "output", "images", "delta_v")
mkpath(out_dir)
_sched_suffix = isempty(SCHEDULE) ? "" : "_$(SCHEDULE)"
_ecc_suffix   = @sprintf("_e%.4f", TARGET_ECC)
_nu_suffix    = @sprintf("_nu%.4f", TARGET_NU_DEG)
out_png  = joinpath(out_dir, "dv_RTN_from_feather$(_sched_suffix)$(_ecc_suffix)$(_nu_suffix)_h_cols.png")
out_pdf  = joinpath(out_dir, "dv_RTN_from_feather$(_sched_suffix)$(_ecc_suffix)$(_nu_suffix)_h_cols.pdf")
savefig(fig, out_png)
savefig(fig, out_pdf)
display(fig)
println("Saved → ", out_png)
println("Saved → ", out_pdf)
