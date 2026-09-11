#!/usr/bin/env julia
# Plots one selected ΔV component (R, T, or N) as a 1×5 row across all
# target altitudes — equivalent to extracting one row from plot_dv_from_feather.jl.
#
# Columns = target altitude h_t  [km]
# X-axis  = N_helpers  (one curve per inclination delta)
#
# Run from the repository root:
#   julia --project=. ORACLE/post_processing/plot_dv_row.jl \
#       --schedule gve_inc --target-ecc 0.0 --target-nu-deg 1.0 --component N

using Arrow
get!(ENV, "GKSwstype", "100")
using Plots
using LaTeXStrings
using Printf

gr()

const REPO_ROOT  = normpath(joinpath(@__DIR__, "..", ".."))
const OUTPUT_DIR = joinpath(REPO_ROOT, "output", "paper_plot_mode")

# ── CLI ───────────────────────────────────────────────────────────────────────
_sched_idx = findfirst(==("--schedule"), ARGS)
const SCHEDULE = _sched_idx !== nothing ? ARGS[_sched_idx + 1] : "naive_next_entering"

_ecc_idx = findfirst(==("--target-ecc"), ARGS)
const TARGET_ECC = _ecc_idx !== nothing ? parse(Float64, ARGS[_ecc_idx + 1]) : 0.0

_nu_idx = findfirst(==("--target-nu-deg"), ARGS)
const TARGET_NU_DEG = _nu_idx !== nothing ? parse(Float64, ARGS[_nu_idx + 1]) : 0.0

_comp_idx = findfirst(==("--component"), ARGS)
const COMPONENT = _comp_idx !== nothing ? uppercase(ARGS[_comp_idx + 1]) : "T"
COMPONENT in ("R", "T", "N") || error("--component must be R, T, or N")

_w_idx = findfirst(==("--width"), ARGS)
_h_idx = findfirst(==("--height"), ARGS)
const FIG_HEIGHT_OVERRIDE = _h_idx !== nothing ? parse(Int, ARGS[_h_idx + 1]) : nothing
const FIG_WIDTH_OVERRIDE  = _w_idx !== nothing ? parse(Int, ARGS[_w_idx + 1]) : nothing

# ── Grid parameters ───────────────────────────────────────────────────────────
const HELPER_ALT_KM          = 1000.0
const TARGET_INCLINATION_DEG = 0.0
const TARGET_ALTITUDES_KM    = [1150.0, 1050.0, 1000.0, 950.0, 850.0]
const INCLINATION_DELTAS_DEG = [0.0, 0.5, 1.0]
const HELPER_COUNTS          = [1, 50, 100, 200, 300]

# ── Feather reader (identical to plot_dv_from_feather.jl) ────────────────────
function read_final_dv(output_dir, helper_alt_km, target_alt_km,
                        helper_inc_deg, target_inc_deg, n_helpers;
                        schedule::String="", target_ecc::Float64=0.0,
                        target_nu_deg::Float64=0.0)
    alt_folder = @sprintf("h%.0fkm_t%.0fkm",    helper_alt_km,  target_alt_km)
    inc_folder = @sprintf("ih%.1fdeg_it%.1fdeg", helper_inc_deg, target_inc_deg)
    n_folder   = @sprintf("N%d",                 n_helpers)
    base_path  = joinpath(output_dir, alt_folder, inc_folder, n_folder)
    isdir(base_path) || return nothing

    t_folders = filter(d -> startswith(d, "T") && endswith(d, "s"), readdir(base_path))
    isempty(t_folders) && return nothing

    sched_folder = isempty(schedule) ? "" :
        @sprintf("%s_e%.4f_nu%.4f", schedule, target_ecc, target_nu_deg)
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

# ── Styling ───────────────────────────────────────────────────────────────────
gap_inner           = 0
tick_fontsize       = 12
main_label_fontsize = 12
col_title_fontsize  = 12
row_label_fontsize  = 12
legend_fontsize     = 12
marker_stroke_width = 0.1
main_xlabel         = L"N_{\mathrm{helpers}}"

comp_label = COMPONENT == "R" ? L"\Delta v_R\ \mathrm{(m/s)}" :
             COMPONENT == "T" ? L"\Delta v_T\ \mathrm{(m/s)}" :
                                L"\Delta v_N\ \mathrm{(m/s)}"

inc_colors = [:blue, :red, :green]
inc_labels = [latexstring("\\Delta i=$(round(d, digits=1))^{\\circ}")
              for d in INCLINATION_DELTAS_DEG]

n_cols    = length(TARGET_ALTITUDES_KM)
fig_width = 220 * n_cols + 80

# ── nice_ticks (same as parent script) ───────────────────────────────────────
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

# ── Build the 1×5 subplots ───────────────────────────────────────────────────
subplots = []

for (col_idx, t_alt) in enumerate(TARGET_ALTITUDES_KM)
    is_left = col_idx == 1

    all_xs  = Vector{Vector{Float64}}()
    all_dvs = Vector{Vector{Float64}}()

    for delta_i in INCLINATION_DELTAS_DEG
        t_inc = TARGET_INCLINATION_DEG + delta_i
        xs  = Float64[]
        dvs = Float64[]
        for n in HELPER_COUNTS
            res = read_final_dv(OUTPUT_DIR, HELPER_ALT_KM, t_alt,
                                TARGET_INCLINATION_DEG, t_inc, n;
                                schedule=SCHEDULE, target_ecc=TARGET_ECC,
                                target_nu_deg=TARGET_NU_DEG)
            res === nothing && continue
            push!(xs,  Float64(n))
            push!(dvs, COMPONENT == "R" ? res.dv_r :
                       COMPONENT == "T" ? res.dv_t : res.dv_n)
        end
        push!(all_xs,  xs)
        push!(all_dvs, dvs)
    end

    all_vals  = vcat(all_dvs...)
    col_title = latexstring("h_{\\mathrm{t}}=$(round(Int,t_alt))\\,\\mathrm{km}")
    row_ylab  = is_left ? comp_label : ""

    if isempty(all_vals)
        sp = plot(legend=false, grid=false, framestyle=:box,
                  title=col_title, titlefontsize=col_title_fontsize,
                  ylabel=row_ylab, yguidefontsize=row_label_fontsize,
                  xlabel=main_xlabel, xguidefontsize=main_label_fontsize,
                  left_margin=gap_inner*Plots.mm, right_margin=gap_inner*Plots.mm,
                  top_margin=gap_inner*Plots.mm,  bottom_margin=gap_inner*Plots.mm)
        annotate!(sp, 0.5, 0.5, text("No Data", :grey, :center, 8))
        push!(subplots, sp)
        continue
    end

    _vmin  = minimum(all_vals); _vmax = maximum(all_vals)
    _vctr  = (_vmin + _vmax) / 2
    _vhalf = max((_vmax - _vmin) / 2 * 1.1, abs(_vctr) * 0.5, 1e-9)
    ytk_vals, ytk_labs = nice_ticks(_vctr - _vhalf, _vctr + _vhalf; n=5)

    ref_xs = first(filter(!isempty, all_xs))

    sp = plot(legend=false, grid=true,
              title=col_title, titlefontsize=col_title_fontsize,
              ylabel=row_ylab, yguidefontsize=row_label_fontsize,
              xlabel=main_xlabel, xguidefontsize=main_label_fontsize,
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

# ── Title strip ─────────────────────────────────────────────────────────────
sched_display = replace(SCHEDULE, "_" => "\\_")
comp_tex = COMPONENT == "R" ? "\\Delta v_R" :
           COMPONENT == "T" ? "\\Delta v_T" : "\\Delta v_N"
title_str = latexstring(
    "\\mathrm{$sched_display}\\ \\ $comp_tex\\ \\ " *
    "e=$(TARGET_ECC)\\ \\ \\nu=$(TARGET_NU_DEG)^{\\circ}"
)
title_sp = plot(framestyle=:none, ticks=nothing, legend=false,
                grid=false, xlims=(0,1), ylims=(0,1))
annotate!(title_sp, 0.5, 0.4, text(title_str, :black, :center, 12))

# ── Legend strip ──────────────────────────────────────────────────────────────
legend_sp = plot(framestyle=:none, ticks=nothing, legend=false,
                 grid=false, xlims=(0,1), ylims=(0,1))

legend_xs = [0.15, 0.45, 0.75]
for (xc, col, lab) in zip(legend_xs, inc_colors, inc_labels)
    plot!(legend_sp, [xc, xc+0.06], [0.6, 0.6]; color=col, lw=1.5, linestyle=:solid)
    scatter!(legend_sp, [xc+0.03], [0.6]; color=col, marker=:circle, ms=3,
             markerstrokewidth=marker_stroke_width, label=false)
    annotate!(legend_sp, xc+0.07, 0.6, text(lab, :black, :left, legend_fontsize))
end

# ── Assemble ──────────────────────────────────────────────────────────────────
l = @layout [a{0.08h}
             grid(1, n_cols)
             b{0.12h}]

fig = plot(title_sp, subplots..., legend_sp,
           layout=l,
           size=(FIG_WIDTH_OVERRIDE !== nothing ? FIG_WIDTH_OVERRIDE : fig_width,
                 FIG_HEIGHT_OVERRIDE !== nothing ? FIG_HEIGHT_OVERRIDE : 370),
           left_margin=5*Plots.mm,  right_margin=2*Plots.mm,
           top_margin=2*Plots.mm,   bottom_margin=8*Plots.mm)

# ── Save ──────────────────────────────────────────────────────────────────────
out_dir = joinpath(REPO_ROOT, "output", "images", "delta_v")
mkpath(out_dir)

_sched_suffix = isempty(SCHEDULE) ? "" : "_$(SCHEDULE)"
_ecc_suffix   = @sprintf("_e%.4f", TARGET_ECC)
_nu_suffix    = @sprintf("_nu%.4f", TARGET_NU_DEG)
stem = "dv_row$(COMPONENT)$(_sched_suffix)$(_ecc_suffix)$(_nu_suffix)_h_cols"

savefig(fig, joinpath(out_dir, stem * ".png"))
savefig(fig, joinpath(out_dir, stem * ".pdf"))
display(fig)
println("Saved → ", joinpath(out_dir, stem * ".pdf"))
