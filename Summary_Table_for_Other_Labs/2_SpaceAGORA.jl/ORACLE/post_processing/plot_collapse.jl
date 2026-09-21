#!/usr/bin/env julia
# Step 1 of the pattern-finding methodology:
#
# Tests the dimensionless collapse hypothesis:
#   κ = N × |T_h − T_t| / T_h
#
# If κ is the true governing variable (not N and h_h separately), then
# ΔvT vs κ curves from all helper altitudes should overlap onto a single
# master curve for each inclination delta.
#
# Outputs:
#   collapse_kappa_dvt.png    — ΔvT vs κ, all altitudes overlaid, coloured by Δi
#   collapse_heatmap.png      — ΔvT heatmap over (κ, Δi) from available data
#   collapse_report.txt       — RMS spread at each κ bin (quantifies collapse quality)
#
# Run:
#   julia --project=. ORACLE/post_processing/plot_collapse.jl

using Arrow
using DataFrames
using Statistics
get!(ENV, "GKSwstype", "100")
using Plots
using LaTeXStrings
using Printf

gr()

const REPO_ROOT  = normpath(joinpath(@__DIR__, "..", ".."))
const OUTPUT_DIR = joinpath(REPO_ROOT, "output", "paper_plot_mode")
const SAVE_DIR   = joinpath(REPO_ROOT, "output", "images")
mkpath(SAVE_DIR)

const _MU      = 3.986004418e14
const _R_EARTH = 6_378_137.0

const TARGET_ALT_KM          = 1000.0
const TARGET_INCLINATION_DEG = 0.0
const HELPER_ALTITUDES_KM    = [1150.0, 1050.0, 950.0, 850.0]
const INCLINATION_DELTAS_DEG = [0.0, 0.5, 1.0]

T_orbital(h_km) = 2π * sqrt((_R_EARTH + h_km * 1e3)^3 / _MU)

# κ = N × |T_h − T_t| / T_h  (crossings per target orbit)
kappa(n, h_h, h_t) = n * abs(T_orbital(h_h) - T_orbital(h_t)) / T_orbital(h_h)

# ── Feather reader ────────────────────────────────────────────────────────────
function read_all_scalars(output_dir, helper_alt_km, target_alt_km,
                           helper_inc_deg, target_inc_deg)
    alt_folder = @sprintf("h%.0fkm_t%.0fkm",    helper_alt_km,  target_alt_km)
    inc_folder = @sprintf("ih%.1fdeg_it%.1fdeg", helper_inc_deg, target_inc_deg)
    inc_path   = joinpath(output_dir, alt_folder, inc_folder)
    isdir(inc_path) || return NamedTuple[]

    rows = NamedTuple[]
    for n_folder in readdir(inc_path)
        startswith(n_folder, "N") || continue
        n = tryparse(Int, n_folder[2:end])
        n === nothing && continue

        t_dirs = filter(d -> startswith(d, "T") && endswith(d, "s"),
                        readdir(joinpath(inc_path, n_folder)))
        isempty(t_dirs) && continue
        fpath = joinpath(inc_path, n_folder, first(t_dirs), "simulation_results.feather")
        isfile(fpath) || continue

        tbl  = Arrow.Table(fpath)
        nrow = length(tbl.time)
        push!(rows, (
            n             = n,
            helper_alt_km = helper_alt_km,
            delta_i       = helper_inc_deg - target_inc_deg,
            kappa         = kappa(n, helper_alt_km, target_alt_km),
            dv_r          = Float64(tbl.dv_r_accumulated[nrow]),
            dv_t          = Float64(tbl.dv_t_accumulated[nrow]),
            dv_n          = Float64(tbl.dv_n_accumulated[nrow]),
        ))
    end
    return rows
end

# ── Gather all data ───────────────────────────────────────────────────────────
println("Loading feather data …")
all_rows = NamedTuple[]
for h_h in HELPER_ALTITUDES_KM
    for δi in INCLINATION_DELTAS_DEG
        rows = read_all_scalars(OUTPUT_DIR, h_h, TARGET_ALT_KM,
                                TARGET_INCLINATION_DEG + δi, TARGET_INCLINATION_DEG)
        append!(all_rows, rows)
    end
end
df = DataFrame(all_rows)
println("  $(nrow(df)) data points loaded across $(length(HELPER_ALTITUDES_KM)) altitudes × $(length(INCLINATION_DELTAS_DEG)) inclination deltas")

sort!(df, :kappa)

# ── Figure 1: Collapse plot — ΔvT vs κ ───────────────────────────────────────
# One marker per (N, h_h, Δi) combination.
# Marker shape  → helper altitude
# Colour        → inclination delta
# If curves collapse: all same-Δi markers form a single curve regardless of h_h

INC_COLORS  = Dict(0.0 => :royalblue, 0.5 => :firebrick, 1.0 => :forestgreen)
ALT_MARKERS = Dict(1150.0 => :circle,   1050.0 => :diamond,
                    950.0 => :utriangle, 850.0  => :square)
ALT_LABELS  = Dict(h => latexstring("h_{\\rm h}=$(round(Int,h))\\,{\\rm km}")
                   for h in HELPER_ALTITUDES_KM)
INC_LABELS  = Dict(δi => latexstring("\\Delta i=$(δi)^\\circ")
                   for δi in INCLINATION_DELTAS_DEG)

fig1 = plot(
    xlabel = L"\kappa = N \cdot |T_h - T_t|\,/\,T_h \quad (\mathrm{crossings/orbit})",
    ylabel = L"\Delta v_T \ \mathrm{(m/s)}",
    title  = "Dimensionless collapse: " *
             latexstring("\\Delta v_T\\ \\mathrm{vs}\\ \\kappa") *
             "\n(shape = altitude, colour = Δi — collapse if same-colour shapes overlap)",
    titlefontsize = 11,
    legend     = :outertopright,
    size       = (980, 520),
    grid       = true,
    left_margin  = 10Plots.mm, right_margin  = 2Plots.mm,
    bottom_margin = 8Plots.mm, top_margin    = 6Plots.mm,
)

# Vertical lines at integer κ (resonance positions — same for all altitudes)
for k in 1:ceil(Int, maximum(df.kappa))
    vline!(fig1, [Float64(k)]; color=:gray, lw=0.6, linestyle=:dot, alpha=0.5, label=false)
end
annotate!(fig1, 0.5, maximum(df.dv_t) * 0.92,
          text("integer κ = resonance", :grey, :left, 8))

# Plot each (altitude × inclination) series
first_alt_plotted = Set{Float64}()
first_inc_plotted = Set{Float64}()
for h_h in HELPER_ALTITUDES_KM
    for δi in INCLINATION_DELTAS_DEG
        sub = sort(filter(r -> r.helper_alt_km == h_h && r.delta_i == δi, eachrow(df)),
                   by=r -> r.kappa)
        isempty(sub) && continue
        ks  = [r.kappa for r in sub]
        dvs = [r.dv_t  for r in sub]

        # Legend entries: first time we see this altitude or inclination
        alt_lab = h_h ∉ first_alt_plotted ? ALT_LABELS[h_h] : false
        inc_lab = δi ∉ first_inc_plotted  ? INC_LABELS[δi]  : false
        push!(first_alt_plotted, h_h); push!(first_inc_plotted, δi)

        scatter!(fig1, ks, dvs;
                 color=INC_COLORS[δi], marker=ALT_MARKERS[h_h],
                 ms=5, markerstrokewidth=0.3, alpha=0.85,
                 label=false)
        plot!(fig1, ks, dvs;
              color=INC_COLORS[δi], lw=0.9, alpha=0.5, linestyle=:solid,
              label=false)
    end
end

# Manual legend: shape = altitude
for h_h in HELPER_ALTITUDES_KM
    scatter!(fig1, [NaN], [NaN];
             color=:black, marker=ALT_MARKERS[h_h], ms=5,
             label=ALT_LABELS[h_h])
end
# Manual legend: colour = inclination
for δi in INCLINATION_DELTAS_DEG
    plot!(fig1, [NaN, NaN], [NaN, NaN];
          color=INC_COLORS[δi], lw=2.0,
          label=INC_LABELS[δi])
end

out1 = joinpath(SAVE_DIR, "collapse_kappa_dvt.png")
savefig(fig1, out1)
println("Saved → $out1")

# ── Figure 2: Heatmap of ΔvT over (κ, Δi) ────────────────────────────────────
# Bin κ into uniform cells and average ΔvT within each cell.
κ_bins  = range(0.0, ceil(maximum(df.kappa)) + 0.5; step=0.25)
Δi_vals = sort(unique(df.delta_i))
n_κ     = length(κ_bins) - 1
n_Δi    = length(Δi_vals)

Z = fill(NaN, n_κ, n_Δi)
for (ji, δi) in enumerate(Δi_vals)
    sub = filter(r -> r.delta_i == δi, eachrow(df))
    for (ki, κ_lo) in enumerate(κ_bins[1:end-1])
        κ_hi  = κ_bins[ki+1]
        in_bin = filter(r -> κ_lo <= r.kappa < κ_hi, sub)
        isempty(in_bin) || (Z[ki, ji] = mean(r.dv_t for r in in_bin))
    end
end

κ_centres = [(κ_bins[i] + κ_bins[i+1]) / 2 for i in 1:n_κ]

fig2 = heatmap(
    Δi_vals, κ_centres, Z;
    xlabel = L"\Delta i\ \mathrm{(deg)}",
    ylabel = L"\kappa = N \cdot |T_h - T_t|\,/\,T_h",
    title  = latexstring("\\langle\\Delta v_T\\rangle\\ \\mathrm{(m/s)}\\ \\mathrm{binned\\ over\\ }(\\kappa,\\ \\Delta i)"),
    titlefontsize = 11,
    colorbar_title = "ΔvT (m/s)",
    c      = :RdBu,
    clims  = let m = maximum(abs.(filter(!isnan, Z)));(-m, m) end,
    size   = (560, 680),
    left_margin  = 12Plots.mm, right_margin  = 12Plots.mm,
    bottom_margin = 8Plots.mm, top_margin    = 6Plots.mm,
    xtickfontsize=11, ytickfontsize=11,
    xguidefontsize=12, yguidefontsize=12,
)
# Integer κ lines = resonances
for k in 1:ceil(Int, maximum(df.kappa))
    hline!(fig2, [Float64(k)]; color=:black, lw=0.8, linestyle=:dash,
           alpha=0.5, label=false)
end

out2 = joinpath(SAVE_DIR, "collapse_heatmap.png")
savefig(fig2, out2)
println("Saved → $out2")

# ── Collapse quality report ───────────────────────────────────────────────────
# For each (κ-bin, Δi), compute RMS spread across altitudes.
# Low RMS = good collapse (altitude doesn't matter beyond κ).
report_path = joinpath(SAVE_DIR, "collapse_report.txt")
open(report_path, "w") do io
    println(io, "Collapse quality report")
    println(io, "  If RMS spread is small relative to the mean, κ alone captures the altitude effect.")
    println(io, "  Large RMS spread → altitude enters through a variable OTHER than κ.\n")
    println(io, @sprintf("%-12s %-10s %-10s %-10s %-8s", "κ-bin centre", "Δi", "mean ΔvT", "RMS spread", "n_pts"))
    println(io, "-"^58)
    for (ji, δi) in enumerate(Δi_vals)
        for (ki, κ_lo) in enumerate(κ_bins[1:end-1])
            κ_hi  = κ_bins[ki+1]
            sub   = filter(r -> r.delta_i == δi && κ_lo <= r.kappa < κ_hi, eachrow(df))
            length(sub) < 2 && continue
            dvs   = [r.dv_t for r in sub]
            μ     = mean(dvs); σ = std(dvs)
            println(io, @sprintf("%-12.2f %-10.1f %-10.4f %-10.4f %-8d",
                                  (κ_lo+κ_hi)/2, δi, μ, σ, length(sub)))
        end
    end
end
println("Saved → $report_path")

# ── Console summary ───────────────────────────────────────────────────────────
println("\n── Collapse quality (RMS spread / |mean|, lower = better collapse) ──")
println(rpad("κ range", 16), rpad("Δi=0.0°", 14), rpad("Δi=0.5°", 14), "Δi=1.0°")
println("-"^58)
for (ki, κ_lo) in enumerate(κ_bins[1:end-1])
    κ_hi = κ_bins[ki+1]
    parts = String[]
    any_data = false
    for δi in Δi_vals
        sub = filter(r -> r.delta_i == δi && κ_lo <= r.kappa < κ_hi, eachrow(df))
        if length(sub) < 2
            push!(parts, rpad("—", 14))
        else
            dvs = [r.dv_t for r in sub]
            μ   = mean(dvs); σ = std(dvs)
            rel = abs(μ) > 1e-9 ? σ/abs(μ) : NaN
            push!(parts, rpad(@sprintf("%.2f (n=%d)", isnan(rel) ? 0.0 : rel, length(sub)), 14))
            any_data = true
        end
    end
    any_data || continue
    println(rpad(@sprintf("[%.1f, %.1f)", κ_lo, κ_hi), 16), join(parts))
end
println("\nInterpretation:")
println("  RMS/|mean| < 0.1  → excellent collapse; κ fully captures altitude effect")
println("  RMS/|mean| 0.1–0.5 → partial collapse; altitude has secondary role")
println("  RMS/|mean| > 0.5  → poor collapse; need additional dimensionless variable")
println("\nNext steps based on collapse quality:")
println("  Good:  Fit surrogate as f(κ, Δi) only — 2D problem")
println("  Poor:  Add altitude ratio h_h/h_t as third axis — 3D problem")
println("         Then run: Δi sweep at fixed (κ, h_h/h_t) to map optimal inclination")
