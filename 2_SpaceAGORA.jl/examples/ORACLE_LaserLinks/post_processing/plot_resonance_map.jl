#!/usr/bin/env julia

const _PP_REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
import Pkg; Pkg.activate(_PP_REPO_ROOT; io=devnull)
# Resonance prediction map for the ORACLE Case-2 laser-link scenario.
#
# For a given target altitude, computes the constellation-drift resonance
# condition for each helper altitude and plots two figures:
#
#   Figure 1 — Resonance map:
#     X-axis : number of helpers N
#     Y-axis : helper altitude
#     Vertical dashes : predicted resonant N*(k) = k × T_h/(T_h − T_t)
#     Overlaid scatter : actual ΔvT from feather files (marker size ∝ |ΔvT|)
#
#   Figure 2 — ΔvT vs N for all helper altitudes on one axis, with
#               resonance lines predicted by the formula overlaid
#
# Run:
#   julia --project=. ORACLE/post_processing/plot_resonance_map.jl
#
# The resonance formula (derived in docs/anomaly_notes.md):
#
#   Per-orbit drift of target through helper constellation:
#       δ [rad] = 2π (1 − T_t/T_h)
#   Helper angular spacing:
#       Δφ = 2π / N
#   Crossings per orbit:
#       f = N (T_h − T_t) / T_h
#   Anomaly condition (f ≈ integer k):
#       N*(k) = k × T_h / (T_h − T_t)

using Arrow
using DataFrames
get!(ENV, "GKSwstype", "100")
using Plots
using LaTeXStrings
using Printf

gr()

# ── Paths ─────────────────────────────────────────────────────────────────────
const REPO_ROOT  = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_DIR = joinpath(REPO_ROOT, "output", "paper_plot_mode")
const SAVE_DIR   = joinpath(REPO_ROOT, "output", "images")
mkpath(SAVE_DIR)

# ── Constants ─────────────────────────────────────────────────────────────────
const _MU      = 3.986004418e14
const _R_EARTH = 6_378_137.0

const TARGET_ALT_KM          = 1000.0
const TARGET_INCLINATION_DEG = 0.0
const HELPER_ALTITUDES_KM    = [1150.0, 1050.0, 950.0, 850.0]  # skip same-orbit
const INCLINATION_DELTAS_DEG = [0.0, 0.5, 1.0]
const N_MAX                  = 400

T_orbital(h_km) = 2π * sqrt((_R_EARTH + h_km * 1e3)^3 / _MU)

# ── Resonance formula ─────────────────────────────────────────────────────────
"""
    resonant_N(helper_alt_km, target_alt_km; N_max=400)

Return all resonant N*(k) values ≤ N_max, along with their k indices.
Returns `nothing` if helper and target are at the same altitude.

Formula: N*(k) = k × T_h / |T_h − T_t|
"""
function resonant_N(helper_alt_km, target_alt_km; N_max=400)
    Tt = T_orbital(target_alt_km)
    Th = T_orbital(helper_alt_km)
    abs(Th - Tt) < 1.0 && return nothing
    base = Th / abs(Th - Tt)
    ks   = Float64[]
    ns   = Float64[]
    for k in 1:10_000
        nstar = k * base
        nstar > N_max && break
        push!(ks, Float64(k))
        push!(ns, nstar)
    end
    return (base=base, k=ks, N_star=ns)
end

# ── Feather reader ────────────────────────────────────────────────────────────
function _feather_path(helper_alt_km, target_alt_km, helper_inc_deg, target_inc_deg, n)
    base_dir = joinpath(OUTPUT_DIR,
        @sprintf("h%.0fkm_t%.0fkm",    helper_alt_km,  target_alt_km),
        @sprintf("ih%.1fdeg_it%.1fdeg", helper_inc_deg, target_inc_deg),
        @sprintf("N%d", n))
    isdir(base_dir) || return nothing
    t_dirs = filter(d -> startswith(d, "T") && endswith(d, "s"), readdir(base_dir))
    isempty(t_dirs) && return nothing
    p = joinpath(base_dir, first(t_dirs), "simulation_results.feather")
    isfile(p) ? p : nothing
end

function read_dv_t(helper_alt_km, target_alt_km, helper_inc_deg, target_inc_deg, n)
    p = _feather_path(helper_alt_km, target_alt_km, helper_inc_deg, target_inc_deg, n)
    p === nothing && return nothing
    tbl = Arrow.Table(p)
    return Float64(tbl.dv_t_accumulated[end])
end

# ── Gather all available N values for each helper altitude ─────────────────────
println("Scanning output directory …")

# For each helper altitude and inclination delta: collect (N, ΔvT) pairs
data = Dict{Float64, Vector{NamedTuple}}()

for h_h in HELPER_ALTITUDES_KM
    alt_folder = @sprintf("h%.0fkm_t%.0fkm", h_h, TARGET_ALT_KM)
    alt_path   = joinpath(OUTPUT_DIR, alt_folder)
    isdir(alt_path) || continue

    for δi in INCLINATION_DELTAS_DEG
        inc_folder = @sprintf("ih%.1fdeg_it%.1fdeg",
                              TARGET_INCLINATION_DEG + δi, TARGET_INCLINATION_DEG)
        inc_path = joinpath(alt_path, inc_folder)
        isdir(inc_path) || continue

        for n_folder in readdir(inc_path)
            startswith(n_folder, "N") || continue
            n = tryparse(Int, n_folder[2:end])
            n === nothing && continue
            dv = read_dv_t(h_h, TARGET_ALT_KM,
                           TARGET_INCLINATION_DEG + δi, TARGET_INCLINATION_DEG, n)
            dv === nothing && continue
            key = h_h
            haskey(data, key) || (data[key] = NamedTuple[])
            push!(data[key], (n=n, delta_i=δi, dv_t=dv))
        end
    end
end

# ── Figure 1: Resonance map ────────────────────────────────────────────────────
# One row per helper altitude.  Grey dashes = predicted N*(k).
# Scatter = actual data, colour = inclination delta, size ∝ |ΔvT|.

INC_COLORS = Dict(0.0 => :blue, 0.5 => :red, 1.0 => :green)
n_rows     = length(HELPER_ALTITUDES_KM)
sp_map     = []

for (row, h_h) in enumerate(HELPER_ALTITUDES_KM)
    is_bottom = row == n_rows
    res = resonant_N(h_h, TARGET_ALT_KM; N_max=N_MAX)
    label_str = latexstring("h_{\\rm h}=$(round(Int,h_h))\\,{\\rm km}")

    sp = plot(
        ylabel     = label_str,
        xlabel     = is_bottom ? L"N_{\rm helpers}" : "",
        xlims      = (0, N_MAX * 1.05),
        ylims      = (-1, 1),           # normalised y (scatter overlaid later)
        yticks     = false,
        legend     = false,
        grid       = false,
        framestyle = :box,
        yguidefontsize = 11,
        xguidefontsize = 11,
        xtickfontsize  = 10,
        bottom_margin  = 2Plots.mm,
    )

    # Predicted resonance lines
    if res !== nothing
        for n_star in res.N_star
            vline!(sp, [n_star];
                   color=:black, lw=1.0, linestyle=:dash, alpha=0.4, label=false)
        end
    end

    # Actual data as scatter: y = dv_t normalised to [-1,1] within this altitude
    pts = get(data, h_h, NamedTuple[])
    if !isempty(pts)
        all_dv = [p.dv_t for p in pts]
        dv_max  = maximum(abs.(all_dv))
        dv_max  < 1e-12 && (dv_max = 1.0)

        for δi in INCLINATION_DELTAS_DEG
            subset = filter(p -> p.delta_i == δi, pts)
            isempty(subset) && continue
            sort!(subset, by=p -> p.n)
            ns  = [p.n     for p in subset]
            dvs = [p.dv_t ./ dv_max for p in subset]
            scatter!(sp, ns, dvs;
                     color=INC_COLORS[δi], ms=4, markerstrokewidth=0.2,
                     alpha=0.8, label=false)
            plot!(sp, ns, dvs;
                  color=INC_COLORS[δi], lw=1.2, alpha=0.5, label=false)
        end
        annotate!(sp, N_MAX * 0.97, 0.85,
                  text("ΔvT/max", :grey, :right, 8))
    end

    push!(sp_map, sp)
end

fig1 = plot(sp_map...;
            layout = grid(n_rows, 1),
            size   = (860, 220 * n_rows + 60),
            title  = [latexstring("Resonance map: h_{\\rm target}=$(round(Int,TARGET_ALT_KM))\\,{\\rm km}") "" "" ""],
            titlefontsize = 12,
            left_margin  = 14Plots.mm, right_margin = 4Plots.mm,
            top_margin   = 8Plots.mm, bottom_margin = 2Plots.mm)

out1 = joinpath(SAVE_DIR, "resonance_map.png")
savefig(fig1, out1)
println("Saved → $out1")

# ── Figure 2: ΔvT vs N with resonance overlay (Δi=0 only, all altitudes) ────
fig2 = plot(
    xlabel = L"N_{\rm helpers}",
    ylabel = L"\Delta v_T\ \mathrm{(m/s)}",
    title  = latexstring("\\Delta v_T\\ vs\\ N\\ with\\ predicted\\ resonances\\ " *
                         "(\\Delta i=0^\\circ,\\ h_{\\rm target}=$(round(Int,TARGET_ALT_KM))\\,{\\rm km})"),
    titlefontsize = 12,
    xlims  = (0, N_MAX * 1.05),
    legend = :outertopright,
    grid   = true,
    size   = (900, 450),
    left_margin = 10Plots.mm, right_margin = 2Plots.mm,
    bottom_margin = 6Plots.mm, top_margin = 6Plots.mm,
)

alt_colors = [:dodgerblue, :crimson, :darkorange, :purple]

for (ci, h_h) in enumerate(HELPER_ALTITUDES_KM)
    pts = filter(p -> p.delta_i == 0.0, get(data, h_h, NamedTuple[]))
    isempty(pts) && continue
    sort!(pts, by=p -> p.n)
    ns  = [p.n     for p in pts]
    dvs = [p.dv_t  for p in pts]
    col = alt_colors[ci]

    plot!(fig2, ns, dvs;
          color=col, lw=1.8, marker=:circle, ms=4, markerstrokewidth=0.2,
          label=latexstring("h_{\\rm h}=$(round(Int,h_h))\\,{\\rm km}"))

    # Resonance tick marks along the top
    res = resonant_N(h_h, TARGET_ALT_KM; N_max=N_MAX)
    if res !== nothing
        for (k_idx, n_star) in zip(res.k, res.N_star)
            vline!(fig2, [n_star];
                   color=col, lw=0.7, linestyle=:dash, alpha=0.35, label=false)
        end
    end
end

out2 = joinpath(SAVE_DIR, "resonance_overlay.png")
savefig(fig2, out2)
println("Saved → $out2")

# ── Print prediction table ─────────────────────────────────────────────────────
println("\n── Resonance prediction table (h_target=$(TARGET_ALT_KM) km) ──")
println(rpad("h_helper", 14), rpad("base T_h/(T_h-T_t)", 22), "Resonant N*(k) ≤ $(N_MAX)")
println("-"^80)
for h_h in HELPER_ALTITUDES_KM
    res = resonant_N(h_h, TARGET_ALT_KM; N_max=N_MAX)
    res === nothing && continue
    ns_str = join([@sprintf("%.0f(k=%d)", n, k) for (n, k) in zip(res.N_star, res.k)], "  ")
    println(rpad("$(h_h) km", 14), rpad(@sprintf("%.2f", res.base), 22), ns_str)
end
println("\nFormula: N*(k) = k × T_h / |T_h − T_t|")
println("Anomaly bandwidth ≈ ±(base × range_window/helper_spacing) N units")
