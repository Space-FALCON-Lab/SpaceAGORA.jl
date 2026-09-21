#!/usr/bin/env julia
# Verifies the stroboscopic resonance mechanism using actual simulation feather data.
#
# For N=150, 200, 250 at the anomalous helper altitude, plots:
#
#   Figure 1 — Orbital phase scatter:
#     X-axis : orbit number  (time / T_target)
#     Y-axis : target's true anomaly at each activation switch
#     If resonant: points form near-vertical stripes (same phase repeating)
#     If non-resonant: diagonal bands (phase slowly drifts through full circle)
#
#   Figure 2 — Polar phase histogram:
#     Angular position = true anomaly at activation
#     Radius          = count of activations at that phase
#     Resonant: 6 distinct spokes.  Non-resonant: uniform ring.
#
#   Figure 3 — Per-activation ΔvT vs orbital phase:
#     Each point = one detected activation; y = ΔvT increment during that activation
#     Shows whether the phases locked by resonance have consistent sign
#
# Run:
#   julia --project=. ORACLE/post_processing/plot_phase_scatter.jl

using Arrow
using DataFrames
using LinearAlgebra
using StaticArrays
using Statistics
get!(ENV, "GKSwstype", "100")
using Plots
using LaTeXStrings
using Printf

gr()

# ── Config ─────────────────────────────────────────────────────────────────────
const REPO_ROOT  = normpath(joinpath(@__DIR__, "..", ".."))
const OUTPUT_DIR = joinpath(REPO_ROOT, "output", "paper_plot_mode")
const SAVE_DIR   = joinpath(REPO_ROOT, "output", "images")
mkpath(SAVE_DIR)

const HELPER_ALT_KM    = 1050.0
const TARGET_ALT_KM    = 1000.0
const DELTA_I_SWEEP    = [0.0, 0.5, 1.0]   # generate one set of plots per Δi
const N_COMPARE        = [50, 100, 150]

const _MU      = 3.986004418e14
const _R_EARTH = 6_378_137.0

T_orbital(h_km) = 2π * sqrt((_R_EARTH + h_km * 1e3)^3 / _MU)
const T_TARGET  = T_orbital(TARGET_ALT_KM)
const T_HELPER  = T_orbital(HELPER_ALT_KM)
const BASE      = T_HELPER / abs(T_HELPER - T_TARGET)

# ── Orbital mechanics helpers ─────────────────────────────────────────────────

# True anomaly [rad, 0–2π] from ECI position and velocity
function true_anomaly(r_vec, v_vec; mu=_MU)
    r  = norm(r_vec)
    v  = norm(v_vec)
    h  = cross(r_vec, v_vec)         # specific angular momentum
    n  = cross(SVector(0.,0.,1.), h) # node vector
    e_vec = ((v^2 - mu/r) .* r_vec .- dot(r_vec, v_vec) .* v_vec) ./ mu
    e  = norm(e_vec)
    if e < 1e-10                     # circular orbit: use argument of latitude
        cosν = dot(n, r_vec) / (norm(n) * r)
        cosν = clamp(cosν, -1.0, 1.0)
        ν    = acos(cosν)
        return dot(n, v_vec) >= 0 ? ν : (2π - ν)
    end
    cosν = clamp(dot(e_vec, r_vec) / (e * r), -1.0, 1.0)
    ν    = acos(cosν)
    return dot(r_vec, v_vec) >= 0 ? ν : (2π - ν)
end

# ── Feather loader and activation extractor ───────────────────────────────────
function _feather_path(h_h, h_t, helper_inc, target_inc, n)
    base_dir = joinpath(OUTPUT_DIR,
        @sprintf("h%.0fkm_t%.0fkm",    h_h,        h_t),
        @sprintf("ih%.1fdeg_it%.1fdeg", helper_inc, target_inc),
        @sprintf("N%d", n))
    isdir(base_dir) || return nothing
    dirs = filter(d -> startswith(d,"T") && endswith(d,"s"), readdir(base_dir))
    isempty(dirs) && return nothing
    p = joinpath(base_dir, first(dirs), "simulation_results.feather")
    isfile(p) ? p : nothing
end

"""
Extract activation events from a feather table.
Returns a NamedTuple vector with fields:
  t_start, t_end   — simulation time window [s]
  orbit_start      — orbit number at t_start
  true_anom_deg    — target's true anomaly at activation midpoint [deg]
  dv_t_increment   — ΔvT accumulated during this activation [m/s]
  helper_idx       — which helper was active
"""
function extract_activations(tbl)
    n       = length(tbl.time)
    helper  = Int.(tbl.laser_active_helper)
    t       = Float64.(tbl.time)
    dv_t    = Float64.(tbl.dv_t_accumulated)

    events = NamedTuple[]
    i = 1
    while i <= n
        # Find the start of an activation window (helper > 0)
        helper[i] == 0 && (i += 1; continue)
        j_start = i
        h_id    = helper[i]

        # Advance while same helper is active
        j_end = j_start
        while j_end + 1 <= n && helper[j_end+1] == h_id
            j_end += 1
        end

        # Midpoint index (for orbital phase)
        j_mid   = (j_start + j_end) ÷ 2
        t_mid   = t[j_mid]

        # True anomaly at midpoint (target = spacecraft 1)
        r_vec = SVector{3,Float64}(
            tbl.sc1_pos_1[j_mid], tbl.sc1_pos_2[j_mid], tbl.sc1_pos_3[j_mid])
        v_vec = SVector{3,Float64}(
            tbl.sc1_vel_1[j_mid], tbl.sc1_vel_2[j_mid], tbl.sc1_vel_3[j_mid])
        ν_deg = rad2deg(true_anomaly(r_vec, v_vec))

        # ΔvT increment across this activation
        dv_inc = dv_t[j_end] - dv_t[j_start]

        push!(events, (
            t_start       = t[j_start],
            t_end         = t[j_end],
            orbit_start   = t[j_start] / T_TARGET,
            true_anom_deg = ν_deg,
            dv_t_increment= dv_inc,
            helper_idx    = h_id,
        ))

        i = j_end + 1
    end
    return events
end

_PALETTE  = [:dodgerblue, :crimson, :darkorange, :purple, :green, :brown]
TS_COLORS = Dict(n => _PALETTE[mod1(i, length(_PALETTE))]
                 for (i, n) in enumerate(N_COMPARE))

for delta_i in DELTA_I_SWEEP
    # Δi label used in titles and filenames (e.g. "di0.5")
    di_tag   = @sprintf("di%.1f", delta_i)
    di_label = latexstring("\\Delta i=$(delta_i)^\\circ")

    # ── Load data ─────────────────────────────────────────────────────────────
    println("\nLoading feather data for N = $N_COMPARE at h_h=$(HELPER_ALT_KM) km, Δi=$(delta_i)° …")

    all_data = Dict{Int, Vector}()
    for n in N_COMPARE
        p = _feather_path(HELPER_ALT_KM, TARGET_ALT_KM, delta_i, 0.0, n)
        if p === nothing
            @warn "No feather found for N=$n Δi=$(delta_i)° — skipping"
            continue
        end
        tbl    = Arrow.Table(p)
        events = extract_activations(tbl)
        all_data[n] = events
        κ = n / BASE
        println("  N=$n (κ=$(round(κ, digits=2))): $(length(events)) activation events detected")
    end
    isempty(all_data) && continue

    # ── Figure 1: Phase vs orbit scatter ──────────────────────────────────────
    panels_fig1 = []
    for (idx, n) in enumerate(N_COMPARE)
        haskey(all_data, n) || continue
        events = all_data[n]
        κ      = n / BASE

        orbits = [e.orbit_start   for e in events]
        phases = [e.true_anom_deg for e in events]

        sp = scatter(orbits, phases;
            xlabel     = idx == length(N_COMPARE) ? "Orbit number" : "",
            ylabel     = L"\nu_{\rm target}\ (^\circ)",
            title      = latexstring("N=$n,\\ \\kappa=$(round(κ,digits=2)),\\ ") * di_label,
            titlefontsize = 10,
            xlims      = (0, 80),
            ylims      = (0, 360),
            yticks     = 0:60:360,
            ms         = 1.5, alpha = 0.6,
            markerstrokewidth = 0,
            color      = TS_COLORS[n],
            legend     = false,
            grid       = true,
            xtickfontsize=10, ytickfontsize=10,
            xguidefontsize=11, yguidefontsize=11,
            bottom_margin = 2Plots.mm,
        )
        # Draw horizontal lines at predicted resonant phases for resonant N
        k_approx = round(Int, κ)
        if abs(κ - k_approx) < 0.15
            phase0 = phases[1]
            for ki in 0:(k_approx-1)
                predicted = mod(phase0 + ki * 360.0/k_approx, 360.0)
                hline!(sp, [predicted]; color=:black, lw=0.8, linestyle=:dash,
                       alpha=0.5, label=false)
            end
        end
        push!(panels_fig1, sp)
    end

    fig1 = plot(panels_fig1...;
        layout = grid(length(panels_fig1), 1),
        size   = (820, 280 * length(panels_fig1) + 40),
        left_margin  = 12Plots.mm, right_margin = 4Plots.mm,
        top_margin   = 8Plots.mm,  bottom_margin = 2Plots.mm,
    )
    out1 = joinpath(SAVE_DIR, "phase_scatter_vs_orbit_$(di_tag).png")
    savefig(fig1, out1)
    println("Saved → $out1")

    # ── Figure 2: Polar phase histograms ──────────────────────────────────────
    n_bins    = 36
    panels_fig2 = []
    for n in N_COMPARE
        haskey(all_data, n) || continue
        events = all_data[n]
        κ      = n / BASE

        phases    = [e.true_anom_deg for e in events]
        bin_edges = range(0.0, 360.0; length=n_bins+1)
        counts    = zeros(n_bins)
        for ν in phases
            bi = clamp(searchsortedlast(bin_edges, ν), 1, n_bins)
            counts[bi] += 1
        end
        θ_centres = [deg2rad((bin_edges[i] + bin_edges[i+1])/2) for i in 1:n_bins]
        r_max     = maximum(counts) * 1.2

        sp = plot(θ_centres, counts;
            proj       = :polar,
            seriestype = :bar,
            color      = TS_COLORS[n], alpha = 0.7,
            linewidth  = 0.3,
            title      = latexstring("N=$n\\ (\\kappa=$(round(κ,digits=2)))"),
            titlefontsize = 10,
            legend     = false,
            rlims      = (0, r_max),
        )
        push!(panels_fig2, sp)
    end

    fig2 = plot(panels_fig2...;
        layout = (1, length(panels_fig2)),
        size   = (300 * length(panels_fig2), 340),
        left_margin  = 4Plots.mm, right_margin  = 4Plots.mm,
        top_margin   = 10Plots.mm, bottom_margin = 4Plots.mm,
    )
    out2 = joinpath(SAVE_DIR, "phase_polar_histogram_$(di_tag).png")
    savefig(fig2, out2)
    println("Saved → $out2")

    # ── Figure 3: Per-activation ΔvT vs orbital phase ─────────────────────────
    panels_fig3 = []
    for (idx, n) in enumerate(N_COMPARE)
        haskey(all_data, n) || continue
        events  = all_data[n]
        κ       = n / BASE
        phases  = [e.true_anom_deg  for e in events]
        dv_incs = [e.dv_t_increment for e in events]
        orbits  = [e.orbit_start    for e in events]

        sp = scatter(phases, dv_incs;
            xlabel     = L"\nu_{\rm target}\ \mathrm{at\ activation\ centre}\ (^\circ)",
            ylabel     = idx == 1 ? L"\Delta v_T\ \mathrm{increment\ (m/s)}" : "",
            title      = latexstring("N=$n,\\ \\kappa=$(round(κ,digits=2)),\\ ") * di_label,
            titlefontsize = 10,
            xlims      = (0, 360), xticks = 0:60:360,
            ms         = 3, alpha = 0.5,
            markerstrokewidth = 0,
            marker_z   = orbits,
            color      = :plasma,
            colorbar   = idx == length(N_COMPARE),
            colorbar_title = "orbit",
            legend     = false,
            grid       = true,
            xtickfontsize=10, ytickfontsize=10,
            xguidefontsize=11, yguidefontsize=11,
            bottom_margin = 4Plots.mm,
        )
        hline!(sp, [0.0]; color=:black, lw=0.8, linestyle=:dash, label=false)
        push!(panels_fig3, sp)
    end

    fig3 = plot(panels_fig3...;
        layout = (1, length(panels_fig3)),
        size   = (320 * length(panels_fig3), 380),
        left_margin  = 12Plots.mm, right_margin  = 14Plots.mm,
        top_margin   = 8Plots.mm,  bottom_margin = 8Plots.mm,
    )
    out3 = joinpath(SAVE_DIR, "phase_dvt_per_activation_$(di_tag).png")
    savefig(fig3, out3)
    println("Saved → $out3")

    # ── Console report ─────────────────────────────────────────────────────────
    println("\n── Phase clustering report (Δi=$(delta_i)°) ──")
    println(@sprintf("%-6s  %-8s  %-12s  %-14s  %-14s",
                     "N", "κ", "n_activations", "phase_spread°", "sum_dv_t m/s"))
    println("-"^60)
    for n in N_COMPARE
        haskey(all_data, n) || continue
        events = all_data[n]
        κ      = n / BASE
        phases = [e.true_anom_deg for e in events]
        spread = length(phases) > 1 ? std(phases) : 0.0
        sum_dv = sum(e.dv_t_increment for e in events)
        println(@sprintf("%-6d  %-8.2f  %-12d  %-14.1f  %-14.5f",
                         n, κ, length(events), spread, sum_dv))
    end
end  # delta_i loop

println("\nDone. Interpretation:")
println("  phase_spread° — std of activation phases (lower = more clustered).")
println("  Resonant N: smaller spread; non-resonant: larger spread.")
