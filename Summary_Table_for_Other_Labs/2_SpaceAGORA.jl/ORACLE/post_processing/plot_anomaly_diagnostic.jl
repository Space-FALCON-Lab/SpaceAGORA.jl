#!/usr/bin/env julia
# Diagnostic 3-panel figure for the N≈200 outlier observed at h_helper=1150 km.
#
# Panel 1 (top):    ΔvT vs N — one curve per inclination delta, anomaly band shaded
# Panel 2 (middle): activations vs N — same layout
# Panel 3 (bottom): normalised active-helper index vs time (orbits) for
#                   N=150, 200, 250 overlaid — reveals when/why they diverge
#
# Run from the repository root:
#   julia --project=. ORACLE/post_processing/plot_anomaly_diagnostic.jl
#
# Optional CLI overrides (positional):
#   julia --project=. ORACLE/post_processing/plot_anomaly_diagnostic.jl \
#         1150.0   # helper_alt_km  (default 1150.0)

using Arrow
using DataFrames
get!(ENV, "GKSwstype", "100")
using Plots
using LaTeXStrings
using Printf

gr()

# ── Paths ─────────────────────────────────────────────────────────────────────
const REPO_ROOT  = normpath(joinpath(@__DIR__, "..", ".."))
const OUTPUT_DIR = joinpath(REPO_ROOT, "output", "paper_plot_mode")
const SAVE_DIR   = joinpath(REPO_ROOT, "output", "images")
mkpath(SAVE_DIR)

# ── Configuration ─────────────────────────────────────────────────────────────
const HELPER_ALT_KM          = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1150.0
const TARGET_ALT_KM          = 1000.0
const TARGET_INCLINATION_DEG = 0.0
const INCLINATION_DELTAS_DEG = [0.0, 0.5, 1.0]

# All N values to include in the sweep panels (panels 1 & 2).
# The script auto-skips any N whose feather file is absent.
const N_SWEEP = [1, 2, 3, 50, 100, 150, 197, 198, 199, 200, 201, 202, 250, 299, 300]

# N values shown in the time-series panel (panel 3). One curve each.
const N_TIMESERIES = [150, 200, 250]

# Target orbital period [s] — used to convert time axis to orbits.
# Approximate: T = 2π √(a³/μ),  a = R_earth + target_alt,  μ = 3.986e14
const _MU        = 3.986004418e14
const _R_EARTH   = 6_378_137.0
const TARGET_PERIOD_S = 2π * sqrt((_R_EARTH + TARGET_ALT_KM * 1e3)^3 / _MU)

# Anomaly band drawn on panels 1 & 2
const ANOMALY_LO = 197
const ANOMALY_HI = 202

# ── Resonance prediction ──────────────────────────────────────────────────────
# N*(k) = k × T_h / (T_h - T_t)  — integer crossings per target orbit.
# Returns all resonant N values in [N_lo, N_hi].
function resonant_N_values(helper_alt_km, target_alt_km; N_lo=1, N_hi=400)
    T(h) = 2π * sqrt((_R_EARTH + h*1e3)^3 / _MU)
    Tt = T(target_alt_km); Th = T(helper_alt_km)
    abs(Th - Tt) < 1.0 && return Float64[]   # same orbit → no drift
    base = Th / abs(Th - Tt)
    return [k * base for k in 1:10_000 if N_lo <= k * base <= N_hi]
end

const RESONANT_N = resonant_N_values(HELPER_ALT_KM, TARGET_ALT_KM; N_lo=1,
                                     N_hi=maximum(N_SWEEP)*1.1)

# ── Colours & labels ──────────────────────────────────────────────────────────
const INC_COLORS  = [:blue, :red, :green]
const INC_LABELS  = [latexstring("\\Delta i=$(round(d, digits=1))^\\circ")
                     for d in INCLINATION_DELTAS_DEG]
const TS_COLORS   = [:dodgerblue, :crimson, :darkorange]

# ── Directory helper ──────────────────────────────────────────────────────────
function _feather_path(helper_alt_km, target_alt_km, helper_inc_deg, target_inc_deg, n)
    base = joinpath(OUTPUT_DIR,
        @sprintf("h%.0fkm_t%.0fkm",    helper_alt_km,  target_alt_km),
        @sprintf("ih%.1fdeg_it%.1fdeg", helper_inc_deg, target_inc_deg),
        @sprintf("N%d", n))
    isdir(base) || return nothing
    t_dirs = filter(d -> startswith(d, "T") && endswith(d, "s"), readdir(base))
    isempty(t_dirs) && return nothing
    p = joinpath(base, first(t_dirs), "simulation_results.feather")
    isfile(p) ? p : nothing
end

# ── Per-scenario scalar extraction ───────────────────────────────────────────
# Returns (dv_t, approx_activations) or nothing.
function read_scalars(helper_alt_km, target_alt_km, helper_inc_deg, target_inc_deg, n)
    p = _feather_path(helper_alt_km, target_alt_km, helper_inc_deg, target_inc_deg, n)
    p === nothing && return nothing
    tbl  = Arrow.Table(p)
    nrow = length(tbl.time)
    dv_t = Float64(tbl.dv_t_accumulated[nrow])
    # Activations: count transitions in laser_active_helper (piecewise-constant at
    # the 1001 saved points — a lower bound on the true ODE-step count, but
    # proportional across scenarios so it's a valid relative metric).
    helper_series = tbl.laser_active_helper
    activations = count(i -> helper_series[i] != helper_series[i-1], 2:nrow)
    return (dv_t=dv_t, activations=activations)
end

# ── Full time-series for the bottom panel ─────────────────────────────────────
# Returns (time_orbits, normalised_helper_index) where helper index is divided
# by N so all three N values share the same [0, 1] y-axis.
function read_timeseries(helper_alt_km, target_alt_km, helper_inc_deg, target_inc_deg, n)
    p = _feather_path(helper_alt_km, target_alt_km, helper_inc_deg, target_inc_deg, n)
    p === nothing && return nothing
    tbl            = Arrow.Table(p)
    time_orbits    = Float64.(tbl.time) ./ TARGET_PERIOD_S
    helper_norm    = Float64.(tbl.laser_active_helper) ./ n   # 0 = off, ~1 = last helper
    return (time_orbits=time_orbits, helper_norm=helper_norm)
end

# ── Build sweep data (panels 1 & 2) ──────────────────────────────────────────
println("Reading sweep data for h_helper=$(HELPER_ALT_KM) km …")

# For each inclination delta: vectors of (N, dv_t, activations) for found files
sweep = [begin
    inc_deg = TARGET_INCLINATION_DEG + δi
    ns, dvts, acts = Float64[], Float64[], Float64[]
    for n in N_SWEEP
        r = read_scalars(HELPER_ALT_KM, TARGET_ALT_KM, inc_deg, TARGET_INCLINATION_DEG, n)
        r === nothing && continue
        push!(ns,   Float64(n))
        push!(dvts, r.dv_t)
        push!(acts, Float64(r.activations))
    end
    (ns=ns, dvts=dvts, acts=acts)
end for δi in INCLINATION_DELTAS_DEG]

# ── Build time-series data (panel 3) — use Δi = 0 for clarity ────────────────
println("Reading time-series data for N = $(N_TIMESERIES) …")

ts_data = [begin
    r = read_timeseries(HELPER_ALT_KM, TARGET_ALT_KM,
                        TARGET_INCLINATION_DEG, TARGET_INCLINATION_DEG, n)
    (n=n, data=r)
end for n in N_TIMESERIES]

# ── Plotting helpers ──────────────────────────────────────────────────────────
function add_anomaly_band!(sp, y_lo, y_hi)
    vspan!(sp, [ANOMALY_LO, ANOMALY_HI];
           fillcolor=:orange, fillalpha=0.15,
           linecolor=:orange, linealpha=0.4, lw=0.5,
           label=false)
end

tick_fs = 11
label_fs = 12

# ── Panel 1: ΔvT vs N ─────────────────────────────────────────────────────────
all_dvts = vcat([s.dvts for s in sweep]...)
y_lo1 = isempty(all_dvts) ? 0.0 : minimum(all_dvts)
y_hi1 = isempty(all_dvts) ? 1.0 : maximum(all_dvts)
y_pad1 = max((y_hi1 - y_lo1) * 0.12, abs(y_hi1) * 0.05, 1e-9)

sp1 = plot(
    title  = latexstring("h_{\\rm helper}=$(round(Int, HELPER_ALT_KM))\\,{\\rm km},\\;" *
                         "h_{\\rm target}=$(round(Int, TARGET_ALT_KM))\\,{\\rm km}"),
    ylabel = L"\Delta v_T\ \mathrm{(m/s)}",
    xlabel = "",
    xlims  = (0, maximum(N_SWEEP) * 1.05),
    ylims  = (y_lo1 - y_pad1, y_hi1 + y_pad1),
    legend = :topright,
    xtickfontsize=tick_fs, ytickfontsize=tick_fs,
    xguidefontsize=label_fs, yguidefontsize=label_fs,
    tick_direction=:out, grid=true,
    bottom_margin=2Plots.mm,
)
add_anomaly_band!(sp1, y_lo1 - y_pad1, y_hi1 + y_pad1)
annotate!(sp1, (ANOMALY_LO + ANOMALY_HI)/2, y_hi1 + y_pad1 * 0.6,
          text("anomaly\nregion", :orange, :center, 8))
# Predicted resonance lines
for N_r in RESONANT_N
    vline!(sp1, [N_r]; color=:gray, lw=0.8, linestyle=:dash, alpha=0.6, label=false)
end

for (k, s) in enumerate(sweep)
    isempty(s.ns) && continue
    plot!(sp1, s.ns, s.dvts;
          color=INC_COLORS[k], lw=1.8, marker=:circle, ms=4,
          markerstrokewidth=0.2, label=INC_LABELS[k])
end

# ── Panel 2: activations vs N ─────────────────────────────────────────────────
all_acts = vcat([s.acts for s in sweep]...)
y_lo2 = isempty(all_acts) ? 0.0 : minimum(all_acts)
y_hi2 = isempty(all_acts) ? 1.0 : maximum(all_acts)
y_pad2 = max((y_hi2 - y_lo2) * 0.12, abs(y_hi2) * 0.05, 1e-9)

sp2 = plot(
    ylabel = L"\mathrm{Activations}\ (\approx)",
    xlabel = "",
    xlims  = (0, maximum(N_SWEEP) * 1.05),
    ylims  = (y_lo2 - y_pad2, y_hi2 + y_pad2),
    legend = false,
    xtickfontsize=tick_fs, ytickfontsize=tick_fs,
    xguidefontsize=label_fs, yguidefontsize=label_fs,
    tick_direction=:out, grid=true,
    bottom_margin=2Plots.mm,
)
add_anomaly_band!(sp2, y_lo2 - y_pad2, y_hi2 + y_pad2)
for N_r in RESONANT_N
    vline!(sp2, [N_r]; color=:gray, lw=0.8, linestyle=:dash, alpha=0.6, label=false)
end

for (k, s) in enumerate(sweep)
    isempty(s.ns) && continue
    plot!(sp2, s.ns, s.acts;
          color=INC_COLORS[k], lw=1.8, marker=:circle, ms=4,
          markerstrokewidth=0.2, label=false)
end

# ── Panel 3: active-helper / N vs time (orbits) ───────────────────────────────
sp3 = plot(
    ylabel = L"\mathrm{Active\ helper\ index}\ /\ N",
    xlabel = L"\mathrm{Simulation\ time\ (orbits)}",
    xlims  = (0, 80),
    ylims  = (-0.05, 1.1),
    legend = :topright,
    xtickfontsize=tick_fs, ytickfontsize=tick_fs,
    xguidefontsize=label_fs, yguidefontsize=label_fs,
    tick_direction=:out, grid=true,
    bottom_margin=4Plots.mm,
)
annotate!(sp3, 1.0, 1.05,
          text("0 = laser off,  1 = highest-index helper firing", :grey, :left, 8))

for (k, entry) in enumerate(ts_data)
    entry.data === nothing && continue
    plot!(sp3, entry.data.time_orbits, entry.data.helper_norm;
          color=TS_COLORS[k], lw=0.8, alpha=0.75,
          linetype=:steppost,
          label=latexstring("N=$(entry.n)"))
end

# ── Combine & save ────────────────────────────────────────────────────────────
fig = plot(sp1, sp2, sp3;
           layout=grid(3, 1; heights=[0.33, 0.33, 0.34]),
           size=(820, 960),
           left_margin=12Plots.mm, right_margin=4Plots.mm,
           top_margin=6Plots.mm)

out_path = joinpath(SAVE_DIR, @sprintf("anomaly_diagnostic_h%.0fkm.png",
                                       HELPER_ALT_KM))
savefig(fig, out_path)
println("Saved → $out_path")
