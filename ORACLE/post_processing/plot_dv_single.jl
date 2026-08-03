#!/usr/bin/env julia
# Plots one ΔV component vs N_helpers for a single target altitude,
# with one curve per inclination delta.
#
# Run from the repository root:
#   julia --project=. ORACLE/post_processing/plot_dv_single.jl [options]
#
# Options:
#   --schedule         gve_sma          (required)
#   --target-ecc       0.0
#   --target-nu-deg    0.0
#   --target-altitude-km  850.0
#   --component        T                (R | T | N)

using Arrow
get!(ENV, "GKSwstype", "100")
using Plots
using LaTeXStrings
using Printf

gr()

const REPO_ROOT  = normpath(joinpath(@__DIR__, "..", ".."))
const OUTPUT_DIR = joinpath(REPO_ROOT, "output", "paper_plot_mode")

# ── CLI ───────────────────────────────────────────────────────────────────────
function _get(args, flag, default)
    idx = findfirst(==(flag), args)
    idx === nothing ? default : args[idx + 1]
end

const SCHEDULE       = _get(ARGS, "--schedule",          "naive_next_entering")
const TARGET_ECC     = parse(Float64, _get(ARGS, "--target-ecc",       "0.0"))
const TARGET_NU_DEG  = parse(Float64, _get(ARGS, "--target-nu-deg",    "0.0"))
const TARGET_ALT_KM  = parse(Float64, _get(ARGS, "--target-altitude-km", "850.0"))
const COMPONENT      = uppercase(_get(ARGS, "--component", "T"))

COMPONENT in ("R", "T", "N") || error("--component must be R, T, or N")

# ── Grid parameters ───────────────────────────────────────────────────────────
const HELPER_ALT_KM          = 1000.0
const HELPER_INCLINATION_DEG = 0.0
const INCLINATION_DELTAS_DEG = [0.0, 0.5, 1.0]
const HELPER_COUNTS          = [1, 50, 100, 150, 200, 250, 300]

# ── Feather reader ────────────────────────────────────────────────────────────
function read_final_dv(output_dir, helper_alt_km, target_alt_km,
                        helper_inc_deg, target_inc_deg, n_helpers;
                        schedule::String, target_ecc::Float64, target_nu_deg::Float64)
    alt_folder = @sprintf("h%.0fkm_t%.0fkm",    helper_alt_km,  target_alt_km)
    inc_folder = @sprintf("ih%.1fdeg_it%.1fdeg", helper_inc_deg, target_inc_deg)
    n_folder   = @sprintf("N%d",                 n_helpers)
    base_path  = joinpath(output_dir, alt_folder, inc_folder, n_folder)
    isdir(base_path) || return nothing

    t_folders = filter(d -> startswith(d, "T") && endswith(d, "s"), readdir(base_path))
    isempty(t_folders) && return nothing

    sched_folder = @sprintf("%s_e%.4f_nu%.4f", schedule, target_ecc, target_nu_deg)
    feather_path = joinpath(base_path, first(t_folders), sched_folder, "simulation_results.feather")
    isfile(feather_path) || return nothing

    tbl = Arrow.Table(feather_path)
    n   = length(tbl.time)
    return (
        dv_r = Float64(tbl.dv_r_accumulated[n]),
        dv_t = Float64(tbl.dv_t_accumulated[n]),
        dv_n = Float64(tbl.dv_n_accumulated[n]),
    )
end

# ── Collect data ──────────────────────────────────────────────────────────────
inc_colors = [:blue, :red, :green]
inc_labels = [latexstring("\\Delta i=$(d)^{\\circ}") for d in INCLINATION_DELTAS_DEG]

all_xs  = Vector{Vector{Float64}}()
all_dvs = Vector{Vector{Float64}}()

for delta_i in INCLINATION_DELTAS_DEG
    t_inc = HELPER_INCLINATION_DEG + delta_i
    xs  = Float64[]
    dvs = Float64[]
    for n in HELPER_COUNTS
        res = read_final_dv(OUTPUT_DIR, HELPER_ALT_KM, TARGET_ALT_KM,
                            HELPER_INCLINATION_DEG, t_inc, n;
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

all_vals = vcat(all_dvs...)
isempty(all_vals) && error("No data found for the requested parameters.")

# ── Plot ──────────────────────────────────────────────────────────────────────
comp_label = COMPONENT == "R" ? L"\Delta v_R\ \mathrm{(m/s)}" :
             COMPONENT == "T" ? L"\Delta v_T\ \mathrm{(m/s)}" :
                                L"\Delta v_N\ \mathrm{(m/s)}"

title_str = latexstring(
    "$(replace(SCHEDULE, "_" => "\\_"))\\ \\ h_{\\mathrm{t}}=$(round(Int, TARGET_ALT_KM))\\,\\mathrm{km}\\ \\ " *
    "e=$(TARGET_ECC)\\ \\ \\nu=$(TARGET_NU_DEG)^{\\circ}"
)

ref_xs = first(filter(!isempty, all_xs))

fig = plot(
    title=title_str, titlefontsize=11,
    xlabel=L"N_{\mathrm{helpers}}", xguidefontsize=12,
    ylabel=comp_label,              yguidefontsize=12,
    xticks=(ref_xs, string.(round.(Int, ref_xs))),
    xrotation=45,
    legend=:outertopright,
    grid=true, tick_direction=:out,
    size=(520, 380),
    left_margin=6Plots.mm, bottom_margin=6Plots.mm,
)

for (k, (xs, dvs)) in enumerate(zip(all_xs, all_dvs))
    isempty(xs) && continue
    plot!(fig, xs, dvs;
          color=inc_colors[k], lw=1.5, linestyle=:solid,
          marker=:circle, ms=4, markerstrokewidth=0.2,
          label=inc_labels[k])
end

# ── Save ──────────────────────────────────────────────────────────────────────
out_dir = joinpath(REPO_ROOT, "output", "images", "delta_v")
mkpath(out_dir)

stem = @sprintf("dv_%s_%s_e%.4f_nu%.4f_t%.0fkm",
    COMPONENT, SCHEDULE, TARGET_ECC, TARGET_NU_DEG, TARGET_ALT_KM)
savefig(fig, joinpath(out_dir, stem * ".png"))
savefig(fig, joinpath(out_dir, stem * ".pdf"))
display(fig)
println("Saved → ", joinpath(out_dir, stem * ".pdf"))
