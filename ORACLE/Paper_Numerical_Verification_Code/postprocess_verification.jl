#!/usr/bin/env julia
# =============================================================================
# Section IV.F — Post-Processing Driver
# =============================================================================
# Reads the Feather files produced by run_verification.jl, computes all
# conservation-law residuals (Sec. IV.F, Eqs. 21–23), and produces the
# validation plots equivalent to Figures 2 and 3 in the paper.
#
# Usage (after running run_verification.jl):
#
#   # Default — expects both j2_off/ and j2_on/ under the default output dir:
#   julia --project=. ORACLE/Paper_Numerical_Verification_Code/postprocess_verification.jl
#
#   # Custom output directory:
#   julia --project=. ORACLE/Paper_Numerical_Verification_Code/postprocess_verification.jl \
#         --output-dir /path/to/Paper_Numerical_Verification_Code
#
#   # Single-run mode (only J2 on or only J2 off was simulated):
#   julia --project=. ORACLE/Paper_Numerical_Verification_Code/postprocess_verification.jl \
#         --j2-on-dir  path/to/j2_on  --j2-off-dir path/to/j2_off
# =============================================================================

# ── 1. Bootstrap SpaceAGORA (loads common.jl → activates project) ────────────
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(REPO_ROOT, "examples", "common.jl"))

# ── 2. Dependencies ──────────────────────────────────────────────────────────
using Arrow
using DataFrames
using Plots
using Printf
using LinearAlgebra
using StaticArrays

# ── 3. Load compute and plot layers ──────────────────────────────────────────
include(joinpath(@__DIR__, "verification_compute.jl"))
include(joinpath(@__DIR__, "verification_plots.jl"))

# ── 4. Default paths ─────────────────────────────────────────────────────────
const DEFAULT_BASE_DIR = joinpath(REPO_ROOT, "output", "Paper_Numerical_Verification_Code")

# ─────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────────────────────
function _pp_usage()
    return """
    Usage:
      julia --project=. ORACLE/Paper_Numerical_Verification_Code/postprocess_verification.jl [options]

    Options:
      --output-dir PATH     Base directory that contains j2_off/ and j2_on/
                            subdirectories written by run_verification.jl.
                            [default: output/Paper_Numerical_Verification_Code]
      --j2-off-dir PATH     Explicit path to the J2-inactive results directory.
      --j2-on-dir  PATH     Explicit path to the J2-active  results directory.
      --no-j2               Skip the J2-inactive case (plot only J2-active).
      --j2-only             Skip the J2-active  case (plot only J2-inactive).
      --laser-range-km KM   Laser range for force recomputation [default: 200.0]
      --laser-power-w W     Laser power   [default: 10000.0]
      --magnification B     OCL magnification [default: 100.0]
      --beta VALUE          Geometric loss [default: 1.0]
      --eta VALUE           Efficiency    [default: 1.0]
      --mass-kg KG          Spacecraft mass [default: 227.0]
    """
end

function _pp_parse_options(argv)
    opts = Dict{Symbol,Any}(
        :output_dir        => DEFAULT_BASE_DIR,
        :j2_off_dir        => nothing,
        :j2_on_dir         => nothing,
        :skip_no_j2        => false,
        :skip_with_j2      => false,
        :laser_range_km    => 200.0,
        :laser_power_w     => 10_000.0,
        :magnification     => 100.0,
        :beta              => 1.0,
        :eta               => 1.0,
        :mass_kg           => 227.0,
    )
    float_keys = (:laser_range_km, :laser_power_w, :magnification, :beta, :eta, :mass_kg)
    i = 1
    while i <= length(argv)
        arg = argv[i]
        arg in ("--help", "-h") && (println(_pp_usage()); exit(0))
        if arg == "--no-j2";   opts[:skip_no_j2]   = true; i += 1; continue; end
        if arg == "--j2-only"; opts[:skip_with_j2] = true; i += 1; continue; end
        startswith(arg, "--") || throw(ArgumentError("Unexpected argument '$arg'."))
        key = Symbol(replace(arg[3:end], '-' => '_'))
        i < length(argv) || throw(ArgumentError("Missing value for $arg."))
        val = argv[i + 1]
        if key in float_keys
            opts[key] = parse(Float64, val)
        elseif key in (:output_dir, :j2_off_dir, :j2_on_dir)
            opts[key] = abspath(val)
        else
            throw(ArgumentError("Unknown option: $arg"))
        end
        i += 2
    end
    return opts
end

# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────
function main(argv=ARGS)
    opts = _pp_parse_options(argv)

    # ── Resolve result directories ────────────────────────────────────────────
    base_dir   = opts[:output_dir]
    j2_off_dir = something(opts[:j2_off_dir], joinpath(base_dir, "j2_off"))
    j2_on_dir  = something(opts[:j2_on_dir],  joinpath(base_dir, "j2_on"))

    # ── Plots output directory ────────────────────────────────────────────────
    plots_dir = joinpath(base_dir, "plots")
    mkpath(plots_dir)

    # ── Laser parameters (for OCL force recomputation) ───────────────────────
    lp = LaserParams(
        range_km      = opts[:laser_range_km],
        power_w       = opts[:laser_power_w],
        magnification = opts[:magnification],
        beta          = opts[:beta],
        eta           = opts[:eta],
    )
    mass_kg = opts[:mass_kg]

    println("\n=== Section IV.F Conservation-Law Validation — Post-Processing ===\n")

    # ── Load and compute diagnostics ──────────────────────────────────────────
    diag_no_j2   = nothing
    diag_with_j2 = nothing

    if !opts[:skip_no_j2]
        if isdir(j2_off_dir) && isfile(joinpath(j2_off_dir, "simulation_results.feather"))
            diag_no_j2 = load_and_compute(j2_off_dir, lp; mass_kg=mass_kg)
            print_diagnostics_summary(diag_no_j2; label="J2 inactive")
        else
            @warn "J2-inactive results not found at: $j2_off_dir\n" *
                  "Run run_verification.jl --no-j2 first, or pass --j2-only to skip."
        end
    end

    if !opts[:skip_with_j2]
        if isdir(j2_on_dir) && isfile(joinpath(j2_on_dir, "simulation_results.feather"))
            diag_with_j2 = load_and_compute(j2_on_dir, lp; mass_kg=mass_kg)
            print_diagnostics_summary(diag_with_j2; label="J2 active")
        else
            @warn "J2-active results not found at: $j2_on_dir\n" *
                  "Run run_verification.jl --j2-only first, or pass --no-j2 to skip."
        end
    end

    # ── Select plotting backend ───────────────────────────────────────────────
    pyplot_available = try; @eval using PyPlot; true; catch; false; end
    if pyplot_available
        pyplot()
    else
        try; @eval using GR; gr(); catch; end
    end

    # ── Produce plots ─────────────────────────────────────────────────────────
    if diag_no_j2 !== nothing && diag_with_j2 !== nothing
        # Side-by-side panels (mirrors paper Figure 2 & 3)
        println("Generating paired Figure-2 style plot (force / torque / angular momentum)...")
        plot_figure2_style(diag_no_j2, diag_with_j2; save_dir=plots_dir)

        println("Generating paired Figure-3 style plot (energy / work / residual)...")
        plot_figure3_style(diag_no_j2, diag_with_j2; save_dir=plots_dir)

        println("Generating per-satellite ΔH_z plot (J2 inactive)...")
        plot_per_satellite_Hz_change(diag_no_j2; save_dir=plots_dir, label="J2 inactive")

        println("Generating per-satellite ΔH_z plot (J2 active)...")
        plot_per_satellite_Hz_change(diag_with_j2; save_dir=plots_dir, label="J2 active")

    elseif diag_no_j2 !== nothing
        println("Generating single-case Figure-2 style plot (J2 inactive)...")
        plot_single_case_figure2(diag_no_j2; save_dir=plots_dir, label="J2 inactive")

        println("Generating single-case Figure-3 style plot (J2 inactive)...")
        plot_single_case_figure3(diag_no_j2; save_dir=plots_dir, label="J2 inactive")

        println("Generating per-satellite ΔH_z plot (J2 inactive)...")
        plot_per_satellite_Hz_change(diag_no_j2; save_dir=plots_dir, label="J2 inactive")

    elseif diag_with_j2 !== nothing
        println("Generating single-case Figure-2 style plot (J2 active)...")
        plot_single_case_figure2(diag_with_j2; save_dir=plots_dir, label="J2 active")

        println("Generating single-case Figure-3 style plot (J2 active)...")
        plot_single_case_figure3(diag_with_j2; save_dir=plots_dir, label="J2 active")

        println("Generating per-satellite ΔH_z plot (J2 active)...")
        plot_per_satellite_Hz_change(diag_with_j2; save_dir=plots_dir, label="J2 active")

    else
        println("No valid result directories found. Run run_verification.jl first.")
        return nothing
    end

    println("\nPlots saved to: $plots_dir")

    # ── Print final validation verdict ────────────────────────────────────────
    println("\n── Validation Verdict ──────────────────────────────────────────────")
    _verdict(diag, tag) = begin
        diag === nothing && return
        peak_Fnet    = maximum(diag.F_net_mag)
        peak_tau     = maximum(diag.tau_net_mag)
        peak_eps_E   = maximum(abs.(diag.eps_E))
        peak_eps_H   = maximum(diag.eps_H_mag)
        W_ocl_final  = diag.W_ocl[end]
        dE_final     = diag.delta_Eorb[end]
        println("$tag:")
        @printf("  |F_net^OCL| ≤ %.2e N  (machine precision = OK if < 1e-10 N)\n", peak_Fnet)
        @printf("  |τ_net^OCL| ≤ %.2e N·m\n",                                    peak_tau)
        @printf("  |ε_E| ≤ %.2e J  (should be << |W_OCL| = %.2e J)\n", peak_eps_E, abs(W_ocl_final))
        @printf("  |ε_H| ≤ %.2e kg·m²/s\n",                                      peak_eps_H)
        rel_eps_E = abs(W_ocl_final) > 0 ? peak_eps_E / abs(W_ocl_final) : NaN
        @printf("  Relative energy residual |ε_E|/|W_OCL| = %.2e\n", rel_eps_E)
        println()
    end
    _verdict(diag_no_j2,   "J2 inactive")
    _verdict(diag_with_j2, "J2 active  ")

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
