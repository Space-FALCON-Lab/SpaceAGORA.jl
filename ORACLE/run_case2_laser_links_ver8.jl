#!/usr/bin/env julia

# 1. load common.jl
const REPO_ROOT = normpath(joinpath(@__DIR__, "..")) # find path to the repository root
include(joinpath(REPO_ROOT, "examples", "common.jl")) # load common.jl for utility functions and types

# 2. load dependencies
using Arrow
using CSV
using DataFrames
using LinearAlgebra
using Printf
using DiffEqBase
using StaticArrays
using .SimulationModel  # NOTE: must come AFTER include(common.jl) — SpaceAGORA defines this submodule

# 3. -animate flag & import "10_Animation_ver2.jl"
const _HAS_GLMAKIE = "--animate" in ARGS && (try; @eval using GLMakie; true; catch; false; end)
_HAS_GLMAKIE && include(joinpath(@__DIR__, "10_Animation_ver2.jl"))

# 4. define output path
const DEFAULT_OUTPUT_DIR = joinpath(REPO_ROOT, "output")

# 5. define paper grid parameters (overridden by --paper-grid flag)
const PAPER_HELPER_ALTITUDES_KM       = (1150.0, 1050.0, 1000.0, 950.0, 850.0)
const PAPER_HELPER_INCLINATION_DELTAS_DEG = (0.0, 0.5, 1.0)
const PAPER_HELPER_COUNTS             = (200, 250, 300) #(50, 100) #(150, 200, 250, 300) #(1, 2, 3)

# 6. settings container & default values
Base.@kwdef struct OracleCase2Options
    helpers::Int = 200
    helper_altitude_km::Float64 = 1050.0
    target_altitude_km::Float64 = 1000.0
    target_inclination_deg::Float64 = 0.0
    helper_inclination_deg::Float64 = 0.0
    orbits::Float64 = 80.0
    schedule::Symbol = :naive_next_entering
    laser_range_km::Float64 = 200.0
    laser_power_w::Float64 = 10_000.0
    magnification::Float64 = 100.0
    beta::Float64 = 1.0
    eta::Float64 = 1.0
    mass_kg::Float64 = 227.0
    dt_max_s::Float64 = 10.0
    paper_grid::Bool = false
    feather_only::Bool = false  # skip CSV, summary.csv, and plots — write feather+manifest only
    output_dir::String = DEFAULT_OUTPUT_DIR
    timeseries_points::Int = 1001
    animate::Bool = false
end

# --------- Functions ---------
#include("functions/0_Module_Setup.jl")
include("functions/1_LOS_Metrics.jl")
include("functions/2_Laser_Forces_ver2.jl")
include("functions/3_Dynamics.jl")
include("functions/4_Diagnostics.jl")
include("functions/5_OE_Converters.jl")
include("functions/6_OE_and_dv_in_RTN.jl")
include("functions/7_Plots.jl")
include("functions/8_LoS_time_series.jl")
include("functions/9_Runners.jl")
#include("functions/10_Animation_ver2.jl")  # loaded conditionally above (requires --animate)

# 6. help text printed when you run the script with --help.
function _usage()
    return """
    Usage:
      julia --project=. ORACLE/run_case2_laser_links.jl [options]

    Options:
      --helpers N
      --helper-altitude-km KM
      --target-altitude-km KM
      --target-inclination-deg DEG
      --helper-inclination-deg DEG
      --orbits N
      --schedule naive_next_entering|positive_along_track
      --laser-range-km KM
      --laser-power-w W
      --magnification B
      --beta VALUE
      --eta VALUE
      --mass-kg KG
      --dt-max-s SEC
      --paper-grid
      --feather-only           Skip CSV, summary.csv, and plots; write feather+manifest only (faster I/O, less disk)
      --output-dir PATH        Directory for per-scenario output files (default: output/oracle_case2_laser_links/)
      --timeseries-points N
      --animate             Show 3D animation after the simulation (requires GLMakie)
    """
end

# 7. group the command-line option names by types
const _FLAG_OPTS   = (:paper_grid, :feather_only, :animate)
const _INT_OPTS    = (:helpers, :timeseries_points)
const _SYMBOL_OPTS = (:schedule,)
const _PATH_OPTS   = (:output_dir,)
const _FLOAT_OPTS  = (
    :helper_altitude_km, :target_altitude_km, :target_inclination_deg,
    :helper_inclination_deg, :orbits,
    :laser_range_km, :laser_power_w, :magnification, :beta, :eta,
    :mass_kg, :dt_max_s,
)

# 8. Reads what the user typed on the command line and turns it into an OracleCase2Options struct.
function _parse_options(argv)::OracleCase2Options
    opts = Dict{Symbol, Any}()
    i = 1
    while i <= length(argv) # example of argv: ["--helpers", "200", "--laser-range-km", "150.0", "--schedule", "positive_along_track"]
        arg = argv[i] # Get the current argument string, e.g. "--helpers" or "--laser-range-km"
        arg in ("--help", "-h") && (println(_usage()); exit(0)) # If the user typed --help or -h, print the usage text and quit
        startswith(arg, "--") || throw(ArgumentError("Unexpected argument '$arg'.\n$(_usage())")) # Every argument must start with "--", otherwise crash with an error
        key = Symbol(replace(arg[3:end], '-' => '_')) # Strip the leading "--" and convert dashes to underscores
        if key in _FLAG_OPTS
            opts[key] = true; i += 1; continue # For flags like --animate where there is no value after them, just mark them as true, step forward by 1, and go to the next argument
        end
        i < length(argv) || throw(ArgumentError("Missing value for $arg.")) # Make sure there IS a next input to read as the value for this arg
        val = argv[i + 1] # Grab the next word as the value, e.g. "200" for --helpers 200
        if     key in _INT_OPTS    opts[key] = parse(Int, val) # if the key is in the integer options, parse the value as an integer
        elseif key in _SYMBOL_OPTS opts[key] = Symbol(val) # if the key is in the symbol options, convert the value to a symbol
        elseif key in _PATH_OPTS   opts[key] = abspath(val) # if the key is in the path options, convert the value to an absolute path
        elseif key in _FLOAT_OPTS  opts[key] = parse(Float64, val) # if the key is in the float options, parse the value as a float
        else   throw(ArgumentError("Unknown option $arg.\n$(_usage())"))
        end
        i += 2 # Step forward by 2 to move past the key and its value
    end
    return OracleCase2Options(; opts...)
end

# 9. Validate that the inputs are within reasonable bounds, and throw an error if not.
function _validate_options(opts::OracleCase2Options)
    opts.helpers >= 1 || throw(ArgumentError("--helpers must be >= 1."))
    opts.helper_altitude_km > 0.0 || throw(ArgumentError("--helper-altitude-km must be positive."))
    opts.target_altitude_km > 0.0 || throw(ArgumentError("--target-altitude-km must be positive."))
    opts.orbits > 0.0 || throw(ArgumentError("--orbits must be positive."))
    opts.schedule in (:naive_next_entering, :positive_along_track) ||
        throw(ArgumentError("--schedule must be naive_next_entering or positive_along_track."))
    opts.laser_range_km >= 0.0 || throw(ArgumentError("--laser-range-km must be nonnegative."))
    opts.laser_power_w >= 0.0 || throw(ArgumentError("--laser-power-w must be nonnegative."))
    opts.magnification >= 0.0 || throw(ArgumentError("--magnification must be nonnegative."))
    opts.beta >= 0.0 || throw(ArgumentError("--beta must be nonnegative."))
    opts.eta >= 0.0 || throw(ArgumentError("--eta must be nonnegative."))
    opts.mass_kg > 0.0 || throw(ArgumentError("--mass-kg must be positive."))
    opts.dt_max_s > 0.0 || throw(ArgumentError("--dt-max-s must be positive."))
    opts.timeseries_points >= 2 || throw(ArgumentError("--timeseries-points must be >= 2."))
    return nothing
end

# 10. Copy opts, assist publication grid search
_with(opts::OracleCase2Options; kwargs...) =
    OracleCase2Options(; (f => getfield(opts, f) for f in fieldnames(OracleCase2Options))..., kwargs...)

function main(argv=ARGS)
    opts = _parse_options(argv)  # parse command-line arguments into an OracleCase2Options struct

    if opts.paper_grid
        # --- Paper-grid mode: sweep over all (helper_altitude, inclination_delta) combinations ---

        for n_helpers in PAPER_HELPER_COUNTS                        # outermost loop: helper counts
            for helper_alt in PAPER_HELPER_ALTITUDES_KM               # middle loop: 5 altitudes
                for helper_inclination_delta in PAPER_HELPER_INCLINATION_DELTAS_DEG  # inner loop: 3 inclination deltas
                    case_opts = _with(opts;
                        helpers=n_helpers,
                        helper_altitude_km=helper_alt,
                        helper_inclination_deg=opts.target_inclination_deg + helper_inclination_delta,
                        output_dir=joinpath(opts.output_dir, "paper_plot_mode"),
                        animate=false,
                    )
                    elapsed = @elapsed begin
                        result = run_open_cavity_case_native(case_opts)
                    end
                    s = result.summary
                    println("  → $(result.results_dir)")
                    @printf(
                        "helpers=%d helper_alt_km=%.1f helper_inc_deg=%.1f dv_R=%.6e dv_T=%.6e dv_N=%.6e activations=%d  [%.1f s]\n",
                        s.helpers, s.helper_altitude_km, s.helper_inclination_deg,
                        s.dv_r_mps, s.dv_t_mps, s.dv_n_mps, s.activations, elapsed
                    )
                end
            end
        end
        println("Output directory: $(opts.output_dir)")

    else
        # --- Single-run mode ---
        elapsed = @elapsed begin
            result = run_open_cavity_case_native(_with(opts; output_dir=joinpath(opts.output_dir, "single_case_mode")))
        end
        s = result.summary
        _print_summary(s)
        @printf("  run time: %.1f s\n", elapsed)
        println("Output directory: $(result.results_dir)")
        img_dir = joinpath(result.results_dir, "images")

        # --- Diagnostic plots (skipped when --feather-only) ---
        opts.feather_only || plot_open_cavity_results(result, opts; IMG_DIR=img_dir, target_only=false)

        # --- Optional 3D animation (only if --animate was passed and GLMakie loaded) ---
        if opts.animate
            if _HAS_GLMAKIE
                p_anim = Dict{Symbol, Any}(:N => opts.helpers + 1)  # total spacecraft count
                _, anim_controls = animate_all_satellites_3d_smooth_helper_target(
                    result.sol, p_anim, opts.helpers)
                println("Animation window opened. Close the window to exit.")
                wait(anim_controls.screen)  # block until the user closes the window
            else
                @warn "GLMakie is not installed — animation skipped. " *
                      "Install it with: using Pkg; Pkg.add(\"GLMakie\")"
            end
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
