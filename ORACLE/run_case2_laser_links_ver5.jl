#!/usr/bin/env julia

# 1. load common.jl
const REPO_ROOT = normpath(joinpath(@__DIR__, "..")) # find path to the repository root
include(joinpath(REPO_ROOT, "examples", "common.jl")) # load common.jl for utility functions and types

# 2. load dependencies
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
const DEFAULT_SUMMARY_CSV = joinpath(REPO_ROOT, "output", "oracle_case2_laser_links", "case2_laser_summary.csv")
const DEFAULT_TIMESERIES_CSV = joinpath(REPO_ROOT, "output", "oracle_case2_laser_links", "case2_laser_timeseries.csv")

# 5. define paper grid parameters (overridden by --paper-grid flag)
const PAPER_HELPER_ALTITUDES_KM = (1150.0, 1050.0, 1000.0, 950.0, 850.0)
const PAPER_HELPER_INCLINATION_DELTAS_DEG = (0.0, 0.5, 1.0)

# 6. settings container & default values
Base.@kwdef struct OracleCase2Options
    helpers::Int = 200
    helper_altitude_km::Float64 = 1050.0
    target_altitude_km::Float64 = 1000.0
    inclination_deg::Float64 = 0.0
    helper_inclination_delta_deg::Float64 = 0.0
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
    output_csv::String = DEFAULT_SUMMARY_CSV
    timeseries_csv::String = DEFAULT_TIMESERIES_CSV
    append_output::Bool = false
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
      --inclination-deg DEG
      --helper-inclination-delta-deg DEG
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
      --output-csv PATH
      --timeseries-csv PATH
      --append-output
      --timeseries-points N
      --animate             Show 3D animation after the simulation (requires GLMakie)
    """
end

# 7. group the command-line option names by types
const _FLAG_OPTS   = (:paper_grid, :append_output, :animate)
const _INT_OPTS    = (:helpers, :timeseries_points)
const _SYMBOL_OPTS = (:schedule,)
const _PATH_OPTS   = (:output_csv, :timeseries_csv)
const _FLOAT_OPTS  = (
    :helper_altitude_km, :target_altitude_km, :inclination_deg,
    :helper_inclination_delta_deg, :orbits,
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

        # Delete old output files unless the user asked to append
        opts.append_output || (isfile(opts.output_csv)      && rm(opts.output_csv))
        opts.append_output || (isfile(opts.timeseries_csv)  && rm(opts.timeseries_csv))

        for helper_alt in PAPER_HELPER_ALTITUDES_KM               # outer loop: 5 altitudes
            for helper_inclination_delta in PAPER_HELPER_INCLINATION_DELTAS_DEG  # inner loop: 3 inclination deltas
                # Run one simulation with this (altitude, inclination_delta) combination,
                # keeping all other opts unchanged
                result = run_open_cavity_case(_with(opts;
                    helper_altitude_km=helper_alt,
                    helper_inclination_delta_deg=helper_inclination_delta,
                    append_output=true,   # always append inside the grid loop
                    animate=false,        # no animation during batch runs
                ))
                s = result.summary

                # Append this run's one-row summary and full timeseries to the CSVs
                _write_csv!(opts.output_csv,     [s];              append=true)
                _write_csv!(opts.timeseries_csv, result.timeseries; append=true)

                # Print a one-line progress update to the terminal
                @printf(
                    "helpers=%d helper_alt_km=%.1f helper_dinc_deg=%.1f dv_R=%.6e dv_T=%.6e dv_N=%.6e activations=%d\n",
                    s.helpers, s.helper_altitude_km, s.helper_inclination_delta_deg,
                    s.dv_r_mps, s.dv_t_mps, s.dv_n_mps, s.activations
                )
            end
        end
        println("Summary CSV: $(opts.output_csv)")
        println("Time-series CSV: $(opts.timeseries_csv)")

    else
        # --- Single-run mode: run once with exactly what the user typed ---
        result = run_open_cavity_case(opts)
        s = result.summary

        _print_summary(s)  # print human-readable results to the terminal

        # Write summary and timeseries CSVs (returns "" if path is empty, so we skip the print)
        out    = _write_csv!(opts.output_csv,     [s];              append=opts.append_output)
        ts_out = _write_csv!(opts.timeseries_csv, result.timeseries; append=opts.append_output)
        isempty(out)    || println("Summary CSV: $out")
        isempty(ts_out) || println("Time-series CSV: $ts_out")

        # --- Diagnostic plots (mirrors prototype run_open_cavity_multi plot block) ---
        img_dir = joinpath(normpath(joinpath(REPO_ROOT, "output",
                                             "oracle_case2_laser_links")), "images")
        plot_open_cavity_results(result, opts; IMG_DIR=img_dir, target_only=false)

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
