#!/usr/bin/env julia
# ORACLE open-cavity laser-link example.
# Run: julia --project=<repo_root> examples/ORACLE_LaserLinks/run_oracle_laser_links.jl [options]
#
# Options:
#   --helpers N
#   --helper-altitude-km KM
#   --target-altitude-km KM
#   --target-inclination-deg DEG
#   --helper-inclination-deg DEG
#   --target-nu-deg DEG
#   --target-ecc VALUE
#   --orbits N
#   --schedule naive_next_entering|positive_along_track|gve_sma|gve_ecc|gve_inc|gve_raan|gve_argp
#   --laser-range-km KM
#   --laser-power-w W
#   --magnification B
#   --beta VALUE
#   --eta VALUE
#   --mass-kg KG
#   --dt-max-s SEC
#   --paper-grid
#   --feather-only
#   --output-dir PATH
#   --timeseries-points N
#   --animate            (requires GLMakie)

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const REPO_PROJECT = joinpath(REPO_ROOT, "Project.toml")

if something(Base.active_project(), "") != REPO_PROJECT
    import Pkg
    Pkg.activate(REPO_ROOT; io=devnull)
end

using SpaceAGORA
using SpaceAGORA.OracleAnalysis

using Printf

const DEFAULT_OUTPUT_DIR = joinpath(REPO_ROOT, "output", "oracle_laser_links")

# Scenario-specific runner (build_oracle_case_config, run_oracle_open_cavity_case,
# plot_open_cavity_results, _build_laser_save_fields).
include(joinpath(@__DIR__, "oracle_runner.jl"))

# Load GLMakie and the animation module only when --animate is passed.
const _HAS_GLMAKIE = if "--animate" in ARGS
    try
        @eval using GLMakie
        include(joinpath(@__DIR__, "..", "..", "src", "analysis", "visualization",
                         "laser_links", "laser_link_animation.jl"))
        true
    catch
        false
    end
else
    false
end

# =============================================================================
# CLI parsing
# =============================================================================

const _FLAG_OPTS   = (:paper_grid, :feather_only, :animate)
const _INT_OPTS    = (:helpers, :timeseries_points)
const _SYMBOL_OPTS = (:schedule, :planet)
const _PATH_OPTS   = (:output_dir,)
const _FLOAT_OPTS  = (
    :helper_altitude_km, :target_altitude_km, :target_inclination_deg,
    :helper_inclination_deg, :target_nu_deg, :target_ecc, :orbits,
    :laser_range_km, :laser_power_w, :magnification, :beta, :eta,
    :mass_kg, :dt_max_s,
)

function _parse_options(argv)::OracleOptions
    opts = Dict{Symbol, Any}(:output_dir => DEFAULT_OUTPUT_DIR)
    i = 1
    while i <= length(argv)
        arg = argv[i]
        arg in ("--help", "-h") && (println(_usage()); exit(0))
        startswith(arg, "--") || throw(ArgumentError("Unexpected argument '$arg'.\n$(_usage())"))
        key = Symbol(replace(arg[3:end], '-' => '_'))
        if key in _FLAG_OPTS
            opts[key] = true; i += 1; continue
        end
        i < length(argv) || throw(ArgumentError("Missing value for $arg."))
        val = argv[i + 1]
        if     key in _INT_OPTS    opts[key] = parse(Int, val)
        elseif key in _SYMBOL_OPTS opts[key] = Symbol(val)
        elseif key in _PATH_OPTS   opts[key] = abspath(val)
        elseif key in _FLOAT_OPTS  opts[key] = parse(Float64, val)
        else   throw(ArgumentError("Unknown option $arg.\n$(_usage())"))
        end
        i += 2
    end
    return OracleOptions(; opts...)
end

function _usage()
    return """
    Usage:
      julia --project=<repo> examples/ORACLE_LaserLinks/run_oracle_laser_links.jl [options]

    Options:
      --helpers N                    Number of helper satellites (default: 200)
      --helper-altitude-km KM
      --target-altitude-km KM
      --target-inclination-deg DEG
      --helper-inclination-deg DEG
      --target-nu-deg DEG
      --target-ecc VALUE
      --orbits N
      --schedule naive_next_entering|positive_along_track|gve_sma|gve_ecc|gve_inc|gve_raan|gve_argp
      --planet earth|mars|venus|titan|moon  (default: earth)
      --laser-range-km KM
      --laser-power-w W
      --magnification B
      --beta VALUE
      --eta VALUE
      --mass-kg KG
      --dt-max-s SEC
      --paper-grid               Sweep paper-grid parameter space
      --feather-only             Write feather+manifest only (skip CSV and plots)
      --output-dir PATH
      --timeseries-points N
      --animate                  Show 3D animation after simulation (requires GLMakie)
    """
end

# =============================================================================
# Main
# =============================================================================

function main(argv=ARGS)
    opts = _parse_options(argv)

    if opts.paper_grid
        for n_helpers in ORACLE_PAPER_HELPER_COUNTS
            for target_alt in ORACLE_PAPER_TARGET_ALTITUDES_KM
                for target_inc in ORACLE_PAPER_TARGET_INCLINATIONS_DEG
                    case_opts = _with(opts;
                        helpers=n_helpers,
                        target_altitude_km=target_alt,
                        target_inclination_deg=target_inc,
                        helper_altitude_km=ORACLE_PAPER_FIXED_HELPER_ALTITUDE_KM,
                        helper_inclination_deg=ORACLE_PAPER_FIXED_HELPER_INCLINATION_DEG,
                        output_dir=joinpath(opts.output_dir, "paper_plot_mode"),
                        animate=false,
                    )
                    elapsed = @elapsed begin
                        result = run_oracle_open_cavity_case(case_opts)
                    end
                    s = result.summary
                    @printf(
                        "helpers=%d target_alt_km=%.1f target_inc_deg=%.1f dv_R=%.6e dv_T=%.6e dv_N=%.6e activations=%d  [%.1f s]\n",
                        s.helpers, s.target_altitude_km, s.target_inclination_deg,
                        s.dv_r_mps, s.dv_t_mps, s.dv_n_mps, s.activations, elapsed
                    )
                    println("  → $(result.results_dir)")
                    result = nothing; GC.gc()
                end
            end
        end
        println("Output directory: $(opts.output_dir)")
    else
        single_opts = _with(opts; output_dir=joinpath(opts.output_dir, "single_case_mode"))
        elapsed = @elapsed begin
            result = run_oracle_open_cavity_case(single_opts)
        end
        s = result.summary
        _print_summary(s)
        @printf("  run time: %.1f s\n", elapsed)
        println("Output directory: $(result.results_dir)")

        img_dir = joinpath(result.results_dir, "images")
        opts.feather_only || plot_open_cavity_results(result, single_opts; IMG_DIR=img_dir, target_only=false)

        # Optional 3D animation.
        if opts.animate
            if _HAS_GLMAKIE
                _, anim_controls = LaserLinkAnimation.animate_oracle_satellites_3d(result.sol, opts.helpers)
                println("Animation window opened. Close the window to exit.")
                wait(anim_controls.screen)
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
