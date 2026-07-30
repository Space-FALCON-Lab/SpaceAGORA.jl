const EDITOR_ANALYSIS_MODE = "monte-carlo"

const EDITOR_ANALYSIS_MODES = (
    "monte-carlo",
    "sensitivity"
)

if !isdefined(@__MODULE__, :NavigationPaths)
    include(joinpath(@__DIR__, "navigation", "paths.jl"))
end
using .NavigationPaths: navigation_output_path, resolve_repo_path

const LEGACY_NOMINAL_OUTPUTS = (
    joinpath(
        NavigationPaths.REPO_ROOT,
        "output_navigation_pairwise_monte_carlo_legacy"
    ),
    joinpath(
        homedir(),
        "SpaceAGORA_outputs",
        "output_navigation_pairwise_monte_carlo_legacy"
    )
)
const LEGACY_SENSITIVITY_OUTPUTS = (
    joinpath(
        NavigationPaths.REPO_ROOT,
        "output_navigation_stress_monte_carlo"
    ),
    joinpath(
        homedir(),
        "SpaceAGORA_outputs",
        "output_navigation_stress_monte_carlo"
    )
)

const PYTHON_PLOT_OPTIONS = (
    "case",
    "input",
    "csv",
    "output",
    "max-target-plots"
)
const PYTHON_ANIMATION_OPTIONS = (
    "case",
    "input",
    "csv",
    "log",
    "preset",
    "show",
    "fps",
    "stride",
    "trail",
    "max-links-per-observer",
    "show-target-trails",
    "hide-target-trails",
    "hide-est-trails",
    "speed-factor",
    "interp-substeps",
    "camera-smooth",
    "event-focus-frames",
    "nav-dt",
    "nav-tol",
    "aligned-only",
    "earth",
    "visibility-km",
    "focus-scale",
    "near-margin",
    "hide-visibility-spheres"
)
const PYTHON_FLAG_OPTIONS = Set((
    "show",
    "show-target-trails",
    "hide-target-trails",
    "hide-est-trails",
    "aligned-only",
    "earth",
    "hide-visibility-spheres"
))

const TERMINAL_UTILITY_MODES = (
    "comparison",
    "iod-geometry",
    "iod-pairwise",
    "plots",
    "animation"
)
const INTERNAL_ANALYSIS_MODES = (
    "monte-carlo",
    "sensitivity-monte-carlo",
    TERMINAL_UTILITY_MODES...
)

function _editor_analysis_options()::Dict{String, String}
    selected = replace(
        lowercase(strip(string(EDITOR_ANALYSIS_MODE))),
        '_' => '-',
        ' ' => '-'
    )
    selected in EDITOR_ANALYSIS_MODES ||
        error("EDITOR_ANALYSIS_MODE must be one of $(EDITOR_ANALYSIS_MODES).")
    return Dict("mode" => selected)
end

function _analysis_options(args::Vector{String})::Dict{String, String}
    isempty(args) && return _editor_analysis_options()
    options = Dict{String, String}()
    for arg in args
        if startswith(arg, "--")
            key_value = split(arg[3:end], '='; limit=2)
            options[key_value[1]] = length(key_value) == 2 ? key_value[2] : "true"
        elseif !haskey(options, "mode")
            options["mode"] = arg
        else
            error("Unexpected positional argument: $(arg)")
        end
    end
    return options
end

function _print_analysis_help()::Nothing
    println("""
    Usage:
      julia --project=. examples/analyze_navigation.jl --mode=<mode> [options]

    With no command-line arguments, set EDITOR_ANALYSIS_MODE at the top of the
    file to "monte-carlo" or "sensitivity", then use the editor Play button.
    These choices analyze output/navigation/nominal and
    output/navigation/sensitivity, respectively.

    Campaign analysis modes (editor Play button or command line):
      monte-carlo
      sensitivity

    Terminal-only utilities:
      comparison
      iod-geometry
      iod-pairwise
      plots                       Static plots from one saved run
      animation                   3D animation from one saved run

    Options:
      --input=<path>       Campaign or single-run directory
      --output=<path>      Plot directory, pairwise output, or animation file
      --case=<name>        Navigation case; defaults to no-da for plots
      --csv=<path>         simulation_results.csv override for plots/animation
      --log=<path>         Optional run log for animation
      --max-target-plots=<n>
                          Representative single-target plots (default: 3)
      --show=<bool>        Show the animation while rendering
      --replay=<bool>      Replay missing truth geometry

    Examples:
      julia --project=. examples/analyze_navigation.jl --mode=monte-carlo
      julia --project=. examples/analyze_navigation.jl --mode=sensitivity
      julia --project=. examples/analyze_navigation.jl --mode=iod-geometry
      julia --project=. examples/analyze_navigation.jl --mode=iod-pairwise
      julia --project=. examples/analyze_navigation.jl --mode=plots --case=no-da
      julia --project=. examples/analyze_navigation.jl --mode=animation --case=no-da --output=/tmp/no_da.mp4

    The remaining modes are command-line utilities used by the commands printed
    after single runs or for optional diagnostics.
    """)
    return nothing
end

function _set_option!(
    options::Dict{String, String},
    option::String,
    environment_name::String
)::Nothing
    haskey(options, option) && (ENV[environment_name] = options[option])
    return nothing
end

function _campaign_input_root(
    input::String,
    environment_name::String,
    marker_file::String,
    default_output::String,
    legacy_outputs::Tuple
)::String
    explicit_input = strip(input)
    if !isempty(explicit_input)
        selected = resolve_repo_path(explicit_input)
        isfile(joinpath(selected, marker_file)) || error(
            "Campaign input does not contain $(marker_file): $(selected)"
        )
        return selected
    end

    environment_input = strip(get(ENV, environment_name, ""))
    if !isempty(environment_input)
        selected = resolve_repo_path(environment_input)
        isfile(joinpath(selected, marker_file)) || error(
            "$(environment_name) does not contain $(marker_file): $(selected)"
        )
        return selected
    end

    candidates = (default_output, legacy_outputs...)
    selected_index = findfirst(
        candidate -> isfile(joinpath(candidate, marker_file)),
        candidates
    )
    selected_index !== nothing && return normpath(candidates[selected_index])

    checked = join(["  - $(candidate)" for candidate in candidates], "\n")
    error(
        "Could not find $(marker_file). Checked:\n$(checked)\n" *
        "Download the campaign or pass --input=<campaign-directory>."
    )
end

@inline function _option_is_true(value::String)::Bool
    return lowercase(strip(value)) in ("1", "true", "yes", "y", "on")
end

@inline function _option_is_false(value::String)::Bool
    return lowercase(strip(value)) in ("0", "false", "no", "n", "off")
end

function _clear_editor_input_override!(mode::String)::Nothing
    mode == "monte-carlo" && pop!(ENV, "SPACEAGORA_MC_OUTPUT", nothing)
    mode == "sensitivity-monte-carlo" &&
        pop!(ENV, "SPACEAGORA_STRESS_OUTPUT", nothing)
    return nothing
end

function _run_python_analysis(
    mode::String,
    options::Dict{String, String}
)::Nothing
    script_name, supported_options = if mode == "plots"
        ("plot_results.py", PYTHON_PLOT_OPTIONS)
    else
        ("animate_results_3d.py", PYTHON_ANIMATION_OPTIONS)
    end
    script_path = normpath(
        joinpath(@__DIR__, "..", "scripts", "plotting", "navigation", script_name)
    )
    isfile(script_path) || error("Navigation plotting script not found: $(script_path)")

    python_override = strip(get(ENV, "PYTHON", ""))
    repo_root = normpath(joinpath(@__DIR__, ".."))
    venv_python = Sys.iswindows() ?
        joinpath(repo_root, ".venv-navigation", "Scripts", "python.exe") :
        joinpath(repo_root, ".venv-navigation", "bin", "python3")
    python_path = if !isempty(python_override)
        isfile(python_override) ? python_override : Sys.which(python_override)
    elseif isfile(venv_python)
        venv_python
    else
        Sys.which("python3")
    end
    python_path === nothing && error(
        "Python executable not found. Create .venv-navigation as described " *
        "in examples/navigation/README.md."
    )
    command_parts = String[python_path, script_path]

    for option in supported_options
        haskey(options, option) || continue
        value = options[option]
        if option in PYTHON_FLAG_OPTIONS
            _option_is_true(value) && push!(command_parts, "--$(option)")
            (_option_is_true(value) || _option_is_false(value)) ||
                error("--$(option) must be true or false")
        else
            push!(command_parts, "--$(option)=$(value)")
        end
    end

    if mode == "animation" && haskey(options, "output")
        push!(command_parts, "--save=$(options["output"])")
    end

    allowed_options = Set(("mode", "help", "h", supported_options...))
    mode == "animation" && push!(allowed_options, "output")
    unsupported = sort([key for key in keys(options) if !(key in allowed_options)])
    isempty(unsupported) ||
        error("Unsupported option(s) for --mode=$(mode): $(join(unsupported, ", "))")

    println("  python script: $(script_path)")
    run(Cmd(command_parts))
    return nothing
end

function _select_analysis!(
    mode::String,
    options::Dict{String, String}
)::Union{Nothing, String}
    input = get(options, "input", "")
    output = get(options, "output", "")

    if mode == "comparison"
        !isempty(input) && (ENV["SPACEAGORA_COMPARISON_OUTPUT"] = input)
        !isempty(output) && (ENV["SPACEAGORA_PLOT_OUTPUT"] = output)
        _set_option!(options, "case", "SPACEAGORA_CASES")
        return joinpath(@__DIR__, "navigation", "analysis", "comparison.jl")
    elseif mode == "monte-carlo"
        campaign_root = _campaign_input_root(
            input,
            "SPACEAGORA_MC_OUTPUT",
            "mc_run_status.csv",
            navigation_output_path("nominal"),
            LEGACY_NOMINAL_OUTPUTS
        )
        ENV["SPACEAGORA_MC_OUTPUT"] = campaign_root
        println("  input: $(campaign_root)")
        return joinpath(@__DIR__, "navigation", "analysis", "monte_carlo.jl")
    elseif mode == "sensitivity-monte-carlo"
        campaign_root = _campaign_input_root(
            input,
            "SPACEAGORA_STRESS_OUTPUT",
            "stress_run_status.csv",
            navigation_output_path("sensitivity"),
            LEGACY_SENSITIVITY_OUTPUTS
        )
        ENV["SPACEAGORA_STRESS_OUTPUT"] = campaign_root
        ENV["SPACEAGORA_STRESS_ANALYZE_ONLY"] = "true"
        println("  input: $(campaign_root)")
        return joinpath(
            @__DIR__,
            "navigation",
            "run",
            "sensitivity_monte_carlo.jl"
        )
    elseif mode == "iod-geometry"
        !isempty(input) && (ENV["SPACEAGORA_MC_OUTPUT"] = input)
        _set_option!(options, "case", "SPACEAGORA_MC_IOD_GEOMETRY_CASE")
        _set_option!(options, "replay", "SPACEAGORA_MC_IOD_GEOMETRY_REPLAY")
        return joinpath(@__DIR__, "navigation", "analysis", "iod_geometry.jl")
    elseif mode == "iod-pairwise"
        empty!(ARGS)
        if !isempty(input)
            push!(ARGS, input)
        elseif !isempty(output)
            push!(ARGS, navigation_output_path("iod_pairwise"))
        end
        if !isempty(output)
            push!(ARGS, output)
        end
        return joinpath(
            @__DIR__,
            "navigation",
            "analysis",
            "iod_pairwise.jl"
        )
    elseif mode == "plots" || mode == "animation"
        _run_python_analysis(mode, options)
        return nothing
    end
    error(
        "Unknown navigation analysis mode: $(mode). " *
        "Valid internal modes: $(INTERNAL_ANALYSIS_MODES)"
    )
end

function main()::Nothing
    using_editor_config = isempty(ARGS)
    options = _analysis_options(ARGS)
    if haskey(options, "help") || haskey(options, "h")
        _print_analysis_help()
        return nothing
    end
    mode = replace(
        get(options, "mode", get(ENV, "SPACEAGORA_NAV_ANALYSIS_MODE", "monte-carlo")),
        '_' => '-',
        ' ' => '-'
    )
    mode == "sensitivity" && (mode = "sensitivity-monte-carlo")
    mode in INTERNAL_ANALYSIS_MODES ||
        error(
            "Unknown navigation analysis mode: $(mode). " *
            "Campaign modes: $(EDITOR_ANALYSIS_MODES); " *
            "terminal utilities: $(TERMINAL_UTILITY_MODES)"
        )
    using_editor_config && _clear_editor_input_override!(mode)
    display_mode = mode == "sensitivity-monte-carlo" ? "sensitivity" : mode
    println("Navigation analysis mode: $(display_mode)")
    using_editor_config && println("  selection: EDITOR_ANALYSIS_MODE")
    selected_file = _select_analysis!(mode, options)
    selected_file === nothing || include(selected_file)
    if mode == "iod-geometry"
        geometry_module = Base.invokelatest(
            getfield,
            @__MODULE__,
            :MonteCarloIODGeometryAnalysis
        )
        geometry_main = Base.invokelatest(getfield, geometry_module, :main)
        Base.invokelatest(geometry_main)
    end
    return nothing
end

main()
