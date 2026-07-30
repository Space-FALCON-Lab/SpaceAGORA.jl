const EDITOR_RUN_CONFIG = (
    mode="single",
    case="independent_local_da",
    output="",
    mission_time=600,
    targets=300,
    runs=1,
    resume=true,
    sweeps="false_alarm",
)

const RUN_MODES = (
    "single",
    "comparison",
    "monte-carlo",
    "sensitivity-sweep",
    "sensitivity-monte-carlo"
)

function _editor_run_options()::Dict{String, String}
    config = EDITOR_RUN_CONFIG
    options = Dict(
        "mode" => string(config.mode),
        "case" => string(config.case),
        "mission-time" => string(config.mission_time),
        "targets" => string(config.targets),
        "runs" => string(config.runs),
        "resume" => string(config.resume),
        "sweeps" => string(config.sweeps),
    )
    !isempty(strip(string(config.output))) &&
        (options["output"] = string(config.output))
    return options
end

function _run_options(args::Vector{String})::Dict{String, String}
    isempty(args) && return _editor_run_options()
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

function _print_run_help()::Nothing
    println("""
    Usage:
      julia --project=. examples/run_navigation.jl --mode=<mode> [options]

    With no command-line arguments, EDITOR_RUN_CONFIG at the top of this file
    is used. Configure it and launch the file with the editor Play button.

    Modes:
      single                     One navigation case; output is temporary by default
      comparison                 One realization for multiple methods; temporary by default
      monte-carlo                Nominal Monte Carlo campaign
      sensitivity-sweep          One-realization sensitivity sweep; temporary by default
      sensitivity-monte-carlo    Monte Carlo measurement-fault campaign

    Single runs print their metrics and plot/animation commands. Comparison and
    sensitivity-sweep print compact results for every executed configuration.
    Complete percentile analyses are produced from the two Monte Carlo outputs.

    Common options:
      --case=<name>              One case or a comma-separated case list
      --output=<path>            Override the mode output directory
      --runs=<count>             Number of Monte Carlo realizations
      --first-run=<index>        First Monte Carlo realization index
      --resume=<true|false>      Resume completed Monte Carlo runs
      --mission-time=<seconds>   Mission duration
      --targets=<count>          Number of clustered targets
      --iod-pairwise=<bool>      Use the dedicated nominal pairwise IOD campaign
      --iod-geometry=<bool>      Save IOD geometry under the nominal campaign
      --sweeps=<list>            Sensitivity sweeps: bias, misdetection, false_alarm

    Examples:
      julia --project=. examples/run_navigation.jl --mode=monte-carlo --case=proposed --runs=100
      julia --project=. examples/run_navigation.jl --mode=sensitivity-monte-carlo --runs=100
      julia --project=. examples/run_navigation.jl --mode=single --case=proposed
      julia --project=. examples/run_navigation.jl --mode=single --case=no-da

    Default output roots:
      nominal MC                 output/navigation/nominal
      sensitivity MC             output/navigation/sensitivity
      single/debug runs          system temporary directory

    Existing SPACEAGORA_* variables remain valid for advanced configuration.
    """)
    return nothing
end

function _absolute_output_path(value::AbstractString)::String
    path = expanduser(String(value))
    return normpath(isabspath(path) ? path : joinpath(@__DIR__, "..", path))
end

@inline _truthy_option(value::AbstractString)::Bool =
    lowercase(strip(value)) in ("1", "true", "yes", "y", "on")

function _output_environment_name(mode::String)::String
    if mode in ("single", "comparison")
        return "SPACEAGORA_COMPARISON_OUTPUT"
    elseif mode == "sensitivity-sweep"
        return "SPACEAGORA_SENSITIVITY_OUTPUT"
    elseif mode == "monte-carlo"
        return "SPACEAGORA_MC_OUTPUT"
    elseif mode == "sensitivity-monte-carlo"
        return "SPACEAGORA_STRESS_OUTPUT"
    end
    error("No output environment variable for navigation mode $(mode)")
end

function _default_output_root(mode::String)::String
    repo_root = normpath(joinpath(@__DIR__, ".."))
    if mode == "single"
        return joinpath(tempdir(), "spaceagora_navigation", "runs", "single")
    elseif mode == "comparison"
        return joinpath(tempdir(), "spaceagora_navigation", "runs", "comparison")
    elseif mode == "sensitivity-sweep"
        return joinpath(
            tempdir(),
            "spaceagora_navigation",
            "runs",
            "sensitivity_sweep"
        )
    elseif mode == "monte-carlo"
        return joinpath(repo_root, "output", "navigation", "nominal")
    elseif mode == "sensitivity-monte-carlo"
        return joinpath(repo_root, "output", "navigation", "sensitivity")
    end
    error("No default output root for navigation mode $(mode)")
end

function _run_output_root(
    mode::String,
    options::Dict{String, String};
    use_environment::Bool=true
)::String
    haskey(options, "output") &&
        return _absolute_output_path(options["output"])
    if mode == "monte-carlo" &&
       _truthy_option(get(options, "iod-pairwise", "false"))
        return normpath(joinpath(
            @__DIR__,
            "..",
            "output",
            "navigation",
            "iod_pairwise"
        ))
    end
    if use_environment
        environment_output = strip(get(
            ENV,
            _output_environment_name(mode),
            ""
        ))
        !isempty(environment_output) &&
            return _absolute_output_path(environment_output)
    end
    return _default_output_root(mode)
end

function _print_editor_follow_up(
    mode::String,
    options::Dict{String, String}
)::Nothing
    output_root = _run_output_root(mode, options; use_environment=false)
    repo_root = normpath(joinpath(@__DIR__, ".."))
    analysis_script = joinpath(repo_root, "examples", "analyze_navigation.jl")

    println()
    println("Configured editor run completed")
    println("  output: $(output_root)")

    if mode == "single"
        case_name = replace(options["case"], '-' => '_')
        results_directory = joinpath(output_root, case_name)
        plots_directory = joinpath(output_root, "plots")
        animation_path = joinpath(output_root, "animation", "$(case_name).gif")
        if isfile(joinpath(results_directory, "simulation_results.csv"))
            println()
            println("Copy and paste to generate the plots:")
            println(
                "julia --project=\"$(repo_root)\" \"$(analysis_script)\" " *
                "--mode=plots --case=$(case_name) " *
                "--input=\"$(results_directory)\" --output=\"$(plots_directory)\""
            )
            println()
            println("After generation, copy and paste to view them:")
            println(
                "code \"$(joinpath(plots_directory, "target_scenario.png"))\" " *
                "\"$(joinpath(plots_directory, "observer_initial_configuration_reference_orbits.png"))\" " *
                "\"$(joinpath(plots_directory, "rmse_all_tracks.png"))\" " *
                "\"$(joinpath(plots_directory, "target_errors"))\"/*.png"
            )
            println()
            println("Copy and paste to generate the animation:")
            println(
                "julia --project=\"$(repo_root)\" \"$(analysis_script)\" " *
                "--mode=animation --case=$(case_name) " *
                "--input=\"$(results_directory)\" --output=\"$(animation_path)\" " *
                "--interp-substeps=2 --hide-target-trails=true " *
                "--hide-est-trails=true --hide-visibility-spheres=true"
            )
            println()
            println("After rendering, copy and paste to view it:")
            println("code \"$(animation_path)\"")
        end
    elseif mode == "monte-carlo"
        println()
        println("Copy and paste to analyze the Monte Carlo results:")
        println(
            "julia --project=\"$(repo_root)\" \"$(analysis_script)\" " *
            "--mode=monte-carlo --input=\"$(output_root)\""
        )
    elseif mode == "sensitivity-monte-carlo"
        println()
        println("Copy and paste to analyze the sensitivity results:")
        println(
            "julia --project=\"$(repo_root)\" \"$(analysis_script)\" " *
            "--mode=sensitivity --input=\"$(output_root)\""
        )
    end
    return nothing
end

function _set_if_present!(
    options::Dict{String, String},
    option::String,
    environment_name::String
)::Nothing
    haskey(options, option) && (ENV[environment_name] = options[option])
    return nothing
end

function _configure_common!(options::Dict{String, String})::Nothing
    _set_if_present!(options, "mission-time", "SPACEAGORA_MISSION_TIME_SEC")
    _set_if_present!(options, "targets", "SPACEAGORA_N_CLUSTER_TARGETS")
    return nothing
end

function _configure_saved_outputs!(
    mode::String,
    options::Dict{String, String}
)::Nothing
    save_timeseries = string(mode == "single")
    ENV["SPACEAGORA_SAVE_SIMULATION_RESULTS"] = save_timeseries
    ENV["SPACEAGORA_SAVE_TARGET_ESTIMATE_FIELDS"] = save_timeseries

    if mode == "monte-carlo"
        iod_pairwise = get(options, "iod-pairwise", "false")
        iod_geometry = get(options, "iod-geometry", "false")
        ENV["SPACEAGORA_SAVE_IOD_PAIRWISE_DIAGNOSTICS"] = iod_pairwise
        ENV["SPACEAGORA_SAVE_IOD_EVENT_GEOMETRY"] = iod_geometry
        ENV["SPACEAGORA_MC_IOD_GEOMETRY_ANALYSIS"] = iod_geometry
    elseif mode == "sensitivity-monte-carlo"
        requested_pairwise = _truthy_option(get(options, "iod-pairwise", "false"))
        requested_geometry = _truthy_option(get(options, "iod-geometry", "false"))
        (requested_pairwise || requested_geometry) && error(
            "Detailed IOD geometry and pairwise diagnostics use the dedicated " *
            "nominal campaigns"
        )
        ENV["SPACEAGORA_SAVE_IOD_PAIRWISE_DIAGNOSTICS"] = "false"
        ENV["SPACEAGORA_SAVE_IOD_EVENT_GEOMETRY"] = "false"
    else
        ENV["SPACEAGORA_SAVE_IOD_PAIRWISE_DIAGNOSTICS"] = "false"
        ENV["SPACEAGORA_SAVE_IOD_EVENT_GEOMETRY"] = "false"
    end
    return nothing
end

function _configure_mode!(
    mode::String,
    options::Dict{String, String};
    use_environment_output::Bool=true
)::String
    _configure_common!(options)
    _configure_saved_outputs!(mode, options)
    case_value = replace(get(options, "case", ""), '-' => '_')
    output_root = _run_output_root(
        mode,
        options;
        use_environment=use_environment_output
    )
    ENV[_output_environment_name(mode)] = output_root

    if mode == "single"
        !isempty(case_value) && (ENV["SPACEAGORA_NAV_CASE"] = case_value)
        return joinpath(@__DIR__, "navigation.jl")
    elseif mode == "comparison"
        !isempty(case_value) && (ENV["SPACEAGORA_CASES"] = case_value)
        return joinpath(@__DIR__, "navigation", "run", "comparison.jl")
    elseif mode == "monte-carlo"
        if _truthy_option(get(options, "iod-pairwise", "false"))
            isempty(case_value) && (case_value = "proposed")
            case_value == "proposed" || error(
                "The dedicated pairwise IOD campaign supports case=proposed"
            )
        end
        !isempty(case_value) && (ENV["SPACEAGORA_MC_CASES"] = case_value)
        _set_if_present!(options, "runs", "SPACEAGORA_MC_RUNS")
        _set_if_present!(options, "first-run", "SPACEAGORA_MC_FIRST_RUN")
        _set_if_present!(options, "resume", "SPACEAGORA_MC_RESUME")
        return joinpath(@__DIR__, "navigation", "run", "monte_carlo.jl")
    elseif mode == "sensitivity-sweep"
        return joinpath(@__DIR__, "navigation", "run", "sensitivity_sweep.jl")
    elseif mode == "sensitivity-monte-carlo"
        !isempty(case_value) && (ENV["SPACEAGORA_STRESS_CASES"] = case_value)
        ENV["SPACEAGORA_STRESS_ANALYZE_ONLY"] = "false"
        _set_if_present!(options, "sweeps", "SPACEAGORA_STRESS_SWEEPS")
        _set_if_present!(options, "runs", "SPACEAGORA_STRESS_RUNS")
        _set_if_present!(options, "first-run", "SPACEAGORA_STRESS_FIRST_RUN")
        _set_if_present!(options, "resume", "SPACEAGORA_STRESS_RESUME")
        return joinpath(
            @__DIR__,
            "navigation",
            "run",
            "sensitivity_monte_carlo.jl"
        )
    end
    error("Unknown navigation run mode: $(mode). Valid modes: $(RUN_MODES)")
end

function main()::Nothing
    using_editor_config = isempty(ARGS)
    options = _run_options(ARGS)
    if haskey(options, "help") || haskey(options, "h")
        _print_run_help()
        return nothing
    end
    mode = replace(
        get(options, "mode", get(ENV, "SPACEAGORA_NAV_RUN_MODE", "monte-carlo")),
        '_' => '-'
    )
    mode in RUN_MODES ||
        error("Unknown navigation run mode: $(mode). Valid modes: $(RUN_MODES)")
    selected_file = _configure_mode!(
        mode,
        options;
        use_environment_output=!using_editor_config
    )
    println("Navigation run mode: $(mode)")
    using_editor_config && println("  configuration: EDITOR_RUN_CONFIG")
    include(selected_file)
    using_editor_config && _print_editor_follow_up(mode, options)
    return nothing
end

main()
