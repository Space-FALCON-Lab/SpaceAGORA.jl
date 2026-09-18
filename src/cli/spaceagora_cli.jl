module SpaceAGORACLI
using ..SimulationModel: SimulationModel

export AssetCheckItem, AssetCheckReport, check_assets, render_asset_report, run_cli

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
# Every child process runs under the repository project, the environment whose
# Project.toml is SpaceAGORA itself and the one `julia --project=.` selects.
# The children used to be launched with `--project=<repo>/.AGORA`, an
# untracked, machine-local directory that no setup step creates any more; a
# child started with a missing project directory gets an empty project and
# cannot load SpaceAGORA. Examples that include examples/common.jl hid the
# defect by re-activating the repository project themselves; the telemetry and
# benchmark launchers, and any example that loads the package directly, did not.
const CHILD_PROJECT = REPO_ROOT
const EXAMPLES_DIR = joinpath(REPO_ROOT, "examples")
const TELEMETRY_LAUNCHER = joinpath(REPO_ROOT, "benchmarks", "studies", "telemetry_orbit_accuracy_study.jl")
const PERF_RUNTIME_LAUNCHER = joinpath(REPO_ROOT, "benchmarks", "studies", "performance_runtime_analysis.jl")
const SMART_LADDER_LAUNCHER = joinpath(REPO_ROOT, "benchmarks", "studies", "performance_smart_parallel_ladder.jl")

include(joinpath(@__DIR__, "assets.jl"))

@inline function _starts_with(arg::String, prefix::String)
    return startswith(arg, prefix)
end

@inline function _value_after_equals(arg::String, prefix::String)
    return split(arg, "=", limit=2)[2]
end

@inline function _normalize_example_name(name::String)::String
    token = strip(name)
    endswith(token, ".jl") && return token
    return string(token, ".jl")
end
@inline _normalize_example_name(name::AbstractString) = _normalize_example_name(String(name))

function _resolve_example_path(name::String)::String
    direct = abspath(name)
    if isfile(direct)
        return direct
    end
    candidate = joinpath(EXAMPLES_DIR, _normalize_example_name(name))
    isfile(candidate) || throw(ArgumentError("Example not found: $(name). Checked $(candidate)."))
    return candidate
end
@inline _resolve_example_path(name::AbstractString) = _resolve_example_path(String(name))

function _print_usage(io::IO=stdout)
    println(io, "Usage:")
    println(io, "  spaceagora run --example=<file> [--output-dir=<dir>] [--smoke] [--visualize] [--print-only]")
    println(io, "  spaceagora visualize --run=<results dir or bundle prefix> [--out=<html>] [--max-frames=<n>] [--frame=inertial|planet_fixed] [--texture=best|4k|8k] [--trail-orbits=<n>] [--model=<id>=<stl|obj|glb>] [--model-scale=<m per unit>] [--model-center=0|1] [--ensemble]")
    println(io, "  spaceagora telemetry [quick|full|smoke] [--output-dir=<dir>] [--enforce=0|1] [--plots=0|1] [--print-only]")
    println(io, "  spaceagora benchmark runtime-analysis [quick|full|smoke] [--output-dir=<dir>] [--print-only]")
    println(io, "  spaceagora benchmark smart-parallel-ladder [quick|full|smoke] [--output-dir=<dir>] [--print-only]")
    println(io, "  spaceagora assets check")
    println(io, "  spaceagora assets manifest")
    println(io, "  spaceagora assets setup-open")
    return 0
end

# The one place the child command is built. `--print-only` renders exactly this
# object, so what is printed is what launches, and the CLI tests compare its
# arguments directly rather than parsing the printed line (a path with spaces
# is one argument here and a broken one in any whitespace-split rendering).
function _child_command(script::String, script_args::Vector{String})::Cmd
    return `$(Base.julia_cmd()) --project=$CHILD_PROJECT $script $script_args`
end

function _run_subprocess(script::String, script_args::Vector{String}; env_pairs::Vector{Pair{String,String}}=Pair{String,String}[], print_only::Bool=false, io::IO=stdout, errio::IO=stderr)::Int
    full = _child_command(script, script_args)
    if print_only
        println(io, "project=$(CHILD_PROJECT)")
        println(io, "script=$(script)")
        if !isempty(env_pairs)
            println(io, "env:")
            for (k, v) in env_pairs
                println(io, "  $(k)=$(v)")
            end
        end
        println(io, "cmd=$(full)")
        return 0
    end
    mkpath(joinpath(REPO_ROOT, "output"))
    cmd_env = isempty(env_pairs) ? full : addenv(full, env_pairs...)
    process = run(pipeline(ignorestatus(cmd_env); stdout=io, stderr=errio), wait=true)
    return process.exitcode
end

function _run_example(args::Vector{String}; io::IO=stdout, errio::IO=stderr)::Int
    example = ""
    output_dir = ""
    smoke = false
    visualize = false
    print_only = false
    for arg in args
        if _starts_with(arg, "--example=")
            example = _value_after_equals(arg, "--example=")
        elseif _starts_with(arg, "--output-dir=")
            output_dir = abspath(_value_after_equals(arg, "--output-dir="))
        elseif arg == "--smoke"
            smoke = true
        elseif arg == "--visualize"
            visualize = true
        elseif arg == "--print-only"
            print_only = true
        else
            throw(ArgumentError("Unknown run argument '$arg'."))
        end
    end
    isempty(example) && throw(ArgumentError("run requires --example=<file>."))
    script = _resolve_example_path(example)
    env_pairs = Pair{String,String}[]
    !isempty(output_dir) && push!(env_pairs, "SPACEAGORA_CLI_OUTPUT_DIR" => output_dir)
    if smoke
        push!(env_pairs, "SPACEAGORA_EXAMPLE_SMOKE" => "1")
        push!(env_pairs, "SPACEAGORA_EXAMPLE_SMOKE_RESULTS" => "1")
    end
    # The example scripts call run_simulation themselves; the environment switch
    # turns on the scene sidecar and builds the viewer page after the run.
    visualize && push!(env_pairs, "SPACEAGORA_VISUALIZATION" => "1")
    return _run_subprocess(script, String[]; env_pairs=env_pairs, print_only=print_only, io=io, errio=errio)
end

# `--run` accepts a results directory (holding simulation_results.feather and
# its scene sidecar), the bundle prefix itself, or, with --ensemble, a campaign
# directory of sample subdirectories.
function _run_visualize(args::Vector{String}; io::IO=stdout, errio::IO=stderr)::Int
    run = ""
    out = nothing
    ensemble = false
    models = Dict{Int, String}()
    kwargs = Pair{Symbol, Any}[]
    for arg in args
        if _starts_with(arg, "--run=")
            run = _value_after_equals(arg, "--run=")
        elseif _starts_with(arg, "--out=")
            out = abspath(_value_after_equals(arg, "--out="))
        elseif _starts_with(arg, "--max-frames=")
            push!(kwargs, :max_frames => parse(Int, _value_after_equals(arg, "--max-frames=")))
        elseif _starts_with(arg, "--frame=")
            push!(kwargs, :frame => Symbol(_value_after_equals(arg, "--frame=")))
        elseif _starts_with(arg, "--texture=")
            value = lowercase(strip(_value_after_equals(arg, "--texture=")))
            push!(kwargs, :texture_resolution => (value == "best" ? :best : value))
        elseif _starts_with(arg, "--trail-orbits=")
            push!(kwargs, :trail_orbits => parse(Float64, _value_after_equals(arg, "--trail-orbits=")))
        elseif _starts_with(arg, "--title=")
            push!(kwargs, :title => _value_after_equals(arg, "--title="))
        elseif arg == "--no-textures"
            push!(kwargs, :textures => false)
        elseif _starts_with(arg, "--model=")
            spec = _value_after_equals(arg, "--model=")
            occursin("=", spec) || throw(ArgumentError("--model expects <spacecraft id>=<file>, got '$spec'."))
            id_text, path = split(spec, "="; limit=2)
            models[parse(Int, id_text)] = String(path)
        elseif _starts_with(arg, "--model-scale=")
            push!(kwargs, :model_scale => parse(Float64, _value_after_equals(arg, "--model-scale=")))
        elseif _starts_with(arg, "--model-center=")
            push!(kwargs, :model_center => lowercase(strip(_value_after_equals(arg, "--model-center="))) in ("1", "true", "yes", "on"))
        elseif arg == "--ensemble"
            ensemble = true
        else
            throw(ArgumentError("Unknown visualize argument '$arg'."))
        end
    end
    isempty(run) && throw(ArgumentError("visualize requires --run=<results dir or bundle prefix>."))
    run = abspath(run)
    out === nothing || push!(kwargs, :out => out)
    if !isempty(models)
        ensemble && throw(ArgumentError("--model applies to single-run pages, not --ensemble."))
        push!(kwargs, :models => models)
    end
    page = if ensemble
        SimulationModel.SceneVisualization.export_ensemble_visualization(run; kwargs...)
    else
        prefix = isdir(run) ? joinpath(run, "simulation_results") : run
        SimulationModel.SceneVisualization.export_visualization(prefix; kwargs...)
    end
    println(io, "viewer=$(page)")
    return 0
end

function _run_telemetry(args::Vector{String}; io::IO=stdout, errio::IO=stderr)::Int
    profile = "quick"
    output_dir = joinpath(REPO_ROOT, "output", "telemetry")
    enforce = false
    generate_plots = false
    print_only = false
    scenarios = ""
    for arg in args
        if arg in ("quick", "full")
            profile = arg
        elseif arg == "smoke"
            profile = "quick"
        elseif _starts_with(arg, "--profile=")
            profile = lowercase(strip(_value_after_equals(arg, "--profile=")))
            profile == "smoke" && (profile = "quick")
        elseif _starts_with(arg, "--output-dir=")
            output_dir = abspath(_value_after_equals(arg, "--output-dir="))
        elseif _starts_with(arg, "--enforce=")
            enforce = lowercase(strip(_value_after_equals(arg, "--enforce="))) in ("1", "true", "yes", "on")
        elseif _starts_with(arg, "--plots=")
            generate_plots = lowercase(strip(_value_after_equals(arg, "--plots="))) in ("1", "true", "yes", "on")
        elseif _starts_with(arg, "--scenarios=")
            scenarios = strip(_value_after_equals(arg, "--scenarios="))
        elseif arg == "--print-only"
            print_only = true
        else
            throw(ArgumentError("Unknown telemetry argument '$arg'."))
        end
    end
    mkpath(output_dir)
    env_pairs = Pair{String,String}[
        "SPACEAGORA_TELEMETRY_OUT_SUMMARY" => joinpath(output_dir, "telemetry_orbit_accuracy_summary.csv"),
        "SPACEAGORA_TELEMETRY_OUT_ERRORS" => joinpath(output_dir, "telemetry_orbit_accuracy_errors.csv"),
        "SPACEAGORA_TELEMETRY_PLOTS" => (generate_plots ? "1" : "0"),
        "SPACEAGORA_TELEMETRY_ENFORCE" => (enforce ? "1" : "0"),
        "SPACEAGORA_TELEMETRY_SCENARIOS" => String(scenarios),
    ]
    return _run_subprocess(TELEMETRY_LAUNCHER, [profile]; env_pairs=env_pairs, print_only=print_only, io=io, errio=errio)
end

function _run_benchmark(args::Vector{String}; io::IO=stdout, errio::IO=stderr)::Int
    isempty(args) && throw(ArgumentError("benchmark requires a mode: runtime-analysis or smart-parallel-ladder."))
    mode = first(args)
    tail = args[2:end]
    output_dir = ""
    print_only = false
    script_args = String[]
    launcher = ""
    env_pairs = Pair{String,String}[]
    if mode == "runtime-analysis"
        launcher = PERF_RUNTIME_LAUNCHER
        output_dir = joinpath(REPO_ROOT, "output", "performance")
        for arg in tail
            if _starts_with(arg, "--output-dir=")
                output_dir = abspath(_value_after_equals(arg, "--output-dir="))
            elseif arg == "--print-only"
                print_only = true
            else
                push!(script_args, arg)
            end
        end
        push!(env_pairs, "SPACEAGORA_PERF_OUTDIR" => output_dir)
    elseif mode == "smart-parallel-ladder"
        launcher = SMART_LADDER_LAUNCHER
        output_dir = joinpath(REPO_ROOT, "output", "performance", "smart_parallel_ladder")
        for arg in tail
            if _starts_with(arg, "--output-dir=")
                output_dir = abspath(_value_after_equals(arg, "--output-dir="))
            elseif arg == "--print-only"
                print_only = true
            else
                push!(script_args, arg)
            end
        end
        push!(env_pairs, "SPACEAGORA_SMART_LADDER_OUTDIR" => output_dir)
    else
        throw(ArgumentError("Unsupported benchmark mode '$mode'. Use runtime-analysis or smart-parallel-ladder."))
    end
    mkpath(output_dir)
    return _run_subprocess(launcher, script_args; env_pairs=env_pairs, print_only=print_only, io=io, errio=errio)
end

"""
    run_cli([args=ARGS]; io=stdout, errio=stderr) -> Int

Stable CLI entrypoint for SpaceAGORA operational commands:

- `run` (add `--visualize` for the 3D viewer page)
- `visualize`
- `telemetry`
- `benchmark`
- `assets check`

This is the package-owned command surface used by the `bin/spaceagora` wrapper.
"""
function run_cli(args::Vector{String}=copy(ARGS); io::IO=stdout, errio::IO=stderr)::Int
    isempty(args) && return _print_usage(io)
    cmd = first(args)
    tail = args[2:end]
    if cmd in ("help", "--help", "-h")
        return _print_usage(io)
    elseif cmd == "assets"
        isempty(tail) && throw(ArgumentError("assets requires a subcommand: check, manifest, or setup-open"))
        subcmd = first(tail)
        if subcmd == "check"
            render_asset_report(check_assets(); io=io)
            return 0
        elseif subcmd == "manifest"
            render_asset_manifest(load_asset_manifest(); io=io)
            return 0
        elseif subcmd == "setup-open"
            setup_open_assets(; io=io)
            return 0
        end
        throw(ArgumentError("assets supports only: check, manifest, or setup-open"))
    elseif cmd == "run"
        return _run_example(tail; io=io, errio=errio)
    elseif cmd == "visualize"
        return _run_visualize(tail; io=io, errio=errio)
    elseif cmd == "telemetry"
        return _run_telemetry(tail; io=io, errio=errio)
    elseif cmd == "benchmark"
        return _run_benchmark(tail; io=io, errio=errio)
    end
    throw(ArgumentError("Unknown SpaceAGORA CLI command '$cmd'."))
end

end # module SpaceAGORACLI
