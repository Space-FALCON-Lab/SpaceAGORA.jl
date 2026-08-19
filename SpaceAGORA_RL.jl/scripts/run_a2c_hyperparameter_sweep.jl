#!/usr/bin/env julia

using TOML

const PACKAGE_DIR = normpath(joinpath(@__DIR__, ".."))
const REPOSITORY_DIR = dirname(PACKAGE_DIR)
const BASE_CONFIG = joinpath(
    PACKAGE_DIR,
    "configs",
    "aerobraking",
    "a2c_spaceagora_physics_marsgram.toml",
)
const SWEEP_ROOT = joinpath(
    PACKAGE_DIR,
    "outputs",
    "hyperparameter_sweeps",
    "a2c_one_factor_20260814",
)
const CONFIG_DIR = joinpath(SWEEP_ROOT, "configs")
const RUNS_DIR = joinpath(SWEEP_ROOT, "runs")
const TRAIN_SCRIPT = joinpath(PACKAGE_DIR, "scripts", "train.jl")

const NOMINAL = (
    entropy_coef=0.01,
    successful_case_repetitions=9,
    learning_rate=1.0e-4,
    segment_length=10,
)

const VARIANTS = [
    (
        name="nominal",
        label="Nominal",
        parameter="nominal",
        changed_value="none",
        values=NOMINAL,
    ),
    (
        name="entropy_0p02",
        label="Changed (entropy coefficient = 0.02)",
        parameter="entropy_coefficient",
        changed_value="0.02",
        values=merge(NOMINAL, (entropy_coef=0.02,)),
    ),
    (
        name="successful_repetitions_3",
        label="Changed (successful case repetitions = 3)",
        parameter="successful_case_repetitions",
        changed_value="3",
        values=merge(NOMINAL, (successful_case_repetitions=3,)),
    ),
    (
        name="learning_rate_7p5e5",
        label="Changed (learning rate = 7.5e-5)",
        parameter="learning_rate",
        changed_value="7.5e-5",
        values=merge(NOMINAL, (learning_rate=7.5e-5,)),
    ),
    (
        name="rollout_length_16",
        label="Changed (rollout length = 16)",
        parameter="rollout_length",
        changed_value="16",
        values=merge(NOMINAL, (segment_length=16,)),
    ),
]

function _usage(io::IO=stdout)
    println(io, """
Usage:
  julia --project=SpaceAGORA_RL.jl \\
    SpaceAGORA_RL.jl/scripts/run_a2c_hyperparameter_sweep.jl [--generate-only]

Generates the five controlled A2C configurations and, unless --generate-only
is supplied, trains and checkpoint-validates each incomplete sweep arm.
""")
end

function _write_config(variant)
    config = TOML.parsefile(BASE_CONFIG)

    config["a2c"]["entropy_coef"] = variant.values.entropy_coef
    config["a2c"]["learning_rate"] = variant.values.learning_rate
    config["a2c"]["segment_length"] = variant.values.segment_length

    config["training"]["successful_case_repetitions"] =
        variant.values.successful_case_repetitions
    config["training"]["global_steps"] = 450_000
    config["training"]["validation_episodes"] = 20
    config["training"]["output_dir"] = RUNS_DIR
    config["reports"]["output_dir"] = joinpath(SWEEP_ROOT, "reports")

    path = joinpath(CONFIG_DIR, "a2c_sweep_$(variant.name).toml")
    open(path, "w") do io
        TOML.print(io, config; sorted=true)
    end
    return path
end

function _write_definition(config_paths)
    definition = Dict(
        "base_config" => BASE_CONFIG,
        "design" => "one_factor_at_a_time",
        "global_steps" => 450_000,
        "validation_episodes" => 20,
        "validation_seed" => 1,
        "training_seed" => 42,
        "variants" => [
            Dict(
                "name" => variant.name,
                "label" => variant.label,
                "parameter" => variant.parameter,
                "changed_value" => variant.changed_value,
                "config_path" => config_paths[index],
                "entropy_coef" => variant.values.entropy_coef,
                "successful_case_repetitions" =>
                    variant.values.successful_case_repetitions,
                "learning_rate" => variant.values.learning_rate,
                "segment_length" => variant.values.segment_length,
            )
            for (index, variant) in enumerate(VARIANTS)
        ],
    )
    open(joinpath(SWEEP_ROOT, "sweep_definition.toml"), "w") do io
        TOML.print(io, definition; sorted=true)
    end
end

function _matching_run_dirs(config_path)
    isdir(RUNS_DIR) || return String[]
    matches = String[]
    for entry in readdir(RUNS_DIR; join=true)
        manifest_path = joinpath(entry, "manifest.toml")
        isfile(manifest_path) || continue
        manifest = TOML.parsefile(manifest_path)
        manifest_config = get(manifest, "config_path", "")
        try
            same_config = samefile(manifest_config, config_path)
            same_config && push!(matches, entry)
        catch
            abspath(manifest_config) == abspath(config_path) && push!(matches, entry)
        end
    end
    return sort(matches)
end

function _completed_run(config_path)
    for run_dir in reverse(_matching_run_dirs(config_path))
        best_path = joinpath(run_dir, "checkpoint_validation", "best_validation_checkpoint.txt")
        isfile(best_path) && return run_dir
    end
    return nothing
end

function main(args=ARGS)
    any(==("--help"), args) && return _usage()
    all(arg -> arg == "--generate-only", args) ||
        throw(ArgumentError("unknown option; use --help for usage"))
    generate_only = any(==("--generate-only"), args)

    mkpath(CONFIG_DIR)
    mkpath(RUNS_DIR)
    config_paths = [_write_config(variant) for variant in VARIANTS]
    _write_definition(config_paths)
    println("wrote sweep definition and configurations to ", SWEEP_ROOT)
    generate_only && return nothing

    for (variant, config_path) in zip(VARIANTS, config_paths)
        completed = _completed_run(config_path)
        if completed !== nothing
            println("skipping completed variant ", variant.name, ": ", completed)
            continue
        end
        println("starting sweep variant ", variant.name, ": ", variant.label)
        command = `$(Base.julia_cmd()) --project=$(PACKAGE_DIR) $(TRAIN_SCRIPT) $(config_path)`
        cd(REPOSITORY_DIR) do
            run(setenv(command, "JULIA_NUM_THREADS" => "1"))
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main()
    catch error
        showerror(stderr, error, catch_backtrace())
        println(stderr)
        exit(1)
    end
end
