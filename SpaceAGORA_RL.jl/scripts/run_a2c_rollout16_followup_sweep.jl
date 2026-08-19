#!/usr/bin/env julia

using TOML

const PACKAGE_DIR = normpath(joinpath(@__DIR__, ".."))
const REPOSITORY_DIR = dirname(PACKAGE_DIR)
const BASELINE_SWEEP_ROOT = joinpath(
    PACKAGE_DIR,
    "outputs",
    "hyperparameter_sweeps",
    "a2c_one_factor_20260814",
)
const BASE_CONFIG = joinpath(
    BASELINE_SWEEP_ROOT,
    "configs",
    "a2c_sweep_rollout_length_16.toml",
)
const BASELINE_RUNS_DIR = joinpath(BASELINE_SWEEP_ROOT, "runs")
const SWEEP_ROOT = joinpath(
    PACKAGE_DIR,
    "outputs",
    "hyperparameter_sweeps",
    "a2c_rollout16_followup_20260816",
)
const CONFIG_DIR = joinpath(SWEEP_ROOT, "configs")
const RUNS_DIR = joinpath(SWEEP_ROOT, "runs")
const TRAIN_SCRIPT = joinpath(PACKAGE_DIR, "scripts", "train.jl")

const BASELINE = (
    entropy_coef=0.01,
    successful_case_repetitions=9,
    learning_rate=1.0e-4,
    segment_length=16,
)

const VARIANTS = [
    (
        name="rollout16_entropy_0p0075",
        label="Changed (rollout 16, entropy coefficient = 0.0075)",
        parameter="entropy_coefficient",
        nominal_value="0.01",
        changed_value="0.0075",
        values=merge(BASELINE, (entropy_coef=0.0075,)),
    ),
    (
        name="rollout16_learning_rate_1p25e4",
        label="Changed (rollout 16, learning rate = 1.25e-4)",
        parameter="learning_rate",
        nominal_value="1e-4",
        changed_value="1.25e-4",
        values=merge(BASELINE, (learning_rate=1.25e-4,)),
    ),
    (
        name="rollout16_successful_repetitions_12",
        label="Changed (rollout 16, successful case repetitions = 12)",
        parameter="successful_case_repetitions",
        nominal_value="9",
        changed_value="12",
        values=merge(BASELINE, (successful_case_repetitions=12,)),
    ),
    (
        name="rollout_length_20",
        label="Changed (rollout length = 20)",
        parameter="rollout_length",
        nominal_value="16",
        changed_value="20",
        values=merge(BASELINE, (segment_length=20,)),
    ),
]

function _usage(io::IO=stdout)
    println(io, """
Usage:
  julia --project=SpaceAGORA_RL.jl \\
    SpaceAGORA_RL.jl/scripts/run_a2c_rollout16_followup_sweep.jl [--generate-only]

Reuses the completed rollout-length-16 run from a2c_one_factor_20260814 as the
baseline. Generates four controlled follow-up configurations and, unless
--generate-only is supplied, trains and checkpoint-validates each incomplete
arm. Re-running the command skips arms with completed best-checkpoint validation.
""")
end

function _matching_run_dirs(runs_dir, config_path)
    isdir(runs_dir) || return String[]
    matches = String[]
    for entry in readdir(runs_dir; join=true)
        manifest_path = joinpath(entry, "manifest.toml")
        isfile(manifest_path) || continue
        manifest = TOML.parsefile(manifest_path)
        manifest_config = String(get(manifest, "config_path", ""))
        try
            samefile(manifest_config, config_path) && push!(matches, entry)
        catch
            abspath(manifest_config) == abspath(config_path) && push!(matches, entry)
        end
    end
    return sort(matches)
end

function _completed_run(runs_dir, config_path)
    for run_dir in reverse(_matching_run_dirs(runs_dir, config_path))
        best_path = joinpath(
            run_dir,
            "checkpoint_validation",
            "best_validation_checkpoint.txt",
        )
        isfile(best_path) && return run_dir
    end
    return nothing
end

function _baseline_run()
    isfile(BASE_CONFIG) || error("rollout-16 baseline config does not exist: $BASE_CONFIG")
    run_dir = _completed_run(BASELINE_RUNS_DIR, BASE_CONFIG)
    run_dir === nothing && error(
        "no completed rollout-16 baseline run found under $BASELINE_RUNS_DIR",
    )
    return run_dir
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

    path = joinpath(CONFIG_DIR, "a2c_followup_$(variant.name).toml")
    open(path, "w") do io
        TOML.print(io, config; sorted=true)
    end
    return path
end

function _write_definition(config_paths, baseline_run)
    baseline_record = Dict(
        "name" => "rollout_length_16_baseline",
        "label" => "Baseline (rollout length = 16)",
        "parameter" => "baseline",
        "changed_value" => "none",
        "config_path" => BASE_CONFIG,
        "run_dir" => baseline_run,
        "entropy_coef" => BASELINE.entropy_coef,
        "successful_case_repetitions" => BASELINE.successful_case_repetitions,
        "learning_rate" => BASELINE.learning_rate,
        "segment_length" => BASELINE.segment_length,
    )
    changed_records = [
        Dict(
            "name" => variant.name,
            "label" => variant.label,
            "parameter" => variant.parameter,
            "nominal_value" => variant.nominal_value,
            "changed_value" => variant.changed_value,
            "config_path" => config_paths[index],
            "entropy_coef" => variant.values.entropy_coef,
            "successful_case_repetitions" =>
                variant.values.successful_case_repetitions,
            "learning_rate" => variant.values.learning_rate,
            "segment_length" => variant.values.segment_length,
        )
        for (index, variant) in enumerate(VARIANTS)
    ]
    definition = Dict(
        "base_config" => BASE_CONFIG,
        "design" => "one_factor_at_a_time_from_rollout_16_baseline",
        "baseline_variant" => "rollout_length_16_baseline",
        "comparison_variants" => getfield.(VARIANTS, :name),
        "global_steps" => 450_000,
        "validation_episodes" => 20,
        "validation_seed" => 1,
        "training_seed" => 42,
        "variants" => vcat([baseline_record], changed_records),
    )
    open(joinpath(SWEEP_ROOT, "sweep_definition.toml"), "w") do io
        TOML.print(io, definition; sorted=true)
    end
end

function main(args=ARGS)
    any(==("--help"), args) && return _usage()
    all(arg -> arg == "--generate-only", args) ||
        throw(ArgumentError("unknown option; use --help for usage"))
    generate_only = any(==("--generate-only"), args)

    baseline_run = _baseline_run()
    mkpath(CONFIG_DIR)
    mkpath(RUNS_DIR)
    config_paths = [_write_config(variant) for variant in VARIANTS]
    _write_definition(config_paths, baseline_run)
    println("reusing rollout-16 baseline run: ", baseline_run)
    println("wrote sweep definition and configurations to ", SWEEP_ROOT)
    generate_only && return nothing

    for (variant, config_path) in zip(VARIANTS, config_paths)
        completed = _completed_run(RUNS_DIR, config_path)
        if completed !== nothing
            println("skipping completed variant ", variant.name, ": ", completed)
            continue
        end
        println("starting follow-up variant ", variant.name, ": ", variant.label)
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
