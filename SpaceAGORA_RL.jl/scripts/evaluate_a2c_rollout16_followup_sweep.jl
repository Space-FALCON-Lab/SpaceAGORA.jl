#!/usr/bin/env julia

using TOML

const PACKAGE_DIR = normpath(joinpath(@__DIR__, ".."))
const REPOSITORY_DIR = dirname(PACKAGE_DIR)
const SWEEP_ROOT = joinpath(
    PACKAGE_DIR,
    "outputs",
    "hyperparameter_sweeps",
    "a2c_rollout16_followup_20260816",
)
const PR_DRL_RUN = joinpath(
    PACKAGE_DIR,
    "outputs",
    "runs",
    "20260811T171808_pr_drl_spaceagora_physics_marsgram-pr_drl",
)
const EVALUATION_SCRIPT = joinpath(PACKAGE_DIR, "scripts", "evaluate_rl_run.jl")
const ANALYSIS_SCRIPT = joinpath(PACKAGE_DIR, "scripts", "analyze_a2c_hyperparameter_sweep.jl")

function _usage(io::IO=stdout)
    println(io, """
Usage:
  julia --project=SpaceAGORA_RL.jl \\
    SpaceAGORA_RL.jl/scripts/evaluate_a2c_rollout16_followup_sweep.jl

Evaluates the rollout-16 baseline and all four completed follow-up arms using
their best validation checkpoints. Runs 100 shared perturbed-wind campaigns per
policy with 16 workers, compares against the same best PR-DRL checkpoint, and
writes the common evaluation plus labeled per-parameter analysis plots.
""")
end

function _completed_run(runs_dir, config_path)
    isdir(runs_dir) || error("sweep runs directory does not exist: $runs_dir")
    matches = String[]
    for run_dir in readdir(runs_dir; join=true)
        manifest_path = joinpath(run_dir, "manifest.toml")
        isfile(manifest_path) || continue
        manifest = TOML.parsefile(manifest_path)
        abspath(String(get(manifest, "config_path", ""))) == abspath(config_path) ||
            continue
        isfile(joinpath(run_dir, "checkpoint_validation", "best_validation_checkpoint.txt")) ||
            continue
        push!(matches, run_dir)
    end
    isempty(matches) && error("no completed validated run for $config_path")
    return last(sort(matches))
end

function _run_for_record(record, runs_dir)
    if haskey(record, "run_dir")
        run_dir = abspath(String(record["run_dir"]))
        isfile(joinpath(run_dir, "checkpoint_validation", "best_validation_checkpoint.txt")) ||
            error("recorded baseline run has no best validation checkpoint: $run_dir")
        return run_dir
    end
    return _completed_run(runs_dir, String(record["config_path"]))
end

function main(args=ARGS)
    any(==("--help"), args) && return _usage()
    isempty(args) || throw(ArgumentError("unknown option; use --help for usage"))

    definition_path = joinpath(SWEEP_ROOT, "sweep_definition.toml")
    isfile(definition_path) || error(
        "sweep definition does not exist; run run_a2c_rollout16_followup_sweep.jl first",
    )
    definition = TOML.parsefile(definition_path)
    variants = Dict(String(record["name"]) => record for record in definition["variants"])
    baseline_name = String(definition["baseline_variant"])
    comparison_names = String.(definition["comparison_variants"])
    variant_order = vcat([baseline_name], comparison_names)
    runs_dir = joinpath(SWEEP_ROOT, "runs")
    run_dirs = [_run_for_record(variants[name], runs_dir) for name in variant_order]

    output_dir = joinpath(SWEEP_ROOT, "common_evaluation")
    command = `$(Base.julia_cmd()) --project=$(PACKAGE_DIR) $(EVALUATION_SCRIPT) $(PR_DRL_RUN) --checkpoint best`
    for run_dir in run_dirs
        command = `$command --compare-run $run_dir`
    end
    command = `$command --wind-mode perturbed --episodes 100 --processes 16 --threads-per-process 1 --progress-every 10 --output $output_dir`
    cd(REPOSITORY_DIR) do
        run(setenv(command, "JULIA_NUM_THREADS" => "1"))
        run(`$(Base.julia_cmd()) --project=$(PACKAGE_DIR) $(ANALYSIS_SCRIPT) $(SWEEP_ROOT) $(output_dir)`)
    end
    return output_dir
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
