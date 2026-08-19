#!/usr/bin/env julia

using TOML

const PACKAGE_DIR = normpath(joinpath(@__DIR__, ".."))
const REPOSITORY_DIR = dirname(PACKAGE_DIR)
const SWEEP_ROOT = joinpath(
    PACKAGE_DIR,
    "outputs",
    "hyperparameter_sweeps",
    "a2c_one_factor_20260814",
)
const PR_DRL_RUN = joinpath(
    PACKAGE_DIR,
    "outputs",
    "runs",
    "20260811T171808_pr_drl_spaceagora_physics_marsgram-pr_drl",
)
const EVALUATION_SCRIPT = joinpath(PACKAGE_DIR, "scripts", "evaluate_rl_run.jl")
const ANALYSIS_SCRIPT = joinpath(PACKAGE_DIR, "scripts", "analyze_a2c_hyperparameter_sweep.jl")
const VARIANT_ORDER = [
    "nominal",
    "entropy_0p02",
    "successful_repetitions_3",
    "learning_rate_7p5e5",
    "rollout_length_16",
]

function _completed_run(runs_dir, config_path)
    matches = String[]
    for run_dir in readdir(runs_dir; join=true)
        manifest_path = joinpath(run_dir, "manifest.toml")
        isfile(manifest_path) || continue
        manifest = TOML.parsefile(manifest_path)
        abspath(String(get(manifest, "config_path", ""))) == abspath(config_path) || continue
        isfile(joinpath(run_dir, "checkpoint_validation", "best_validation_checkpoint.txt")) ||
            continue
        push!(matches, run_dir)
    end
    isempty(matches) && error("no completed validated run for $config_path")
    return last(sort(matches))
end

function main()
    definition = TOML.parsefile(joinpath(SWEEP_ROOT, "sweep_definition.toml"))
    variants = Dict(String(record["name"]) => record for record in definition["variants"])
    runs_dir = joinpath(SWEEP_ROOT, "runs")
    run_dirs = [
        _completed_run(runs_dir, String(variants[name]["config_path"]))
        for name in VARIANT_ORDER
    ]

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
