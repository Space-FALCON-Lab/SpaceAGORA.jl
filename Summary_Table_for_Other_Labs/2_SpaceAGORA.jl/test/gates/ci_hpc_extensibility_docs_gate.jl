const REPO_ROOT = dirname(dirname(dirname(@__FILE__)))

required_paths = [
    joinpath(REPO_ROOT, "docs", "src", "public_api_policy.md"),
    joinpath(REPO_ROOT, "docs", "src", "distributed_hpc.md"),
    joinpath(REPO_ROOT, "docs", "src", "extensibility.md"),
    joinpath(REPO_ROOT, "templates", "force_torque_model_template.jl"),
    joinpath(REPO_ROOT, "templates", "density_model_template.jl"),
    joinpath(REPO_ROOT, "templates", "control_hook_template.jl"),
    joinpath(REPO_ROOT, "scripts", "hpc", "local_process_runtime_analysis.sh"),
    joinpath(REPO_ROOT, "scripts", "hpc", "slurm_process_runtime_analysis.sh"),
]

missing = filter(path -> !isfile(path), required_paths)
isempty(missing) || error("Missing Wave 3/Wave 5 documentation assets: $(missing)")

hpc_doc = read(joinpath(REPO_ROOT, "docs", "src", "distributed_hpc.md"), String)
ext_doc = read(joinpath(REPO_ROOT, "docs", "src", "extensibility.md"), String)
api_doc = read(joinpath(REPO_ROOT, "docs", "src", "public_api_policy.md"), String)

for marker in (
    "SPACEAGORA_PERF_PROCS",
    "SPACEAGORA_PERF_WORKER_PROJECT",
    "SPACEAGORA_PARALLEL_POLICY_STATE_PATH",
    "shared depot",
    "scripts/hpc/slurm_process_runtime_analysis.sh",
)
    occursin(marker, hpc_doc) || error("distributed_hpc.md missing marker: $(marker)")
end

for marker in (
    "calcForceTorque",
    "getDensity",
    "calcControlEffect!",
    "templates/force_torque_model_template.jl",
    "templates/control_hook_template.jl",
)
    occursin(marker, ext_doc) || error("extensibility.md missing marker: $(marker)")
end

for marker in (
    "root `SpaceAGORA` module",
    "Public API",
    "internal",
    "generated public API",
)
    occursin(marker, api_doc) || error("public_api_policy.md missing marker: $(marker)")
end

println("hpc_extensibility_docs_gate_ok")
