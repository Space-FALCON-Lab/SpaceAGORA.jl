const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

Base.include(@__MODULE__, joinpath(@__DIR__, "runtime_any_hotpath_checks.jl"))

const HOTPATH_FILES = String[
    joinpath(REPO_ROOT, "src", "orchestration", "simulation", "state", "runtime_types.jl"),
    joinpath(REPO_ROOT, "src", "orchestration", "simulation", "engine", "setup.jl"),
    joinpath(REPO_ROOT, "src", "models", "dynamics", "effectors", "aerodynamic_wrench_models.jl"),
    joinpath(REPO_ROOT, "src", "models", "dynamics", "effectors", "perturbations.jl"),
]

# Split source keeps exactly the original hot-path scan coverage. New children
# are covered recursively so later source splits cannot bypass the contract.
const SPLIT_HOTPATH_DIR = joinpath(REPO_ROOT, "src", "models", "dynamics", "effectors", "perturbations")
const MULTIBODY_POLICY_FILE = joinpath(REPO_ROOT, "src", "orchestration", "parallel", "policy", "multibody_routing.jl")

const CALLBACK_DIR = joinpath(REPO_ROOT, "src", "orchestration", "simulation", "callbacks")

targets = [HOTPATH_FILES; MULTIBODY_POLICY_FILE]

for scan_root in (CALLBACK_DIR, SPLIT_HOTPATH_DIR)
    isdir(scan_root) || error("Missing runtime hot-path scan root: $(replace(relpath(scan_root, REPO_ROOT), '\\' => '/'))")
    for (root, _, files) in walkdir(scan_root)
        for file in files
            endswith(file, ".jl") || continue
            push!(targets, joinpath(root, file))
        end
    end
end

violations = RuntimeAnyHotpathChecks.scan(REPO_ROOT, targets, CALLBACK_DIR;
    runtime_types_path=first(HOTPATH_FILES))

if !isempty(violations)
    println("runtime_any_hotpath_violations:")
    for violation in sort(unique(violations))
        println("  - $violation")
    end
    error("Runtime Any hot-path gate failed")
end

println("runtime_any_hotpath_gate_ok")
