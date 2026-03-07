const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const MODELS_DIR = joinpath(REPO_ROOT, "src", "dynamics", "models")

isdir(MODELS_DIR) &&
    error("Canonical topology violation: src/dynamics/models must not exist.")

const SCAN_ROOTS = (
    joinpath(REPO_ROOT, "src"),
    joinpath(REPO_ROOT, "test"),
)

const FORBIDDEN_TOKENS = (
    "src/dynamics/models/",
    "\"dynamics\", \"models\"",
    "\"dynamics\",\"models\"",
)

violations = String[]
for root in SCAN_ROOTS
    isdir(root) || continue
    for (dir, _, files) in walkdir(root)
        for file in files
            path = joinpath(dir, file)
            (endswith(path, ".jl") || endswith(path, ".md") || endswith(path, ".yml") || endswith(path, ".yaml")) || continue
            rel = relpath(path, REPO_ROOT)
            startswith(rel, joinpath("test", "ci_")) && endswith(rel, "_gate.jl") && continue
            src = read(path, String)
            for tok in FORBIDDEN_TOKENS
                occursin(tok, src) || continue
                push!(violations, "$rel: contains forbidden dynamics/models token '$tok'")
            end
        end
    end
end

if !isempty(violations)
    println("no_dynamics_models_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("No-dynamics-models gate failed")
end

println("no_dynamics_models_gate_ok")
