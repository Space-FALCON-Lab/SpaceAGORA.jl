const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const SRC_ROOT = joinpath(REPO_ROOT, "src")

const RETIRED_DIRS = (
    joinpath(SRC_ROOT, "control"),
    joinpath(SRC_ROOT, "physical_models"),
    joinpath(SRC_ROOT, "guidance"),
    joinpath(SRC_ROOT, "integrator"),
    joinpath(SRC_ROOT, "utils"),
    joinpath(SRC_ROOT, "simulation_model")
)

for dir in RETIRED_DIRS
    isdir(dir) && error("Retired legacy source tree still exists: $(relpath(dir, REPO_ROOT))")
end

const FORBIDDEN_SRC_TOKENS = (
    "src/control/",
    "src/physical_models/",
    "src/guidance/",
    "src/integrator/",
    "src/utils/",
    "src/simulation_model/",
    "Compatibility wrapper: canonical path forwarding to legacy implementation."
)

violations = String[]
for (root, _, files) in walkdir(SRC_ROOT)
    for file in files
        endswith(file, ".jl") || continue
        path = joinpath(root, file)
        rel = relpath(path, REPO_ROOT)
        src = read(path, String)
        for tok in FORBIDDEN_SRC_TOKENS
            occursin(tok, src) || continue
            push!(violations, "$(rel): contains forbidden token '$tok'")
        end
    end
end

if !isempty(violations)
    println("legacy_ownership_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Legacy ownership gate failed")
end

println("no_legacy_ownership_gate_ok")
