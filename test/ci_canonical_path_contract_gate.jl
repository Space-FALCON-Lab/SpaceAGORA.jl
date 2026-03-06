const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const SCAN_ROOTS = (
    joinpath(REPO_ROOT, "test"),
    joinpath(REPO_ROOT, "docs"),
    joinpath(REPO_ROOT, ".github", "workflows")
)

const FORBIDDEN_PATH_TOKENS = (
    "src/control/",
    "src/physical_models/",
    "src/guidance/",
    "src/integrator/",
    "src/utils/",
    "src/simulation_model/",
    "src/simulation/execution/simulation_elements.jl",
    "src/simulation/execution/simulation_execution.jl",
    "src/simulation/execution/run_simulation.jl"
)

@inline function _should_scan(path::String)::Bool
    endswith(path, ".jl") && return true
    endswith(path, ".md") && return true
    endswith(path, ".yml") && return true
    endswith(path, ".yaml") && return true
    return false
end

@inline function _is_excluded(path::String)::Bool
    rel = relpath(path, REPO_ROOT)
    # Historical AI notes are archival and excluded from canonical-path enforcement.
    startswith(rel, joinpath("test", "ai_")) && return true
    startswith(rel, joinpath("test", "ci_")) && endswith(rel, "_gate.jl") && return true
    return false
end

violations = String[]
for root in SCAN_ROOTS
    isdir(root) || continue
    for (dir, _, files) in walkdir(root)
        for file in files
            path = joinpath(dir, file)
            _should_scan(path) || continue
            _is_excluded(path) && continue
            rel = relpath(path, REPO_ROOT)
            src = read(path, String)
            for tok in FORBIDDEN_PATH_TOKENS
                occursin(tok, src) || continue
                push!(violations, "$(rel): contains forbidden canonical-path token '$tok'")
            end
        end
    end
end

if !isempty(violations)
    println("canonical_path_contract_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Canonical path contract gate failed")
end

println("canonical_path_contract_gate_ok")
