const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const CANONICAL_DIRS = [
    joinpath(REPO_ROOT, "src", "simulation", "engine"),
    joinpath(REPO_ROOT, "src", "simulation", "callbacks"),
    joinpath(REPO_ROOT, "src", "parallel", "routing"),
    joinpath(REPO_ROOT, "benchmarks", "studies"),
    joinpath(REPO_ROOT, "benchmarks", "scripts")
]

const FILE_RE = r"^[a-z0-9]+(_[a-z0-9]+)*\.jl$"

violations = String[]

for dir in CANONICAL_DIRS
    isdir(dir) || continue
    for (root, _, files) in walkdir(dir)
        for file in files
            endswith(file, ".jl") || continue
            file == ".gitkeep" && continue
            if isnothing(match(FILE_RE, file))
                push!(violations, relpath(joinpath(root, file), REPO_ROOT))
            end
        end
    end
end

if !isempty(violations)
    println("naming_contract_violations:")
    for v in sort(violations)
        println("  - $v")
    end
    error("Naming contract gate failed")
end

println("naming_contract_gate_ok")
