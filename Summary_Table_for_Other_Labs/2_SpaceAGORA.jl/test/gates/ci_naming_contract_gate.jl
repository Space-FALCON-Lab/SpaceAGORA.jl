const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const CANONICAL_DIRS = [
    joinpath(REPO_ROOT, "src", "simulation", "engine"),
    joinpath(REPO_ROOT, "src", "simulation", "callbacks"),
    joinpath(REPO_ROOT, "src", "parallel", "routing"),
    joinpath(REPO_ROOT, "benchmarks", "studies"),
    joinpath(REPO_ROOT, "benchmarks", "scripts")
]

const FILE_RE = r"^[a-z0-9]+(_[a-z0-9]+)*\.jl$"
const MODULE_NAMING_TARGETS = [
    joinpath(REPO_ROOT, "src", "core", "state", "reference_system_config.jl"),
    joinpath(REPO_ROOT, "src", "vehicle", "spacecraft", "model.jl"),
]

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

for path in MODULE_NAMING_TARGETS
    isfile(path) || error("Missing module naming target: $(relpath(path, REPO_ROOT))")
    src = read(path, String)
    occursin("module ref_sys", src) &&
        push!(violations, "$(relpath(path, REPO_ROOT)) contains retired snake_case module name 'ref_sys'")
    occursin("module PhysicalModel", src) &&
        push!(violations, "$(relpath(path, REPO_ROOT)) contains retired module name 'PhysicalModel'")
end

if !isempty(violations)
    println("naming_contract_violations:")
    for v in sort(violations)
        println("  - $v")
    end
    error("Naming contract gate failed")
end

println("naming_contract_gate_ok")
