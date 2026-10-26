include(joinpath(@__DIR__, "..", "contracts", "edg_ownership_checks.jl"))
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

isdir(joinpath(REPO_ROOT, "src", "models", "gnc", "guidance", "targeting_control")) &&
    error("Retired guidance targeting_control path still exists.")

violations = String[]

# Guidance must not include control source files.
guidance_root = joinpath(REPO_ROOT, "src", "models", "gnc", "guidance", "aerobraking")
for (root, _, files) in walkdir(guidance_root)
    for file in files
        endswith(file, ".jl") || continue
        path = joinpath(root, file)
        rel = relpath(path, REPO_ROOT)
        src = read(path, String)
        EDGOwnershipChecks.has_control_include(src) && push!(violations, "$rel: guidance includes control source directly")
    end
end

for root_rel in (
    joinpath("src", "models", "gnc", "guidance"),
    joinpath("src", "models", "gnc", "planning"),
    joinpath("src", "models", "gnc", "companions"),
    joinpath("src", "models", "gnc", "contracts"),
    joinpath("src", "interfaces", "analysis", "planning"),
    joinpath("src", "models", "gnc", "control"),
)
    root = joinpath(REPO_ROOT, root_rel)
    isdir(root) || error("Missing GNC boundary scan root: $root_rel")
    for (dir, _, files) in walkdir(root)
        for file in files
            endswith(file, ".jl") || continue
            path = joinpath(dir, file)
            rel = relpath(path, REPO_ROOT)
            src = read(path, String)
            occursin("DynamicEffectors", src) && push!(violations, "$rel: GNC source still depends on DynamicEffectors directly")
        end
    end
end

append!(violations, EDGOwnershipChecks.violations(EDGOwnershipChecks.source_map(REPO_ROOT)))

if !isempty(violations)
    println("gnc_aerobraking_boundary_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("GNC aerobraking boundary gate failed")
end

println("gnc_aerobraking_boundary_gate_ok")
