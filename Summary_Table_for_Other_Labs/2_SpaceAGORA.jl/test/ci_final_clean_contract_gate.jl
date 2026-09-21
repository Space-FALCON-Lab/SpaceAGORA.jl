const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const SRC_ROOT = joinpath(REPO_ROOT, "src")

violations = String[]

# 1) No `legacy_` in src
for (root, _, files) in walkdir(SRC_ROOT)
    for file in files
        endswith(file, ".jl") || continue
        path = joinpath(root, file)
        rel = relpath(path, REPO_ROOT)
        src = read(path, String)
        occursin("legacy_", src) && push!(violations, "$rel: contains forbidden token legacy_")
    end
end

# 2) Retired GNC effectors hook file paths must not exist
for rel in (
    joinpath("src", "gnc", "control", "effectors.jl"),
    joinpath("src", "gnc", "guidance", "effectors.jl"),
    joinpath("src", "gnc", "navigation", "effectors.jl"),
    joinpath("src", "gnc", "control", "bridge_helpers.jl"),
    joinpath("src", "gnc", "control", "bridge_helpers.jl"),
    joinpath("src", "gnc", "guidance", "targeting_control"),
    joinpath("src", "gnc", "control", "heatload_control"),
)
    ispath(joinpath(REPO_ROOT, rel)) && push!(violations, "$rel: retired file or directory still exists")
end

# 3) Retired hook module names must not appear in src
for (root, _, files) in walkdir(SRC_ROOT)
    for file in files
        endswith(file, ".jl") || continue
        path = joinpath(root, file)
        rel = relpath(path, REPO_ROOT)
        src = read(path, String)
        for old_name in ("ControlEffectors", "GuidanceEffectors", "NavigationEffectors")
            occursin(old_name, src) && push!(violations, "$rel: contains retired module name $old_name")
        end
    end
end

if !isempty(violations)
    println("final_clean_contract_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Final clean contract gate failed")
end

println("final_clean_contract_gate_ok")
