const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const GUIDANCE_ROOTS = (joinpath(REPO_ROOT, "packages", "SpaceAGORAHYPR", "src"),
    (joinpath(REPO_ROOT, "src", "gnc", owner) for owner in ("guidance", "shared", "hypr", "rrt"))...)

violations = String[]
for owner_root in GUIDANCE_ROOTS, (root, _, files) in walkdir(owner_root)
    for file in files
        endswith(file, ".jl") || continue
        path = joinpath(root, file)
        rel = relpath(path, REPO_ROOT)
        src = read(path, String)
        for line in split(src, '\n')
            stripped = strip(first(split(line, '#', limit=2)))
            isempty(stripped) && continue
            if startswith(stripped, "include(") && occursin("control", stripped)
                push!(violations, "$rel: $stripped")
            end
        end
    end
end

if !isempty(violations)
    println("guidance_control_cross_include_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Guidance/control cross-include gate failed")
end

println("no_guidance_control_cross_include_gate_ok")
