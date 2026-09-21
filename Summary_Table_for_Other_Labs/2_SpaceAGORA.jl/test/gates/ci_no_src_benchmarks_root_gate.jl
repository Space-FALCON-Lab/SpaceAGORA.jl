const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SRC_BENCHMARKS_ROOT = joinpath(REPO_ROOT, "src", "benchmarks")

isdir(SRC_BENCHMARKS_ROOT) &&
    error("Canonical topology violation: src/benchmarks must not exist. Use top-level benchmarks/ instead.")

const SCAN_ROOTS = (
    joinpath(REPO_ROOT, "src"),
    joinpath(REPO_ROOT, "test"),
)

violations = String[]
for root in SCAN_ROOTS
    isdir(root) || continue
    for (dir, _, files) in walkdir(root)
        for file in files
            path = joinpath(dir, file)
            (endswith(path, ".jl") || endswith(path, ".md") || endswith(path, ".yml") || endswith(path, ".yaml")) || continue
            rel = relpath(path, REPO_ROOT)
            startswith(rel, joinpath("test", "gates", "ci_")) && endswith(rel, "_gate.jl") && continue
            src = read(path, String)
            occursin("src/benchmarks/", src) && push!(violations, "$rel: contains forbidden token 'src/benchmarks/'")
        end
    end
end

if !isempty(violations)
    println("no_src_benchmarks_root_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("No-src-benchmarks-root gate failed")
end

println("no_src_benchmarks_root_gate_ok")
