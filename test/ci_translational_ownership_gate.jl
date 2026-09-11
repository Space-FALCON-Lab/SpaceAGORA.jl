const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const SRC_ROOT = joinpath(REPO_ROOT, "src")
const TRANSLATIONAL_ROOT = joinpath(SRC_ROOT, "dynamics", "translational")

# Translational equations must be owned by src/dynamics/translational/* only.
const FORBIDDEN_PATTERNS = (
    r"du_view\.pos\s*\.=\s*sc_view\.vel",
    r"du_view\.pos\s*\.=\s*0\.0",
    r"du_view\.vel\s*\.=\s*forces\s*/\s*sc_view\.mass",
    r"du_view\.vel\s*\.=\s*net_force\s*/\s*sc_view\.mass",
    r"du_view\.mass\s*=\s*mass_rate",
    r"du_view\.mass\s*=\s*0\.0",
)

violations = String[]

for (root, _, files) in walkdir(SRC_ROOT)
    for file in files
        endswith(file, ".jl") || continue
        path = joinpath(root, file)
        startswith(path, TRANSLATIONAL_ROOT) && continue

        rel = relpath(path, REPO_ROOT)
        src = read(path, String)
        active_src = join(
            (strip(first(split(line, '#', limit=2))) for line in split(src, '\n') if !isempty(strip(first(split(line, '#', limit=2))))),
            '\n'
        )

        for rx in FORBIDDEN_PATTERNS
            occursin(rx, active_src) || continue
            push!(violations, "$rel: contains forbidden inline translational dynamics pattern '$rx'.")
        end
    end
end

if !isempty(violations)
    println("translational_ownership_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Translational ownership gate failed")
end

println("translational_ownership_gate_ok")
