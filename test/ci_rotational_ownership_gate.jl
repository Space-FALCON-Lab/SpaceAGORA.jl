const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const SRC_ROOT = joinpath(REPO_ROOT, "src")
const ROTATIONAL_ROOT = joinpath(SRC_ROOT, "dynamics", "rotational")

const FORBIDDEN_PATTERNS = (
    r"0\.5\s*\*\s*quat_mult\(",
    r"cross\(\s*ω_body\s*,\s*inertia_tensor\s*\*\s*ω_body\s*\)",
    r"inertia_tensor\s*\\\s*\(\s*τ_body\s*-\s*cross\(",
)

violations = String[]

for (root, _, files) in walkdir(SRC_ROOT)
    for file in files
        endswith(file, ".jl") || continue
        path = joinpath(root, file)
        startswith(path, ROTATIONAL_ROOT) && continue

        rel = relpath(path, REPO_ROOT)
        src = read(path, String)
        active_src = join(
            (strip(first(split(line, '#', limit=2))) for line in split(src, '\n') if !isempty(strip(first(split(line, '#', limit=2))))),
            '\n'
        )

        for rx in FORBIDDEN_PATTERNS
            occursin(rx, active_src) || continue
            push!(violations, "$rel: contains forbidden inline rotational dynamics pattern '$rx'.")
        end
    end
end

if !isempty(violations)
    println("rotational_ownership_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Rotational ownership gate failed")
end

println("rotational_ownership_gate_ok")
