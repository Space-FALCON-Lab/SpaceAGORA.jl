const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

function _active_source(path::String)
    src = read(path, String)
    return join(
        (
            strip(first(split(line, '#', limit=2)))
            for line in split(src, '\n')
            if !isempty(strip(first(split(line, '#', limit=2))))
        ),
        '\n',
    )
end

const FILE_TARGETS = [
    joinpath(REPO_ROOT, "src", "mission", "operations", "aerobraking_policy", "selector_stub.jl"),
    joinpath(REPO_ROOT, "src", "gnc", "internal", "bridge_helpers.jl"),
    joinpath(REPO_ROOT, "src", "simulation", "solver_orchestration", "implicit_midpoint_jacobian.jl"),
]

const DIR_TARGETS = [
    joinpath(REPO_ROOT, "src", "gnc", "control", "aerobraking"),
    joinpath(REPO_ROOT, "src", "gnc", "guidance", "aerobraking", "t_edg"),
]

const FORBIDDEN_REGEXES = (
    r"\b(?:input\.)?args\s+isa\s+AbstractDict\b",
    r"\bgravity_const\s*\(",
    r"\bgravity_invsquared\s*\(",
    r"\bgravity_invsquared_J2\s*\(",
    r"\bgravity_GRAM\s*\(",
)

violations = String[]
targets = String[]
append!(targets, FILE_TARGETS)
for dir in DIR_TARGETS
    isdir(dir) || continue
    for (root, _, files) in walkdir(dir)
        for file in files
            endswith(file, ".jl") || continue
            push!(targets, joinpath(root, file))
        end
    end
end

for path in sort(unique(targets))
    isfile(path) || error("Missing runtime legacy-surface target: $(relpath(path, REPO_ROOT))")
    active_src = _active_source(path)
    rel = relpath(path, REPO_ROOT)
    for rx in FORBIDDEN_REGEXES
        occursin(rx, active_src) || continue
        push!(violations, "$rel: contains forbidden runtime legacy token '$rx'")
    end
end

if !isempty(violations)
    println("runtime_legacy_surface_violations:")
    for violation in sort(unique(violations))
        println("  - $violation")
    end
    error("Runtime legacy surface gate failed")
end

println("runtime_legacy_surface_gate_ok")
