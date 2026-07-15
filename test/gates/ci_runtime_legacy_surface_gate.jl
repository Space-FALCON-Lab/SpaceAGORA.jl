const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

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
    joinpath(REPO_ROOT, "src", "gnc", "control", "aerobraking", "control_commands.jl"),
    joinpath(REPO_ROOT, "src", "gnc", "control", "aerobraking", "constraint_tracking.jl"),
    joinpath(REPO_ROOT, "src", "gnc", "control", "aerobraking", "tracking_executor.jl"),
    joinpath(REPO_ROOT, "src", "gnc", "guidance", "aerobraking", "t_edg", "trajectory_predictor.jl"),
    joinpath(REPO_ROOT, "src", "gnc", "guidance", "aerobraking", "t_edg", "eom_predictor.jl"),
    joinpath(REPO_ROOT, "src", "gnc", "guidance", "aerobraking", "t_edg", "targeting_solver.jl"),
    joinpath(REPO_ROOT, "src", "io", "config", "io_config.jl"),
    joinpath(REPO_ROOT, "src", "io", "outputs", "io_outputs.jl"),
    joinpath(REPO_ROOT, "src", "simulation", "engine", "persistence.jl"),
    joinpath(REPO_ROOT, "src", "simulation", "engine", "reporting.jl"),
    joinpath(REPO_ROOT, "src", "simulation", "engine", "setup.jl"),
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

const FILE_SPECIFIC_FORBIDDEN = Dict(
    joinpath(REPO_ROOT, "src", "mission", "operations", "aerobraking_policy", "selector_stub.jl") => (
        r"args\s*\[",
        r"_compat_",
    ),
    joinpath(REPO_ROOT, "src", "gnc", "internal", "bridge_helpers.jl") => (
        r"_compat_",
    ),
    joinpath(REPO_ROOT, "src", "gnc", "control", "aerobraking", "control_commands.jl") => (
        r"args\s*\[",
        r"param\[[0-9]+\]",
        r"_compat_",
    ),
    joinpath(REPO_ROOT, "src", "gnc", "control", "aerobraking", "constraint_tracking.jl") => (
        r"args\s*\[",
        r"param\[[0-9]+\]",
        r"_compat_",
    ),
    joinpath(REPO_ROOT, "src", "gnc", "control", "aerobraking", "tracking_executor.jl") => (
        r"args\s*\[",
        r"_compat_",
        r"monte_carlo_guidance_environment\s*\(",
    ),
    joinpath(REPO_ROOT, "src", "gnc", "guidance", "aerobraking", "t_edg", "trajectory_predictor.jl") => (
        r"args\s*\[",
        r"param\[[0-9]+\]",
        r"_compat_",
    ),
    joinpath(REPO_ROOT, "src", "gnc", "guidance", "aerobraking", "t_edg", "eom_predictor.jl") => (
        r"args\s*\[",
        r"param\[[0-9]+\]",
        r"_compat_",
    ),
    joinpath(REPO_ROOT, "src", "gnc", "guidance", "aerobraking", "t_edg", "targeting_solver.jl") => (
        r"args\s*\[",
        r"_compat_",
    ),
    joinpath(REPO_ROOT, "src", "io", "config", "io_config.jl") => (
        r"_compat_results_csv_path",
    ),
    joinpath(REPO_ROOT, "src", "io", "outputs", "io_outputs.jl") => (
        r"_write_compat_results_csv!",
    ),
    joinpath(REPO_ROOT, "src", "simulation", "engine", "persistence.jl") => (
        r"_compat_results_csv_path",
        r"_write_compat_results_csv!",
    ),
    joinpath(REPO_ROOT, "src", "simulation", "engine", "reporting.jl") => (
        r"_compat_",
    ),
    joinpath(REPO_ROOT, "src", "simulation", "engine", "setup.jl") => (
        r"_compat_",
    ),
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
    for rx in get(FILE_SPECIFIC_FORBIDDEN, path, ())
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
