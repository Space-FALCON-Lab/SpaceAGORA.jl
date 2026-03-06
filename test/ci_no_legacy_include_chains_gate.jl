const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const SCAN_ROOTS = (
    joinpath(REPO_ROOT, "src", "analysis", "reports"),
    joinpath(REPO_ROOT, "src", "core"),
    joinpath(REPO_ROOT, "src", "environment"),
    joinpath(REPO_ROOT, "src", "gnc", "control"),
    joinpath(REPO_ROOT, "src", "gnc", "guidance", "targeting_control"),
    joinpath(REPO_ROOT, "src", "vehicle", "resources"),
    joinpath(REPO_ROOT, "src", "vehicle", "actuators", "laser_terminal"),
    joinpath(REPO_ROOT, "src", "mission", "constellation", "network"),
    joinpath(REPO_ROOT, "src", "gnc", "estimation"),
    joinpath(REPO_ROOT, "src", "simulation", "engine"),
    joinpath(REPO_ROOT, "src", "simulation", "callbacks"),
    joinpath(REPO_ROOT, "src", "simulation", "solver_orchestration"),
    joinpath(REPO_ROOT, "src", "parallel")
)

const FORBIDDEN_TOKENS = (
    "__legacy_",
)

const FORBIDDEN_REGEXES = (
    r"include\(joinpath\(@__DIR__,\s*\"\.\.\",\s*\"\.\.\",\s*\"\.\.\",",
)

const ALLOWED_RAW_INCLUDE_FILES = Set([
    joinpath("src", "core", "simulation_model.jl"),
    joinpath("src", "environment", "physical_models.jl"),
    joinpath("src", "environment", "ephemerides", "planet_data.jl"),
    joinpath("src", "environment", "ephemerides", "planets.jl"),
    joinpath("src", "vehicle", "resources", "resources.jl"),
    joinpath("src", "gnc", "control", "control.jl"),
    joinpath("src", "gnc", "control", "eoms.jl"),
    joinpath("src", "gnc", "control", "eom_ctrl.jl"),
    joinpath("src", "gnc", "control", "effectors.jl"),
    joinpath("src", "gnc", "control", "propulsive_maneuvers.jl"),
    joinpath("src", "gnc", "control", "closed_form_solution.jl"),
    joinpath("src", "gnc", "control", "legacy_include_helpers.jl"),
    joinpath("src", "gnc", "control", "heatload_control", "utils_timeswitch.jl"),
    joinpath("src", "gnc", "control", "heatload_control", "time_switch_calcs.jl"),
    joinpath("src", "gnc", "control", "heatload_control", "second_tsw_calcs.jl"),
    joinpath("src", "gnc", "control", "heatload_control", "security_mode.jl"),
    joinpath("src", "gnc", "guidance", "targeting_control", "sim_targeting.jl"),
    joinpath("src", "gnc", "guidance", "targeting_control", "eom_targeting.jl"),
    joinpath("src", "gnc", "guidance", "targeting_control", "targeting.jl"),
    joinpath("src", "simulation", "engine", "simulation_engine.jl"),
    joinpath("src", "simulation", "engine", "setup.jl"),
    joinpath("src", "simulation", "callbacks", "callbacks.jl"),
    joinpath("src", "simulation", "callbacks", "registry.jl"),
    joinpath("src", "simulation", "solver_orchestration", "integrators.jl"),
    joinpath("src", "simulation", "solver_orchestration", "implicit_midpoint_jacobian.jl"),
    joinpath("src", "simulation", "solver_orchestration", "include_helpers.jl"),
    joinpath("src", "parallel", "routing", "parallel_profiles.jl"),
])

violations = String[]
for root in SCAN_ROOTS
    isdir(root) || continue
    for (dir, _, files) in walkdir(root)
        for file in files
            endswith(file, ".jl") || continue
            path = joinpath(dir, file)
            rel = relpath(path, REPO_ROOT)
            src = read(path, String)
            active_src = join(
                (strip(first(split(line, '#', limit=2))) for line in split(src, '\n') if !isempty(strip(first(split(line, '#', limit=2))))),
                '\n'
            )

            for token in FORBIDDEN_TOKENS
                occursin(token, active_src) || continue
                push!(violations, "$rel: contains forbidden legacy token '$token'.")
            end
            for rx in FORBIDDEN_REGEXES
                occursin(rx, active_src) || continue
                push!(violations, "$rel: contains forbidden cross-layer include chain pattern '$rx'.")
            end

            has_raw_include = any(occursin(r"^\s*include\(", line) for line in split(active_src, '\n'))
            if has_raw_include && !(rel in ALLOWED_RAW_INCLUDE_FILES)
                push!(violations, "$rel: raw include(...) is not permitted here; move include orchestration into an approved helper/aggregator file.")
            end
        end
    end
end

if !isempty(violations)
    println("no_legacy_include_chains_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Legacy include-chain gate failed")
end

println("no_legacy_include_chains_gate_ok")
