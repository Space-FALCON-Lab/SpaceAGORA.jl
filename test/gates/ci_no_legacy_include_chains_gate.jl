const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const SCAN_ROOTS = (
    joinpath(REPO_ROOT, "src", "analysis", "reports"),
    joinpath(REPO_ROOT, "src", "core"),
    joinpath(REPO_ROOT, "src", "environment"),
    joinpath(REPO_ROOT, "src", "gnc", "control"),
    joinpath(REPO_ROOT, "src", "gnc", "guidance"),
    joinpath(REPO_ROOT, "src", "simulation", "engine"),
    joinpath(REPO_ROOT, "src", "simulation", "callbacks"),
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
    joinpath("src", "environment", "ephemerides", "ephemerides_models.jl"),
    joinpath("src", "environment", "ephemerides", "planet_data.jl"),
    joinpath("src", "environment", "ephemerides", "planets.jl"),
    joinpath("src", "gnc", "control", "control_hooks.jl"),
    joinpath("src", "gnc", "control", "propulsive_maneuvers.jl"),
    joinpath("src", "gnc", "guidance", "guidance_models.jl"),
    joinpath("src", "gnc", "guidance", "guidance_hooks.jl"),
    joinpath("src", "gnc", "navigation", "navigation_hooks.jl"),
    joinpath("src", "simulation", "engine", "simulation_engine.jl"),
    joinpath("src", "simulation", "engine", "setup.jl"),
    joinpath("src", "simulation", "callbacks", "callbacks.jl"),
    joinpath("src", "simulation", "callbacks", "registry.jl"),
    joinpath("src", "simulation", "callbacks", "density_callbacks.jl"),
    joinpath("src", "simulation", "callbacks", "gram_track_cache.jl"),
    joinpath("src", "parallel", "policy", "parallel_policy.jl"),
    joinpath("src", "parallel", "routing", "parallel_profiles.jl"),
    joinpath("src", "parallel", "process", "parallel_process.jl"),
    # Fourth src/parallel/<area>/<area>.jl module aggregator, same shape as the
    # three above: defines its module, includes its own files, exports its
    # surface, owns no behavior. Added with the cost model and omitted here by
    # oversight, so this gate has failed since that commit.
    #
    # Deliberately NOT added to CANONICAL_AGGREGATOR_FILES in
    # ci_src_completeness_contract_gate.jl, which is a different allowlist with
    # a different meaning: that one requires the
    # "# Canonical aggregator: no behavior ownership." header and forbids any
    # function or type definition. parallel_policy.jl sits in this list and not
    # in that one for the same reason.
    joinpath("src", "parallel", "cost", "parallel_cost.jl"),
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

            if startswith(rel, joinpath("src", "gnc", "guidance", "aerobraking"))
                if occursin("include(", active_src) && occursin("control", active_src)
                    push!(violations, "$rel: guidance aerobraking file includes control source directly.")
                end
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
