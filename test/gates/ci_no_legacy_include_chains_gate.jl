include(joinpath(@__DIR__, "..", "contracts", "edg_ownership_checks.jl"))
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const SCAN_ROOTS = (
    joinpath(REPO_ROOT, "src", "foundation", "core"),
    joinpath(REPO_ROOT, "src", "models", "environment"),
    joinpath(REPO_ROOT, "src", "models", "gnc", "control"),
    joinpath(REPO_ROOT, "src", "models", "gnc", "guidance"),
    joinpath(REPO_ROOT, "src", "models", "gnc", "planning"),
    joinpath(REPO_ROOT, "src", "models", "gnc", "companions"),
    joinpath(REPO_ROOT, "src", "models", "gnc", "contracts"),
    joinpath(REPO_ROOT, "src", "interfaces", "analysis", "planning"),
    joinpath(REPO_ROOT, "src", "orchestration", "simulation", "composition"),
    joinpath(REPO_ROOT, "src", "orchestration", "simulation", "state"),
    joinpath(REPO_ROOT, "src", "orchestration", "simulation", "config", "presets"),
    joinpath(REPO_ROOT, "src", "orchestration", "simulation", "engine"),
    joinpath(REPO_ROOT, "src", "orchestration", "simulation", "callbacks"),
    joinpath(REPO_ROOT, "src", "orchestration", "parallel")
)

const FORBIDDEN_TOKENS = (
    "__legacy_",
)

const FORBIDDEN_REGEXES = (
    r"include\(joinpath\(@__DIR__,\s*\"\.\.\",\s*\"\.\.\",\s*\"\.\.\",",
)

# These two explicit composition seams cross physical groups while retaining
# their established Julia ownership. Other files still obey the depth check.
const CROSS_GROUP_COMPOSITION_OWNERS = Set([
    joinpath("src", "orchestration", "simulation", "composition", "simulation_model.jl"),
    joinpath("src", "models", "gnc", "planning", "rpo", "comparison", "reporting_compat.jl"),
])

const ALLOWED_RAW_INCLUDE_FILES = Set([
    # Ordered model manifests and planning composition, not duplicate owners.
    joinpath("src", "models", "environment", "gravity", "gravity_models.jl"),
    joinpath("src", "models", "gnc", "companions", "hypr", "robot_arm.jl"),
    joinpath("src", "models", "gnc", "planning", "robot_arm", "planning.jl"),
    joinpath("src", "models", "gnc", "planning", "rpo", "comparison", "planner_comparison.jl"),
    joinpath("src", "models", "gnc", "planning", "rpo", "comparison", "reporting_compat.jl"),
    # The comparison report module includes its figure functions, kept in one file of their own.
    joinpath("src", "interfaces", "analysis", "planning", "rpo_comparison_reports.jl"),
    joinpath("src", "models", "gnc", "planning", "rpo", "lifecycle.jl"),
    # The RPO planner contract loads its validation and headroom parts; it sat in
    # gnc/interfaces/, outside these scan roots, until the GNC regroup.
    joinpath("src", "models", "gnc", "planning", "rpo", "contract", "rpo_planner.jl"),
    joinpath("src", "orchestration", "simulation", "composition", "simulation_model.jl"),
    joinpath("src", "models", "environment", "physical_models.jl"),
    joinpath("src", "models", "environment", "ephemerides", "ephemerides_models.jl"),
    joinpath("src", "models", "environment", "ephemerides", "planets.jl"),
    joinpath("src", "models", "gnc", "control", "hooks.jl"),
    joinpath("src", "models", "gnc", "control", "maneuvers", "propulsive_maneuvers.jl"),
    joinpath("src", "models", "gnc", "guidance", "models.jl"),
    joinpath("src", "models", "gnc", "guidance", "hooks.jl"),
    # EDG physical services include only their reviewed capability wrapper file.
    # The ownership gate checks that exact include; algorithms.jl has none.
    joinpath("src", "models", "gnc", "guidance", "aerobraking", "typed_edg", "services.jl"),
    # Existing helper module now loads its shared and RRT source owners.
    joinpath("src", "models", "gnc", "companions", "hypr", "hypr_utils.jl"),
    joinpath("src", "models", "gnc", "navigation", "hooks.jl"),
    joinpath("src", "orchestration", "simulation", "engine", "simulation_engine.jl"),
    joinpath("src", "orchestration", "simulation", "engine", "setup.jl"),
    joinpath("src", "orchestration", "simulation", "callbacks", "callbacks.jl"),
    joinpath("src", "orchestration", "simulation", "callbacks", "registry.jl"),
    joinpath("src", "orchestration", "simulation", "callbacks", "density_callbacks.jl"),
    joinpath("src", "orchestration", "simulation", "callbacks", "gram_track_cache.jl"),
    joinpath("src", "orchestration", "parallel", "policy", "parallel_policy.jl"),
    joinpath("src", "orchestration", "parallel", "routing", "parallel_profiles.jl"),
    joinpath("src", "orchestration", "parallel", "process", "parallel_process.jl"),
    # Fourth src/orchestration/parallel/<area>/<area>.jl module aggregator, same shape as the
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
    joinpath("src", "orchestration", "parallel", "cost", "parallel_cost.jl"),
])

violations = String[]
for root in SCAN_ROOTS
    isdir(root) || error("Missing include ownership scan root: $(relpath(root, REPO_ROOT))")
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
                rel in CROSS_GROUP_COMPOSITION_OWNERS && continue
                occursin(rx, active_src) || continue
                push!(violations, "$rel: contains forbidden cross-layer include chain pattern '$rx'.")
            end

            if startswith(rel, joinpath("src", "models", "gnc", "guidance", "aerobraking"))
                if EDGOwnershipChecks.has_control_include(src)
                    push!(violations, "$rel: guidance aerobraking file includes control source directly.")
                end
            end

            if rel == joinpath("src", "models", "gnc", "guidance", "aerobraking", "typed_edg", "algorithms.jl")
                append!(violations, EDGOwnershipChecks.aggregator_violations(src))
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
