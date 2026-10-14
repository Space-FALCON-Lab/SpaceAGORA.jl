# Single path-policy gate for the repository topology contracts:
#   1. retired paths that must stay absent,
#   2. canonical owner paths that must exist,
#   3. forbidden path tokens, scanned per root and file type.
#
# It replaces the former ci_no_dynamics_models_gate, ci_no_src_benchmarks_root_gate,
# ci_no_legacy_ownership_gate, ci_final_clean_contract_gate,
# ci_io_surface_nonempty_gate and ci_canonical_path_contract_gate, the path half
# of ci_benchmark_wrapper_parity_gate, and the retired-path loops that
# ci_architecture_contract_gate and ci_src_completeness_contract_gate used to
# repeat. Record new retirements here instead of adding a gate.
using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Physical grouping only: these folders do not add Julia namespace modules.
const SOURCE_GROUPS = (
    "interfaces" => ("cli", "io", "analysis"),
    "orchestration" => ("simulation", "mission", "parallel"),
    "models" => ("vehicle", "environment", "dynamics", "gnc"),
    "foundation" => ("core", "assets"),
)

function source_layout_violations(root)
    errors = String[]
    source = joinpath(root, "src")
    isdir(source) || return ["src: source root is missing"]
    actual_groups = Set(name for name in readdir(source) if isdir(joinpath(source, name)))
    expected_groups = Set(first.(SOURCE_GROUPS))
    actual_groups == expected_groups || push!(errors,
        "src: expected exactly the four source groups; found $(sort!(collect(actual_groups)))")
    for (group, domains) in SOURCE_GROUPS
        parent = joinpath(source, group)
        if !isdir(parent) || islink(parent)
            push!(errors, "src/$group: source group must be a real directory")
            continue
        end
        actual_domains = Set(name for name in readdir(parent) if isdir(joinpath(parent, name)))
        actual_domains == Set(domains) || push!(errors,
            "src/$group: domain directories do not match the declared layout")
        for domain in domains
            relative = joinpath("src", group, domain)
            path = joinpath(root, relative)
            if !isdir(path) || islink(path)
                push!(errors, "$relative: required source domain must be a real directory")
            elseif !any(endswith(file, ".jl") for (_, _, files) in walkdir(path) for file in files)
                push!(errors, "$relative: source domain has no Julia implementation files")
            end
            old_path = joinpath(source, domain)
            (ispath(old_path) || islink(old_path)) &&
                push!(errors, "src/$domain: former domain root must remain absent")
        end
    end
    for name in ("SpaceAGORA.jl", "precompile_workload.jl")
        isfile(joinpath(source, name)) || push!(errors, "src/$name: package entry/helper must remain at the source root")
    end
    extensions = joinpath(root, "ext")
    isdir(extensions) && !islink(extensions) || push!(errors, "ext: Julia extensions must remain at the repository root")
    project_path = joinpath(root, "Project.toml")
    if isfile(project_path)
        for name in keys(get(TOML.parsefile(project_path), "extensions", Dict()))
            (isfile(joinpath(extensions, name * ".jl")) ||
             isfile(joinpath(extensions, name, name * ".jl"))) ||
                push!(errors, "ext/$name: registered Julia extension entry is missing")
        end
    else
        push!(errors, "Project.toml: package metadata is missing")
    end
    return errors
end

# Retired by the topology cleanup; must not reappear (relative to REPO_ROOT).
const RETIRED_PATHS = (
    # retired source trees
    joinpath("src", "control"),
    joinpath("src", "physical_models"),
    joinpath("src", "guidance"),
    joinpath("src", "integrator"),
    joinpath("src", "utils"),
    joinpath("src", "simulation_model"),
    joinpath("src", "benchmarks"),
    joinpath("src", "examples"),
    joinpath("src", "models", "dynamics", "models"),
    joinpath("src", "models", "dynamics", "multibody_cloth"),
    joinpath("src", "interfaces", "analysis", "plotting"),
    joinpath("src", "interfaces", "analysis", "reports"),
    joinpath("src", "models", "gnc", "estimation"),
    joinpath("src", "orchestration", "mission", "constellation"),
    joinpath("src", "orchestration", "parallel", "state"),
    joinpath("src", "orchestration", "parallel", "tuning"),
    joinpath("src", "orchestration", "simulation", "events"),
    joinpath("src", "orchestration", "simulation", "execution"),
    joinpath("src", "orchestration", "simulation", "solver_orchestration"),
    joinpath("src", "models", "vehicle", "actuators", "laser_terminal"),
    joinpath("src", "models", "vehicle", "resources"),
    joinpath("src", "models", "vehicle", "sensors"),
    joinpath("src", "models", "gnc", "guidance", "targeting_control"),
    joinpath("src", "models", "gnc", "control", "heatload_control"),
    # browser rendering belongs to the optional SpaceAGORAVisualizer companion
    joinpath("src", "viewer"),
    # retired source files
    joinpath("src", "foundation", "core", "state", "simulation_configuration.jl"),
    joinpath("src", "foundation", "core", "utils", "typed_example_utils.jl"),
    joinpath("src", "models", "gnc", "control", "effectors.jl"),
    joinpath("src", "models", "gnc", "guidance", "effectors.jl"),
    joinpath("src", "models", "gnc", "navigation", "effectors.jl"),
    joinpath("src", "models", "gnc", "control", "bridge_helpers.jl"),
    joinpath("src", "orchestration", "mission", "campaigns", "run.jl"),
    joinpath("src", "orchestration", "mission", "campaigns", "set_and_run.jl"),
    joinpath("src", "orchestration", "mission", "campaigns", "define_mission.jl"),
    joinpath("src", "orchestration", "mission", "campaigns", "mission_model.jl"),
    joinpath("src", "orchestration", "mission", "campaigns", "initial_cond_calc.jl"),
    joinpath("src", "models", "vehicle", "actuators", "thruster_effectors.jl"),
    joinpath("src", "models", "vehicle", "spacecraft", "spacecraft_analysis.jl"),
    # benchmark/study wrappers retired from test/ (canonical homes under benchmarks/)
    joinpath("test", "performance_runtime_analysis.jl"),
    joinpath("test", "performance_smart_parallel_ladder.jl"),
    joinpath("test", "performance_smart_parallel_ladder_cross_machine.jl"),
    joinpath("test", "performance_split_imex_compare.jl"),
    joinpath("test", "performance_static_vs_parallel.jl"),
    joinpath("test", "performance_effector_reduction_microbench.jl"),
    joinpath("test", "performance_paper_pipeline.jl"),
    joinpath("test", "gram_interpolation_vs_point_to_point_analysis.jl"),
    joinpath("test", "gram_single_call_vs_point_to_point_analysis.jl"),
    joinpath("test", "gram_real_sim_runtime_compare.jl"),
    joinpath("test", "gram_real_sim_surrogate_matrix.jl"),
    joinpath("test", "gram_real_sim_surrogate_decision_table.jl"),
    joinpath("test", "gram_offline_db_accuracy_study.jl"),
    joinpath("test", "gram_planet_rho_altitude_sweep.jl"),
    joinpath("test", "telemetry_hybrid_tuner.jl"),
    joinpath("test", "telemetry_odyssey_tuner.jl"),
    joinpath("test", "telemetry_orbit_accuracy_study.jl"),
    joinpath("test", "telemetry_orbit_accuracy_plots.jl"),
    # completed one-time migrations
    joinpath("scripts", "dev", "test_reorg_b1.sh"),
)

# Canonical owner paths that must exist.
const REQUIRED_PATHS = (
    joinpath("src", "models", "dynamics", "compliant_multibody", "compliant_multibody.jl"),
    joinpath("src", "models", "dynamics", "compliant_multibody", "compliant_robot_arm_dynamics.jl"),
    joinpath("src", "models", "dynamics", "compliant_multibody", "compliant_attachment.jl"),
    joinpath("src", "models", "dynamics", "compliant_multibody", "compliant_attachment_dynamics.jl"),
    joinpath("src", "orchestration", "simulation", "config", "configuration.jl"),
    joinpath("src", "orchestration", "simulation", "config", "run_settings.jl"),
    joinpath("src", "orchestration", "simulation", "config", "solver_settings.jl"),
    joinpath("src", "orchestration", "simulation", "config", "environment_settings.jl"),
    joinpath("src", "orchestration", "simulation", "config", "constellation_configuration.jl"),
    joinpath("src", "orchestration", "simulation", "config", "simulation_configuration.jl"),
    joinpath("src", "interfaces", "io", "paths", "io_paths.jl"),
    joinpath("src", "interfaces", "io", "checkpoints", "checkpoint_storage.jl"),
    joinpath("src", "interfaces", "io", "outputs", "io_outputs.jl"),
    joinpath("src", "interfaces", "analysis", "visualization", "scene", "scene_visualization.jl"),
    joinpath("benchmarks", "studies", "performance_runtime_analysis.jl"),
    joinpath("benchmarks", "studies", "performance_smart_parallel_ladder.jl"),
    joinpath("benchmarks", "studies", "parallelization_performance.jl"),
    joinpath("benchmarks", "studies", "performance_smart_parallel_ladder_cross_machine.jl"),
    joinpath("benchmarks", "studies", "performance_split_imex_compare.jl"),
    joinpath("benchmarks", "studies", "performance_static_vs_parallel.jl"),
    joinpath("benchmarks", "studies", "performance_effector_reduction_microbench.jl"),
    joinpath("benchmarks", "scripts", "performance_paper_pipeline.jl"),
    joinpath("benchmarks", "studies", "gram_interpolation_vs_point_to_point_analysis.jl"),
    joinpath("benchmarks", "studies", "gram_single_call_vs_point_to_point_analysis.jl"),
    joinpath("benchmarks", "studies", "gram_real_sim_runtime_compare.jl"),
    joinpath("benchmarks", "studies", "gram_real_sim_surrogate_matrix.jl"),
    joinpath("benchmarks", "studies", "gram_real_sim_surrogate_decision_table.jl"),
    joinpath("benchmarks", "studies", "gram_offline_db_accuracy_study.jl"),
    joinpath("benchmarks", "studies", "gram_planet_rho_altitude_sweep.jl"),
    joinpath("benchmarks", "studies", "telemetry_hybrid_tuner.jl"),
    joinpath("benchmarks", "studies", "telemetry_odyssey_tuner.jl"),
    joinpath("benchmarks", "studies", "telemetry_orbit_accuracy_study.jl"),
    joinpath("benchmarks", "studies", "telemetry_orbit_accuracy_plots.jl"),
)

const TEXT_EXTS = (".jl", ".md", ".yml", ".yaml")

# Forbidden tokens: (rule name, scan roots, file extensions, tokens).
const TOKEN_RULES = (
    (
        name = "retired compliant multibody source tree",
        roots = ("src", "test", "docs"),
        exts = TEXT_EXTS,
        tokens = ("multibody_cloth",),
    ),
    (
        name = "retired dynamics/models tree",
        roots = ("src", "test"),
        exts = TEXT_EXTS,
        tokens = ("src/models/dynamics/models/", "\"dynamics\", \"models\"", "\"dynamics\",\"models\""),
    ),
    (
        name = "src/benchmarks root",
        roots = ("src", "test"),
        exts = TEXT_EXTS,
        tokens = ("src/benchmarks/",),
    ),
    (
        name = "legacy ownership path",
        roots = ("src",),
        exts = (".jl",),
        tokens = (
            "src/control/",
            "src/physical_models/",
            "src/guidance/",
            "src/integrator/",
            "src/orchestration/simulation/events/",
            "src/orchestration/simulation/execution/",
            "src/orchestration/simulation/solver_orchestration/",
            "src/utils/",
            "src/simulation_model/",
            "Compatibility wrapper: canonical path forwarding to legacy implementation.",
        ),
    ),
    (
        name = "retired identifier",
        roots = ("src",),
        exts = (".jl",),
        tokens = ("legacy_", "ControlEffectors", "GuidanceEffectors", "NavigationEffectors"),
    ),
    (
        name = "canonical path",
        roots = ("examples", "test", "docs", joinpath(".github", "workflows")),
        exts = TEXT_EXTS,
        tokens = (
            "SimulationEngine.SimulationModel",
            "SimulationModel.SPICE_LOCK",
            "SimulationModel.GRAM_LOCK",
            "src/benchmarks/",
            "src/control/",
            "src/physical_models/",
            "src/guidance/",
            "src/models/dynamics/models/",
            "src/integrator/",
            "src/utils/",
            "src/simulation_model/",
            "src/orchestration/mission/campaigns/define_mission.jl",
            "src/orchestration/mission/campaigns/mission_model.jl",
            "src/orchestration/mission/campaigns/initial_cond_calc.jl",
            "src/orchestration/mission/campaigns/run.jl",
            "src/orchestration/mission/campaigns/set_and_run.jl",
            "src/models/vehicle/actuators/thruster_effectors.jl",
            "src/models/vehicle/spacecraft/spacecraft_analysis.jl",
            "src/models/gnc/guidance/targeting_control/",
            "src/models/gnc/control/heatload_control/",
            "src/orchestration/simulation/events/",
            "src/orchestration/simulation/execution/",
            "src/orchestration/simulation/solver_orchestration/",
            "src/examples/",
            "src/interfaces/analysis/plotting/",
            "src/foundation/core/utils/typed_example_utils.jl",
        ),
    ),
)

# Files exempt from token scans: the gates themselves (they carry the token
# lists), archival AI notes, and the contract documents that describe the
# retired paths on purpose (including their docs-site mirrors).
@inline function _scan_exempt(rel::String)::Bool
    startswith(rel, joinpath("test", "gates", "ci_")) && endswith(rel, "_gate.jl") && return true
    startswith(rel, joinpath("test", "ai_")) && return true
    for doc in (
        joinpath("docs", "architecture", "canonical_topology_contract.md"),
        joinpath("docs", "architecture", "src_completeness_contract.md"),
        joinpath("docs", "architecture", "src_canonical_owner_audit.md"),
    )
        rel == doc && return true
        rel == joinpath("docs", "src", "generated", "contracts", relpath(doc, "docs")) && return true
    end
    return false
end

violations = source_layout_violations(REPO_ROOT)

for rel in RETIRED_PATHS
    ispath(joinpath(REPO_ROOT, rel)) && push!(violations, "$rel: retired path must not exist")
end

for rel in REQUIRED_PATHS
    isfile(joinpath(REPO_ROOT, rel)) || push!(violations, "$rel: canonical owner file is missing")
end

scan_roots = unique(vcat([collect(rule.roots) for rule in TOKEN_RULES]...))
for root_rel in scan_roots
    root = joinpath(REPO_ROOT, root_rel)
    isdir(root) || error("Missing path-policy scan root: $root_rel")
    rules = [rule for rule in TOKEN_RULES if root_rel in rule.roots]
    for (dir, _, files) in walkdir(root)
        for file in files
            path = joinpath(dir, file)
            rel = relpath(path, REPO_ROOT)
            _scan_exempt(rel) && continue
            ext = splitext(file)[2]
            applicable = [rule for rule in rules if ext in rule.exts]
            isempty(applicable) && continue
            src = read(path, String)
            for rule in applicable, tok in rule.tokens
                occursin(tok, src) || continue
                push!(violations, "$rel: contains forbidden $(rule.name) token '$tok'")
            end
        end
    end
end

if !isempty(violations)
    println("path_policy_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Path policy gate failed")
end

println("path_policy_gate_ok")
