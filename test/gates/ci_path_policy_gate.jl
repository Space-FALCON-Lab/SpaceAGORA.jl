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
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

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
    joinpath("src", "dynamics", "models"),
    joinpath("src", "analysis", "plotting"),
    joinpath("src", "analysis", "reports"),
    joinpath("src", "gnc", "estimation"),
    joinpath("src", "mission", "constellation"),
    joinpath("src", "parallel", "state"),
    joinpath("src", "parallel", "tuning"),
    joinpath("src", "simulation", "events"),
    joinpath("src", "simulation", "execution"),
    joinpath("src", "simulation", "solver_orchestration"),
    joinpath("src", "vehicle", "actuators", "laser_terminal"),
    joinpath("src", "vehicle", "resources"),
    joinpath("src", "vehicle", "sensors"),
    joinpath("src", "gnc", "guidance", "targeting_control"),
    joinpath("src", "gnc", "control", "heatload_control"),
    # retired source files
    joinpath("src", "core", "utils", "typed_example_utils.jl"),
    joinpath("src", "gnc", "control", "effectors.jl"),
    joinpath("src", "gnc", "guidance", "effectors.jl"),
    joinpath("src", "gnc", "navigation", "effectors.jl"),
    joinpath("src", "gnc", "control", "bridge_helpers.jl"),
    joinpath("src", "mission", "campaigns", "run.jl"),
    joinpath("src", "mission", "campaigns", "set_and_run.jl"),
    joinpath("src", "mission", "campaigns", "define_mission.jl"),
    joinpath("src", "mission", "campaigns", "mission_model.jl"),
    joinpath("src", "mission", "campaigns", "initial_cond_calc.jl"),
    joinpath("src", "vehicle", "actuators", "thruster_effectors.jl"),
    joinpath("src", "vehicle", "spacecraft", "spacecraft_analysis.jl"),
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
)

# Canonical owner paths that must exist.
const REQUIRED_PATHS = (
    joinpath("src", "io", "config", "io_config.jl"),
    joinpath("src", "io", "serialization", "io_serialization.jl"),
    joinpath("src", "io", "logging", "io_logging.jl"),
    joinpath("src", "io", "outputs", "io_outputs.jl"),
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
        name = "retired dynamics/models tree",
        roots = ("src", "test"),
        exts = TEXT_EXTS,
        tokens = ("src/dynamics/models/", "\"dynamics\", \"models\"", "\"dynamics\",\"models\""),
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
            "src/simulation/events/",
            "src/simulation/execution/",
            "src/simulation/solver_orchestration/",
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
            "src/dynamics/models/",
            "src/integrator/",
            "src/utils/",
            "src/simulation_model/",
            "src/mission/campaigns/define_mission.jl",
            "src/mission/campaigns/mission_model.jl",
            "src/mission/campaigns/initial_cond_calc.jl",
            "src/mission/campaigns/run.jl",
            "src/mission/campaigns/set_and_run.jl",
            "src/vehicle/actuators/thruster_effectors.jl",
            "src/vehicle/spacecraft/spacecraft_analysis.jl",
            "src/gnc/guidance/targeting_control/",
            "src/gnc/control/heatload_control/",
            "src/simulation/events/",
            "src/simulation/execution/",
            "src/simulation/solver_orchestration/",
            "src/examples/",
            "src/analysis/plotting/",
            "src/core/utils/typed_example_utils.jl",
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

violations = String[]

for rel in RETIRED_PATHS
    ispath(joinpath(REPO_ROOT, rel)) && push!(violations, "$rel: retired path must not exist")
end

for rel in REQUIRED_PATHS
    isfile(joinpath(REPO_ROOT, rel)) || push!(violations, "$rel: canonical owner file is missing")
end

scan_roots = unique(vcat([collect(rule.roots) for rule in TOKEN_RULES]...))
for root_rel in scan_roots
    root = joinpath(REPO_ROOT, root_rel)
    isdir(root) || continue
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
