const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const SCAN_ROOTS = (
    joinpath(REPO_ROOT, "examples"),
    joinpath(REPO_ROOT, "test"),
    joinpath(REPO_ROOT, "docs"),
    joinpath(REPO_ROOT, ".github", "workflows")
)

const FORBIDDEN_PATH_TOKENS = (
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
    "src/core/utils/typed_example_utils.jl"
)

@inline function _should_scan(path::String)::Bool
    endswith(path, ".jl") && return true
    endswith(path, ".md") && return true
    endswith(path, ".yml") && return true
    endswith(path, ".yaml") && return true
    return false
end

@inline function _is_excluded(path::String)::Bool
    rel = relpath(path, REPO_ROOT)
    # Historical AI notes are archival and excluded from canonical-path enforcement.
    startswith(rel, joinpath("test", "ai_")) && return true
    startswith(rel, joinpath("test", "ci_")) && endswith(rel, "_gate.jl") && return true
    rel == joinpath("docs", "architecture", "canonical_topology_contract.md") && return true
    rel == joinpath("docs", "architecture", "src_completeness_contract.md") && return true
    rel == joinpath("docs", "architecture", "src_canonical_owner_audit.md") && return true
    rel == joinpath("docs", "src", "generated", "contracts", "architecture", "canonical_topology_contract.md") && return true
    rel == joinpath("docs", "src", "generated", "contracts", "architecture", "src_completeness_contract.md") && return true
    rel == joinpath("docs", "src", "generated", "contracts", "architecture", "src_canonical_owner_audit.md") && return true
    return false
end

violations = String[]
for root in SCAN_ROOTS
    isdir(root) || continue
    for (dir, _, files) in walkdir(root)
        for file in files
            path = joinpath(dir, file)
            _should_scan(path) || continue
            _is_excluded(path) && continue
            rel = relpath(path, REPO_ROOT)
            src = read(path, String)
            for tok in FORBIDDEN_PATH_TOKENS
                occursin(tok, src) || continue
                push!(violations, "$(rel): contains forbidden canonical-path token '$tok'")
            end
        end
    end
end

if !isempty(violations)
    println("canonical_path_contract_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Canonical path contract gate failed")
end

println("canonical_path_contract_gate_ok")
