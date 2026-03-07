const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

function _read(rel)
    path = joinpath(REPO_ROOT, rel)
    isfile(path) || error("Missing required file: $rel")
    return read(path, String)
end

spaceagora_src = _read(joinpath("src", "SpaceAGORA.jl"))
engine_src = _read(joinpath("src", "simulation", "engine", "simulation_engine.jl"))
telemetry_src = _read(joinpath("src", "analysis", "verification", "telemetry_verification.jl"))
rhs_src = _read(joinpath("src", "simulation", "engine", "dynamics_rhs.jl"))
force_torque_src = _read(joinpath("src", "dynamics", "coupled", "force_torque_models.jl"))
guidance_hooks_src = _read(joinpath("src", "gnc", "guidance", "guidance_hooks.jl"))
control_hooks_src = _read(joinpath("src", "gnc", "control", "control_hooks.jl"))
callback_registry_src = _read(joinpath("src", "simulation", "callbacks", "registry.jl"))

occursin("run_simulation(args...; kwargs...) = SimulationEngine.run_simulation(args...; kwargs...)", spaceagora_src) ||
    error("SpaceAGORA.run_simulation is not forwarded to SimulationEngine.")
occursin("include(joinpath(@__DIR__, \"core\", \"simulation_model.jl\"))", spaceagora_src) ||
    error("SpaceAGORA does not include src/core/simulation_model.jl as the canonical sibling module.")
findfirst("include(joinpath(@__DIR__, \"core\", \"simulation_model.jl\"))", spaceagora_src) <
    findfirst("include(joinpath(@__DIR__, \"simulation\", \"engine\", \"simulation_engine.jl\"))", spaceagora_src) ||
    error("SpaceAGORA is not loading SimulationModel before SimulationEngine.")
occursin("using ..SimulationModel", engine_src) ||
    error("SimulationEngine is not loading SimulationModel as a sibling module.")
!occursin("include(joinpath(@__DIR__, \"..\", \"..\", \"core\", \"simulation_model.jl\"))", engine_src) ||
    error("SimulationEngine still directly includes core/simulation_model.jl.")

!occursin("simulation\", \"execution\", \"run_simulation.jl", telemetry_src) ||
    error("TelemetryVerification still directly includes simulation/execution/run_simulation.jl.")
!occursin("simulation\", \"engine\", \"simulation_engine.jl", telemetry_src) ||
    error("TelemetryVerification still directly includes simulation/engine/simulation_engine.jl.")
!occursin("SimulationEngine.SimulationModel", telemetry_src) ||
    error("TelemetryVerification still references SimulationEngine.SimulationModel.")

occursin("SimulationEngine.run_simulation", telemetry_src) ||
    error("TelemetryVerification is not calling SimulationEngine.run_simulation.")

for child_module in (
    "module GravityEffectors",
    "module AerodynamicEffectors",
    "module PerturbationEffectors",
    "module ThrusterModels",
    "module GuidanceModels",
)
    occursin(child_module, force_torque_src) ||
        error("DynamicEffectors split contract violation: missing child module declaration '$child_module'.")
end

occursin("using ..DynamicEffectors.GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel", guidance_hooks_src) ||
    error("GuidanceHooks is not importing the guidance model from DynamicEffectors.GuidanceModels.")
occursin("using ..DynamicEffectors.ThrusterModels: BaseThrusterModel", guidance_hooks_src) ||
    error("GuidanceHooks is not importing BaseThrusterModel from DynamicEffectors.ThrusterModels.")
occursin("using ..DynamicEffectors.GravityEffectors: aerobraking_gravity_force_ii", guidance_hooks_src) ||
    error("GuidanceHooks is not importing aerobraking_gravity_force_ii from DynamicEffectors.GravityEffectors.")
occursin("using ..DynamicEffectors.GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel", control_hooks_src) ||
    error("ControlHooks is not importing the guidance model from DynamicEffectors.GuidanceModels.")
occursin("using ..DynamicEffectors.ThrusterModels: BaseThrusterModel", control_hooks_src) ||
    error("ControlHooks is not importing BaseThrusterModel from DynamicEffectors.ThrusterModels.")
occursin("using ..DynamicEffectors.GravityEffectors: aerobraking_gravity_force_ii", control_hooks_src) ||
    error("ControlHooks is not importing aerobraking_gravity_force_ii from DynamicEffectors.GravityEffectors.")
occursin("using ..DynamicEffectors.ThrusterModels: BaseThrusterModel", callback_registry_src) ||
    error("Callback registry is not importing BaseThrusterModel from DynamicEffectors.ThrusterModels.")
occursin("using ..DynamicEffectors.AerodynamicEffectors: AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight", callback_registry_src) ||
    error("Callback registry is not importing aerodynamic coefficient models from DynamicEffectors.AerodynamicEffectors.")
occursin("using ..DynamicEffectors.GravityEffectors: InverseSquaredJ2GravityModel", callback_registry_src) ||
    error("Callback registry is not importing gravity models from DynamicEffectors.GravityEffectors.")

isfile(joinpath(REPO_ROOT, "src", "simulation", "execution", "run_simulation.jl")) &&
    error("Wave 2 contract violation: legacy simulation execution wrapper still exists at src/simulation/execution/run_simulation.jl")
isfile(joinpath(REPO_ROOT, "src", "simulation", "execution", "simulation_elements.jl")) &&
    error("Canonical-ownership contract violation: legacy simulation elements file still exists at src/simulation/execution/simulation_elements.jl")
isfile(joinpath(REPO_ROOT, "src", "simulation", "execution", "simulation_execution.jl")) &&
    error("Canonical-ownership contract violation: legacy simulation execution file still exists at src/simulation/execution/simulation_execution.jl")
isfile(joinpath(REPO_ROOT, "src", "simulation", "execution", "run.jl")) ||
    error("Canonical-ownership contract violation: missing simulation execution dispatcher at src/simulation/execution/run.jl")
isfile(joinpath(REPO_ROOT, "src", "simulation", "execution", "dispatch.jl")) ||
    error("Canonical-ownership contract violation: missing simulation execution dispatcher at src/simulation/execution/dispatch.jl")

run_dispatch_src = _read(joinpath("src", "simulation", "execution", "run.jl"))
occursin("dispatch.jl", run_dispatch_src) ||
    error("Canonical-ownership contract violation: simulation/execution/run.jl is not dispatching through dispatch.jl.")

for retired in (
    joinpath("src", "mission", "campaigns", "run.jl"),
    joinpath("src", "mission", "campaigns", "set_and_run.jl"),
    joinpath("src", "mission", "campaigns", "define_mission.jl"),
    joinpath("src", "mission", "campaigns", "mission_model.jl"),
    joinpath("src", "mission", "campaigns", "initial_cond_calc.jl"),
    joinpath("src", "vehicle", "actuators", "thruster_effectors.jl"),
    joinpath("src", "vehicle", "spacecraft", "spacecraft_analysis.jl")
)
    isfile(joinpath(REPO_ROOT, retired)) &&
        error("Canonical-ownership contract violation: retired file still exists at $retired")
end

isdir(joinpath(REPO_ROOT, "src", "benchmarks")) &&
    error("Canonical-ownership contract violation: src/benchmarks must not exist (use top-level benchmarks/).")
isdir(joinpath(REPO_ROOT, "src", "dynamics", "models")) &&
    error("Canonical-ownership contract violation: src/dynamics/models must not exist.")
for retired_dir in (
    joinpath("src", "examples"),
    joinpath("src", "analysis", "plotting"),
    joinpath("src", "analysis", "reports"),
    joinpath("src", "parallel", "state"),
    joinpath("src", "parallel", "tuning"),
    joinpath("src", "vehicle", "sensors"),
)
    isdir(joinpath(REPO_ROOT, retired_dir)) &&
        error("Canonical-ownership contract violation: retired source directory still exists at $retired_dir")
end
isfile(joinpath(REPO_ROOT, "src", "core", "utils", "typed_example_utils.jl")) &&
    error("Canonical-ownership contract violation: retired example helper still exists at src/core/utils/typed_example_utils.jl")

for rel in (
    joinpath("src", "mission", "initial_conditions.jl"),
    joinpath("src", "vehicle", "actuators", "actuator_hooks.jl"),
    joinpath("src", "vehicle", "actuators", "thruster", "thruster_hooks.jl"),
    joinpath("src", "vehicle", "structure", "structure_models.jl"),
    joinpath("src", "vehicle", "structure", "assembly_graph.jl"),
    joinpath("src", "vehicle", "structure", "mass_properties.jl"),
    joinpath("src", "vehicle", "structure", "geometry_properties.jl"),
    joinpath("src", "dynamics", "coupled", "force_torque_models.jl"),
    joinpath("src", "gnc", "guidance", "aerobraking", "interfaces.jl"),
    joinpath("src", "gnc", "guidance", "aerobraking", "dispatcher.jl"),
    joinpath("src", "gnc", "guidance", "aerobraking", "e_edg", "e_edg_strategy.jl"),
    joinpath("src", "gnc", "guidance", "aerobraking", "t_edg", "t_edg_strategy.jl"),
    joinpath("src", "mission", "operations", "aerobraking_policy", "policy_types.jl"),
)
    isfile(joinpath(REPO_ROOT, rel)) ||
        error("Canonical-ownership contract violation: missing canonical owner file: $rel")
end

isdir(joinpath(REPO_ROOT, "src", "gnc", "guidance", "targeting_control")) &&
    error("Canonical-ownership contract violation: retired targeting_control path still exists.")

for retired in ("control", "physical_models", "guidance", "integrator", "utils", "simulation_model")
    isdir(joinpath(REPO_ROOT, "src", retired)) &&
        error("Canonical-ownership contract violation: retired source tree still exists at src/$retired")
end

# Gravity ownership contract:
# - canonical gravity owner is environment/gravity
# - translational folder must not own gravity models
isfile(joinpath(REPO_ROOT, "src", "environment", "gravity", "gravity_models.jl")) ||
    error("Canonical-ownership contract violation: missing gravity owner at src/environment/gravity/gravity_models.jl")
isfile(joinpath(REPO_ROOT, "src", "dynamics", "translational", "gravity_models.jl")) &&
    error("Canonical-ownership contract violation: gravity model must not be owned by src/dynamics/translational")

# Translational ownership contract:
# - translational RHS equations are owned by dynamics/translational and consumed by engine RHS.
for rel in (
    joinpath("src", "dynamics", "translational", "translational_models.jl"),
    joinpath("src", "dynamics", "translational", "position_kinematics.jl"),
    joinpath("src", "dynamics", "translational", "point_mass_dynamics.jl"),
)
    isfile(joinpath(REPO_ROOT, rel)) ||
        error("Canonical-ownership contract violation: missing translational owner file: $rel")
end

occursin("DynamicsTranslational.assign_full_translational_rhs!", rhs_src) ||
    error("Engine RHS is not using canonical translational owner for full dynamics path.")
occursin("DynamicsTranslational.assign_slow_translational_rhs!", rhs_src) ||
    error("Engine RHS is not using canonical translational owner for slow dynamics path.")
occursin("DynamicsTranslational.assign_control_only_translational_rhs!", rhs_src) ||
    error("Engine RHS is not using canonical translational owner for control-only path.")

engine_dir = joinpath(REPO_ROOT, "src", "simulation", "engine")
for (root, _, files) in walkdir(engine_dir)
    for file in files
        endswith(file, ".jl") || continue
        rel = relpath(joinpath(root, file), REPO_ROOT)
        rel in (
            joinpath("src", "simulation", "engine", "adapters", "from_env.jl"),
            joinpath("src", "simulation", "engine", "adapters", "from_simulation_configuration.jl")
        ) && continue
        src = read(joinpath(root, file), String)
        if occursin("get(ENV", src) || occursin("ENV[", src)
            error("ENV usage found outside engine adapters: $rel")
        end
    end
end

perf_launcher = _read(joinpath("benchmarks", "studies", "performance_runtime_analysis.jl"))
(
    occursin(joinpath("performance_runtime_analysis", "main.jl"), perf_launcher) ||
    occursin("\"performance_runtime_analysis\", \"main.jl\"", perf_launcher)
) ||
    error("Canonical runtime-analysis launcher is not forwarding to split main.jl.")
for rel in (
    joinpath("benchmarks", "studies", "performance_runtime_analysis", "main.jl"),
    joinpath("benchmarks", "studies", "performance_runtime_analysis", "case_catalog.jl"),
    joinpath("benchmarks", "studies", "performance_runtime_analysis", "measurement.jl"),
    joinpath("benchmarks", "studies", "performance_runtime_analysis", "reporting.jl"),
    joinpath("benchmarks", "studies", "performance_runtime_analysis", "cli.jl")
)
    isfile(joinpath(REPO_ROOT, rel)) || error("Missing split runtime-analysis file: $rel")
end

wrapper_files = (
    "test/performance_runtime_analysis.jl",
    "test/performance_smart_parallel_ladder.jl",
    "test/performance_smart_parallel_ladder_cross_machine.jl",
    "test/performance_split_imex_compare.jl",
    "test/performance_static_vs_parallel.jl",
    "test/performance_effector_reduction_microbench.jl",
    "test/performance_paper_pipeline.jl",
    "test/gram_interpolation_vs_point_to_point_analysis.jl",
    "test/gram_single_call_vs_point_to_point_analysis.jl",
    "test/gram_real_sim_runtime_compare.jl",
    "test/gram_real_sim_surrogate_matrix.jl",
    "test/gram_real_sim_surrogate_decision_table.jl",
    "test/gram_offline_db_accuracy_study.jl",
    "test/gram_planet_rho_altitude_sweep.jl",
    "test/telemetry_hybrid_tuner.jl",
    "test/telemetry_odyssey_tuner.jl",
    "test/telemetry_orbit_accuracy_study.jl",
    "test/telemetry_orbit_accuracy_plots.jl"
)
for rel in wrapper_files
    isfile(joinpath(REPO_ROOT, rel)) &&
        error("Wave 2 contract violation: benchmark wrapper still exists: $rel")
end

println("architecture_contract_gate_ok")
