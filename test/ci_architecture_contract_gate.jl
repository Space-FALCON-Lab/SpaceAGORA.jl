const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

function _read(rel)
    path = joinpath(REPO_ROOT, rel)
    isfile(path) || error("Missing required file: $rel")
    return read(path, String)
end

spaceagora_src = _read(joinpath("src", "SpaceAGORA.jl"))
runtime_services_src = _read(joinpath("src", "simulation", "runtime_services.jl"))
simulation_model_src = _read(joinpath("src", "core", "simulation_model.jl"))
reference_system_config_src = _read(joinpath("src", "core", "state", "reference_system_config.jl"))
spacecraft_model_src = _read(joinpath("src", "vehicle", "spacecraft", "model.jl"))
legacy_model_codes_src = _read(joinpath("src", "core", "types", "legacy_model_codes.jl"))
runtime_types_src = _read(joinpath("src", "core", "types", "runtime_types.jl"))
command_types_src = _read(joinpath("src", "gnc", "command_types.jl"))
top_level_thruster_models_src = _read(joinpath("src", "vehicle", "actuators", "thruster", "thruster_models_module.jl"))
top_level_guidance_models_src = _read(joinpath("src", "gnc", "guidance", "guidance_models.jl"))
top_level_gravity_effectors_src = _read(joinpath("src", "environment", "gravity", "gravity_effectors.jl"))
engine_src = _read(joinpath("src", "simulation", "engine", "simulation_engine.jl"))
engine_public_api_src = _read(joinpath("src", "simulation", "engine", "public_api.jl"))
engine_execution_src = _read(joinpath("src", "simulation", "engine", "execution.jl"))
telemetry_src = _read(joinpath("src", "analysis", "verification", "telemetry_verification.jl"))
rhs_src = _read(joinpath("src", "simulation", "engine", "dynamics_rhs.jl"))
force_torque_src = _read(joinpath("src", "dynamics", "coupled", "force_torque_models.jl"))
gravity_effectors_src = _read(joinpath("src", "dynamics", "coupled", "force_torque_models", "gravity_effectors.jl"))
aerodynamic_effectors_src = _read(joinpath("src", "dynamics", "coupled", "force_torque_models", "aerodynamic_effectors.jl"))
perturbation_effectors_src = _read(joinpath("src", "dynamics", "coupled", "force_torque_models", "perturbation_effectors.jl"))
thruster_models_src = _read(joinpath("src", "dynamics", "coupled", "force_torque_models", "thruster_models.jl"))
guidance_models_src = _read(joinpath("src", "dynamics", "coupled", "force_torque_models", "guidance_models.jl"))
guidance_hooks_src = _read(joinpath("src", "gnc", "guidance", "guidance_hooks.jl"))
control_hooks_src = _read(joinpath("src", "gnc", "control", "control_hooks.jl"))
callback_registry_src = _read(joinpath("src", "simulation", "callbacks", "registry.jl"))
density_callback_assembly_src = _read(joinpath("src", "simulation", "callbacks", "density_callbacks", "assembly.jl"))
density_models_src = _read(joinpath("src", "environment", "atmosphere", "density_models.jl"))
ephemerides_models_src = _read(joinpath("src", "environment", "ephemerides", "ephemerides_models.jl"))
setup_src = _read(joinpath("src", "simulation", "engine", "setup.jl"))
control_commands_src = _read(joinpath("src", "gnc", "control", "aerobraking", "control_commands.jl"))
constraint_tracking_src = _read(joinpath("src", "gnc", "control", "aerobraking", "constraint_tracking.jl"))
trajectory_predictor_src = _read(joinpath("src", "gnc", "guidance", "aerobraking", "t_edg", "trajectory_predictor.jl"))
eom_predictor_src = _read(joinpath("src", "gnc", "guidance", "aerobraking", "t_edg", "eom_predictor.jl"))

occursin("run_simulation(args...; kwargs...) = SimulationEngine.run_simulation(args...; kwargs...)", spaceagora_src) ||
    error("SpaceAGORA.run_simulation is not forwarded to SimulationEngine.")
occursin("include(joinpath(@__DIR__, \"simulation\", \"runtime_services.jl\"))", spaceagora_src) ||
    error("SpaceAGORA does not include src/simulation/runtime_services.jl as the canonical runtime-services sibling module.")
occursin("include(joinpath(@__DIR__, \"core\", \"simulation_model.jl\"))", spaceagora_src) ||
    error("SpaceAGORA does not include src/core/simulation_model.jl as the canonical sibling module.")
findfirst("include(joinpath(@__DIR__, \"simulation\", \"runtime_services.jl\"))", spaceagora_src) <
    findfirst("include(joinpath(@__DIR__, \"core\", \"simulation_model.jl\"))", spaceagora_src) ||
    error("SpaceAGORA is not loading RuntimeServices before SimulationModel.")
findfirst("include(joinpath(@__DIR__, \"core\", \"simulation_model.jl\"))", spaceagora_src) <
    findfirst("include(joinpath(@__DIR__, \"simulation\", \"engine\", \"simulation_engine.jl\"))", spaceagora_src) ||
    error("SpaceAGORA is not loading SimulationModel before SimulationEngine.")
occursin("module RuntimeServices", runtime_services_src) ||
    error("RuntimeServices owner module is missing.")
occursin("const SPICE_LOCK = ReentrantLock()", runtime_services_src) ||
    error("RuntimeServices is not the canonical owner of SPICE_LOCK.")
occursin("const GRAM_LOCK = ReentrantLock()", runtime_services_src) ||
    error("RuntimeServices is not the canonical owner of GRAM_LOCK.")
!occursin("const SPICE_LOCK", simulation_model_src) ||
    error("SimulationModel still owns SPICE_LOCK.")
!occursin("const GRAM_LOCK", simulation_model_src) ||
    error("SimulationModel still owns GRAM_LOCK.")
occursin("runtime_services.jl", simulation_model_src) ||
    error("SimulationModel standalone include contract is not bootstrapping RuntimeServices.")
occursin("include(joinpath(@__DIR__, \"..\", \"vehicle\", \"spacecraft\", \"model.jl\"))", simulation_model_src) ||
    error("SimulationModel is not including the canonical spacecraft model owner.")
occursin("@reexport using .SpacecraftModels", simulation_model_src) ||
    error("SimulationModel is not reexporting the canonical SpacecraftModels owner.")
!occursin("@reexport using .PhysicalModel", simulation_model_src) ||
    error("SimulationModel still reexports the retired PhysicalModel owner.")
occursin("include(joinpath(@__DIR__, \"..\", \"core\", \"types\", \"legacy_model_codes.jl\"))", simulation_model_src) ||
    error("SimulationModel does not include src/core/types/legacy_model_codes.jl before runtime types.")
findfirst("include(joinpath(@__DIR__, \"..\", \"core\", \"types\", \"legacy_model_codes.jl\"))", simulation_model_src) <
    findfirst("include(joinpath(@__DIR__, \"..\", \"core\", \"types\", \"runtime_types.jl\"))", simulation_model_src) ||
    error("SimulationModel is not loading LegacyModelCodes before ConfigTypes.")
occursin("module ReferenceSystems", reference_system_config_src) ||
    error("Reference system config is not using the CamelCase ReferenceSystems module.")
!occursin("module ref_sys", reference_system_config_src) ||
    error("Reference system config still uses the retired ref_sys module name.")
occursin("module SpacecraftModels", spacecraft_model_src) ||
    error("Spacecraft model owner is not using the canonical SpacecraftModels module name.")
!occursin("module PhysicalModel", spacecraft_model_src) ||
    error("Spacecraft model owner still uses the retired PhysicalModel module name.")
occursin("module LegacyModelCodes", legacy_model_codes_src) ||
    error("Legacy model selector enums were not split into src/core/types/legacy_model_codes.jl.")
occursin("using ..LegacyModelCodes:", runtime_types_src) ||
    error("ConfigTypes is not importing legacy model codes from LegacyModelCodes.")
occursin("using ..SpacecraftModels: SpacecraftModel", runtime_types_src) ||
    error("ConfigTypes is not importing SpacecraftModel from SpacecraftModels.")
!occursin("@enum Legacy", runtime_types_src) ||
    error("ConfigTypes still defines Legacy*Code enums directly in runtime_types.jl.")
occursin("include(joinpath(@__DIR__, \"..\", \"gnc\", \"command_types.jl\"))", simulation_model_src) ||
    error("SimulationModel does not include src/gnc/command_types.jl before runtime types.")
findfirst("include(joinpath(@__DIR__, \"..\", \"gnc\", \"command_types.jl\"))", simulation_model_src) <
    findfirst("include(joinpath(@__DIR__, \"..\", \"core\", \"types\", \"runtime_types.jl\"))", simulation_model_src) ||
    error("SimulationModel is not loading CommandTypes before ConfigTypes.")
occursin("module CommandTypes", command_types_src) ||
    error("CommandTypes owner module is missing.")
occursin("struct PropulsiveManeuverCommand", command_types_src) ||
    error("CommandTypes is missing PropulsiveManeuverCommand.")
occursin("struct AerobrakingControlCommand", command_types_src) ||
    error("CommandTypes is missing AerobrakingControlCommand.")
occursin("using ..CommandTypes: PropulsiveManeuverCommand", runtime_types_src) ||
    error("ConfigTypes is not importing PropulsiveManeuverCommand from CommandTypes.")
occursin("maneuver_commands::Vector{PropulsiveManeuverCommand}", runtime_types_src) ||
    error("SharedBuffers is missing the typed maneuver_commands buffer.")
occursin("include(joinpath(@__DIR__, \"..\", \"vehicle\", \"actuators\", \"thruster\", \"thruster_models_module.jl\"))", simulation_model_src) ||
    error("SimulationModel does not include the top-level ThrusterModels owner.")
occursin("module ThrusterModels", top_level_thruster_models_src) ||
    error("Top-level ThrusterModels owner module is missing.")
occursin("include(joinpath(@__DIR__, \"thruster_models.jl\"))", top_level_thruster_models_src) ||
    error("Top-level ThrusterModels owner is not including vehicle/actuators/thruster/thruster_models.jl.")
occursin("include(joinpath(@__DIR__, \"..\", \"gnc\", \"guidance\", \"guidance_models.jl\"))", simulation_model_src) ||
    error("SimulationModel does not include the top-level GuidanceModels owner.")
occursin("module GuidanceModels", top_level_guidance_models_src) ||
    error("Top-level GuidanceModels owner module is missing.")
occursin("include(joinpath(@__DIR__, \"thruster_guidance\", \"thruster_guidance_models.jl\"))", top_level_guidance_models_src) ||
    error("Top-level GuidanceModels owner is not including the thruster guidance model definitions.")
findfirst("include(joinpath(@__DIR__, \"..\", \"vehicle\", \"actuators\", \"thruster\", \"thruster_models_module.jl\"))", simulation_model_src) <
    findfirst("include(joinpath(@__DIR__, \"..\", \"dynamics\", \"coupled\", \"force_torque_models.jl\"))", simulation_model_src) ||
    error("SimulationModel is not loading ThrusterModels before DynamicEffectors.")
findfirst("include(joinpath(@__DIR__, \"..\", \"gnc\", \"guidance\", \"guidance_models.jl\"))", simulation_model_src) <
    findfirst("include(joinpath(@__DIR__, \"..\", \"dynamics\", \"coupled\", \"force_torque_models.jl\"))", simulation_model_src) ||
    error("SimulationModel is not loading GuidanceModels before DynamicEffectors.")
occursin("include(joinpath(@__DIR__, \"..\", \"environment\", \"gravity\", \"gravity_effectors.jl\"))", simulation_model_src) ||
    error("SimulationModel does not include the top-level GravityEffectors boundary module.")
occursin("module GravityEffectors", top_level_gravity_effectors_src) ||
    error("Top-level GravityEffectors boundary module is missing.")
findfirst("include(joinpath(@__DIR__, \"..\", \"dynamics\", \"coupled\", \"force_torque_models.jl\"))", simulation_model_src) <
    findfirst("include(joinpath(@__DIR__, \"..\", \"environment\", \"gravity\", \"gravity_effectors.jl\"))", simulation_model_src) ||
    error("SimulationModel is not loading the top-level GravityEffectors boundary after DynamicEffectors.")
occursin("using ..SimulationModel", engine_src) ||
    error("SimulationEngine is not loading SimulationModel as a sibling module.")
occursin("import ..RuntimeServices", engine_src) ||
    error("SimulationEngine is not loading RuntimeServices as a sibling module.")
!occursin("isdefined(parentmodule(@__MODULE__), :SimulationModel)", engine_src) ||
    error("SimulationEngine still carries the legacy parent-module include-order guard for SimulationModel.")
!occursin("isdefined(parentmodule(@__MODULE__), :RuntimeServices)", engine_src) ||
    error("SimulationEngine still carries the legacy parent-module include-order guard for RuntimeServices.")
!occursin("include(joinpath(@__DIR__, \"..\", \"..\", \"core\", \"simulation_model.jl\"))", engine_src) ||
    error("SimulationEngine still directly includes core/simulation_model.jl.")
!occursin("include(joinpath(@__DIR__, \"..\", \"runtime_services.jl\"))", engine_src) ||
    error("SimulationEngine still directly includes runtime_services.jl.")
!occursin("include(joinpath(@__DIR__, \"..\", \"..\", \"simulation\", \"runtime_services.jl\"))", engine_src) ||
    error("SimulationEngine still directly includes runtime_services.jl.")

occursin("using ...RuntimeServices: GRAM_LOCK", density_models_src) ||
    error("Environment density models are not importing GRAM_LOCK from RuntimeServices.")
occursin("using ...RuntimeServices: SPICE_LOCK", ephemerides_models_src) ||
    error("Ephemerides models are not importing SPICE_LOCK from RuntimeServices.")
occursin("using ...RuntimeServices: SPICE_LOCK, GRAM_LOCK", callback_registry_src) ||
    error("Callback registry is not importing SPICE_LOCK and GRAM_LOCK from RuntimeServices.")
occursin("RuntimeServices.SPICE_LOCK", setup_src) ||
    error("Simulation engine setup is not using RuntimeServices.SPICE_LOCK.")
!occursin("SimulationModel.SPICE_LOCK", density_models_src * ephemerides_models_src * callback_registry_src * setup_src) ||
    error("Package-owned source still references SimulationModel.SPICE_LOCK.")
!occursin("SimulationModel.GRAM_LOCK", density_models_src * ephemerides_models_src * callback_registry_src * setup_src) ||
    error("Package-owned source still references SimulationModel.GRAM_LOCK.")
occursin("using ..ReferenceSystems", guidance_hooks_src) ||
    error("GuidanceHooks is not importing the CamelCase ReferenceSystems module.")
occursin("using ..ReferenceSystems", control_hooks_src) ||
    error("ControlHooks is not importing the CamelCase ReferenceSystems module.")
!occursin("ref_sys", guidance_hooks_src * control_hooks_src * control_commands_src * constraint_tracking_src * trajectory_predictor_src * eom_predictor_src) ||
    error("Package-owned source still references the retired ref_sys module name.")
!occursin("Any[", density_callback_assembly_src) ||
    error("Density callback assembly still builds callback collections through Any[].")
occursin("_append_callback", density_callback_assembly_src) ||
    error("Density callback assembly is not using typed callback tuple composition.")
occursin("function run_simulation(\n    args::SimulationConfiguration;", engine_execution_src) ||
    error("SimulationEngine execution boundary is not typed on SimulationConfiguration.")
occursin("_require_simulation_configuration", engine_public_api_src) ||
    error("SimulationEngine public API is not validating untyped run_simulation inputs.")
occursin("function run_simulation(args; kwargs...)", engine_public_api_src) ||
    error("SimulationEngine public API no longer provides the generic validation wrapper.")
occursin("Base.depwarn", engine_public_api_src) ||
    error("SimulationEngine public API is missing the Base.depwarn shim for untyped run_simulation inputs.")
!occursin("function run_simulation(\n    args;", engine_execution_src) ||
    error("SimulationEngine execution.jl still owns the generic run_simulation(args; ...) boundary.")

!occursin("simulation\", \"execution\", \"run_simulation.jl", telemetry_src) ||
    error("TelemetryVerification still directly includes simulation/execution/run_simulation.jl.")
!occursin("simulation\", \"engine\", \"simulation_engine.jl", telemetry_src) ||
    error("TelemetryVerification still directly includes simulation/engine/simulation_engine.jl.")
!occursin("SimulationEngine.SimulationModel", telemetry_src) ||
    error("TelemetryVerification still references SimulationEngine.SimulationModel.")
!occursin("isdefined(parentmodule(@__MODULE__), :SimulationModel)", telemetry_src) ||
    error("TelemetryVerification still carries the legacy parent-module include-order guard for SimulationModel.")
!occursin("isdefined(parentmodule(@__MODULE__), :SimulationEngine)", telemetry_src) ||
    error("TelemetryVerification still carries the legacy parent-module include-order guard for SimulationEngine.")

occursin("SimulationEngine.run_simulation", telemetry_src) ||
    error("TelemetryVerification is not calling SimulationEngine.run_simulation.")

for include_name in (
    "gravity_effectors.jl",
    "aerodynamic_effectors.jl",
    "perturbation_effectors.jl",
    "thruster_models.jl",
    "guidance_models.jl",
)
    occursin("\"$include_name\"", force_torque_src) ||
        error("DynamicEffectors split contract violation: missing include for '$include_name'.")
end
occursin("module GravityEffectors", gravity_effectors_src) ||
    error("DynamicEffectors split contract violation: missing GravityEffectors module owner file.")
occursin("module AerodynamicEffectors", aerodynamic_effectors_src) ||
    error("DynamicEffectors split contract violation: missing AerodynamicEffectors module owner file.")
occursin("module PerturbationEffectors", perturbation_effectors_src) ||
    error("DynamicEffectors split contract violation: missing PerturbationEffectors module owner file.")
occursin("module ThrusterModels", thruster_models_src) ||
    error("DynamicEffectors split contract violation: missing ThrusterModels module owner file.")
occursin("module GuidanceModels", guidance_models_src) ||
    error("DynamicEffectors split contract violation: missing GuidanceModels module owner file.")
occursin("using ...ThrusterModels: BaseThrusterModel", thruster_models_src) ||
    error("DynamicEffectors split contract violation: ThrusterModels shim is not forwarding from the top-level ThrusterModels owner.")
occursin("using ...GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel", guidance_models_src) ||
    error("DynamicEffectors split contract violation: GuidanceModels shim is not forwarding from the top-level GuidanceModels owner.")
!occursin(r"^\s+module\s+(GravityEffectors|AerodynamicEffectors|PerturbationEffectors|ThrusterModels|GuidanceModels)"m, force_torque_src) ||
    error("DynamicEffectors split contract violation: child modules still live inline in force_torque_models.jl.")
occursin("export aerobraking_gravity_force_ii, srp, srp_cannonball_accel", force_torque_src) ||
    error("DynamicEffectors is not re-exporting the flat helper symbols expected by SimulationModel.")
occursin("using .AerodynamicEffectors: _threadid_capacity, _multibody_use_threads, _multibody_thread_decision", force_torque_src) ||
    error("DynamicEffectors is not forwarding the multibody parallel helpers needed by package-owned callers.")
occursin("using .PerturbationEffectors: srp, srp_cannonball_accel, _spice_query_name", force_torque_src) ||
    error("DynamicEffectors is not forwarding the perturbation helper symbols needed by package-owned callers.")

occursin("using ..GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel", guidance_hooks_src) ||
    error("GuidanceHooks is not importing the guidance model from the top-level GuidanceModels owner.")
occursin("using ..CommandTypes: PropulsiveManeuverCommand, AerobrakingControlCommand", guidance_hooks_src) ||
    error("GuidanceHooks is not importing command types from CommandTypes.")
occursin("using ..GravityEffectors: aerobraking_gravity_force_ii", guidance_hooks_src) ||
    error("GuidanceHooks is not importing aerobraking_gravity_force_ii from the top-level GravityEffectors boundary.")
!occursin("DynamicEffectors", guidance_hooks_src) ||
    error("GuidanceHooks still depends on DynamicEffectors directly.")
occursin("using ..ThrusterModels: BaseThrusterModel", control_hooks_src) ||
    error("ControlHooks is not importing BaseThrusterModel from the top-level ThrusterModels owner.")
occursin("using ..CommandTypes: PropulsiveManeuverCommand", control_hooks_src) ||
    error("ControlHooks is not importing PropulsiveManeuverCommand from CommandTypes.")
occursin("using ..GravityEffectors: aerobraking_gravity_force_ii", control_hooks_src) ||
    error("ControlHooks is not importing aerobraking_gravity_force_ii from the top-level GravityEffectors boundary.")
!occursin("DynamicEffectors", control_hooks_src) ||
    error("ControlHooks still depends on DynamicEffectors directly.")
occursin("using ..ThrusterModels: BaseThrusterModel", callback_registry_src) ||
    error("Callback registry is not importing BaseThrusterModel from the top-level ThrusterModels owner.")
!occursin("using ..DynamicEffectors.ThrusterModels: BaseThrusterModel", callback_registry_src) ||
    error("Callback registry still imports BaseThrusterModel from DynamicEffectors.ThrusterModels.")
occursin("using ..DynamicEffectors.AerodynamicEffectors: AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight", callback_registry_src) ||
    error("Callback registry is not importing aerodynamic coefficient models from DynamicEffectors.AerodynamicEffectors.")
occursin("using ..GravityEffectors: InverseSquaredJ2GravityModel", callback_registry_src) ||
    error("Callback registry is not importing gravity models from the top-level GravityEffectors boundary.")
!occursin("using ..DynamicEffectors.GravityEffectors: InverseSquaredJ2GravityModel", callback_registry_src) ||
    error("Callback registry still imports gravity models from DynamicEffectors.GravityEffectors.")

for retired_simulation_dir in (
    joinpath("src", "simulation", "events"),
    joinpath("src", "simulation", "execution"),
    joinpath("src", "simulation", "solver_orchestration"),
)
    isdir(joinpath(REPO_ROOT, retired_simulation_dir)) &&
        error("Canonical-ownership contract violation: retired simulation tree still exists at $retired_simulation_dir")
end

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
    joinpath("src", "simulation", "events"),
    joinpath("src", "simulation", "execution"),
    joinpath("src", "simulation", "solver_orchestration"),
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
