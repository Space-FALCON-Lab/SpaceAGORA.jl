# Migrated out of test/suites/01_contract_and_api_tests.jl (see
# test/TEST_RESTRUCTURE_PLAN.md, Phase 2). Unlike the other ci_*_gate.jl
# files in this directory, these check raw module/export/include-order
# structure via the sandbox modules and raw-include mechanics that
# test/helpers/bootstrap.jl sets up (EXPORT_IMPORT_SANDBOX, INCLUDE_ORDER_SANDBOX,
# REPO_ROOT, the standalone SimulationModel/SimulationEngine copies) --
# rewriting them to the self-contained `using SpaceAGORA` style the other
# gates use would change what's actually being verified, so this one
# includes bootstrap.jl instead.
include(joinpath(@__DIR__, "..", "helpers", "bootstrap.jl"))

@testset "SimulationModel Export Contract" begin
    sandbox = EXPORT_IMPORT_SANDBOX
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
    @test Core.eval(sandbox, :(isdefined(@__MODULE__, :RuntimeServices)))
    @test_nowarn Core.eval(sandbox, :(using .SimulationModel))

    required_public_names = [
        :SimulationConfiguration,
        :InitialCondition,
        :InitialTime,
        :MissionConfiguration,
        :EnvironmentModel,
        :SimpleEphemeridesModel,
        :PiecewiseExponentialAtmosphereModel,
        :NRLMSISE00AtmosphereModel,
        :SpiceEphemeridesModel,
        :make_no_gram_environment,
        :IntegrationTolerances,
        :ControlModel,
        :GuidanceModel,
        :NavigationModel,
        :DynamicsModel
    ]

    for sym in required_public_names
        @test Core.eval(sandbox, :(isdefined(@__MODULE__, $(QuoteNode(sym)))))
    end

    for sym in (
        :asim_ctrl,
        :asim_ctrl_plot,
        :control_solarpanels_heatrate,
        :control_solarpanels_heatload,
        :control_solarpanels_openloop,
    )
        @test !Core.eval(sandbox, :(isdefined(@__MODULE__, $(QuoteNode(sym)))))
        @test !Core.eval(sandbox, :(Base.isexported(SimulationModel.ControlHooks, $(QuoteNode(sym)))))
    end

    @test !Core.eval(sandbox, :(isdefined(SimulationModel, :SPICE_LOCK)))
    @test !Core.eval(sandbox, :(isdefined(SimulationModel, :GRAM_LOCK)))
    @test !Core.eval(sandbox, :(isdefined(SimulationModel, :ref_sys)))
    @test !Core.eval(sandbox, :(isdefined(SimulationModel, :PhysicalModel)))
    @test Core.eval(sandbox, :(isdefined(SimulationModel, :ReferenceSystems)))
    @test Core.eval(sandbox, :(isdefined(SimulationModel, :SpacecraftModels)))
    @test Core.eval(sandbox, :(isdefined(SimulationModel, :CommandTypes)))
    @test Core.eval(sandbox, :(isdefined(SimulationModel, :PropulsiveManeuverCommand)))
    @test Core.eval(sandbox, :(isdefined(SimulationModel, :LegacyGravityModelCode)))
    @test Core.eval(sandbox, :(isdefined(RuntimeServices, :SPICE_LOCK)))
    @test Core.eval(sandbox, :(isdefined(RuntimeServices, :GRAM_LOCK)))
end

@testset "Include-Order + Name Ambiguity Smoke" begin
    sandbox = INCLUDE_ORDER_SANDBOX

    Base.include_string(sandbox, """
    module ConflictingExports
        export SimulationConfiguration
        struct SimulationConfiguration end
    end
    using .ConflictingExports
    """)

    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
    Core.eval(sandbox, :(const quat_mult = SimulationModel.quat_mult))
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
    @test isdefined(sandbox, :RuntimeServices)
    @test isdefined(sandbox, :SimulationEngine)
    @test Core.eval(sandbox, :(isdefined(SimulationEngine, :run_simulation)))
end

@testset "SimulationModel Sibling-Module Topology Contract" begin
    spaceagora_src = read(joinpath(REPO_ROOT, "src", "SpaceAGORA.jl"), String)
    runtime_services_src = read(joinpath(REPO_ROOT, "src", "simulation", "runtime_services.jl"), String)
    simulation_model_src = read(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"), String)
    reference_system_config_src = read(joinpath(REPO_ROOT, "src", "core", "state", "reference_system_config.jl"), String)
    spacecraft_model_src = read(joinpath(REPO_ROOT, "src", "vehicle", "spacecraft", "model.jl"), String)
    legacy_model_codes_src = read(joinpath(REPO_ROOT, "src", "core", "types", "compat_model_codes.jl"), String)
    runtime_types_src = read(joinpath(REPO_ROOT, "src", "core", "types", "runtime_types.jl"), String)
    command_types_src = read(joinpath(REPO_ROOT, "src", "gnc", "command_types.jl"), String)
    engine_src = read(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"), String)
    engine_public_api_src = read(joinpath(REPO_ROOT, "src", "simulation", "engine", "public_api.jl"), String)
    engine_execution_src = read(joinpath(REPO_ROOT, "src", "simulation", "engine", "execution.jl"), String)
    telemetry_src = read(joinpath(REPO_ROOT, "src", "analysis", "verification", "telemetry_verification.jl"), String)
    density_callback_assembly_src = read(joinpath(REPO_ROOT, "src", "simulation", "callbacks", "density_callbacks", "assembly.jl"), String)
    legacy_nested_model_path = "SimulationEngine" * ".SimulationModel"

    @test occursin("include(joinpath(@__DIR__, \"simulation\", \"runtime_services.jl\"))", spaceagora_src)
    @test occursin("include(joinpath(@__DIR__, \"core\", \"simulation_model.jl\"))", spaceagora_src)
    @test findfirst("include(joinpath(@__DIR__, \"simulation\", \"runtime_services.jl\"))", spaceagora_src) <
        findfirst("include(joinpath(@__DIR__, \"core\", \"simulation_model.jl\"))", spaceagora_src)
    @test findfirst("include(joinpath(@__DIR__, \"core\", \"simulation_model.jl\"))", spaceagora_src) <
        findfirst("include(joinpath(@__DIR__, \"simulation\", \"engine\", \"simulation_engine.jl\"))", spaceagora_src)
    @test occursin("module RuntimeServices", runtime_services_src)
    @test occursin("const SPICE_LOCK = ReentrantLock()", runtime_services_src)
    # GRAM_LOCK must alias SPICE_LOCK, not be an independent ReentrantLock():
    # libGRAM.dylib's statically-linked CSPICE globally exports the same
    # internal symbol names as SpaceAGORA's own SPICE.jl bindings, so a
    # native-GRAM call and a SpaceAGORA ephemerides call must never run
    # concurrently under two different locks (see the comment on
    # RuntimeServices.GRAM_LOCK/SPICE_LOCK).
    @test occursin("const GRAM_LOCK = SPICE_LOCK", runtime_services_src)
    @test !occursin("const SPICE_LOCK", simulation_model_src)
    @test !occursin("const GRAM_LOCK", simulation_model_src)
    @test occursin("runtime_services.jl", simulation_model_src)
    @test occursin("include(joinpath(@__DIR__, \"..\", \"vehicle\", \"spacecraft\", \"model.jl\"))", simulation_model_src)
    @test occursin("@reexport using .SpacecraftModels", simulation_model_src)
    @test !occursin("@reexport using .PhysicalModel", simulation_model_src)
    @test occursin("include(joinpath(@__DIR__, \"..\", \"core\", \"types\", \"compat_model_codes.jl\"))", simulation_model_src)
    @test findfirst("include(joinpath(@__DIR__, \"..\", \"core\", \"types\", \"compat_model_codes.jl\"))", simulation_model_src) <
        findfirst("include(joinpath(@__DIR__, \"..\", \"core\", \"types\", \"runtime_types.jl\"))", simulation_model_src)
    @test occursin("module ReferenceSystems", reference_system_config_src)
    @test !occursin("module ref_sys", reference_system_config_src)
    @test occursin("module SpacecraftModels", spacecraft_model_src)
    @test !occursin("module PhysicalModel", spacecraft_model_src)
    @test occursin("module LegacyModelCodes", legacy_model_codes_src)
    @test occursin("using ..LegacyModelCodes:", runtime_types_src)
    @test occursin("using ..SpacecraftModels: SpacecraftModel", runtime_types_src)
    @test !occursin("@enum Legacy", runtime_types_src)
    @test occursin("include(joinpath(@__DIR__, \"..\", \"gnc\", \"command_types.jl\"))", simulation_model_src)
    @test findfirst("include(joinpath(@__DIR__, \"..\", \"gnc\", \"command_types.jl\"))", simulation_model_src) <
        findfirst("include(joinpath(@__DIR__, \"..\", \"core\", \"types\", \"runtime_types.jl\"))", simulation_model_src)
    @test occursin("module CommandTypes", command_types_src)
    @test occursin("using ..CommandTypes: PropulsiveManeuverCommand", runtime_types_src)
    @test occursin("maneuver_commands::Vector{PropulsiveManeuverCommand}", runtime_types_src)
    @test occursin("using ..SimulationModel", engine_src)
    @test occursin("import ..RuntimeServices", engine_src)
    @test !occursin("isdefined(parentmodule(@__MODULE__), :SimulationModel)", engine_src)
    @test !occursin("isdefined(parentmodule(@__MODULE__), :RuntimeServices)", engine_src)
    @test !occursin("include(joinpath(@__DIR__, \"..\", \"..\", \"core\", \"simulation_model.jl\"))", engine_src)
    @test !occursin("include(joinpath(@__DIR__, \"..\", \"runtime_services.jl\"))", engine_src)
    @test !occursin("include(joinpath(@__DIR__, \"..\", \"..\", \"simulation\", \"runtime_services.jl\"))", engine_src)
    @test occursin("_require_simulation_configuration", engine_public_api_src)
    @test occursin("function run_simulation(args; kwargs...)", engine_public_api_src)
    @test occursin("Base.depwarn", engine_public_api_src)
    @test occursin("function run_simulation(\n    args::SimulationConfiguration;", engine_execution_src)
    @test !occursin("function run_simulation(\n    args;", engine_execution_src)
    @test !occursin("Any[", density_callback_assembly_src)
    @test occursin("_append_callback", density_callback_assembly_src)
    @test !occursin(legacy_nested_model_path, telemetry_src)
    @test !occursin("isdefined(parentmodule(@__MODULE__), :SimulationModel)", telemetry_src)
    @test !occursin("isdefined(parentmodule(@__MODULE__), :SimulationEngine)", telemetry_src)
end

@testset "Typed Engine Boundary Contract" begin
    sandbox = Module(:SpaceAGORATypedEngineSandbox)
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "SpaceAGORA.jl"))
    @test Core.eval(sandbox, :(hasmethod(SpaceAGORA.SimulationEngine.run_simulation, Tuple{SpaceAGORA.SimulationModel.SimulationConfiguration})))
    @test_logs (:warn, r"typed SimulationConfiguration") @test_throws ArgumentError Core.eval(
        sandbox,
        :(SpaceAGORA.SimulationEngine.run_simulation(Dict(:bad => 1)))
    )
    @test_logs (:warn, r"typed SimulationConfiguration") @test_throws ArgumentError Core.eval(
        sandbox,
        :(SpaceAGORA.run_simulation(Dict(:bad => 1)))
    )
end

@testset "DynamicEffectors Internal Split Contract" begin
    force_torque_src = read(joinpath(REPO_ROOT, "src", "dynamics", "coupled", "force_torque_models.jl"), String)
    gravity_effectors_src = read(joinpath(REPO_ROOT, "src", "dynamics", "coupled", "force_torque_models", "gravity_effectors.jl"), String)
    aerodynamic_effectors_src = read(joinpath(REPO_ROOT, "src", "dynamics", "coupled", "force_torque_models", "aerodynamic_effectors.jl"), String)
    perturbation_effectors_src = read(joinpath(REPO_ROOT, "src", "dynamics", "coupled", "force_torque_models", "perturbation_effectors.jl"), String)
    thruster_models_src = read(joinpath(REPO_ROOT, "src", "dynamics", "coupled", "force_torque_models", "thruster_models.jl"), String)
    guidance_models_src = read(joinpath(REPO_ROOT, "src", "dynamics", "coupled", "force_torque_models", "guidance_models.jl"), String)
    top_level_thruster_models_src = read(joinpath(REPO_ROOT, "src", "vehicle", "actuators", "thruster", "thruster_models_module.jl"), String)
    top_level_guidance_models_src = read(joinpath(REPO_ROOT, "src", "gnc", "guidance", "guidance_models.jl"), String)
    top_level_gravity_effectors_src = read(joinpath(REPO_ROOT, "src", "environment", "gravity", "gravity_effectors.jl"), String)
    command_types_src = read(joinpath(REPO_ROOT, "src", "gnc", "command_types.jl"), String)
    guidance_hooks_src = read(joinpath(REPO_ROOT, "src", "gnc", "guidance", "guidance_hooks.jl"), String)
    control_hooks_src = read(joinpath(REPO_ROOT, "src", "gnc", "control", "control_hooks.jl"), String)
    callback_registry_src = read(joinpath(REPO_ROOT, "src", "simulation", "callbacks", "registry.jl"), String)

    for include_name in (
        "gravity_effectors.jl",
        "aerodynamic_effectors.jl",
        "perturbation_effectors.jl",
        "thruster_models.jl",
        "guidance_models.jl",
    )
        @test occursin("\"$include_name\"", force_torque_src)
    end
    @test occursin("module GravityEffectors", gravity_effectors_src)
    @test occursin("module AerodynamicEffectors", aerodynamic_effectors_src)
    @test occursin("module PerturbationEffectors", perturbation_effectors_src)
    @test occursin("module ThrusterModels", thruster_models_src)
    @test occursin("module GuidanceModels", guidance_models_src)
    @test occursin("module ThrusterModels", top_level_thruster_models_src)
    @test occursin("module GuidanceModels", top_level_guidance_models_src)
    @test occursin("module GravityEffectors", top_level_gravity_effectors_src)
    @test occursin("module CommandTypes", command_types_src)
    @test !occursin(r"^\s+module\s+(GravityEffectors|AerodynamicEffectors|PerturbationEffectors|ThrusterModels|GuidanceModels)"m, force_torque_src)

    @test occursin("function calcForceTorque end", force_torque_src)
    @test occursin("using ..EffectorSampling: wrench, wrench_caching!, environment_requirements, solver_partition", force_torque_src)
    @test occursin("using .GravityEffectors: ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel", force_torque_src)
    @test occursin("using .AerodynamicEffectors: AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight", force_torque_src)
    @test occursin("using .PerturbationEffectors: NBodyGravityModel, GravitationalHarmonicsModel, SolarRadiationPressureModel", force_torque_src)
    @test occursin("using .ThrusterModels: BaseThrusterModel", force_torque_src)
    @test occursin("using .GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel", force_torque_src)
    @test occursin("using ...ThrusterModels: BaseThrusterModel", thruster_models_src)
    @test occursin("using ...GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel", guidance_models_src)
    @test occursin("export aerobraking_gravity_force_ii, srp, srp_cannonball_accel", force_torque_src)

    @test occursin("using ..GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel", guidance_hooks_src)
    @test occursin("using ..CommandTypes: PropulsiveManeuverCommand, AerobrakingControlCommand", guidance_hooks_src)
    @test occursin("using ..GravityEffectors: aerobraking_gravity_force_ii", guidance_hooks_src)
    @test !occursin("DynamicEffectors", guidance_hooks_src)
    @test occursin("using ..ThrusterModels: BaseThrusterModel", control_hooks_src)
    @test occursin("using ..CommandTypes: PropulsiveManeuverCommand", control_hooks_src)
    @test occursin("using ..GravityEffectors: aerobraking_gravity_force_ii", control_hooks_src)
    @test !occursin("DynamicEffectors", control_hooks_src)
    @test occursin("using ..ThrusterModels: BaseThrusterModel", callback_registry_src)
    @test occursin("using ..DynamicEffectors.AerodynamicEffectors: AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight", callback_registry_src)
    @test occursin("using ..GravityEffectors: InverseSquaredJ2GravityModel", callback_registry_src)

    @test isdefined(SimulationModel, :ConstantGravityModel)
    @test isdefined(SimulationModel, :NBodyGravityModel)
    @test isdefined(SimulationModel, :AerodynamicCoefficientConstant)
    @test isdefined(SimulationModel, :BaseThrusterModel)
    @test isdefined(SimulationModel, :AerobrakingCampaignPropulsiveManeuverGuidanceModel)
    @test isdefined(SimulationModel, :PropulsiveManeuverCommand)
    @test isdefined(SimulationModel, :StateSample)
    @test isdefined(SimulationModel, :EnvironmentSample)
    @test isdefined(SimulationModel, :wrench)
    @test isdefined(SimulationModel, :environment_requirements)
    @test isdefined(SimulationModel, :solver_partition)
    @test isdefined(SimulationModel, :gravity_backbone_kick_structure)
    @test isdefined(SimulationModel, :gravity_backbone_kick_acceleration_ii)
    @test isdefined(SimulationModel, :aerobraking_gravity_force_ii)
    @test isdefined(SimulationModel.DynamicEffectors, :_spice_query_name)
    @test isdefined(SimulationModel.DynamicEffectors, :_parse_bool_env)
    @test isdefined(SimulationModel.DynamicEffectors, :_multibody_outer_parallel_hint)
    @test isdefined(SimulationModel.DynamicEffectors, :_multibody_use_threads)
    @test isdefined(SimulationModel.DynamicEffectors, :collect_and_reset_link_wrenches!)
    @test isdefined(SimulationModel.DynamicEffectors, :_aero_workspace_for_sat!)
    @test isdefined(SimulationModel.DynamicEffectors, :_nbody_workspace_for_sat!)
    @test isdefined(SimulationModel.DynamicEffectors, :_harmonics_workspace_for_sat!)
    @test isdefined(SimulationModel.DynamicEffectors, :_nbody_body_position_from_cache_j2000_m)
    @test isdefined(SimulationModel.DynamicEffectors, :_srp_sun_position_from_cache_j2000_m)
    @test isdefined(SimulationModel.DynamicEffectors, :eclipse_area_calc)
end

@testset "Simulation Filename Canonical Contract" begin
    engine_path = joinpath(REPO_ROOT, "src", "simulation", "engine", "execution.jl")
    legacy_events_dir = joinpath(REPO_ROOT, "src", "simulation", "events")
    legacy_exec_dir = joinpath(REPO_ROOT, "src", "simulation", "execution")
    legacy_solver_dir = joinpath(REPO_ROOT, "src", "simulation", "solver_orchestration")
    legacy_aerobraking_path = joinpath(REPO_ROOT, "src", "simulation", "Aerobraking.jl")
    legacy_complete_path = joinpath(REPO_ROOT, "src", "simulation", "Complete_passage.jl")
    legacy_mission_model_path = joinpath(REPO_ROOT, "src", "mission", "mission_model.jl")
    legacy_define_mission_path = joinpath(REPO_ROOT, "src", "mission", "define_mission.jl")
    legacy_save_results_path = joinpath(REPO_ROOT, "src", "io", "outputs", "save_results.jl")
    legacy_mc_set_path = joinpath(REPO_ROOT, "src", "mission", "campaigns", "montecarlo_set.jl")
    legacy_mc_perturbations_path = joinpath(REPO_ROOT, "src", "mission", "campaigns", "montecarlo_perturbations.jl")

    @test isfile(engine_path)
    @test !isdir(legacy_events_dir)
    @test !isdir(legacy_exec_dir)
    @test !isdir(legacy_solver_dir)
    @test !isfile(legacy_aerobraking_path)
    @test !isfile(legacy_complete_path)
    @test !isfile(legacy_mission_model_path)
    @test !isfile(legacy_define_mission_path)
    @test !isfile(legacy_save_results_path)
    @test !isfile(legacy_mc_set_path)
    @test !isfile(legacy_mc_perturbations_path)
end

@testset "SharedBuffers Type Contract" begin
    shared_buffers_type = SimulationModel.ConfigTypes.SharedBuffers

    @test fieldtype(shared_buffers_type, :density_models) ==
        Vector{Union{SimulationModel.GRAMAtmosphereModel, SimulationModel.GRAMAtmosphereModelSurrogate}}
    @test fieldtype(shared_buffers_type, :gram_density_cache) ==
        Vector{Union{Nothing, SimulationModel.GramTrackCache}}
    @test fieldtype(shared_buffers_type, :gram_isolated_pool_models) ==
        Vector{SimulationModel.GRAMAtmosphereModel}
    @test fieldtype(shared_buffers_type, :harmonics_workspaces) ==
        Vector{Union{Nothing, Dict{UInt, SimulationModel.HarmonicsScratchWorkspace}}}
    @test fieldtype(shared_buffers_type, :nbody_workspaces) ==
        Vector{Union{Nothing, SimulationModel.NBodyScratchWorkspace}}
    @test fieldtype(shared_buffers_type, :aero_workspaces) ==
        Vector{Union{Nothing, SimulationModel.AeroScratchWorkspace}}
    @test fieldtype(shared_buffers_type, :maneuver_commands) ==
        Vector{SimulationModel.PropulsiveManeuverCommand}
    @test SimulationModel.SaveData == Dict{Symbol, Any}
end

@testset "SpaceAGORA run_simulation Doc Contract" begin
    sandbox = Module(:SpaceAGORADocSandbox)
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "SpaceAGORA.jl"))
    doc = Core.eval(sandbox, :(Base.Docs.doc(SpaceAGORA.run_simulation)))
    doc_text = sprint(show, MIME"text/plain"(), doc)

    @test occursin("isolate_state", doc_text)
    @test occursin("deep-copies the simulation configuration", doc_text)
    @test occursin("concurrent", doc_text)
    @test occursin("mutate shared state", doc_text)
end

@testset "SpaceAGORA Public Abstract Interface Contract" begin
    sandbox = Module(:SpaceAGORAAbstractSandbox)
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "SpaceAGORA.jl"))

    for sym in (
        :AbstractForceTorqueModel,
        :AbstractDensityModel,
        :AbstractControlEffectorModel,
        :AbstractEphemeridesModel,
        :AbstractThermalModel,
        :AbstractThrusterModel,
        :AbstractGuidanceModel,
    )
        @test Core.eval(sandbox, :(Base.isexported(SpaceAGORA, $(QuoteNode(sym)))))
        doc = Core.eval(sandbox, :(Base.Docs.doc(getproperty(SpaceAGORA, $(QuoteNode(sym))))))
        @test doc !== nothing
    end

    @test !Core.eval(sandbox, :(Base.isexported(SpaceAGORA, :SimulationModel)))
    @test !Core.eval(sandbox, :(Base.isexported(SpaceAGORA, :RuntimeServices)))
end

@testset "Typed Effector Sampling Public Contract" begin
    sandbox = Module(:SpaceAGORAWrenchSandbox)
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "SpaceAGORA.jl"))

    for sym in (
        :StateSample,
        :PlanetFrameSample,
        :AtmosphereSample,
        :SolarEphemerisSample,
        :ThirdBodyEphemerisSample,
        :EnvironmentSample,
        :EffectorEnvironmentRequirements,
        :wrench,
        :environment_requirements,
        :solver_partition,
        :gravity_backbone_structure,
        :gravity_backbone_acceleration_ii,
        :gravity_backbone_kick_structure,
        :gravity_backbone_kick_acceleration_ii,
    )
        @test Core.eval(sandbox, :(Base.isexported(SpaceAGORA, $(QuoteNode(sym)))))
        doc = Core.eval(sandbox, :(Base.Docs.doc(getproperty(SpaceAGORA, $(QuoteNode(sym))))))
        @test doc !== nothing
    end
end

@testset "Documenter Strictness Contract" begin
    docs_make = read(joinpath(REPO_ROOT, "docs", "make.jl"), String)
    getting_started = read(joinpath(REPO_ROOT, "docs", "src", "getting_started.md"), String)
    from_env_src = read(joinpath(REPO_ROOT, "src", "simulation", "engine", "adapters", "from_env.jl"), String)

    @test occursin("modules = [SpaceAGORA]", docs_make)
    @test occursin("doctest = true", docs_make)
    @test occursin("checkdocs = :exports", docs_make)
    @test occursin("checkdocs_ignored_modules = Module[", docs_make)
    @test occursin("SpaceAGORA.RuntimeServices", docs_make)
    @test occursin("SpaceAGORA.SimulationEngine", docs_make)
    @test occursin("SpaceAGORA.SimulationCampaigns", docs_make)
    @test occursin("SpaceAGORA.SimulationModel", docs_make)
    @test occursin("SpaceAGORA.ParallelProfiles", docs_make)
    @test occursin("SpaceAGORA.TelemetryVerification", docs_make)
    @test occursin("SpaceAGORA.SpaceAGORACLI", docs_make)
    @test occursin("warnonly = false", docs_make)
    @test occursin("spaceagora_no_gram_example_args", docs_make)
    @test occursin("```jldoctest", getting_started)
    @test occursin("```jldoctest", from_env_src)
end

@testset "Monte Carlo Campaign Runner Contract" begin
    sandbox = Module(:SpaceAGORAMonteCarloSandbox)
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "SpaceAGORA.jl"))

    @test Core.eval(sandbox, :(Base.isexported(SpaceAGORA, :MonteCarloSpec)))
    @test Core.eval(sandbox, :(Base.isexported(SpaceAGORA, :MonteCarloSampleResult)))
    @test Core.eval(sandbox, :(Base.isexported(SpaceAGORA, :MonteCarloResult)))
    @test Core.eval(sandbox, :(Base.isexported(SpaceAGORA, :run_monte_carlo)))

    result = Core.eval(sandbox, quote
        SpaceAGORA.run_monte_carlo(10:12; threads=1) do seed
            seed + 1
        end
    end)
    @test result.threads == 1
    @test length(result.samples) == 3
    @test length(result.successful) == 3
    @test isempty(result.failed)
    @test [sample.seed for sample in result.samples] == [10, 11, 12]
    @test [sample.value for sample in result.samples] == [11, 12, 13]

    mixed = Core.eval(sandbox, quote
        spec = SpaceAGORA.MonteCarloSpec(seeds=1:3, threads=1)
        SpaceAGORA.run_monte_carlo(spec) do seed
            seed == 2 && error("seed two failed")
            seed
        end
    end)
    @test length(mixed.samples) == 3
    @test length(mixed.successful) == 2
    @test length(mixed.failed) == 1
    @test mixed.failed[1].seed == 2

    @test_throws ArgumentError Core.eval(sandbox, :(
        SpaceAGORA.run_monte_carlo(identity, 1:1; threads=Threads.nthreads() + 1)
    ))
end

@testset "Aerobraking Selector Contract" begin
    selector = SimulationModel.DefaultAerobrakingPolicySelector()
    cfg_default = SimulationModel.AerobrakingPolicyConfig()
    cfg_targeting = SimulationModel.AerobrakingPolicyConfig(default_strategy=SimulationModel.T_EDG)
    input = SimulationModel.GuidanceHooks.AerobrakingGuidanceInput(
        ip=nothing,
        mission=nothing,
        args=Dict{Symbol, Any}(:aerobraking_strategy => "e_edg"),
    )

    @test SimulationModel.select_strategy(selector, cfg_default, input) == SimulationModel.E_EDG
    @test SimulationModel.select_strategy(selector, cfg_targeting, input) == SimulationModel.T_EDG
end

@testset "Root Environment Contract" begin
    agora_project = joinpath(REPO_ROOT, "Project.toml")
    agora_manifest = joinpath(REPO_ROOT, "Manifest.toml")

    @test isfile(agora_project)
    @test isfile(agora_manifest)

    tracked = read(`sh -lc "cd '$REPO_ROOT' && git ls-files Project.toml Manifest.toml"`, String)
    @test occursin("Project.toml", tracked)
    @test occursin("Manifest.toml", tracked)

    readme = read(joinpath(REPO_ROOT, "README.md"), String)
    getting_started = read(joinpath(REPO_ROOT, "docs", "src", "getting_started.md"), String)
    @test occursin("canonical committed execution environment", readme)
    @test occursin("canonical committed execution environment", getting_started)
    @test occursin("there is no bootstrap step", lowercase(readme))
end

@testset "Project Compat Coverage Contract" begin
    stdlibs = Set([
        "Artifacts",
        "Base64",
        "Dates",
        "DelimitedFiles",
        "Distributed",
        "FileWatching",
        "InteractiveUtils",
        "LibCURL",
        "LibGit2",
        "Libdl",
        "LinearAlgebra",
        "Logging",
        "Markdown",
        "Mmap",
        "Pkg",
        "Printf",
        "Profile",
        "Random",
        "REPL",
        "SHA",
        "Serialization",
        "SharedArrays",
        "Sockets",
        "SparseArrays",
        "Statistics",
        "SuiteSparse",
        "TOML",
        "Tar",
        "Test",
        "UUIDs",
        "Unicode"
    ])

    function nonstdlib_deps_without_compat(project_path::String)
        project = TOML.parsefile(project_path)
        deps = Set(String.(keys(get(project, "deps", Dict()))))
        compat = Set(String.(keys(get(project, "compat", Dict()))))
        delete!(compat, "julia")
        return sort(collect(setdiff(setdiff(deps, stdlibs), compat)))
    end

    @test isempty(nonstdlib_deps_without_compat(joinpath(REPO_ROOT, "Project.toml")))
end

@testset "Precompilation Contract" begin
    root_project = TOML.parsefile(joinpath(REPO_ROOT, "Project.toml"))
    spaceagora_src = read(joinpath(REPO_ROOT, "src", "SpaceAGORA.jl"), String)
    precompile_src = read(joinpath(REPO_ROOT, "src", "precompile_workload.jl"), String)
    clean_depot_smoke = read(joinpath(REPO_ROOT, "test", "ci_clean_depot_smoke.jl"), String)
    legacy_nested_model_path = "SimulationEngine" * ".SimulationModel"

    @test get(get(root_project, "deps", Dict()), "PrecompileTools", nothing) == "aea7be01-6a6a-4083-8856-8a6e6704d82a"
    @test get(get(root_project, "compat", Dict()), "PrecompileTools", nothing) == "1"

    @test occursin("__precompile__(true)", spaceagora_src)
    @test !occursin("__precompile__(false)", spaceagora_src)
    @test occursin("using PrecompileTools: @compile_workload, @setup_workload", spaceagora_src)
    @test occursin("include(joinpath(@__DIR__, \"precompile_workload.jl\"))", spaceagora_src)

    @test occursin("@compile_workload", precompile_src)
    @test occursin("parse_parallel_profile(\"R2\")", precompile_src)
    @test occursin("simulation_engine_config_from_env", precompile_src)
    @test occursin("run_simulation(engine_config, args; return_solution=true)", precompile_src)
    @test occursin("SimpleEphemeridesModel()", precompile_src)
    @test occursin("ExponentialAtmosphereModel(planet)", precompile_src)
    @test !occursin(legacy_nested_model_path, precompile_src)
    @test isdefined(SpaceAGORA, :_SPACEAGORA_PRECOMPILE_ENV)
    @test isdefined(SpaceAGORA, :_spaceagora_precompile_args)
    @test isdefined(SpaceAGORA, :_run_spaceagora_precompile_workload)

    precompile_args = SpaceAGORA._spaceagora_precompile_args()
    @test nameof(typeof(precompile_args.environment_model.density_model)) == :ExponentialAtmosphereModel
    @test nameof(typeof(precompile_args.environment_model.ephemerides_model)) == :SimpleEphemeridesModel
    @test precompile_args.simulation_settings.results === false
    @test SpaceAGORA._SPACEAGORA_PRECOMPILE_ENV["SPACEAGORA_PARALLEL_PROFILE"] == "R2"

    precompile_solution = SpaceAGORA._run_spaceagora_precompile_workload(workspace=tempdir())
    @test precompile_solution !== nothing
    @test last(precompile_solution.t) > first(precompile_solution.t)

    @test occursin("using SpaceAGORA", clean_depot_smoke)
    @test !occursin("include(joinpath(REPO_ROOT, \"src\", \"core\", \"simulation_model.jl\"))", clean_depot_smoke)
    @test !occursin("include(joinpath(REPO_ROOT, \"src\", \"simulation\", \"engine\", \"simulation_engine.jl\"))", clean_depot_smoke)
end

@testset "Dependabot Contract" begin
    dependabot_path = joinpath(REPO_ROOT, ".github", "dependabot.yml")
    @test isfile(dependabot_path)
    dependabot_src = read(dependabot_path, String)
    @test occursin("version: 2", dependabot_src)
    @test occursin("package-ecosystem: \"github-actions\"", dependabot_src)
    @test occursin("package-ecosystem: \"julia\"", dependabot_src)
    @test occursin("directory: \"/\"", dependabot_src)
    @test !occursin("directory: \"/.AGORA\"", dependabot_src)
    @test occursin("directory: \"/docs\"", dependabot_src)
end











# ensure_guidance_sandbox_loaded! moved to test/helpers/sandbox_modules.jl so
# it's available regardless of include order (test/unit/mission/ also needs
# it, and used to run before this file did -- see test/TEST_RESTRUCTURE_PLAN.md
# Phase 2 notes).









