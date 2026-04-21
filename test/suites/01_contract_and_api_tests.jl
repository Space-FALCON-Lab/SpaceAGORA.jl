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
    @test occursin("const GRAM_LOCK = ReentrantLock()", runtime_services_src)
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
    @test occursin("using ..EffectorSampling: wrench, environment_requirements, solver_partition", force_torque_src)
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
    shared_buffers_type = SimulationModel.ConfigTypes.SharedBuffers{2}

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
    @test occursin("SpaceAGORA.SimulationModel", docs_make)
    @test occursin("SpaceAGORA.ParallelProfiles", docs_make)
    @test occursin("SpaceAGORA.TelemetryVerification", docs_make)
    @test occursin("SpaceAGORA.SpaceAGORACLI", docs_make)
    @test occursin("warnonly = false", docs_make)
    @test occursin("spaceagora_no_gram_example_args", docs_make)
    @test occursin("```jldoctest", getting_started)
    @test occursin("```jldoctest", from_env_src)
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

@testset ".AGORA Environment Contract" begin
    agora_project = joinpath(REPO_ROOT, ".AGORA", "Project.toml")
    agora_manifest = joinpath(REPO_ROOT, ".AGORA", "Manifest.toml")

    @test isfile(agora_project)
    @test isfile(agora_manifest)

    tracked = read(`sh -lc "cd '$REPO_ROOT' && git ls-files .AGORA/Project.toml .AGORA/Manifest.toml"`, String)
    @test occursin(".AGORA/Project.toml", tracked)
    @test occursin(".AGORA/Manifest.toml", tracked)

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
    @test isempty(nonstdlib_deps_without_compat(joinpath(REPO_ROOT, ".AGORA", "Project.toml")))
end

@testset "Precompilation Contract" begin
    root_project = TOML.parsefile(joinpath(REPO_ROOT, "Project.toml"))
    agora_project = TOML.parsefile(joinpath(REPO_ROOT, ".AGORA", "Project.toml"))
    spaceagora_src = read(joinpath(REPO_ROOT, "src", "SpaceAGORA.jl"), String)
    precompile_src = read(joinpath(REPO_ROOT, "src", "precompile_workload.jl"), String)
    clean_depot_smoke = read(joinpath(REPO_ROOT, "test", "ci_clean_depot_smoke.jl"), String)
    legacy_nested_model_path = "SimulationEngine" * ".SimulationModel"

    @test get(get(root_project, "deps", Dict()), "PrecompileTools", nothing) == "aea7be01-6a6a-4083-8856-8a6e6704d82a"
    @test get(get(root_project, "compat", Dict()), "PrecompileTools", nothing) == "1"
    @test get(get(agora_project, "deps", Dict()), "PrecompileTools", nothing) == "aea7be01-6a6a-4083-8856-8a6e6704d82a"
    @test get(get(agora_project, "compat", Dict()), "PrecompileTools", nothing) == "1"

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
    @test occursin("directory: \"/.AGORA\"", dependabot_src)
    @test occursin("directory: \"/docs\"", dependabot_src)
end











function ensure_guidance_sandbox_loaded!()
    if GUIDANCE_SANDBOX_LOADED[]
        return
    end

    Core.eval(GUIDANCE_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "gnc", "guidance", "thruster_guidance", "thruster_guidance_models.jl"))))
    Core.eval(GUIDANCE_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "gnc", "guidance", "thruster_guidance", "thruster_guidance_functions.jl"))))
    GUIDANCE_SANDBOX_LOADED[] = true
end









@testset "API Convenience Constructors" begin
    @testset "Mission/Environment Config Validation" begin
        mc_default = MissionConfiguration()
        @test mc_default.mission_type == MissionTime

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            SimulationModel.SimConfig._deprecated_mission_type_input_warned[] = false
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            mc_str = MissionConfiguration(mission_type="Time", mission_time=600.0, number_of_orbits=1, num_steps_to_save=10)
            @test mc_str.mission_type == MissionTime
            @test mc_str.mission_type == "Time"
            @test "Time" == mc_str.mission_type
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            mc_sym = MissionConfiguration(mission_type=:orbits, mission_time=600.0, number_of_orbits=2, num_steps_to_save=10)
            @test mc_sym.mission_type == MissionOrbits
            @test mc_sym.mission_type == :orbits
            @test :orbits == mc_sym.mission_type
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            mc_enum = MissionConfiguration(mission_type=MissionOrbits, mission_time=600.0, number_of_orbits=3, num_steps_to_save=10)
            @test mc_enum.mission_type == MissionOrbits
            @test mc_enum.mission_type == "Orbits"
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "1") do
            SimulationModel.SimConfig._deprecated_mission_type_input_warned[] = false
            @test_logs (:warn, r"mission_type=.*deprecated") MissionConfiguration(mission_type="Time")
            @test SimulationModel.SimConfig._deprecated_mission_type_input_warned[] == true
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            SimulationModel.SimConfig._deprecated_mission_type_input_warned[] = false
            @test_logs MissionConfiguration(mission_type="Time")
            @test SimulationModel.SimConfig._deprecated_mission_type_input_warned[] == false
        end

        @test_throws ArgumentError MissionConfiguration(mission_type="invalid")
        @test_throws ArgumentError MissionConfiguration(mission_time=0.0)
        @test_throws ArgumentError MissionConfiguration(number_of_orbits=0)
        @test_throws ArgumentError MissionConfiguration(num_steps_to_save=0)

        env_ok = EnvironmentModel(
            planet=EARTH,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SpiceEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=8,
            topo_order=8,
            wind=false
        )
        @test env_ok.EI == 120.0
        @test env_ok.ephemerides_model isa SpiceEphemeridesModel

        env_no_gram_default = make_no_gram_environment()
        @test env_no_gram_default.ephemerides_model isa SimpleEphemeridesModel
        @test env_no_gram_default.density_model isa NoAtmosphereModel
        @test env_no_gram_default.planet isa Earth

        env_no_gram = make_no_gram_environment(planet=:mars, atmosphere=:exponential, EI_km=140.0)
        @test env_no_gram.ephemerides_model isa SimpleEphemeridesModel
        @test env_no_gram.density_model isa ExponentialAtmosphereModel

        exp_model = ExponentialAtmosphereModel(1.0e-4, 120.0e3, 12.0e3)
        @test exp_model.temperature_k == 200.0
        @test exp_model.valid_min_altitude_m == 120.0e3
        @test exp_model.valid_max_altitude_m == 180.0e3

        nrl_model = NRLMSISE00AtmosphereModel()
        @test nrl_model.f107a == 150.0
        @test nrl_model.f107 == 150.0
        @test nrl_model.ap == 4.0
        @test nrl_model.valid_min_altitude_m == 0.0
        @test nrl_model.valid_max_altitude_m == 1_000.0e3

        dt_nrl = DateTime(2024, 1, 1, 0, 0, 0)
        j2000_dt = DateTime(2000, 1, 1, 12, 0, 0)
        el_time_nrl = Dates.value(dt_nrl - j2000_dt) / 1000.0
        p_nrl = (
            args=(
                initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
                environment_model=(planet=EARTH,),
            ),
        )

        fixed_nrl = NRLMSISE00AtmosphereModel(f107a=120.0, f107=130.0, ap=6.0)
        expected_nrl_fixed = SatelliteToolbox.AtmosphericModels.nrlmsise00(
            dt_nrl, 400.0e3, 0.1, 0.2, 120.0, 130.0, 6.0
        )
        rho_nrl_fixed, T_nrl_fixed, wind_nrl_fixed = getDensity(fixed_nrl, 400.0e3, 0.1, 0.2, el_time_nrl, false)
        @test isapprox(rho_nrl_fixed, expected_nrl_fixed.total_density; atol=0.0, rtol=1e-12)
        @test isapprox(T_nrl_fixed, expected_nrl_fixed.temperature; atol=0.0, rtol=1e-12)
        @test wind_nrl_fixed == SVector{3, Float64}(0.0, 0.0, 0.0)

        rho_nrl_rel, T_nrl_rel, wind_nrl_rel = getDensity(fixed_nrl, 400.0e3, 0.1, 0.2, 0.0, false, p_nrl)
        @test isapprox(rho_nrl_rel, expected_nrl_fixed.total_density; atol=0.0, rtol=1e-12)
        @test isapprox(T_nrl_rel, expected_nrl_fixed.temperature; atol=0.0, rtol=1e-12)
        @test wind_nrl_rel == SVector{3, Float64}(0.0, 0.0, 0.0)

        provider_hits = Ref(0)
        provider_nrl = NRLMSISE00AtmosphereModel(
            index_provider=(instant, h, lat, lon) -> begin
                provider_hits[] += 1
                return (f107a=95.0, f107=105.0, ap=[8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0])
            end
        )
        expected_nrl_provider = SatelliteToolbox.AtmosphericModels.nrlmsise00(
            dt_nrl, 400.0e3, 0.1, 0.2, 95.0, 105.0, [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0]
        )
        rho_nrl_provider, T_nrl_provider, wind_nrl_provider = getDensity(provider_nrl, 400.0e3, 0.1, 0.2, el_time_nrl, false)
        @test provider_hits[] == 1
        @test isapprox(rho_nrl_provider, expected_nrl_provider.total_density; atol=0.0, rtol=1e-12)
        @test isapprox(T_nrl_provider, expected_nrl_provider.temperature; atol=0.0, rtol=1e-12)
        @test wind_nrl_provider == SVector{3, Float64}(0.0, 0.0, 0.0)

        space_indices_nrl = NRLMSISE00AtmosphereModel(use_space_indices=true)
        @test space_indices_nrl.index_provider isa SimulationModel.EnvironmentModels.NRLMSISE00SpaceIndicesProvider
        low_altitude_indices = space_indices_nrl.index_provider(dt_nrl, 70.0e3, 0.1, 0.2)
        @test low_altitude_indices == (f107a=150.0, f107=150.0, ap=4.0)

        helper_calls = Tuple{Symbol, DateTime}[]
        helper_instant = DateTime(2024, 1, 2, 10, 30, 0)
        function fake_space_indices_lookup(index, instant::DateTime)
            index_sym = if index isa Val{:F10adj_avg_center81}
                :F10adj_avg_center81
            elseif index isa Val{:F10adj}
                :F10adj
            elseif index isa Val{:Ap_daily}
                :Ap_daily
            elseif index isa Val{:Ap}
                :Ap
            else
                error("Unexpected index lookup in test: $index")
            end
            push!(helper_calls, (index_sym, instant))
            if index_sym === :F10adj_avg_center81
                return 177.0
            elseif index_sym === :F10adj
                return 166.0
            elseif index_sym === :Ap_daily
                return 99
            end
            day = Date(instant)
            if day == Date(2024, 1, 2)
                return (21, 22, 23, 24, 25, 26, 27, 28)
            elseif day == Date(2024, 1, 1)
                return (11, 12, 13, 14, 15, 16, 17, 18)
            elseif day == Date(2023, 12, 31)
                return (1, 2, 3, 4, 5, 6, 7, 8)
            end
            error("Unexpected Ap lookup date in test: $day")
        end

        helper_indices = SimulationModel.EnvironmentModels._nrlmsise_space_indices_indices(
            fake_space_indices_lookup,
            helper_instant
        )
        @test helper_indices.f107a == 177.0
        @test helper_indices.f107 == 166.0
        @test helper_indices.ap == SVector{7, Float64}(99.0, 24.0, 23.0, 22.0, 21.0, 14.5, 4.5)
        @test (Symbol(:F10adj_avg_center81), helper_instant) in helper_calls
        @test (Symbol(:F10adj), helper_instant - Day(1)) in helper_calls
        @test SimulationModel.EnvironmentModels._nrlmsise_ap_slot_index(helper_instant) == 4

        @test_throws ArgumentError NRLMSISE00AtmosphereModel(ap=[1.0, 2.0])
        @test_throws ArgumentError NRLMSISE00AtmosphereModel(
            use_space_indices=true,
            index_provider=(instant) -> (f107a=1.0, f107=1.0, ap=1.0)
        )
        @test_throws ArgumentError NRLMSISE00AtmosphereModel(space_indices_force_download=true)
        @test_throws ArgumentError NRLMSISE00AtmosphereModel(valid_min_altitude_m=1.0, valid_max_altitude_m=0.0)

        ρ_mid = 1.225 * exp(-20.0e3 / 8.5e3)
        piecewise = PiecewiseExponentialAtmosphereModel(
            [0.0, 20.0e3, 60.0e3],
            [1.225, ρ_mid],
            [8.5e3, 12.0e3];
            temperature_k=210.0
        )
        @test piecewise.valid_min_altitude_m == 0.0
        @test piecewise.valid_max_altitude_m == 60.0e3

        ρ_piece_low, T_piece_low, wind_piece_low = getDensity(piecewise, 10.0e3, 0.0, 0.0, 0.0, false)
        @test isapprox(ρ_piece_low, 1.225 * exp(-10.0e3 / 8.5e3); atol=0.0, rtol=1e-12)
        @test T_piece_low == 210.0
        @test wind_piece_low == SVector{3, Float64}(0.0, 0.0, 0.0)

        ρ_piece_high, T_piece_high, wind_piece_high = getDensity(piecewise, 30.0e3, 0.0, 0.0, 0.0, false)
        @test isapprox(ρ_piece_high, ρ_mid * exp(-(30.0e3 - 20.0e3) / 12.0e3); atol=0.0, rtol=1e-12)
        @test T_piece_high == 210.0
        @test wind_piece_high == SVector{3, Float64}(0.0, 0.0, 0.0)

        rhos_piece = zeros(3)
        Ts_piece = zeros(3)
        winds_piece = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:3]
        getDensityBatch!(
            rhos_piece,
            Ts_piece,
            winds_piece,
            piecewise,
            [10.0e3, 30.0e3, 80.0e3],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.0,
            false,
            nothing
        )
        @test isapprox(rhos_piece[1], ρ_piece_low; atol=0.0, rtol=1e-12)
        @test isapprox(rhos_piece[2], ρ_piece_high; atol=0.0, rtol=1e-12)
        @test isapprox(rhos_piece[3], ρ_mid * exp(-(80.0e3 - 20.0e3) / 12.0e3); atol=0.0, rtol=1e-12)
        @test all(==(210.0), Ts_piece)
        @test all(==(SVector{3, Float64}(0.0, 0.0, 0.0)), winds_piece)

        @test_throws ArgumentError PiecewiseExponentialAtmosphereModel([0.0], [1.0], [8.5e3])
        @test_throws ArgumentError PiecewiseExponentialAtmosphereModel([0.0, 20.0e3, 10.0e3], [1.0, 0.1], [8.5e3, 7.0e3])

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=-1.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SpiceEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=false,
            topo_degree=8,
            topo_order=8,
            wind=false
        )

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SpiceEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=-1,
            topo_order=8,
            wind=false
        )

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SpiceEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=8,
            topo_order=-1,
            wind=false
        )
    end

    ic = InitialCondition()
    @test ic isa InitialCondition
    @test ic.a == 0.0
    @test ic.e == 0.0

    link = Link()
    @test link isa Link{0}

    joint = Joint()
    @test joint isa Joint

    sc = SpacecraftModel()
    @test sc isa SpacecraftModel
    @test sc.root.root

    custom_root = Link{0}(root=true, m=100.0)
    sc_custom = SpacecraftModel(root=custom_root, id=42)
    @test sc_custom.id == 42
    @test sc_custom.root === custom_root
    @test any(link -> link === custom_root, sc_custom.links)

    nbody = NBodyGravityModel(["Sun"], "Earth", SPICE_PATH)
    @test nbody isa NBodyGravityModel
    @test nbody.primary_body_name == "Earth"
    @test nbody.body_names == ("Sun",)

    nbody_mars = NBodyGravityModel(["Sun"], "Mars", SPICE_PATH)
    nbody_venus = NBodyGravityModel(["Sun"], "Venus", SPICE_PATH)
    nbody_titan = NBodyGravityModel(["Sun"], "Titan", SPICE_PATH)
    @test lowercase(nbody_mars.planet.name) == "mars"
    @test lowercase(nbody_venus.planet.name) == "venus"
    @test lowercase(nbody_titan.planet.name) == "titan"
    @test_throws ArgumentError NBodyGravityModel(["Sun"], "Pluto", SPICE_PATH)

    nbody_jupiter = NBodyGravityModel(["Jupiter"], "Earth", SPICE_PATH)
    nbody_state = [EARTH.Rp_e + 500e3, 0.0, 0.0, 0.0, 0.0, 0.0, 500.0]
    sc_nbody = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=175.0)
    args_nbody = build_config(
        spacecraft=sc_nbody,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    force_nbody, torque_nbody = calcForceTorque(nbody_jupiter, nbody_state, ODEParams{1}(args=args_nbody), 1)
    @test all(isfinite, force_nbody)
    @test torque_nbody == SVector{3, Float64}(0.0, 0.0, 0.0)

    harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
    harmonics_l20 = GravitationalHarmonicsModel(20, 20, harmonics_file, EARTH)
    @test size(harmonics_l20.C) == (21, 21)
    @test size(harmonics_l20.S) == (21, 21)
    @test harmonics_l20.coefficient_normalization == :full
    @test_throws ArgumentError GravitationalHarmonicsModel(10, 11, harmonics_file, EARTH)

    child_link = Link(root=false, q=MVector{4, Float64}(sin(pi / 4), 0.0, 0.0, cos(pi / 4)))
    rot_child = rotate_to_body(child_link)
    @test size(rot_child) == (3, 3)
    @test isapprox(det(Matrix(rot_child)), 1.0; atol=1e-12)
    @test norm(Matrix(rot_child) - Matrix{Float64}(I, 3, 3)) > 0.1

    @testset "Quaternion DCM Conversion Negative-Trace Branch" begin
        dcm_180_x = SMatrix{3, 3, Float64}(1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0)
        q_neg_trace = SimulationModel.dcm_to_quaternion(dcm_180_x)
        @test isapprox(norm(q_neg_trace), 1.0; atol=1e-12, rtol=0.0)

        dcm_roundtrip = SimulationModel.rot(q_neg_trace)
        @test isapprox(Matrix(dcm_roundtrip), Matrix(dcm_180_x); atol=1e-12, rtol=0.0)
    end

    @testset "Effector Rate Validation" begin
        @test GuidanceModel((), Float64[]) isa GuidanceModel
        @test NavigationModel((), Float64[]) isa NavigationModel
        @test ControlModel((), Float64[]) isa ControlModel

        @test_throws ArgumentError GuidanceModel((:g1,), Float64[])
        @test_throws ArgumentError NavigationModel((:n1,), Float64[])
        @test_throws ArgumentError ControlModel((:c1,), Float64[])

        @test_throws ArgumentError GuidanceModel((:g1,), [0.0])
        @test_throws ArgumentError NavigationModel((:n1,), [-1.0])
        @test_throws ArgumentError ControlModel((:c1,), [Inf])

        sc1 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
        sc2 = make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0)
        args_bad_slots = build_config_multi(
            spacecraft=[sc1, sc2],
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=60.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            control_effectors=(make_base_thruster_model(thrust=1.0, Δv=1.0, start_burn_time=-1.0, stop_burn_time=-1.0),),
            control_rates=[1.0],
            keplerian=true
        )
        @test_throws ArgumentError run_case_silent(args_bad_slots)
    end
end

@testset "Run Simulation Metadata Return" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=20.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        keplerian=true
    )

    metadata = mktempdir() do tmp
        cd(tmp) do
            run_simulation(args; return_solution=true, return_solver_metadata=true)
        end
    end

    @test metadata isa NamedTuple
    @test hasproperty(metadata, :solution)
    @test hasproperty(metadata, :solver_mode)
    @test hasproperty(metadata, :solver_trace)
    @test hasproperty(metadata, :parallel_policy)
    @test hasproperty(metadata, :spice_counters)
end

@testset "Run Simulation Metadata Return (Gravity Backbone)" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=20.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        keplerian=true
    )

    metadata = mktempdir() do tmp
        cd(tmp) do
            withenv(
                "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
                "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "2.0"
            ) do
                run_simulation(args; return_solution=true, return_solver_metadata=true)
            end
        end
    end

    @test metadata isa NamedTuple
    @test metadata.solver_mode == "gravity_backbone_split"
    @test _is_gravity_backbone_state(metadata.solution.u[end])
    @test isapprox(
        _state_mass_kg(metadata.solution.u[end], args, 1),
        args.dynamics_model.spacecraft[1].dry_mass + args.dynamics_model.spacecraft[1].prop_mass;
        atol=0.0,
        rtol=0.0
    )
    @test isempty(_state_heat_loads(metadata.solution.u[end], args, 1)) == false
end

@testset "RHS Completeness: Mass Derivative" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )

    u = build_initial_conditions(args)
    du = copy(u)
    du.sc[1].mass = 789.0
    p = ODEParams{1}(args=args)
    spacecraft_dynamics!(du, u, p, 0.0)
    @test du.sc[1].mass == 0.0

    du_inactive = copy(u)
    du_inactive.sc[1].mass = 123.0
    p_inactive = ODEParams{1}(args=args, is_active=[false])
    spacecraft_dynamics!(du_inactive, u, p_inactive, 0.0)
    @test du_inactive.sc[1].mass == 0.0
end







@testset "Guidance Thruster Campaign Helpers" begin
    ensure_guidance_sandbox_loaded!()
    sandbox = GUIDANCE_SANDBOX

    guidance_model = sandbox.AerobrakingCampaignPropulsiveManeuverGuidanceModel(
        maneuver_orbit_number=[2, 3],
        maneuver_Δv=[-5.0, 20.0]
    )

    args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p = ODEParams{1}(args=args)
    u = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])

    p.orbit_counter[1] = 3
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test p.shared_buffers.maneuver_commands[1].valid == true
    @test p.shared_buffers.maneuver_commands[1].delta_v_mps == 20.0
    @test isapprox(p.shared_buffers.maneuver_commands[1].direction_rad, 0.0; atol=1e-12, rtol=0.0)

    p.orbit_counter[1] = 2
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test p.shared_buffers.maneuver_commands[1].valid == true
    @test p.shared_buffers.maneuver_commands[1].delta_v_mps == 5.0
    @test isapprox(p.shared_buffers.maneuver_commands[1].direction_rad, π; atol=1e-12, rtol=0.0)

    p.orbit_counter[1] = 1
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test p.shared_buffers.maneuver_commands[1].valid == true
    @test p.shared_buffers.maneuver_commands[1].delta_v_mps == 0.0
end

@testset "Odyssey Maneuver Schedule Bridge" begin
    maneuvers = odyssey_campaign_maneuvers(1:20)

    @test !(1 in maneuvers.maneuver_orbit_number)
    @test 7 in maneuvers.maneuver_orbit_number
    @test 14 in maneuvers.maneuver_orbit_number

    idx7 = findfirst(==(7), maneuvers.maneuver_orbit_number)
    idx14 = findfirst(==(14), maneuvers.maneuver_orbit_number)
    @test idx7 !== nothing
    @test idx14 !== nothing
    @test maneuvers.maneuver_Δv[idx7] > 0.0
    @test maneuvers.maneuver_Δv[idx14] < 0.0

    guidance_model = AerobrakingCampaignPropulsiveManeuverGuidanceModel(
        maneuver_orbit_number=maneuvers.maneuver_orbit_number,
        maneuver_Δv=maneuvers.maneuver_Δv
    )
    thruster = BaseThrusterModel(
        thrust=[4.0],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[-1.0],
        stop_burn_time=[-1.0],
        Isp=[300.0]
    )
    @test length(guidance_model.maneuver_orbit_number) == length(guidance_model.maneuver_Δv)
    @test length(thruster.thrust) == 1

    ensure_guidance_sandbox_loaded!()
    sandbox = GUIDANCE_SANDBOX
    sandbox_guidance = sandbox.AerobrakingCampaignPropulsiveManeuverGuidanceModel(
        maneuver_orbit_number=maneuvers.maneuver_orbit_number,
        maneuver_Δv=maneuvers.maneuver_Δv
    )

    args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p = ODEParams{1}(args=args)
    u = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])

    p.orbit_counter[1] = 7
    sandbox.calcGuidanceEffect!(sandbox_guidance, u, p, 0.0, 1)
    @test p.shared_buffers.maneuver_commands[1].delta_v_mps > 0.0
    @test isapprox(p.shared_buffers.maneuver_commands[1].direction_rad, 0.0; atol=1e-12, rtol=0.0)

    p.orbit_counter[1] = 14
    sandbox.calcGuidanceEffect!(sandbox_guidance, u, p, 0.0, 1)
    @test p.shared_buffers.maneuver_commands[1].delta_v_mps > 0.0
    @test isapprox(p.shared_buffers.maneuver_commands[1].direction_rad, π; atol=1e-12, rtol=0.0)
end
