module PublicAPISymbols

using Base.Docs

export PUBLIC_API_SECTIONS, public_api_specs, render_public_api_markdown, undocumented_public_api_specs

const PUBLIC_API_SECTIONS = [
    (
        title = "Simulation",
        items = [
            (owner = :SpaceAGORA, symbol = :run_simulation, rendered = "SpaceAGORA.run_simulation"),
            (owner = :SpaceAGORA, symbol = :prewarm_nbody_ephemeris_cache, rendered = "SpaceAGORA.prewarm_nbody_ephemeris_cache"),
            (owner = :SpaceAGORA, symbol = :load_nbody_ephemeris_cache!, rendered = "SpaceAGORA.load_nbody_ephemeris_cache!"),
            (owner = :SimulationCampaigns, symbol = :MonteCarloSpec, rendered = "SpaceAGORA.MonteCarloSpec"),
            (owner = :SimulationCampaigns, symbol = :MonteCarloSampleResult, rendered = "SpaceAGORA.MonteCarloSampleResult"),
            (owner = :SimulationCampaigns, symbol = :MonteCarloResult, rendered = "SpaceAGORA.MonteCarloResult"),
            (owner = :SimulationCampaigns, symbol = :run_monte_carlo, rendered = "SpaceAGORA.run_monte_carlo")
        ]
    ),
    (
        title = "Runtime Configuration",
        items = [
            (owner = :SimulationEngine, symbol = :ParallelConfig, rendered = "SpaceAGORA.ParallelConfig"),
            (owner = :SimulationEngine, symbol = :SolverConfig, rendered = "SpaceAGORA.SolverConfig"),
            (owner = :SimulationEngine, symbol = :RuntimePolicyConfig, rendered = "SpaceAGORA.RuntimePolicyConfig"),
            (owner = :SimulationEngine, symbol = :ArtifactConfig, rendered = "SpaceAGORA.ArtifactConfig"),
            (owner = :SimulationEngine, symbol = :SimulationEngineConfig, rendered = "SpaceAGORA.SimulationEngineConfig"),
            (owner = :SimulationEngine, symbol = :simulation_engine_config_from_env, rendered = "SpaceAGORA.simulation_engine_config_from_env")
        ]
    ),
    (
        title = "No-GRAM and Baseline Models",
        items = [
            (owner = :SpaceAGORA, symbol = :NoAtmosphereModel, rendered = "SpaceAGORA.NoAtmosphereModel"),
            (owner = :SpaceAGORA, symbol = :ExponentialAtmosphereModel, rendered = "SpaceAGORA.ExponentialAtmosphereModel"),
            (owner = :SpaceAGORA, symbol = :PiecewiseExponentialAtmosphereModel, rendered = "SpaceAGORA.PiecewiseExponentialAtmosphereModel"),
            (owner = :SpaceAGORA, symbol = :NRLMSISE00AtmosphereModel, rendered = "SpaceAGORA.NRLMSISE00AtmosphereModel"),
            (owner = :SpaceAGORA, symbol = :init_nrlmsise_space_indices!, rendered = "SpaceAGORA.init_nrlmsise_space_indices!"),
            (owner = :SpaceAGORA, symbol = :SimpleEphemeridesModel, rendered = "SpaceAGORA.SimpleEphemeridesModel"),
            (owner = :SpaceAGORA, symbol = :make_no_gram_planet, rendered = "SpaceAGORA.make_no_gram_planet"),
            (owner = :SpaceAGORA, symbol = :make_no_gram_density_model, rendered = "SpaceAGORA.make_no_gram_density_model"),
            (owner = :SpaceAGORA, symbol = :make_no_gram_environment, rendered = "SpaceAGORA.make_no_gram_environment")
        ]
    ),
    (
        title = "Model Extension Interfaces",
        items = [
            (owner = :SpaceAGORA, symbol = :AbstractForceTorqueModel, rendered = "SpaceAGORA.AbstractForceTorqueModel"),
            (owner = :SpaceAGORA, symbol = :AbstractPlanet, rendered = "SpaceAGORA.AbstractPlanet"),
            (owner = :SpaceAGORA, symbol = :AbstractDensityModel, rendered = "SpaceAGORA.AbstractDensityModel"),
            (owner = :SpaceAGORA, symbol = :AbstractControlEffectorModel, rendered = "SpaceAGORA.AbstractControlEffectorModel"),
            (owner = :SpaceAGORA, symbol = :AbstractEphemeridesModel, rendered = "SpaceAGORA.AbstractEphemeridesModel"),
            (owner = :SpaceAGORA, symbol = :AbstractThermalModel, rendered = "SpaceAGORA.AbstractThermalModel"),
            (owner = :SpaceAGORA, symbol = :AbstractThrusterModel, rendered = "SpaceAGORA.AbstractThrusterModel"),
            (owner = :SpaceAGORA, symbol = :AbstractGuidanceModel, rendered = "SpaceAGORA.AbstractGuidanceModel"),
            (owner = :SpaceAGORA, symbol = :StateSample, rendered = "SpaceAGORA.StateSample"),
            (owner = :SpaceAGORA, symbol = :PlanetFrameSample, rendered = "SpaceAGORA.PlanetFrameSample"),
            (owner = :SpaceAGORA, symbol = :AtmosphereSample, rendered = "SpaceAGORA.AtmosphereSample"),
            (owner = :SpaceAGORA, symbol = :SolarEphemerisSample, rendered = "SpaceAGORA.SolarEphemerisSample"),
            (owner = :SpaceAGORA, symbol = :ThirdBodyEphemerisSample, rendered = "SpaceAGORA.ThirdBodyEphemerisSample"),
            (owner = :SpaceAGORA, symbol = :EnvironmentSample, rendered = "SpaceAGORA.EnvironmentSample"),
            (owner = :SpaceAGORA, symbol = :EffectorEnvironmentRequirements, rendered = "SpaceAGORA.EffectorEnvironmentRequirements")
        ]
    ),
    (
        title = "Model Extension Hooks",
        items = [
            (owner = :SpaceAGORA, symbol = :calcForceTorque, rendered = "SpaceAGORA.calcForceTorque"),
            (owner = :SpaceAGORA, symbol = :wrench, rendered = "SpaceAGORA.wrench"),
            (owner = :SpaceAGORA, symbol = :environment_requirements, rendered = "SpaceAGORA.environment_requirements"),
            (owner = :SpaceAGORA, symbol = :solver_partition, rendered = "SpaceAGORA.solver_partition"),
            (owner = :SpaceAGORA, symbol = :gravity_backbone_structure, rendered = "SpaceAGORA.gravity_backbone_structure"),
            (owner = :SpaceAGORA, symbol = :gravity_backbone_acceleration_ii, rendered = "SpaceAGORA.gravity_backbone_acceleration_ii"),
            (owner = :SpaceAGORA, symbol = :gravity_backbone_kick_structure, rendered = "SpaceAGORA.gravity_backbone_kick_structure"),
            (owner = :SpaceAGORA, symbol = :gravity_backbone_kick_acceleration_ii, rendered = "SpaceAGORA.gravity_backbone_kick_acceleration_ii"),
            (owner = :SpaceAGORA, symbol = :getDensity, rendered = "SpaceAGORA.getDensity"),
            (owner = :SpaceAGORA, symbol = :getDensityBatch!, rendered = "SpaceAGORA.getDensityBatch!"),
            (owner = :SpaceAGORA, symbol = :calcControlEffect!, rendered = "SpaceAGORA.calcControlEffect!"),
            (owner = :SpaceAGORA, symbol = :calcControlForceTorque, rendered = "SpaceAGORA.calcControlForceTorque"),
            (owner = :SpaceAGORA, symbol = :calcControlMassFlowRate, rendered = "SpaceAGORA.calcControlMassFlowRate")
        ]
    ),
    (
        title = "Robotics and Planning",
        items = [
            (owner = :SpaceAGORA, symbol = :ClothArmBasePose, rendered = "SpaceAGORA.ClothArmBasePose"),
            (owner = :SpaceAGORA, symbol = :ClothArmModel, rendered = "SpaceAGORA.ClothArmModel"),
            (owner = :SpaceAGORA, symbol = :RobotArmPlan, rendered = "SpaceAGORA.RobotArmPlan")
        ]
    ),
    (
        title = "Parallel Profiles and Routing",
        items = [
            (owner = :ParallelProfiles, symbol = :ParallelProfile, rendered = "SpaceAGORA.ParallelProfile"),
            (owner = :ParallelProfiles, symbol = :ParallelProfileConfig, rendered = "SpaceAGORA.ParallelProfileConfig"),
            (owner = :ParallelProfiles, symbol = :parse_parallel_profile, rendered = "SpaceAGORA.parse_parallel_profile"),
            (owner = :ParallelProfiles, symbol = :parallel_profile_name, rendered = "SpaceAGORA.parallel_profile_name"),
            (owner = :ParallelProfiles, symbol = :profile_config, rendered = "SpaceAGORA.profile_config"),
            (owner = :ParallelProfiles, symbol = :profile_env_pairs, rendered = "SpaceAGORA.profile_env_pairs"),
            (owner = :ParallelProfiles, symbol = :with_parallel_profile, rendered = "SpaceAGORA.with_parallel_profile"),
            (owner = :ParallelProfiles, symbol = :OuterRouteFeatures, rendered = "SpaceAGORA.OuterRouteFeatures"),
            (owner = :ParallelProfiles, symbol = :OuterRouteTuning, rendered = "SpaceAGORA.OuterRouteTuning"),
            (owner = :ParallelProfiles, symbol = :OuterRouteState, rendered = "SpaceAGORA.OuterRouteState"),
            (owner = :ParallelProfiles, symbol = :reset_outer_route_state!, rendered = "SpaceAGORA.reset_outer_route_state!"),
            (owner = :ParallelProfiles, symbol = :outer_route_signature, rendered = "SpaceAGORA.outer_route_signature"),
            (owner = :ParallelProfiles, symbol = :outer_route_stats_snapshot, rendered = "SpaceAGORA.outer_route_stats_snapshot"),
            (owner = :ParallelProfiles, symbol = :default_outer_route, rendered = "SpaceAGORA.default_outer_route"),
            (owner = :ParallelProfiles, symbol = :outer_route_candidates, rendered = "SpaceAGORA.outer_route_candidates"),
            (owner = :ParallelProfiles, symbol = :select_outer_route!, rendered = "SpaceAGORA.select_outer_route!"),
            (owner = :ParallelProfiles, symbol = :record_outer_route_feedback!, rendered = "SpaceAGORA.record_outer_route_feedback!")
        ]
    ),
    (
        title = "Verification",
        items = [
            (owner = :TelemetryVerification, symbol = :VerificationRequest, rendered = "SpaceAGORA.VerificationRequest"),
            (owner = :TelemetryVerification, symbol = :VerificationResult, rendered = "SpaceAGORA.VerificationResult"),
            (owner = :TelemetryVerification, symbol = :run_verification, rendered = "SpaceAGORA.run_verification"),
            (owner = :TelemetryVerification, symbol = :run_verification_cli, rendered = "SpaceAGORA.run_verification_cli"),
            (owner = :TelemetryVerification, symbol = :run_study, rendered = "SpaceAGORA.run_study")
        ]
    ),
    (
        title = "CLI and Assets",
        items = [
            (owner = :SpaceAGORACLI, symbol = :AssetCheckItem, rendered = "SpaceAGORA.AssetCheckItem"),
            (owner = :SpaceAGORACLI, symbol = :AssetCheckReport, rendered = "SpaceAGORA.AssetCheckReport"),
            (owner = :SpaceAGORA, symbol = :check_assets, rendered = "SpaceAGORA.check_assets"),
            (owner = :SpaceAGORA, symbol = :render_asset_report, rendered = "SpaceAGORA.render_asset_report"),
            (owner = :SpaceAGORA, symbol = :run_cli, rendered = "SpaceAGORA.run_cli")
        ]
    )
]

@inline function _owner_module(spaceagora::Module, owner::Symbol)::Module
    if owner === :SpaceAGORA
        return spaceagora
    elseif owner === :SimulationEngine
        return getproperty(spaceagora, :SimulationEngine)
    elseif owner === :SimulationCampaigns
        return getproperty(spaceagora, :SimulationCampaigns)
    elseif owner === :ParallelProfiles
        return getproperty(spaceagora, :ParallelProfiles)
    elseif owner === :TelemetryVerification
        return getproperty(spaceagora, :TelemetryVerification)
    elseif owner === :SpaceAGORACLI
        return getproperty(spaceagora, :SpaceAGORACLI)
    end
    error("Unsupported public API owner: $(owner)")
end

function public_api_specs(spaceagora::Module)
    specs = NamedTuple[]
    for section in PUBLIC_API_SECTIONS
        for item in section.items
            push!(specs, (
                section = section.title,
                owner = item.owner,
                owner_module = _owner_module(spaceagora, item.owner),
                symbol = item.symbol,
                rendered = item.rendered,
            ))
        end
    end
    return specs
end

function render_public_api_markdown(spaceagora::Module)::String
    io = IOBuffer()
    println(io, "# Public API")
    println(io)
    println(io, "This page documents the stable exported interface available from `SpaceAGORA`.")
    println(io)
    for section in PUBLIC_API_SECTIONS
        println(io, "## $(section.title)")
        println(io)
        println(io, "```@docs")
        for item in section.items
            println(io, item.rendered)
        end
        println(io, "```")
        println(io)
    end
    return String(take!(io))
end

function undocumented_public_api_specs(spaceagora::Module)
    missing = NamedTuple[]
    for spec in public_api_specs(spaceagora)
        binding = Docs.Binding(spec.owner_module, spec.symbol)
        doc = Docs.doc(binding)
        if doc === nothing
            push!(missing, spec)
        end
    end
    return missing
end

end # module PublicAPISymbols
