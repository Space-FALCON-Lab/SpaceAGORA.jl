module PublicAPISymbols

using Base.Docs

export PUBLIC_API_SECTIONS, public_api_specs, render_public_api_markdown, undocumented_public_api_specs

const PUBLIC_API_SECTIONS = [
    (
        title = "Simulation",
        items = [
            (owner = :SpaceAGORA, symbol = :run_simulation, rendered = "SpaceAGORA.run_simulation")
        ]
    ),
    (
        title = "Runtime Configuration",
        items = [
            (owner = :SimulationEngine, symbol = :ParallelConfig, rendered = "SpaceAGORA.SimulationEngine.ParallelConfig"),
            (owner = :SimulationEngine, symbol = :SolverConfig, rendered = "SpaceAGORA.SimulationEngine.SolverConfig"),
            (owner = :SimulationEngine, symbol = :RuntimePolicyConfig, rendered = "SpaceAGORA.SimulationEngine.RuntimePolicyConfig"),
            (owner = :SimulationEngine, symbol = :ArtifactConfig, rendered = "SpaceAGORA.SimulationEngine.ArtifactConfig"),
            (owner = :SimulationEngine, symbol = :SimulationEngineConfig, rendered = "SpaceAGORA.SimulationEngine.SimulationEngineConfig"),
            (owner = :SimulationEngine, symbol = :simulation_engine_config_from_env, rendered = "SpaceAGORA.SimulationEngine.simulation_engine_config_from_env")
        ]
    ),
    (
        title = "Parallel Profiles and Routing",
        items = [
            (owner = :ParallelProfiles, symbol = :ParallelProfile, rendered = "SpaceAGORA.ParallelProfiles.ParallelProfile"),
            (owner = :ParallelProfiles, symbol = :ParallelProfileConfig, rendered = "SpaceAGORA.ParallelProfiles.ParallelProfileConfig"),
            (owner = :ParallelProfiles, symbol = :parse_parallel_profile, rendered = "SpaceAGORA.ParallelProfiles.parse_parallel_profile"),
            (owner = :ParallelProfiles, symbol = :parallel_profile_name, rendered = "SpaceAGORA.ParallelProfiles.parallel_profile_name"),
            (owner = :ParallelProfiles, symbol = :profile_config, rendered = "SpaceAGORA.ParallelProfiles.profile_config"),
            (owner = :ParallelProfiles, symbol = :profile_env_pairs, rendered = "SpaceAGORA.ParallelProfiles.profile_env_pairs"),
            (owner = :ParallelProfiles, symbol = :with_parallel_profile, rendered = "SpaceAGORA.ParallelProfiles.with_parallel_profile"),
            (owner = :ParallelProfiles, symbol = :OuterRouteFeatures, rendered = "SpaceAGORA.ParallelProfiles.OuterRouteFeatures"),
            (owner = :ParallelProfiles, symbol = :OuterRouteTuning, rendered = "SpaceAGORA.ParallelProfiles.OuterRouteTuning"),
            (owner = :ParallelProfiles, symbol = :OuterRouteState, rendered = "SpaceAGORA.ParallelProfiles.OuterRouteState"),
            (owner = :ParallelProfiles, symbol = :reset_outer_route_state!, rendered = "SpaceAGORA.ParallelProfiles.reset_outer_route_state!"),
            (owner = :ParallelProfiles, symbol = :outer_route_signature, rendered = "SpaceAGORA.ParallelProfiles.outer_route_signature"),
            (owner = :ParallelProfiles, symbol = :outer_route_stats_snapshot, rendered = "SpaceAGORA.ParallelProfiles.outer_route_stats_snapshot"),
            (owner = :ParallelProfiles, symbol = :default_outer_route, rendered = "SpaceAGORA.ParallelProfiles.default_outer_route"),
            (owner = :ParallelProfiles, symbol = :outer_route_candidates, rendered = "SpaceAGORA.ParallelProfiles.outer_route_candidates"),
            (owner = :ParallelProfiles, symbol = :select_outer_route!, rendered = "SpaceAGORA.ParallelProfiles.select_outer_route!"),
            (owner = :ParallelProfiles, symbol = :record_outer_route_feedback!, rendered = "SpaceAGORA.ParallelProfiles.record_outer_route_feedback!")
        ]
    ),
    (
        title = "Verification",
        items = [
            (owner = :TelemetryVerification, symbol = :VerificationRequest, rendered = "SpaceAGORA.TelemetryVerification.VerificationRequest"),
            (owner = :TelemetryVerification, symbol = :VerificationResult, rendered = "SpaceAGORA.TelemetryVerification.VerificationResult"),
            (owner = :TelemetryVerification, symbol = :run_verification, rendered = "SpaceAGORA.TelemetryVerification.run_verification"),
            (owner = :TelemetryVerification, symbol = :run_verification_cli, rendered = "SpaceAGORA.TelemetryVerification.run_verification_cli"),
            (owner = :TelemetryVerification, symbol = :run_study, rendered = "SpaceAGORA.TelemetryVerification.run_study")
        ]
    )
]

@inline function _owner_module(spaceagora::Module, owner::Symbol)::Module
    if owner === :SpaceAGORA
        return spaceagora
    elseif owner === :SimulationEngine
        return getproperty(spaceagora, :SimulationEngine)
    elseif owner === :ParallelProfiles
        return getproperty(spaceagora, :ParallelProfiles)
    elseif owner === :TelemetryVerification
        return getproperty(spaceagora, :TelemetryVerification)
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
    println(io, "This page documents the stable exported interface. Symbols below are available from `SpaceAGORA`; the rendered docstrings come from their canonical owner modules.")
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
