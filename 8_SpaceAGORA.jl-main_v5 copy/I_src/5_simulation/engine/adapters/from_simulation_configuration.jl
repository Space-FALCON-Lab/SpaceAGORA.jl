function SimulationEngineConfig(args::SimulationConfiguration; env=ENV)
    config = simulation_engine_config_from_env(env)

    # Keep this adapter intentionally conservative during compatibility window.
    # We bind typed artifact policy defaults from SimulationConfiguration and
    # let explicit env overrides still win via `env_overrides`.
    artifacts = ArtifactConfig(
        save_bundle=config.artifacts.save_bundle,
        warn_deprecated_config=config.artifacts.warn_deprecated_config
    )

    return SimulationEngineConfig(
        parallel=config.parallel,
        solver=config.solver,
        runtime_policy=config.runtime_policy,
        artifacts=artifacts,
        env_overrides=config.env_overrides
    )
end
