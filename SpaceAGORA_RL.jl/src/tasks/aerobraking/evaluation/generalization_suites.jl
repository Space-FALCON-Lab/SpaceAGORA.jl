function generalization_suite_configs(config::AerobrakingScenarioConfig)
    return Dict(
        "nominal" => config,
        "iid_randomized" => default_aerobraking_config(phase=config.phase, nominal=false,
                                                       max_passes=config.termination_config.max_passes,
                                                       training=config.training),
    )
end
