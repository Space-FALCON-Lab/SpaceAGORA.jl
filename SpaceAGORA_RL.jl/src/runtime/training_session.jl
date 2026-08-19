struct TrainingSession{L}
    config::ResolvedConfig
    backend::SpaceAGORAAerobrakingBackend
    learner::L
    rng::MersenneTwister
    output_dir::String
    manifest::RunManifest
end

function build_training_session(config::ResolvedConfig=resolve_config(); run_id::Union{Nothing,String}=nothing)
    rng = MersenneTwister(config.training.seed)
    backend = SpaceAGORAAerobrakingBackend(config.scenario)
    device = resolve_training_device(config.training.device)
    if device isa CUDATrainingDevice
        return Base.invokelatest(_build_training_session, config, run_id, rng, backend, device)
    end
    return _build_training_session(config, run_id, rng, backend, device)
end

function _build_training_session(config::ResolvedConfig, run_id::Union{Nothing,String},
                                 rng::MersenneTwister,
                                 backend::SpaceAGORAAerobrakingBackend,
                                 device::AbstractTrainingDevice)
    learner = if is_ddqn_family_algorithm(config.training.algorithm)
        DDQNLearner(rng, config.ddqn; schedule=config.epsilon, device=device)
    elseif config.training.algorithm == :a2c
        A2CLearner(rng, config.a2c; device=device)
    else
        throw(ArgumentError("unknown training algorithm: $(config.training.algorithm)"))
    end
    manifest = RunManifest(config; run_id=run_id)
    output_dir = joinpath(config.training.output_dir, manifest.run_id)
    mkpath(output_dir)
    write_manifest(
        joinpath(output_dir, "manifest.toml"),
        manifest;
        gram_wind_mode=config.scenario.spaceagora_gram_wind_mode,
    )
    return TrainingSession(config, backend, learner, rng, output_dir, manifest)
end
