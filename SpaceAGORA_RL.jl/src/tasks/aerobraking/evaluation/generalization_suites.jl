const PAPER_IID_EVALUATION_EPISODES = 40
const PAPER_GENERALIZATION_EVALUATION_EPISODES = 100

function _paper_pr_drl_reward_config()
    return RewardConfig(
        heat_low_w_cm2=0.05,
        heat_high_w_cm2=0.25,
        heat_medium_w_cm2=0.30,
        heat_hard_w_cm2=0.45,
    )
end

function _paper_pr_drl_termination_config(max_passes::Integer)
    return TerminationConfig(
        impact_periapsis_altitude_m=85e3,
        out_of_passage_periapsis_altitude_m=135e3,
        max_passes=Int(max_passes),
        terminal_on_thermal_violation=true,
    )
end

function _paper_pr_drl_randomization_config(; nominal::Bool=false,
                                            process_noise_scale::Real=0.4,
                                            aerodynamic_coefficient_dispersion::Bool=true)
    return AerobrakingRandomizationConfig(
        nominal=nominal,
        apoapsis_jitter_m=nominal ? 0.0 : 2500.0,
        periapsis_jitter_m=nominal ? 0.0 : 2500.0,
        angle_jitter_deg=nominal ? 0.0 : 0.25,
        nonnominal_inclination_low_deg=88.6,
        nonnominal_inclination_high_deg=98.6,
        nonnominal_aop_low_deg=60.0,
        nonnominal_aop_high_deg=90.0,
        nonnominal_raan_low_deg=110.0,
        nonnominal_raan_high_deg=120.0,
        initial_date_start=nominal ? nothing : Date(2001, 12, 1),
        initial_date_days=nominal ? 0 : 31,
        randomize_initial_time_of_day=!nominal,
        initial_true_anomaly_jitter_deg=nominal ? 0.0 : 0.025,
        process_noise=Float64(process_noise_scale) != 0.0,
        process_noise_scale=Float64(process_noise_scale),
        aerodynamic_coefficient_dispersion=aerodynamic_coefficient_dispersion,
        aerodynamic_coefficient_span=0.10,
        aerodynamic_cd_span=0.10,
        aerodynamic_cl_span=0.10,
        marsgram_perturbation_scale=1.0,
        marsgram_seed_base=1001,
    )
end

function paper_pr_drl_evaluation_config(; process_noise_scale::Real=0.4,
                                        training::Bool=false,
                                        backend_mode::Symbol=:paper_surrogate,
                                        phase::AbstractString="Main",
                                        max_passes::Integer=250)
    return default_aerobraking_config(
        phase=phase,
        nominal=false,
        max_passes=Int(max_passes),
        backend_mode=backend_mode,
        training=training,
        reward_config=_paper_pr_drl_reward_config(),
        termination_config=_paper_pr_drl_termination_config(max_passes),
        randomization_config=_paper_pr_drl_randomization_config(
            nominal=false,
            process_noise_scale=process_noise_scale,
            aerodynamic_coefficient_dispersion=true,
        ),
    )
end

paper_pr_drl_marsgram_evaluation_config(; kwargs...) =
    paper_pr_drl_evaluation_config(; backend_mode=:spaceagora_marsgram, kwargs...)

paper_pr_drl_physics_evaluation_config(; kwargs...) =
    paper_pr_drl_evaluation_config(; backend_mode=:spaceagora_physics, kwargs...)

function paper_odyssey_flight_evaluation_config(; process_noise_scale::Real=0.0,
                                                training::Bool=false,
                                                backend_mode::Symbol=:paper_surrogate,
                                                phase::AbstractString="Main",
                                                max_passes::Integer=250)
    return default_aerobraking_config(
        phase=phase,
        nominal=true,
        max_passes=Int(max_passes),
        backend_mode=backend_mode,
        training=training,
        reward_config=_paper_pr_drl_reward_config(),
        termination_config=_paper_pr_drl_termination_config(max_passes),
        randomization_config=_paper_pr_drl_randomization_config(
            nominal=true,
            process_noise_scale=process_noise_scale,
            aerodynamic_coefficient_dispersion=false,
        ),
    )
end

function paper_evaluation_scenario(
    config::AerobrakingScenarioConfig;
    max_passes::Integer=config.termination_config.max_passes,
    randomization_config::AerobrakingRandomizationConfig=config.randomization_config,
)
    term = config.termination_config
    evaluation_termination = TerminationConfig(
        impact_periapsis_altitude_m=term.impact_periapsis_altitude_m,
        out_of_passage_periapsis_altitude_m=term.out_of_passage_periapsis_altitude_m,
        max_passes=Int(max_passes),
        terminal_on_thermal_violation=term.terminal_on_thermal_violation,
    )
    return default_aerobraking_config(
        phase=config.phase,
        nominal=randomization_config.nominal,
        max_passes=Int(max_passes),
        backend_mode=config.backend_mode,
        training=false,
        spaceagora_atmosphere_model=config.spaceagora_atmosphere_model,
        spaceagora_tabulated_flight_file=config.spaceagora_tabulated_flight_file,
        spaceagora_tabulated_flight_sigma=config.spaceagora_tabulated_flight_sigma,
        spaceagora_gravity_harmonics_degree=config.spaceagora_gravity_harmonics_degree,
        spaceagora_gravity_harmonics_order=config.spaceagora_gravity_harmonics_order,
        spaceagora_gravity_harmonics_file=config.spaceagora_gravity_harmonics_file,
        reward_config=config.reward_config,
        termination_config=evaluation_termination,
        randomization_config=randomization_config,
    )
end

function generalization_suite_configs(config::AerobrakingScenarioConfig)
    evaluation_config = paper_evaluation_scenario(config)
    iid_randomization = config.randomization_config.nominal ?
                        _paper_pr_drl_randomization_config(
                            nominal=false,
                            process_noise_scale=0.4,
                            aerodynamic_coefficient_dispersion=true,
                        ) :
                        config.randomization_config
    nominal_randomization = _paper_pr_drl_randomization_config(
        nominal=true,
        process_noise_scale=0.0,
        aerodynamic_coefficient_dispersion=false,
    )
    return Dict(
        "nominal" => paper_evaluation_scenario(
            evaluation_config;
            randomization_config=nominal_randomization,
        ),
        "iid_randomized" => paper_evaluation_scenario(
            evaluation_config;
            randomization_config=iid_randomization,
        ),
    )
end
