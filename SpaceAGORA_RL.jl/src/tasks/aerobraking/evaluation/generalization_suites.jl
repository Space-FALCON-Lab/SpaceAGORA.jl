const PAPER_IID_EVALUATION_EPISODES = 40
const PAPER_GENERALIZATION_EVALUATION_EPISODES = 100
const PAPER_EVALUATION_MODES = ("conservative", "tolerant")
const GENERALIZATION_EVALUATION_REFERENCE_CASE = "iid_reference"
const GENERALIZATION_EVALUATION_CASES = (
    "nominal",
    "exponential_density",
    "aggressive_atmosphere",
    "short_campaign",
    "long_campaign",
    "high_accuracy_spaceagora",
)
const SPACEAGORA_AGGRESSIVE_MGCM_DUST_LEVELS = (0.3, 0.0, 0.0)
const SPACEAGORA_AGGRESSIVE_DUST_STORM = (250.0, 48.0, 3.0, 0.0, 0.0, 0.0)

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
                                        max_passes::Integer=1000)
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
                                                max_passes::Integer=1000)
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
    terminal_on_thermal_violation::Bool=
        config.termination_config.terminal_on_thermal_violation,
)
    term = config.termination_config
    evaluation_termination = TerminationConfig(
        impact_periapsis_altitude_m=term.impact_periapsis_altitude_m,
        out_of_passage_periapsis_altitude_m=term.out_of_passage_periapsis_altitude_m,
        max_passes=Int(max_passes),
        terminal_on_thermal_violation=terminal_on_thermal_violation,
    )
    return default_aerobraking_config(
        phase=config.phase,
        nominal=randomization_config.nominal,
        max_passes=Int(max_passes),
        backend_mode=config.backend_mode,
        training=false,
        spaceagora_atmosphere_model=config.spaceagora_atmosphere_model,
        spaceagora_gram_wind_mode=config.spaceagora_gram_wind_mode,
        spaceagora_gram_once_per_step=config.spaceagora_gram_once_per_step,
        spaceagora_mars_mgcm_dust_levels=config.spaceagora_mars_mgcm_dust_levels,
        spaceagora_mars_dust_storm=config.spaceagora_mars_dust_storm,
        spaceagora_integration_config=config.spaceagora_integration_config,
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

function _generalization_randomization_variant(
    config::AerobrakingRandomizationConfig;
    nominal::Bool=config.nominal,
    marsgram_perturbation_scale::Real=config.marsgram_perturbation_scale,
)
    return AerobrakingRandomizationConfig(
        nominal=nominal,
        apoapsis_jitter_m=nominal ? 0.0 : config.apoapsis_jitter_m,
        periapsis_jitter_m=nominal ? 0.0 : config.periapsis_jitter_m,
        angle_jitter_deg=nominal ? 0.0 : config.angle_jitter_deg,
        nonnominal_inclination_low_deg=config.nonnominal_inclination_low_deg,
        nonnominal_inclination_high_deg=config.nonnominal_inclination_high_deg,
        nonnominal_aop_low_deg=config.nonnominal_aop_low_deg,
        nonnominal_aop_high_deg=config.nonnominal_aop_high_deg,
        nonnominal_raan_low_deg=config.nonnominal_raan_low_deg,
        nonnominal_raan_high_deg=config.nonnominal_raan_high_deg,
        initial_date_start=nominal ? nothing : config.initial_date_start,
        initial_date_days=nominal ? 0 : config.initial_date_days,
        randomize_initial_time_of_day=!nominal && config.randomize_initial_time_of_day,
        initial_true_anomaly_jitter_deg=
            nominal ? 0.0 : config.initial_true_anomaly_jitter_deg,
        process_noise=!nominal && config.process_noise,
        process_noise_scale=nominal ? 0.0 : config.process_noise_scale,
        aerodynamic_coefficient_dispersion=
            !nominal && config.aerodynamic_coefficient_dispersion,
        aerodynamic_coefficient_span=config.aerodynamic_coefficient_span,
        aerodynamic_cd_span=config.aerodynamic_cd_span,
        aerodynamic_cl_span=config.aerodynamic_cl_span,
        marsgram_perturbation_scale=Float64(marsgram_perturbation_scale),
        marsgram_seed_base=config.marsgram_seed_base,
    )
end

function _spaceagora_generalization_variant(
    base::AerobrakingScenarioConfig;
    phase::AbstractString=base.phase,
    atmosphere_model::Symbol=base.spaceagora_atmosphere_model,
    randomization_config::AerobrakingRandomizationConfig=base.randomization_config,
    mars_mgcm_dust_levels=base.spaceagora_mars_mgcm_dust_levels,
    mars_dust_storm=base.spaceagora_mars_dust_storm,
    integration_config::SpaceAGORAIntegrationConfig=base.spaceagora_integration_config,
)
    return default_aerobraking_config(
        phase=phase,
        nominal=randomization_config.nominal,
        max_passes=base.termination_config.max_passes,
        backend_mode=base.backend_mode,
        training=false,
        spaceagora_atmosphere_model=atmosphere_model,
        spaceagora_gram_wind_mode=base.spaceagora_gram_wind_mode,
        spaceagora_gram_once_per_step=base.spaceagora_gram_once_per_step,
        spaceagora_mars_mgcm_dust_levels=mars_mgcm_dust_levels,
        spaceagora_mars_dust_storm=mars_dust_storm,
        spaceagora_integration_config=integration_config,
        spaceagora_tabulated_flight_file=base.spaceagora_tabulated_flight_file,
        spaceagora_tabulated_flight_sigma=base.spaceagora_tabulated_flight_sigma,
        spaceagora_gravity_harmonics_degree=base.spaceagora_gravity_harmonics_degree,
        spaceagora_gravity_harmonics_order=base.spaceagora_gravity_harmonics_order,
        spaceagora_gravity_harmonics_file=base.spaceagora_gravity_harmonics_file,
        reward_config=base.reward_config,
        termination_config=base.termination_config,
        randomization_config=randomization_config,
    )
end

function _high_accuracy_spaceagora_integration(
    base::SpaceAGORAIntegrationConfig,
)
    return SpaceAGORAIntegrationConfig(
        solver_mode=base.solver_mode,
        split_imex_solver=base.split_imex_solver,
        reltol_orbit=min(base.reltol_orbit, 1e-10),
        abstol_orbit=min(base.abstol_orbit, 1e-10),
        dt_max_orbit_s=min(base.dt_max_orbit_s, 5.0),
        reltol_atmosphere=min(base.reltol_atmosphere, 1e-10),
        abstol_atmosphere=min(base.abstol_atmosphere, 1e-12),
        dt_max_atmosphere_s=min(base.dt_max_atmosphere_s, 0.5),
    )
end

"""
Construct the SpaceAGORA-native IID reference and six Table VI-inspired cases.

The supplied configuration is the policy's actual training scenario. The
returned policy-evaluation scenarios change SpaceAGORA model inputs only; they
do not substitute the legacy ABTS simulator or its integrators. Thermal
violations are counted without ending a campaign so Total TV and final mission
outcomes remain observable.
"""
function generalization_evaluation_suite(config::AerobrakingScenarioConfig)
    evaluation = paper_evaluation_scenario(
        config;
        terminal_on_thermal_violation=false,
    )
    _is_spaceagora_live_backend(evaluation.backend_mode) || throw(ArgumentError(
        "generalization_evaluation_suite requires a SpaceAGORA live-physics backend",
    ))
    evaluation.spaceagora_atmosphere_model in (:gram, :marsgram) ||
        throw(ArgumentError(
            "generalization_evaluation_suite requires a native GRAM atmosphere baseline",
        ))

    nominal_randomization = _generalization_randomization_variant(
        evaluation.randomization_config;
        nominal=true,
    )
    nominal = _spaceagora_generalization_variant(
        evaluation;
        randomization_config=nominal_randomization,
        mars_mgcm_dust_levels=nothing,
        mars_dust_storm=nothing,
    )
    aggressive_randomization = _generalization_randomization_variant(
        nominal_randomization;
        nominal=true,
        marsgram_perturbation_scale=
            2.0 * evaluation.randomization_config.marsgram_perturbation_scale,
    )

    return Pair{String,AerobrakingScenarioConfig}[
        GENERALIZATION_EVALUATION_REFERENCE_CASE => evaluation,
        "nominal" => nominal,
        "exponential_density" => _spaceagora_generalization_variant(
            nominal;
            atmosphere_model=:exponential,
            mars_mgcm_dust_levels=nothing,
            mars_dust_storm=nothing,
        ),
        "aggressive_atmosphere" => _spaceagora_generalization_variant(
            nominal;
            randomization_config=aggressive_randomization,
            mars_mgcm_dust_levels=SPACEAGORA_AGGRESSIVE_MGCM_DUST_LEVELS,
            mars_dust_storm=SPACEAGORA_AGGRESSIVE_DUST_STORM,
        ),
        "short_campaign" => _spaceagora_generalization_variant(
            nominal;
            phase="Walkout",
        ),
        "long_campaign" => _spaceagora_generalization_variant(
            nominal;
            phase="Campaign",
        ),
        "high_accuracy_spaceagora" => _spaceagora_generalization_variant(
            nominal;
            integration_config=_high_accuracy_spaceagora_integration(
                nominal.spaceagora_integration_config,
            ),
        ),
    ]
end

function paper_evaluation_mode_scenarios(config::AerobrakingScenarioConfig)
    return Dict(
        "conservative" => paper_evaluation_scenario(
            config;
            terminal_on_thermal_violation=true,
        ),
        "tolerant" => paper_evaluation_scenario(
            config;
            terminal_on_thermal_violation=false,
        ),
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
