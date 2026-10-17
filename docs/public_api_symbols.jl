module PublicAPISymbols

using Base.Docs

export PUBLIC_API_SECTIONS, public_api_specs, render_public_api_markdown, undocumented_public_api_specs,
    KNOWN_PLACEHOLDER_DOCS, stale_placeholder_entries

const PUBLIC_API_SECTIONS = [
    (title = "Constellation laser interlinks", items = [
        (owner = :SpaceAGORA, symbol = :ConstellationInterLinkModel, rendered = "SpaceAGORA.ConstellationInterLinkModel"),
        (owner = :SpaceAGORA, symbol = :InterLinkGraph, rendered = "SpaceAGORA.InterLinkGraph"),
        (owner = :SpaceAGORA, symbol = :InterLinkParameters, rendered = "SpaceAGORA.InterLinkParameters"),
        (owner = :SpaceAGORA, symbol = :InterLinkKey, rendered = "SpaceAGORA.InterLinkKey"),
        (owner = :SpaceAGORA, symbol = :InterLinkHistory, rendered = "SpaceAGORA.InterLinkHistory"),
        (owner = :SpaceAGORA, symbol = :SchedulingPolicyModel, rendered = "SpaceAGORA.SchedulingPolicyModel"),
        (owner = :SpaceAGORA, symbol = :register_allowed_link!, rendered = "SpaceAGORA.register_allowed_link!"),
        (owner = :SpaceAGORA, symbol = :update_availability!, rendered = "SpaceAGORA.update_availability!"),
        (owner = :SpaceAGORA, symbol = :score_candidates!, rendered = "SpaceAGORA.score_candidates!"),
        (owner = :SpaceAGORA, symbol = :select_interlinks!, rendered = "SpaceAGORA.select_interlinks!"),
        (owner = :SpaceAGORA, symbol = :schedule_interlinks!, rendered = "SpaceAGORA.schedule_interlinks!"),
        (owner = :SpaceAGORA, symbol = :interlink_scheduler_callback, rendered = "SpaceAGORA.interlink_scheduler_callback"),
        (owner = :SpaceAGORA, symbol = :force_on_endpoint, rendered = "SpaceAGORA.force_on_endpoint"),
        (owner = :SpaceAGORA, symbol = :laser_force_on_spacecraft, rendered = "SpaceAGORA.laser_force_on_spacecraft"),
    ]),
    (title = "RPO planner pilot", items = [
        (owner = :SpaceAGORA, symbol = :hypr_available, rendered = "SpaceAGORA.hypr_available"),
        (owner = :SpaceAGORA, symbol = :HYPRUnavailableError, rendered = "SpaceAGORA.HYPRUnavailableError"),
        (owner = :SpaceAGORA, symbol = :AbstractRPOPlanner, rendered = "SpaceAGORA.AbstractRPOPlanner"),
        (owner = :SpaceAGORA, symbol = :RPOPlanningConstraints, rendered = "SpaceAGORA.RPOPlanningConstraints"),
        (owner = :SpaceAGORA, symbol = :RPOValidationSettings, rendered = "SpaceAGORA.RPOValidationSettings"),
        (owner = :SpaceAGORA, symbol = :RPOPlanningRequest, rendered = "SpaceAGORA.RPOPlanningRequest"),
        (owner = :SpaceAGORA, symbol = :RPOReference, rendered = "SpaceAGORA.RPOReference"),
        (owner = :SpaceAGORA, symbol = :RPOPlanningResult, rendered = "SpaceAGORA.RPOPlanningResult"),
        (owner = :SpaceAGORA, symbol = :RPOPlannerCapabilities, rendered = "SpaceAGORA.RPOPlannerCapabilities"),
        (owner = :SpaceAGORA, symbol = :RPOValidationResult, rendered = "SpaceAGORA.RPOValidationResult"),
        (owner = :SpaceAGORA, symbol = :RPOPlanningHeadroom, rendered = "SpaceAGORA.RPOPlanningHeadroom"),
        (owner = :SpaceAGORA, symbol = :planner_capabilities, rendered = "SpaceAGORA.planner_capabilities"),
        (owner = :SpaceAGORA, symbol = :initialize_planner, rendered = "SpaceAGORA.initialize_planner"),
        (owner = :SpaceAGORA, symbol = :plan_rpo!, rendered = "SpaceAGORA.plan_rpo!"),
        (owner = :SpaceAGORA, symbol = :retime_rpo!, rendered = "SpaceAGORA.retime_rpo!"),
        (owner = :SpaceAGORA, symbol = :validate_rpo_result, rendered = "SpaceAGORA.validate_rpo_result"),
        (owner = :SpaceAGORA, symbol = :validate_rpo_capabilities, rendered = "SpaceAGORA.validate_rpo_capabilities"),
        (owner = :SpaceAGORA, symbol = :rpo_reference_is_current, rendered = "SpaceAGORA.rpo_reference_is_current"),
        (owner = :SpaceAGORA, symbol = :DirectRPOPlanner, rendered = "SpaceAGORA.DirectRPOPlanner"),
        (owner = :SpaceAGORA, symbol = :HYPRRPOPlanner, rendered = "SpaceAGORA.HYPRRPOPlanner"),
        (owner = :SpaceAGORA, symbol = :RPOPSOConfig, rendered = "SpaceAGORA.RPOPSOConfig"),
        (owner = :SpaceAGORA, symbol = :rpo_pso_config, rendered = "SpaceAGORA.rpo_pso_config"),
        (owner = :SpaceAGORA, symbol = :RPOPlanningEvent, rendered = "SpaceAGORA.RPOPlanningEvent"),
        (owner = :SpaceAGORA, symbol = :RPOPlanningError, rendered = "SpaceAGORA.RPOPlanningError"),
        (owner = :SpaceAGORA, symbol = :make_rpo_configuration, rendered = "SpaceAGORA.make_rpo_configuration"),
        (owner = :SpaceAGORA, symbol = :rpo_run_report, rendered = "SpaceAGORA.rpo_run_report"),
        (owner = :SpaceAGORA, symbol = :AbstractRPOGeometryProvider, rendered = "SpaceAGORA.AbstractRPOGeometryProvider"),
        (owner = :SpaceAGORA, symbol = :rpo_geometry_revision, rendered = "SpaceAGORA.rpo_geometry_revision"),
        (owner = :SpaceAGORA, symbol = :rpo_geometry_snapshot, rendered = "SpaceAGORA.rpo_geometry_snapshot"),
        (owner = :SpaceAGORA, symbol = :RPOReplanningSpheres, rendered = "SpaceAGORA.RPOReplanningSpheres"),
    ]),
    (title = "GNC algorithm contract", items = [
        (owner = :SpaceAGORA, symbol = :AbstractGNCAlgorithm, rendered = "SpaceAGORA.AbstractGNCAlgorithm"),
        (owner = :SpaceAGORA, symbol = :AbstractAlgorithmFamily, rendered = "SpaceAGORA.AbstractAlgorithmFamily"),
        (owner = :SpaceAGORA, symbol = :AbstractAlgorithmRequest, rendered = "SpaceAGORA.AbstractAlgorithmRequest"),
        (owner = :SpaceAGORA, symbol = :AbstractAlgorithmResult, rendered = "SpaceAGORA.AbstractAlgorithmResult"),
        (owner = :SpaceAGORA, symbol = :AbstractValidationResult, rendered = "SpaceAGORA.AbstractValidationResult"),
        (owner = :SpaceAGORA, symbol = :GNCAlgorithmError, rendered = "SpaceAGORA.GNCAlgorithmError"),
        (owner = :SpaceAGORA, symbol = :AlgorithmResult, rendered = "SpaceAGORA.AlgorithmResult"),
        (owner = :SpaceAGORA, symbol = :ValidationResult, rendered = "SpaceAGORA.ValidationResult"),
        (owner = :SpaceAGORA, symbol = :RPOReferenceFamily, rendered = "SpaceAGORA.RPOReferenceFamily"),
        (owner = :SpaceAGORA, symbol = :algorithm_family, rendered = "SpaceAGORA.algorithm_family"),
        (owner = :SpaceAGORA, symbol = :default_capabilities, rendered = "SpaceAGORA.default_capabilities"),
        (owner = :SpaceAGORA, symbol = :algorithm_capabilities, rendered = "SpaceAGORA.algorithm_capabilities"),
        (owner = :SpaceAGORA, symbol = :require_algorithm_support, rendered = "SpaceAGORA.require_algorithm_support"),
        (owner = :SpaceAGORA, symbol = :initialize_algorithm, rendered = "SpaceAGORA.initialize_algorithm"),
        (owner = :SpaceAGORA, symbol = :propose!, rendered = "SpaceAGORA.propose!"),
        (owner = :SpaceAGORA, symbol = :propose_revision!, rendered = "SpaceAGORA.propose_revision!"),
        (owner = :SpaceAGORA, symbol = :consumes_rng, rendered = "SpaceAGORA.consumes_rng"),
        (owner = :SpaceAGORA, symbol = :validate_capabilities, rendered = "SpaceAGORA.validate_capabilities"),
        (owner = :SpaceAGORA, symbol = :validate_proposal, rendered = "SpaceAGORA.validate_proposal"),
        (owner = :SpaceAGORA, symbol = :install!, rendered = "SpaceAGORA.install!"),
        (owner = :SpaceAGORA, symbol = :is_current, rendered = "SpaceAGORA.is_current"),
        (owner = :SpaceAGORA, symbol = :on_rejection, rendered = "SpaceAGORA.on_rejection"),
        (owner = :SpaceAGORA, symbol = :fallback_product, rendered = "SpaceAGORA.fallback_product"),
        (owner = :SpaceAGORA, symbol = :family_contract_version, rendered = "SpaceAGORA.family_contract_version"),
        (owner = :SpaceAGORA, symbol = :gnc_stream_words, rendered = "SpaceAGORA.gnc_stream_words"),
        (owner = :SpaceAGORA, symbol = :AlgorithmRuntime, rendered = "SpaceAGORA.AlgorithmRuntime"),
        (owner = :SpaceAGORA, symbol = :run_algorithm_request!, rendered = "SpaceAGORA.run_algorithm_request!"),
        (owner = :SpaceAGORA, symbol = :algorithm_runtime, rendered = "SpaceAGORA.algorithm_runtime"),
        (owner = :SpaceAGORA, symbol = :bound_algorithm, rendered = "SpaceAGORA.bound_algorithm"),
        (owner = :SpaceAGORA, symbol = :bound_family, rendered = "SpaceAGORA.bound_family"),
        (owner = :SpaceAGORA, symbol = :bound_link, rendered = "SpaceAGORA.bound_link"),
        (owner = :SpaceAGORA, symbol = :build_algorithm_request, rendered = "SpaceAGORA.build_algorithm_request"),
        (owner = :SpaceAGORA, symbol = :trusted_queries, rendered = "SpaceAGORA.trusted_queries"),
        (owner = :SpaceAGORA, symbol = :fail_algorithm_request!, rendered = "SpaceAGORA.fail_algorithm_request!"),
        (owner = :SpaceAGORA, symbol = :algorithm_result_type, rendered = "SpaceAGORA.algorithm_result_type"),
        (owner = :SpaceAGORA, symbol = :is_revision_request, rendered = "SpaceAGORA.is_revision_request"),
        (owner = :SpaceAGORA, symbol = :proposed_product, rendered = "SpaceAGORA.proposed_product"),
        (owner = :SpaceAGORA, symbol = :installed_record, rendered = "SpaceAGORA.installed_record"),
        (owner = :SpaceAGORA, symbol = :rejected_record, rendered = "SpaceAGORA.rejected_record"),
    ]),
    (title = "GNC runtime", items = [
        (owner = :SpaceAGORA, symbol = :tick_schedule, rendered = "SpaceAGORA.tick_schedule"),
        (owner = :SpaceAGORA, symbol = :gnc_callbacks, rendered = "SpaceAGORA.gnc_callbacks"),
        (owner = :SpaceAGORA, symbol = :register_tstops!, rendered = "SpaceAGORA.register_tstops!"),
        (owner = :SpaceAGORA, symbol = :gnc_clock_started!, rendered = "SpaceAGORA.gnc_clock_started!"),
        (owner = :SpaceAGORA, symbol = :gnc_period_s, rendered = "SpaceAGORA.gnc_period_s"),
    ]),
    (
        title = "Simulation Setup",
        items = [
            (owner = :SpaceAGORA, symbol = :SimulationConfiguration, rendered = "SpaceAGORA.SimulationConfiguration"),
            (owner = :SpaceAGORA, symbol = :MissionConfiguration, rendered = "SpaceAGORA.MissionConfiguration"),
            (owner = :SpaceAGORA, symbol = :MissionType, rendered = "SpaceAGORA.MissionType"),
            (owner = :SpaceAGORA, symbol = :MissionTime, rendered = "SpaceAGORA.MissionTime"),
            (owner = :SpaceAGORA, symbol = :MissionOrbits, rendered = "SpaceAGORA.MissionOrbits"),
            (owner = :SpaceAGORA, symbol = :EnvironmentModel, rendered = "SpaceAGORA.EnvironmentModel"),
            (owner = :SpaceAGORA, symbol = :DynamicsModel, rendered = "SpaceAGORA.DynamicsModel"),
            (owner = :SpaceAGORA, symbol = :GuidanceModel, rendered = "SpaceAGORA.GuidanceModel"),
            (owner = :SpaceAGORA, symbol = :NavigationModel, rendered = "SpaceAGORA.NavigationModel"),
            (owner = :SpaceAGORA, symbol = :ControlModel, rendered = "SpaceAGORA.ControlModel"),
            (owner = :SpaceAGORA, symbol = :SpacecraftModel, rendered = "SpaceAGORA.SpacecraftModel"),
            (owner = :SpaceAGORA, symbol = :Link, rendered = "SpaceAGORA.Link"),
            (owner = :SpaceAGORA, symbol = :Joint, rendered = "SpaceAGORA.Joint"),
            (owner = :SpaceAGORA, symbol = :InitialCondition, rendered = "SpaceAGORA.InitialCondition"),
            (owner = :SpaceAGORA, symbol = :CartesianInitialCondition, rendered = "SpaceAGORA.CartesianInitialCondition"),
            (owner = :SpaceAGORA, symbol = :InitialTime, rendered = "SpaceAGORA.InitialTime"),
            (owner = :SpaceAGORA, symbol = :IntegrationTolerances, rendered = "SpaceAGORA.IntegrationTolerances"),
            (owner = :SpaceAGORA, symbol = :SaveField, rendered = "SpaceAGORA.SaveField"),
            (owner = :SpaceAGORA, symbol = :default_save_fields, rendered = "SpaceAGORA.default_save_fields"),
            (owner = :SpaceAGORA, symbol = :MaxwellianHeat, rendered = "SpaceAGORA.MaxwellianHeat"),
            (owner = :SpaceAGORA, symbol = :AerodynamicCoefficientfM, rendered = "SpaceAGORA.AerodynamicCoefficientfM"),
            (owner = :SpaceAGORA, symbol = :InverseSquaredGravityModel, rendered = "SpaceAGORA.InverseSquaredGravityModel"),
            (owner = :SpaceAGORA, symbol = :InverseSquaredJ2GravityModel, rendered = "SpaceAGORA.InverseSquaredJ2GravityModel"),
            (owner = :SpaceAGORA, symbol = :NBodyGravityModel, rendered = "SpaceAGORA.NBodyGravityModel"),
            (owner = :SpaceAGORA, symbol = :GravitationalHarmonicsModel, rendered = "SpaceAGORA.GravitationalHarmonicsModel"),
            (owner = :SpaceAGORA, symbol = :make_example_config, rendered = "SpaceAGORA.make_example_config"),
            (owner = :SpaceAGORA, symbol = :make_three_body_spacecraft, rendered = "SpaceAGORA.make_three_body_spacecraft")
        ]
    ),
    (
        title = "Simulation",
        items = [
            (owner = :SpaceAGORA, symbol = :run_simulation, rendered = "SpaceAGORA.run_simulation"),
            (owner = :SpaceAGORA, symbol = :SimulationResults, rendered = "SpaceAGORA.SimulationResults"),
            (owner = :SpaceAGORA, symbol = :MissionPlan, rendered = "SpaceAGORA.MissionPlan"),
            (owner = :SpaceAGORA, symbol = :MissionPhase, rendered = "SpaceAGORA.MissionPhase"),
            (owner = :SpaceAGORA, symbol = :MissionSegment, rendered = "SpaceAGORA.MissionSegment"),
            (owner = :SpaceAGORA, symbol = :MissionDeadline, rendered = "SpaceAGORA.MissionDeadline"),
            (owner = :SpaceAGORA, symbol = :MissionDuration, rendered = "SpaceAGORA.MissionDuration"),
            (owner = :SpaceAGORA, symbol = :MissionCondition, rendered = "SpaceAGORA.MissionCondition"),
            (owner = :SpaceAGORA, symbol = :MissionEvent, rendered = "SpaceAGORA.MissionEvent"),
            (owner = :SpaceAGORA, symbol = :MissionRun, rendered = "SpaceAGORA.MissionRun"),
            (owner = :SpaceAGORA, symbol = :run_mission, rendered = "SpaceAGORA.run_mission"),
            (owner = :SpaceAGORA, symbol = :StateAnchor, rendered = "SpaceAGORA.StateAnchor"),
            (owner = :SpaceAGORA, symbol = :get_state_anchor_callback, rendered = "SpaceAGORA.get_state_anchor_callback"),
            (owner = :SpaceAGORA, symbol = :prewarm_nbody_ephemeris_cache, rendered = "SpaceAGORA.prewarm_nbody_ephemeris_cache"),
            (owner = :SpaceAGORA, symbol = :load_nbody_ephemeris_cache!, rendered = "SpaceAGORA.load_nbody_ephemeris_cache!"),
            (owner = :SimulationCampaigns, symbol = :MonteCarloSpec, rendered = "SpaceAGORA.MonteCarloSpec"),
            (owner = :SimulationCampaigns, symbol = :MonteCarloSampleResult, rendered = "SpaceAGORA.MonteCarloSampleResult"),
            (owner = :SimulationCampaigns, symbol = :MonteCarloResult, rendered = "SpaceAGORA.MonteCarloResult"),
            (owner = :SimulationCampaigns, symbol = :run_monte_carlo, rendered = "SpaceAGORA.run_monte_carlo"),
            (owner = :SimulationCampaigns, symbol = :run_constellation_ensemble, rendered = "SpaceAGORA.run_constellation_ensemble")
        ]
    ),
    (
        title = "Optional Adapter Services",
        items = [
            (owner = :SpaceAGORA, symbol = :visualizer_available, rendered = "SpaceAGORA.visualizer_available"),
            (owner = :SpaceAGORA, symbol = :VisualizerUnavailableError, rendered = "SpaceAGORA.VisualizerUnavailableError"),
            (owner = :SpaceAGORA, symbol = :export_viewer_payload, rendered = "SpaceAGORA.export_viewer_payload"),
            (owner = :SpaceAGORA, symbol = :ExternalPropagationServices, rendered = "SpaceAGORA.ExternalPropagationServices"),
            (owner = :SpaceAGORA, symbol = :EDGServices, rendered = "SpaceAGORA.EDGServices")
        ]
    ),
    (
        title = "Runtime Configuration",
        items = [
            (owner = :SpaceAGORA, symbol = :SimulationSettings, rendered = "SpaceAGORA.SimulationSettings"),
            (owner = :SimulationEngine, symbol = :ParallelConfig, rendered = "SpaceAGORA.ParallelConfig"),
            (owner = :SimulationEngine, symbol = :SolverConfig, rendered = "SpaceAGORA.SolverConfig"),
            (owner = :SimulationEngine, symbol = :RuntimePolicyConfig, rendered = "SpaceAGORA.RuntimePolicyConfig"),
            (owner = :SimulationEngine, symbol = :ArtifactConfig, rendered = "SpaceAGORA.ArtifactConfig"),
            (owner = :SimulationEngine, symbol = :SimulationEngineConfig, rendered = "SpaceAGORA.SimulationEngineConfig"),
            (owner = :SimulationEngine, symbol = :simulation_engine_config_from_env, rendered = "SpaceAGORA.simulation_engine_config_from_env")
        ]
    ),
    (
        title = "Mesh Aerodynamics",
        items = [
            (owner = :SpaceAGORA, symbol = :MeshAeroPanels, rendered = "SpaceAGORA.MeshAeroPanels"),
            (owner = :SpaceAGORA, symbol = :MeshAeroSurrogate, rendered = "SpaceAGORA.MeshAeroSurrogate"),
            (owner = :SpaceAGORA, symbol = :AerodynamicCoefficientMeshSurrogate, rendered = "SpaceAGORA.AerodynamicCoefficientMeshSurrogate"),
            (owner = :SpaceAGORA, symbol = :MESH_AERO_MAX_DEGREE, rendered = "SpaceAGORA.MESH_AERO_MAX_DEGREE"),
            (owner = :SpaceAGORA, symbol = :mesh_aero_panels, rendered = "SpaceAGORA.mesh_aero_panels"),
            (owner = :SpaceAGORA, symbol = :panel_aero_coefficients, rendered = "SpaceAGORA.panel_aero_coefficients"),
            (owner = :SpaceAGORA, symbol = :panel_aero_coefficients_split, rendered = "SpaceAGORA.panel_aero_coefficients_split"),
            (owner = :SpaceAGORA, symbol = :panel_shadow_mask, rendered = "SpaceAGORA.panel_shadow_mask"),
            (owner = :SpaceAGORA, symbol = :panel_projected_area, rendered = "SpaceAGORA.panel_projected_area"),
            (owner = :SpaceAGORA, symbol = :fit_mesh_aero_surrogate, rendered = "SpaceAGORA.fit_mesh_aero_surrogate"),
            (owner = :SpaceAGORA, symbol = :mesh_aero_coefficients, rendered = "SpaceAGORA.mesh_aero_coefficients"),
            (owner = :SpaceAGORA, symbol = :write_mesh_aero_surrogate, rendered = "SpaceAGORA.write_mesh_aero_surrogate"),
            (owner = :SpaceAGORA, symbol = :read_mesh_aero_surrogate, rendered = "SpaceAGORA.read_mesh_aero_surrogate")
        ]
    ),
    (
        title = "No-GRAM and Baseline Models",
        items = [
            (owner = :SpaceAGORA, symbol = :NoAtmosphereModel, rendered = "SpaceAGORA.NoAtmosphereModel"),
            (owner = :SpaceAGORA, symbol = :ExponentialAtmosphereModel, rendered = "SpaceAGORA.ExponentialAtmosphereModel"),
            (owner = :SpaceAGORA, symbol = :PiecewiseExponentialAtmosphereModel, rendered = "SpaceAGORA.PiecewiseExponentialAtmosphereModel"),
            (owner = :SpaceAGORA, symbol = :GRAMGridAtmosphereModel, rendered = "SpaceAGORA.GRAMGridAtmosphereModel"),
            (owner = :SpaceAGORA, symbol = :GRAMNearSurfaceAtmosphereModel, rendered = "SpaceAGORA.GRAMNearSurfaceAtmosphereModel"),
            (owner = :SpaceAGORA, symbol = :CombinedAtmosphereModel, rendered = "SpaceAGORA.CombinedAtmosphereModel"),
            (owner = :SpaceAGORA, symbol = :SurrogatePresetResolution, rendered = "SpaceAGORA.SurrogatePresetResolution"),
            (owner = :SpaceAGORA, symbol = :available_surrogate_presets, rendered = "SpaceAGORA.available_surrogate_presets"),
            (owner = :SpaceAGORA, symbol = :resolve_surrogate_preset, rendered = "SpaceAGORA.resolve_surrogate_preset"),
            (owner = :SpaceAGORA, symbol = :surrogate_preset_model, rendered = "SpaceAGORA.surrogate_preset_model"),
            (owner = :SpaceAGORA, symbol = :atmosphere_provenance, rendered = "SpaceAGORA.atmosphere_provenance"),
            (owner = :SpaceAGORA, symbol = :odyssey_surrogate_assets, rendered = "SpaceAGORA.odyssey_surrogate_assets"),
            (owner = :SpaceAGORA, symbol = :NRLMSISE00AtmosphereModel, rendered = "SpaceAGORA.NRLMSISE00AtmosphereModel"),
            (owner = :SpaceAGORA, symbol = :init_nrlmsise_space_indices!, rendered = "SpaceAGORA.init_nrlmsise_space_indices!"),
            (owner = :SpaceAGORA, symbol = :SimpleEphemeridesModel, rendered = "SpaceAGORA.SimpleEphemeridesModel"),
            (owner = :SpaceAGORA, symbol = :make_no_gram_planet, rendered = "SpaceAGORA.make_no_gram_planet"),
            (owner = :SpaceAGORA, symbol = :make_no_gram_density_model, rendered = "SpaceAGORA.make_no_gram_density_model"),
            (owner = :SpaceAGORA, symbol = :make_no_gram_environment, rendered = "SpaceAGORA.make_no_gram_environment")
        ]
    ),
    (
        title = "Plume Interaction",
        items = [
            (owner = :SpaceAGORA, symbol = :PlumeSurfaceConfig, rendered = "SpaceAGORA.PlumeSurfaceConfig"),
            (owner = :SpaceAGORA, symbol = :PlumeSurfaceInteractionModel, rendered = "SpaceAGORA.PlumeSurfaceInteractionModel"),
            (owner = :SpaceAGORA, symbol = :PlumeSurfaceState, rendered = "SpaceAGORA.PlumeSurfaceState"),
            (owner = :SpaceAGORA, symbol = :plume_surface_footprint, rendered = "SpaceAGORA.plume_surface_footprint"),
            (owner = :SpaceAGORA, symbol = :plume_erosion_onset_height, rendered = "SpaceAGORA.plume_erosion_onset_height"),
            (owner = :SpaceAGORA, symbol = :plume_ground_effect_force, rendered = "SpaceAGORA.plume_ground_effect_force"),
            (owner = :SpaceAGORA, symbol = :plume_quantities, rendered = "SpaceAGORA.plume_quantities")
        ]
    ),
    (
        title = "Lunar Landing",
        items = [
            (owner = :SpaceAGORA, symbol = :DescentPhaseTargets, rendered = "SpaceAGORA.DescentPhaseTargets"),
            (owner = :SpaceAGORA, symbol = :apollo11_descent_targets, rendered = "SpaceAGORA.apollo11_descent_targets"),
            (owner = :SpaceAGORA, symbol = :ApolloDescentConfig, rendered = "SpaceAGORA.ApolloDescentConfig"),
            (owner = :SpaceAGORA, symbol = :ApolloDescentState, rendered = "SpaceAGORA.ApolloDescentState"),
            (owner = :SpaceAGORA, symbol = :ApolloDescentGuidanceModel, rendered = "SpaceAGORA.ApolloDescentGuidanceModel"),
            (owner = :SpaceAGORA, symbol = :ApolloDescentControlConfig, rendered = "SpaceAGORA.ApolloDescentControlConfig"),
            (owner = :SpaceAGORA, symbol = :ApolloDescentControlModel, rendered = "SpaceAGORA.ApolloDescentControlModel"),
            (owner = :SpaceAGORA, symbol = :descent_attitude_command, rendered = "SpaceAGORA.descent_attitude_command")
        ]
    ),
    (
        title = "Terrain Queries",
        items = [
            (owner = :SpaceAGORA, symbol = :AbstractTerrainModel, rendered = "SpaceAGORA.AbstractTerrainModel"),
            (owner = :SpaceAGORA, symbol = :NoTerrainModel, rendered = "SpaceAGORA.NoTerrainModel"),
            (owner = :SpaceAGORA, symbol = :DEMGrid, rendered = "SpaceAGORA.DEMGrid"),
            (owner = :SpaceAGORA, symbol = :DEMTerrainModel, rendered = "SpaceAGORA.DEMTerrainModel"),
            (owner = :SpaceAGORA, symbol = :terrain_height, rendered = "SpaceAGORA.terrain_height"),
            (owner = :SpaceAGORA, symbol = :terrain_radius, rendered = "SpaceAGORA.terrain_radius"),
            (owner = :SpaceAGORA, symbol = :load_dem_grid, rendered = "SpaceAGORA.load_dem_grid"),
            (owner = :SpaceAGORA, symbol = :load_site_terrain, rendered = "SpaceAGORA.load_site_terrain"),
            (owner = :SpaceAGORA, symbol = :dem_grid_covers, rendered = "SpaceAGORA.dem_grid_covers")
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
            (owner = :SpaceAGORA, symbol = :AbstractNavigationModel, rendered = "SpaceAGORA.AbstractNavigationModel"),
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
            (owner = :SpaceAGORA, symbol = :wrench, rendered = "SpaceAGORA.wrench"),
            (owner = :SpaceAGORA, symbol = :environment_requirements, rendered = "SpaceAGORA.environment_requirements"),
            (owner = :SpaceAGORA, symbol = :solver_partition, rendered = "SpaceAGORA.solver_partition"),
            (owner = :SpaceAGORA, symbol = :gravity_backbone_structure, rendered = "SpaceAGORA.gravity_backbone_structure"),
            (owner = :SpaceAGORA, symbol = :gravity_backbone_acceleration_ii, rendered = "SpaceAGORA.gravity_backbone_acceleration_ii"),
            (owner = :SpaceAGORA, symbol = :gravity_backbone_kick_structure, rendered = "SpaceAGORA.gravity_backbone_kick_structure"),
            (owner = :SpaceAGORA, symbol = :gravity_backbone_kick_acceleration_ii, rendered = "SpaceAGORA.gravity_backbone_kick_acceleration_ii"),
            (owner = :SpaceAGORA, symbol = :getDensity, rendered = "SpaceAGORA.getDensity"),
            (owner = :SpaceAGORA, symbol = :getDensityBatch!, rendered = "SpaceAGORA.getDensityBatch!"),
            (owner = :SpaceAGORA, symbol = :with_density_model_epoch, rendered = "SpaceAGORA.with_density_model_epoch"),
            (owner = :SpaceAGORA, symbol = :calcControlEffect!, rendered = "SpaceAGORA.calcControlEffect!"),
            (owner = :SpaceAGORA, symbol = :control_wrench, rendered = "SpaceAGORA.control_wrench"),
            (owner = :SpaceAGORA, symbol = :control_propellant_rate, rendered = "SpaceAGORA.control_propellant_rate"),
            (owner = :SpaceAGORA, symbol = :control_wheel_torque, rendered = "SpaceAGORA.control_wheel_torque"),
            (owner = :SpaceAGORA, symbol = :control_thruster_levels, rendered = "SpaceAGORA.control_thruster_levels"),
            (owner = :SpaceAGORA, symbol = :touchdown_spec, rendered = "SpaceAGORA.touchdown_spec"),
            (owner = :SpaceAGORA, symbol = :calcGuidanceEffect!, rendered = "SpaceAGORA.calcGuidanceEffect!"),
            (owner = :SpaceAGORA, symbol = :calcNavigationEffect!, rendered = "SpaceAGORA.calcNavigationEffect!"),
            (owner = :SpaceAGORA, symbol = :getHeatRate, rendered = "SpaceAGORA.getHeatRate"),
            (owner = :SpaceAGORA, symbol = :bind_spacecraft, rendered = "SpaceAGORA.bind_spacecraft")
        ]
    ),
    (
        title = "Deprecated hooks",
        items = [
            (owner = :SpaceAGORA, symbol = :calcForceTorque, rendered = "SpaceAGORA.calcForceTorque"),
            (owner = :SpaceAGORA, symbol = :calcControlForceTorque, rendered = "SpaceAGORA.calcControlForceTorque"),
            (owner = :SpaceAGORA, symbol = :calcControlMassFlowRate, rendered = "SpaceAGORA.calcControlMassFlowRate"),
            (owner = :SpaceAGORA, symbol = :calcReactionWheelTorque, rendered = "SpaceAGORA.calcReactionWheelTorque")
        ]
    ),
    (
        title = "Actuators",
        items = [
            (owner = :SpaceAGORA, symbol = :AbstractActuator, rendered = "SpaceAGORA.AbstractActuator"),
            (owner = :SpaceAGORA, symbol = :Thruster, rendered = "SpaceAGORA.Thruster"),
            (owner = :SpaceAGORA, symbol = :ReactionWheel, rendered = "SpaceAGORA.ReactionWheel"),
            (owner = :SpaceAGORA, symbol = :Magnetorquer, rendered = "SpaceAGORA.Magnetorquer"),
            (owner = :SpaceAGORA, symbol = :LaserTerminal, rendered = "SpaceAGORA.LaserTerminal"),
            (owner = :SpaceAGORA, symbol = :PanelHinge, rendered = "SpaceAGORA.PanelHinge"),
            (owner = :SpaceAGORA, symbol = :ThrusterGimbal, rendered = "SpaceAGORA.ThrusterGimbal"),
            (owner = :SpaceAGORA, symbol = :thrust_direction, rendered = "SpaceAGORA.thrust_direction"),
            (owner = :SpaceAGORA, symbol = :ControlMomentGyro, rendered = "SpaceAGORA.ControlMomentGyro"),
            (owner = :SpaceAGORA, symbol = :cmg_momentum, rendered = "SpaceAGORA.cmg_momentum"),
            (owner = :SpaceAGORA, symbol = :cmg_torque_axis, rendered = "SpaceAGORA.cmg_torque_axis"),
            (owner = :SpaceAGORA, symbol = :cmg_angle_saturate, rendered = "SpaceAGORA.cmg_angle_saturate"),
            (owner = :SpaceAGORA, symbol = :cmg_steering, rendered = "SpaceAGORA.cmg_steering"),
            (owner = :SpaceAGORA, symbol = :ActuatorContext, rendered = "SpaceAGORA.ActuatorContext"),
            (owner = :SpaceAGORA, symbol = :actuator_wrench, rendered = "SpaceAGORA.actuator_wrench"),
            (owner = :SpaceAGORA, symbol = :actuator_limits, rendered = "SpaceAGORA.actuator_limits"),
            (owner = :SpaceAGORA, symbol = :actuator_saturate, rendered = "SpaceAGORA.actuator_saturate"),
            (owner = :SpaceAGORA, symbol = :actuator_consumption, rendered = "SpaceAGORA.actuator_consumption"),
            (owner = :SpaceAGORA, symbol = :actuator_fraction, rendered = "SpaceAGORA.actuator_fraction"),
            (owner = :SpaceAGORA, symbol = :actuator_unit_wrench, rendered = "SpaceAGORA.actuator_unit_wrench"),
            (owner = :SpaceAGORA, symbol = :actuation_mode, rendered = "SpaceAGORA.actuation_mode"),
            (owner = :SpaceAGORA, symbol = :actuator_levels, rendered = "SpaceAGORA.actuator_levels"),
            (owner = :SpaceAGORA, symbol = :add_actuator!, rendered = "SpaceAGORA.add_actuator!"),
            (owner = :SpaceAGORA, symbol = :mounted_actuators, rendered = "SpaceAGORA.mounted_actuators"),
            (owner = :SpaceAGORA, symbol = :actuator_context, rendered = "SpaceAGORA.actuator_context"),
            (owner = :SpaceAGORA, symbol = :allocation_matrix, rendered = "SpaceAGORA.allocation_matrix"),
            (owner = :SpaceAGORA, symbol = :mounted_wrench, rendered = "SpaceAGORA.mounted_wrench"),
            (owner = :SpaceAGORA, symbol = :allocate_actuators, rendered = "SpaceAGORA.allocate_actuators"),
            (owner = :SpaceAGORA, symbol = :six_axis_thrusters, rendered = "SpaceAGORA.six_axis_thrusters"),
            (owner = :SpaceAGORA, symbol = :propellant_mass_rate, rendered = "SpaceAGORA.propellant_mass_rate"),
            (owner = :SpaceAGORA, symbol = :ActuatorRequirements, rendered = "SpaceAGORA.ActuatorRequirements"),
            (owner = :SpaceAGORA, symbol = :actuator_requirements, rendered = "SpaceAGORA.actuator_requirements"),
            (owner = :SpaceAGORA, symbol = :applied_actuator_levels, rendered = "SpaceAGORA.applied_actuator_levels"),
            (owner = :SpaceAGORA, symbol = :mounted_power, rendered = "SpaceAGORA.mounted_power"),
            (owner = :SpaceAGORA, symbol = :SolarArray, rendered = "SpaceAGORA.SolarArray"),
            (owner = :SpaceAGORA, symbol = :Battery, rendered = "SpaceAGORA.Battery"),
            (owner = :SpaceAGORA, symbol = :PowerLoad, rendered = "SpaceAGORA.PowerLoad"),
            (owner = :SpaceAGORA, symbol = :add_power!, rendered = "SpaceAGORA.add_power!"),
            (owner = :SpaceAGORA, symbol = :ActuatorFault, rendered = "SpaceAGORA.ActuatorFault"),
            (owner = :SpaceAGORA, symbol = :add_fault!, rendered = "SpaceAGORA.add_fault!")
        ]
    ),
    (
        title = "Sensors",
        items = [
            (owner = :SpaceAGORA, symbol = :AbstractSensor, rendered = "SpaceAGORA.AbstractSensor"),
            (owner = :SpaceAGORA, symbol = :RateGyro, rendered = "SpaceAGORA.RateGyro"),
            (owner = :SpaceAGORA, symbol = :StarTracker, rendered = "SpaceAGORA.StarTracker"),
            (owner = :SpaceAGORA, symbol = :Magnetometer, rendered = "SpaceAGORA.Magnetometer"),
            (owner = :SpaceAGORA, symbol = :SensorTruth, rendered = "SpaceAGORA.SensorTruth"),
            (owner = :SpaceAGORA, symbol = :sensor_measure, rendered = "SpaceAGORA.sensor_measure"),
            (owner = :SpaceAGORA, symbol = :sensor_dimension, rendered = "SpaceAGORA.sensor_dimension"),
            (owner = :SpaceAGORA, symbol = :sensor_requirements, rendered = "SpaceAGORA.sensor_requirements"),
            (owner = :SpaceAGORA, symbol = :sensor_period_s, rendered = "SpaceAGORA.sensor_period_s"),
            (owner = :SpaceAGORA, symbol = :initialize_sensor!, rendered = "SpaceAGORA.initialize_sensor!"),
            (owner = :SpaceAGORA, symbol = :add_sensor!, rendered = "SpaceAGORA.add_sensor!"),
            (owner = :SpaceAGORA, symbol = :mounted_sensors, rendered = "SpaceAGORA.mounted_sensors"),
            (owner = :SpaceAGORA, symbol = :sensor_measurement, rendered = "SpaceAGORA.sensor_measurement"),
            (owner = :SpaceAGORA, symbol = :sensor_measurement_time, rendered = "SpaceAGORA.sensor_measurement_time")
        ]
    ),
    (
        title = "Robotics and Planning",
        items = [
            (owner = :SpaceAGORA, symbol = :RobotArmBasePose, rendered = "SpaceAGORA.RobotArmBasePose"),
            (owner = :SpaceAGORA, symbol = :RobotArmLink, rendered = "SpaceAGORA.RobotArmLink"),
            (owner = :SpaceAGORA, symbol = :RobotArmJoint, rendered = "SpaceAGORA.RobotArmJoint"),
            (owner = :SpaceAGORA, symbol = :RobotArmModel, rendered = "SpaceAGORA.RobotArmModel"),
            (owner = :SpaceAGORA, symbol = :RobotArmPose, rendered = "SpaceAGORA.RobotArmPose"),
            (owner = :SpaceAGORA, symbol = :RobotArmState, rendered = "SpaceAGORA.RobotArmState"),
            (owner = :SpaceAGORA, symbol = :default_robot_arm_model, rendered = "SpaceAGORA.default_robot_arm_model"),
            (owner = :SpaceAGORA, symbol = :robot_arm_fk, rendered = "SpaceAGORA.robot_arm_fk"),
            (owner = :SpaceAGORA, symbol = :robot_arm_fk_state, rendered = "SpaceAGORA.robot_arm_fk_state"),
            (owner = :SpaceAGORA, symbol = :robot_arm_end_effector_pose, rendered = "SpaceAGORA.robot_arm_end_effector_pose"),
            (owner = :SpaceAGORA, symbol = :robot_arm_ik, rendered = "SpaceAGORA.robot_arm_ik"),
            (owner = :SpaceAGORA, symbol = :robot_arm_total_reach, rendered = "SpaceAGORA.robot_arm_total_reach"),
            (owner = :SpaceAGORA, symbol = :ClothArmBasePose, rendered = "SpaceAGORA.ClothArmBasePose"),
            (owner = :SpaceAGORA, symbol = :ClothArmLink, rendered = "SpaceAGORA.ClothArmLink"),
            (owner = :SpaceAGORA, symbol = :ClothArmJoint, rendered = "SpaceAGORA.ClothArmJoint"),
            (owner = :SpaceAGORA, symbol = :ClothArmModel, rendered = "SpaceAGORA.ClothArmModel"),
            (owner = :SpaceAGORA, symbol = :ClothArmPose, rendered = "SpaceAGORA.ClothArmPose"),
            (owner = :SpaceAGORA, symbol = :ClothArmState, rendered = "SpaceAGORA.ClothArmState"),
            (owner = :SpaceAGORA, symbol = :default_cloth_arm_model, rendered = "SpaceAGORA.default_cloth_arm_model"),
            (owner = :SpaceAGORA, symbol = :cloth_fk, rendered = "SpaceAGORA.cloth_fk"),
            (owner = :SpaceAGORA, symbol = :cloth_fk_state, rendered = "SpaceAGORA.cloth_fk_state"),
            (owner = :SpaceAGORA, symbol = :cloth_end_effector_pose, rendered = "SpaceAGORA.cloth_end_effector_pose"),
            (owner = :SpaceAGORA, symbol = :cloth_ik, rendered = "SpaceAGORA.cloth_ik"),
            (owner = :SpaceAGORA, symbol = :cloth_total_reach, rendered = "SpaceAGORA.cloth_total_reach"),
            (owner = :SpaceAGORA, symbol = :closest_surface_target, rendered = "SpaceAGORA.closest_surface_target"),
            (owner = :SpaceAGORA, symbol = :RobotArmPlannerConfig, rendered = "SpaceAGORA.RobotArmPlannerConfig"),
            (owner = :SpaceAGORA, symbol = :RobotArmPlan, rendered = "SpaceAGORA.RobotArmPlan"),
            (owner = :SpaceAGORA, symbol = :robot_arm_plan_sample, rendered = "SpaceAGORA.robot_arm_plan_sample"),
            (owner = :SpaceAGORA, symbol = :RobotArmSphereObstacle, rendered = "SpaceAGORA.RobotArmSphereObstacle"),
            (owner = :SpaceAGORA, symbol = :RobotArmHYPRConfig, rendered = "SpaceAGORA.RobotArmHYPRConfig"),
            (owner = :SpaceAGORA, symbol = :RobotArmHYPRResult, rendered = "SpaceAGORA.RobotArmHYPRResult"),
            (owner = :SpaceAGORA, symbol = :plan_robot_arm_motion, rendered = "SpaceAGORA.plan_robot_arm_motion"),
            (owner = :SpaceAGORA, symbol = :plan_robot_arm_motion_hypr, rendered = "SpaceAGORA.plan_robot_arm_motion_hypr"),
            (owner = :SpaceAGORA, symbol = :robot_arm_sample_hypr_path, rendered = "SpaceAGORA.robot_arm_sample_hypr_path"),
            (owner = :SpaceAGORA, symbol = :robot_arm_clearance_stats_from_samples, rendered = "SpaceAGORA.robot_arm_clearance_stats_from_samples"),
            (owner = :SpaceAGORA, symbol = :robot_arm_hypr_path_cost_components, rendered = "SpaceAGORA.robot_arm_hypr_path_cost_components")
        ]
    ),
    (
        title = "Compliant Multibody Dynamics",
        items = [
            (owner = :SpaceAGORA, symbol = :CompliantBody, rendered = "SpaceAGORA.CompliantBody"),
            (owner = :SpaceAGORA, symbol = :CompliantJoint, rendered = "SpaceAGORA.CompliantJoint"),
            (owner = :SpaceAGORA, symbol = :CompliantMultibodyModel, rendered = "SpaceAGORA.CompliantMultibodyModel"),
            (owner = :SpaceAGORA, symbol = :CompliantMultibodyTrajectory, rendered = "SpaceAGORA.CompliantMultibodyTrajectory"),
            (owner = :SpaceAGORA, symbol = :CompliantTopologyNode, rendered = "SpaceAGORA.CompliantTopologyNode"),
            (owner = :SpaceAGORA, symbol = :CompliantTopologyEdge, rendered = "SpaceAGORA.CompliantTopologyEdge"),
            (owner = :SpaceAGORA, symbol = :CompliantTopologyBuild, rendered = "SpaceAGORA.CompliantTopologyBuild"),
            (owner = :SpaceAGORA, symbol = :CompliantJointActuator, rendered = "SpaceAGORA.CompliantJointActuator"),
            (owner = :SpaceAGORA, symbol = :CompliantJointLoad, rendered = "SpaceAGORA.CompliantJointLoad"),
            (owner = :SpaceAGORA, symbol = :CompliantAttachment, rendered = "SpaceAGORA.CompliantAttachment"),
            (owner = :SpaceAGORA, symbol = :rectangular_prism_inertia, rendered = "SpaceAGORA.rectangular_prism_inertia"),
            (owner = :SpaceAGORA, symbol = :thin_panel_inertia, rendered = "SpaceAGORA.thin_panel_inertia"),
            (owner = :SpaceAGORA, symbol = :build_compliant_topology, rendered = "SpaceAGORA.build_compliant_topology"),
            (owner = :SpaceAGORA, symbol = :build_rectangular_compliant_grid, rendered = "SpaceAGORA.build_rectangular_compliant_grid"),
            (owner = :SpaceAGORA, symbol = :compliant_state_vector, rendered = "SpaceAGORA.compliant_state_vector"),
            (owner = :SpaceAGORA, symbol = :compliant_state_parts, rendered = "SpaceAGORA.compliant_state_parts"),
            (owner = :SpaceAGORA, symbol = :compliant_multibody_dynamics, rendered = "SpaceAGORA.compliant_multibody_dynamics"),
            (owner = :SpaceAGORA, symbol = :compliant_joint_loads, rendered = "SpaceAGORA.compliant_joint_loads"),
            (owner = :SpaceAGORA, symbol = :step_compliant_multibody_rk4, rendered = "SpaceAGORA.step_compliant_multibody_rk4"),
            (owner = :SpaceAGORA, symbol = :step_compliant_multibody_implicit_midpoint, rendered = "SpaceAGORA.step_compliant_multibody_implicit_midpoint"),
            (owner = :SpaceAGORA, symbol = :simulate_compliant_multibody, rendered = "SpaceAGORA.simulate_compliant_multibody"),
            (owner = :SpaceAGORA, symbol = :CompliantMultibody, rendered = "SpaceAGORA.CompliantMultibody"),
            (owner = :SpaceAGORA, symbol = :CompliantRobotArmDynamics, rendered = "SpaceAGORA.CompliantRobotArmDynamics"),
            (owner = :SpaceAGORA, symbol = :CompliantRobotArmReferenceState, rendered = "SpaceAGORA.CompliantRobotArmReferenceState"),
            (owner = :SpaceAGORA, symbol = :CompliantRobotArmSimulation, rendered = "SpaceAGORA.CompliantRobotArmSimulation"),
            (owner = :SpaceAGORA, symbol = :compliant_robot_arm_reference_state, rendered = "SpaceAGORA.compliant_robot_arm_reference_state"),
            (owner = :SpaceAGORA, symbol = :compliant_robot_arm_multibody, rendered = "SpaceAGORA.compliant_robot_arm_multibody"),
            (owner = :SpaceAGORA, symbol = :compliant_robot_arm_initial_state, rendered = "SpaceAGORA.compliant_robot_arm_initial_state"),
            (owner = :SpaceAGORA, symbol = :compliant_robot_arm_rest_quaternions, rendered = "SpaceAGORA.compliant_robot_arm_rest_quaternions"),
            (owner = :SpaceAGORA, symbol = :compliant_robot_arm_end_effector, rendered = "SpaceAGORA.compliant_robot_arm_end_effector"),
            (owner = :SpaceAGORA, symbol = :compliant_robot_arm_actuators, rendered = "SpaceAGORA.compliant_robot_arm_actuators"),
            (owner = :SpaceAGORA, symbol = :coupled_compliant_robot_arm_state_shape, rendered = "SpaceAGORA.coupled_compliant_robot_arm_state_shape"),
            (owner = :SpaceAGORA, symbol = :initialize_coupled_compliant_robot_arm_state!, rendered = "SpaceAGORA.initialize_coupled_compliant_robot_arm_state!"),
            (owner = :SpaceAGORA, symbol = :assign_coupled_compliant_robot_arm_rhs!, rendered = "SpaceAGORA.assign_coupled_compliant_robot_arm_rhs!"),
            (owner = :SpaceAGORA, symbol = :simulate_compliant_robot_arm_plan, rendered = "SpaceAGORA.simulate_compliant_robot_arm_plan"),
            (owner = :SpaceAGORA, symbol = :ClothRobotArmSimulation, rendered = "SpaceAGORA.ClothRobotArmSimulation"),
            (owner = :SpaceAGORA, symbol = :cloth_robot_arm_multibody, rendered = "SpaceAGORA.cloth_robot_arm_multibody"),
            (owner = :SpaceAGORA, symbol = :cloth_robot_arm_initial_state, rendered = "SpaceAGORA.cloth_robot_arm_initial_state"),
            (owner = :SpaceAGORA, symbol = :cloth_robot_arm_rest_quaternions, rendered = "SpaceAGORA.cloth_robot_arm_rest_quaternions"),
            (owner = :SpaceAGORA, symbol = :cloth_robot_arm_end_effector, rendered = "SpaceAGORA.cloth_robot_arm_end_effector"),
            (owner = :SpaceAGORA, symbol = :simulate_cloth_robot_arm_plan, rendered = "SpaceAGORA.simulate_cloth_robot_arm_plan"),
            (owner = :SpaceAGORA, symbol = :cloth_robot_arm_actuators, rendered = "SpaceAGORA.cloth_robot_arm_actuators"),
            (owner = :SpaceAGORA, symbol = :coupled_cloth_robot_arm_state_shape, rendered = "SpaceAGORA.coupled_cloth_robot_arm_state_shape"),
            (owner = :SpaceAGORA, symbol = :initialize_coupled_cloth_robot_arm_state!, rendered = "SpaceAGORA.initialize_coupled_cloth_robot_arm_state!"),
            (owner = :SpaceAGORA, symbol = :assign_coupled_cloth_robot_arm_rhs!, rendered = "SpaceAGORA.assign_coupled_cloth_robot_arm_rhs!")
        ]
    ),
    (
        title = "Robot Arm Control",
        items = [
            (owner = :SpaceAGORA, symbol = :RobotArmHeldActuation, rendered = "SpaceAGORA.RobotArmHeldActuation"),
            (owner = :SpaceAGORA, symbol = :RobotArmJointMPCController, rendered = "SpaceAGORA.RobotArmJointMPCController"),
            (owner = :SpaceAGORA, symbol = :RobotArmControlEffector, rendered = "SpaceAGORA.RobotArmControlEffector"),
            (owner = :SpaceAGORA, symbol = :RobotArmReactionEffector, rendered = "SpaceAGORA.RobotArmReactionEffector"),
            (owner = :SpaceAGORA, symbol = :init_robot_arm_joint_mpc, rendered = "SpaceAGORA.init_robot_arm_joint_mpc"),
            (owner = :SpaceAGORA, symbol = :robot_arm_joint_mpc_reference_preview, rendered = "SpaceAGORA.robot_arm_joint_mpc_reference_preview"),
            (owner = :SpaceAGORA, symbol = :robot_arm_joint_mpc_control, rendered = "SpaceAGORA.robot_arm_joint_mpc_control"),
            (owner = :SpaceAGORA, symbol = :robot_arm_measured_joint_state, rendered = "SpaceAGORA.robot_arm_measured_joint_state")
        ]
    ),
    (
        title = "Guidance and RPO Assets",
        items = [
            (owner = :SpaceAGORA, symbol = :ApoapsisTargetPeriapsisRaiseGuidanceModel, rendered = "SpaceAGORA.ApoapsisTargetPeriapsisRaiseGuidanceModel"),
            (owner = :SpaceAGORA, symbol = :AerobrakingEnergyDepletionConfig, rendered = "SpaceAGORA.AerobrakingEnergyDepletionConfig"),
            (owner = :SpaceAGORA, symbol = :AerobrakingEnergyDepletionControlModel, rendered = "SpaceAGORA.AerobrakingEnergyDepletionControlModel"),
            (owner = :SpaceAGORA, symbol = :AerobrakingEnergyDepletionGuidanceModel, rendered = "SpaceAGORA.AerobrakingEnergyDepletionGuidanceModel"),
            (owner = :SpaceAGORA, symbol = :AerobrakingEnergyDepletionState, rendered = "SpaceAGORA.AerobrakingEnergyDepletionState"),
            (owner = :SpaceAGORA, symbol = :SolarPanelAngleOfAttackControlModel, rendered = "SpaceAGORA.SolarPanelAngleOfAttackControlModel"),
            (owner = :SpaceAGORA, symbol = :station_geometry_path, rendered = "SpaceAGORA.station_geometry_path"),
            (owner = :SpaceAGORA, symbol = :station_cad_path, rendered = "SpaceAGORA.station_cad_path"),
            (owner = :SpaceAGORA, symbol = :load_rpo_station_pointcloud, rendered = "SpaceAGORA.load_rpo_station_pointcloud"),
            (owner = :SpaceAGORA, symbol = :load_rpo_station_cad_triangles, rendered = "SpaceAGORA.load_rpo_station_cad_triangles"),
            (owner = :SpaceAGORA, symbol = :load_rpo_station_cad_pointcloud, rendered = "SpaceAGORA.load_rpo_station_cad_pointcloud")
        ]
    ),
    (
        title = "Visualization Scene",
        items = [
            (owner = :SpaceAGORA, symbol = :VisualizationScene, rendered = "SpaceAGORA.VisualizationScene"),
            (owner = :SpaceAGORA, symbol = :PlanetSpec, rendered = "SpaceAGORA.PlanetSpec"),
            (owner = :SpaceAGORA, symbol = :SpacecraftGeometry, rendered = "SpaceAGORA.SpacecraftGeometry"),
            (owner = :SpaceAGORA, symbol = :LinkBox, rendered = "SpaceAGORA.LinkBox"),
            (owner = :SpaceAGORA, symbol = :AtmosphereSpec, rendered = "SpaceAGORA.AtmosphereSpec"),
            (owner = :SpaceAGORA, symbol = :atmosphere_spec, rendered = "SpaceAGORA.atmosphere_spec"),
            (owner = :SpaceAGORA, symbol = :ArmGeometry, rendered = "SpaceAGORA.ArmGeometry"),
            (owner = :SpaceAGORA, symbol = :arm_geometry, rendered = "SpaceAGORA.arm_geometry"),
            (owner = :SpaceAGORA, symbol = :load_model_triangles, rendered = "SpaceAGORA.load_model_triangles"),
            (owner = :SpaceAGORA, symbol = :model_bounding_box, rendered = "SpaceAGORA.model_bounding_box"),
            (owner = :SpaceAGORA, symbol = :sample_model_pointcloud, rendered = "SpaceAGORA.sample_model_pointcloud"),
            (owner = :SpaceAGORA, symbol = :articulate_triangles, rendered = "SpaceAGORA.articulate_triangles"),
            (owner = :SpaceAGORA, symbol = :articulation_payload, rendered = "SpaceAGORA.articulation_payload"),
            (owner = :SpaceAGORA, symbol = :spacecraft_geometry, rendered = "SpaceAGORA.spacecraft_geometry"),
            (owner = :SpaceAGORA, symbol = :planet_spec, rendered = "SpaceAGORA.planet_spec"),
            (owner = :SpaceAGORA, symbol = :planet_rotation_table, rendered = "SpaceAGORA.planet_rotation_table"),
            (owner = :SpaceAGORA, symbol = :build_visualization_scene, rendered = "SpaceAGORA.build_visualization_scene"),
            (owner = :SpaceAGORA, symbol = :visualization_scene_path, rendered = "SpaceAGORA.visualization_scene_path"),
            (owner = :SpaceAGORA, symbol = :write_visualization_scene, rendered = "SpaceAGORA.write_visualization_scene"),
            (owner = :SpaceAGORA, symbol = :read_visualization_scene, rendered = "SpaceAGORA.read_visualization_scene"),
            (owner = :SpaceAGORA, symbol = :velocity_aligned_quaternion, rendered = "SpaceAGORA.velocity_aligned_quaternion"),
            (owner = :SpaceAGORA, symbol = :visualization_frame_budget, rendered = "SpaceAGORA.visualization_frame_budget"),
            (owner = :SpaceAGORA, symbol = :export_visualization, rendered = "SpaceAGORA.export_visualization"),
            (owner = :SpaceAGORA, symbol = :with_visualization_scene, rendered = "SpaceAGORA.with_visualization_scene"),
            (owner = :SpaceAGORA, symbol = :write_viewer_dev_payload, rendered = "SpaceAGORA.write_viewer_dev_payload"),
            (owner = :SpaceAGORA, symbol = :EnsembleSample, rendered = "SpaceAGORA.EnsembleSample"),
            (owner = :SpaceAGORA, symbol = :sample_results_directory, rendered = "SpaceAGORA.sample_results_directory"),
            (owner = :SpaceAGORA, symbol = :with_results_directory, rendered = "SpaceAGORA.with_results_directory"),
            (owner = :SpaceAGORA, symbol = :write_ensemble_manifest, rendered = "SpaceAGORA.write_ensemble_manifest"),
            (owner = :SpaceAGORA, symbol = :export_ensemble_visualization, rendered = "SpaceAGORA.export_ensemble_visualization"),
            (owner = :SimulationCampaigns, symbol = :run_monte_carlo_visualization, rendered = "SpaceAGORA.run_monte_carlo_visualization")
        ]
    ),
    (
        title = "Process Workers",
        items = [
            (owner = :ParallelProcess, symbol = :ProcessPool, rendered = "SpaceAGORA.ProcessPool"),
            (owner = :ParallelProcess, symbol = :PROCESS_WARMUP, rendered = "SpaceAGORA.PROCESS_WARMUP"),
            (owner = :ParallelProcess, symbol = :campaign_process_pool, rendered = "SpaceAGORA.campaign_process_pool"),
            (owner = :ParallelProcess, symbol = :ensure_process_workers!, rendered = "SpaceAGORA.ensure_process_workers!"),
            (owner = :ParallelProcess, symbol = :shutdown_process_pool!, rendered = "SpaceAGORA.shutdown_process_pool!"),
            (owner = :ParallelProcess, symbol = :adopt_process_workers!, rendered = "SpaceAGORA.adopt_process_workers!")
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
    ),
    # The Conformance module's own exports; they are not exported from the root.
    (
        title = "Contributor conformance",
        items = [
            (owner = :SpaceAGORA, symbol = :Conformance, rendered = "SpaceAGORA.Conformance"),
            (owner = :Conformance, symbol = :ConformanceCheck, rendered = "SpaceAGORA.Conformance.ConformanceCheck"),
            (owner = :Conformance, symbol = :ConformanceReport, rendered = "SpaceAGORA.Conformance.ConformanceReport"),
            (owner = :Conformance, symbol = :passed, rendered = "SpaceAGORA.Conformance.passed"),
            (owner = :Conformance, symbol = :ConformanceFixture, rendered = "SpaceAGORA.Conformance.ConformanceFixture"),
            (owner = :Conformance, symbol = :conformance_fixture, rendered = "SpaceAGORA.Conformance.conformance_fixture"),
            (owner = :Conformance, symbol = :check_force_torque_model, rendered = "SpaceAGORA.Conformance.check_force_torque_model"),
            (owner = :Conformance, symbol = :check_service_facade, rendered = "SpaceAGORA.Conformance.check_service_facade"),
            (owner = :Conformance, symbol = :test_report, rendered = "SpaceAGORA.Conformance.test_report"),
            (owner = :Conformance, symbol = :test_force_torque_model, rendered = "SpaceAGORA.Conformance.test_force_torque_model"),
            (owner = :Conformance, symbol = :test_service_facade, rendered = "SpaceAGORA.Conformance.test_service_facade")
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
    elseif owner === :ParallelProcess
        return getproperty(spaceagora, :ParallelProcess)
    elseif owner === :TelemetryVerification
        return getproperty(spaceagora, :TelemetryVerification)
    elseif owner === :SpaceAGORACLI
        return getproperty(spaceagora, :SpaceAGORACLI)
    elseif owner === :Conformance
        return getproperty(spaceagora, :Conformance)
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

function render_public_api_markdown(spaceagora::Module;
                                    sections=PUBLIC_API_SECTIONS,
                                    title::AbstractString="Public API")::String
    io = IOBuffer()
    println(io, "# $(title)")
    println(io)
    println(io, "This page documents the supported public interface available from `SpaceAGORA`.")
    println(io)
    for section in sections
        println(io, "## $(section.title)")
        println(io)
        # Julia stores a module alias's docstring in the target module. An alias
        # of an internal module therefore cannot use Documenter's @docs lookup,
        # which intentionally excludes those implementation modules.
        module_aliases = filter(section.items) do item
            owner = _owner_module(spaceagora, item.owner)
            value = getproperty(owner, item.symbol)
            value isa Module && parentmodule(value) !== owner
        end
        println(io, "```@docs")
        for item in section.items
            item in module_aliases || println(io, item.rendered)
        end
        println(io, "```")
        println(io)
        for item in module_aliases
            binding = Docs.Binding(_owner_module(spaceagora, item.owner), item.symbol)
            doc = Docs.doc(binding)
            doc === nothing && error("Undocumented module alias: $(item.rendered)")
            println(io, "### `$(item.rendered)`")
            println(io)
            show(io, MIME"text/markdown"(), doc)
            println(io)
        end
    end
    return String(take!(io))
end

# Public names, as (owner module, symbol), whose rendered documentation is still
# the "No documentation found" placeholder. Each is listed so that any new gap
# fails the build. The build also fails when a listed name no longer renders the
# placeholder, so documenting one requires removing it here.
const KNOWN_PLACEHOLDER_DOCS = Set([
    (:SpaceAGORA, :thin_panel_inertia),
    (:SpaceAGORA, :ApoapsisTargetPeriapsisRaiseGuidanceModel),
    (:SpaceAGORA, :station_geometry_path),
    (:SpaceAGORA, :station_cad_path),
    (:SpaceAGORA, :load_rpo_station_pointcloud),
    (:SpaceAGORA, :load_rpo_station_cad_triangles),
    (:SpaceAGORA, :load_rpo_station_cad_pointcloud),
])

_spec_key(spec) = (nameof(spec.owner_module), spec.symbol)

_is_placeholder_doc(doc) = doc === nothing || occursin("No documentation found", string(doc))

function undocumented_public_api_specs(spaceagora::Module)
    missing = NamedTuple[]
    for spec in public_api_specs(spaceagora)
        binding = Docs.Binding(spec.owner_module, spec.symbol)
        # `Docs.doc` returns a placeholder rather than `nothing` for an
        # undocumented name, and a re-export written as `@doc (@doc inner) name`
        # copies that placeholder into a real docstring, which `hasdoc` then
        # counts. Check the text as well.
        placeholder = _is_placeholder_doc(Docs.doc(binding))
        if !Docs.hasdoc(binding.mod, binding.var) ||
                (placeholder && !(_spec_key(spec) in KNOWN_PLACEHOLDER_DOCS))
            push!(missing, spec)
        end
    end
    return missing
end

"""
Entries of `KNOWN_PLACEHOLDER_DOCS` that no longer render the placeholder or are
no longer public names; each must be removed from the list.
"""
function stale_placeholder_entries(spaceagora::Module)
    placeholders = Set(
        _spec_key(spec) for spec in public_api_specs(spaceagora)
        if _is_placeholder_doc(Docs.doc(Docs.Binding(spec.owner_module, spec.symbol)))
    )
    return sort!(collect(setdiff(KNOWN_PLACEHOLDER_DOCS, placeholders)))
end

end # module PublicAPISymbols
