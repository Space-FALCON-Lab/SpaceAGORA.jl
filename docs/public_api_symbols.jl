module PublicAPISymbols

using Base.Docs

export PUBLIC_API_SECTIONS, public_api_specs, render_public_api_markdown, undocumented_public_api_specs

const PUBLIC_API_SECTIONS = [
    (
        title = "Simulation",
        items = [
            (owner = :SpaceAGORA, symbol = :run_simulation, rendered = "SpaceAGORA.run_simulation"),
            (owner = :SpaceAGORA, symbol = :StateAnchor, rendered = "SpaceAGORA.StateAnchor"),
            (owner = :SpaceAGORA, symbol = :get_state_anchor_callback, rendered = "SpaceAGORA.get_state_anchor_callback"),
            (owner = :SpaceAGORA, symbol = :prewarm_nbody_ephemeris_cache, rendered = "SpaceAGORA.prewarm_nbody_ephemeris_cache"),
            (owner = :SpaceAGORA, symbol = :load_nbody_ephemeris_cache!, rendered = "SpaceAGORA.load_nbody_ephemeris_cache!"),
            (owner = :SimulationCampaigns, symbol = :MonteCarloSpec, rendered = "SpaceAGORA.MonteCarloSpec"),
            (owner = :SimulationCampaigns, symbol = :MonteCarloSampleResult, rendered = "SpaceAGORA.MonteCarloSampleResult"),
            (owner = :SimulationCampaigns, symbol = :MonteCarloResult, rendered = "SpaceAGORA.MonteCarloResult"),
            (owner = :SimulationCampaigns, symbol = :run_monte_carlo, rendered = "SpaceAGORA.run_monte_carlo"),
            (owner = :SimulationCampaigns, symbol = :run_constellation_ensemble, rendered = "SpaceAGORA.run_constellation_ensemble"),
            (owner = :SimulationCampaigns, symbol = :campaign_route_features, rendered = "SpaceAGORA.campaign_route_features"),
            (owner = :SimulationCampaigns, symbol = :campaign_outer_route_state, rendered = "SpaceAGORA.campaign_outer_route_state")
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
            (owner = :SpaceAGORA, symbol = :calcControlMassFlowRate, rendered = "SpaceAGORA.calcControlMassFlowRate"),
            (owner = :SpaceAGORA, symbol = :control_thruster_levels, rendered = "SpaceAGORA.control_thruster_levels")
        ]
    ),
    (
        title = "Robotics and Planning",
        items = [
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
        title = "Compliant Multibody and Cloth Coupling",
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
        title = "Terrain and Powered Descent",
        items = [
            (owner = :SpaceAGORA, symbol = :NoTerrainModel, rendered = "SpaceAGORA.NoTerrainModel"),
            (owner = :SpaceAGORA, symbol = :DEMGrid, rendered = "SpaceAGORA.DEMGrid"),
            (owner = :SpaceAGORA, symbol = :DEMTerrainModel, rendered = "SpaceAGORA.DEMTerrainModel"),
            (owner = :SpaceAGORA, symbol = :terrain_height, rendered = "SpaceAGORA.terrain_height"),
            (owner = :SpaceAGORA, symbol = :terrain_radius, rendered = "SpaceAGORA.terrain_radius"),
            (owner = :SpaceAGORA, symbol = :load_dem_grid, rendered = "SpaceAGORA.load_dem_grid"),
            (owner = :SpaceAGORA, symbol = :load_site_terrain, rendered = "SpaceAGORA.load_site_terrain"),
            (owner = :SpaceAGORA, symbol = :DescentPhaseTargets, rendered = "SpaceAGORA.DescentPhaseTargets"),
            (owner = :SpaceAGORA, symbol = :ApolloDescentConfig, rendered = "SpaceAGORA.ApolloDescentConfig"),
            (owner = :SpaceAGORA, symbol = :ApolloDescentState, rendered = "SpaceAGORA.ApolloDescentState"),
            (owner = :SpaceAGORA, symbol = :ApolloDescentGuidanceModel, rendered = "SpaceAGORA.ApolloDescentGuidanceModel"),
            (owner = :SpaceAGORA, symbol = :apollo11_descent_targets, rendered = "SpaceAGORA.apollo11_descent_targets"),
            (owner = :SpaceAGORA, symbol = :ApolloDescentControlConfig, rendered = "SpaceAGORA.ApolloDescentControlConfig"),
            (owner = :SpaceAGORA, symbol = :ApolloDescentControlModel, rendered = "SpaceAGORA.ApolloDescentControlModel"),
            (owner = :SpaceAGORA, symbol = :descent_attitude_command, rendered = "SpaceAGORA.descent_attitude_command"),
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
        title = "Regolith Erosion Regimes",
        items = [
            (owner = :SpaceAGORA, symbol = :RegolithProperties, rendered = "SpaceAGORA.RegolithProperties"),
            (owner = :SpaceAGORA, symbol = :lunar_mare_regolith, rendered = "SpaceAGORA.lunar_mare_regolith"),
            (owner = :SpaceAGORA, symbol = :ErosionEnvironment, rendered = "SpaceAGORA.ErosionEnvironment"),
            (owner = :SpaceAGORA, symbol = :erosion_environment, rendered = "SpaceAGORA.erosion_environment"),
            (owner = :SpaceAGORA, symbol = :ErosionRegimeKind, rendered = "SpaceAGORA.ErosionRegimeKind"),
            (owner = :SpaceAGORA, symbol = :NoErosion, rendered = "SpaceAGORA.NoErosion"),
            (owner = :SpaceAGORA, symbol = :ViscousErosion, rendered = "SpaceAGORA.ViscousErosion"),
            (owner = :SpaceAGORA, symbol = :DiffusionDrivenFlowRegime, rendered = "SpaceAGORA.DiffusionDrivenFlowRegime"),
            (owner = :SpaceAGORA, symbol = :BearingCapacityFailureRegime, rendered = "SpaceAGORA.BearingCapacityFailureRegime"),
            (owner = :SpaceAGORA, symbol = :AbstractErosionRegime, rendered = "SpaceAGORA.AbstractErosionRegime"),
            (owner = :SpaceAGORA, symbol = :ViscousErosionRoberts, rendered = "SpaceAGORA.ViscousErosionRoberts"),
            (owner = :SpaceAGORA, symbol = :ViscousErosionEnergyFlux, rendered = "SpaceAGORA.ViscousErosionEnergyFlux"),
            (owner = :SpaceAGORA, symbol = :DiffusionDrivenFlow, rendered = "SpaceAGORA.DiffusionDrivenFlow"),
            (owner = :SpaceAGORA, symbol = :BearingCapacityFailure, rendered = "SpaceAGORA.BearingCapacityFailure"),
            (owner = :SpaceAGORA, symbol = :erosion_rate, rendered = "SpaceAGORA.erosion_rate"),
            (owner = :SpaceAGORA, symbol = :erosion_onset, rendered = "SpaceAGORA.erosion_onset"),
            (owner = :SpaceAGORA, symbol = :regime_kind, rendered = "SpaceAGORA.regime_kind"),
            (owner = :SpaceAGORA, symbol = :default_erosion_regimes, rendered = "SpaceAGORA.default_erosion_regimes"),
            (owner = :SpaceAGORA, symbol = :regolith_erosion_rate, rendered = "SpaceAGORA.regolith_erosion_rate"),
            (owner = :SpaceAGORA, symbol = :shields_threshold_shear_pa, rendered = "SpaceAGORA.shields_threshold_shear_pa"),
            (owner = :SpaceAGORA, symbol = :energy_flux_threshold_shear_pa, rendered = "SpaceAGORA.energy_flux_threshold_shear_pa"),
            (owner = :SpaceAGORA, symbol = :soil_bearing_capacity_pa, rendered = "SpaceAGORA.soil_bearing_capacity_pa"),
            (owner = :SpaceAGORA, symbol = :mean_thermal_speed_mps, rendered = "SpaceAGORA.mean_thermal_speed_mps"),
            (owner = :SpaceAGORA, symbol = :mean_lift_height_m, rendered = "SpaceAGORA.mean_lift_height_m"),
            (owner = :SpaceAGORA, symbol = :gas_dynamic_viscosity_pa_s, rendered = "SpaceAGORA.gas_dynamic_viscosity_pa_s"),
            (owner = :SpaceAGORA, symbol = :pressure_diffusion_depth_m, rendered = "SpaceAGORA.pressure_diffusion_depth_m"),
            (owner = :SpaceAGORA, symbol = :soil_tensile_strength_pa, rendered = "SpaceAGORA.soil_tensile_strength_pa")
        ]
    ),
    (
        title = "Mesh Aerodynamics",
        items = [
            (owner = :SpaceAGORA, symbol = :MeshAeroPanels, rendered = "SpaceAGORA.MeshAeroPanels"),
            (owner = :SpaceAGORA, symbol = :MeshAeroSurrogate, rendered = "SpaceAGORA.MeshAeroSurrogate"),
            (owner = :SpaceAGORA, symbol = :AerodynamicCoefficientMeshSurrogate, rendered = "SpaceAGORA.AerodynamicCoefficientMeshSurrogate"),
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
            (owner = :ParallelProfiles, symbol = :record_outer_route_feedback!, rendered = "SpaceAGORA.record_outer_route_feedback!"),
            (owner = :ParallelProcess, symbol = :ProcessPool, rendered = "SpaceAGORA.ProcessPool"),
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
