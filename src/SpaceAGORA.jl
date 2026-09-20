__precompile__(true)

module SpaceAGORA

## 1. Include package modules
include(joinpath(@__DIR__, "parallel", "routing", "parallel_profiles.jl"))
include(joinpath(@__DIR__, "parallel", "process", "parallel_process.jl"))
include(joinpath(@__DIR__, "simulation", "runtime_services.jl"))
include(joinpath(@__DIR__, "core", "simulation_model.jl"))
include(joinpath(@__DIR__, "simulation", "engine", "simulation_engine.jl"))
include(joinpath(@__DIR__, "simulation", "campaigns", "simulation_campaigns.jl"))
include(joinpath(@__DIR__, "analysis", "verification", "telemetry_verification.jl"))
include(joinpath(@__DIR__, "assets", "rpo_station_assets.jl"))
include(joinpath(@__DIR__, "assets", "odyssey_surrogate_assets.jl"))
include(joinpath(@__DIR__, "analysis", "visualization", "rpo", "rpo_visualization.jl"))
include(joinpath(@__DIR__, "cli", "spaceagora_cli.jl"))


## 2. Bring needed names from package modules into the scope of SpaceAGORA.jl
# 2.1. Parallel Profiles
using .ParallelProfiles: ParallelProfile, ParallelProfileConfig
using .ParallelProfiles: parse_parallel_profile, parallel_profile_name, profile_config, profile_env_pairs, with_parallel_profile
using .ParallelProfiles: OuterRouteFeatures, OuterRouteTuning, OuterRouteState
using .ParallelProfiles: reset_outer_route_state!, outer_route_signature, outer_route_stats_snapshot
using .ParallelProfiles: default_outer_route, outer_route_candidates, select_outer_route!, record_outer_route_feedback!

# 2.2. Parallel Process
using .ParallelProcess: ProcessPool, campaign_process_pool, ensure_process_workers!, shutdown_process_pool!, adopt_process_workers!

## 2.3. Simulation Engine
using .SimulationEngine: ParallelConfig, SolverConfig, RuntimePolicyConfig, ArtifactConfig, SimulationEngineConfig
using .SimulationEngine: simulation_engine_config_from_env
using .SimulationEngine: prewarm_nbody_ephemeris_cache, load_nbody_ephemeris_cache!
run_simulation(args...; kwargs...) = SimulationEngine.run_simulation(args...; kwargs...)

## 2.4. Simulation Campaigns
using .SimulationCampaigns: MonteCarloSpec, MonteCarloSampleResult, MonteCarloResult, run_monte_carlo
using .SimulationCampaigns: run_constellation_ensemble
using .SimulationCampaigns: run_monte_carlo_visualization
using .SimulationCampaigns: campaign_route_features, campaign_outer_route_state

## 2.5. Simulation Model
using .SimulationModel: StateAnchor, get_state_anchor_callback
using .SimulationModel.AbstractTypes: AbstractForceTorqueModel, AbstractPlanet, AbstractDensityModel
using .SimulationModel.AbstractTypes: AbstractControlEffectorModel, AbstractEphemeridesModel
using .SimulationModel.AbstractTypes: AbstractThermalModel, AbstractThrusterModel, AbstractGuidanceModel
using .SimulationModel: StateSample, PlanetFrameSample, AtmosphereSample, SolarEphemerisSample
using .SimulationModel: ThirdBodyEphemerisSample, EnvironmentSample, EffectorEnvironmentRequirements
using .SimulationModel: ClothArmBasePose, ClothArmLink, ClothArmJoint, ClothArmModel, ClothArmPose, ClothArmState
using .SimulationModel: default_cloth_arm_model, cloth_fk, cloth_fk_state, cloth_end_effector_pose
using .SimulationModel: cloth_ik, cloth_total_reach, closest_surface_target
using .SimulationModel: RobotArmPlannerConfig, RobotArmPlan, plan_robot_arm_motion, robot_arm_plan_sample
using .SimulationModel: RobotArmSphereObstacle, RobotArmHYPRConfig, RobotArmHYPRResult
using .SimulationModel: plan_robot_arm_motion_hypr, robot_arm_sample_hypr_path
using .SimulationModel: robot_arm_clearance_stats_from_samples, robot_arm_hypr_path_cost_components
using .SimulationModel: ClothRobotArmReferenceState, cloth_reference_state
using .SimulationModel: CompliantBody, CompliantJoint, CompliantMultibodyModel, CompliantMultibodyTrajectory
using .SimulationModel: CompliantTopologyNode, CompliantTopologyEdge, CompliantTopologyBuild
using .SimulationModel: CompliantJointActuator, CompliantJointLoad
using .SimulationModel: rectangular_prism_inertia, thin_panel_inertia
using .SimulationModel: build_compliant_topology, build_rectangular_compliant_grid
using .SimulationModel: compliant_state_vector, compliant_state_parts, compliant_multibody_dynamics
using .SimulationModel: compliant_joint_loads
using .SimulationModel: step_compliant_multibody_rk4, step_compliant_multibody_implicit_midpoint
using .SimulationModel: simulate_compliant_multibody
using .SimulationModel: ClothRobotArmSimulation, cloth_robot_arm_multibody, cloth_robot_arm_initial_state
using .SimulationModel: cloth_robot_arm_rest_quaternions, cloth_robot_arm_end_effector
using .SimulationModel: cloth_robot_arm_actuators
using .SimulationModel: coupled_cloth_robot_arm_state_shape, initialize_coupled_cloth_robot_arm_state!
using .SimulationModel: assign_coupled_cloth_robot_arm_rhs!
using .SimulationModel: simulate_cloth_robot_arm_plan
using .SimulationModel: RobotArmHeldActuation, RobotArmJointMPCController, RobotArmControlEffector, RobotArmReactionEffector
using .SimulationModel: init_robot_arm_joint_mpc, robot_arm_joint_mpc_reference_preview
using .SimulationModel: robot_arm_joint_mpc_control, robot_arm_measured_joint_state
using .SimulationModel: NoAtmosphereModel, ExponentialAtmosphereModel, PiecewiseExponentialAtmosphereModel
using .SimulationModel: GRAMGridAtmosphereModel
using .SimulationModel.EnvironmentModels: SurrogatePresetResolution, available_surrogate_presets, resolve_surrogate_preset, surrogate_preset_model, atmosphere_provenance
using .SimulationModel: NRLMSISE00AtmosphereModel, init_nrlmsise_space_indices!
using .SimulationModel: SimpleEphemeridesModel
using .SimulationModel.TerrainModels: AbstractTerrainModel, NoTerrainModel, DEMGrid, DEMTerrainModel
using .SimulationModel.TerrainModels: terrain_height, terrain_radius, load_dem_grid, load_site_terrain, dem_grid_covers
using .SimulationModel: make_no_gram_planet, make_no_gram_density_model, make_no_gram_environment
using .SimulationModel: calcForceTorque, wrench, environment_requirements, solver_partition
using .SimulationModel: gravity_backbone_structure, gravity_backbone_acceleration_ii
using .SimulationModel: gravity_backbone_kick_structure, gravity_backbone_kick_acceleration_ii
using .SimulationModel: PlumeSurfaceConfig, PlumeSurfaceInteractionModel, PlumeSurfaceState, plume_surface_footprint, plume_erosion_onset_height, plume_ground_effect_force, plume_quantities
using .SimulationModel: getDensity, getDensityBatch!, with_density_model_epoch
using .SimulationModel: DescentPhaseTargets, apollo11_descent_targets, ApolloDescentConfig, ApolloDescentState, ApolloDescentGuidanceModel
using .SimulationModel: ApolloDescentControlConfig, ApolloDescentControlModel, descent_attitude_command
using .SimulationModel: calcControlEffect!, calcControlForceTorque, calcControlMassFlowRate
using .SimulationModel: control_thruster_levels
using .SimulationModel: AerobrakingEnergyDepletionConfig, AerobrakingEnergyDepletionState
using .SimulationModel: AerobrakingEnergyDepletionGuidanceModel, AerobrakingEnergyDepletionControlModel
using .SimulationModel: SolarPanelAngleOfAttackControlModel
using .SimulationModel: ApoapsisTargetPeriapsisRaiseGuidanceModel
using .SimulationModel: VisualizationScene, PlanetSpec, SpacecraftGeometry, LinkBox, AtmosphereSpec, atmosphere_spec
using .SimulationModel: ArmGeometry, arm_geometry
using .SimulationModel: load_model_triangles, model_bounding_box, sample_model_pointcloud, articulate_triangles, articulation_payload
using .SimulationModel: spacecraft_geometry, planet_spec, planet_rotation_table, build_visualization_scene
using .SimulationModel: visualization_scene_path, write_visualization_scene, read_visualization_scene
using .SimulationModel: velocity_aligned_quaternion, visualization_frame_budget
using .SimulationModel: MeshAeroPanels, MeshAeroSurrogate, AerodynamicCoefficientMeshSurrogate, MESH_AERO_MAX_DEGREE
using .SimulationModel: mesh_aero_panels, panel_aero_coefficients, panel_aero_coefficients_split, panel_shadow_mask, panel_projected_area
using .SimulationModel: fit_mesh_aero_surrogate, mesh_aero_coefficients, write_mesh_aero_surrogate, read_mesh_aero_surrogate
using .SimulationModel: export_visualization, with_visualization_scene, write_viewer_dev_payload
using .SimulationModel: EnsembleSample, sample_results_directory, with_results_directory, write_ensemble_manifest, export_ensemble_visualization

## 2.6. Telemetry Verification
using .TelemetryVerification: VerificationRequest, VerificationResult
using .TelemetryVerification: run_verification, run_verification_cli, run_study

using .OdysseySurrogateAssets: odyssey_surrogate_assets
export odyssey_surrogate_assets

## 2.7. RPO Station Assets
using .RPOStationAssets: station_geometry_path, station_cad_path, load_rpo_station_pointcloud, load_rpo_station_cad_triangles, load_rpo_station_cad_pointcloud

## 2.8. RPO Visualization
using .RPOVisualization: rpo_path_plot, rpo_tracking_plot

## 2.9. SpaceAGORA CLI
using .SpaceAGORACLI: AssetCheckItem, AssetCheckReport
using .SpaceAGORACLI: check_assets, render_asset_report, run_cli


@doc (@doc SimulationModel.EnvironmentModels.SurrogatePresetResolution) SurrogatePresetResolution
@doc (@doc SimulationModel.EnvironmentModels.available_surrogate_presets) available_surrogate_presets
@doc (@doc SimulationModel.EnvironmentModels.resolve_surrogate_preset) resolve_surrogate_preset
@doc (@doc SimulationModel.EnvironmentModels.surrogate_preset_model) surrogate_preset_model
@doc (@doc SimulationModel.EnvironmentModels.atmosphere_provenance) atmosphere_provenance
@doc (@doc OdysseySurrogateAssets.odyssey_surrogate_assets) odyssey_surrogate_assets

## 3. Attach description of each function / model to the SpaceAGORA module's bindings
# Forward the docstrings onto this module's bindings: the docs build resolves
# `@docs SpaceAGORA.X` blocks against SpaceAGORA's own doc metadata, and the
# CI environment does not follow the explicit-import alias for these.
# 3.1. Simulation Engine
@doc (@doc SimulationEngine.ParallelConfig) ParallelConfig
@doc (@doc SimulationEngine.SolverConfig) SolverConfig
@doc (@doc SimulationEngine.RuntimePolicyConfig) RuntimePolicyConfig
@doc (@doc SimulationEngine.ArtifactConfig) ArtifactConfig
@doc (@doc SimulationEngine.SimulationEngineConfig) SimulationEngineConfig
@doc (@doc SimulationEngine.simulation_engine_config_from_env) simulation_engine_config_from_env
@doc (@doc SimulationEngine.run_simulation) run_simulation
@doc (@doc SimulationEngine.prewarm_nbody_ephemeris_cache) prewarm_nbody_ephemeris_cache
@doc (@doc SimulationEngine.load_nbody_ephemeris_cache!) load_nbody_ephemeris_cache!

# 3.2. Simulation Campaigns
@doc (@doc SimulationCampaigns.MonteCarloSpec) MonteCarloSpec
@doc (@doc SimulationCampaigns.MonteCarloSampleResult) MonteCarloSampleResult
@doc (@doc SimulationCampaigns.MonteCarloResult) MonteCarloResult
@doc (@doc SimulationCampaigns.run_monte_carlo) run_monte_carlo
@doc (@doc SimulationCampaigns.run_constellation_ensemble) run_constellation_ensemble
@doc (@doc SimulationCampaigns.campaign_route_features) campaign_route_features
@doc (@doc SimulationCampaigns.campaign_outer_route_state) campaign_outer_route_state

# 3.3. Simulation Model
@doc (@doc SimulationModel.AerobrakingEnergyDepletionConfig) AerobrakingEnergyDepletionConfig
@doc (@doc SimulationModel.AerobrakingEnergyDepletionState) AerobrakingEnergyDepletionState
@doc (@doc SimulationModel.AerobrakingEnergyDepletionGuidanceModel) AerobrakingEnergyDepletionGuidanceModel
@doc (@doc SimulationModel.AerobrakingEnergyDepletionControlModel) AerobrakingEnergyDepletionControlModel
@doc (@doc SimulationModel.SolarPanelAngleOfAttackControlModel) SolarPanelAngleOfAttackControlModel
@doc (@doc SimulationModel.StateAnchor) StateAnchor
@doc (@doc SimulationModel.get_state_anchor_callback) get_state_anchor_callback
@doc (@doc SimulationModel.AbstractTypes.AbstractForceTorqueModel) AbstractForceTorqueModel
@doc (@doc SimulationModel.AbstractTypes.AbstractPlanet) AbstractPlanet
@doc (@doc SimulationModel.AbstractTypes.AbstractDensityModel) AbstractDensityModel
@doc (@doc SimulationModel.AbstractTypes.AbstractControlEffectorModel) AbstractControlEffectorModel
@doc (@doc SimulationModel.AbstractTypes.AbstractEphemeridesModel) AbstractEphemeridesModel
@doc (@doc SimulationModel.AbstractTypes.AbstractThermalModel) AbstractThermalModel
@doc (@doc SimulationModel.AbstractTypes.AbstractThrusterModel) AbstractThrusterModel
@doc (@doc SimulationModel.AbstractTypes.AbstractGuidanceModel) AbstractGuidanceModel
@doc (@doc SimulationModel.StateSample) StateSample
@doc (@doc SimulationModel.PlanetFrameSample) PlanetFrameSample
@doc (@doc SimulationModel.AtmosphereSample) AtmosphereSample
@doc (@doc SimulationModel.SolarEphemerisSample) SolarEphemerisSample
@doc (@doc SimulationModel.ThirdBodyEphemerisSample) ThirdBodyEphemerisSample
@doc (@doc SimulationModel.EnvironmentSample) EnvironmentSample
@doc (@doc SimulationModel.EffectorEnvironmentRequirements) EffectorEnvironmentRequirements
@doc (@doc SimulationModel.ClothArmModel) ClothArmModel
@doc (@doc SimulationModel.ClothArmBasePose) ClothArmBasePose
@doc (@doc SimulationModel.ClothArmLink) ClothArmLink
@doc (@doc SimulationModel.ClothArmJoint) ClothArmJoint
@doc (@doc SimulationModel.ClothArmPose) ClothArmPose
@doc (@doc SimulationModel.ClothArmState) ClothArmState
@doc (@doc SimulationModel.default_cloth_arm_model) default_cloth_arm_model
@doc (@doc SimulationModel.cloth_fk) cloth_fk
@doc (@doc SimulationModel.cloth_fk_state) cloth_fk_state
@doc (@doc SimulationModel.cloth_end_effector_pose) cloth_end_effector_pose
@doc (@doc SimulationModel.cloth_ik) cloth_ik
@doc (@doc SimulationModel.cloth_total_reach) cloth_total_reach
@doc (@doc SimulationModel.closest_surface_target) closest_surface_target
@doc (@doc SimulationModel.RobotArmPlannerConfig) RobotArmPlannerConfig
@doc (@doc SimulationModel.RobotArmPlan) RobotArmPlan
@doc (@doc SimulationModel.robot_arm_plan_sample) robot_arm_plan_sample
@doc (@doc SimulationModel.RobotArmSphereObstacle) RobotArmSphereObstacle
@doc (@doc SimulationModel.RobotArmHYPRConfig) RobotArmHYPRConfig
@doc (@doc SimulationModel.RobotArmHYPRResult) RobotArmHYPRResult
@doc (@doc SimulationModel.plan_robot_arm_motion) plan_robot_arm_motion
@doc (@doc SimulationModel.plan_robot_arm_motion_hypr) plan_robot_arm_motion_hypr
@doc (@doc SimulationModel.robot_arm_sample_hypr_path) robot_arm_sample_hypr_path
@doc (@doc SimulationModel.robot_arm_clearance_stats_from_samples) robot_arm_clearance_stats_from_samples
@doc (@doc SimulationModel.robot_arm_hypr_path_cost_components) robot_arm_hypr_path_cost_components
@doc (@doc SimulationModel.CompliantBody) CompliantBody
@doc (@doc SimulationModel.CompliantJoint) CompliantJoint
@doc (@doc SimulationModel.CompliantMultibodyModel) CompliantMultibodyModel
@doc (@doc SimulationModel.CompliantMultibodyTrajectory) CompliantMultibodyTrajectory
@doc (@doc SimulationModel.CompliantTopologyNode) CompliantTopologyNode
@doc (@doc SimulationModel.CompliantTopologyEdge) CompliantTopologyEdge
@doc (@doc SimulationModel.CompliantTopologyBuild) CompliantTopologyBuild
@doc (@doc SimulationModel.CompliantJointActuator) CompliantJointActuator
@doc (@doc SimulationModel.CompliantJointLoad) CompliantJointLoad
@doc (@doc SimulationModel.rectangular_prism_inertia) rectangular_prism_inertia
@doc (@doc SimulationModel.thin_panel_inertia) thin_panel_inertia
@doc (@doc SimulationModel.build_compliant_topology) build_compliant_topology
@doc (@doc SimulationModel.build_rectangular_compliant_grid) build_rectangular_compliant_grid
@doc (@doc SimulationModel.compliant_state_vector) compliant_state_vector
@doc (@doc SimulationModel.compliant_state_parts) compliant_state_parts
@doc (@doc SimulationModel.compliant_multibody_dynamics) compliant_multibody_dynamics
@doc (@doc SimulationModel.compliant_joint_loads) compliant_joint_loads
@doc (@doc SimulationModel.step_compliant_multibody_rk4) step_compliant_multibody_rk4
@doc (@doc SimulationModel.step_compliant_multibody_implicit_midpoint) step_compliant_multibody_implicit_midpoint
@doc (@doc SimulationModel.simulate_compliant_multibody) simulate_compliant_multibody
@doc (@doc SimulationModel.ClothRobotArmSimulation) ClothRobotArmSimulation
@doc (@doc SimulationModel.cloth_robot_arm_multibody) cloth_robot_arm_multibody
@doc (@doc SimulationModel.cloth_robot_arm_initial_state) cloth_robot_arm_initial_state
@doc (@doc SimulationModel.cloth_robot_arm_rest_quaternions) cloth_robot_arm_rest_quaternions
@doc (@doc SimulationModel.cloth_robot_arm_end_effector) cloth_robot_arm_end_effector
@doc (@doc SimulationModel.simulate_cloth_robot_arm_plan) simulate_cloth_robot_arm_plan
@doc (@doc SimulationModel.cloth_robot_arm_actuators) cloth_robot_arm_actuators
@doc (@doc SimulationModel.coupled_cloth_robot_arm_state_shape) coupled_cloth_robot_arm_state_shape
@doc (@doc SimulationModel.initialize_coupled_cloth_robot_arm_state!) initialize_coupled_cloth_robot_arm_state!
@doc (@doc SimulationModel.assign_coupled_cloth_robot_arm_rhs!) assign_coupled_cloth_robot_arm_rhs!
@doc (@doc SimulationModel.RobotArmHeldActuation) RobotArmHeldActuation
@doc (@doc SimulationModel.RobotArmJointMPCController) RobotArmJointMPCController
@doc (@doc SimulationModel.RobotArmControlEffector) RobotArmControlEffector
@doc (@doc SimulationModel.RobotArmReactionEffector) RobotArmReactionEffector
@doc (@doc SimulationModel.init_robot_arm_joint_mpc) init_robot_arm_joint_mpc
@doc (@doc SimulationModel.robot_arm_joint_mpc_reference_preview) robot_arm_joint_mpc_reference_preview
@doc (@doc SimulationModel.robot_arm_joint_mpc_control) robot_arm_joint_mpc_control
@doc (@doc SimulationModel.robot_arm_measured_joint_state) robot_arm_measured_joint_state
@doc (@doc SimulationModel.ApoapsisTargetPeriapsisRaiseGuidanceModel) ApoapsisTargetPeriapsisRaiseGuidanceModel
@doc (@doc SimulationModel.AbstractTypes.AbstractTerrainModel) AbstractTerrainModel
@doc (@doc SimulationModel.TerrainModels.NoTerrainModel) NoTerrainModel
@doc (@doc SimulationModel.TerrainModels.DEMGrid) DEMGrid
@doc (@doc SimulationModel.TerrainModels.DEMTerrainModel) DEMTerrainModel
@doc (@doc SimulationModel.TerrainModels.terrain_height) terrain_height
@doc (@doc SimulationModel.TerrainModels.terrain_radius) terrain_radius
@doc (@doc SimulationModel.TerrainModels.load_dem_grid) load_dem_grid
@doc (@doc SimulationModel.TerrainModels.load_site_terrain) load_site_terrain
@doc (@doc SimulationModel.TerrainModels.dem_grid_covers) dem_grid_covers
@doc (@doc SimulationModel.NoAtmosphereModel) NoAtmosphereModel
@doc (@doc SimulationModel.ExponentialAtmosphereModel) ExponentialAtmosphereModel
@doc (@doc SimulationModel.PiecewiseExponentialAtmosphereModel) PiecewiseExponentialAtmosphereModel
@doc (@doc SimulationModel.GRAMGridAtmosphereModel) GRAMGridAtmosphereModel
@doc (@doc SimulationModel.NRLMSISE00AtmosphereModel) NRLMSISE00AtmosphereModel
@doc (@doc SimulationModel.DescentPhaseTargets) DescentPhaseTargets
@doc (@doc SimulationModel.apollo11_descent_targets) apollo11_descent_targets
@doc (@doc SimulationModel.ApolloDescentConfig) ApolloDescentConfig
@doc (@doc SimulationModel.ApolloDescentState) ApolloDescentState
@doc (@doc SimulationModel.ApolloDescentGuidanceModel) ApolloDescentGuidanceModel
@doc (@doc SimulationModel.ApolloDescentControlConfig) ApolloDescentControlConfig
@doc (@doc SimulationModel.ApolloDescentControlModel) ApolloDescentControlModel
@doc (@doc SimulationModel.descent_attitude_command) descent_attitude_command
@doc (@doc SimulationModel.init_nrlmsise_space_indices!) init_nrlmsise_space_indices!
@doc (@doc SimulationModel.SimpleEphemeridesModel) SimpleEphemeridesModel
@doc (@doc SimulationModel.make_no_gram_planet) make_no_gram_planet
@doc (@doc SimulationModel.make_no_gram_density_model) make_no_gram_density_model
@doc (@doc SimulationModel.make_no_gram_environment) make_no_gram_environment
# Copy canonical hook text into fresh documentation metadata for the root binding.
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.DynamicEffectors.calcForceTorque)).text) calcForceTorque
@doc (@doc SimulationModel.wrench) wrench
@doc (@doc SimulationModel.environment_requirements) environment_requirements
@doc (@doc SimulationModel.solver_partition) solver_partition
@doc (@doc SimulationModel.gravity_backbone_structure) gravity_backbone_structure
@doc (@doc SimulationModel.gravity_backbone_acceleration_ii) gravity_backbone_acceleration_ii
@doc (@doc SimulationModel.gravity_backbone_kick_structure) gravity_backbone_kick_structure
@doc (@doc SimulationModel.gravity_backbone_kick_acceleration_ii) gravity_backbone_kick_acceleration_ii
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.DynamicEffectors.PlumeSurfaceInteraction.PlumeSurfaceConfig)).text) PlumeSurfaceConfig
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.DynamicEffectors.PlumeSurfaceInteraction.PlumeSurfaceInteractionModel)).text) PlumeSurfaceInteractionModel
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.DynamicEffectors.PlumeSurfaceInteraction.PlumeSurfaceState)).text) PlumeSurfaceState
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.DynamicEffectors.PlumeSurfaceInteraction.plume_surface_footprint(::SimulationModel.DynamicEffectors.PlumeSurfaceInteraction.PlumeSurfaceConfig, ::Real, ::Real))).text) plume_surface_footprint
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.DynamicEffectors.PlumeSurfaceInteraction.plume_erosion_onset_height(::SimulationModel.DynamicEffectors.PlumeSurfaceInteraction.PlumeSurfaceConfig, ::Real))).text) plume_erosion_onset_height
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.DynamicEffectors.PlumeSurfaceInteraction.plume_ground_effect_force(::SimulationModel.DynamicEffectors.PlumeSurfaceInteraction.PlumeSurfaceConfig, ::Real, ::Real))).text) plume_ground_effect_force
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.DynamicEffectors.PlumeSurfaceInteraction.plume_quantities(::SimulationModel.DynamicEffectors.PlumeSurfaceInteraction.PlumeSurfaceConfig, ::Real, ::Real))).text) plume_quantities
@doc (@doc SimulationModel.getDensity) getDensity
@doc (@doc SimulationModel.getDensityBatch!) getDensityBatch!
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.EnvironmentModels.with_density_model_epoch(
    ::SimulationModel.AbstractDensityModel, ::Any,
))).text) with_density_model_epoch
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.ControlHooks.calcControlEffect!(
    ::SimulationModel.ControlHooks.BaseThrusterModel,
    ::SimulationModel.ControlHooks.ComponentVector,
    ::SimulationModel.ODEParams, ::Float64, ::Int64,
))).text) calcControlEffect!
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.ControlHooks.calcControlForceTorque(
    ::SimulationModel.ControlHooks.BaseThrusterModel,
    ::AbstractVector, ::SimulationModel.ODEParams, ::Int64, ::Float64,
))).text) calcControlForceTorque
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.ControlHooks.calcControlMassFlowRate(
    ::SimulationModel.AbstractTypes.AbstractControlEffectorModel,
    ::AbstractVector, ::SimulationModel.ODEParams, ::Int64, ::Float64,
))).text) calcControlMassFlowRate
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.ControlHooks.control_thruster_levels(::Any, ::Int))).text) control_thruster_levels

# 3.4. RPO Station Assets
@doc (@doc RPOStationAssets.station_geometry_path) station_geometry_path
@doc (@doc RPOStationAssets.station_cad_path) station_cad_path
@doc (@doc RPOStationAssets.load_rpo_station_pointcloud) load_rpo_station_pointcloud
@doc (@doc RPOStationAssets.load_rpo_station_cad_triangles) load_rpo_station_cad_triangles
@doc (@doc RPOStationAssets.load_rpo_station_cad_pointcloud) load_rpo_station_cad_pointcloud

# Scene visualization and display geometry
@doc (@doc SimulationModel.SceneVisualization.VisualizationScene) VisualizationScene
@doc (@doc SimulationModel.SceneVisualization.PlanetSpec) PlanetSpec
@doc (@doc SimulationModel.SceneVisualization.SpacecraftGeometry) SpacecraftGeometry
@doc (@doc SimulationModel.SceneVisualization.LinkBox) LinkBox
@doc (@doc SimulationModel.SceneVisualization.AtmosphereSpec) AtmosphereSpec
@doc (@doc SimulationModel.SceneVisualization.atmosphere_spec) atmosphere_spec
@doc (@doc SimulationModel.SceneVisualization.ArmGeometry) ArmGeometry
@doc (@doc SimulationModel.SceneVisualization.arm_geometry) arm_geometry
@doc (@doc SimulationModel.Structure.load_model_triangles) load_model_triangles
@doc (@doc SimulationModel.Structure.model_bounding_box) model_bounding_box
@doc (@doc SimulationModel.Structure.sample_model_pointcloud) sample_model_pointcloud
@doc (@doc SimulationModel.Structure.articulate_triangles) articulate_triangles
@doc (@doc SimulationModel.Structure.articulation_payload) articulation_payload
@doc (@doc SimulationModel.SceneVisualization.spacecraft_geometry) spacecraft_geometry
@doc (@doc SimulationModel.SceneVisualization.planet_spec) planet_spec
@doc (@doc SimulationModel.SceneVisualization.planet_rotation_table) planet_rotation_table
@doc (@doc SimulationModel.SceneVisualization.build_visualization_scene) build_visualization_scene
@doc (@doc SimulationModel.SceneVisualization.visualization_scene_path) visualization_scene_path
@doc (@doc SimulationModel.SceneVisualization.write_visualization_scene) write_visualization_scene
@doc (@doc SimulationModel.SceneVisualization.read_visualization_scene) read_visualization_scene
@doc (@doc SimulationModel.SceneVisualization.velocity_aligned_quaternion) velocity_aligned_quaternion
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.MeshAeroPanels) MeshAeroPanels
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.MeshAeroSurrogate) MeshAeroSurrogate
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.AerodynamicCoefficientMeshSurrogate) AerodynamicCoefficientMeshSurrogate
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.MESH_AERO_MAX_DEGREE) MESH_AERO_MAX_DEGREE
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.mesh_aero_panels) mesh_aero_panels
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.panel_aero_coefficients) panel_aero_coefficients
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.panel_aero_coefficients_split) panel_aero_coefficients_split
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.panel_shadow_mask) panel_shadow_mask
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.panel_projected_area) panel_projected_area
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.fit_mesh_aero_surrogate) fit_mesh_aero_surrogate
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.mesh_aero_coefficients) mesh_aero_coefficients
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.write_mesh_aero_surrogate) write_mesh_aero_surrogate
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.read_mesh_aero_surrogate) read_mesh_aero_surrogate
@doc (@doc SimulationModel.SceneVisualization.visualization_frame_budget) visualization_frame_budget
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.SceneVisualization.export_visualization(::AbstractString))).text) export_visualization
@doc (@doc SimulationModel.SceneVisualization.with_visualization_scene) with_visualization_scene
@doc (@doc SimulationModel.SceneVisualization.write_viewer_dev_payload) write_viewer_dev_payload
@doc (@doc SimulationModel.SceneVisualization.EnsembleSample) EnsembleSample
@doc (@doc SimulationModel.SceneVisualization.sample_results_directory) sample_results_directory
@doc (@doc SimulationModel.SceneVisualization.with_results_directory) with_results_directory
@doc (@doc SimulationModel.SceneVisualization.write_ensemble_manifest) write_ensemble_manifest
@doc (@doc SimulationModel.SceneVisualization.export_ensemble_visualization) export_ensemble_visualization
@doc (@doc SimulationCampaigns.run_monte_carlo_visualization) run_monte_carlo_visualization

# 3.5. Parallel Profiles
@doc (@doc ParallelProfiles.ParallelProfile) ParallelProfile
@doc (@doc ParallelProfiles.ParallelProfileConfig) ParallelProfileConfig
@doc (@doc ParallelProfiles.parse_parallel_profile) parse_parallel_profile
@doc (@doc ParallelProfiles.parallel_profile_name) parallel_profile_name
@doc (@doc ParallelProfiles.profile_config) profile_config
@doc (@doc ParallelProfiles.profile_env_pairs) profile_env_pairs
@doc (@doc ParallelProfiles.with_parallel_profile) with_parallel_profile
@doc (@doc ParallelProfiles.OuterRouteFeatures) OuterRouteFeatures
@doc (@doc ParallelProfiles.OuterRouteTuning) OuterRouteTuning
@doc (@doc ParallelProfiles.OuterRouteState) OuterRouteState
@doc (@doc ParallelProfiles.reset_outer_route_state!) reset_outer_route_state!
@doc (@doc ParallelProfiles.outer_route_signature) outer_route_signature
@doc (@doc ParallelProfiles.outer_route_stats_snapshot) outer_route_stats_snapshot
@doc (@doc ParallelProfiles.default_outer_route) default_outer_route
@doc (@doc ParallelProfiles.outer_route_candidates) outer_route_candidates
@doc (@doc ParallelProfiles.select_outer_route!) select_outer_route!
@doc (@doc ParallelProfiles.record_outer_route_feedback!) record_outer_route_feedback!

# 3.6. Parallel Process
@doc (@doc ParallelProcess.ProcessPool) ProcessPool
@doc (@doc ParallelProcess.campaign_process_pool) campaign_process_pool
@doc (@doc ParallelProcess.ensure_process_workers!) ensure_process_workers!
@doc (@doc ParallelProcess.shutdown_process_pool!) shutdown_process_pool!
@doc (@doc ParallelProcess.adopt_process_workers!) adopt_process_workers!

# 3.7. Telemetry Verification
@doc (@doc TelemetryVerification.VerificationRequest) VerificationRequest
@doc (@doc TelemetryVerification.VerificationResult) VerificationResult
@doc (@doc TelemetryVerification.run_verification) run_verification
@doc (@doc TelemetryVerification.run_verification_cli) run_verification_cli
@doc (@doc TelemetryVerification.run_study) run_study

# 3.8. SpaceAGORA CLI
@doc (@doc SpaceAGORACLI.AssetCheckItem) AssetCheckItem
@doc (@doc SpaceAGORACLI.AssetCheckReport) AssetCheckReport
@doc (@doc SpaceAGORACLI.check_assets) check_assets
@doc (@doc SpaceAGORACLI.render_asset_report) render_asset_report
@doc (@doc SpaceAGORACLI.run_cli) run_cli


## 4. Declare exports
export ParallelProfile, ParallelProfileConfig
export parse_parallel_profile, parallel_profile_name, profile_config, profile_env_pairs, with_parallel_profile
export OuterRouteFeatures, OuterRouteTuning, OuterRouteState
export reset_outer_route_state!, outer_route_signature, outer_route_stats_snapshot
export default_outer_route, outer_route_candidates, select_outer_route!, record_outer_route_feedback!
export ProcessPool, campaign_process_pool, ensure_process_workers!, shutdown_process_pool!, adopt_process_workers!
export ParallelConfig, SolverConfig, RuntimePolicyConfig, ArtifactConfig, SimulationEngineConfig
export simulation_engine_config_from_env
export prewarm_nbody_ephemeris_cache, load_nbody_ephemeris_cache!
export MonteCarloSpec, MonteCarloSampleResult, MonteCarloResult, run_monte_carlo
export run_constellation_ensemble
export campaign_route_features, campaign_outer_route_state
export StateAnchor, get_state_anchor_callback
export AbstractForceTorqueModel, AbstractPlanet, AbstractDensityModel, AbstractControlEffectorModel
export AbstractEphemeridesModel, AbstractThermalModel, AbstractThrusterModel, AbstractGuidanceModel
export StateSample, PlanetFrameSample, AtmosphereSample, SolarEphemerisSample
export ThirdBodyEphemerisSample, EnvironmentSample, EffectorEnvironmentRequirements
export ClothArmBasePose, ClothArmLink, ClothArmJoint, ClothArmModel, ClothArmPose, ClothArmState
export default_cloth_arm_model, cloth_fk, cloth_fk_state, cloth_end_effector_pose
export cloth_ik, cloth_total_reach, closest_surface_target
export RobotArmPlannerConfig, RobotArmPlan, plan_robot_arm_motion, robot_arm_plan_sample
export RobotArmSphereObstacle, RobotArmHYPRConfig, RobotArmHYPRResult
export plan_robot_arm_motion_hypr, robot_arm_sample_hypr_path
export robot_arm_clearance_stats_from_samples, robot_arm_hypr_path_cost_components
export ClothRobotArmSimulation, cloth_robot_arm_multibody, cloth_robot_arm_initial_state
export cloth_robot_arm_rest_quaternions, cloth_robot_arm_end_effector, simulate_cloth_robot_arm_plan
export cloth_robot_arm_actuators
export CompliantBody, CompliantJoint, CompliantMultibodyModel, CompliantMultibodyTrajectory
export CompliantTopologyNode, CompliantTopologyEdge, CompliantTopologyBuild
export CompliantJointActuator, CompliantJointLoad
export rectangular_prism_inertia, thin_panel_inertia
export build_compliant_topology, build_rectangular_compliant_grid
export compliant_state_vector, compliant_state_parts, compliant_multibody_dynamics, compliant_joint_loads
export step_compliant_multibody_rk4, step_compliant_multibody_implicit_midpoint, simulate_compliant_multibody
export coupled_cloth_robot_arm_state_shape, initialize_coupled_cloth_robot_arm_state!
export assign_coupled_cloth_robot_arm_rhs!
export RobotArmHeldActuation, RobotArmJointMPCController, RobotArmControlEffector, RobotArmReactionEffector
export init_robot_arm_joint_mpc, robot_arm_joint_mpc_reference_preview
export robot_arm_joint_mpc_control, robot_arm_measured_joint_state
export AbstractTerrainModel, NoTerrainModel, DEMGrid, DEMTerrainModel
export terrain_height, terrain_radius, load_dem_grid, load_site_terrain, dem_grid_covers
export NoAtmosphereModel, ExponentialAtmosphereModel, PiecewiseExponentialAtmosphereModel
export GRAMGridAtmosphereModel
export SurrogatePresetResolution, available_surrogate_presets, resolve_surrogate_preset, surrogate_preset_model, atmosphere_provenance
export NRLMSISE00AtmosphereModel, init_nrlmsise_space_indices!
export SimpleEphemeridesModel
export make_no_gram_planet, make_no_gram_density_model, make_no_gram_environment
export calcForceTorque, wrench, environment_requirements, solver_partition
export gravity_backbone_structure, gravity_backbone_acceleration_ii
export gravity_backbone_kick_structure, gravity_backbone_kick_acceleration_ii
export PlumeSurfaceConfig, PlumeSurfaceInteractionModel, PlumeSurfaceState, plume_surface_footprint, plume_erosion_onset_height, plume_ground_effect_force, plume_quantities
export DescentPhaseTargets, apollo11_descent_targets, ApolloDescentConfig, ApolloDescentState, ApolloDescentGuidanceModel
export ApolloDescentControlConfig, ApolloDescentControlModel, descent_attitude_command
export getDensity, getDensityBatch!, with_density_model_epoch
export calcControlEffect!, calcControlForceTorque, calcControlMassFlowRate
export control_thruster_levels
export AerobrakingEnergyDepletionConfig, AerobrakingEnergyDepletionState
export AerobrakingEnergyDepletionGuidanceModel, AerobrakingEnergyDepletionControlModel
export SolarPanelAngleOfAttackControlModel
export ApoapsisTargetPeriapsisRaiseGuidanceModel
export VerificationRequest, VerificationResult
export run_verification, run_verification_cli, run_study, run_simulation
export station_geometry_path, station_cad_path, load_rpo_station_pointcloud, load_rpo_station_cad_triangles, load_rpo_station_cad_pointcloud
export VisualizationScene, PlanetSpec, SpacecraftGeometry, LinkBox, AtmosphereSpec, atmosphere_spec, ArmGeometry, arm_geometry
export load_model_triangles, model_bounding_box, sample_model_pointcloud, articulate_triangles, articulation_payload
export spacecraft_geometry, planet_spec, planet_rotation_table, build_visualization_scene
export visualization_scene_path, write_visualization_scene, read_visualization_scene
export velocity_aligned_quaternion, visualization_frame_budget
export MeshAeroPanels, MeshAeroSurrogate, AerodynamicCoefficientMeshSurrogate, MESH_AERO_MAX_DEGREE
export mesh_aero_panels, panel_aero_coefficients, panel_aero_coefficients_split, panel_shadow_mask, panel_projected_area
export fit_mesh_aero_surrogate, mesh_aero_coefficients, write_mesh_aero_surrogate, read_mesh_aero_surrogate
export export_visualization, with_visualization_scene, write_viewer_dev_payload
export EnsembleSample, sample_results_directory, with_results_directory, write_ensemble_manifest
export export_ensemble_visualization, run_monte_carlo_visualization
export AssetCheckItem, AssetCheckReport, check_assets, render_asset_report, run_cli


## 5. Precompile Workload
using PrecompileTools: @compile_workload, @setup_workload
include(joinpath(@__DIR__, "precompile_workload.jl"))

# The Monte Carlo dispatchers compile on their first campaign in a process --
# the job channel, the feeders and local consumers of the mixed dispatcher, the
# sample wrapper, the steady-cost estimator. Measured on the paper harness
# (L12, independent_1sat_1hr, 64 samples): the runner's first pool campaign
# cost 3.1-3.2 s against 1.8-2.2 s for the static pool path's own cold start
# on both machines, and 0.2-0.6 s warm. A production process pays that once;
# the harness pays it on the first repeat of every point. Exercised with a
# trivial sample so the generic machinery is in the pkgimage; the user's sample
# closure itself still specialises on first call. The body lives in
# `SimulationCampaigns._warm_campaign_dispatchers` so the test suite can run
# the same code at run time.
@compile_workload begin
	SimulationCampaigns._warm_campaign_dispatchers()
end

## 6. Runtime Initialization
# Runtime wiring that must not be baked into the precompiled image: these Refs
# hold closures over EnvironmentModels functions, so assigning them at include
# time would serialize a closure from an earlier world age. __init__ runs on
# every load of the cached image, which is what this needs.
function __init__()
	try
		SimulationModel.SimulationCallbacks._install_density_service_hooks!()
	catch err
		@warn "Could not install distributed density service hooks; the service will be unavailable." exception=(err, catch_backtrace())
	end
	return nothing
end


end # module SpaceAGORA
