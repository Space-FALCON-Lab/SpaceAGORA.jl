__precompile__(true)

module SpaceAGORA

## 1. Include package modules
include(joinpath(@__DIR__, "models", "gnc", "contracts", "gnc_algorithms.jl"))
include(joinpath(@__DIR__, "models", "gnc", "planning", "rpo", "contract", "rpo_planner.jl"))
include(joinpath(@__DIR__, "models", "gnc", "contracts", "gnc_runtime.jl"))
include(joinpath(@__DIR__, "models", "gnc", "contracts", "simulation_lifecycle.jl"))
include(joinpath(@__DIR__, "orchestration", "parallel", "routing", "parallel_profiles.jl"))
include(joinpath(@__DIR__, "orchestration", "parallel", "process", "parallel_process.jl"))
include(joinpath(@__DIR__, "interfaces", "analysis", "visualization", "visualizer_services.jl"))
include(joinpath(@__DIR__, "orchestration", "simulation", "runtime_services.jl"))
include(joinpath(@__DIR__, "orchestration", "simulation", "callback_services.jl"))
include(joinpath(@__DIR__, "orchestration", "simulation", "composition", "simulation_model.jl"))
include(joinpath(@__DIR__, "models", "gnc", "planning", "rpo", "direct_planner.jl"))
include(joinpath(@__DIR__, "models", "gnc", "companions", "hypr", "rpo_planner_adapter.jl"))
include(joinpath(@__DIR__, "models", "gnc", "planning", "rpo", "comparison", "planner_contract.jl"))
include(joinpath(@__DIR__, "models", "gnc", "planning", "rpo", "lifecycle.jl"))
include(joinpath(@__DIR__, "models", "gnc", "contracts", "algorithm_lifecycle.jl"))
include(joinpath(@__DIR__, "models", "gnc", "control", "maneuvers", "burn_planning.jl"))
include(joinpath(@__DIR__, "models", "gnc", "planning", "robot_arm", "joint_reference_family.jl"))
include(joinpath(@__DIR__, "models", "gnc", "guidance", "aerobraking", "edg_guidance_family.jl"))
include(joinpath(@__DIR__, "models", "gnc", "companions", "hypr", "services.jl"))
include(joinpath(@__DIR__, "orchestration", "simulation", "engine", "simulation_engine.jl"))
include(joinpath(@__DIR__, "orchestration", "simulation", "campaigns", "simulation_campaigns.jl"))
include(joinpath(@__DIR__, "orchestration", "mission", "mission_workflow.jl"))
include(joinpath(@__DIR__, "interfaces", "analysis", "verification", "telemetry_verification.jl"))
include(joinpath(@__DIR__, "foundation", "assets", "rpo_station_assets.jl"))
include(joinpath(@__DIR__, "foundation", "assets", "odyssey_surrogate_assets.jl"))
include(joinpath(@__DIR__, "interfaces", "analysis", "visualization", "rpo", "rpo_visualization.jl"))
include(joinpath(@__DIR__, "foundation", "assets", "asset_catalog.jl"))
include(joinpath(@__DIR__, "interfaces", "cli", "spaceagora_cli.jl"))
include(joinpath(@__DIR__, "orchestration", "simulation", "external_propagation_services.jl"))
public ExternalPropagationServices
include(joinpath(@__DIR__, "models", "gnc", "companions", "edg", "services.jl"))
public EDGServices
include(joinpath(@__DIR__, "interfaces", "analysis", "conformance", "conformance.jl"))
public Conformance



using .MissionWorkflow: MissionPlan, MissionPhase, MissionEvent, MissionRun, run_mission, MissionSegment, MissionDeadline, MissionDuration, MissionCondition
export MissionPlan, MissionPhase, MissionEvent, MissionRun, run_mission, MissionSegment, MissionDeadline, MissionDuration, MissionCondition
for name in (:MissionPlan, :MissionPhase, :MissionEvent, :MissionRun, :run_mission, :MissionSegment, :MissionDeadline, :MissionDuration, :MissionCondition)
    @eval @doc (@doc MissionWorkflow.$name) $name
end
using .VisualizerServices: VisualizerUnavailableError
const visualizer_available = VisualizerServices.available
export VisualizerUnavailableError, visualizer_available
@doc (@doc VisualizerServices.VisualizerUnavailableError) VisualizerUnavailableError
@doc "Return whether a compatible optional SpaceAGORAVisualizer is loaded in this environment." visualizer_available
using .SimulationModel.SceneVisualization: export_viewer_payload
export export_viewer_payload
@doc (@doc SimulationModel.SceneVisualization.export_viewer_payload) export_viewer_payload

# Opt-in public RPO planner pilot.
using .RPOPlannerInterfaces: AbstractRPOPlanner, RPOPlanningConstraints, RPOValidationSettings, RPOPlanningRequest, RPOReference, RPOPlanningResult, RPOPlannerCapabilities, RPOValidationResult, RPOPlanningHeadroom, planner_capabilities, initialize_planner, plan_rpo!, retime_rpo!, validate_rpo_result, validate_rpo_capabilities, rpo_reference_is_current
using .DirectRPOPlanning: DirectRPOPlanner
using .RPOPlannerInterfaces: RPOReferenceFamily
using .RPOPlannerInterfaces: AbstractRPOGeometryProvider, rpo_geometry_revision, rpo_geometry_snapshot
using .GNCAlgorithms: bound_family
using .GNCAlgorithms: AbstractGNCAlgorithm, AbstractAlgorithmFamily, AbstractAlgorithmRequest, AbstractAlgorithmResult, AbstractValidationResult, GNCAlgorithmError, AlgorithmResult, ValidationResult
using .GNCAlgorithms: algorithm_family, default_capabilities, algorithm_capabilities, require_algorithm_support, initialize_algorithm, propose!, propose_revision!, consumes_rng
using .GNCAlgorithms: validate_capabilities, validate_proposal, install!, is_current, on_rejection, fallback_product, family_contract_version
using .GNCAlgorithms: gnc_stream_words, AlgorithmRuntime, run_algorithm_request!, algorithm_runtime, bound_algorithm, bound_link, build_algorithm_request, trusted_queries, fail_algorithm_request!, algorithm_result_type, is_revision_request, proposed_product, installed_record, rejected_record
using .GNCRuntime: tick_schedule, gnc_callbacks, register_tstops!, gnc_clock_started!, gnc_period_s
using .GNCRuntime: Produces, Consumes, gnc_links, link_cardinality, check_link, link_initialize!
using .GNCRuntime: SpacecraftScope, AllSpacecraft, BoundSpacecraft, CoupledSpacecraft, spacecraft_scope, skip_inactive, held_command, held_command_time, release_command!
using .HYPRRPOPlanning: HYPRRPOPlanner
using .HYPRRPOPlanning: RetainedPathRPOPlanner
using .SimulationModel.HYPRSupport: HYPRUnavailableError, hypr_available
export HYPRUnavailableError, hypr_available
@doc (@doc SimulationModel.HYPRSupport.hypr_available) hypr_available
@doc (@doc SimulationModel.HYPRSupport.HYPRUnavailableError) HYPRUnavailableError
using .SimulationModel: RPOPSOConfig, rpo_pso_config, SimulationSettings
export SimulationSettings
# --- Simulation setup API: `using SpaceAGORA` alone is enough to configure a run ---
using .SimulationModel: AerodynamicCoefficientfM
export AerodynamicCoefficientfM
using .SimulationModel: SolarRadiationPressureModel
export SolarRadiationPressureModel
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.AerodynamicCoefficientfM) AerodynamicCoefficientfM
using .SimulationModel: SimulationConfiguration, MissionConfiguration, MissionType, MissionTime, MissionOrbits, EnvironmentModel, DynamicsModel, GuidanceModel, NavigationModel, ControlModel, SpacecraftModel, Link, Joint, InitialCondition, CartesianInitialCondition, InitialTime, IntegrationTolerances, SaveField, default_save_fields, MaxwellianHeat, InverseSquaredGravityModel, InverseSquaredJ2GravityModel, NBodyGravityModel, GravitationalHarmonicsModel, make_example_config, make_three_body_spacecraft
export SimulationConfiguration, MissionConfiguration, MissionType, MissionTime, MissionOrbits, EnvironmentModel, DynamicsModel, GuidanceModel, NavigationModel, ControlModel, SpacecraftModel, Link, Joint, InitialCondition, CartesianInitialCondition, InitialTime, IntegrationTolerances, SaveField, default_save_fields, MaxwellianHeat, InverseSquaredGravityModel, InverseSquaredJ2GravityModel, NBodyGravityModel, GravitationalHarmonicsModel, make_example_config, make_three_body_spacecraft
using .SimulationModel: AbstractInterlink, LinkContext, LaserMomentumLink, link_range_m, link_available, link_force
using .SimulationModel: interlink_wrench, link_applies_torque
using .SimulationModel: InterlinkModel, InterLinkGraph, InterLinkKey, InterLinkHistory, register_link!, interlink_force_on_spacecraft
using .SimulationModel: link_directed, link_terminal_levels, link_terminal, terminal_supports, terminal_count
using .SimulationModel: register_links!, all_pairs, mounted_link_terminals, terminal_number, terminal_pose
using .SimulationModel: same_plane_pairs, ring_pairs, grid_pairs
using .SimulationModel: GroundStation, add_ground_station!, station_pairs
using .SimulationModel: PowerBeamLink, link_received_power_w, beam_spot_radius_m, beam_capture_fraction
using .SimulationModel: LaserBeamLink, AbstractLaserBeam, GaussianBeam, TerminalPointing, beam_force_factor
using .SimulationModel: DataLink, AbstractDataBudget, FixedRateBudget, ShannonBudget, link_data_rate_bps, budget_rate_bps
using .SimulationModel: terminal_shed_tier, contact_windows, contact_time_s
using .SimulationModel: SchedulingPolicyModel, update_availability!, score_candidates!, select_interlinks!, schedule_interlinks!, interlink_scheduler_callback
using .SimulationModel: AbstractSchedulingPolicy, LinkCandidate, SchedulingContext, candidate_score
using .SimulationModel: SemimajorAxisRatePolicy, EccentricityRatePolicy, PowerDeliveryPolicy, DataThroughputPolicy
using .SimulationModel: EnergyBalancingPolicy, WeightedPolicy, PolicyTerm, state_of_charge
using .SimulationModel: DataBacklogPolicy, stored_bits, store_room_bits
using .SimulationModel: ConstellationInterLinkModel, InterLinkParameters, register_allowed_link!, force_on_endpoint, laser_force_on_spacecraft
export AbstractInterlink, LinkContext, LaserMomentumLink, link_range_m, link_available, link_force
export interlink_wrench, link_applies_torque
export InterlinkModel, InterLinkGraph, InterLinkKey, InterLinkHistory, register_link!, interlink_force_on_spacecraft
export link_directed, link_terminal_levels, link_terminal, terminal_supports, terminal_count
export register_links!, all_pairs, mounted_link_terminals, terminal_number, terminal_pose
export same_plane_pairs, ring_pairs, grid_pairs
export GroundStation, add_ground_station!, station_pairs
export PowerBeamLink, link_received_power_w, beam_spot_radius_m, beam_capture_fraction
export LaserBeamLink, AbstractLaserBeam, GaussianBeam, TerminalPointing, beam_force_factor
export DataLink, AbstractDataBudget, FixedRateBudget, ShannonBudget, link_data_rate_bps, budget_rate_bps
export terminal_shed_tier, contact_windows, contact_time_s
export SchedulingPolicyModel, update_availability!, score_candidates!, select_interlinks!, schedule_interlinks!, interlink_scheduler_callback
export AbstractSchedulingPolicy, LinkCandidate, SchedulingContext, candidate_score
export SemimajorAxisRatePolicy, EccentricityRatePolicy, PowerDeliveryPolicy, DataThroughputPolicy
export EnergyBalancingPolicy, WeightedPolicy, PolicyTerm, state_of_charge
export DataBacklogPolicy, stored_bits, store_room_bits
export ConstellationInterLinkModel, InterLinkParameters, register_allowed_link!, force_on_endpoint, laser_force_on_spacecraft
@doc (@doc SimulationModel.InterLinkModels.AbstractInterlink) AbstractInterlink
@doc (@doc SimulationModel.InterLinkModels.LinkContext) LinkContext
@doc (@doc SimulationModel.InterLinkModels.LaserMomentumLink) LaserMomentumLink
@doc (@doc SimulationModel.InterLinkModels.link_range_m) link_range_m
@doc (@doc SimulationModel.InterLinkModels.link_available) link_available
@doc (@doc SimulationModel.InterLinkModels.link_force) link_force
@doc (@doc SimulationModel.InterLinkModels.interlink_wrench) interlink_wrench
@doc (@doc SimulationModel.InterLinkModels.link_applies_torque) link_applies_torque
@doc (@doc SimulationModel.InterLinkModels.InterlinkModel) InterlinkModel
@doc (@doc SimulationModel.InterLinkModels.register_link!) register_link!
@doc (@doc SimulationModel.InterLinkModels.interlink_force_on_spacecraft) interlink_force_on_spacecraft
@doc (@doc SimulationModel.InterLinkModels.update_availability!) update_availability!
@doc (@doc SimulationModel.InterLinkModels.ConstellationInterLinkModel) ConstellationInterLinkModel
@doc (@doc SimulationModel.InterLinkModels.InterLinkParameters) InterLinkParameters
@doc (@doc SimulationModel.InterLinkModels.register_allowed_link!) register_allowed_link!
@doc (@doc SimulationModel.InterLinkModels.link_directed) link_directed
@doc (@doc SimulationModel.InterLinkModels.link_terminal_levels) link_terminal_levels
@doc (@doc SimulationModel.InterLinkModels.link_terminal) link_terminal
@doc (@doc SimulationModel.InterLinkModels.terminal_supports) terminal_supports
@doc (@doc SimulationModel.InterLinkModels.terminal_count) terminal_count
@doc (@doc SimulationModel.InterLinkModels.register_links!) register_links!
@doc (@doc SimulationModel.InterLinkModels.all_pairs) all_pairs
@doc (@doc SimulationModel.InterLinkModels.same_plane_pairs) same_plane_pairs
@doc (@doc SimulationModel.InterLinkModels.PowerBeamLink) PowerBeamLink
@doc (@doc SimulationModel.InterLinkModels.link_received_power_w) link_received_power_w
@doc (@doc SimulationModel.InterLinkModels.beam_spot_radius_m) beam_spot_radius_m
@doc (@doc SimulationModel.InterLinkModels.beam_capture_fraction) beam_capture_fraction
@doc (@doc SimulationModel.InterLinkModels.LaserBeamLink) LaserBeamLink
@doc (@doc SimulationModel.InterLinkModels.AbstractLaserBeam) AbstractLaserBeam
@doc (@doc SimulationModel.InterLinkModels.GaussianBeam) GaussianBeam
@doc (@doc SimulationModel.InterLinkModels.TerminalPointing) TerminalPointing
@doc (@doc SimulationModel.InterLinkModels.beam_force_factor) beam_force_factor
@doc (@doc SimulationModel.InterLinkModels.DataLink) DataLink
@doc (@doc SimulationModel.InterLinkModels.AbstractDataBudget) AbstractDataBudget
@doc (@doc SimulationModel.InterLinkModels.FixedRateBudget) FixedRateBudget
@doc (@doc SimulationModel.InterLinkModels.ShannonBudget) ShannonBudget
@doc (@doc SimulationModel.InterLinkModels.link_data_rate_bps) link_data_rate_bps
@doc (@doc SimulationModel.InterLinkModels.budget_rate_bps) budget_rate_bps
@doc (@doc SimulationModel.InterLinkModels.terminal_shed_tier) terminal_shed_tier
@doc (@doc SimulationModel.InterLinkModels.contact_windows) contact_windows
@doc (@doc SimulationModel.InterLinkModels.contact_time_s) contact_time_s
@doc (@doc SimulationModel.InterLinkModels.ring_pairs) ring_pairs
@doc (@doc SimulationModel.InterLinkModels.grid_pairs) grid_pairs
@doc (@doc SimulationModel.InterLinkModels.GroundStation) GroundStation
@doc (@doc SimulationModel.InterLinkModels.add_ground_station!) add_ground_station!
@doc (@doc SimulationModel.InterLinkModels.station_pairs) station_pairs
@doc (@doc SimulationModel.InterlinkTerminals.mounted_link_terminals) mounted_link_terminals
@doc (@doc SimulationModel.InterlinkTerminals.terminal_number) terminal_number
@doc (@doc SimulationModel.InterLinkModels.terminal_pose) terminal_pose
@doc """Indexed terminal vertices, candidate edges, the link of each edge (`edge_properties`) and availability/selection flags. Register candidates with `register_link!`.""" InterLinkGraph
@doc """Canonically ordered pair of `(spacecraft array index, terminal number)` endpoints; a ground station end is `(-s, antenna number)`.""" InterLinkKey
@doc """One accepted-step schedule sample: time in seconds, the selected `InterLinkKey` values (`active`), the available ones (`available`, recorded in a run with a data link or with `InterlinkModel(; record_available=true)`), the score that selected each active link (`scores`), and, with `InterlinkModel(; record_terms=true)` and a `WeightedPolicy`, each top-level term's score for each active link (`terms`, one row per term, one column per active link). A sample built from a time and the selected keys records neither scores nor available links. `contact_windows` and `link_schedule` read a run's samples.""" InterLinkHistory
@doc (@doc SimulationModel.SchedulingPolicies.SchedulingPolicyModel) SchedulingPolicyModel
@doc (@doc SimulationModel.SchedulingPolicies.score_candidates!) score_candidates!
@doc (@doc SimulationModel.SchedulingPolicies.AbstractSchedulingPolicy) AbstractSchedulingPolicy
@doc (@doc SimulationModel.SchedulingPolicies.LinkCandidate) LinkCandidate
@doc (@doc SimulationModel.SchedulingPolicies.SchedulingContext) SchedulingContext
@doc (@doc SimulationModel.SchedulingPolicies.candidate_score) candidate_score
@doc (@doc SimulationModel.SchedulingPolicies.SemimajorAxisRatePolicy) SemimajorAxisRatePolicy
@doc (@doc SimulationModel.SchedulingPolicies.EccentricityRatePolicy) EccentricityRatePolicy
@doc (@doc SimulationModel.SchedulingPolicies.PowerDeliveryPolicy) PowerDeliveryPolicy
@doc (@doc SimulationModel.SchedulingPolicies.DataThroughputPolicy) DataThroughputPolicy
@doc (@doc SimulationModel.SchedulingPolicies.EnergyBalancingPolicy) EnergyBalancingPolicy
@doc (@doc SimulationModel.SchedulingPolicies.WeightedPolicy) WeightedPolicy
@doc (@doc SimulationModel.SchedulingPolicies.PolicyTerm) PolicyTerm
@doc (@doc SimulationModel.SchedulingPolicies.state_of_charge) state_of_charge
@doc (@doc SimulationModel.SchedulingPolicies.DataBacklogPolicy) DataBacklogPolicy
@doc (@doc SimulationModel.SchedulingPolicies.stored_bits) stored_bits
@doc (@doc SimulationModel.SchedulingPolicies.store_room_bits) store_room_bits
@doc (@doc SimulationModel.SchedulingPolicies.select_interlinks!) select_interlinks!
@doc """Refresh link availability, score candidates, and select the maximum-total-score terminal matching.""" schedule_interlinks!
@doc """Create the accepted-step scheduler callback used automatically by `run_simulation` when interlinks are enabled.""" interlink_scheduler_callback
@doc (@doc SimulationModel.DynamicEffectors.LaserForceEffectors.force_on_endpoint) force_on_endpoint
@doc (@doc SimulationModel.DynamicEffectors.LaserForceEffectors.laser_force_on_spacecraft) laser_force_on_spacecraft
@doc (@doc SimulationModel.ExampleConfiguration.make_example_config) make_example_config
@doc (@doc SimulationModel.ExampleConfiguration.make_three_body_spacecraft) make_three_body_spacecraft
@doc """Dynamics model of a run: the spacecraft list and the tuple of force/torque effectors, `DynamicsModel([spacecraft], effectors)`.""" DynamicsModel
@doc """Initial orbit of a spacecraft from Keplerian elements (`a` and `e`, or `ra` and `rp`; `i`, `ω`, `Ω`, and the true anomaly `ν` or the mean anomaly `M`, in degrees), with an optional initial attitude `q` (x, y, z, w, inertial to body) and body rate `ang_vel` (rad/s), each any vector or tuple of 4 or 3 numbers.""" InitialCondition
@doc """Complete one-run scenario: simulation settings, mission configuration, environment, dynamics, guidance, navigation and control models, initial time, and integration tolerances. Pass it to `run_simulation`.""" SimulationConfiguration
@doc """Mission length and bookkeeping: `mission_type` (`MissionTime` or `MissionOrbits`), `mission_time`, `number_of_orbits`, `keplerian`, `orientation_sim`, and `num_steps_to_save`.""" MissionConfiguration
@doc """Enum selecting how a mission ends: `MissionTime` (fixed duration) or `MissionOrbits` (fixed orbit count).""" MissionType
@doc """`MissionType` value: the mission ends after `MissionConfiguration.mission_time` seconds.""" MissionTime
@doc """`MissionType` value: the mission ends after `MissionConfiguration.number_of_orbits` orbits.""" MissionOrbits
@doc """Planet, atmosphere density model, ephemerides model, thermal model, entry-interface altitude `EI` (km), and topography/wind switches for a run.""" EnvironmentModel
@doc """Guidance effectors and their update rates for a run (`guidance_effectors`, `guidance_rates`).""" GuidanceModel
@doc """Navigation effectors and their update rates for a run (`navigation_effectors`, `navigation_rates`).""" NavigationModel
@doc """Control effectors and their update rates for a run (`control_effectors`, `control_rates`).""" ControlModel
@doc """Multi-link spacecraft: joints, links, root link, mass properties, propellant mass, initial condition, id, and its `constellation_slot` (a `ConstellationSlot`, set by `build_constellation`, or `nothing`). Build the links with `Link` and `Joint`, or use `make_three_body_spacecraft`.""" SpacecraftModel
@doc """One rigid body of a `SpacecraftModel`: mass, dimensions, reference area, offset `r`, attitude `q`, and reflection coefficient. Exported at the root; note that other packages may also export a `Link` name.""" Link
@doc """Articulated connection between two `Link`s of a `SpacecraftModel`. Exported at the root; note that other packages may also export a `Joint` name.""" Joint
@doc """Initial spacecraft state given as a Cartesian position and velocity, with an optional initial attitude `q` (x, y, z, w, inertial to body) and body rate `ang_vel` (rad/s), each any vector or tuple of 4 or 3 numbers.""" CartesianInitialCondition
@doc """Epoch at which a simulation starts (calendar date and time).""" InitialTime
@doc """ODE solver tolerances and maximum step sizes for the orbit and atmosphere phases (`reltol_*`, `abstol_*`, `dt_max_*`).""" IntegrationTolerances
@doc """Declaration of one quantity written to the saved results.""" SaveField
@doc """    default_save_fields(config::SimulationConfiguration)

Return the default `SaveField`s for `config`.""" default_save_fields
@doc """Free-molecular heating model with a thermal accommodation factor (`thermal_accomodation_factor`) for a given planet.""" MaxwellianHeat
@doc """Point-mass (inverse-square) gravity effector.""" InverseSquaredGravityModel
@doc """Point-mass gravity with the J2 oblateness term.""" InverseSquaredJ2GravityModel
@doc """Third-body gravity effector using ephemerides for the listed bodies.""" NBodyGravityModel
@doc """
    SolarRadiationPressureModel(Cr, A, AU_m=149_597_870_700.0; kwargs...)

Solar radiation pressure for reflectivity coefficient `Cr` and reference area
`A` in square metres. Direct sunlight is enabled by default (`direct=true`);
planetary albedo and infrared contributions are opt-in (`albedo=false`,
`ir=false`). Their parameters are `planet_albedo=0.3` and
`planet_ir_flux_w_m2=237.0`. Use this effector in `DynamicsModel` or the
`dynamic_effectors` argument of `make_example_config`. Active direct sunlight
or albedo requires `SpiceEphemeridesModel` and kernels for the selected planet
and epoch; a run using `SimpleEphemeridesModel` rejects that combination.
""" SolarRadiationPressureModel
@doc """
    GravitationalHarmonicsModel(L, M, coefficients_file, planet; kwargs...)

Spherical-harmonic gravity to degree `L` and order `M` (`M <= L`) from a
coefficients file with columns `degree,order,C,S`. The central point-mass term
is included by default (`include_central=true`), so this model replaces
`InverseSquaredGravityModel` or `InverseSquaredJ2GravityModel` rather than being
added to them. `coefficient_normalization` names the file's convention (`:full`,
the default, `:schmidt` or `:unnormalized`). `coefficients_normalized=false`
also converts unnormalized coefficients, so use one of the two, not both, or
the coefficients are converted twice. The other keywords are `j2_source`,
`reference_radius_m`, `gm_m3s2` and `include_central`. A degree beyond the
file's contents is refused.
""" GravitationalHarmonicsModel
@doc "Output and checkpoint settings for a simulation. The RPO pilot runs with checkpoints. What a checkpoint supports, what it refuses and what lies outside its guarantee is listed in the checkpoint table of the RPO planner guide, \"Run an RPO case with interchangeable planners\"." SimulationSettings
@doc """
    RPOPSOConfig(; kwargs...)

Existing HYPR configuration, exported for explicit `HYPRRPOPlanner` assembly.
It controls swarm counts, geometric cost/search, adaptation and reference retiming.
The adapter copies it, maps the request's clearance/interval, applies the declared
planning reserve, and retains optimizer-returned settings. The public pilot example
shows a deterministic bounded configuration; its numerical limits are not defaults
for arbitrary missions. Existing internal constructors retain their behavior.
""" RPOPSOConfig
@doc """
    rpo_pso_config(config; overrides...)

Copy the existing HYPR configuration with explicit field overrides. This root alias
preserves the existing configuration owner and calculations; pass the result to
`HYPRRPOPlanner`. Request-owned limits and interval still apply at planning time.
""" rpo_pso_config

# The generic algorithm binding and the burn-plan family.
using .GNCAlgorithms: family_request, family_link, preflight_family
using .AlgorithmLifecycle: AlgorithmEvent, AlgorithmGuidance, AlgorithmGuidanceError, algorithm_run_report
using .BurnPlanning: AbstractBurnPlanner, BurnPlanFamily, BurnPlanRequest, BurnPlanSlot, ApoapsisBurnPlanner, validate_burn_plan, burn_plan_is_current, burn_plan_sizing
using .SimulationModel.ControlHooks: burn_plan_rejections
@doc (@doc SimulationModel.ControlHooks.burn_plan_rejections) burn_plan_rejections
# The burn scheduler and its burn plan, which a burn-plan run configures and a planner proposes.
using .SimulationModel: BaseThrusterModel, PropulsiveBurnPlan
@doc (@doc SimulationModel.BaseThrusterModel) BaseThrusterModel
@doc (@doc SimulationModel.PropulsiveBurnPlan) PropulsiveBurnPlan
export family_request, family_link, preflight_family
export AlgorithmEvent, AlgorithmGuidance, AlgorithmGuidanceError, algorithm_run_report
export AbstractBurnPlanner, BurnPlanFamily, BurnPlanRequest, BurnPlanSlot, ApoapsisBurnPlanner, validate_burn_plan, burn_plan_is_current, burn_plan_rejections, burn_plan_sizing
export BaseThrusterModel, PropulsiveBurnPlan

# The joint-reference family and its link to the robot-arm controller.
using .JointReferencePlanning: AbstractJointReferencePlanner, JointReferenceFamily, JointReferenceRequest, JointReferenceSettings, QuinticJointPlanner, validate_joint_reference, joint_reference_is_current
using .SimulationModel: RobotArmPlanBuffer
export AbstractJointReferencePlanner, JointReferenceFamily, JointReferenceRequest, JointReferenceSettings, QuinticJointPlanner, validate_joint_reference, joint_reference_is_current, RobotArmPlanBuffer

# The EDG guidance family, served by the EDGServices facade.
using .EDGGuidanceFamilies: EDGGuidanceFamily, EDGGuidanceCapabilities, EDGGuidanceAlgorithm, EDGGuidanceSettings, EDGGuidanceRequest, EDGGuidanceDecision, validate_edg_guidance_decision
export EDGGuidanceFamily, EDGGuidanceCapabilities, EDGGuidanceAlgorithm, EDGGuidanceSettings, EDGGuidanceRequest, EDGGuidanceDecision, validate_edg_guidance_decision

using .RPOPlannerLifecycle: RPOPlanningEvent, RPOPlanningError, make_rpo_configuration, rpo_run_report
using .RPOPlannerLifecycle: RPOReplanningSpheres
using .SimulationModel: HYPRComparisonPlanner, UnrefinedPSOPlanner, RRTConnectPlanner, RRTConnectBezierPlanner, RRTStarPlanner, CHOMPPlanner, STOMPPlanner
export HYPRComparisonPlanner, UnrefinedPSOPlanner, RRTConnectPlanner, RRTConnectBezierPlanner, RRTStarPlanner, CHOMPPlanner, STOMPPlanner
for name in (:HYPRComparisonPlanner, :UnrefinedPSOPlanner, :RRTConnectPlanner, :RRTConnectBezierPlanner, :RRTStarPlanner, :CHOMPPlanner, :STOMPPlanner)
    @eval @doc (@doc SimulationModel.GuidanceHooks.$name) $name
end
export AbstractGNCAlgorithm, AbstractAlgorithmFamily, AbstractAlgorithmRequest, AbstractAlgorithmResult, AbstractValidationResult, GNCAlgorithmError, AlgorithmResult, ValidationResult, RPOReferenceFamily
export algorithm_family, default_capabilities, algorithm_capabilities, require_algorithm_support, initialize_algorithm, propose!, propose_revision!, consumes_rng
export validate_capabilities, validate_proposal, install!, is_current, on_rejection, fallback_product, family_contract_version
export gnc_stream_words, AlgorithmRuntime, run_algorithm_request!, algorithm_runtime, bound_algorithm, bound_link, build_algorithm_request, trusted_queries, fail_algorithm_request!, algorithm_result_type, is_revision_request, proposed_product, installed_record, rejected_record
export tick_schedule, gnc_callbacks, register_tstops!, gnc_clock_started!, gnc_period_s
export Produces, Consumes, gnc_links, link_cardinality, check_link, link_initialize!
export SpacecraftScope, AllSpacecraft, BoundSpacecraft, CoupledSpacecraft, spacecraft_scope, skip_inactive, held_command, held_command_time, release_command!
export AbstractRPOPlanner, RPOPlanningConstraints, RPOValidationSettings, RPOPlanningRequest, RPOReference, RPOPlanningResult, RPOPlannerCapabilities, RPOValidationResult, RPOPlanningHeadroom, planner_capabilities, initialize_planner, plan_rpo!, retime_rpo!, validate_rpo_result, validate_rpo_capabilities, rpo_reference_is_current, DirectRPOPlanner, HYPRRPOPlanner, RPOPSOConfig, rpo_pso_config, RPOPlanningEvent, RPOPlanningError, make_rpo_configuration, rpo_run_report
export bound_family, AbstractRPOGeometryProvider, rpo_geometry_revision, rpo_geometry_snapshot, RPOReplanningSpheres
export RetainedPathRPOPlanner

## 2. Bring needed names from package modules into the scope of SpaceAGORA.jl
# 2.1. Parallel Profiles
# The R0-R7 profiles and the outer-route bandit are internal and benchmark
# machinery, reachable qualified as `SpaceAGORA.ParallelProfiles.<name>`. The
# supported way to parallelize is `SolverConfig(parallel=true)`.

# 2.2. Parallel Process
using .ParallelProcess: ProcessPool, campaign_process_pool, ensure_process_workers!, shutdown_process_pool!, adopt_process_workers!, PROCESS_WARMUP

## 2.3. Simulation Engine
using .SimulationEngine: ParallelConfig, SolverConfig, RuntimePolicyConfig, ArtifactConfig, SimulationEngineConfig
using .SimulationEngine: simulation_engine_config_from_env
using .SimulationEngine: SimulationResults, link_schedule, run_table, long_table
using .SimulationEngine: ControlForecastSnapshot, capture_control_forecast, run_control_forecast, with_control_forecast_worker
using .SimulationEngine: prewarm_nbody_ephemeris_cache, load_nbody_ephemeris_cache!
run_simulation(args...; kwargs...) = SimulationEngine.run_simulation(args...; kwargs...)

## 2.4. Simulation Campaigns
using .SimulationCampaigns: MonteCarloSpec, MonteCarloSampleResult, MonteCarloResult, run_monte_carlo
using .SimulationCampaigns: run_constellation_ensemble, ConstellationEnsembleResult, read_constellation_results
using .SimulationCampaigns: run_monte_carlo_visualization

## 2.5. Simulation Model
using .SimulationModel.FrameTransforms: orbital_elements_to_cartesian, planet_fixed_to_geodetic
using .SimulationModel: PolynomialFitAtmosphereModel, density_altitude_derivative
using .SimulationModel: AerodynamicCommandedAreaDragModel, apply_panel_angle_of_attack!
export orbital_elements_to_cartesian, planet_fixed_to_geodetic
export PolynomialFitAtmosphereModel, density_altitude_derivative
export AerodynamicCommandedAreaDragModel, apply_panel_angle_of_attack!
@doc (@doc SimulationModel.FrameTransforms.orbital_elements_to_cartesian) orbital_elements_to_cartesian
@doc (@doc SimulationModel.FrameTransforms.planet_fixed_to_geodetic) planet_fixed_to_geodetic
@doc (@doc SimulationModel.EnvironmentModels.PolynomialFitAtmosphereModel) PolynomialFitAtmosphereModel
@doc (@doc SimulationModel.EnvironmentModels.density_altitude_derivative) density_altitude_derivative
@doc (@doc SimulationModel.DynamicEffectors.AerodynamicEffectors.AerodynamicCommandedAreaDragModel) AerodynamicCommandedAreaDragModel
@doc (@doc SimulationModel.ControlHooks.apply_panel_angle_of_attack!) apply_panel_angle_of_attack!
using .SimulationModel: StateAnchor, get_state_anchor_callback
using .SimulationModel.AbstractTypes: AbstractForceTorqueModel, AbstractPlanet, AbstractDensityModel
using .SimulationModel.AbstractTypes: AbstractControlEffectorModel, AbstractEphemeridesModel
using .SimulationModel.AbstractTypes: AbstractThermalModel, AbstractThrusterModel, AbstractGuidanceModel
using .SimulationModel.AbstractTypes: AbstractNavigationModel
using .SimulationModel.AbstractTypes: AbstractActuator, AbstractSensor
using .SimulationModel: StateSample, PlanetFrameSample, AtmosphereSample, SolarEphemerisSample
using .SimulationModel: ThirdBodyEphemerisSample, EnvironmentSample, EffectorEnvironmentRequirements
# Compliant multibody dynamics and robot kinematics. Historical Cloth names below
# remain aliases to the same objects for script, HYPR, and serialization compatibility.
using .SimulationModel: CompliantMultibody, CompliantRobotArmDynamics
public CompliantMultibody, CompliantRobotArmDynamics
@doc "Compliant lumped-body dynamics for bodies connected by joints with stiffness and damping. The historical module alias remains available as `SimulationModel.ClothMultibody`." CompliantMultibody
@doc "Compliant robot-arm realization, simulation, and spacecraft coupling. The historical module alias remains available as `SimulationModel.ClothRobotArmDynamics`." CompliantRobotArmDynamics
using .SimulationModel: RobotArmBasePose, RobotArmLink, RobotArmJoint, RobotArmModel, RobotArmPose, RobotArmState, default_robot_arm_model, robot_arm_fk, robot_arm_fk_state, robot_arm_end_effector_pose, robot_arm_ik, robot_arm_total_reach
export RobotArmBasePose, RobotArmLink, RobotArmJoint, RobotArmModel, RobotArmPose, RobotArmState, default_robot_arm_model, robot_arm_fk, robot_arm_fk_state, robot_arm_end_effector_pose, robot_arm_ik, robot_arm_total_reach
@doc (@doc SimulationModel.RobotArmBasePose) RobotArmBasePose
@doc (@doc SimulationModel.RobotArmLink) RobotArmLink
@doc (@doc SimulationModel.RobotArmJoint) RobotArmJoint
@doc (@doc SimulationModel.RobotArmModel) RobotArmModel
@doc (@doc SimulationModel.RobotArmPose) RobotArmPose
@doc (@doc SimulationModel.RobotArmState) RobotArmState
@doc (@doc SimulationModel.default_robot_arm_model) default_robot_arm_model
@doc (@doc SimulationModel.robot_arm_fk) robot_arm_fk
@doc (@doc SimulationModel.robot_arm_fk_state) robot_arm_fk_state
@doc (@doc SimulationModel.robot_arm_end_effector_pose) robot_arm_end_effector_pose
@doc (@doc SimulationModel.robot_arm_ik) robot_arm_ik
@doc (@doc SimulationModel.robot_arm_total_reach) robot_arm_total_reach
using .SimulationModel: CompliantRobotArmReferenceState, CompliantRobotArmSimulation, compliant_robot_arm_reference_state, compliant_robot_arm_multibody, compliant_robot_arm_initial_state, compliant_robot_arm_rest_quaternions, compliant_robot_arm_end_effector, compliant_robot_arm_actuators, coupled_compliant_robot_arm_state_shape, initialize_coupled_compliant_robot_arm_state!, assign_coupled_compliant_robot_arm_rhs!, simulate_compliant_robot_arm_plan
export CompliantRobotArmReferenceState, CompliantRobotArmSimulation, compliant_robot_arm_reference_state, compliant_robot_arm_multibody, compliant_robot_arm_initial_state, compliant_robot_arm_rest_quaternions, compliant_robot_arm_end_effector, compliant_robot_arm_actuators, coupled_compliant_robot_arm_state_shape, initialize_coupled_compliant_robot_arm_state!, assign_coupled_compliant_robot_arm_rhs!, simulate_compliant_robot_arm_plan
@doc (@doc SimulationModel.CompliantRobotArmReferenceState) CompliantRobotArmReferenceState
@doc (@doc SimulationModel.CompliantRobotArmSimulation) CompliantRobotArmSimulation
@doc (@doc SimulationModel.compliant_robot_arm_reference_state) compliant_robot_arm_reference_state
@doc (@doc SimulationModel.compliant_robot_arm_multibody) compliant_robot_arm_multibody
@doc (@doc SimulationModel.compliant_robot_arm_initial_state) compliant_robot_arm_initial_state
@doc (@doc SimulationModel.compliant_robot_arm_rest_quaternions) compliant_robot_arm_rest_quaternions
@doc (@doc SimulationModel.compliant_robot_arm_end_effector) compliant_robot_arm_end_effector
@doc (@doc SimulationModel.compliant_robot_arm_actuators) compliant_robot_arm_actuators
@doc (@doc SimulationModel.coupled_compliant_robot_arm_state_shape) coupled_compliant_robot_arm_state_shape
@doc (@doc SimulationModel.initialize_coupled_compliant_robot_arm_state!) initialize_coupled_compliant_robot_arm_state!
@doc (@doc SimulationModel.assign_coupled_compliant_robot_arm_rhs!) assign_coupled_compliant_robot_arm_rhs!
@doc (@doc SimulationModel.simulate_compliant_robot_arm_plan) simulate_compliant_robot_arm_plan

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
using .SimulationModel: CompliantJointActuator, CompliantJointLoad, CompliantAttachment
using .SimulationModel: rectangular_prism_inertia, thin_panel_inertia
using .SimulationModel: build_compliant_topology, build_rectangular_compliant_grid
using .SimulationModel: ConstellationSlot, OrbitalPlane, build_constellation, walker_planes, walker_constellation
using .SimulationModel: ScopedEffector, scoped, SpacecraftSelection
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
using .SimulationModel: GRAMGridAtmosphereModel, GRAMNearSurfaceAtmosphereModel, CombinedAtmosphereModel
using .SimulationModel.EnvironmentModels: SurrogatePresetResolution, available_surrogate_presets, resolve_surrogate_preset, surrogate_preset_model, atmosphere_provenance
using .SimulationModel: NRLMSISE00AtmosphereModel, init_nrlmsise_space_indices!
using .SimulationModel: SimpleEphemeridesModel, SpiceEphemeridesModel
using .SimulationModel.TerrainModels: AbstractTerrainModel, NoTerrainModel, DEMGrid, DEMTerrainModel
using .SimulationModel.TerrainModels: terrain_height, terrain_radius, load_dem_grid, load_site_terrain, dem_grid_covers
using .SimulationModel: make_no_gram_planet, make_no_gram_density_model, make_no_gram_environment
using .SimulationModel: wrench, environment_requirements, solver_partition
using .SimulationModel: gravity_backbone_structure, gravity_backbone_acceleration_ii
using .SimulationModel: gravity_backbone_kick_structure, gravity_backbone_kick_acceleration_ii
using .SimulationModel: PlumeSurfaceConfig, PlumeSurfaceInteractionModel, PlumeSurfaceState, plume_surface_footprint, plume_erosion_onset_height, plume_ground_effect_force, plume_quantities
using .SimulationModel: getDensity, getDensityBatch!, with_density_model_epoch
using .SimulationModel: ephemerides_frame_lpi, ephemerides_sun_ii, ephemerides_epoch_seconds
using .SimulationModel: ephemerides_requires_spice, ephemerides_cache_key
using .SimulationModel: DescentPhaseTargets, apollo11_descent_targets, ApolloDescentConfig, ApolloDescentState, ApolloDescentGuidanceModel
using .SimulationModel: ApolloDescentControlConfig, ApolloDescentControlModel, descent_attitude_command
using .SimulationModel: control_wrench, control_propellant_rate, control_wheel_torque
using .SimulationModel: control_thruster_levels, touchdown_spec
using .SimulationModel: getHeatRate
using .SimulationModel: update_control!, update_guidance!, update_navigation!
# Actuators: the contract (Actuators), the thruster kind (Components) and mounting (ActuatorMounting).
using .SimulationModel: Thruster, ReactionWheel, Magnetorquer, LaserTerminal, PowerReceiver, RadioTerminal, PanelHinge, ActuatorContext
using .SimulationModel: ThrusterGimbal, thrust_direction
using .SimulationModel: ControlMomentGyro, cmg_momentum, cmg_torque_axis, cmg_angle_saturate, cmg_steering
using .SimulationModel: CMGInnerGimbal, cmg_torque_axes
using .SimulationModel: actuator_wrench, actuator_limits, actuator_saturate, actuator_consumption, actuator_fraction, actuator_unit_wrench, actuation_mode, actuator_levels, allocate_actuators
using .SimulationModel: add_actuator!, mounted_actuators, actuator_context, allocation_matrix, mounted_wrench, six_axis_thrusters
using .SimulationModel: propellant_mass_rate, ActuatorRequirements, actuator_requirements, applied_actuator_levels, mounted_power
using .SimulationModel: SolarArray, Battery, PowerLoad, add_power!
using .SimulationModel: AbstractDataSource, DataStore, DataSource, add_data!, data_production_bps
using .SimulationModel: ActuatorFault, add_fault!
using .SimulationModel: RateGyro, StarTracker, Magnetometer, SensorTruth, sensor_measure, sensor_dimension, sensor_requirements, sensor_period_s, initialize_sensor!
using .SimulationModel: GNSSReceiver, SunSensor, Accelerometer, sensor_needs_specific_force, is_gravitational
using .SimulationModel: add_sensor!, mounted_sensors, sensor_measurement, sensor_measurement_time
using .SimulationLifecycle: bind_spacecraft
using .SimulationLifecycle: preflight_guidance, preflight_control, initialize_navigation!, initialize_guidance!, initialize_control!, checkpoint_state, restore_checkpoint_state!, requires_atmosphere_events, atmosphere_transition!, before_reference_control!
export preflight_guidance, preflight_control, initialize_navigation!, initialize_guidance!, initialize_control!, checkpoint_state, restore_checkpoint_state!, requires_atmosphere_events, atmosphere_transition!, before_reference_control!
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
using .SimulationModel.NavigationHooks: RPOStationGeometry, RPOCubeSatGeometry

## 2.8. RPO Visualization
using .RPOVisualization: rpo_path_plot, rpo_tracking_plot
# The figures come from the PlotlyJS package extension.
using .SimulationModel.GuidanceHooks.RPOComparisonReports: PlotlyJSUnavailableError
export PlotlyJSUnavailableError
@doc (@doc SimulationModel.GuidanceHooks.RPOComparisonReports.PlotlyJSUnavailableError) PlotlyJSUnavailableError

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
# `@doc X` returns one nested Markdown.MD per docstring of X. Documenter unwraps a
# lone nested MD, but it renders each of several as plain text and warns about Julia
# interpolation. A name with more than one docstring (a type and its constructor,
# several methods) is therefore forwarded through `_flat_doc`, which splices the
# blocks into one MD in the same order, with the rule between docstrings that the
# REPL shows. The result carries no `:results` metadata, which would make the REPL
# draw a rule between every block.
function _flat_doc(md::T) where {T}
    # `@doc X` gives a Markdown.MD only while the REPL docsystem is loaded (today through
    # DataFrames -> PrettyTables); otherwise it gives the first DocStr, forwarded as is.
    (hasfield(T, :content) && hasfield(T, :meta)) || return md
    blocks = Any[]
    function add_blocks!(content)
        for block in content
            block isa T ? add_blocks!(block.content) : push!(blocks, block)
        end
    end
    for (i, doc) in enumerate(md.content)
        i > 1 && push!(blocks, parentmodule(T).HorizontalRule())
        add_blocks!((doc,))
    end
    return T(blocks)
end

# 3.1. Simulation Engine
@doc (@doc SimulationEngine.ParallelConfig) ParallelConfig
@doc (@doc SimulationEngine.SolverConfig) SolverConfig
@doc (@doc SimulationEngine.RuntimePolicyConfig) RuntimePolicyConfig
@doc (@doc SimulationEngine.ArtifactConfig) ArtifactConfig
@doc (@doc SimulationEngine.SimulationEngineConfig) SimulationEngineConfig
@doc (@doc SimulationEngine.simulation_engine_config_from_env) simulation_engine_config_from_env
@doc (@doc SimulationEngine.SimulationResults) SimulationResults
@doc (@doc SimulationEngine.ControlForecastSnapshot) ControlForecastSnapshot
@doc (@doc SimulationEngine.capture_control_forecast) capture_control_forecast
@doc (@doc SimulationEngine.run_control_forecast) run_control_forecast
@doc (@doc SimulationEngine.with_control_forecast_worker) with_control_forecast_worker
@doc (@doc SimulationEngine.link_schedule) link_schedule
@doc (@doc SimulationEngine.run_table) run_table
@doc (@doc SimulationEngine.long_table) long_table
@doc (@doc SimulationEngine.run_simulation) run_simulation
@doc (@doc SimulationEngine.prewarm_nbody_ephemeris_cache) prewarm_nbody_ephemeris_cache
@doc (@doc SimulationEngine.load_nbody_ephemeris_cache!) load_nbody_ephemeris_cache!

# 3.2. Simulation Campaigns
@doc (@doc SimulationCampaigns.MonteCarloSpec) MonteCarloSpec
@doc (@doc SimulationCampaigns.MonteCarloSampleResult) MonteCarloSampleResult
@doc (@doc SimulationCampaigns.MonteCarloResult) MonteCarloResult
@doc (@doc SimulationCampaigns.run_monte_carlo) run_monte_carlo
@doc (@doc SimulationCampaigns.run_constellation_ensemble) run_constellation_ensemble
@doc (@doc SimulationCampaigns.ConstellationEnsembleResult) ConstellationEnsembleResult
@doc (@doc SimulationCampaigns.read_constellation_results) read_constellation_results

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
@doc (@doc SimulationModel.AbstractTypes.AbstractActuator) AbstractActuator
@doc (@doc SimulationModel.AbstractTypes.AbstractSensor) AbstractSensor
@doc (@doc SimulationModel.Thruster) Thruster
@doc (@doc SimulationModel.ReactionWheel) ReactionWheel
@doc (@doc SimulationModel.Magnetorquer) Magnetorquer
@doc (@doc SimulationModel.LaserTerminal) LaserTerminal
@doc (@doc SimulationModel.PowerReceiver) PowerReceiver
@doc (@doc SimulationModel.RadioTerminal) RadioTerminal
@doc (@doc SimulationModel.PanelHinge) PanelHinge
@doc (@doc SimulationModel.ThrusterGimbal) ThrusterGimbal
@doc (@doc SimulationModel.thrust_direction) thrust_direction
@doc (@doc SimulationModel.ControlMomentGyro) ControlMomentGyro
@doc _flat_doc(@doc SimulationModel.cmg_momentum) cmg_momentum
@doc (@doc SimulationModel.cmg_torque_axis) cmg_torque_axis
@doc (@doc SimulationModel.cmg_angle_saturate) cmg_angle_saturate
@doc (@doc SimulationModel.CMGInnerGimbal) CMGInnerGimbal
@doc (@doc SimulationModel.cmg_torque_axes) cmg_torque_axes
@doc (@doc SimulationModel.cmg_steering) cmg_steering
@doc (@doc SimulationModel.ActuatorContext) ActuatorContext
@doc _flat_doc(@doc SimulationModel.actuator_wrench) actuator_wrench
@doc (@doc SimulationModel.actuator_limits) actuator_limits
@doc (@doc SimulationModel.actuator_saturate) actuator_saturate
@doc _flat_doc(@doc SimulationModel.actuator_consumption) actuator_consumption
@doc (@doc SimulationModel.actuator_fraction) actuator_fraction
@doc (@doc SimulationModel.actuator_unit_wrench) actuator_unit_wrench
@doc (@doc SimulationModel.actuation_mode) actuation_mode
@doc _flat_doc(@doc SimulationModel.actuator_levels) actuator_levels
@doc (@doc SimulationModel.allocate_actuators) allocate_actuators
@doc (@doc SimulationModel.add_actuator!) add_actuator!
@doc (@doc SimulationModel.mounted_actuators) mounted_actuators
@doc (@doc SimulationModel.actuator_context) actuator_context
@doc (@doc SimulationModel.allocation_matrix) allocation_matrix
@doc (@doc SimulationModel.mounted_wrench) mounted_wrench
@doc (@doc SimulationModel.six_axis_thrusters) six_axis_thrusters
@doc (@doc SimulationModel.propellant_mass_rate) propellant_mass_rate
@doc (@doc SimulationModel.ActuatorRequirements) ActuatorRequirements
@doc (@doc SimulationModel.actuator_requirements) actuator_requirements
@doc (@doc SimulationModel.applied_actuator_levels) applied_actuator_levels
@doc (@doc SimulationModel.mounted_power) mounted_power
@doc (@doc SimulationModel.SolarArray) SolarArray
@doc (@doc SimulationModel.Battery) Battery
@doc (@doc SimulationModel.PowerLoad) PowerLoad
@doc (@doc SimulationModel.add_power!) add_power!
@doc (@doc SimulationModel.AbstractDataSource) AbstractDataSource
@doc (@doc SimulationModel.DataStore) DataStore
@doc (@doc SimulationModel.DataSource) DataSource
@doc (@doc SimulationModel.add_data!) add_data!
@doc (@doc SimulationModel.data_production_bps) data_production_bps
@doc (@doc SimulationModel.ActuatorFault) ActuatorFault
@doc (@doc SimulationModel.add_fault!) add_fault!
@doc (@doc SimulationModel.RateGyro) RateGyro
@doc (@doc SimulationModel.GNSSReceiver) GNSSReceiver
@doc (@doc SimulationModel.SunSensor) SunSensor
@doc (@doc SimulationModel.Accelerometer) Accelerometer
@doc (@doc SimulationModel.sensor_needs_specific_force) sensor_needs_specific_force
@doc (@doc SimulationModel.is_gravitational) is_gravitational
@doc (@doc SimulationModel.StarTracker) StarTracker
@doc (@doc SimulationModel.Magnetometer) Magnetometer
@doc (@doc SimulationModel.SensorTruth) SensorTruth
@doc (@doc SimulationModel.sensor_measure) sensor_measure
@doc (@doc SimulationModel.sensor_dimension) sensor_dimension
@doc (@doc SimulationModel.sensor_requirements) sensor_requirements
@doc (@doc SimulationModel.sensor_period_s) sensor_period_s
@doc (@doc SimulationModel.initialize_sensor!) initialize_sensor!
@doc (@doc SimulationModel.add_sensor!) add_sensor!
@doc (@doc SimulationModel.mounted_sensors) mounted_sensors
@doc (@doc SimulationModel.sensor_measurement) sensor_measurement
@doc (@doc SimulationModel.sensor_measurement_time) sensor_measurement_time
@doc (@doc SimulationModel.AbstractTypes.AbstractGuidanceModel) AbstractGuidanceModel
@doc (@doc SimulationModel.AbstractTypes.AbstractNavigationModel) AbstractNavigationModel
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
@doc (@doc SimulationModel.RobotArmPlanBuffer) RobotArmPlanBuffer
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
@doc (@doc SimulationModel.CompliantAttachment) CompliantAttachment
@doc (@doc SimulationModel.rectangular_prism_inertia) rectangular_prism_inertia
@doc (@doc SimulationModel.thin_panel_inertia) thin_panel_inertia
@doc (@doc SimulationModel.build_compliant_topology) build_compliant_topology
@doc (@doc SimulationModel.build_rectangular_compliant_grid) build_rectangular_compliant_grid
@doc (@doc SimulationModel.ConstellationSlot) ConstellationSlot
@doc (@doc SimulationModel.OrbitalPlane) OrbitalPlane
@doc (@doc SimulationModel.build_constellation) build_constellation
@doc (@doc SimulationModel.walker_planes) walker_planes
@doc (@doc SimulationModel.walker_constellation) walker_constellation
@doc (@doc SimulationModel.scoped) scoped
@doc (@doc SimulationModel.ScopedEffector) ScopedEffector
@doc (@doc SimulationModel.SpacecraftSelection) SpacecraftSelection
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
@doc (@doc SimulationModel.GRAMNearSurfaceAtmosphereModel) GRAMNearSurfaceAtmosphereModel
@doc (@doc SimulationModel.CombinedAtmosphereModel) CombinedAtmosphereModel
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
@doc """
    SpiceEphemeridesModel()

Use SPICE kernels for planetary frames, the Sun direction and the simulation
start epoch. Furnish kernels covering the chosen body and epoch before running
(for example through a SPICE-backed planet constructor). This constructor selects
the existing backend; it does not download or load kernels itself.
""" SpiceEphemeridesModel
@doc (@doc SimulationModel.make_no_gram_planet) make_no_gram_planet
@doc (@doc SimulationModel.make_no_gram_density_model) make_no_gram_density_model
@doc (@doc SimulationModel.make_no_gram_environment) make_no_gram_environment
@doc (@doc SimulationModel.wrench) wrench
@doc (@doc SimulationModel.environment_requirements) environment_requirements
@doc (@doc SimulationModel.solver_partition) solver_partition
@doc (@doc SimulationModel.gravity_backbone_structure) gravity_backbone_structure
@doc (@doc SimulationModel.gravity_backbone_acceleration_ii) gravity_backbone_acceleration_ii
@doc (@doc SimulationModel.gravity_backbone_kick_structure) gravity_backbone_kick_structure
@doc (@doc SimulationModel.gravity_backbone_kick_acceleration_ii) gravity_backbone_kick_acceleration_ii
# Copy canonical hook text into fresh documentation metadata for the root binding.
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
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.EphemeridesModels.ephemerides_frame_lpi)).text) ephemerides_frame_lpi
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.EphemeridesModels.ephemerides_sun_ii)).text) ephemerides_sun_ii
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.EphemeridesModels.ephemerides_epoch_seconds)).text) ephemerides_epoch_seconds
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.EphemeridesModels.ephemerides_requires_spice)).text) ephemerides_requires_spice
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.EphemeridesModels.ephemerides_cache_key)).text) ephemerides_cache_key
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.ControlHooks.update_control!)).text) update_control!
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.ControlHooks.control_wrench)).text) control_wrench
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.ControlHooks.control_propellant_rate)).text) control_propellant_rate
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.ControlHooks.control_wheel_torque)).text) control_wheel_torque
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.ControlHooks.control_thruster_levels(::Any, ::Int))).text) control_thruster_levels
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.ControlHooks.touchdown_spec(::Any, ::Int))).text) touchdown_spec
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.GuidanceHooks.update_guidance!)).text) update_guidance!
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.NavigationHooks.update_navigation!)).text) update_navigation!
@doc Base.Docs.docstr((Base.Docs.@ref(SimulationModel.VehicleThermalModels.getHeatRate)).text) getHeatRate

# 3.4. RPO Station Assets
@doc (@doc RPOStationAssets.station_geometry_path) station_geometry_path
@doc (@doc RPOStationAssets.station_cad_path) station_cad_path
@doc (@doc RPOStationAssets.load_rpo_station_pointcloud) load_rpo_station_pointcloud
@doc (@doc RPOStationAssets.load_rpo_station_cad_triangles) load_rpo_station_cad_triangles
@doc (@doc RPOStationAssets.load_rpo_station_cad_pointcloud) load_rpo_station_cad_pointcloud
@doc _flat_doc(@doc SimulationModel.NavigationHooks.RPOStationGeometry) RPOStationGeometry
@doc _flat_doc(@doc SimulationModel.NavigationHooks.RPOCubeSatGeometry) RPOCubeSatGeometry

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

# 3.6. Parallel Process
@doc (@doc ParallelProcess.ProcessPool) ProcessPool
@doc (@doc ParallelProcess.campaign_process_pool) campaign_process_pool
@doc (@doc ParallelProcess.ensure_process_workers!) ensure_process_workers!
@doc (@doc ParallelProcess.shutdown_process_pool!) shutdown_process_pool!
@doc (@doc ParallelProcess.adopt_process_workers!) adopt_process_workers!
@doc (@doc ParallelProcess.PROCESS_WARMUP) PROCESS_WARMUP

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
export ProcessPool, campaign_process_pool, ensure_process_workers!, shutdown_process_pool!, adopt_process_workers!, PROCESS_WARMUP
export ParallelConfig, SolverConfig, RuntimePolicyConfig, ArtifactConfig, SimulationEngineConfig
export simulation_engine_config_from_env
export prewarm_nbody_ephemeris_cache, load_nbody_ephemeris_cache!
export MonteCarloSpec, MonteCarloSampleResult, MonteCarloResult, run_monte_carlo
export run_constellation_ensemble, ConstellationEnsembleResult, read_constellation_results
export StateAnchor, get_state_anchor_callback
export AbstractForceTorqueModel, AbstractPlanet, AbstractDensityModel, AbstractControlEffectorModel
export AbstractEphemeridesModel, AbstractThermalModel, AbstractThrusterModel, AbstractGuidanceModel
export AbstractNavigationModel
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
export CompliantJointActuator, CompliantJointLoad, CompliantAttachment
export rectangular_prism_inertia, thin_panel_inertia
export build_compliant_topology, build_rectangular_compliant_grid
export ConstellationSlot, OrbitalPlane, build_constellation, walker_planes, walker_constellation
export ScopedEffector, scoped, SpacecraftSelection
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
export GRAMGridAtmosphereModel, GRAMNearSurfaceAtmosphereModel, CombinedAtmosphereModel
export SurrogatePresetResolution, available_surrogate_presets, resolve_surrogate_preset, surrogate_preset_model, atmosphere_provenance
export NRLMSISE00AtmosphereModel, init_nrlmsise_space_indices!
export SimpleEphemeridesModel, SpiceEphemeridesModel
export make_no_gram_planet, make_no_gram_density_model, make_no_gram_environment
export wrench, environment_requirements, solver_partition
export gravity_backbone_structure, gravity_backbone_acceleration_ii
export gravity_backbone_kick_structure, gravity_backbone_kick_acceleration_ii
export PlumeSurfaceConfig, PlumeSurfaceInteractionModel, PlumeSurfaceState, plume_surface_footprint, plume_erosion_onset_height, plume_ground_effect_force, plume_quantities
export DescentPhaseTargets, apollo11_descent_targets, ApolloDescentConfig, ApolloDescentState, ApolloDescentGuidanceModel
export ApolloDescentControlConfig, ApolloDescentControlModel, descent_attitude_command
export getDensity, getDensityBatch!, with_density_model_epoch
export ephemerides_frame_lpi, ephemerides_sun_ii, ephemerides_epoch_seconds
export ephemerides_requires_spice, ephemerides_cache_key
export control_wrench, control_propellant_rate, control_wheel_torque
export control_thruster_levels, touchdown_spec, bind_spacecraft
export getHeatRate
export update_control!, update_guidance!, update_navigation!
export AbstractActuator, Thruster, ReactionWheel, Magnetorquer, LaserTerminal, PowerReceiver, RadioTerminal, PanelHinge, ActuatorContext
export ThrusterGimbal, thrust_direction
export ControlMomentGyro, cmg_momentum, cmg_torque_axis, cmg_angle_saturate, cmg_steering
export CMGInnerGimbal, cmg_torque_axes
export actuator_wrench, actuator_limits, actuator_saturate, actuator_consumption, actuator_fraction, actuator_unit_wrench, actuation_mode, actuator_levels, allocate_actuators
export add_actuator!, mounted_actuators, actuator_context, allocation_matrix, mounted_wrench, six_axis_thrusters
export propellant_mass_rate, ActuatorRequirements, actuator_requirements, applied_actuator_levels, mounted_power
export SolarArray, Battery, PowerLoad, add_power!
export AbstractDataSource, DataStore, DataSource, add_data!, data_production_bps
export ActuatorFault, add_fault!
export AbstractSensor, RateGyro, StarTracker, Magnetometer, SensorTruth, sensor_measure, sensor_dimension, sensor_requirements, sensor_period_s, initialize_sensor!
export GNSSReceiver, SunSensor, Accelerometer, sensor_needs_specific_force, is_gravitational
export add_sensor!, mounted_sensors, sensor_measurement, sensor_measurement_time
export AerobrakingEnergyDepletionConfig, AerobrakingEnergyDepletionState
export AerobrakingEnergyDepletionGuidanceModel, AerobrakingEnergyDepletionControlModel
export SolarPanelAngleOfAttackControlModel
export ApoapsisTargetPeriapsisRaiseGuidanceModel
export VerificationRequest, VerificationResult
export run_verification, run_verification_cli, run_study, run_simulation, SimulationResults, link_schedule, run_table, long_table
export ControlForecastSnapshot, capture_control_forecast, run_control_forecast, with_control_forecast_worker
export station_geometry_path, station_cad_path, load_rpo_station_pointcloud, load_rpo_station_cad_triangles, load_rpo_station_cad_pointcloud
export RPOStationGeometry, RPOCubeSatGeometry
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

## 6. Runtime Initialization
# Runtime wiring that must not be baked into the precompiled image: these Refs
# hold closures over EnvironmentModels functions, so assigning them at include
# time would serialize a closure from an earlier world age. __init__ runs on
# every load of the cached image, which is what this needs.
function __init__()
	# Nothing measured about the precompiling host may stand in for this one
	# (see `_reset_process_local_state!` in precompile_workload.jl).
	_reset_process_local_state!()
	try
		SimulationModel.SimulationCallbacks._install_density_service_hooks!()
	catch err
		@warn "Could not install distributed density service hooks; the service will be unavailable." exception=(err, catch_backtrace())
	end
	return nothing
end



# Core-owned trial RHS contract for optional coordinate solvers.
const trial_rhs_safe = SimulationEngine.trial_rhs_safe
@doc (@doc SimulationEngine.trial_rhs_safe) trial_rhs_safe
export trial_rhs_safe
end # module SpaceAGORA
