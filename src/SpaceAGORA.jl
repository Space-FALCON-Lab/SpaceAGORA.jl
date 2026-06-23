__precompile__(true)

module SpaceAGORA

using PrecompileTools: @compile_workload, @setup_workload

include(joinpath(@__DIR__, "parallel", "routing", "parallel_profiles.jl"))
include(joinpath(@__DIR__, "simulation", "runtime_services.jl"))
include(joinpath(@__DIR__, "core", "simulation_model.jl"))
include(joinpath(@__DIR__, "simulation", "engine", "simulation_engine.jl"))
include(joinpath(@__DIR__, "simulation", "campaigns", "simulation_campaigns.jl"))
include(joinpath(@__DIR__, "analysis", "verification", "telemetry_verification.jl"))
include(joinpath(@__DIR__, "assets", "rpo_station_assets.jl"))
include(joinpath(@__DIR__, "analysis", "visualization", "rpo", "rpo_visualization.jl"))
include(joinpath(@__DIR__, "cli", "spaceagora_cli.jl"))
include(joinpath(@__DIR__, "dynamics", "keplerian", "keplerian_propagator.jl"))
include(joinpath(@__DIR__, "mission", "constellation_design", "constellation_design.jl"))

using .ParallelProfiles: ParallelProfile, ParallelProfileConfig
using .ParallelProfiles: parse_parallel_profile, parallel_profile_name, profile_config, profile_env_pairs, with_parallel_profile
using .ParallelProfiles: OuterRouteFeatures, OuterRouteTuning, OuterRouteState
using .ParallelProfiles: reset_outer_route_state!, outer_route_signature, outer_route_stats_snapshot
using .ParallelProfiles: default_outer_route, outer_route_candidates, select_outer_route!, record_outer_route_feedback!
using .SimulationEngine: ParallelConfig, SolverConfig, RuntimePolicyConfig, ArtifactConfig, SimulationEngineConfig
using .SimulationEngine: simulation_engine_config_from_env
import .SimulationEngine: prewarm_nbody_ephemeris_cache, load_nbody_ephemeris_cache!
using .SimulationCampaigns: MonteCarloSpec, MonteCarloSampleResult, MonteCarloResult, run_monte_carlo
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
using .SimulationModel: NRLMSISE00AtmosphereModel, init_nrlmsise_space_indices!
using .SimulationModel: SimpleEphemeridesModel
using .SimulationModel: make_no_gram_planet, make_no_gram_density_model, make_no_gram_environment
using .SimulationModel: calcForceTorque, wrench, environment_requirements, solver_partition
using .SimulationModel: gravity_backbone_structure, gravity_backbone_acceleration_ii
using .SimulationModel: gravity_backbone_kick_structure, gravity_backbone_kick_acceleration_ii
using .SimulationModel: getDensity, getDensityBatch!
using .SimulationModel: calcControlEffect!, calcControlForceTorque, calcControlMassFlowRate
using .SimulationModel: ApoapsisTargetPeriapsisRaiseGuidanceModel
using .TelemetryVerification: VerificationRequest, VerificationResult
using .TelemetryVerification: run_verification, run_verification_cli, run_study
using .RPOStationAssets: station_geometry_path, station_cad_path, load_rpo_station_pointcloud, load_rpo_station_cad_triangles, load_rpo_station_cad_pointcloud
using .RPOVisualization: rpo_path_plot, rpo_tracking_plot
using .SpaceAGORACLI: AssetCheckItem, AssetCheckReport
using .KeplerianPropagator: oe2cart, kepler_satellite_state, kepler_prop_state, kepler_client_state
using .ConstellationDesign: run_stage0_seeding, run_constellation_design, run_stage2_verification, run_capo_pipeline

@doc (@doc SimulationEngine.ParallelConfig) ParallelConfig
@doc (@doc SimulationEngine.SolverConfig) SolverConfig
@doc (@doc SimulationEngine.RuntimePolicyConfig) RuntimePolicyConfig
@doc (@doc SimulationEngine.ArtifactConfig) ArtifactConfig
@doc (@doc SimulationEngine.SimulationEngineConfig) SimulationEngineConfig
@doc (@doc SimulationEngine.simulation_engine_config_from_env) simulation_engine_config_from_env
@doc (@doc SimulationCampaigns.MonteCarloSpec) MonteCarloSpec
@doc (@doc SimulationCampaigns.MonteCarloSampleResult) MonteCarloSampleResult
@doc (@doc SimulationCampaigns.MonteCarloResult) MonteCarloResult
@doc (@doc SimulationCampaigns.run_monte_carlo) run_monte_carlo
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
@doc (@doc SimulationModel.RobotArmPlan) RobotArmPlan

"""
    NoAtmosphereModel()

Density-model constructor for no-atmosphere baseline runs and no-GRAM onboarding
scenarios.
"""
NoAtmosphereModel

"""
    ExponentialAtmosphereModel(planet)
    ExponentialAtmosphereModel(rho_ref, h_ref, H; temperature_k=200.0, valid_min_altitude_m=h_ref, valid_max_altitude_m=h_ref + 5H)

Single-scale-height analytic density-model constructor for baseline runs that
should not depend on GRAM assets. The `planet` convenience form uses the
built-in reference density and scale-height constants for the chosen body. This
model assumes zero winds and constant temperature, and its validity-range
keywords are advisory only; evaluation extrapolates the same exponential
outside the documented band.
"""
ExponentialAtmosphereModel

"""
    PiecewiseExponentialAtmosphereModel(h_breaks_m, rho_refs, Hs; h_refs=h_breaks_m[1:end-1], temperature_k=200.0, valid_min_altitude_m=first(h_breaks_m), valid_max_altitude_m=last(h_breaks_m))

Multi-layer analytic density-model constructor for bounded no-GRAM studies that
need more shape than a single scale height. The model assumes zero winds and
constant temperature and uses the nearest configured layer to extrapolate
outside the advisory validity band.
"""
PiecewiseExponentialAtmosphereModel

"""
    NRLMSISE00AtmosphereModel(; f107a=150.0, f107=150.0, ap=4.0, index_provider=nothing, use_space_indices=false, space_indices_force_download=false, include_anomalous_oxygen=true, valid_min_altitude_m=0.0, valid_max_altitude_m=1000e3)

NRLMSISE-00 atmosphere-model constructor for runs that need empirical
thermospheric density without GRAM. Use fixed `f107a`, `f107`, and `ap` values
or provide `index_provider`, a callable that returns `(f107a, f107, ap)` for a
requested instant. Set `use_space_indices=true` to use the built-in
CelesTrak-backed F10.7/Ap dataset path, optionally prewarmed through
[`init_nrlmsise_space_indices!`](@ref). The standard NRLMSISE-00 validity band
is approximately `0 m` to `1000 km`; the validity fields document that range
but do not clamp evaluation.
"""
NRLMSISE00AtmosphereModel

"""
    init_nrlmsise_space_indices!(; force_download=false)

Initialize or refresh the CelesTrak space-weather dataset used by
`NRLMSISE00AtmosphereModel(use_space_indices=true)`.

Call this before long runs if you want any dataset download or refresh to
happen before the first atmosphere evaluation.
"""
init_nrlmsise_space_indices!

"""
    SimpleEphemeridesModel(; reference_epoch_seconds=0.0, prime_meridian_at_reference_rad=0.0)

Analytic ephemerides/frame backend for onboarding and open-data runs that should
not depend on local SPICE kernels.
"""
SimpleEphemeridesModel

@doc (@doc SimulationModel.make_no_gram_planet) make_no_gram_planet
@doc (@doc SimulationModel.make_no_gram_density_model) make_no_gram_density_model
@doc (@doc SimulationModel.make_no_gram_environment) make_no_gram_environment

@doc (@doc KeplerianPropagator.oe2cart) oe2cart
@doc (@doc KeplerianPropagator.kepler_satellite_state) kepler_satellite_state
@doc (@doc KeplerianPropagator.kepler_prop_state) kepler_prop_state
@doc (@doc KeplerianPropagator.kepler_client_state) kepler_client_state

"""
    calcForceTorque(model, x, p, i) -> (force_n, torque_n_m)

Stable extension hook for custom [`AbstractForceTorqueModel`](@ref)
implementations. Extend this method for package or user models that contribute
translational and rotational wrench terms to the simulation RHS.
"""
calcForceTorque

"""
    wrench(model, x::StateSample, env::EnvironmentSample, t::Float64) -> (force_ii, torque_body)

Preferred additive extension hook for custom [`AbstractForceTorqueModel`](@ref)
implementations. The engine owns stage-consistent sampling and caching, then
passes a typed state/environment bundle into `wrench`.

Return inertial-frame force and body-frame torque in SI units. Implementations
should behave as pure functions of `(model, x, env, t)`.
"""
wrench

"""
    environment_requirements(model) -> EffectorEnvironmentRequirements

Preferred additive declaration hook for the sampled environment capabilities a
[`wrench`](@ref) implementation requires. The default requests no sampled
environment fields.
"""
environment_requirements

"""
    solver_partition(model) -> Symbol

Optional additive declaration hook for `split_imex` solver partitioning of
dynamic effectors.

Return `:implicit` to place the effector on the atmosphere-implicit IMEX side,
or `:explicit` to keep it on the non-stiff explicit side. The default is
`:explicit`.
"""
solver_partition

"""
    gravity_backbone_structure(model) -> Symbol

Optional additive declaration hook for the `gravity_backbone_split` solver
mode.

Return `:position_only_static_gravity` for effectors that can participate in
the gravity-only translational backbone, or `:unsupported` otherwise. The
default is `:unsupported`.
"""
gravity_backbone_structure

"""
    gravity_backbone_acceleration_ii(model, x::StateSample, env::EnvironmentSample, t::Float64) -> accel_ii

Optional additive acceleration hook for `gravity_backbone_split`.

Implementations must return inertial-frame translational acceleration in SI
units for effectors that declare
[`gravity_backbone_structure`](@ref) == `:position_only_static_gravity`.
"""
gravity_backbone_acceleration_ii

"""
    gravity_backbone_kick_structure(model) -> Symbol

Optional additive declaration hook for explicit translational perturbation kicks
in `gravity_backbone_split`.

Return `:velocity_kick_explicit` for effectors that should be applied as
explicit velocity kicks around the gravity core, or `:unsupported` otherwise.
The default is `:unsupported`.
"""
gravity_backbone_kick_structure

"""
    gravity_backbone_kick_acceleration_ii(model, x::StateSample, env::EnvironmentSample, t::Float64) -> accel_ii

Optional additive acceleration hook for explicit velocity kicks in
`gravity_backbone_split`.

Implementations must return inertial-frame translational acceleration in SI
units for effectors that declare
[`gravity_backbone_kick_structure`](@ref) == `:velocity_kick_explicit`.
"""
gravity_backbone_kick_acceleration_ii

"""
    getDensity(model, h, lat, lon, el_time, wind[, p]) -> (rho, temperature, wind_vec)

Stable extension hook for custom [`AbstractDensityModel`](@ref)
implementations. The scalar form returns density, temperature, and wind for a
single atmosphere query.
"""
getDensity

"""
    getDensityBatch!(rhos, Ts, winds, model, hs, lats, lons, el_time, wind, p)

Optional batch extension hook for [`AbstractDensityModel`](@ref)
implementations that can answer many atmosphere queries more efficiently than
repeated scalar `getDensity` dispatch.
"""
getDensityBatch!

"""
    calcControlEffect!(model, u, p, t, i)

Stable extension hook for [`AbstractControlEffectorModel`](@ref)
implementations that update control-related shared state during the simulation
loop.
"""
calcControlEffect!

"""
    calcControlForceTorque(model, u, p, i, t)

Stable extension hook for [`AbstractControlEffectorModel`](@ref)
implementations that contribute force and torque terms to the spacecraft
dynamics.
"""
calcControlForceTorque

"""
    calcControlMassFlowRate(model, u, p, i, t)

Stable extension hook for [`AbstractControlEffectorModel`](@ref)
implementations that consume propellant. Effectors that do not model propellant
consumption should return `0.0`.
"""
calcControlMassFlowRate

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

@doc (@doc TelemetryVerification.VerificationRequest) VerificationRequest
@doc (@doc TelemetryVerification.VerificationResult) VerificationResult
@doc (@doc TelemetryVerification.run_verification) run_verification
@doc (@doc TelemetryVerification.run_verification_cli) run_verification_cli
@doc (@doc TelemetryVerification.run_study) run_study

@doc (@doc SpaceAGORACLI.AssetCheckItem) AssetCheckItem
@doc (@doc SpaceAGORACLI.AssetCheckReport) AssetCheckReport

export ParallelProfile, ParallelProfileConfig
export parse_parallel_profile, parallel_profile_name, profile_config, profile_env_pairs, with_parallel_profile
export OuterRouteFeatures, OuterRouteTuning, OuterRouteState
export reset_outer_route_state!, outer_route_signature, outer_route_stats_snapshot
export default_outer_route, outer_route_candidates, select_outer_route!, record_outer_route_feedback!
export ParallelConfig, SolverConfig, RuntimePolicyConfig, ArtifactConfig, SimulationEngineConfig
export simulation_engine_config_from_env
export prewarm_nbody_ephemeris_cache, load_nbody_ephemeris_cache!
export MonteCarloSpec, MonteCarloSampleResult, MonteCarloResult, run_monte_carlo
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
export NoAtmosphereModel, ExponentialAtmosphereModel, PiecewiseExponentialAtmosphereModel
export NRLMSISE00AtmosphereModel, init_nrlmsise_space_indices!
export SimpleEphemeridesModel
export make_no_gram_planet, make_no_gram_density_model, make_no_gram_environment
export calcForceTorque, wrench, environment_requirements, solver_partition
export gravity_backbone_structure, gravity_backbone_acceleration_ii
export gravity_backbone_kick_structure, gravity_backbone_kick_acceleration_ii
export getDensity, getDensityBatch!
export calcControlEffect!, calcControlForceTorque, calcControlMassFlowRate
export ApoapsisTargetPeriapsisRaiseGuidanceModel
export VerificationRequest, VerificationResult
export run_verification, run_verification_cli, run_study, run_simulation
export station_geometry_path, station_cad_path, load_rpo_station_pointcloud, load_rpo_station_cad_triangles, load_rpo_station_cad_pointcloud
export AssetCheckItem, AssetCheckReport, check_assets, render_asset_report, run_cli

"""
    run_simulation(args...; isolate_state=true, kwargs...)

Stable package entrypoint for simulation execution used by calibration integrations.

By default, `isolate_state=true` deep-copies the simulation configuration before execution.
This preserves correctness and reentrancy across repeated runs and concurrent callers by
preventing one run from mutating shared campaign or model state that another run still
references.

Set `isolate_state=false` only as an advanced performance lever when the caller owns the
configuration instance and will not reuse it concurrently or across runs that may mutate
shared state. This can reduce setup cost for large mission definitions or many short runs,
but it trades away the default isolation guarantee.

# Examples
```jldoctest
julia> args = spaceagora_no_gram_example_args();

julia> sol = run_simulation(args; return_solution=true);

julia> length(sol.t) > 1
true
```
"""
run_simulation(args...; kwargs...) = SimulationEngine.run_simulation(args...; kwargs...)
run_simulation(config::SimulationEngineConfig, args...; kwargs...) = SimulationEngine.run_simulation(config, args...; kwargs...)

"""
    prewarm_nbody_ephemeris_cache(args; dt_s=nothing, mission_end_s=nothing, save_path=nothing) -> cache
    prewarm_nbody_ephemeris_cache(config, args; dt_s=nothing, mission_end_s=nothing, save_path=nothing) -> cache

Precompute and register a process-local N-body SPICE ephemeris cache for later
[`run_simulation`](@ref) calls. This is intended for Monte Carlo campaigns that
reuse the same third-body set, start epoch, mission span, and cache sample
spacing across many runs. The returned cache is keyed by the same deterministic
boundary that the runtime setup already uses.

If `save_path` is provided, the cache is also serialized to disk so other Julia
worker processes can call [`load_nbody_ephemeris_cache!`](@ref) and reuse the
same precomputed ephemeris without rebuilding it from SPICE.
"""
prewarm_nbody_ephemeris_cache(args...; kwargs...) = SimulationEngine.prewarm_nbody_ephemeris_cache(args...; kwargs...)

"""
    load_nbody_ephemeris_cache!(path; replace=true) -> cache

Load a serialized N-body ephemeris cache created by
[`prewarm_nbody_ephemeris_cache`](@ref) and register it in the current Julia
process so later [`run_simulation`](@ref) calls can reuse it. This is intended
for multi-process Monte Carlo campaigns where each worker should load the same
precomputed SPICE cache once before running many trajectories.
"""
load_nbody_ephemeris_cache!(args...; kwargs...) = SimulationEngine.load_nbody_ephemeris_cache!(args...; kwargs...)

"""
    check_assets(; repo_root=pwd()) -> AssetCheckReport

Inspect the current repository asset layout and report which baseline, optional,
and high-fidelity asset roots are available.
"""
check_assets(args...; kwargs...) = SpaceAGORACLI.check_assets(args...; kwargs...)

"""
    render_asset_report(report; io=stdout)

Render a human-readable asset status report.
"""
render_asset_report(args...; kwargs...) = SpaceAGORACLI.render_asset_report(args...; kwargs...)

"""
    run_cli([args=ARGS]; io=stdout, errio=stderr) -> Int

Stable CLI entrypoint for SpaceAGORA operational commands:

- `run`
- `telemetry`
- `benchmark`
- `assets check`

This is the package-owned command surface used by the `bin/spaceagora` wrapper.
"""
run_cli(args...; kwargs...) = SpaceAGORACLI.run_cli(args...; kwargs...)

include(joinpath(@__DIR__, "precompile_workload.jl"))

end # module SpaceAGORA
