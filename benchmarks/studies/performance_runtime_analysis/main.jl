const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_OUTPUT_DIR = joinpath(REPO_ROOT, "output", "performance")
const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH_HARMONICS_FILE = joinpath(REPO_ROOT, "data/planets/gravity", "EarthGGM05C.csv")
const EARTH_GRAM_SURROGATE_FILE = joinpath(REPO_ROOT, "data", "planets", "gram_surrogate", "earth_surrogate.jls")
const PERF_BASELINE_SCENARIO = "single_j2"
const _PERF_POLICY_ENV_NAMES = (
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE",
    "SPACEAGORA_DENSITY_CALLBACK_PARALLEL",
    "SPACEAGORA_CONTROL_CALLBACK_PARALLEL",
    "SPACEAGORA_MULTIBODY_PARALLEL",
)
const _PERF_POLICY_ENV_BASELINE = Dict{String, Union{Nothing, String}}(
    name => (haskey(ENV, name) ? String(ENV[name]) : nothing)
    for name in _PERF_POLICY_ENV_NAMES
)
const _PERF_THREADS_BACKEND_WARNING_EMITTED = Ref(false)

# Use the prepared caller environment and the package's canonical module identities.
import SpaceAGORA
const ParallelProfiles = SpaceAGORA.ParallelProfiles
const RuntimeServices = SpaceAGORA.RuntimeServices
using .ParallelProfiles

using CSV
using DataFrames
using Arrow
using Dates
using Distributed
using LinearAlgebra
using Random
using Sockets
using SPICE
using StaticArrays
using Statistics

if myid() == 1
    ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
    using Plots
end

const SimulationModel = SpaceAGORA.SimulationModel
using .SimulationModel
const MissionType = SimulationModel.SimConfig.MissionType
const MissionTime = SimulationModel.SimConfig.MissionTime
const MissionOrbits = SimulationModel.SimConfig.MissionOrbits
const SimulationConfiguration = SimulationModel.SimConfig.SimulationConfiguration
const InitialTime = SimulationModel.SimConfig.InitialTime
const IntegrationTolerances = SimulationModel.SimConfig.IntegrationTolerances
const SimulationSettings = SimulationModel.SimConfig.SimulationSettings
const MissionConfiguration = SimulationModel.SimConfig.MissionConfiguration
const EnvironmentModel = SimulationModel.SimConfig.EnvironmentModel
const Earth = SimulationModel.Planets.Earth
const Mars = SimulationModel.Planets.Mars
const Link = SimulationModel.SpacecraftModels.Link
const Joint = SimulationModel.SpacecraftModels.Joint
const SpacecraftModel = SimulationModel.SpacecraftModels.SpacecraftModel
const DynamicsModel = SimulationModel.SpacecraftModels.DynamicsModel
const InitialCondition = SimulationModel.SpacecraftModels.InitialCondition
const GuidanceModel = SimulationModel.SpacecraftModels.GuidanceModel
const NavigationModel = SimulationModel.SpacecraftModels.NavigationModel
const ControlModel = SimulationModel.SpacecraftModels.ControlModel
const AbstractPlanet = SimulationModel.AbstractTypes.AbstractPlanet
const NoAtmosphereModel = SimulationModel.EnvironmentModels.NoAtmosphereModel
const GRAMAtmosphereModel = SimulationModel.EnvironmentModels.GRAMAtmosphereModel
const GRAMAtmosphereModelSurrogate = SimulationModel.EnvironmentModels.GRAMAtmosphereModelSurrogate
const MaxwellianHeat = SimulationModel.VehicleThermalModels.MaxwellianHeat
const InverseSquaredGravityModel = SimulationModel.DynamicEffectors.InverseSquaredGravityModel
const InverseSquaredJ2GravityModel = SimulationModel.DynamicEffectors.InverseSquaredJ2GravityModel
const GravitationalHarmonicsModel = SimulationModel.DynamicEffectors.GravitationalHarmonicsModel
const NBodyGravityModel = SimulationModel.DynamicEffectors.NBodyGravityModel
const SolarRadiationPressureModel = SimulationModel.DynamicEffectors.SolarRadiationPressureModel
const AerodynamicCoefficientfM = SimulationModel.DynamicEffectors.AerodynamicCoefficientfM
const BaseThrusterModel = SimulationModel.ThrusterModels.BaseThrusterModel
const AerobrakingCampaignPropulsiveManeuverGuidanceModel =
    SimulationModel.GuidanceModels.AerobrakingCampaignPropulsiveManeuverGuidanceModel

function _ensure_perf_gramsuite_loaded!()
    vendored_gramsuite = joinpath(REPO_ROOT, "data", "GRAMSuite.jl")
    RuntimeServices._prepare_gramsuite_load_path!(vendored_gramsuite; prepend=false)
    @eval import GRAMSuite
    return nothing
end

_ensure_perf_gramsuite_loaded!()

# GRAM constructors, copies, serialization, locks and native ownership belong
# to the package extension; this benchmark defines no competing adapters.
Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) === nothing &&
    error("Runtime benchmark requires the canonical SpaceAGORAGRAMSuiteExt extension.")
const SimulationEngine = SpaceAGORA.SimulationEngine
const run_simulation = SpaceAGORA.run_simulation


include(joinpath(@__DIR__, "reporting", "runtime_reports.jl"))
using .RuntimeReports
using .RuntimeReports: _parse_positive_int_env, _parse_nonnegative_int_env, _parse_nonnegative_float_env, _parse_bool_env,
    _perf_smoke_mode, _perf_case_name_filter, _exclude_entry_scenarios, _split_rollout_enabled,
    _split_rollout_enforce, _split_rollout_case_names, _split_rollout_solver_variants, _split_rollout_max_slowdown_ratio,
    _split_rollout_pos_rel_tol, _split_rollout_vel_rel_tol, _split_rollout_q_angle_tol_rad, _split_rollout_omega_rel_tol,
    _split_rollout_sample_count, _multirate_rollout_enabled, _multirate_rollout_enforce, _multirate_rollout_case_names,
    _multirate_rollout_max_slowdown_ratio, _multirate_rollout_pos_rel_tol, _multirate_rollout_vel_rel_tol, _multirate_rollout_q_angle_tol_rad,
    _multirate_rollout_omega_rel_tol, _multirate_rollout_sample_count, _multirate_env_float_or_missing, _multirate_env_int_or_missing,
    _multirate_rollout_setting_snapshot, _safe_unique_join, _fmt, _scenario_metric,
    _safe_stat, _ci95_bounds, _ci95_low, _ci95_high,
    _sem, _cv_pct, _sum_nonmissing, _safe_share,
    _summary_group_keys, _case_with_solver, _solution_state_at, _relative_vector_delta,
    _quaternion_angle_delta_rad, _trajectory_delta_metrics

include(joinpath(@__DIR__, "case_catalog.jl"))
include(joinpath(@__DIR__, "measurement.jl"))
include(joinpath(@__DIR__, "reporting.jl"))
include(joinpath(@__DIR__, "cli.jl"))
