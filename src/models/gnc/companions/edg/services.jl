"""
Versioned value and prediction services for the optional experimental EDG package,
contract 3.0.0. The provider implements the two decision services with the EDG
run state first:

- `guidance_decision!(state, config, u, p, t, spacecraft_index)`;
- `control_decision!(state, config, controlled_panel_links, u, p, t, spacecraft_index)`.

`state` is the run's `AerobrakingEnergyDepletionState`, which the provider
updates, `config` its `AerobrakingEnergyDepletionConfig`, `u` a copied
observation, `p` a borrowed prediction context, `t` elapsed seconds and
`spacecraft_index` the spacecraft's position in the run, from 1. The queries take
the context first: `_edg_environment_state(p, u, t, spacecraft_index)` and
`_edg_targeting_control_clock(p)`, which returns the control clock the engine
reported to that state when the run started (and again at each solve's start).
"""
module EDGServices
using ..SimulationModel.EDGServices: CONTRACT_VERSION, CompatibilityError, EDGUnavailableError,
    require_version, require_provider, available,
    AbstractPredictionContext, EDGControlDecision, planet_constants, thermal_parameters,
    entry_altitude_m, spacecraft_geometry, density_scale_height, frame_state, gravity_acceleration,
    _edg_control_sat_state, _edg_control_pos_vel_mass
# Only the capability signatures cross the public facade. Raw engine overloads
# remain private, even though they implement the same physical queries.
import ..SimulationModel.EDGServices: _edg_environment_state as _environment,
    _edg_in_drag_passage as _in_passage, _edg_ephemeris_time as _epoch,
    _edg_targeting_prediction_environment as _trajectory, _edg_targeting_control_clock as _clock
_edg_environment_state(p::AbstractPredictionContext, u, t::Float64, spacecraft_index::Int) =
    _environment(p, u, t, spacecraft_index)
_edg_in_drag_passage(p::AbstractPredictionContext, env) = _in_passage(p, env)
_edg_ephemeris_time(p::AbstractPredictionContext, t::Float64) = _epoch(p, t)
_edg_targeting_prediction_environment(p::AbstractPredictionContext, r, v, t::Float64) = _trajectory(p, r, v, t)
_edg_targeting_control_clock(p::AbstractPredictionContext) = _clock(p)
using ..SimulationModel.GuidanceHooks: AerobrakingEnergyDepletionConfig, AerobrakingEnergyDepletionState
using ..SimulationModel.CommandTypes: AerobrakingControlCommand
using ..SimulationModel.AerodynamicEffectors: aerodynamic_coefficient_fM
using ..SimulationModel.FrameTransforms: rvtoorbitalelement, rtolatlong
using ..SimulationModel.EDGAlgorithms: guidance_decision!, control_decision!
const IMPLEMENTED_FUNCTIONS = (:guidance_decision!, :control_decision!)
const CONSUMED_SERVICES = (
    :AbstractPredictionContext,
    :EDGControlDecision,
    :AerobrakingEnergyDepletionConfig,
    :AerobrakingEnergyDepletionState,
    :AerobrakingControlCommand,
    :planet_constants,
    :thermal_parameters,
    :entry_altitude_m,
    :spacecraft_geometry,
    :density_scale_height,
    :frame_state,
    :gravity_acceleration,
    :_edg_control_sat_state,
    :_edg_control_pos_vel_mass,
    :_edg_environment_state,
    :_edg_in_drag_passage,
    :_edg_ephemeris_time,
    :_edg_targeting_prediction_environment,
    :_edg_targeting_control_clock,
    :aerodynamic_coefficient_fM,
    :rvtoorbitalelement,
    :rtolatlong,
)
end # module EDGServices
