# Canonical aggregator: no behavior ownership.
# Typed EDG owns prediction and decisions; public model types retain their owners.
module EDGAlgorithms
using ..GuidanceHooks: AerobrakingEnergyDepletionConfig, AerobrakingEnergyDepletionState
using ..ConfigTypes: ODEParams
using ..CommandTypes: AerobrakingControlCommand
using ..AerodynamicEffectors: aerodynamic_coefficient_fM
using ..GravityEffectors: aerobraking_gravity_force_ii
using ..FrameTransforms: r_intor_p!, rtolatlong, latlongtoNED, rvtoorbitalelement
using LinearAlgebra, StaticArrays, Roots, SpecialFunctions
using ..EDGServices: _edg_control_sat_state, _edg_control_pos_vel_mass, _edg_environment_state, _edg_in_drag_passage, _edg_ephemeris_time, _edg_planet_frame_lpi, _edg_targeting_prediction_environment, _edg_sample_prediction_atmosphere

include(joinpath(@__DIR__, "heat_rate.jl"))
include(joinpath(@__DIR__, "heat_load.jl"))
include(joinpath(@__DIR__, "structural_load.jl"))
include(joinpath(@__DIR__, "targeting.jl"))
include(joinpath(@__DIR__, "guidance_decision.jl"))
include(joinpath(@__DIR__, "angle_decision.jl"))

end # module EDGAlgorithms
