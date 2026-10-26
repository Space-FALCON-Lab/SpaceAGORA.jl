"""Robot-arm reaction effector module that couples joint-load feedback into spacecraft force and torque models."""
module RobotArmReactionEffectors

using StaticArrays
using ...AbstractTypes: AbstractForceTorqueModel
using ...ClothMultibody: CompliantJointActuator
using ...RobotArmPlanning: RobotArmPlan, robot_arm_plan_sample
using ...EffectorSampling: StateSample, EnvironmentSample
import ..DynamicEffectors: wrench

export RobotArmReactionEffector

"""
Force-torque effector that feeds robot-arm joint reaction loads back to the base of
the spacecraft at run position `spacecraft_index`.
"""
Base.@kwdef mutable struct RobotArmReactionEffector <: AbstractForceTorqueModel
    spacecraft_index::Int = 1
    plan::Union{Nothing, RobotArmPlan} = nothing
    updated_at_s::Float64 = 0.0
    force_scale::Float64 = 0.0
    torque_scale::Float64 = 0.0
    k_translation_n_m::Any = 5.0e3
    c_translation_n_s_m::Any = 30.0
    k_rotation_n_m_rad::Any = 15.0
    c_rotation_n_m_s_rad::Any = 0.5
    joint_actuators::Vector{CompliantJointActuator} = CompliantJointActuator[]
end

# The robot-arm reaction force and torque on the spacecraft at `x.spacecraft_index`.
function wrench(model::RobotArmReactionEffector, x::StateSample, env::EnvironmentSample, t::Float64)
    zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    x.spacecraft_index == model.spacecraft_index || return zero3, zero3
    model.plan === nothing && return zero3, zero3
    # Placeholder hook for the coupled backend: the feedforward path is present,
    # but reaction wrench estimation belongs with compliant multibody/RNEA dynamics.
    return zero3, zero3
end


end # module RobotArmReactionEffectors
