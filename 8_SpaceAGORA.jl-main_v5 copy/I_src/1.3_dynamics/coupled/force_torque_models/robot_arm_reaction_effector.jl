"""Robot-arm reaction effector module that couples joint-load feedback into spacecraft force and torque models."""
module RobotArmReactionEffectors

using StaticArrays
using ...AbstractTypes: AbstractForceTorqueModel
using ...ClothMultibody: CompliantJointActuator
using ...RobotArmPlanning: RobotArmPlan, robot_arm_plan_sample
import ..DynamicEffectors: calcForceTorque

export RobotArmReactionEffector

"""Force-torque effector that feeds robot-arm joint reaction loads back to the spacecraft base."""
Base.@kwdef mutable struct RobotArmReactionEffector <: AbstractForceTorqueModel
    spacecraft_idx::Int = 1
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

"""Compute the robot-arm reaction force and torque applied to the spacecraft."""
function calcForceTorque(model::RobotArmReactionEffector, sc_view, p, sat_idx::Int)
    sat_idx == model.spacecraft_idx || return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    model.plan === nothing && return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    # Placeholder hook for the coupled backend: the feedforward path is present,
    # but reaction wrench estimation belongs with Cloth/RNEA dynamics.
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

end # module RobotArmReactionEffectors
