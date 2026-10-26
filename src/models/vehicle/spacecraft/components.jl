module Components

using StaticArrays
using LinearAlgebra
using ..AbstractTypes: AbstractActuator, AbstractSensor
using ..QuaternionMath: rot, dcm_to_quaternion
using Random: randn
using ..Sensors: SensorTruth
import ..Sensors: sensor_measure, sensor_dimension, sensor_requirements, sensor_needs_specific_force
using ..ThrusterModels: SixAxisThrusterModel
using ..Actuators: ActuatorContext, ActuatorRequirements, STANDARD_GRAVITY_MPS2, actuator_saturate, propellant_mass_rate
import ..Actuators: actuator_wrench, actuator_limits, actuator_consumption, actuator_requirements
export Facet, Thruster, Magnet, ReactionWheelAssembly, create_facet_list
export six_axis_thrusters, ReactionWheel, wheel_momentum_saturate, wheel_power, Magnetorquer, LaserTerminal, PowerReceiver, RadioTerminal, PanelHinge, hinge_step
export ThrusterGimbal, gimbal_step, thrust_direction
export ControlMomentGyro, cmg_momentum, cmg_torque_axis, cmg_angle_saturate, CMGInnerGimbal, cmg_torque_axes
export SolarArray, Battery, PowerLoad
export AbstractDataSource, DataStore, DataSource
export ActuatorFault, FaultSet
export RateGyro, StarTracker, Magnetometer, GNSSReceiver, SunSensor, Accelerometer

const I3 = SMatrix{3, 3, Float64}(diagm(ones(3)))

# Actuator hardware records live beside the actuator layer.
include(joinpath(@__DIR__, "..", "actuators", "thruster.jl"))
# Control moment gyros, held by the reaction-wheel assembly below.
include(joinpath(@__DIR__, "..", "actuators", "cmg.jl"))
include(joinpath(@__DIR__, "..", "actuators", "reaction_wheel.jl"))
include(joinpath(@__DIR__, "..", "actuators", "magnetorquer.jl"))
include(joinpath(@__DIR__, "..", "actuators", "laser.jl"))
# The receiving end of a power beam (an interlink terminal; see interlinks/power_links.jl).
include(joinpath(@__DIR__, "..", "actuators", "power_receiver.jl"))
# A radio, either end of a data link (an interlink terminal; see interlinks/data_links.jl).
include(joinpath(@__DIR__, "..", "actuators", "radio_terminal.jl"))
include(joinpath(@__DIR__, "..", "actuators", "hinge.jl"))
include(joinpath(@__DIR__, "..", "actuators", "gimbal.jl"))
# Power-system records (mounted with PowerSystems.add_power!).
include(joinpath(@__DIR__, "..", "power", "power_components.jl"))
# Onboard-data records (mounted with DataStorage.add_data!).
include(joinpath(@__DIR__, "..", "data", "data_components.jl"))
# The actuator-fault record (mounted with ActuatorFaults.add_fault!).
include(joinpath(@__DIR__, "..", "actuators", "fault.jl"))
# The built-in sensor kinds (mounted with SensorMounting.add_sensor!).
include(joinpath(@__DIR__, "..", "sensors", "sensor_kinds.jl"))

@kwdef mutable struct Facet
    area::Float64 = 0.0 # Area of the facet
    attitude::MVector{4, Float64} = @MVector [0.0, 0.0, 0.0, 1.0] # Attitude of the facet relative to the body frame
    normal_vector::MVector{3, Float64} = @MVector [1.0, 0.0, 0.0] # Normal vector in the facet body frame, for flat plate this is in the x-direction
    cp::MVector{3, Float64} = @MVector [0.0, 0.0, 0.0] # Center of pressure with respect to the center of mass of the Link containing the facet, in the Link frame
    ρ::Float64 = 0.0 # Diffuse coefficient
    δ::Float64 = 0.0 # Specular coefficient
    name::String = "" # Identifying name
end

function create_facet_list(area_list::Vector{Float64}, attitude_list::Vector{SVector{4, Float64}}, normal_vector_list::Vector{SVector{3, Float64}},
    cp_loc_list::Vector{SVector{3, Float64}},
    diffuse_coeffs_list::Vector{Float64}, specular_coeffs_list::Vector{Float64}, facet_names_list::Vector{String})
    list_length = length(area_list)
    lists = Vector[area_list, attitude_list, normal_vector_list, cp_loc_list, specular_coeffs_list, diffuse_coeffs_list]
    @assert all(l -> length(l) == list_length, lists) "Not all lists are the same length"
    facet_vector = Vector{Facet}(undef, list_length)
    for i in eachindex(area_list)
        facet_vector[i] = Facet(area_list[i], attitude_list[i], normal_vector_list[i], cp_loc_list[i], diffuse_coeffs_list[i], specular_coeffs_list[i], facet_names_list[i])
    end
    return facet_vector
end

end # module Components