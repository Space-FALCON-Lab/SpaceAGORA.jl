struct ConstantDensityModel <: SimulationModel.AbstractDensityModel
    rho::Float64
    temp::Float64
end

struct TimedTangentialThrusterModel <: SimulationModel.AbstractControlEffectorModel
    thrust::Float64
    direction_sign::Float64 # +1.0 prograde, -1.0 retrograde
    start_time::Float64
    stop_time::Float64
end

mutable struct CountingGuidanceModel <: SimulationModel.AbstractTypes.AbstractGuidanceModel
    hits::Vector{Int}
end

mutable struct CountingNavigationModel
    hits::Vector{Int}
end

mutable struct CountingControlModel <: SimulationModel.AbstractControlEffectorModel
    hits::Vector{Int}
end

struct ThrowingForceModel <: SimulationModel.AbstractForceTorqueModel
end

struct NaNForceModel <: SimulationModel.AbstractForceTorqueModel
end

struct NaNParamForceModel <: SimulationModel.AbstractForceTorqueModel
end

struct ConstantTorqueModel <: SimulationModel.AbstractForceTorqueModel
    torque::SVector{3, Float64}
end

struct ConstantForceModel <: SimulationModel.AbstractForceTorqueModel
    force::SVector{3, Float64}
end

struct WrenchOnlyForceModel <: SimulationModel.AbstractForceTorqueModel
    force::SVector{3, Float64}
    torque::SVector{3, Float64}
end

struct AtmosphereProbeWrenchModel <: SimulationModel.AbstractForceTorqueModel
end

struct ImplicitLegacyForceModel <: SimulationModel.AbstractForceTorqueModel
    force::SVector{3, Float64}
    torque::SVector{3, Float64}
end

struct InvalidPartitionForceModel <: SimulationModel.AbstractForceTorqueModel
end

struct BackboneCustomGravityModel <: SimulationModel.AbstractForceTorqueModel
    accel::SVector{3, Float64}
end

struct InvalidBackboneStructureModel <: SimulationModel.AbstractForceTorqueModel
end

struct SolarBackboneModel <: SimulationModel.AbstractForceTorqueModel
end

struct InvalidBackboneKickStructureModel <: SimulationModel.AbstractForceTorqueModel
end

struct PlanetFrameKickModel <: SimulationModel.AbstractForceTorqueModel
end

struct ThrowingOrbitPlanet <: SimulationModel.AbstractPlanet
    Rp_e::Float64
end

function SimulationModel.EnvironmentModels.getDensity(
    model::ConstantDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)
    return model.rho, model.temp, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcControlForceTorque(
    model::TimedTangentialThrusterModel,
    u::AbstractVector{Float64},
    p::ODEParams,
    i::Int64,
    t::Float64
)
    if t < model.start_time || t > model.stop_time
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    v = SVector{3, Float64}(u.vel)
    vm = norm(v)
    if vm == 0.0 || !isfinite(vm)
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    force = model.thrust * model.direction_sign * (v / vm)
    return force, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcControlEffect!(
    model::TimedTangentialThrusterModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    return nothing
end

function SimulationModel.calcGuidanceEffect!(
    model::CountingGuidanceModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    model.hits[i] += 1
    return nothing
end

function SimulationModel.calcNavigationEffect!(
    model::CountingNavigationModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    model.hits[i] += 1
    return nothing
end

function SimulationModel.calcControlEffect!(
    model::CountingControlModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    model.hits[i] += 1
    return nothing
end

function SimulationModel.calcForceTorque(
    model::ThrowingForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    error("intentional derivative failure")
end

function SimulationModel.calcForceTorque(
    model::NaNForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return SVector{3, Float64}(NaN, NaN, NaN), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcForceTorque(
    model::NaNParamForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    p.shared_buffers.current_time[] = NaN
    return SVector{3, Float64}(NaN, NaN, NaN), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcForceTorque(
    model::ConstantTorqueModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return SVector{3, Float64}(0.0, 0.0, 0.0), model.torque
end

function SimulationModel.calcForceTorque(
    model::ConstantForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return model.force, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.wrench(
    model::WrenchOnlyForceModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64
)
    return model.force, model.torque
end

SimulationModel.solver_partition(::ImplicitLegacyForceModel) = :implicit
SimulationModel.solver_partition(::InvalidPartitionForceModel) = :bad_partition
SimulationModel.gravity_backbone_structure(::BackboneCustomGravityModel) = :position_only_static_gravity
SimulationModel.gravity_backbone_structure(::InvalidBackboneStructureModel) = :bad_structure
SimulationModel.gravity_backbone_structure(::SolarBackboneModel) = :position_only_static_gravity
SimulationModel.gravity_backbone_kick_structure(::InvalidBackboneKickStructureModel) = :bad_structure
SimulationModel.gravity_backbone_kick_structure(::PlanetFrameKickModel) = :velocity_kick_explicit

function SimulationModel.gravity_backbone_acceleration_ii(
    model::BackboneCustomGravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64
)
    return model.accel
end

function SimulationModel.gravity_backbone_acceleration_ii(
    model::SolarBackboneModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64
)
    return SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.environment_requirements(::SolarBackboneModel)
    return EffectorEnvironmentRequirements(solar=true)
end

function SimulationModel.environment_requirements(::PlanetFrameKickModel)
    return EffectorEnvironmentRequirements(planet_frame=true)
end

function SimulationModel.gravity_backbone_kick_acceleration_ii(
    model::PlanetFrameKickModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64
)
    return SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcForceTorque(
    model::ImplicitLegacyForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return model.force, model.torque
end

function SimulationModel.calcForceTorque(
    model::InvalidPartitionForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.environment_requirements(::AtmosphereProbeWrenchModel)
    return EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
end

function SimulationModel.wrench(
    model::AtmosphereProbeWrenchModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64
)
    env.atmosphere === nothing && error("expected atmosphere sample")
    env.planet_frame === nothing && error("expected planet-frame sample")
    return SVector{3, Float64}(env.atmosphere.rho_kg_m3, env.planet_frame.alt_m, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.ControlHooks.rvtoorbitalelement(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet::ThrowingOrbitPlanet
)
    throw(ErrorException("forced-orbital-element-conversion-failure"))
end
