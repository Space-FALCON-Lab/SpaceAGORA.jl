module EffectorSampling

using StaticArrays

export StateSample,
    PlanetFrameSample,
    AtmosphereSample,
    SolarEphemerisSample,
    ThirdBodyEphemerisSample,
    EnvironmentSample,
    EffectorEnvironmentRequirements,
    wrench,
    environment_requirements

"""
    StateSample

Typed state view passed to the additive [`wrench`](@ref) force/torque hook.

`pos_ii` and `vel_ii` are inertial-frame SI vectors, `mass_kg` is spacecraft
mass, `q_ib` is the inertial-to-body attitude quaternion when available, and
`ω_body` is the body-frame angular velocity when available. `spacecraft` is the
typed spacecraft model handle for effectors that need static geometry or
inertia data without reaching back into the integrator state.
"""
struct StateSample{S}
    pos_ii::SVector{3, Float64}
    vel_ii::SVector{3, Float64}
    mass_kg::Float64
    q_ib::Union{Nothing, SVector{4, Float64}}
    ω_body::Union{Nothing, SVector{3, Float64}}
    spacecraft::S
end

StateSample(
    pos_ii::SVector{3, Float64},
    vel_ii::SVector{3, Float64},
    mass_kg::Real;
    q_ib::Union{Nothing, AbstractVector{<:Real}}=nothing,
    ω_body::Union{Nothing, AbstractVector{<:Real}}=nothing,
    spacecraft=nothing,
) = StateSample(
    pos_ii,
    vel_ii,
    Float64(mass_kg),
    q_ib === nothing ? nothing : SVector{4, Float64}(q_ib),
    ω_body === nothing ? nothing : SVector{3, Float64}(ω_body),
    spacecraft,
)

"""
    PlanetFrameSample

Stage-consistent planet-relative kinematics for the current ODE evaluation.
"""
struct PlanetFrameSample
    l_pi::SMatrix{3, 3, Float64, 9}
    pos_pp::SVector{3, Float64}
    vel_pp::SVector{3, Float64}
    alt_m::Float64
    lat_rad::Float64
    lon_rad::Float64
end

"""
    AtmosphereSample

Stage-consistent atmospheric properties in SI units.
"""
struct AtmosphereSample
    rho_kg_m3::Float64
    temperature_k::Float64
    wind_pp::SVector{3, Float64}
end

"""
    SolarEphemerisSample

Stage-consistent Sun position expressed in the inertial frame.
"""
struct SolarEphemerisSample
    sun_pos_ii::SVector{3, Float64}
end

"""
    ThirdBodyEphemerisSample{N}

Stage-consistent third-body inertial positions for an N-body effector.
"""
struct ThirdBodyEphemerisSample{N}
    names::NTuple{N, String}
    positions_ii::NTuple{N, SVector{3, Float64}}
end

"""
    EnvironmentSample

Typed environment bundle passed to [`wrench`](@ref). `planet` carries the
typed static planet model used by the current simulation; the remaining fields
are optional sampled capabilities requested by the effector.
"""
struct EnvironmentSample{P, PF, AT, SE, TB}
    planet::P
    planet_frame::PF
    atmosphere::AT
    solar::SE
    third_bodies::TB
end

EnvironmentSample(
    planet=nothing;
    planet_frame=nothing,
    atmosphere=nothing,
    solar=nothing,
    third_bodies=nothing,
) = EnvironmentSample(planet, planet_frame, atmosphere, solar, third_bodies)

"""
    EffectorEnvironmentRequirements

Capability request describing which sampled environment fields should be built
for a [`wrench`](@ref) evaluation.
"""
struct EffectorEnvironmentRequirements{TB <: Tuple{Vararg{String}}}
    planet_frame::Bool
    atmosphere::Bool
    solar::Bool
    third_body_names::TB
end

EffectorEnvironmentRequirements(;
    planet_frame::Bool=false,
    atmosphere::Bool=false,
    solar::Bool=false,
    third_body_names::Tuple=(),
) = EffectorEnvironmentRequirements(planet_frame, atmosphere, solar, third_body_names)

"""
    environment_requirements(model) -> EffectorEnvironmentRequirements

Preferred additive declaration for the sampled environment data required by a
[`wrench`](@ref) implementation. The default requests no additional sampled
fields.
"""
@inline environment_requirements(::Any) = EffectorEnvironmentRequirements()

"""
    wrench(model, x::StateSample, env::EnvironmentSample, t::Float64)

Preferred additive force/torque extension hook for dynamic effectors.

Returns `(force_ii, torque_body)` in SI units, where force is expressed in the
inertial frame and torque is expressed in the body frame. Implementations should
behave as pure functions of `(model, x, env, t)`.
"""
function wrench end

end # module EffectorSampling
